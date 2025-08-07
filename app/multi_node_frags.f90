program hierarchical_mpi_mbe
   use mpi_f08
   use pic_blas_interfaces, only: pic_gemm
   use pic_timer, only: pic_timer_type
   use pic_types, only: dp, default_int
   use pic_string_utils, only: to_string
   use pic_mbe, only: get_nfrags, create_monomer_list, generate_fragment_list
   use pic_fragment, only: pic_fragment_block, count_nonzeros
   implicit none

   ! Fragment generation parameters
   integer(default_int), parameter :: n_monomers = 15
   integer(default_int), parameter :: max_level = 3
   integer(default_int), parameter :: n = 128  ! monomer matrix size
   
   ! MPI and timing variables
   integer :: rank, size, ierr
   integer :: node_rank, node_size
   type(MPI_Comm) :: comm_node
   type(pic_timer_type) :: timer
   real(dp) :: elapsed_time, flops
   
   ! Fragment data structures
   integer(default_int), allocatable :: monomers(:)
   integer(default_int), allocatable :: polymers(:, :)
   integer(default_int) :: n_expected_fragments, fragment_count
   type(pic_fragment_block), allocatable :: fragments(:)
   
   ! Variables for node leader detection
   integer :: global_node_rank
   integer, allocatable :: all_node_leader_ranks(:), node_leader_ranks(:)
   integer :: i, num_nodes, j

   ! MPI Initialization
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

   call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_node, ierr)
   call MPI_Comm_rank(comm_node, node_rank, ierr)
   call MPI_Comm_size(comm_node, node_size, ierr)

   if (size < 3) then
      print *, "This program requires at least 3 processes."
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
   end if

   ! Determine global ranks of node leaders
   global_node_rank = -1
   if (node_rank == 0) global_node_rank = rank

   allocate (all_node_leader_ranks(size))
   call MPI_Allgather(global_node_rank, 1, MPI_INTEGER, all_node_leader_ranks, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

   ! Extract valid node leader ranks
   num_nodes = count(all_node_leader_ranks /= -1)
   allocate (node_leader_ranks(num_nodes))
   i = 0
   do concurrent(j=1:size)
      if (all_node_leader_ranks(j) /= -1) then
         i = i + 1
         node_leader_ranks(i) = all_node_leader_ranks(j)
      end if
   end do

   ! Generate fragments only on rank 0
   if (rank == 0) then
      n_expected_fragments = get_nfrags(n_monomers, max_level)
      allocate (monomers(n_monomers))
      allocate (polymers(n_expected_fragments, max_level))
      allocate (fragments(n_expected_fragments))
      
      print *, "Expected fragments = ", n_expected_fragments
      
      ! Initialize polymers array
      polymers = 0_default_int
      call create_monomer_list(monomers)
      polymers(:n_monomers, 1) = monomers(:n_monomers)
      fragment_count = n_monomers
      
      ! Generate all fragment combinations
      call generate_fragment_list(monomers, max_level, polymers, fragment_count)
      
      if(fragment_count /= n_expected_fragments) then
         print *, "Error: fragment_count does not match expected fragments!"
         call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if
      
      print *, "Generated ", fragment_count, " fragments for distribution"
      call timer%start()
   end if

   ! Role assignment
   if (rank == 0 .and. node_rank == 0) then
      ! Global coordinator (first node leader)
      call global_coordinator(rank, comm_node, fragment_count, polymers, max_level, node_leader_ranks, num_nodes)
   else if (node_rank == 0) then
      ! Node coordinator (other nodes)
      call node_coordinator(rank, comm_node, max_level)
   else
      ! Worker
      call node_worker(rank, comm_node, n, max_level)
   end if

   ! Timing and Flops Report
   if (rank == 0) then
      call timer%stop()
      elapsed_time = timer%get_elapsed_time()
      ! Estimate flops based on matrix operations per fragment
      !flops = real(fragment_count, dp) * 2.0_dp * (n**3)  ! Rough estimate
      flops = 0.0_dp
      call calculate_exact_flops(polymers, fragment_count, max_level, n, flops) 
      print *, "Total elapsed time for all fragments:", elapsed_time, "seconds"
      print *, "Total slops ", flops
      print *, "Estimated flop rate: ", flops / elapsed_time / 1.0e9_dp, " GFLOP/s"
   end if

   call MPI_Finalize(ierr)

contains
   subroutine calculate_exact_flops(polymers, fragment_count, max_level, matrix_size, total_flops)
      integer, intent(in) :: polymers(:,:), fragment_count, max_level, matrix_size
      real(dp), intent(out) :: total_flops
      integer :: i, fragment_size, n_monomers, n_dimers, n_trimers
      real(dp) :: monomer_flops, dimer_flops, trimer_flops
      real(dp) :: monomer_size, dimer_size, trimer_size 
      
      n_monomers = 0
      n_dimers = 0
      n_trimers = 0
     
      monomer_size = real(matrix_size,dp) 
      dimer_size = real(2 * matrix_size,dp) 
      trimer_size = real(3 * matrix_size,dp)
      ! Count fragments by size
      do i = 1, fragment_count
         fragment_size = count(polymers(i, :) > 0)
         select case (fragment_size)
         case (1)
            n_monomers = n_monomers + 1
         case (2)
            n_dimers = n_dimers + 1
         case (3)
            n_trimers = n_trimers + 1
         end select
      end do
      
      ! Calculate flops for each fragment type (assuming GEMM: 2*n^3 flops)
      monomer_flops = real(n_monomers, dp) * 2.0_dp * real(matrix_size*matrix_size*matrix_size,dp)
      dimer_flops = real(n_dimers, dp) * 2.0_dp * real(dimer_size * dimer_size * dimer_size,dp)
      trimer_flops = real(n_trimers, dp) * 2.0_dp * real(trimer_size * trimer_size * trimer_size, dp)
      
      total_flops = monomer_flops + dimer_flops + trimer_flops
      
      print *, "Fragment breakdown:"
      print *, "  Monomers: ", n_monomers, " (", monomer_flops/1.0e9_dp, " GFLOP)"
      print *, "  Dimers:   ", n_dimers, " (", dimer_flops/1.0e9_dp, " GFLOP)"
      print *, "  Trimers:  ", n_trimers, " (", trimer_flops/1.0e9_dp, " GFLOP)"
      print *, "  Total:    ", fragment_count, " (", total_flops/1.0e9_dp, " GFLOP)"
   end subroutine calculate_exact_flops

   ! Fragment work routine
   subroutine process_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size)
      integer, intent(in) :: fragment_idx, fragment_size, matrix_size
      integer, intent(in) :: fragment_indices(fragment_size)
      real(dp), allocatable :: A(:,:), B(:,:), C(:,:)
      integer :: i
      
      ! Allocate and initialize fragment matrix
      allocate(A(fragment_size * matrix_size, fragment_size * matrix_size))
      allocate(B(fragment_size * matrix_size, fragment_size * matrix_size))
      allocate(C(fragment_size * matrix_size, fragment_size * matrix_size))

      A = real(fragment_size * fragment_idx, dp)
      B = real(fragment_size * fragment_idx, dp)
      C = 0.0_dp
      
      call pic_gemm(A,B,C)
      
      deallocate(A, B, C)
   end subroutine process_fragment

   ! Global coordinator that distributes fragment indices to node leaders
   subroutine global_coordinator(rank, comm_node, total_fragments, polymers, max_level, node_leader_ranks, num_nodes)
      integer, intent(in) :: rank, total_fragments, max_level, num_nodes
      integer, intent(in) :: polymers(:,:), node_leader_ranks(:)
      type(MPI_Comm), intent(in) :: comm_node
      integer :: current_fragment, finished_nodes, ierr
      integer :: request_source, dummy_msg
      type(MPI_Status) :: status
      integer :: node_rank, node_size
      logical :: handling_local_workers
      
      ! For handling local workers
      integer :: local_finished_workers
      integer :: local_dummy, fragment_size
      integer, allocatable :: fragment_indices(:)
      type(MPI_Status) :: local_status
      logical :: local_message_pending

      call MPI_Comm_rank(comm_node, node_rank, ierr)
      call MPI_Comm_size(comm_node, node_size, ierr)

      current_fragment = 1  ! Start from fragment 1
      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (node_size > 1)

      print *, "Global coordinator starting with", total_fragments, "fragments for", num_nodes, "nodes"

      do while (finished_nodes < num_nodes)
         ! Check for requests from other node leaders
         call MPI_Iprobe(MPI_ANY_SOURCE, 300, MPI_COMM_WORLD, local_message_pending, status, ierr)

         if (local_message_pending) then
            call MPI_Recv(dummy_msg, 1, MPI_INTEGER, MPI_ANY_SOURCE, 300, MPI_COMM_WORLD, status, ierr)
            request_source = status%MPI_SOURCE

            if (current_fragment <= total_fragments) then
               ! Send fragment info to requesting node leader
               call send_fragment_to_node(current_fragment, polymers, max_level, request_source)
               current_fragment = current_fragment + 1
            else
               ! No more fragments - signal node to finish
               call MPI_Send(-1, 1, MPI_INTEGER, request_source, 302, MPI_COMM_WORLD, ierr)
               finished_nodes = finished_nodes + 1
            end if
         end if

         ! Handle local workers if this node has any
         if (handling_local_workers .and. local_finished_workers < node_size - 1) then
            call MPI_Iprobe(MPI_ANY_SOURCE, 200, comm_node, local_message_pending, local_status, ierr)

            if (local_message_pending) then
               call MPI_Recv(local_dummy, 1, MPI_INTEGER, MPI_ANY_SOURCE, 200, comm_node, local_status, ierr)

               if (current_fragment <= total_fragments) then
                  call send_fragment_to_worker(current_fragment, polymers, max_level, local_status%MPI_SOURCE, comm_node)
                  current_fragment = current_fragment + 1
               else
                  call MPI_Send(-1, 1, MPI_INTEGER, local_status%MPI_SOURCE, 202, comm_node, ierr)
                  local_finished_workers = local_finished_workers + 1
               end if
            end if
         end if

         ! If we're done with local workers, count this node as finished
         if (handling_local_workers .and. local_finished_workers >= node_size - 1) then
            if (finished_nodes < num_nodes) then
               finished_nodes = finished_nodes + 1
               handling_local_workers = .false.
               print *, "Global coordinator finished local workers"
            end if
         end if

         ! Small delay to prevent busy waiting
         call sleep(0)
      end do

      print *, "Global coordinator finished all fragments"
   end subroutine global_coordinator

   ! Helper subroutine to send fragment data to node coordinator
   subroutine send_fragment_to_node(fragment_idx, polymers, max_level, dest_rank)
      integer, intent(in) :: fragment_idx, max_level, dest_rank
      integer, intent(in) :: polymers(:,:)
      integer :: fragment_size, ierr
      integer, allocatable :: fragment_indices(:)
      
      ! Determine fragment size
      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate(fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)
      
      ! Send fragment index first
      call MPI_Send(fragment_idx, 1, MPI_INTEGER, dest_rank, 301, MPI_COMM_WORLD, ierr)
      ! Send fragment size
      call MPI_Send(fragment_size, 1, MPI_INTEGER, dest_rank, 301, MPI_COMM_WORLD, ierr)
      ! Send fragment indices
      call MPI_Send(fragment_indices, fragment_size, MPI_INTEGER, dest_rank, 301, MPI_COMM_WORLD, ierr)
      
      deallocate(fragment_indices)
   end subroutine send_fragment_to_node

   ! Helper subroutine to send fragment data to local worker
   subroutine send_fragment_to_worker(fragment_idx, polymers, max_level, dest_rank, comm)
      integer, intent(in) :: fragment_idx, max_level, dest_rank
      integer, intent(in) :: polymers(:,:)
      type(MPI_Comm), intent(in) :: comm
      integer :: fragment_size, ierr
      integer, allocatable :: fragment_indices(:)
      
      ! Determine fragment size
      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate(fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)
      
      ! Send fragment index first
      call MPI_Send(fragment_idx, 1, MPI_INTEGER, dest_rank, 201, comm, ierr)
      ! Send fragment size
      call MPI_Send(fragment_size, 1, MPI_INTEGER, dest_rank, 201, comm, ierr)
      ! Send fragment indices
      call MPI_Send(fragment_indices, fragment_size, MPI_INTEGER, dest_rank, 201, comm, ierr)
      
      deallocate(fragment_indices)
   end subroutine send_fragment_to_worker

   ! Node coordinator that requests fragments from global coordinator
   subroutine node_coordinator(rank, comm_node, max_level)
      integer, intent(in) :: rank, max_level
      type(MPI_Comm), intent(in) :: comm_node
      integer :: fragment_idx, fragment_size, ierr, dummy_msg
      integer :: node_rank, node_size, finished_workers
      integer, allocatable :: fragment_indices(:)
      type(MPI_Status) :: status, global_status
      logical :: local_message_pending, more_fragments
      integer :: local_dummy

      call MPI_Comm_rank(comm_node, node_rank, ierr)
      call MPI_Comm_size(comm_node, node_size, ierr)

      finished_workers = 0
      more_fragments = .true.
      dummy_msg = 0

      do while (finished_workers < node_size - 1)
         ! Check for local worker requests
         call MPI_Iprobe(MPI_ANY_SOURCE, 200, comm_node, local_message_pending, status, ierr)

         if (local_message_pending) then
            call MPI_Recv(local_dummy, 1, MPI_INTEGER, MPI_ANY_SOURCE, 200, comm_node, status, ierr)

            if (more_fragments) then
               ! Request fragment from global coordinator
               call MPI_Send(dummy_msg, 1, MPI_INTEGER, 0, 300, MPI_COMM_WORLD, ierr)
               call MPI_Recv(fragment_idx, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, global_status, ierr)

               if (global_status%MPI_TAG == 301) then
                  ! Got a fragment - receive the rest of the data
                  call MPI_Recv(fragment_size, 1, MPI_INTEGER, 0, 301, MPI_COMM_WORLD, global_status, ierr)
                  allocate(fragment_indices(fragment_size))
                  call MPI_Recv(fragment_indices, fragment_size, MPI_INTEGER, 0, 301, MPI_COMM_WORLD, global_status, ierr)
                  
                  ! Forward to worker
                  call MPI_Send(fragment_idx, 1, MPI_INTEGER, status%MPI_SOURCE, 201, comm_node, ierr)
                  call MPI_Send(fragment_size, 1, MPI_INTEGER, status%MPI_SOURCE, 201, comm_node, ierr)
                  call MPI_Send(fragment_indices, fragment_size, MPI_INTEGER, status%MPI_SOURCE, 201, comm_node, ierr)
                  
                  deallocate(fragment_indices)
               else
                  ! No more fragments
                  call MPI_Send(-1, 1, MPI_INTEGER, status%MPI_SOURCE, 202, comm_node, ierr)
                  finished_workers = finished_workers + 1
                  more_fragments = .false.
               end if
            else
               ! No more fragments available
               call MPI_Send(-1, 1, MPI_INTEGER, status%MPI_SOURCE, 202, comm_node, ierr)
               finished_workers = finished_workers + 1
            end if
         end if

         ! Small delay to prevent busy waiting
         call sleep(0)
      end do
   end subroutine node_coordinator

   ! Worker that processes fragments
   subroutine node_worker(rank, comm_node, matrix_size, max_level)
      integer, intent(in) :: rank, matrix_size, max_level
      type(MPI_Comm), intent(in) :: comm_node
      integer :: fragment_idx, fragment_size, ierr, dummy_msg
      integer, allocatable :: fragment_indices(:)
      type(MPI_Status) :: status

      dummy_msg = 0

      do
         call MPI_Send(dummy_msg, 1, MPI_INTEGER, 0, 200, comm_node, ierr)
         call MPI_Recv(fragment_idx, 1, MPI_INTEGER, 0, MPI_ANY_TAG, comm_node, status, ierr)

         select case (status%MPI_TAG)
         case (201)
            ! Receive fragment data
            call MPI_Recv(fragment_size, 1, MPI_INTEGER, 0, 201, comm_node, status, ierr)
            allocate(fragment_indices(fragment_size))
            call MPI_Recv(fragment_indices, fragment_size, MPI_INTEGER, 0, 201, comm_node, status, ierr)
            
            ! Process the fragment
            call process_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size)
            
            deallocate(fragment_indices)
         case (202)
            exit
         end select
      end do
   end subroutine node_worker

end program hierarchical_mpi_mbe
