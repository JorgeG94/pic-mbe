program hierarchical_mpi_mbe
   use mpi_f08
   use mpi_comm_simple  ! Your wrapper
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
   integer(default_int), parameter :: n = 1024  ! monomer matrix size
   
   ! MPI and timing variables using wrapper
   type(comm_t) :: world_comm, node_comm
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

   ! MPI Initialization using wrapper
   call MPI_Init()
   world_comm = comm_world()
   
   ! Create node communicator
   node_comm = world_comm%split()  ! Use shared memory split
   
   if (world_comm%size() < 3) then
      print *, "This program requires at least 3 processes."
      call world_comm%abort(1)
   end if

   ! Determine global ranks of node leaders using wrapper
   global_node_rank = -1
   if (node_comm%leader()) global_node_rank = world_comm%rank()

   allocate (all_node_leader_ranks(world_comm%size()))
   call MPI_Allgather(global_node_rank, 1, MPI_INTEGER, all_node_leader_ranks, 1, MPI_INTEGER, world_comm%get())

   ! Extract valid node leader ranks
   num_nodes = count(all_node_leader_ranks /= -1)
   allocate (node_leader_ranks(num_nodes))
   i = 0
   do concurrent(j=1:world_comm%size())
      if (all_node_leader_ranks(j) /= -1) then
         i = i + 1
         node_leader_ranks(i) = all_node_leader_ranks(j)
      end if
   end do

   ! Generate fragments only on rank 0
   if (world_comm%leader()) then
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
         call world_comm%abort(1)
      end if
      
      print *, "Generated ", fragment_count, " fragments for distribution"
      call timer%start()
   end if

   ! Role assignment using wrapper methods
   if (world_comm%leader() .and. node_comm%leader()) then
      ! Global coordinator (first node leader)
      call global_coordinator(world_comm, node_comm, fragment_count, polymers, max_level, node_leader_ranks, num_nodes)
   else if (node_comm%leader()) then
      ! Node coordinator (other nodes)
      call node_coordinator(world_comm, node_comm, max_level)
   else
      ! Worker
      call node_worker(world_comm, node_comm, n, max_level)
   end if

   ! Timing and Flops Report
   if (world_comm%leader()) then
      call timer%stop()
      elapsed_time = timer%get_elapsed_time()
      flops = 0.0_dp
      call calculate_exact_flops(polymers, fragment_count, max_level, n, flops) 
      print *, "Total elapsed time for all fragments:", elapsed_time, "seconds"
      print *, "Total flops ", flops
      print *, "Estimated flop rate: ", flops / elapsed_time / 1.0e9_dp, " GFLOP/s"
   end if

   ! Cleanup using wrapper
   call node_comm%finalize()
   call world_comm%finalize()
   call MPI_Finalize()

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

      subroutine global_coordinator(world_comm, node_comm, total_fragments, polymers, max_level, node_leader_ranks, num_nodes)
      type(comm_t), intent(in) :: world_comm, node_comm
      integer, intent(in) :: total_fragments, max_level, num_nodes
      integer, intent(in) :: polymers(:,:), node_leader_ranks(:)
      integer :: current_fragment, finished_nodes
      integer :: request_source, dummy_msg
      type(MPI_Status) :: status
      logical :: handling_local_workers
      
      ! For handling local workers
      integer :: local_finished_workers
      integer :: local_dummy, fragment_size
      integer, allocatable :: fragment_indices(:)
      type(MPI_Status) :: local_status
      logical :: local_message_pending

      current_fragment = 1
      finished_nodes = 0
      local_finished_workers = 0
      handling_local_workers = (node_comm%size() > 1)

      print *, "Global coordinator starting with", total_fragments, "fragments for", num_nodes, "nodes"

      do while (finished_nodes < num_nodes)
         ! Check for requests from other node leaders using wrapper
         call iprobe(world_comm, MPI_ANY_SOURCE, 300, local_message_pending, status)

         if (local_message_pending) then
            call recv(world_comm, dummy_msg, MPI_ANY_SOURCE, 300)
            request_source = status%MPI_SOURCE

            if (current_fragment <= total_fragments) then
               call send_fragment_to_node(world_comm, current_fragment, polymers, max_level, request_source)
               current_fragment = current_fragment + 1
            else
               call send(world_comm, -1, request_source, 302)
               finished_nodes = finished_nodes + 1
            end if
         end if

         ! Handle local workers
         if (handling_local_workers .and. local_finished_workers < node_comm%size() - 1) then
            call iprobe(node_comm, MPI_ANY_SOURCE, 200, local_message_pending, local_status)

            if (local_message_pending) then
               call recv(node_comm, local_dummy, MPI_ANY_SOURCE, 200)

               if (current_fragment <= total_fragments) then
                  call send_fragment_to_worker(node_comm, current_fragment, polymers, max_level, local_status%MPI_SOURCE)
                  current_fragment = current_fragment + 1
               else
                  call send(node_comm, -1, local_status%MPI_SOURCE, 202)
                  local_finished_workers = local_finished_workers + 1
               end if
            end if
         end if

         ! Check if done with local workers
         if (handling_local_workers .and. local_finished_workers >= node_comm%size() - 1) then
            if (finished_nodes < num_nodes) then
               finished_nodes = finished_nodes + 1
               handling_local_workers = .false.
               print *, "Global coordinator finished local workers"
            end if
         end if

         call sleep(0)
      end do

      print *, "Global coordinator finished all fragments"
   end subroutine global_coordinator

     subroutine send_fragment_to_node(world_comm, fragment_idx, polymers, max_level, dest_rank)
      type(comm_t), intent(in) :: world_comm
      integer, intent(in) :: fragment_idx, max_level, dest_rank
      integer, intent(in) :: polymers(:,:)
      integer :: fragment_size
      integer, allocatable :: fragment_indices(:)
      
      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate(fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)
      
      call send(world_comm, fragment_idx, dest_rank, 301)
      call send(world_comm, fragment_size, dest_rank, 301)
      call send(world_comm, fragment_indices, dest_rank, 301)
      
      deallocate(fragment_indices)
   end subroutine send_fragment_to_node

   subroutine send_fragment_to_worker(node_comm, fragment_idx, polymers, max_level, dest_rank)
      type(comm_t), intent(in) :: node_comm
      integer, intent(in) :: fragment_idx, max_level, dest_rank
      integer, intent(in) :: polymers(:,:)
      integer :: fragment_size
      integer, allocatable :: fragment_indices(:)
      
      fragment_size = count(polymers(fragment_idx, :) > 0)
      allocate(fragment_indices(fragment_size))
      fragment_indices = polymers(fragment_idx, 1:fragment_size)
      
      call send(node_comm, fragment_idx, dest_rank, 201)
      call send(node_comm, fragment_size, dest_rank, 201)  
      call send(node_comm, fragment_indices, dest_rank, 201)
      
      deallocate(fragment_indices)
   end subroutine send_fragment_to_worker

end program hierarchical_mpi_mbe
