module pic_mpi_algorithms
  use mpi_f08
  use mpi_comm_simple
  use pic_types
   use pic_blas_interfaces, only: pic_gemm
   use pic_timer
  implicit none 

  contains 

     subroutine process_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size)
      integer, intent(in) :: fragment_idx, fragment_size, matrix_size
      integer, intent(in) :: fragment_indices(fragment_size)
      real(dp), allocatable :: A(:,:), B(:,:), C(:,:)
      integer :: i
      type(pic_timer_type) :: gemm_timer
      real(dp) :: elapsed_time
      
      ! Allocate and initialize fragment matrix
      allocate(A(fragment_size * matrix_size, fragment_size * matrix_size))
      allocate(B(fragment_size * matrix_size, fragment_size * matrix_size))
      allocate(C(fragment_size * matrix_size, fragment_size * matrix_size))


      A = real(fragment_size * fragment_idx, dp)
      B = real(fragment_size * fragment_idx, dp)
      C = 0.0_dp
      
      call gemm_timer%start()
      call pic_gemm(A,B,C)
      call gemm_timer%stop()
      elapsed_time = gemm_timer%get_elapsed_time()

      print *, "Gemm for fragment", fragment_indices, " was ", elapsed_time, " seconds"
      
      deallocate(A, B, C)
   end subroutine process_fragment

  subroutine calculate_exact_flops(polymers, fragment_count, max_level, matrix_size, total_flops)
   use pic_types, only: dp
   implicit none
   integer, intent(in) :: polymers(:,:), fragment_count, max_level, matrix_size
   real(dp), intent(out) :: total_flops
   integer :: i, fragment_size, n_monomers, n_dimers, n_trimers
   real(dp) :: monomer_flops, dimer_flops, trimer_flops
   real(dp) :: monomer_size, dimer_size, trimer_size 

   n_monomers = 0
   n_dimers   = 0
   n_trimers  = 0

   monomer_size = real(matrix_size, dp)
   dimer_size   = real(2 * matrix_size, dp)
   trimer_size  = real(3 * matrix_size, dp)

   do i = 1, fragment_count
      fragment_size = count(polymers(i, :) > 0)
      select case (fragment_size)
      case (1); n_monomers = n_monomers + 1
      case (2); n_dimers   = n_dimers   + 1
      case (3); n_trimers  = n_trimers  + 1
      end select
   end do

   monomer_flops = real(n_monomers, dp) * 2.0_dp * monomer_size**3
   dimer_flops   = real(n_dimers, dp)   * 2.0_dp * dimer_size**3
   trimer_flops  = real(n_trimers, dp)  * 2.0_dp * trimer_size**3

   total_flops = monomer_flops + dimer_flops + trimer_flops

   print *, "Fragment breakdown:"
   print *, "  Monomers: ", n_monomers, " (", monomer_flops / 1.0e9_dp, " GFLOP)"
   print *, "  Dimers:   ", n_dimers,   " (", dimer_flops   / 1.0e9_dp, " GFLOP)"
   print *, "  Trimers:  ", n_trimers,  " (", trimer_flops  / 1.0e9_dp, " GFLOP)"
   print *, "  Total:    ", fragment_count, " (", total_flops / 1.0e9_dp, " GFLOP)"
end subroutine calculate_exact_flops

subroutine global_coordinator(world_comm, node_comm, total_fragments, polymers, max_level, node_leader_ranks, num_nodes)
   type(comm_t), intent(in) :: world_comm, node_comm
   integer, intent(in) :: total_fragments, max_level, num_nodes
   integer, intent(in) :: polymers(:,:), node_leader_ranks(:)

   integer :: current_fragment, finished_nodes
   integer :: request_source, dummy_msg
   type(MPI_Status) :: status, local_status
   logical :: handling_local_workers
   logical :: has_pending

   ! For local workers
   integer :: local_finished_workers, fragment_size, local_dummy
   integer, allocatable :: fragment_indices(:)

   current_fragment = 1
   finished_nodes = 0
   local_finished_workers = 0
   handling_local_workers = (node_comm%size() > 1)

   print *, "Global coordinator starting with", total_fragments, "fragments for", num_nodes, "nodes"

   do while (finished_nodes < num_nodes)

      ! Remote node coordinator requests
      call iprobe(world_comm, MPI_ANY_SOURCE, 300, has_pending, status)
      if (has_pending) then
         call recv(world_comm, dummy_msg, status%MPI_SOURCE, 300)
         request_source = status%MPI_SOURCE

         if (current_fragment <= total_fragments) then
            call send_fragment_to_node(world_comm, current_fragment, polymers, max_level, request_source)
            current_fragment = current_fragment + 1
         else
            call send(world_comm, -1, request_source, 302)
            finished_nodes = finished_nodes + 1
         end if
      end if

      ! Local workers (shared memory)
      if (handling_local_workers .and. local_finished_workers < node_comm%size() - 1) then
         call iprobe(node_comm, MPI_ANY_SOURCE, 200, has_pending, local_status)
         if (has_pending) then
            call recv(node_comm, local_dummy, local_status%MPI_SOURCE, 200)

            if (current_fragment <= total_fragments) then
               call send_fragment_to_worker(node_comm, current_fragment, polymers, max_level, local_status%MPI_SOURCE)
               current_fragment = current_fragment + 1
            else
               call send(node_comm, -1, local_status%MPI_SOURCE, 202)
               local_finished_workers = local_finished_workers + 1
            end if
         end if
      end if

      ! Finalize local worker completion
      if (handling_local_workers .and. local_finished_workers >= node_comm%size() - 1) then
         handling_local_workers = .false.
         if (num_nodes == 1) then
            finished_nodes = finished_nodes + 1
            print *, "Manually incremented finished_nodes for self"
         else
            finished_nodes = finished_nodes + 1
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


subroutine node_coordinator(world_comm, node_comm, max_level)
   class(comm_t), intent(in) :: world_comm, node_comm
   integer(int32), intent(in) :: max_level

   integer(int32) :: fragment_idx, fragment_size, ierr, dummy_msg
   integer(int32) :: finished_workers
   integer(int32), allocatable :: fragment_indices(:)
   type(MPI_Status) :: status, global_status
   logical :: local_message_pending, more_fragments
   integer(int32) :: local_dummy

   finished_workers = 0
   more_fragments = .true.
   dummy_msg = 0

   do while (finished_workers < node_comm%size() - 1)
      call iprobe(node_comm, MPI_ANY_SOURCE, 200, local_message_pending, status)

      if (local_message_pending) then
         call recv(node_comm, local_dummy, MPI_ANY_SOURCE, 200)

         if (more_fragments) then
            call send(world_comm, dummy_msg, 0, 300)
            call recv(world_comm, fragment_idx, 0, MPI_ANY_TAG, global_status)

            if (global_status%MPI_TAG == 301) then
               call recv(world_comm, fragment_size, 0, 301, global_status)
               allocate(fragment_indices(fragment_size))
               call recv(world_comm, fragment_indices, 0, 301, global_status)

               call send(node_comm, fragment_idx, status%MPI_SOURCE, 201)
               call send(node_comm, fragment_size, status%MPI_SOURCE, 201)
               call send(node_comm, fragment_indices, status%MPI_SOURCE, 201)

               deallocate(fragment_indices)
            else
               call send(node_comm, -1, status%MPI_SOURCE, 202)
               finished_workers = finished_workers + 1
               more_fragments = .false.
            end if
         else
            call send(node_comm, -1, status%MPI_SOURCE, 202)
            finished_workers = finished_workers + 1
         end if
      end if

      call sleep(0)
   end do
end subroutine node_coordinator


subroutine node_worker(world_comm, node_comm, matrix_size, max_level)

   class(comm_t), intent(in) :: world_comm, node_comm
   integer, intent(in) :: matrix_size, max_level

   integer(int32) :: fragment_idx, fragment_size, ierr, dummy_msg
   integer(int32), allocatable :: fragment_indices(:)
   type(MPI_Status) :: status

   dummy_msg = 0

   do
      call send(node_comm, dummy_msg, 0, 200)
      call recv(node_comm, fragment_idx, 0, MPI_ANY_TAG, status)

      select case (status%MPI_TAG)
      case (201)
         call recv(node_comm, fragment_size, 0, 201, status)
         allocate(fragment_indices(fragment_size))
         call recv(node_comm, fragment_indices, 0, 201, status)

         call process_fragment(fragment_idx, fragment_indices, fragment_size, matrix_size)

         deallocate(fragment_indices)
      case (202)
         exit
      end select
   end do
end subroutine node_worker



end module pic_mpi_algorithms
