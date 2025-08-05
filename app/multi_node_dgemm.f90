program hierarchical_mpi
   use mpi_f08
   use pic_blas_interfaces
   use pic_timer
   use pic_types, only: dp
   implicit none

   abstract interface
      subroutine work(task_id)
         integer, intent(in) :: task_id
      end subroutine work
   end interface

   procedure(work), pointer :: pwork => null()
   integer, parameter :: total_tasks = 4096
   integer, parameter :: m = 5192
   real(dp) :: flops
   integer :: rank, size, ierr
   integer :: node_rank, node_size
   type(MPI_Comm) :: comm_node
   type(pic_timer_type) :: timer
   real(dp) :: elapsed_time

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

   pwork => dummy_work

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

   if (rank == 0) then
      call timer%start()
   end if

   ! Role assignment
   if (rank == 0 .and. node_rank == 0) then
      ! Global coordinator (first node leader)
      call global_coordinator(rank, comm_node, total_tasks, node_leader_ranks, num_nodes)
   else if (node_rank == 0) then
      ! Node coordinator (other nodes)
      call node_coordinator(rank, comm_node)
   else
      ! Worker
      call node_worker(rank, comm_node, pwork)
   end if

   ! Timing and Flops Report
   if (rank == 0) then
      call timer%stop()
      elapsed_time = timer%get_elapsed_time()
      flops = total_tasks*2*m*m*m
      print *, "Total elapsed time for all tasks:", elapsed_time, "seconds"
      print *, "Flop rate is ", real(flops, dp)/elapsed_time/1.0e9_dp, " GFLOP/s"
   end if

   call MPI_Finalize(ierr)

contains

   subroutine dummy_work(task_id)
      integer, intent(in) :: task_id
      integer :: i, j, k
      real(dp), allocatable :: A(:, :), B(:, :), C(:, :)

      allocate (A(m, m), B(m, m), C(m, m))
      A = 1.0_dp
      B = 1.0_dp
      C = 0.0_dp

      call pic_gemm(A, B, C)

      deallocate (A, B, C)
   end subroutine dummy_work

! Global coordinator that distributes tasks dynamically to node leaders
   subroutine global_coordinator(rank, comm_node, total_tasks, node_leader_ranks, num_nodes)
      use mpi_f08
      implicit none
      integer, intent(in) :: rank, total_tasks
      integer, intent(in) :: node_leader_ranks(:), num_nodes
      type(MPI_Comm), intent(in) :: comm_node
      integer :: current_task, finished_nodes, ierr
      integer :: request_source, dummy_msg
      type(MPI_Status) :: status
      integer :: node_rank, node_size
      logical :: handling_local_workers

      ! For handling local workers
      integer :: local_task_queue, local_finished_workers
      integer :: local_request_source, local_dummy
      type(MPI_Status) :: local_status
      logical :: local_message_pending

      call MPI_Comm_rank(comm_node, node_rank, ierr)
      call MPI_Comm_size(comm_node, node_size, ierr)

      current_task = 0
      finished_nodes = 0
      local_task_queue = 0
      local_finished_workers = 0
      handling_local_workers = (node_size > 1)

      print *, "Global coordinator starting with", total_tasks, "tasks for", num_nodes, "nodes"

      do while (finished_nodes < num_nodes)
         ! Check for requests from other node leaders
         call MPI_Iprobe(MPI_ANY_SOURCE, 300, MPI_COMM_WORLD, local_message_pending, status, ierr)

         if (local_message_pending) then
            call MPI_Recv(dummy_msg, 1, MPI_INTEGER, MPI_ANY_SOURCE, 300, MPI_COMM_WORLD, status, ierr)
            request_source = status%MPI_SOURCE

            if (current_task < total_tasks) then
               ! Send next task to requesting node leader
               call MPI_Send(current_task, 1, MPI_INTEGER, request_source, 301, MPI_COMM_WORLD, ierr)
               current_task = current_task + 1
            else
               ! No more tasks - signal node to finish
               call MPI_Send(-1, 1, MPI_INTEGER, request_source, 302, MPI_COMM_WORLD, ierr)
               finished_nodes = finished_nodes + 1
            end if
         end if

         ! Handle local workers if this node has any
         if (handling_local_workers .and. local_finished_workers < node_size - 1) then
            call MPI_Iprobe(MPI_ANY_SOURCE, 200, comm_node, local_message_pending, local_status, ierr)

            if (local_message_pending) then
               call MPI_Recv(local_dummy, 1, MPI_INTEGER, MPI_ANY_SOURCE, 200, comm_node, local_status, ierr)

               if (current_task < total_tasks) then
                  call MPI_Send(current_task, 1, MPI_INTEGER, local_status%MPI_SOURCE, 201, comm_node, ierr)
                  current_task = current_task + 1
               else
                  call MPI_Send(-1, 1, MPI_INTEGER, local_status%MPI_SOURCE, 202, comm_node, ierr)
                  local_finished_workers = local_finished_workers + 1
               end if
            end if
         end if

         ! If we're done with local workers, count this node as finished
         if (handling_local_workers .and. local_finished_workers >= node_size - 1 .and. finished_nodes < num_nodes) then
            finished_nodes = finished_nodes + 1
            handling_local_workers = .false.
            print *, "Global coordinator finished local workers"
         end if

         ! Small delay to prevent busy waiting
         call sleep(0)
      end do

      print *, "Global coordinator finished all tasks"
   end subroutine global_coordinator

! Node coordinator that requests tasks from global coordinator
   subroutine node_coordinator(rank, comm_node)
      use mpi_f08
      implicit none
      integer, intent(in) :: rank
      type(MPI_Comm), intent(in) :: comm_node
      integer :: task_id, ierr, dummy_msg
      integer :: node_rank, node_size
      integer :: finished_workers
      type(MPI_Status) :: status, global_status
      logical :: local_message_pending, more_tasks
      integer :: local_dummy

      call MPI_Comm_rank(comm_node, node_rank, ierr)
      call MPI_Comm_size(comm_node, node_size, ierr)

      finished_workers = 0
      more_tasks = .true.
      dummy_msg = 0

      !print *, "Node coordinator rank", rank, "starting"

      do while (finished_workers < node_size - 1)
         ! Check for local worker requests
         call MPI_Iprobe(MPI_ANY_SOURCE, 200, comm_node, local_message_pending, status, ierr)

         if (local_message_pending) then
            call MPI_Recv(local_dummy, 1, MPI_INTEGER, MPI_ANY_SOURCE, 200, comm_node, status, ierr)

            if (more_tasks) then
               ! Request task from global coordinator
               call MPI_Send(dummy_msg, 1, MPI_INTEGER, 0, 300, MPI_COMM_WORLD, ierr)
               call MPI_Recv(task_id, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, global_status, ierr)

               if (global_status%MPI_TAG == 301) then
                  ! Got a task - forward to worker
                  call MPI_Send(task_id, 1, MPI_INTEGER, status%MPI_SOURCE, 201, comm_node, ierr)
               else
                  ! No more tasks
                  call MPI_Send(-1, 1, MPI_INTEGER, status%MPI_SOURCE, 202, comm_node, ierr)
                  finished_workers = finished_workers + 1
                  more_tasks = .false.
               end if
            else
               ! No more tasks available
               call MPI_Send(-1, 1, MPI_INTEGER, status%MPI_SOURCE, 202, comm_node, ierr)
               finished_workers = finished_workers + 1
            end if
         end if

         ! Small delay to prevent busy waiting
         call sleep(0)
      end do

      !print *, "Node coordinator rank", rank, "finished"
   end subroutine node_coordinator

   subroutine node_worker(rank, comm_node, work_routine)
      use mpi_f08
      implicit none
      integer, intent(in) :: rank
      type(MPI_Comm), intent(in) :: comm_node
      integer :: task_id, ierr, dummy_msg
      type(MPI_Status) :: status
      procedure(dummy_work) :: work_routine

      !print *, "Worker rank", rank, "starting"
      dummy_msg = 0

      do
         call MPI_Send(dummy_msg, 1, MPI_INTEGER, 0, 200, comm_node, ierr)
         call MPI_Recv(task_id, 1, MPI_INTEGER, 0, MPI_ANY_TAG, comm_node, status, ierr)

         select case (status%MPI_TAG)
         case (201)
            !print *, "Worker rank", rank, "processing task", task_id
            !call dummy_work(task_id)
            call work_routine(task_id)
         case (202)
            !print *, "Worker rank", rank, "done"
            exit
         end select
      end do

      !print *, "Worker rank", rank, "finished"
   end subroutine node_worker

end program hierarchical_mpi
