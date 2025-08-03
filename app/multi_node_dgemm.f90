program hierarchical_mpi
    use mpi_f08
    use pic_blas_interfaces
    use pic_timer
    use pic_types, only: dp
    implicit none

    integer, parameter :: total_tasks = 128
    integer :: rank, size, ierr
    integer :: node_rank, node_size
    type(MPI_Comm) :: comm_node
    type(pic_timer_type) :: timer
    real(dp) :: elapsed_time

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_node, ierr)
    call MPI_Comm_rank(comm_node, node_rank, ierr)
    call MPI_Comm_size(comm_node, node_size, ierr)
    if(size < 3) then
        print *, "This program requires at least 3 processes."
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if 

    if(rank == 0) then
        call timer%start()
    end if

if (rank == 0) then
    if (node_rank == 0) then
        ! Rank 0 is global and node coordinator (first rank on this node)
        call global_and_node_coordinator(rank, comm_node, total_tasks)
    else
        call program_coordinator(rank, total_tasks)
    end if
else if (node_rank == 0) then
    call node_coordinator(rank, comm_node)
else
    call node_worker(rank, comm_node)
end if

    if(rank == 0) then
        call timer%stop()
        elapsed_time = timer%get_elapsed_time()
        print *, "Total elapsed time for all tasks:", elapsed_time, "seconds"
    end if

    call MPI_Finalize(ierr)

    contains 

        subroutine dummy_work(task_id)
        integer, intent(in) :: task_id
        integer :: i, j, k
        integer :: m
        real(dp), allocatable :: A(:,:), B(:,:), C(:,:)

        m = 4096
        allocate(A(m, m), B(m, m), C(m, m))
        A = 1.0_dp
        B = 1.0_dp
        C = 0.0_dp


        call pic_gemm(A,B,C)

        deallocate(A, B, C)

    end subroutine dummy_work

    subroutine program_coordinator(rank, total_tasks)
    use mpi_f08
    implicit none
    integer, intent(in) :: rank, total_tasks
    integer :: size, ierr, i, task_range(2)
    integer :: tasks_per_node, remaining, start_task, end_task
    integer :: node_count

    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    node_count = size - 1  ! exclude rank 0

    tasks_per_node = total_tasks / node_count
    remaining = mod(total_tasks, node_count)
    start_task = 0

    do i = 1, node_count
        end_task = start_task + tasks_per_node - 1
        if (i <= remaining) end_task = end_task + 1

        task_range = [start_task, end_task]
        call MPI_Send(task_range, 2, MPI_INTEGER, i, 100, MPI_COMM_WORLD, ierr)
        start_task = end_task + 1
    end do

    !print *, "Global coordinator: all tasks distributed."
end subroutine program_coordinator

subroutine global_and_node_coordinator(rank, comm_node, total_tasks)
    use mpi_f08
    implicit none
    integer, intent(in) :: rank, total_tasks
    type(MPI_Comm), intent(in) :: comm_node
    integer :: task_range(2)

    task_range = [0, total_tasks - 1]
    !print *, "Rank", rank, "is global AND node coordinator."

    call node_coordinator_worker_logic(rank, comm_node, task_range)
end subroutine global_and_node_coordinator

subroutine node_coordinator(rank, comm_node)
    use mpi_f08
    implicit none
    integer, intent(in) :: rank
    type(MPI_Comm), intent(in) :: comm_node
    integer :: task_range(2), ierr
    type(MPI_Status) :: status

    call MPI_Recv(task_range, 2, MPI_INTEGER, 0, 100, MPI_COMM_WORLD, status, ierr)

    !print *, "Node coordinator rank", rank, "got task range", task_range(1), "-", task_range(2)

    call node_coordinator_worker_logic(rank, comm_node, task_range)
end subroutine node_coordinator

subroutine node_coordinator_worker_logic(rank, comm_node, task_range)
    use mpi_f08
    implicit none
    integer, intent(in) :: rank
    integer, intent(in) :: task_range(2)
    type(MPI_Comm), intent(in) :: comm_node
    integer :: node_rank, node_size, ierr
    integer :: task_id, task_idx, finished_workers
    type(MPI_Status) :: status

    call MPI_Comm_rank(comm_node, node_rank, ierr)
    call MPI_Comm_size(comm_node, node_size, ierr)

    task_idx = task_range(1)
    finished_workers = 0

    do while (finished_workers < node_size - 1)
        call MPI_Recv(task_id, 1, MPI_INTEGER, MPI_ANY_SOURCE, 200, comm_node, status, ierr)

        if (task_idx <= task_range(2)) then
            call MPI_Send(task_idx, 1, MPI_INTEGER, status%MPI_SOURCE, 201, comm_node, ierr)
            task_idx = task_idx + 1
        else
            call MPI_Send(-1, 1, MPI_INTEGER, status%MPI_SOURCE, 202, comm_node, ierr)
            finished_workers = finished_workers + 1
        end if
    end do

    !print *, "Node coordinator rank", rank, "finished all local tasks."
end subroutine node_coordinator_worker_logic

subroutine node_worker(rank, comm_node)
    use mpi_f08
    implicit none
    integer, intent(in) :: rank
    type(MPI_Comm), intent(in) :: comm_node
    integer :: task_id, ierr
    type(MPI_Status) :: status

    print *, "Worker rank", rank, "starting."

    do
        call MPI_Send(0, 1, MPI_INTEGER, 0, 200, comm_node, ierr)
        call MPI_Recv(task_id, 1, MPI_INTEGER, 0, MPI_ANY_TAG, comm_node, status, ierr)

        select case (status%MPI_TAG)
        case (201)
            !print *, "Worker rank", rank, "processing task", task_id
            call dummy_work(task_id)
        case (202)
            !print *, "Worker rank", rank, "done."
            exit
        end select
    end do
end subroutine node_worker

end program hierarchical_mpi