program mpi_dgemm_task_distributor
    use mpi_f08
    use iso_fortran_env, only: dp => real64, int64
    use pic_timer
    use pic_blas_interfaces
    use pic_flop_rate
    implicit none

    ! Parameters
    type(pic_timer_type) :: rank_timer, coord_timer
    type(flop_rate_type) :: my_flop_rate
    ! size of maitrx
    integer, parameter :: m = 4096
    integer(int64) :: flops
    real(dp) :: elapsed_time
    integer, parameter :: total_tasks = 1024
    integer, parameter :: num_initial = 4
    integer, parameter :: tag_request = 1, tag_work = 2, tag_done = 3

    ! MPI variables
    integer :: rank, num_procs, ierr

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)

    if (num_procs < 2) then
        if (rank == 0) print *, "Need at least 2 MPI processes."
        call MPI_Finalize(ierr)
        stop
    end if

    if(rank == 0) then 
    call my_flop_rate%start_time()
    endif
    
    if (rank == 0) then
        
        call coord_timer%start()
        call static_coordinator(num_procs)
        call coord_timer%stop()
        print *, "Coordinator finished in", coord_timer%get_elapsed_time(), "seconds"
    else
        call rank_timer%start()
        call static_worker(rank, num_procs)
        call rank_timer%stop()
        print *, "Rank", rank, "finished in", rank_timer%get_elapsed_time(), "seconds"
    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if(rank == 0) then 
    flops = total_tasks * 2 * m * m * m  
    call my_flop_rate%add_flops(flops)
    call my_flop_rate%stop_time()
    elapsed_time = my_flop_rate%get_time()
    print *, "Elapsed time ", elapsed_time
    call my_flop_rate%report()
    endif

    call MPI_Finalize(ierr)

contains

subroutine static_coordinator(num_procs)
    integer, intent(in) :: num_procs
    integer :: task_id, i, start_task, end_task
    integer :: tasks_per_worker, extra_tasks
    integer :: rank

    tasks_per_worker = total_tasks / (num_procs - 1)
    extra_tasks = mod(total_tasks, num_procs - 1)
    task_id = 0

    do rank = 1, num_procs - 1
        ! Compute range of tasks for this rank
        start_task = task_id
        end_task = start_task + tasks_per_worker - 1
        if (rank <= extra_tasks) end_task = end_task + 1

        ! Send tasks [start_task .. end_task] to this rank
        do i = start_task, end_task
            call MPI_Send(i, 1, MPI_INTEGER, rank, tag_work, MPI_COMM_WORLD)
        end do

        task_id = end_task + 1
    end do

    print *, "Coordinator: all tasks statically assigned."
end subroutine static_coordinator


    subroutine static_dynamic_coordinator(num_procs)
        integer, intent(in) :: num_procs
        integer :: task_id, i, j, tasks_assigned, workers_done
        integer :: source_rank
        logical :: done(num_procs)
        type(MPI_Status) :: status

        done = .false.
        done(1) = .true.  ! Coordinator itself is done
        tasks_assigned = 0
        workers_done = 1  ! Start from 1 to include coordinator

        ! Static distribution
        do i = 1, num_procs - 1
            do j = 1, num_initial
                if (tasks_assigned < total_tasks) then
                    call MPI_Send(tasks_assigned, 1, MPI_INTEGER, i, tag_work, MPI_COMM_WORLD)
                    tasks_assigned = tasks_assigned + 1
                end if
            end do
        end do

        ! Dynamic distribution
        do while (workers_done < num_procs)
            call MPI_Recv(task_id, 1, MPI_INTEGER, MPI_ANY_SOURCE, tag_request, MPI_COMM_WORLD, status)
            source_rank = status%MPI_SOURCE

            if (tasks_assigned < total_tasks) then
                call MPI_Send(tasks_assigned, 1, MPI_INTEGER, source_rank, tag_work, MPI_COMM_WORLD)
                tasks_assigned = tasks_assigned + 1
            else
                call MPI_Send(-1, 1, MPI_INTEGER, source_rank, tag_done, MPI_COMM_WORLD)
                if (.not. done(source_rank + 1)) then
                    done(source_rank + 1) = .true.
                    workers_done = count(done)
                end if
            end if
        end do

        print *, "Coordinator: all tasks completed and workers done."
    end subroutine static_dynamic_coordinator

subroutine static_worker(rank, num_procs)
    integer, intent(in) :: rank, num_procs
    integer :: task_id, num_tasks, i
        type(MPI_Status) :: status

    ! Compute my share of tasks
    num_tasks = total_tasks / (num_procs - 1)
    if (rank <= mod(total_tasks, num_procs - 1)) num_tasks = num_tasks + 1

    do i = 1, num_tasks
        call MPI_Recv(task_id, 1, MPI_INTEGER, 0, tag_work, MPI_COMM_WORLD, status)
        call dummy_work(task_id)
    end do
end subroutine static_worker


    subroutine worker(rank)
        integer, intent(in) :: rank
        integer :: task_id, ierr
        type(MPI_Status) :: status

        do
            call MPI_Recv(task_id, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status)

            select case (status%MPI_TAG)
            case (tag_work)
            !    print *, "Rank", rank, "processing task", task_id
                call dummy_work(task_id)
                call MPI_Send(0, 1, MPI_INTEGER, 0, tag_request, MPI_COMM_WORLD)
            case (tag_done)
                print *, "Rank", rank, "done."
                exit
            end select
        end do
    end subroutine worker

    subroutine dummy_work(task_id)
        integer, intent(in) :: task_id
        integer(int64) :: i, j, k
        real(dp), allocatable :: A(:,:), B(:,:), C(:,:)
        type(pic_timer_type) :: timer
        real(dp) :: elapsed_time

        !call timer%start()
        allocate(A(m, m), B(m, m), C(m, m))
        A = 1.0_dp
        B = 1.0_dp
        C = 0.0_dp
        !print *, "Time to allocate and initialize matrices was: ", elapsed_time, "seconds"

        call pic_gemm(A,B,C)
        !call timer%stop()
        !elapsed_time = timer%get_elapsed_time()
        !print *, "Time to perform DGEMM for task", task_id, "was: ", elapsed_time, "seconds"


        deallocate(A, B, C)
    end subroutine dummy_work

end program mpi_dgemm_task_distributor

