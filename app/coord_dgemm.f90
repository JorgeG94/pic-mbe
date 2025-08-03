program mpi_dgemm_task_distributor
    use mpi_f08
    use iso_fortran_env, only: dp => real64, int64
    use pic_timer
    use pic_blas_interfaces
    implicit none

    ! Parameters
    type(pic_timer_type) :: my_timer, rank_timer
    integer, parameter :: total_tasks = 32
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
    call my_timer%start()
    endif
    
    if (rank == 0) then
        call static_coordinator(num_procs)
    else
        call rank_timer%start()
        call worker(rank)
        call rank_timer%stop()
        print *, "Rank", rank, "finished in", rank_timer%get_elapsed_time(), "seconds"
    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if(rank == 0) then 
    call my_timer%stop()
    call my_timer%print_time()
    endif

    call MPI_Finalize(ierr)

contains

subroutine static_coordinator(num_procs)
    integer, intent(in) :: num_procs
    integer :: task_id, i, tasks_assigned, target_rank
    type(MPI_Status) :: status

    tasks_assigned = 0

    ! Distribute all tasks statically
    do task_id = 0, total_tasks - 1
        ! Round-robin or block distribution: skip coordinator (rank 0)
        target_rank = mod(task_id, num_procs - 1) + 1  ! Ranks 1..(num_procs-1)
        call MPI_Send(task_id, 1, MPI_INTEGER, target_rank, tag_work, MPI_COMM_WORLD)
        tasks_assigned = tasks_assigned + 1
    end do

    ! Notify all workers that no more work is coming
    do i = 1, num_procs - 1
        call MPI_Send(-1, 1, MPI_INTEGER, i, tag_done, MPI_COMM_WORLD)
    end do

    print *, "Coordinator: all tasks assigned statically."
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
        integer :: m
        real(dp), allocatable :: A(:,:), B(:,:), C(:,:)
        type(pic_timer_type) :: timer
        real(dp) :: elapsed_time

        m = 4096
        call timer%start()
        allocate(A(m, m), B(m, m), C(m, m))
        A = 1.0_dp
        B = 1.0_dp
        C = 0.0_dp
        !print *, "Time to allocate and initialize matrices was: ", elapsed_time, "seconds"

        call pic_gemm(A,B,C)
        call timer%stop()
        elapsed_time = timer%get_elapsed_time()
        print *, "Time to perform DGEMM for task", task_id, "was: ", elapsed_time, "seconds"


        deallocate(A, B, C)
    end subroutine dummy_work

end program mpi_dgemm_task_distributor

