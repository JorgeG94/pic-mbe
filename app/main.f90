program hierarchical_mpi_mbe
   use mpi_f08
   use mpi_comm_simple        
   use pic_blas_interfaces, only: pic_gemm
   use pic_timer, only: timer_type
   use pic_types, only: dp, default_int
   use pic_io, only: to_char
   use pic_mbe, only: get_nfrags, create_monomer_list, generate_fragment_list
   use pic_fragment, only: pic_fragment_block, count_nonzeros
   use pic_mpi_algorithms
   implicit none

   ! Fragment generation parameters
   integer(default_int), parameter :: n_monomers = 50
   integer(default_int), parameter :: max_level = 3
   integer(default_int), parameter :: n = 256  ! monomer matrix size

   ! MPI wrappers
   type(comm_t) :: world_comm, node_comm
   integer :: ierr

   ! Timing
   type(timer_type) :: timer
   real(dp) :: elapsed_time, flops

   ! Fragment data structures
   integer(default_int), allocatable :: monomers(:)
   integer(default_int), allocatable :: polymers(:, :)
   integer(default_int) :: n_expected_fragments, fragment_count
   type(pic_fragment_block), allocatable :: fragments(:)

   ! Node leader detection
   integer :: global_node_rank
   integer, allocatable :: all_node_leader_ranks(:), node_leader_ranks(:)
   integer :: i, j, num_nodes

   !==============================
   ! MPI Initialization
   !==============================
   call MPI_Init(ierr)
   world_comm = comm_world()
   node_comm  = world_comm%split()   ! shared memory communicator

   if (world_comm%size() < 3) then
      if (world_comm%leader()) print *, "This program requires at least 3 processes."
      call MPI_Abort(world_comm%get(), 1, ierr)
   end if

   !==============================
   ! Determine node leaders
   !==============================
   global_node_rank = -1
   if (node_comm%leader()) global_node_rank = world_comm%rank()

   allocate(all_node_leader_ranks(world_comm%size()))
   call MPI_Allgather(global_node_rank, 1, MPI_INTEGER, all_node_leader_ranks, 1, MPI_INTEGER, world_comm%get(), ierr)

   num_nodes = count(all_node_leader_ranks /= -1)
   allocate(node_leader_ranks(num_nodes))
   i = 0
   do j = 1,world_comm%size()
      if (all_node_leader_ranks(j) /= -1) then
         i = i + 1
         node_leader_ranks(i) = all_node_leader_ranks(j)
      end if
   end do

   !==============================
   ! Generate fragments (on rank 0)
   !==============================
   if (world_comm%leader()) then
      n_expected_fragments = get_nfrags(n_monomers, max_level)
      allocate(monomers(n_monomers))
      allocate(polymers(n_expected_fragments, max_level))
      allocate(fragments(n_expected_fragments))

      print *, "Expected fragments = ", n_expected_fragments

      polymers = 0_default_int
      call create_monomer_list(monomers)
      polymers(:n_monomers, 1) = monomers(:n_monomers)
      fragment_count = n_monomers

      call generate_fragment_list(monomers, max_level, polymers, fragment_count)

      if (fragment_count /= n_expected_fragments) then
         print *, "Error: fragment_count does not match expected fragments!"
         call MPI_Abort(world_comm%get(), 1, ierr)
      end if

      print *, "Generated ", fragment_count, " fragments for distribution"
      call timer%start()
   end if

   !==============================
   ! Role assignment
   !==============================
   if (world_comm%leader() .and. node_comm%leader()) then
      call global_coordinator(world_comm, node_comm, fragment_count, polymers, max_level, node_leader_ranks, num_nodes)
   else if (node_comm%leader()) then
      call node_coordinator(world_comm, node_comm, max_level)
   else
      call node_worker(world_comm, node_comm, n, max_level)
   end if

   !==============================
   ! Final timing and flops
   !==============================
   call world_comm%barrier()

   if (world_comm%leader()) then
      call timer%stop()
      elapsed_time = timer%get_elapsed_time()
      flops = 0.0_dp
      call calculate_exact_flops(polymers, fragment_count, max_level, n, flops)

      print *, "Total elapsed time for all fragments:", elapsed_time, "seconds"
      print *, "Total flops ", flops
      print *, "Estimated flop rate: ", flops / elapsed_time / 1.0e9_dp, " GFLOP/s"
   end if

   !==============================
   ! Finalization
   !==============================
   call node_comm%finalize()
   call world_comm%finalize()
   call MPI_Finalize()

end program hierarchical_mpi_mbe
