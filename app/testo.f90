program hierarchical_mpi_mbe
   use mpi_f08
   use mpi_comm_simple        
   use pic_blas_interfaces, only: pic_gemm
   use pic_timer, only: pic_timer_type
   use pic_types, only: dp, default_int
   use pic_string, only: to_string
   use pic_mbe, only: get_nfrags, create_monomer_list, generate_fragment_list, print_combos
   use pic_fragment, only: pic_fragment_type,new_fragment, count_nonzeros, to_string
   use pic_mpi_algorithms
   implicit none

   ! Fragment generation parameters
   integer(default_int), parameter :: n_monomers = 10
   integer(default_int), parameter :: max_level = 3
   integer(default_int), parameter :: n = 256  ! monomer matrix size

   ! MPI wrappers
   type(comm_t) :: world_comm, node_comm
   integer :: ierr

   ! Timing
   type(pic_timer_type) :: timer
   real(dp) :: elapsed_time, flops

   ! Fragment data structures
   integer(default_int), allocatable :: monomers(:)
   integer(default_int), allocatable :: polymers(:,:)
   integer(default_int) :: n_expected_fragments, fragment_count
   type(pic_fragment_type), allocatable :: fragments(:)

   ! Node leader detection
   integer :: global_node_rank
   integer, allocatable :: all_node_leader_ranks(:), node_leader_ranks(:)
   integer :: i, j, num_nodes


      n_expected_fragments = get_nfrags(n_monomers, max_level)
      allocate(monomers(n_monomers))
      allocate(polymers(n_expected_fragments, max_level))
      allocate(fragments(n_expected_fragments))

      print *, "Expected fragments = ", n_expected_fragments

      polymers = 0_default_int
      call create_monomer_list(monomers)
      polymers(:n_monomers,1) = monomers(:n_monomers)
      fragment_count = n_monomers

      call generate_fragment_list(monomers, max_level, polymers, fragment_count)

      print *, polymers(1,:)


      if (fragment_count /= n_expected_fragments) then
         print *, "Error: fragment_count does not match expected fragments!"
      end if

      print *, "Generated ", fragment_count, " fragments for distribution"

      block 
      integer(default_int), allocatable :: frag_indices(:)
      integer(default_int) :: frag_level, matrix_size

      matrix_size = n 

      do i = 1, fragment_count 
        fragments(i) = new_fragment(polymers(i,:),i, n)
        print *, to_string(fragments(i))
      end do 

      end block 

      call timer%start()
!      do i = fragment_count, 1, -1 
!        call process_fragment(i,fragments(i)%indices, fragments(i)%n_monomers, fragments(i)%n_basis_functions)
!      end do 

!      call print_combos(polymers)

      !do i = fragment_count, 1, -1
      !  frag_indices = polymers(i,:)
      !  frag_level = size(frag_indices)
      !  call process_fragment(i, frag_indices, frag_level, matrix_size)
      !end do 



      call timer%stop()

      elapsed_time = timer%get_elapsed_time()
      flops = 0.0_dp
      call calculate_exact_flops(polymers, fragment_count, max_level, n, flops)

      print *, "Total elapsed time for all fragments:", elapsed_time, "seconds"
      print *, "Total flops ", flops
      print *, "Estimated flop rate: ", flops / elapsed_time / 1.0e9_dp, " GFLOP/s"


end program hierarchical_mpi_mbe
