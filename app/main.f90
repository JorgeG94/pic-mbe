program main
   use sample, only: hello_there
   use second_sample, only: answer
   use pic, only: pic_print_banner
   use pic_matrix_printer_v2, only: print_array_v2 
   use pic_mbe
   use pic_timer
   use pic_fragment
   implicit none

   integer(default_int), parameter :: n_monomers = 15
   type(pic_timer_type) :: my_timer
   integer(default_int), parameter :: max_level = 3
   integer(default_int), allocatable :: monomers(:)
   integer(default_int), allocatable :: polymers(:,:)
   integer(default_int), allocatable :: fragment_count(:,:)
   integer(default_int) :: n_expected_fragments
   integer(default_int) :: i
   integer(default_int) :: count, n_bytes
   integer(default_int), parameter :: n = 256
   type(pic_fragment_block), allocatable :: fragments(:)
   real(dp) :: elapsed_time

   n_expected_fragments = n_monomers

   allocate(monomers(n_monomers))
   allocate(fragment_count(1,max_level))
   fragment_count(1,1) = n_monomers

   do i = 2, max_level 
    fragment_count(1,i) = binomial(n_monomers,i)
    n_expected_fragments = n_expected_fragments + binomial(n_monomers, i)
   end do 
   allocate(polymers(n_expected_fragments, max_level))
   print *, "frags = ", n_expected_fragments

   ! fragment block type allocated to the total number of frags
   allocate(fragments(n_expected_fragments))

   ! initialize the array to 0 because life is good eh
   polymers = 0_default_int
   call create_monomer_list(monomers)
   do i = 1, n_monomers
    polymers(i,1) = monomers(i)
   end do
   count = n_monomers

   n_bytes = 0
   do i = 1, max_level 
    n_bytes = n_bytes + i*n*n
   end do 
   print *, " bytes to be allocated: ", real(n_bytes * 8,dp)/1e9_dp

   call my_timer%start()
   call generate_combinations(monomers, max_level, polymers, count)
   call my_timer%stop()
   elapsed_time = my_timer%get_elapsed_time()
   print *, "Time for combinations is " // to_string(elapsed_time)

   call my_timer%start()
   make_frags: block 
    integer(default_int) :: k 

    do i = 1, n_monomers 
      allocate(fragments(i)%matrix(n,n))
      fragments(i)%matrix = real(i,dp)
    end do 

    !$omp parallel do private(i,k)
    do i = 1, n_expected_fragments 

      k = count_nonzeros(polymers(i,:))

      ! only for polymers
      if(k > 1) then
        allocate(fragments(i)%indices(k))
        fragments(i)%indices = polymers(i,1:k)

        allocate(fragments(i)%matrix(k*n,k*n))
        fragments(i)%matrix = real(k*i, dp) 
      end if

    end do
    !$omp end parallel do 

   end block make_frags
   call my_timer%stop()
   elapsed_time = my_timer%get_elapsed_time()
   print *, "Time to generate polymers " // to_string(elapsed_time)


end program main
