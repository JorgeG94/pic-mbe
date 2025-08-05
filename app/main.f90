program main
   use sample, only: hello_there
   use second_sample, only: answer
   use pic, only: pic_print_banner
   use pic_matrix_printer_v2, only: print_array_v2
   use pic_mbe, only: get_nfrags, create_monomer_list, generate_fragment_list
   use pic_timer, only: pic_timer_type
   use pic_fragment, only: pic_fragment_block, count_nonzeros
   use pic_types, only: dp, default_int
   use pic_string_utils, only: to_string
   implicit none

   integer(default_int), parameter :: n_monomers = 15
   type(pic_timer_type) :: my_timer
   integer(default_int), parameter :: max_level = 3
   integer(default_int), allocatable :: monomers(:)
   integer(default_int), allocatable :: polymers(:, :)
   integer(default_int) :: n_expected_fragments
   integer(default_int) :: i
   integer(default_int) :: fragment_count, n_bytes
   integer(default_int), parameter :: n = 256
   type(pic_fragment_block), allocatable :: fragments(:)
   real(dp) :: elapsed_time

   n_expected_fragments = get_nfrags(n_monomers, max_level)

   allocate (monomers(n_monomers))

   allocate (polymers(n_expected_fragments, max_level))
   print *, "frags = ", n_expected_fragments

   ! fragment block type allocated to the total number of frags
   allocate (fragments(n_expected_fragments))

   ! initialize the array to 0 because life is good eh
   polymers = 0_default_int
   call create_monomer_list(monomers)

   polymers(:n_monomers, 1) = monomers(:n_monomers)

   fragment_count = n_monomers

   call my_timer%start()
   call generate_fragment_list(monomers, max_level, polymers, fragment_count)
   call my_timer%stop()
   elapsed_time = my_timer%get_elapsed_time()
   print *, "fragment_count of polymer combinations is "//to_string(fragment_count)

   if(fragment_count /= n_expected_fragments) then
      print *, "Error: fragment_count of polymer combinations does not match expected fragments!"
      stop 1
   end if

   print *, "Time for polymer combinations is "//to_string(elapsed_time)


   call my_timer%start()
   make_frags: block
      integer(default_int) :: k

      do i = 1, n_monomers
         allocate (fragments(i)%matrix(n, n))
         fragments(i)%matrix = real(i, dp)
      end do

      !$omp parallel do private(i,k)
      do i = 1, n_expected_fragments

         k = count_nonzeros(polymers(i, :))

         ! only for polymers
         if (k > 1) then
            allocate (fragments(i)%indices(k))
            fragments(i)%indices = polymers(i, 1:k)

            allocate (fragments(i)%matrix(k*n, k*n))
            fragments(i)%matrix = real(k*i, dp)
         end if

      end do
      !$omp end parallel do

   end block make_frags
   call my_timer%stop()
   elapsed_time = my_timer%get_elapsed_time()
   print *, "Time to generate polymers "//to_string(elapsed_time)


end program main
