module pic_fragment
   use pic_types, only: default_int, dp
   use pic_string_utils, only: to_string
   use pic_matrix_printer_v2, only: print_array_v2
   implicit none
   private
   public :: pic_fragment_type
   public :: new_fragment
   public :: count_nonzeros

   type :: pic_fragment_block
      integer(default_int), allocatable :: indices(:)
      real(dp), allocatable :: matrix(:, :)
   contains
      procedure :: print_indices
   end type pic_fragment_block

   type :: pic_fragment_type 
    integer(default_int), allocatable :: indices(:)
    integer(default_int) :: n_monomers
    integer(default_int) :: n_basis_functions
    contains 
      procedure :: get_estimated_cost
      procedure :: print_fragment
   end type pic_fragment_type


contains
  subroutine print_fragment(self)
    class(pic_fragment_type), intent(in) :: self 

    print *, "Indices", self%indices 
    print *, "n monomers", self%n_monomers 
  end subroutine print_fragment
  function get_estimated_cost(self) result(cost)
    class(pic_fragment_type), intent(in) :: self 
    integer(default_int) :: cost 

    cost = self%n_basis_functions * self%n_basis_functions 
  end function get_estimated_cost

  function new_fragment(indices, n_basis_functions) result(fragment)
    integer(default_int), intent(in) :: indices(:)
    integer(default_int), intent(in) :: n_basis_functions 
    type(pic_fragment_type) :: fragment 

    fragment%n_basis_functions = n_basis_functions 
    allocate(fragment%indices, source=indices)
    fragment%n_monomers = size(fragment%indices)
  end function new_fragment

   function count_nonzeros(row) result(count)
      integer(default_int), intent(in) :: row(:)
      integer(default_int) :: i
      integer(default_int) :: count
      count = 0
      do i = 1, size(row)
         if (row(i) /= 0) count = count + 1
      end do
   end function count_nonzeros

   subroutine print_indices(self)
      class(pic_fragment_block), intent(in) :: self

      call print_array_v2(self%indices)

   end subroutine print_indices


end module pic_fragment
