module pic_fragment
   use pic_types, only: default_int, dp
   use pic_array, only: pic_print_array
   use pic_string, only: to_string
   implicit none
   private
   public :: pic_fragment_type
   public :: new_fragment
   public :: count_nonzeros
   public :: to_string

   type :: pic_fragment_block
      integer(default_int), allocatable :: indices(:)
      real(dp), allocatable :: matrix(:, :)
   contains
      procedure :: print_indices
   end type pic_fragment_block

   type :: pic_fragment_type 
    integer(default_int), allocatable :: indices(:)
    integer(default_int) :: fragment_id
    integer(default_int) :: n_monomers
    integer(default_int) :: n_basis_functions
    contains 
      procedure :: get_estimated_cost
   end type pic_fragment_type


   enum, bind(c) 
    enumerator :: MONOMER = 1
    enumerator :: DIMER = 2
    enumerator :: TRIMER = 3 
    enumerator :: TETRAMER = 4
    enumerator :: PENTAMER = 5
   end enum
    
    interface to_string
      procedure :: fragment_to_string
    end interface

contains

function fragment_to_string(fragment) result(trimmed_str)
    class(pic_fragment_type), intent(in) :: fragment
    character(len=:), allocatable :: trimmed_str
    character(len=32) :: temp_str
    character(len=500) :: full_str
    integer :: i

    ! Get fragment type name
    character(len=:), allocatable :: fragment_type_name

    select case (fragment%n_monomers)
        case (MONOMER)
            fragment_type_name = "MONOMER"
        case (DIMER)
            fragment_type_name = "DIMER"
        case (TRIMER)
            fragment_type_name = "TRIMER"
        case (TETRAMER)
            fragment_type_name = "TETRAMER"
        case (PENTAMER)
            fragment_type_name = "PENTAMER"
        case default
            write(temp_str, '(I0,A)') fragment%n_monomers, "-MER"
            fragment_type_name = trim(temp_str)
    end select

    ! Build the complete string
    full_str = "Fragment info for " // fragment_type_name // new_line('a')

    ! Add fragment ID
    write(temp_str, '(I0)') fragment%fragment_id
    full_str = trim(full_str) // " Fragment ID: " // trim(temp_str) // new_line('a')

    write(temp_str, '(I0)') fragment%n_basis_functions 
    full_str = trim(full_str) // " Basis functions: " // trim(temp_str) // new_line('a')
    ! Add indices array
    full_str = trim(full_str) // " Indices: ["
    do i = 1, size(fragment%indices)
        write(temp_str, '(I0)') fragment%indices(i)
        if (i < size(fragment%indices)) then
            full_str = trim(full_str) // trim(temp_str) // ", "
        else
            full_str = trim(full_str) // trim(temp_str)
        end if
    end do
    full_str = trim(full_str) // "]" // new_line('a')

    ! Allocate and return the trimmed result
    trimmed_str = trim(full_str)
end function fragment_to_string 

  function get_estimated_cost(self) result(cost)
    class(pic_fragment_type), intent(in) :: self 
    integer(default_int) :: cost 

    cost = self%n_basis_functions * self%n_basis_functions 
  end function get_estimated_cost

  function new_fragment(indices, id, n_basis_functions) result(fragment)
    integer(default_int), intent(in) :: indices(:)
    integer(default_int), intent(in) :: id
    integer(default_int), intent(in) :: n_basis_functions 
    type(pic_fragment_type) :: fragment 

    fragment%n_basis_functions = n_basis_functions 
    fragment%fragment_id = id
    allocate(fragment%indices, source=indices)
    fragment%n_monomers = count(fragment%indices > 0)
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

      call pic_print_array(self%indices)

   end subroutine print_indices


end module pic_fragment
