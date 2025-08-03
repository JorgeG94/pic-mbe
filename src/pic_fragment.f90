module pic_fragment
use pic_types, only: default_int, dp 
use pic_string_utils, only: to_string
use pic_matrix_printer_v2, only: print_array_v2 
implicit none 

type :: pic_fragment_block 
  integer(default_int), allocatable :: indices(:)
  real(dp), allocatable :: matrix(:,:)

  contains 
  
  procedure :: print_indices

end type pic_fragment_block

contains 

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

!print *, "Number of monomers " // to_string(size(self%indices,1))
call print_array_v2(self%indices)


end subroutine print_indices 


end module pic_fragment
