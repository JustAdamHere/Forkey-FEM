module class_matrix_sparse
!     use common
!     use class_matrix

!     implicit none

!     type, extends(matrix) :: matrix_sparse
!         real(dp), dimension(:), allocatable :: entries
!         integer, dimension(:), allocatable  :: columnNos
!         integer, dimension(:), allocatable  :: rowStarts
!     contains
!         procedure :: constructor
!         procedure :: destructor

!         procedure :: get_entry

!         procedure :: multiply
!     end type

! contains
!     subroutine constructor(this, a_entries, a_columnNos, a_rowStarts)
!         class(matrix_sparse)                :: this
!         real(dp), dimension(:), allocatable :: a_entries
!         integer, dimension(:), allocatable  :: a_columnNos
!         integer, dimension(:), allocatable  :: a_rowStarts

!         this%entries   = a_entries
!         this%columnNos = a_columnNos
!         this%rowStarts = a_rowStarts
!     end subroutine

!     subroutine destructor(this)
!         class(matrix_sparse) :: this
!     end subroutine

!     function get_entry(this, a_i, a_j)
!         class(matrix_sparse) :: this
!         integer              :: a_i
!         integer              :: a_j
!         real(dp)             :: get_entry

!         integer :: i, rowStart, rowEnd, rowCurrent

!         rowStart = this%rowStarts(a_j)
!         rowEnd   = this%rowStarts(a_j+1)

!         i = rowStart
!         do while(i<a_i)
!             i = i+1
!         end do

!         rowCurrent = this%rowStarts(i)
!         if (rowCurrent == a_i) then
!             get_entry = this%entries(i)
!         else
!             get_entry = 0
!         end if
!     end function

!     function multiply(this, a_matrix)
!         class(matrix_sparse) :: this
!         class(matrix)        :: a_matrix

!         class(matrix_sparse), pointer :: multiply
!     end function

end module class_matrix_sparse