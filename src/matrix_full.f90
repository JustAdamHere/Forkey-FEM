module class_matrix_full
!     use common
!     use class_matrix

!     implicit none

!     public

!     type, extends(matrix) :: matrix_full
!         real(dp), dimension(:, :), allocatable :: entries
!     contains
!         procedure :: constructor
!         procedure :: destructor

!         procedure :: get_entry

!         procedure :: multiply
!     end type

! contains
!     subroutine constructor(this, a_matrixEntries)
!         class(matrix_full)                     :: this
!         real(dp), dimension(:, :), allocatable :: a_matrixEntries

!         this%entries = a_matrixEntries
!     end subroutine

!     subroutine destructor(this)
!         class(matrix_full) :: this
!     end subroutine

!     function get_entry(this, a_i, a_j)
!         class(matrix_full) :: this
!         integer            :: a_i
!         integer            :: a_j
!         real(dp)           :: get_entry

!         get_entry = this%entries(a_i, a_j)
!     end function

!     function multiply(this, a_matrix)
!         class(matrix_full) :: this
!         class(matrix)      :: a_matrix

!         class(matrix_full), pointer :: multiply
!     end function

end module class_matrix_full