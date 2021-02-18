module class_matrix
!     use common

!     implicit none

!     type, abstract :: matrix
    
!     contains
!         procedure(interface_get_entry), deferred :: get_entry
!     end type

!     abstract interface
!         function interface_get_entry(this, a_i, a_j)
!             use common
!             import matrix

!             class(matrix) :: this
!             integer       :: a_i
!             integer       :: a_j
!             real(dp)      :: interface_get_entry
!         end function
!     end interface

! contains
    
end module class_matrix