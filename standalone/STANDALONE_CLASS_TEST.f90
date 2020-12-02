module class_element
  implicit none

  type, public :: Element
     integer :: elementNo
     real, dimension(2) :: nodeCoordinates
   contains
     procedure :: Element_cctor
     procedure :: pnt => Element_print
  end type Element

contains
!------------------------------------------------------------------
!>   Prints out the data associated with this element.
!!
!! Author: Adam Blakey
!! 
!! Date Created :: 2020-11-12
!------------------------------------------------------------------
  subroutine Element_print(this)
!------------------------------------------------------------------
    class(Element), intent(in) :: this
    
    print *, 'Element number = ', this%elementNo, ' coord 1 = ', this%nodeCoordinates(1), ' coord 2 = ', this%nodeCoordinates(2)
  end subroutine Element_print

!------------------------------------------------------------------
!>   Constructor for this instance, taking in the the element
!!   number and the node coordinates.
!!
!! Author: Adam Blakey
!! 
!! Date Created :: 2020-11-12
!------------------------------------------------------------------
  subroutine Element_cctor(this, a_elementNo, a_nodeCoordinates)
!------------------------------------------------------------------
    class(Element), intent(out) :: this
    integer, intent(in) :: a_elementNo
    real, dimension(2), intent(in) :: a_nodeCoordinates

    this%elementNo = a_elementNo
    this%nodeCoordinates = a_nodeCoordinates
  end subroutine

end module class_element

program element_test
  use class_element
  implicit none

  real, dimension(2) :: coords
  type(Element) :: ele

  coords(1) = 5.5
  coords(2) = 2.1

  
  ele = Element(1, coords)
  call ele%pnt
end program element_test