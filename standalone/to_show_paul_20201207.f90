!-------------------------------------------------------------------
!> MODULE
!! class_element
!!
!! Abstract class for element, holding only the element number for 
!!   global element numbering for demonstration purposes.
!-------------------------------------------------------------------
module class_element
    implicit none

    public 

    type, abstract :: element
        integer :: elementNo
    end type
end module

!-------------------------------------------------------------------
!> MODULE
!! class_element_interval
!! 
!! Class for an element that is an interval. Holds only the left and
!!   right coordinate of the element for demonstration purposes.
!-------------------------------------------------------------------
module class_element_interval  
    use class_element

    implicit none

    public

    type, extends(element) :: element_interval
        real :: leftSide
        real :: rightSide
    contains
        procedure element_interval_constructor
    end type

contains
    subroutine element_interval_constructor(this, a_elementNo, a_leftSide, a_rightSide)
        class(element_interval) :: this
        integer                 :: a_elementNo
        real                    :: a_leftSide
        real                    :: a_rightSide

        this%elementNo = a_elementNo
        this%leftSide  = a_leftSide
        this%rightSide = a_rightSide
    end subroutine
end module

!-------------------------------------------------------------------
!> MODULE
!! class_mesh
!!
!! Class for meshes, holding elements of any extended class from 
!!   'element' (although the constructor here just sets up element
!!   intervals).
!-------------------------------------------------------------------
module class_mesh
    use class_element
    use class_element_interval

    implicit none

    public

    type :: mesh
        integer                                   :: noElements
        class(element), dimension(:), allocatable :: meshElements
    contains
        procedure :: mesh_constructor
        procedure :: mesh_destructor
    end type

contains   
    subroutine mesh_constructor(this, a_noElements)
        class(mesh) :: this
        integer     :: a_noElements
        integer     :: i

        allocate(element_interval :: this%meshElements(a_noElements))

        do i = 1, size(this%meshElements)
            ! Just to be incredibly simple, we create several of the same element on [0, 1].
            this%meshElements(i)%element_interval_constructor(i, 0.0, 1.0)
        end do
    end subroutine

    subroutine mesh_destructor(this)
        class(mesh) :: this

        deallocate(this%meshElements)
    end subroutine
end module

!-------------------------------------------------------------------
!> PROGRAM
!! test
!!
!! A test program, just instantiating a mesh with 4 elements.
!-------------------------------------------------------------------
program test
    use class_element
    use class_element_interval
    use class_mesh

    implicit none

    type(mesh) :: myMesh

    call myMesh%mesh_constructor(4)

    print *, myMesh%noElements

    call myMesh%mesh_destructor()
end program