module class_mesh
    use common
    !use point
    use class_element
    use class_element_interval

    implicit none

    public

    type mesh
        integer :: noElements
        class(element), dimension(:), allocatable :: elements
    contains
        ! Constructors.
        procedure :: constructor => mesh_constructor
        ! Deconstructors.
        procedure :: destructor => mesh_deconstructor
    end type mesh

contains
    subroutine mesh_constructor(this, a_exampleNo)
        class(mesh) :: this
        integer    :: a_exampleNo

        integer                             :: noNodes
        integer, dimension(2)               :: nodeIndices
        real(dp), dimension(:), allocatable :: nodeCoordinates
        integer                             :: polynomialDegree

        real(dp) :: domainLeft
        real(dp) :: domainRight
        real(dp) :: h
        integer  :: noElements

        integer        :: i

        !type(element_types), pointer :: temp

        noNodes = 2
        nodeIndices(1) = 1
        nodeIndices(2) = 2
        polynomialDegree = 1

        allocate(nodeCoordinates(2))
        nodeCoordinates(1) = 0
        nodeCoordinates(2) = 1

        ! Simple 1D interval with 4 equal-sized elements on [0, 1].
        if (a_exampleNo == 0) then
            noElements  = 4
            domainLeft  = 0
            domainRight = 1
            h = (domainRight - domainLeft)/noElements

            allocate(this%elements(noElements))

            do i = 1, size(this%elements)
                allocate(element_interval :: this%elements(i)%element_type)

                ! nodeIndices(1) = domainLeft + (i-1)*h
                ! nodeIndices(2) = domainLeft + i*h

                nodeIndices(1) = i
                nodeIndices(2) = i+1
                

                call this%elements(i)%element_type%constructor(i, &
                    noNodes, nodeIndices, nodeCoordinates, polynomialDegree)
            end do
        end if

        deallocate(nodeCoordinates)
    end subroutine

    subroutine mesh_deconstructor(this)
        class(mesh) :: this

        integer :: i

        do i = 1, size(this%elements)
            call this%elements(i)%element_type%destructor()
            deallocate(this%elements(i)%element_type)
        end do

        deallocate(this%elements)
    end subroutine

end module class_mesh