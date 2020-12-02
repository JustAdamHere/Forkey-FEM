module class_mesh
    use common
    use class_element
    use class_element_interval

    implicit none

    type mesh
        private
            integer :: noElements
            class(element), dimension(:), allocatable :: elements
    contains
        ! Deconstructors.
        procedure :: mesh_deconstructor
    end type mesh

    interface mesh
        procedure mesh_constructor
    end interface

contains
    function mesh_constructor(a_exampleNo)
        type(mesh) :: mesh_constructor
        integer    :: a_exampleNo

        integer                             :: elementNo
        integer                             :: noNodes
        integer, dimension(2)               :: nodeIndices
        real(dp), dimension(:), allocatable :: nodeCoordinates
        integer                             :: polynomialDegree

        integer        :: i

        elementNo = 1
        noNodes = 2
        nodeIndices(1) = 1
        nodeIndices(2) = 2
        polynomialDegree = 1

        allocate(nodeCoordinates(2))
        nodeCoordinates(1) = 0
        nodeCoordinates(2) = 1

        ! Simple 1D interval with 4 equal-sized elements on [0, 1].
        if (a_exampleNo == 0) then
            !allocate(element_interval(elementNo, noNodes, nodeIndices, &
            !    nodeCoordinates, polynomialDegree) :: mesh_constructor%elements(4))
            allocate(element_interval :: mesh_constructor%elements(4))

            do i = 1, size(mesh_constructor%elements)
                mesh_constructor%elements(i) ! Can't seem to be able to index this.
            end do
        end if

        deallocate(nodeCoordinates)
    end function

    subroutine mesh_deconstructor(this)
        class(mesh) :: this

        deallocate(this%elements)
    end subroutine

end module class_mesh