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
        ! Constructors.
        procedure :: constructor => mesh_constructor
        ! Deconstructors.
        procedure :: destructor => mesh_deconstructor
    end type mesh

contains
    subroutine mesh_constructor(this, a_exampleNo)
        class(mesh) :: this
        integer    :: a_exampleNo

        integer                             :: elementNo
        integer                             :: noNodes
        integer, dimension(2)               :: nodeIndices
        real(dp), dimension(:), allocatable :: nodeCoordinates
        integer                             :: polynomialDegree

        integer        :: i

        !type(element_types), pointer :: temp

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
            !    nodeCoordinates, polynomialDegree) :: this%elements(4))
            allocate(this%elements(4))

            do i = 1, size(this%elements)
                allocate(element_interval :: this%elements(i)%element_type)

                !temp = this%elements(i)
                !temp%element_type%element%constructor(elementNo, &
                !    noNodes, nodeIndices, nodeCoordinates, polynomialDegree)
                !call temp%constructor(elementNo, &
                !    noNodes, nodeIndices, nodeCoordinates, polynomialDegree)

                call this%elements(i)%element_type%constructor(elementNo, &
                    noNodes, nodeIndices, nodeCoordinates, polynomialDegree)

                print *, this%elements(i)%element_type%get_Jacobian()
            end do
        end if
    end subroutine

    subroutine mesh_deconstructor(this)
        class(mesh) :: this

        integer :: i

        do i = 1, size(this%elements)
            call this%elements(i)%element_type%destructor()
            !deallocate(this%elements(i)%element_type)
        end do

        !deallocate(this%elements)
    end subroutine

end module class_mesh