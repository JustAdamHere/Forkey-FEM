module class_mesh
    use common
    !use point
    use class_element
    use class_element_interval

    implicit none

    public

    type mesh
        integer                                   :: noElements
        integer                                   :: noFaces
        class(element), dimension(:), allocatable :: elements
        integer                                   :: problemDimension
        integer, dimension(:, :), allocatable     :: elementConnectivity
        real(dp), dimension(:), allocatable       :: nodeCoordinates
    contains
        ! Constructors.
        procedure :: constructor => mesh_constructor
        procedure :: constructor_ex => mesh_constructor_ex
        procedure :: constructor_eq => mesh_constructor_eq

        ! Deconstructors.
        procedure :: destructor => mesh_deconstructor
    end type mesh

contains
    subroutine mesh_constructor(this, a_problemDimension)
        class(mesh) :: this
        integer     :: a_problemDimension
    end subroutine

    subroutine mesh_constructor_ex(this, a_exampleNo)
        class(mesh) :: this
        integer     :: a_exampleNo

        integer, dimension(2) :: nodeIndices
        integer               :: polynomialDegree

        real(dp) :: domainLeft
        real(dp) :: domainRight
        real(dp) :: h

        integer :: i, j

        !type(element_types), pointer :: temp

        polynomialDegree = 1

        ! Simple 1D interval with 4 equal-sized elements on [0, 1].
        if (a_exampleNo == 0) then
            this%problemDimension = 1

            this%noElements  = 4
            this%noFaces     = 5
            domainLeft  = 0
            domainRight = 1
            h = (domainRight - domainLeft)/this%noElements

            ! Node coordinates.
            allocate(this%nodeCoordinates(this%noElements+1))
            this%nodeCoordinates(1) = domainLeft

            ! Elements storage.
            allocate(this%elements(this%noElements))

            do i = 1, size(this%elements)
                allocate(element_interval :: this%elements(i)%element_type)

                this%nodeCoordinates(i+1) = domainLeft + i*h

                nodeIndices(1) = i
                nodeIndices(2) = i+1

                call this%elements(i)%element_type%constructor(i, &
                    nodeIndices, this%nodeCoordinates, polynomialDegree)
            end do

            ! Element connectivity.
            allocate(this%elementConnectivity(this%noElements, 2))

            do i = 1, this%noElements
                if (i-1 > 0) then
                    this%elementConnectivity(i, 1) = i-1
                end if

                if (i+1 <= this%noElements) then
                    this%elementConnectivity(i, 2) = i+1
                end if
            end do
        end if
    end subroutine

    subroutine mesh_constructor_eq(this, a_noElements, a_polynomialDegree)
        class(mesh) :: this
        integer     :: a_noElements
        integer     :: a_polynomialDegree

        integer, dimension(2) :: nodeIndices

        real(dp) :: domainLeft
        real(dp) :: domainRight
        real(dp) :: h

        integer :: i, j

        !type(element_types), pointer :: temp

        this%problemDimension = 1

        this%noElements  = a_noElements
        this%noFaces     = a_noElements+1
        domainLeft  = 0
        domainRight = 1
        h = (domainRight - domainLeft)/this%noElements

        ! Node coordinates.
        allocate(this%nodeCoordinates(this%noElements+1))
        this%nodeCoordinates(1) = domainLeft

        ! Elements storage.
        allocate(this%elements(this%noElements))

        do i = 1, size(this%elements)
            allocate(element_interval :: this%elements(i)%element_type)

            this%nodeCoordinates(i+1) = domainLeft + i*h

            nodeIndices(1) = i
            nodeIndices(2) = i+1

            call this%elements(i)%element_type%constructor(i, &
                nodeIndices, this%nodeCoordinates, a_polynomialDegree)
        end do

        ! Element connectivity.
        allocate(this%elementConnectivity(this%noElements, 2))

        do i = 1, this%noElements
            if (i-1 > 0) then
                this%elementConnectivity(i, 1) = i-1
            end if

            if (i+1 <= this%noElements) then
                this%elementConnectivity(i, 2) = i+1
            end if
        end do
    end subroutine

    subroutine mesh_deconstructor(this)
        class(mesh) :: this

        integer :: i

        deallocate(this%elementConnectivity)

        do i = 1, size(this%elements)
            call this%elements(i)%element_type%destructor()
            deallocate(this%elements(i)%element_type)
        end do

        deallocate(this%elements)

        deallocate(this%nodeCoordinates)
    end subroutine

end module class_mesh