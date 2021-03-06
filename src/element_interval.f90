module class_element_interval
    use common
    use quadrature
    use class_point
    use class_element

    implicit none

    public

    type, extends(element_types) :: element_interval

    contains
        ! Constructors.
        procedure :: constructor => element_interval_constructor

        ! Destructors.
        procedure :: destructor => element_interval_destructor

        ! Getters.
        procedure :: get_elementQuadrature => element_interval_get_elementQuadrature
        procedure :: get_Jacobian          => element_interval_get_Jacobian
        procedure :: get_localCoordinates  => element_interval_get_localCoordinates
        procedure :: get_globalCoordinates => element_interval_get_globalCoordinates
        procedure :: get_typeName          => element_interval_get_typeName

        ! Maps.
        procedure :: map_localToGlobal => element_interval_map_localToGlobal

        ! Misc methods.
        procedure :: basisLegendre         => element_interval_basisLegendre
        procedure :: basisLobatto          => element_interval_basisLobatto
        procedure :: evaluateBasisLegendre => element_interval_evaluateBasisLegendre
        procedure :: evaluateBasisLobatto  => element_interval_evaluateBasisLobatto
    end type element_interval

contains
    subroutine element_interval_constructor(this, a_elementNo, a_nodeIndices, a_nodeCoordinates, a_polynomialDegree)
        class(element_interval)                         :: this
        integer, intent(in)                             :: a_elementNo
        integer, dimension(2), intent(in)               :: a_nodeIndices
        real(dp), dimension(:), allocatable, intent(in) :: a_nodeCoordinates
        integer, intent(in)                             :: a_polynomialDegree

        this%elementNo        = a_elementNo
        this%nodeIndices      = a_nodeIndices
        this%nodeCoordinates  = a_nodeCoordinates
        this%polynomialDegree = a_polynomialDegree

        this%noNodes = 2
    end subroutine

    subroutine element_interval_destructor(this)
        class(element_interval) :: this

        deallocate(this%nodeCoordinates)
    end subroutine

    subroutine element_interval_get_elementQuadrature(this, a_coordinates, a_weights)
        class(element_interval), intent(inout) :: this
        real(dp), dimension(:), allocatable    :: a_coordinates
        real(dp), dimension(:), allocatable    :: a_weights

        integer :: n
        integer :: i

        n = ceiling(real(2*this%polynomialDegree + 1, dp)/2) + 1

        allocate(a_coordinates(n))
        allocate(a_weights    (n))

        do i = 1, n
            a_coordinates(i) = GaussLegendrePoint (n, i)
            a_weights    (i) = GaussLegendreWeight(n, i)
        end do
    end subroutine element_interval_get_elementQuadrature

    function element_interval_get_Jacobian(this)
        class(element_interval), intent(inout) :: this
        real(dp)                               :: element_interval_get_Jacobian

        element_interval_get_Jacobian = (this%get_globalCoordinates(2) - this%get_globalCoordinates(1))/2
    end function element_interval_get_Jacobian

    function element_interval_get_localCoordinates(this, a_coordinateNumber)
        class(element_interval), intent(inout) :: this
        real(dp)                               :: element_interval_get_localCoordinates
        integer                                :: a_coordinateNumber

        if (a_coordinateNumber == 1) then
            element_interval_get_localCoordinates = -1
        else
            element_interval_get_localCoordinates = 1
        end if
    end function

    function element_interval_get_globalCoordinates(this, a_coordinateNumber)
        class(element_interval), intent(inout) :: this
        real(dp)                               :: element_interval_get_globalCoordinates
        integer                                :: a_coordinateNumber

        integer :: index

        index = this%nodeIndices(a_coordinateNumber)

        element_interval_get_globalCoordinates = this%nodeCoordinates(index)
    end function

    function element_interval_get_typeName(this)
        class(element_interval), intent(inout) :: this
        character(len=32)                      :: element_interval_get_typeName

        element_interval_get_typeName = 'element_interval'
    end function

    function element_interval_map_localToGlobal(this, a_xi)
        class(element_interval), intent(inout) :: this
        real(dp)                               :: a_xi
        real(dp)                               :: element_interval_map_localToGlobal

        element_interval_map_localToGlobal = this%get_globalCoordinates(1) + (a_xi + 1)*this%get_Jacobian()
    end function

    function element_interval_basisLegendre(this, a_degree, a_deriv, a_point)
        use quadrature

        class(element_interval), intent(inout) :: this
        integer                                :: a_degree
        integer                                :: a_deriv
        real(dp)                               :: a_point
        real(dp)                               :: element_interval_basisLegendre
        real(dp)                               :: temp1
        real(dp)                               :: temp2

        element_interval_basisLegendre = legendrePolynomial(a_degree, a_deriv, a_point)
    end function

    recursive function element_interval_basisLobatto(this, a_degree, a_deriv, a_point) result(BL)
        class(element_interval), intent(inout) :: this
        integer                                :: a_degree
        integer                                :: a_deriv
        real(dp)                               :: a_point
        real(dp)                               :: BL
        real(dp)                               :: temp1
        real(dp)                               :: temp2

        if (-1 <= a_point .and. a_point <= 1) then
            if (a_deriv == 0 .and. a_degree <= 1) then
                if (a_degree == 0) then
                    BL = (1-a_point)/2
                else if (a_degree == 1) then
                    BL = (1+a_point)/2
                end if
            else if (a_deriv == 1 .and. a_degree <= 1) then
                if (a_degree == 0) then
                    BL = -0.5
                else if (a_degree == 1) then
                    BL = 0.5
                end if
            else
                temp1 = this%basisLegendre(a_degree,   a_deriv, a_point)
                temp2 = this%basisLegendre(a_degree-2, a_deriv, a_point)
                BL = sqrt(a_degree-1-0.5)*(temp1 - temp2)
            end if
        else
            BL = 0
        end if
    end function

    subroutine element_interval_evaluateBasisLegendre(this, a_basisPoints, a_degree, a_deriv, a_points)
        class(element_interval)             :: this
        real(dp), dimension(:), allocatable :: a_basisPoints
        integer                             :: a_degree
        integer                             :: a_deriv
        real(dp), dimension(:), allocatable :: a_points

        integer :: i, n

        n = size(a_points, 1)

        allocate(a_basisPoints(n))

        do i = 1, n
            a_basisPoints(i) = element_interval_basisLegendre(this, a_degree, a_deriv, a_points(i))
        end do
    end subroutine

    subroutine element_interval_evaluateBasisLobatto(this, a_basisPoints, a_degree, a_deriv, a_points)
        class(element_interval)             :: this
        real(dp), dimension(:), allocatable :: a_basisPoints
        integer                             :: a_degree
        integer                             :: a_deriv
        real(dp), dimension(:), allocatable :: a_points

        integer :: i, n

        n = size(a_points, 1)

        allocate(a_basisPoints(n))

        do i = 1, n
            a_basisPoints(i) = element_interval_basisLobatto(this, a_degree, a_deriv, a_points(i))
        end do
    end subroutine
end module class_element_interval