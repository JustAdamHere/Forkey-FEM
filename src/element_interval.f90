module class_element_interval
    use common
    use class_element
    use quadrature

    implicit none

    private
    public :: element_interval

    type, extends(element_types) :: element_interval

    contains
        ! Constructors.
        procedure :: constructor => element_interval_constructor

        ! Destructors.
        procedure :: destructor => element_interval_destructor

        ! Getters.
        procedure :: get_elementQuadrature
        procedure :: get_Jacobian

        ! Misc methods.
        procedure :: basisLegendre
        procedure :: basisLobatto
        procedure :: mapLocalToGlobal
    end type element_interval

contains
    subroutine element_interval_constructor(this, a_elementNo, a_noNodes, a_nodeIndices, a_nodeCoordinates, a_polynomialDegree)
        class(element_interval)                         :: this
        integer, intent(in)                             :: a_elementNo
        integer, intent(in)                             :: a_noNodes
        integer, dimension(2), intent(in)               :: a_nodeIndices
        real(dp), dimension(:), allocatable, intent(in) :: a_nodeCoordinates
        integer, intent(in)                             :: a_polynomialDegree

        this%elementNo        = a_elementNo
        this%noNodes          = a_noNodes
        this%nodeIndices      = a_nodeIndices
        this%nodeCoordinates  = a_nodeCoordinates
        this%polynomialDegree = a_polynomialDegree
    end subroutine

    subroutine element_interval_destructor(this)
        class(element_interval) :: this

        deallocate(this%nodeCoordinates)
    end subroutine

    subroutine get_elementQuadrature(this, a_coordinates, a_weights)
        class(element_interval), intent(inout) :: this
        real(dp), dimension(10)                :: a_coordinates
        real(dp), dimension(10)                :: a_weights
        
    end subroutine get_elementQuadrature

    function get_Jacobian(this)
        class(element_interval), intent(inout) :: this
        real(dp)                               :: get_Jacobian

        integer :: index1
        integer :: index2
        real(dp), dimension(2) :: nodeCoords 

        index1 = 1
        index2 = 2
        nodeCoords = this%nodeCoordinates

        get_Jacobian = nodeCoords(index2) - nodeCoords(index1)
    end function get_Jacobian

    function basisLegendre(this, a_degree, a_deriv, a_point)
        use quadrature

        class(element_interval), intent(inout) :: this
        integer                                :: a_degree
        integer                                :: a_deriv
        real(dp)                               :: a_point
        real(dp)                               :: basisLegendre
        real(dp)                               :: temp1
        real(dp)                               :: temp2

        basisLegendre = legendrePolynomial(a_degree, a_deriv, a_point)
    end function

    recursive function basisLobatto(this, a_degree, a_deriv, a_point) result(BL)
        class(element_interval), intent(inout) :: this
        integer                                :: a_degree
        integer                                :: a_deriv
        real(dp)                               :: a_point
        real(dp)                               :: BL
        real(dp)                               :: temp1
        real(dp)                               :: temp2

        if (-1 <= a_point .and. a_point <= 1) then
            if (a_degree == 0) then
                if (a_deriv == 0) then
                    BL = (1-a_point)/2
                else if (a_deriv == 1) then
                    BL = (1+a_point)/2
                else
                    BL = 0
                end if
            else if (a_degree == 1) then
                if (a_deriv == 0) then
                    BL = -0.5
                else if (a_deriv == 1) then
                    BL = 0.5
                else
                    BL = 0
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

    function mapLocalToGlobal(this, a_point)
        class(element_interval), intent(inout) :: this
        real(dp)                               :: a_point
        real(dp)                               :: mapLocalToGlobal

        integer :: index
        real(dp), dimension(2) :: nodeCoords 

        index = 1
        nodeCoords = this%nodeCoordinates

        mapLocalToGlobal = nodeCoords(1) + (a_point + 1)*this%get_Jacobian()
    end function
end module class_element_interval