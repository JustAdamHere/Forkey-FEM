module class_element
    use common

    implicit none

    private
    public :: element, element_types

    type element
        class(element_types), pointer :: element_type => null()
    end type

    type, abstract :: element_types
        integer                             :: elementNo
        integer                             :: noNodes
        integer                             :: polynomialDegree
        integer, allocatable, dimension(:)  :: nodeIndices
        real(dp), dimension(:), allocatable :: nodeCoordinates
    contains
        ! Constructors.
        procedure(interface_constructor), deferred :: constructor

        ! Destructors.
        procedure(interface_destructor), deferred :: destructor

        ! Getters.
        procedure(interface_get_elementQuadrature), deferred :: get_elementQuadrature
        procedure(interface_get_Jacobian), deferred          :: get_Jacobian
        procedure(interface_get_localCoordinates), deferred  :: get_localCoordinates
        procedure(interface_get_globalCoordinates), deferred :: get_globalCoordinates

        ! Maps.
        procedure(interface_map_localToGlobal), deferred :: map_localToGlobal

        ! Misc methods.
        procedure(interface_basisLegendre), deferred    :: basisLegendre
        procedure(interface_basisLobatto), deferred     :: basisLobatto
    end type

    abstract interface
        subroutine interface_constructor(this, a_elementNo, a_nodeIndices, a_nodeCoordinates, a_polynomialDegree)
            use common 
            import element_types

            class(element_types)                            :: this
            integer, intent(in)                             :: a_elementNo
            integer, dimension(2), intent(in)               :: a_nodeIndices
            real(dp), dimension(:), allocatable, intent(in) :: a_nodeCoordinates
            integer, intent(in)                             :: a_polynomialDegree
        end subroutine
    end interface

    abstract interface
        subroutine interface_destructor(this)
            use common 
            import element_types

            class(element_types) :: this
        end subroutine
    end interface

    abstract interface
        subroutine interface_get_elementQuadrature(this, a_coordinates, a_weights)
            use common 
            import element_types

            class(element_types), intent(inout) :: this
            real(dp), dimension(:), allocatable :: a_coordinates
            real(dp), dimension(:), allocatable :: a_weights
        end subroutine
    end interface

    abstract interface
        function interface_get_Jacobian(this)
            use common
            import element_types

            class(element_types), intent(inout) :: this
            real(dp)                            :: interface_get_Jacobian
        end function
    end interface

    abstract interface
        function interface_get_localCoordinates(this, a_coordinateNumber)
            use common
            import element_types

            class(element_types), intent(inout) :: this
            real(dp)                            :: interface_get_localCoordinates
            integer                             :: a_coordinateNumber
        end function
    end interface

    abstract interface
        function interface_get_globalCoordinates(this, a_coordinateNumber)
            use common
            import element_types

            class(element_types), intent(inout) :: this
            real(dp)                            :: interface_get_globalCoordinates
            integer                             :: a_coordinateNumber
        end function
    end interface

    abstract interface
        function interface_map_localToGlobal(this, a_xi)
            use common
            import element_types

            class(element_types), intent(inout) :: this
            real(dp)                            :: a_xi
            real(dp)                            :: interface_map_localToGlobal
        end function
    end interface

    abstract interface
        function interface_basisLegendre(this, a_degree, a_deriv, a_point)
            use common
            import element_types

            class(element_types), intent(inout) :: this
            integer                       :: a_degree
            integer                       :: a_deriv
            real(dp)                      :: a_point
            real(dp)                      :: interface_basisLegendre
        end function
    end interface

    abstract interface
        recursive function interface_basisLobatto(this, a_degree, a_deriv, a_point) result(BL)
            use common
            import element_types

            class(element_types), intent(inout) :: this
            integer                       :: a_degree
            integer                       :: a_deriv
            real(dp)                      :: a_point
            real(dp)                      :: BL
        end function
    end interface

contains
    
end module class_element