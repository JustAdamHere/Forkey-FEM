module class_element
    use common

    implicit none

    private
    public :: element

    type, abstract :: element
        integer                             :: elementNo
        integer                             :: noNodes
        integer                             :: polynomialDegree
        integer, allocatable, dimension(:)  :: nodeIndices
        real(dp), dimension(:), allocatable :: nodeCoordinates
    contains
        ! Getters.
        procedure(interface_get_elementQuadrature), deferred :: get_elementQuadrature
        procedure(interface_get_Jacobian), deferred          :: get_Jacobian

        ! Misc methods.
        procedure(interface_basisLegendre), deferred    :: basisLegendre
        procedure(interface_basisLobatto), deferred     :: basisLobatto
        procedure(interface_mapLocalToGlobal), deferred :: mapLocalToGlobal
    end type element

    abstract interface
        subroutine interface_get_elementQuadrature(this, a_coordinates, a_weights)
            use common 
            import element

            class(element), intent(inout) :: this
            real(dp), dimension(10)       :: a_coordinates
            real(dp), dimension(10)       :: a_weights
        end subroutine
    end interface

    abstract interface
        function interface_get_Jacobian(this)
            use common
            import element

            class(element), intent(inout) :: this
            real(dp)                      :: interface_get_Jacobian
        end function
    end interface

    abstract interface
        function interface_basisLegendre(this, a_degree, a_deriv, a_point)
            use common
            import element

            class(element), intent(inout) :: this
            integer                       :: a_degree
            integer                       :: a_deriv
            real(dp)                      :: a_point
            real(dp)                      :: interface_basisLegendre
        end function
    end interface

    abstract interface
        recursive function interface_basisLobatto(this, a_degree, a_deriv, a_point) result(BL)
            use common
            import element

            class(element), intent(inout) :: this
            integer                       :: a_degree
            integer                       :: a_deriv
            real(dp)                      :: a_point
            real(dp)                      :: BL
        end function
    end interface

    abstract interface
        function interface_mapLocalToGlobal(this, a_point)
            use common
            import element

            class(element), intent(inout) :: this
            real(dp)                      :: a_point
            real(dp)                      :: interface_mapLocalToGlobal
        end function
    end interface

contains
    
end module class_element