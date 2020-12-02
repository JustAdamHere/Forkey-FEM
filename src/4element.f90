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
        procedure                                            :: get_elementNo
        procedure(interface_get_elementQuadrature), deferred :: get_elementQuadrature
        procedure(interface_get_Jacobian), deferred          :: get_Jacobian
        procedure                                            :: get_nodeCoordinates
        procedure                                            :: get_nodeIndices
        procedure                                            :: get_noNodes
        procedure                                            :: get_polynomialDegree
        procedure                                            :: get_rawNodeCoordinates

        ! Setters.
        procedure :: set_elementNo
        procedure :: set_nodeCoordinates
        procedure :: set_nodeIndices
        procedure :: set_noNodes
        procedure :: set_polynomialDegree

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
    function get_elementNo(this)
        class(element), intent(inout) :: this
        integer                       :: get_elementNo

        get_elementNo = this%elementNo
    end function get_elementNo

    function get_nodeCoordinates(this)
        class(element), intent(inout) :: this
        real(dp), dimension(2)        :: get_nodeCoordinates

        get_nodeCoordinates = this%nodeCoordinates
    end function get_nodeCoordinates

    function get_nodeIndices(this)
        class(element), intent(inout) :: this
        integer, dimension(2)         :: get_nodeIndices

        get_nodeIndices = this%nodeIndices
    end function get_nodeIndices

    function get_noNodes(this)
        class(element), intent(inout) :: this
        integer                       :: get_noNodes

        get_noNodes = this%noNodes
    end function get_noNodes

    function get_polynomialDegree(this)
        class(element), intent(inout) :: this
        integer                       :: get_polynomialDegree

        get_polynomialDegree = this%polynomialDegree
    end function get_polynomialDegree

    function get_rawNodeCoordinates(this)
        class(element), intent(inout)   :: this
        real(dp), dimension(:), pointer :: get_rawNodeCoordinates

        get_rawNodeCoordinates = this%nodeCoordinates
    end function get_rawNodeCoordinates

    subroutine set_elementNo(this, a_elementNo)
        class(element), intent(inout) :: this
        integer                       :: a_elementNo

        this%elementNo = a_elementNo
    end subroutine set_elementNo

    subroutine set_nodeCoordinates(this, a_nodeCoordinates)
        class(element), intent(inout)       :: this
        real(dp), dimension(:), allocatable :: a_nodeCoordinates

        this%nodeCoordinates = a_nodeCoordinates
    end subroutine set_nodeCoordinates
        
    subroutine set_nodeIndices(this, a_nodeIndices)
        class(element), intent(inout) :: this
        integer, dimension(2)         :: a_nodeIndices

        this%nodeIndices = a_nodeIndices
    end subroutine set_nodeIndices
        
    subroutine set_noNodes(this, a_noNodes)
        class(element), intent(inout) :: this
        integer                       :: a_noNodes

        this%noNodes = a_noNodes
    end subroutine set_noNodes

    subroutine set_polynomialDegree(this, a_polynomialDegree)
        class(element), intent(inout) :: this
        integer                       :: a_polynomialDegree

        this%polynomialDegree = a_polynomialDegree
    end subroutine set_polynomialDegree
end module class_element