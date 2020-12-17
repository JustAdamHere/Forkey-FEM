module class_element
    use common

    implicit none

    ! Make only element visible outside the module.
    private
    public :: element

    ! Abstract element type, with one deferred procedure.
    type, abstract :: element
    contains
        procedure(interface_basisLegendre), deferred :: basisLegendre2
    end type element

    ! Interface for deferred procedure.
    abstract interface
        function interface_basisLegendre(this)
            use common
            import element

            class(element), intent(inout) :: this
            real(dp)                      :: interface_basisLegendre
        end function
    end interface

contains
    
end module class_element