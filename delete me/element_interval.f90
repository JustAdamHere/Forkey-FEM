module class_element_interval
    use common
    use class_element

    implicit none

    ! element_interval type, derived from element, with one implemented procedure.
    type, extends(element) :: element_interval

    contains
        procedure :: basisLegendre2
    end type element_interval

contains
    ! Definition of this one procedure.
    function basisLegendre2(this)
        class(element_interval), intent(inout) :: this
        real(dp)                               :: basisLegendre2

        basisLegendre2 = 1.0_dp ! This gives an error like below:
        !Result mismatch for the overriding procedure ‘basislegendre’ at (1): Type mismatch in function result (REAL(8)/INTEGER(4))
    end function
end module class_element_interval