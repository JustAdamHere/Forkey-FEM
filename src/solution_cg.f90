module class_solution_cg
    use common
    use class_element
    use class_element_interval
    use class_mesh
    use class_solution

    implicit none

    public

    type, extends(solution) :: solution_cg
        procedure(double_double), pointer, nopass :: f
        procedure(double_double), pointer, nopass :: c
    contains
        ! Constructors.
        procedure :: constructor

        ! Solvers.
        procedure :: solve

        ! Internal methods.
        procedure :: a
        procedure :: l
    end type

    interface
        function double_double(x)
            use common

            real(dp) :: x
            real(dp) :: double_double
        end function
    end interface

contains
    subroutine constructor(this)
        class(solution_cg) :: this
    end subroutine

    subroutine solve(this)
        class(solution_cg) :: this
    end subroutine

    function a(this)
        class(solution_cg) :: this
        real(dp)           :: a
    end function

    function l(this)
        class(solution_cg) :: this
        real(dp)           :: l
    end function



end module