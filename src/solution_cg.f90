module class_solution_cg
    use common
    use class_element
    use class_element_interval
    use class_mesh
    use class_solution

    implicit none

    public

    type, extends(solution) :: solution_cg
        procedure(double_double), pointer, nopass :: f => null()
        procedure(double_double), pointer, nopass :: c => null()
    contains
        ! Constructors.
        procedure :: constructor
        procedure :: destructor

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
    subroutine constructor(this, a_mesh, a_f, a_epsilon, a_c)
        class(solution_cg)       :: this
        class(mesh)              :: a_mesh
        procedure(double_double) :: a_f
        real(dp)                 :: a_epsilon
        procedure(double_double) :: a_c

        this%f => a_f
    end subroutine

    subroutine destructor(this)
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