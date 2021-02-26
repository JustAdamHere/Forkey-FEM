module refinement
    use common
    use class_element
    use class_element_interval
    use class_mesh
    use class_solution
    use class_solution_cg
    use class_solution_dg

contains
    ! subroutine refinement_refine(a_mesh, a_meshNew, a_solution, a_solutionNew, a_adaptivityTolerance, a_adaptivityMaxIterations, &
    !     a_refineh, a_refinep, a_refineGlobally, a_output, a_u, a_u_1)
    !     type(mesh), pointer        :: a_mesh
    !     type(mesh), pointer        :: a_meshNew
    !     type(solution_dg), pointer :: a_solution    ! <- Need to change to class(solution).
    !     type(solution_dg), pointer :: a_solutionNew ! <- Need to change to class(solution).
    !     real(dp)                   :: a_adaptivityTolerance
    !     integer                    :: a_adaptivityMaxIterations
    !     logical                    :: a_refineh
    !     logical                    :: a_refinep
    !     logical                    :: a_refineGlobally
    !     logical                    :: a_output
    !     procedure(double_double)   :: a_u
    !     procedure(double_double)   :: a_u_1

    !     type(mesh), pointer        :: oldMesh
    !     type(mesh), pointer        :: newMesh
    !     type(solution_dg), pointer :: oldSolution
    !     type(solution_dg), pointer :: newSolution
    !     integer                    :: noElements
    !     integer                    :: noPolys
    !     integer                    :: i

    !     oldMesh     => a_mesh
    !     oldSolution => a_solution
    !     allocate(newMesh)
    !     allocate(newSolution)

    !     if (a_refineGlobally .and. a_refineh) then
    !         do i = 1, a_adaptivityMaxIterations
    !             noElements = oldMesh%noElements
    !             noPolys    = oldMesh%elements(1)%element_type%polynomialDegree ! <-- Assumes constant polynomial degree.

    !             !call newMesh%constructor_eq(2*noElements, noPolys)
    !             call newMesh%constructor_eq(8, 1)
    !             call newSolution%constructor(newMesh, oldSolution%c, oldSolution%epsilon, oldSolution%f, func_zero)

    !             a_meshNew     => newMesh
    !             a_solutionNew => newSolution
    !         end do
    !     end if
    ! end subroutine

    subroutine refinement_refine(a_mesh, a_meshNew, a_solution, a_solutionNew, a_adaptivityTolerance, a_adaptivityMaxIterations, &
        a_refineh, a_refinep, a_refineGlobally, a_output, a_u, a_u_1)
        type(mesh), target         :: a_mesh
        type(mesh), pointer        :: a_meshNew
        type(solution_dg), target  :: a_solution    ! <- Need to change to class(solution).
        type(solution_dg), pointer :: a_solutionNew ! <- Need to change to class(solution).
        real(dp)                   :: a_adaptivityTolerance
        integer                    :: a_adaptivityMaxIterations
        logical                    :: a_refineh
        logical                    :: a_refinep
        logical                    :: a_refineGlobally
        logical                    :: a_output
        procedure(double_double)   :: a_u
        procedure(double_double)   :: a_u_1

        type(mesh), pointer        :: oldMesh
        type(mesh), pointer        :: newMesh
        type(solution_dg), pointer :: oldSolution
        type(solution_dg), pointer :: newSolution
        integer                    :: noElements
        integer                    :: noPolys
        integer                    :: i

        oldMesh     => a_mesh
        oldSolution => a_solution
        allocate(newMesh)
        allocate(newSolution)

        if (a_refineGlobally .and. a_refineh) then
            !do i = 1, a_adaptivityMaxIterations
                noElements = oldMesh%noElements
                noPolys    = oldMesh%elements(1)%element_type%polynomialDegree ! <-- Assumes constant polynomial degree.

                call newMesh%constructor_eq(2*noElements, noPolys)
                call newSolution%constructor(newMesh, oldSolution%f, oldSolution%epsilon, oldSolution%c, func_sinpi2)

                a_meshNew     => newMesh
                a_solutionNew => newSolution
            !end do
        end if
    end subroutine

    subroutine refinement_refine_h(a_mesh, a_meshNew, a_solution, a_solutionNew, a_errorIndicators)
        type(mesh)                          :: a_mesh
        type(mesh), pointer                 :: a_meshNew
        class(solution)                     :: a_solution
        class(solution), pointer            :: a_solutionNew
        real(dp), dimension(:), allocatable :: a_errorIndicators
    end subroutine

    subroutine refinement_refine_p(a_mesh, a_meshNew, a_solution, a_solutionNew, a_errorIndicators)
        type(mesh)                          :: a_mesh
        type(mesh), pointer                 :: a_meshNew
        class(solution)                     :: a_solution
        class(solution),pointer             :: a_solutionNew
        real(dp), dimension(:), allocatable :: a_errorIndicators
    end subroutine

    subroutine refinement_refine_hp(a_mesh, a_meshNew, a_solution, a_solutionNew, a_errorIndicators)
        type(mesh)                          :: a_mesh
        type(mesh), pointer                 :: a_meshNew
        class(solution)                     :: a_solution
        class(solution),pointer             :: a_solutionNew
        real(dp), dimension(:), allocatable :: a_errorIndicators
    end subroutine
end module