module refinement
    use common
    use class_element
    use class_element_interval
    use class_mesh
    use class_solution
    use class_solution_cg

contains
    subroutine refinement_refine(a_mesh, a_meshNew, a_solution, a_solutionNew, a_adaptivityTolerance, a_adaptivityMaxIterations, &
        a_refineh, a_refinep, a_refineGlobally)
        type(mesh)               :: a_mesh
        type(mesh), pointer      :: a_meshNew
        class(solution)          :: a_solution
        class(solution), pointer :: a_solutionNew
        real(dp)                 :: a_adaptivityTolerance
        integer                  :: a_adaptivityMaxIterations
        logical                  :: a_refineh
        logical                  :: a_refinep
        logical                  :: a_refineGlobally
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