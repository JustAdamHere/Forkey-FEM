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
        a_refineh, a_refinep, a_refineGlobally, a_output)
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

        type(mesh), pointer        :: oldMesh
        type(mesh), pointer        :: newMesh
        type(solution_dg), pointer :: oldSolution
        type(solution_dg), pointer :: newSolution
        integer                    :: noElements
        integer                    :: noPolys
        integer                    :: i
        real(dp)                   :: prevL2Norm, nextL2Norm, ratioL2Norm
        real(dp)                   :: prevH1Norm, nextH1Norm, ratioH1Norm

        integer           :: fileNo
        character(len=32) :: fileName

        oldMesh     => a_mesh
        oldSolution => a_solution

        if (a_refineGlobally .and. a_refineh) then
            if (a_output) then
                fileNo   = 2
                fileName = './data/convergence.dat'

                open(fileNo, file=fileName, status='old')

                print *, "       DoFs", "   L2 Ratio", "                  H1 Ratio"
            end if

            if (a_output) then
                prevL2Norm  = 0
                prevH1Norm  = 0
                nextL2Norm  = sqrt(oldSolution%compute_L2NormDifference2())
                nextH1Norm  = sqrt(oldSolution%compute_H1NormDifference2())
                ratioL2Norm = prevL2Norm/nextL2Norm
                ratioH1Norm = prevH1Norm/nextH1Norm

                write(fileNo, *) oldSolution%DoFs, nextL2Norm, nextH1Norm

                print *, oldSolution%DoFs, ratioL2Norm, ratioH1Norm
            end if

            !-------+-------------------|
            ! START | Refine h uniformly.
            !-------+-------------------|

            allocate(newMesh)
            allocate(newSolution)

            noElements = oldMesh%noElements
            noPolys    = oldMesh%elements(1)%element_type%polynomialDegree ! <-- Assumes constant polynomial degree.

            call newMesh%constructor_eq(2*noElements, noPolys)
            call newSolution%constructor(newMesh, oldSolution%f, oldSolution%epsilon, oldSolution%c, oldSolution%u, oldSolution%u_1)
            call newSolution%solve()

            oldMesh     => newMesh
            oldSolution => newSolution

            !-----+-------------------|
            ! END | Refine h uniformly.
            !-----+-------------------|

            if (a_output) then
                prevL2Norm  = nextL2Norm
                prevH1Norm  = nextH1Norm
                nextL2Norm  = sqrt(newSolution%compute_L2NormDifference2())
                nextH1Norm  = sqrt(newSolution%compute_H1NormDifference2())
                ratioL2Norm = prevL2Norm/nextL2Norm
                ratioH1Norm = prevH1Norm/nextH1Norm

                write(fileNo, *) newSolution%DoFs, nextL2Norm, nextH1Norm

                print *, newSolution%DoFs, ratioL2Norm, ratioH1Norm
            end if

            do i = 2, a_adaptivityMaxIterations
                !-------+-------------------|
                ! START | Refine h uniformly.
                !-------+-------------------|

                allocate(newMesh)
                allocate(newSolution)

                noElements = oldMesh%noElements
                noPolys    = oldMesh%elements(1)%element_type%polynomialDegree ! <-- Assumes constant polynomial degree.

                call newMesh%constructor_eq(2*noElements, noPolys)
                call newSolution%constructor(newMesh, oldSolution%f, oldSolution%epsilon, oldSolution%c, oldSolution%u, &
                    oldSolution%u_1)
                call newSolution%solve()

                deallocate(oldMesh)
                deallocate(oldSolution)

                oldMesh     => newMesh
                oldSolution => newSolution

                !-------+-------------------|
                ! END | Refine h uniformly.
                !-------+-------------------|

                if (a_output) then
                    prevL2Norm  = nextL2Norm
                    prevH1Norm  = nextH1Norm
                    nextL2Norm  = sqrt(newSolution%compute_L2NormDifference2())
                    nextH1Norm  = sqrt(newSolution%compute_H1NormDifference2())
                    ratioL2Norm = prevL2Norm/nextL2Norm
                    ratioH1Norm = prevH1Norm/nextH1Norm

                    write(fileNo, *) newSolution%DoFs, nextL2Norm, nextH1Norm

                    print *, newSolution%DoFs, ratioL2Norm, ratioH1Norm
                end if
            end do

            if (a_output) then
                close(fileNo)
            end if
        end if

        a_meshNew     => newMesh
        a_solutionNew => newSolution
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