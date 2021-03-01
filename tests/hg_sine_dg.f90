program test
    use common
    use quadrature
    use class_element
    use class_element_interval
    use class_mesh
    use solvers_linear
    use class_solution
    use class_solution_dg
    use refinement

    implicit none

    type(mesh), target         :: myMesh
    type(solution_dg), target  :: mySolution
    type(mesh), pointer        :: myNewMesh
    type(solution_dg), pointer :: myNewSolution

    ! allocate(myMesh)
    ! allocate(mySolution)

    call myMesh%constructor_eq(2, 1)
    call mySolution%constructor(myMesh, func_pi2sinpi2, 1.0_dp, func_zero, func_sinpi2, func_sinpi2_)

    call refinement_refine(myMesh, myNewMesh, mySolution, myNewSolution, 1e-15_dp, 5, .true., .false., .true., .true.)

    call myNewSolution%output_solution_u() ! (already solved in refinement procedure)

    call myNewSolution%destructor()
    call myNewMesh%destructor()
    call mySolution%destructor()
    call myMesh%destructor()

    deallocate(myNewSolution)
    deallocate(myNewMesh)
contains

end program test