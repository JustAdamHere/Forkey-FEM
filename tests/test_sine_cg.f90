program test
    use common
    use quadrature
    use class_element
    use class_element_interval
    use class_mesh
    use solvers_linear
    use class_solution
    use class_solution_cg

    implicit none

    type(mesh)        :: myMesh
    type(solution_cg) :: mySolution

    call myMesh%constructor_eq(4, 1)
    call mySolution%constructor(myMesh, func_pi2sinpi2, 1.0_dp, func_zero, func_sinpi2, func_sinpi2_)

    call mySolution%solve()
    call mySolution%output_solution_u()
    print *, mySolution%compute_L2NormDifference2()
    print *, mySolution%compute_H1NormDifference2()

    call mySolution%destructor()
    call myMesh%destructor()

contains

end program test