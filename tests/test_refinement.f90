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
    call mySolution%constructor(myMesh, func_one, 0.0001_dp, func_one, func_boundaryem4)

    call mySolution%solve()
    call mySolution%output_solution()

    call mySolution%destructor()
    call myMesh%destructor()

contains

end program test