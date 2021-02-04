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

    ! type(element_interval)              :: my_element
    ! integer                             :: elementNo
    ! integer                             :: noNodes
    ! integer, dimension(2)               :: nodeIndices
    ! real(dp), dimension(:), allocatable :: nodeCoordinates
    ! integer                             :: polynomialDegree

    ! elementNo = 1
    ! noNodes = 2
    ! nodeIndices(1) = 1
    ! nodeIndices(2) = 2
    ! polynomialDegree = 1

    ! allocate(nodeCoordinates(2))
    ! nodeCoordinates(1) = 0
    ! nodeCoordinates(2) = 1

    ! call my_element%constructor(elementNo, noNodes, nodeIndices, nodeCoordinates, polynomialDegree)
    ! print *, my_element%get_Jacobian()

    ! call my_element%destructor()

    ! deallocate(nodeCoordinates)



    ! type(mesh)        :: myMesh
    ! type(solution_cg) :: mySolution

    ! real(dp) :: epsilon

    ! call myMesh%constructor_ex(0)
    ! call mySolution%constructor(myMesh, f, epsilon, c)

    ! print *, "Adam"
    ! print *, myMesh%elements(2)%element_type%map_localToGlobal(0.0_dp)
    ! print *, myMesh%elementConnectivity(:, :)
    ! print *, ""

    ! call mySolution%destructor()
    ! call myMesh%destructor()



    ! |==========================================================
    ! |                       EXAMPLE                           |
    ! |==========================================================
    
    !     5  0  3  2
    !     0  3  0  1
    !     0  0  0  2
    !    -1  1  0  1
    
    ! |=========================================================|

    real(dp), dimension(:, :), allocatable :: myMatrix
    real(dp), dimension(:), allocatable    :: myVector
    real(dp), dimension(:), allocatable    :: mySolution

    integer n

    n = 4

    allocate(myMatrix(n, n))
    allocate(myVector(n))
    allocate(mySolution(n))

    myMatrix       = 0
    myMatrix(1, 1) = 5
    myMatrix(1, 4) = -1
    myMatrix(2, 2) = 3
    myMatrix(2, 4) = 1
    myMatrix(3, 1) = 3
    myMatrix(3, 3) = 0
    myMatrix(4, 1) = 2
    myMatrix(4, 2) = 1
    myMatrix(4, 3) = 2
    myMatrix(4, 4) = 1

    myVector(1) = 1
    myVector(2) = 2
    myVector(3) = 3
    myVector(4) = 4

    mySolution = 0

    !call GaussJordan(myMatrix, myVector, mySolution)
    !call GMRES(myMatrix, myVector, mySolution, 200, 1e-15_dp)
    call direct(myMatrix, myVector, mySolution)

    print *, myMatrix
    print *, mySolution

    print *, matmul(myMatrix, mySolution)

    deallocate(mySolution)
    deallocate(myVector)
    deallocate(myMatrix)

contains

    function f(x)
        use common

        real(dp) :: x
        real(dp) :: f
        
        f = x
    end function

    function c(x)
        use common

        real(dp) :: x
        real(dp) :: c
        
        c = x
    end function

end program test