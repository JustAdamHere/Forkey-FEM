program test
    use common
    use quadrature
    use class_element
    use class_element_interval
    use class_mesh

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

    type(mesh) :: myMesh

    integer :: i
    real(dp), parameter :: h = 0.25

    call myMesh%constructor_ex(0)

    print *, "Adam"
    print *, myMesh%elements(2)%element_type%map_localToGlobal(0.0_dp)
    print *, myMesh%elementConnectivity(:, :)
    print *, ""

    call myMesh%destructor()
end program test