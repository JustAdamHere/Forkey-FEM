program test
    use common
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

    call myMesh%constructor(1)

    call myMesh%destructor()
end program test