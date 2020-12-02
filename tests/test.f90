program test
    use common
    use class_element
    use class_element_interval

    implicit none

    type(element_interval)              :: my_element
    integer                             :: elementNo
    integer                             :: noNodes
    integer, dimension(2)               :: nodeIndices
    real(dp), dimension(:), allocatable :: nodeCoordinates
    integer                             :: polynomialDegree

    elementNo = 1
    noNodes = 2
    nodeIndices(1) = 1
    nodeIndices(2) = 2
    polynomialDegree = 1

    allocate(nodeCoordinates(2))
    nodeCoordinates(1) = 0
    nodeCoordinates(2) = 1

    my_element = element_interval(elementNo, noNodes, nodeIndices, nodeCoordinates, polynomialDegree)
    print *, my_element%elementNo

    deallocate(nodeCoordinates)
end program test