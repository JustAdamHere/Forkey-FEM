program matrix_test
    use common
    use class_matrix
    use class_matrix_full
    use class_matrix_sparse
    use solvers_linear

    implicit none

    real(dp), dimension(:, :), allocatable :: fulEntries
    real(dp), dimension(:), allocatable    :: spaEntries
    integer, dimension(:), allocatable     :: spaColNo
    integer, dimension(:), allocatable     :: spaRowSt

    integer :: n
    integer :: m

    type(matrix_full)   :: matFul
    type(matrix_sparse) :: matSpa
    type(matrix_sparse) :: matSpa2

    n = 5
    m = 9

    allocate(fulEntries(n, n))
    allocate(spaEntries(m))
    allocate(spaColNo(m))
    allocate(spaRowSt(n+1))

    ! |==========================================================
    ! |                       EXAMPLE                           |
    ! |==========================================================
    ! 
    !     5  0  3  2
    !     0  3  0  1
    !     0  0  0  2
    !    -1  1  0  1
    !
    ! |=========================================================|

    fulEntries       = 0
    fulEntries(1, 1) = 5
    fulEntries(1, 3) = 3
    fulEntries(1, 4) = 2
    fulEntries(2, 2) = 3
    fulEntries(2, 4) = 1
    fulEntries(3, 4) = 2
    fulEntries(4, 1) = -1
    fulEntries(4, 2) = 1
    fulEntries(4, 4) = 1

    spaEntries = (/5, 3, 2, 3, 1, 2, -1, 1, 1 /)
    spaColNo   = (/1, 3, 4, 2, 4, 4,  1, 2, 4 /)
    spaRowSt   = (/1, 4, 6, 7, 10 /)

    call matFul%constructor(fulEntries)
    call matSpa%constructor(spaEntries, spaColNo, spaRowSt)

    call matFul%destructor()
    call matSpa%destructor()

    deallocate(spaRowSt)
    deallocate(spaColNo)
    deallocate(spaEntries)
    deallocate(fulEntries)

end program matrix_test