module solvers_linear
    use common
    use class_matrix
    use class_matrix_full
    use class_matrix_sparse

    implicit none

contains
    ! Algorithm based upon https://en.wikipedia.org/wiki/Generalized_minimal_residual_method.
    subroutine GMRESOLD(a_matrix, a_vector, a_solution, a_maxIterations, a_tolerance)
        real(dp), dimension(:, :), allocatable :: a_matrix
        real(dp), dimension(:), allocatable    :: a_vector
        real(dp), dimension(:), allocatable    :: a_solution
        integer                                :: a_maxIterations
        real(dp)                               :: a_tolerance

        real(dp), dimension(:), allocatable    :: y

        integer :: dim
        integer :: n, i, j, k
        logical :: stopIterating

        real(dp), dimension(:), allocatable :: r, beta
        real(dp)                            :: rNorm, vectorNorm, error

        real(dp), dimension(a_maxIterations)   :: sn
        real(dp), dimension(a_maxIterations)   :: cs
        real(dp), dimension(a_maxIterations+1) :: e
        real(dp), dimension(a_maxIterations+1) :: e1

        real(dp), dimension(:, :), allocatable                    :: Q
        real(dp), dimension(a_maxIterations+1, a_maxIterations)   :: H

        real(dp), dimension(:, :), allocatable :: H_sub
        real(dp), dimension(:, :), allocatable :: Q_sub
        real(dp), dimension(:), allocatable    :: beta_sub

        ! Sets the residual, using a_solution as initial vector.
        r = a_vector - a_vector * a_solution

        ! Dimension of problem.
        dim = size(r)

        ! Allocates Q.
        allocate(Q(dim, a_maxIterations+1))

        ! Initialises norms.
        vectorNorm = norm(a_vector)
        error      = norm(r)/vectorNorm

        ! Gives fake value to initialise the vector of errors.
        e = -1

        ! Initialise the 1D vectors.
        sn      = 0
        cs      = 0
        e1      = 0
        e1(1)   = 1
        e (1)   = error
        rNorm   = norm(r)
        Q(:, 1) = r/rNorm
        beta    = rNorm * e1

        k = 0
        stopIterating = .false.
        do while (k <= a_maxIterations .and. .not. stopIterating)
            k = k + 1

            ! Arnoldi.
            call ArnoldiOLD(a_matrix, Q, k, a_maxIterations, H(1:k+1, k), Q(:, k+1))

            ! Eliminate the last element in H ith row and update the rotation matrix.
            call applyGivensRotation(H(1:k+1, k), cs, sn, k, a_maxIterations)

            ! Update the residual vector.
            beta(k+1) = -sn(k) * beta(k)
            beta(k)   = cs(k)  * beta(k)
            error     = abs(beta(k+1)) / vectorNorm

            ! Save the error.
            e(k+1) = error

            ! Stops if error is under the tolerance.
            if (error <= a_tolerance) then
                stopIterating = .true.
            end if
        end do

        ! Calculates the result.
        allocate(H_sub(k, k))
        H_sub = H(1:k, 1:k)
        allocate(beta_sub(k))
        beta_sub = beta(1:k)
        allocate(y(k))
        call GaussJordan(H_sub, beta_sub, y)
        ! y = H(1:k, 1:k) \ beta(1:k)
        
        !a_solution = a_solution + Q(1:k, 1:k) * y
        a_solution = a_solution + matmul(Q(1:k, 1:k), y)
    end subroutine

    subroutine GMRES(a_matrix, a_vector, a_solution, a_maxIterations, a_tolerance)
        ! Input variables.
        real(dp), dimension(:, :), allocatable :: a_matrix
        real(dp), dimension(:), allocatable    :: a_vector
        real(dp), dimension(:), allocatable    :: a_solution
        integer                                :: a_maxIterations
        real(dp)                               :: a_tolerance

        ! Subroutine variables.
        real(dp), dimension(:), allocatable    :: r0
        real(dp)                               :: beta
        real(dp), dimension(:), allocatable    :: wj
        real(dp), dimension(:, :), allocatable :: V
        real(dp), dimension(:, :), allocatable :: Hm
        real(dp), dimension(:), allocatable    :: ym
        real(dp), dimension(:), allocatable    :: e1
        integer                                :: i, j
        integer                                :: dim
        integer                                :: m
        real(dp), dimension(:), allocatable    :: cs
        real(dp), dimension(:), allocatable    :: sn

        ! Aux.
        dim = size(a_vector)
        m   = a_maxIterations
        allocate(V(m+1, dim)) ! <- Is this right?!
        allocate(e1(m+1))
        e1 = 0
        e1(1) = 1

        ! STEP 1
        r0      = a_vector - matmul(a_matrix, a_solution)
        beta    = norm(r0)
        V(1, :) = r0/beta ! v1

        ! STEP 2
        allocate(Hm(m+1, m))
        Hm = 0

        ! STEPS 3—11
        do j = 1, m
            wj = matmul(a_matrix, V(i, :))
            do i = 1, m
                Hm(i, j) = dot_product(wj, V(i, :))
                wj       = wj - Hm(i, j)*V(i, :)
            end do

            ! STEP 9—10
            Hm(j+1, j) = norm(wj)
            if (Hm(j+1, j) == 0) then
                m = j
            end if
            V(j+1, :) = wj/Hm(j+1, j)
        end do

        ! STEP 12
        !Rotations and stuff


        ! Deallocation.
        deallocate(Hm)
        deallocate(V)
    end subroutine

    ! "Stolen" from AptoFEM.
    subroutine direct(a_matrix, a_vector, a_solution)
        real(dp), dimension(:, :), allocatable :: a_matrix
        real(dp), dimension(:), allocatable    :: a_vector
        real(dp), dimension(:), allocatable    :: a_solution

        integer                                :: n
        real(dp), dimension(:, :), allocatable :: a
        real(dp), dimension(:), allocatable    :: x
        integer,  dimension(:), allocatable    :: permutation
        integer                                :: i

        n = size(a_matrix, 1)
        allocate(permutation(n), x(n), a(n, n))
        a = a_matrix

        ! LU factorisation.
        call LUFact(a, permutation)

        ! Solve ax=b.
        x = a_vector
        call fwd(a, permutation, x)
        call bkd(a, permutation, x)

        do i = 1, n
            a_solution(permutation(i)) = x(i)
        end do

        deallocate(permutation, x, a)
    end subroutine

    ! "Stolen" from AptoFEM.
    subroutine LUFact(a_matrix, a_permutation)
        ! Input variables.
        real(dp), dimension(:, :), allocatable :: a_matrix
        integer,  dimension(:),    allocatable :: a_permutation

        ! Auxiliary variables.
        integer  :: n
        real(dp) :: s, s1, v
        integer  :: k, len, p1, j1, i, j

        n = size(a_matrix, 1)

        do k = 1, n
            a_permutation(k) = k
        end do

        do k = 1, n-1
            len = n-k
            s = abs(a_matrix(k, a_permutation(k)))

            j1 = k
            do j = k+1, n
                s1 = abs(a_matrix(k, a_permutation(j)))
                if (s1 > s) then
                    j1 = j
                    s  = s1
                end if
            end do

            if (s == 0.0_dp) then
                write(1023, *) ' Singular, possibly zero matrix in LUFact.'
                do i = 1, n
                    write(1023, *) (a_matrix(i, j), j=1, n)
                end do
                stop
            end if

            if (j1 /= k) then
                p1      = a_permutation(k)
                a_permutation(k)  = a_permutation(j1)
                a_permutation(j1) = p1
            end if

            v = 1.0_dp/a_matrix(k, a_permutation(k))
            a_matrix(k+1:k+len, a_permutation(k)) = v*a_matrix(k+1:k+len, a_permutation(k))

            do j = k+1, n
                a_matrix(k+1:k+len, a_permutation(j)) = a_matrix(k+1:k+len, a_permutation(j)) &
                    - a_matrix(k, a_permutation(j))*a_matrix(k+1:k+len, a_permutation(k))
            end do
        end do
    end subroutine

    ! "Stolen" from AptoFEM.
    subroutine fwd(a_matrix, a_permutation, a_solution)
        ! Input variables.
        real(dp), dimension(:, :), allocatable, intent(in)    :: a_matrix
        integer,  dimension(:),    allocatable, intent(in)    :: a_permutation
        real(dp), dimension(:),    allocatable, intent(inout) :: a_solution

        ! Auxiliary variables.
        integer :: n, j

        n = size(a_matrix, 1)

        do j = 1, n-1
            a_solution(j+1:n) = a_solution(j+1:n) - a_solution(j)*a_matrix(j+1:n, a_permutation(j))
        end do
    end subroutine

    ! "Stolen" from AptoFEM.
    subroutine bkd(a_matrix, a_permutation, a_solution)
        ! Input variables.
        real(dp), dimension(:, :), allocatable, intent(in)    :: a_matrix
        integer,  dimension(:),    allocatable, intent(in)    :: a_permutation
        real(dp), dimension(:),    allocatable, intent(inout) :: a_solution

        ! Auxiliary variables.
        integer :: n, k

        n = size(a_matrix, 1)

        do k = n, 1, -1
            a_solution(k) = a_solution(k)/a_matrix(k, a_permutation(k))
            a_solution(1:k-1) = a_solution(1:k-1) - a_solution(k)*a_matrix(1:k-1, a_permutation(k))
        end do
    end subroutine

    subroutine GaussJordan(a_matrix, a_vector, a_solution)
        real(dp), dimension(:, :), allocatable :: a_matrix
        real(dp), dimension(:), allocatable    :: a_vector
        real(dp), dimension(:), allocatable    :: a_solution

        integer :: n, i, j, k
        real(dp) :: d
        real(dp), dimension(:, :), allocatable :: a
        real(dp), dimension(:, :), allocatable :: inverse

        n = size(a_matrix, 1) ! Assumes square, as it should be!

        allocate(a(2*n, n))
        a = 0

        ! Sets the coefficients of the matrix correctly.
        do i = 1, n
            do j = 1, n
                a(i, j) = a_matrix(j, i)
            end do
        end do

        ! Sets identity reference matrix.
        do i = 1, n
            a(i+n, i) = 1
        end do

        ! Reducing to diagonal matrix.
        do i = 1, n
            do j = 1, n
                if (j /= i) then
                    d = a(i, j) / a(i, i)
                    !print *, a(i, j), " ", a(i, i), " ", d
                    do k = 1, 2*n
                        a(k, j) = a(k, j) - d*a(k, i)
                    end do
                end if
            end do
        end do

        ! Reducing to unit matrix.
        do i = 1, n
            d = a(i, i)
            do j = 1, 2*n
                a(j, i) = a(j, i)/d
            end do
        end do

        !print *, a

        ! Setting the inverse.
        allocate(inverse(n, n))
        do i = 1, n
            do j = 1, n
                inverse(j, i) = a(j+n, i)
            end do
        end do

        a_solution = matmul(inverse, a_vector)

        ! print *, size(inverse, 1)
        ! print *, size(inverse, 2)
        ! print *, size(a_vector)

        ! do i = 1, n
        !     do j = 1, n
        !         a_solution(i) = a_solution(i) + inverse(j, i)*a_vector(j)
        !     end do    
        ! end do

        !a_solution = inverse*a_vector
        !a_solution = a_vector*inverse
    end subroutine
    
end module solvers_linear