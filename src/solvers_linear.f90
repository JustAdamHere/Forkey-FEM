module solvers_linear
    use common
    use class_matrix
    use class_matrix_full
    use class_matrix_sparse

    implicit none

contains
    ! Algorithm based upon https://en.wikipedia.org/wiki/Generalized_minimal_residual_method.
    subroutine GMRES(a_matrix, a_vector, a_solution, a_maxIterations, a_tolerance)
        real(dp), dimension(:, :), allocatable :: a_matrix
        real(dp), dimension(:), allocatable    :: a_vector
        real(dp), dimension(:), allocatable    :: a_solution
        integer                                :: a_maxIterations
        real(dp)                               :: a_tolerance

        real(dp), dimension(:), allocatable    :: y

        integer :: dim
        integer :: n, k, j
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
        real(dp), dimension(:), allocatable    :: beta_sub

        ! Sets the residual, using a_solution as initial vector.
        r = a_vector - a_vector * a_solution

        ! Dimension of problem.
        dim = size(a_vector)

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
        Q(:, 1) = r/rNorm ! <- Problem here
        beta    = rNorm * e1

        k = 0
        stopIterating = .false.
        do while (k <= a_maxIterations .and. .not. stopIterating)
            k = k + 1

            ! Arnoldi.
            call Arnoldi(a_matrix, Q, k, a_maxIterations, H(1:k+1, k), Q(:, k+1))

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
                a(j, i) = a_matrix(j, i)
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
                    d = a(j, i) / a(i, i)
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

        print *, inverse

        ! print *, size(inverse, 1)
        ! print *, size(inverse, 2)
        ! print *, size(a_vector)
        do i = 1, n
            do j = 1, n
                a_solution(i) = a_solution(i) + inverse(j, i)*a_vector(j)
            end do    
        end do
        !a_solution = inverse*a_vector
        !a_solution = a_vector*inverse
    end subroutine
    
end module solvers_linear