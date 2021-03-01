module common
    implicit none

    integer, parameter  :: dp = selected_real_kind(15)
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

    interface
        function interface_basis(a_degree, a_deriv, a_point)
            import dp

            integer                                :: a_degree
            integer                                :: a_deriv
            real(dp)                               :: a_point
        end function
    end interface

    interface
        function double_double(x)
            import dp

            real(dp) :: x
            real(dp) :: double_double
        end function
    end interface

contains
    function appendToArray(a_array, a_value)
        real(dp), dimension(:), allocatable :: a_array
        real(dp)                            :: a_value
        real(dp), dimension(:), allocatable :: appendToArray

        integer :: n

        n = size(a_array)

        allocate(appendToArray(n+1))

        appendToArray(1:n) = a_array
        appendToArray(n+1) = a_value

        deallocate(a_array)
    end function

    function norm(a_array)
        real(dp), dimension(:) :: a_array
        real(dp)               :: norm

        integer :: i

        ! Calculates Euclidian norm.
        norm = 0
        do i = 1, size(a_array)
            norm = norm + abs(a_array(i)) 
        end do
        norm = sqrt(norm)
    end function

    subroutine ArnoldiOLD(a_A, a_Q, a_k, a_m, h, q)
        real(dp), dimension(:, :)  :: a_A
        real(dp), dimension(:, :)  :: a_Q
        integer                    :: a_k
        integer                    :: a_m
        real(dp), dimension(a_k+1) :: h
        real(dp), dimension(a_m)   :: q

        integer :: i

        q = matmul(a_A, a_Q(:, a_k))
        do i = 1, a_k
            h(i) = dot_product(q, a_Q(:, i))
            q    = q - h(i) * a_Q(:, i)
        end do
        h(a_k+1) = norm(q)
        q        = q / h(a_k+1)
    end subroutine

    ! Taken from page 37:
    !  Tim Kelley,
    !  Iterative Methods for Linear and Nonlinear Equations,
    !  SIAM, 2004,
    !  ISBN: 0898713528,
    !  LC: QA297.8.K45.
    subroutine Arnoldi(A, b, x0, k, V)
        ! Input variables.
        real(dp), dimension(:, :), allocatable :: A
        real(dp), dimension(:),    allocatable :: b
        real(dp), dimension(:),    allocatable :: x0
        integer                                :: k
        real(dp), dimension(:),    allocatable :: V

       
    end subroutine

    subroutine applyGivensRotation(a_h, a_cs, a_sn, a_k, a_m)
        integer :: a_k
        integer :: a_m
        real(dp), dimension(a_k) :: a_h
        real(dp), dimension(a_m) :: a_cs
        real(dp), dimension(a_m) :: a_sn

        real(dp) :: temp
        integer  :: i

        do i = 1, a_k-1
            temp     =  a_cs(i)*a_h(i) + a_sn(i)*a_h(i+1)
            a_h(i+1) = -a_sn(i)*a_h(i) + a_cs(i)*a_h(i+1)
            a_h(i)   = temp
        end do

        call GivensRotation(a_h(a_k), a_h(a_k+1), a_cs(a_k), a_sn(a_k))

        ! Eliminate H(i+1, i)
        a_h(a_k)   = a_cs(a_k)*a_h(a_k) + a_sn(a_k)*a_h(a_k+1)
        a_h(a_k+1) = 0

    end subroutine

    subroutine GivensRotation(a_v1, a_v2, cs, sn)
        real(dp) :: a_v1
        real(dp) :: a_v2
        real(dp) :: cs
        real(dp) :: sn

        real(dp) :: t

        t  = sqrt(a_v1**2 + a_v2**2)
        cs = a_v1/t
        sn = a_v2/t
    end subroutine

    subroutine evaluateBasis(a_basisPoints, a_basis, a_degree, a_deriv, a_points)
        real(dp), dimension(:), allocatable :: a_basisPoints
        procedure(interface_basis), pointer :: a_basis
        integer                             :: a_degree
        integer                             :: a_deriv
        real(dp), dimension(:), allocatable :: a_points

        allocate(a_basisPoints(size(a_points, 1)))
    end subroutine

    function func_sin(x)
        real(dp) :: func_sin
        real(dp) :: x

        func_sin = sin(x)
    end function

    function func_one(x)
        real(dp) :: func_one
        real(dp) :: x

        func_one = 1
    end function

    function func_zero(x)
        real(dp) :: func_zero
        real(dp) :: x

        func_zero = 0
    end function

    function func_boundaryem4(x)
        real(dp) :: func_boundaryem4
        real(dp) :: x

        real(dp) :: a

        a = 1e-4

        func_boundaryem4 = -exp(x/sqrt(a))/(exp(1.0_dp/sqrt(a)) + 1) - (exp(-x/sqrt(a)) &
            * exp(1.0_dp/sqrt(a)))/(exp(1.0_dp/sqrt(a)) + 1) + 1
    end function

    function func_boundaryem4_(x)
        real(dp) :: func_boundaryem4_
        real(dp) :: x

        real(dp) :: a

        a = 1e-4

        func_boundaryem4_ = -exp(x/sqrt(a))/(exp(1.0_dp/sqrt(a)) + 1)/sqrt(a) &
        + (exp(-x/sqrt(a)) * exp(1.0_dp/sqrt(a)))/(exp(1.0_dp/sqrt(a)) + 1)/sqrt(a)
    end function

    function func_sinpi2(x)
        real(dp) :: func_sinpi2
        real(dp) :: x

        func_sinpi2 = sin(2.0_dp*pi*x)
    end function

    function func_sinpi2_(x)
        real(dp) :: func_sinpi2_
        real(dp) :: x

        func_sinpi2_ = 2.0_dp*pi * sin(2.0_dp*pi*x)
    end function

    function func_pi2sinpi2(x)
        real(dp) :: func_pi2sinpi2
        real(dp) :: x

        func_pi2sinpi2 = (2.0_dp*pi)**2 * sin(2.0_dp*pi*x)
    end function
end module common