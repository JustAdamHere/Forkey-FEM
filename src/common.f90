module common
    implicit none

    integer, parameter :: dp = selected_real_kind(15)
    integer, parameter :: pi = 4.0_dp * atan(1.0_dp)

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

    subroutine Arnoldi(a_A, a_Q, a_k, a_m, h, q)
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
end module common