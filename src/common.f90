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
end module common