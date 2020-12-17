module quadrature
    use common

    implicit none

contains
    recursive function legendrePolynomial(a_degree, a_deriv, a_point) result(LP)
        integer  :: a_degree
        integer  :: a_deriv
        real(dp) :: a_point
        real(dp) :: LP

        real(dp) :: temp1
        real(dp) :: temp2
        real(dp) :: temp3

        if (-1 <= a_point .and. a_point <= 1) then
            if (a_degree == 0) then
                if (a_deriv == 0) then
                    LP = 1
                else
                    LP = 0
                end if
            else if (a_degree == 1) then
                if (a_deriv == 0) then
                    LP = a_point
                else if (a_deriv == 1) then
                    LP = 1
                else
                    LP = 0
                end if
            else
                if (a_deriv == 0) then
                    temp1 = (2*a_degree-1)*a_point/a_degree * legendrePolynomial(a_degree-1, a_deriv, a_point)
                    temp2 = (a_degree-1)/a_degree           * legendrePolynomial(a_degree-2, a_deriv, a_point)
                    LP = temp1 - temp2
                else
                    temp1 = (2*a_degree-1)*a_point/(a_degree-1)               * legendrePolynomial(a_degree-1, a_deriv,   a_point)
                    temp2 = real(a_degree, dp)/(a_degree-1)                   * legendrePolynomial(a_degree-2, a_deriv,   a_point)
                    temp3 = (a_deriv-1) * real(2*a_degree-1, dp)/(a_degree-1) * legendrePolynomial(a_degree-1, a_deriv-1, a_point)
                    LP = (temp1 - temp2)*temp3
                end if
            end if
        else
            LP = 0
        end if
    end function
end module quadrature