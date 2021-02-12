module quadrature
    use common

    implicit none

    ! Should add cache for Gauss-Legendre points and weights.

contains
    recursive function LegendrePolynomial(a_degree, a_deriv, a_point) result(LP)
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
                    temp2 = real(a_degree-1, dp)/a_degree   * legendrePolynomial(a_degree-2, a_deriv, a_point)
                    LP = temp1 - temp2
                else
                    temp1 = real(2*a_degree-1, dp)*a_point/(a_degree-1)             * & 
                        legendrePolynomial(a_degree-1, a_deriv,   a_point)
                    temp2 = real(a_degree, dp)/(a_degree-1)                         * & 
                        legendrePolynomial(a_degree-2, a_deriv,   a_point)
                    temp3 = real(a_deriv-1, dp)*real(2*a_degree-1, dp)/(a_degree-1) * & 
                        legendrePolynomial(a_degree-1, a_deriv-1, a_point)
                    LP = temp1 - temp2 + temp3
                end if
            end if
        else
            LP = 0
        end if
    end function

    function LegendrePolynomialRoot(a_degree, a_rootNo)
        integer  :: a_degree
        integer  :: a_rootNo
        real(dp) :: LegendrePolynomialRoot

        integer             :: iterations    = 0
        real(dp), parameter :: tolerance     = 1e-15
        integer, parameter  :: maxIterations = 1e4   ! <-- This might want to be removed? 

        LegendrePolynomialRoot = -cos((2.0_dp*a_rootNo - 1)/(2.0_dp*a_degree)*pi)

        do while (abs(LegendrePolynomial(a_degree, 0, LegendrePolynomialRoot)) >= tolerance .and. iterations < maxIterations)
            LegendrePolynomialRoot = LegendrePolynomialRoot - LegendrePolynomial(a_degree, 0, LegendrePolynomialRoot)/ &
               LegendrePolynomial(a_degree, 1, LegendrePolynomialRoot)
            iterations = iterations + 1
        end do
    end function

    function GaussLegendrePoint(a_degree, a_pointNo)
        integer  :: a_degree
        integer  :: a_pointNo
        real(dp) :: GaussLegendrePoint

        GaussLegendrePoint = LegendrePolynomialRoot(a_degree, a_pointNo)
    end function

    function GaussLegendreWeight(a_degree, a_weightNo)
        integer  :: a_degree
        integer  :: a_weightNo
        real(dp) :: GaussLegendreWeight

        real(dp) :: point
        point = GaussLegendrePoint(a_degree, a_weightNo)

        GaussLegendreWeight = 2.0_dp / ((1-point**2) * LegendrePolynomial(a_degree, 1, point)**2)
    end function
end module quadrature