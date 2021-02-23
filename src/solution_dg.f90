module class_solution_dg
    use common
    use class_element
    use class_element_interval
    use class_mesh
    use class_solution
    use solvers_linear

    implicit none

    public

    type, extends(solution) :: solution_dg
        procedure(double_double), pointer, nopass :: f => null()
        real(dp)                                  :: epsilon
        procedure(double_double), pointer, nopass :: c => null()
        integer, dimension(:), allocatable        :: startDoFs
    contains
        ! Constructors.
        procedure :: constructor => solution_dg_constructor
        procedure :: destructor  => solution_dg_destructor

        ! Solvers.
        procedure :: solve => solution_dg_solve

        ! Calculators.
        procedure :: calculate_DoFs        => solution_dg_calculate_DoFs
        procedure :: calculate_elementDoFs => solution_dg_calculate_elementDoFs

        ! Getters.
        procedure :: get_typeName => solution_dg_get_typeName

        ! Internal methods.
        procedure :: a => solution_dg_a
        procedure :: l => solution_dg_l

        ! Computers.
        procedure :: compute_uh        => solution_dg_compute_uh
        procedure :: compute_uh_single => solution_dg_compute_uh_single

        ! Outputters.
        procedure :: output_solution   => solution_dg_output_solution
        procedure :: output_solution_u => solution_dg_output_solution_u
    end type

contains
    subroutine solution_dg_constructor(this, a_mesh, a_f, a_epsilon, a_c, a_u)
        class(solution_dg)       :: this
        class(mesh)              :: a_mesh
        procedure(double_double) :: a_f
        real(dp)                 :: a_epsilon
        procedure(double_double) :: a_c
        procedure(double_double) :: a_u

        this%solutionMesh =  a_mesh ! Is this making a copy?!
        this%f            => a_f
        this%epsilon      =  a_epsilon
        this%c            => a_c

        call solution_dg_calculate_DoFs(this)

        allocate(this%uh(this%DoFs))
    end subroutine

    subroutine solution_dg_destructor(this)
        class(solution_dg) :: this

        deallocate(this%startDoFs)
        deallocate(this%uh)
    end subroutine

    subroutine solution_dg_solve(this)
        class(solution_dg) :: this

        integer                                :: n, m, o
        real(dp), dimension(:, :), allocatable :: stiffnessMatrix 
        real(dp), dimension(:), allocatable    :: loadVector
        integer                                :: i1, j1
        integer                                :: i, j
        integer                                :: k, f, k1, k2
        type(element)                          :: currentElement, previousElement, nextElement
        integer                                :: p, p1, p2
        real(dp), dimension(:), allocatable    :: basis, basis1, basis2, basis1_, basis2_
        real(dp), dimension(:), allocatable    :: quadPoints, quadWeights
        real(dp), dimension(:), allocatable    :: loadVector_u0, u0
        integer, dimension(:), allocatable     :: elementDoFs
        integer, dimension(:), allocatable     :: previousElementDoFs, nextElementDoFs
        real(dp)                               :: Jacobian, Jacobian1, Jacobian2
        real(dp)                               :: theta, sigma
        real(dp)                               :: h
        real(dp)                               :: u,  u_,  v,  v_
        real(dp)                               :: um, um_, vm, vm_
        real(dp)                               :: up, up_, vp, vp_
        real(dp)                               :: b

        n = this%DoFs
        m = this%solutionMesh%noElements
        o = this%solutionMesh%noFaces

        allocate(stiffnessMatrix(n, n))
        allocate(loadVector(n))

        do k = 1, m
            currentElement = this%solutionMesh%elements(k)
            p              = currentElement%element_type%polynomialDegree

            call currentElement%element_type%get_elementQuadrature(quadPoints, quadWeights)

            call solution_dg_calculate_elementDoFs(this, k, elementDoFs)

            do j1 = 0, p
                j = elementDoFs(j1+1)
                call currentElement%element_type%evaluateBasisLegendre(basis, j1, 0, quadPoints)

                loadVector(j) = loadVector(j) + solution_dg_l(this, currentElement, basis)

                do i1 = 0, p
                    i = elementDoFs(i1+1)

                    call currentElement%element_type%evaluateBasisLegendre(basis1,  i1, 0, quadPoints)
                    call currentElement%element_type%evaluateBasisLegendre(basis2,  j1, 0, quadPoints)
                    call currentElement%element_type%evaluateBasisLegendre(basis1_, i1, 1, quadPoints)
                    call currentElement%element_type%evaluateBasisLegendre(basis2_, j1, 1, quadPoints)

                    stiffnessMatrix(i, j) = stiffnessMatrix(i, j) &
                        + solution_dg_a(this, currentElement, basis1, basis2, basis1_, basis2_)

                    deallocate(basis2_)
                    deallocate(basis1_)
                    deallocate(basis2)
                    deallocate(basis1)
                end do

                deallocate(basis)
            end do

            deallocate(elementDoFs)

            deallocate(quadPoints)
            deallocate(quadWeights)
        end do

        ! Set flux parameter.
        theta = -1

        do f = 1, o
            ! Left face.
            if (f == 1) then
                k = f

                currentElement = this%solutionMesh%elements(k)
                p              = currentElement%element_type%polynomialDegree
                Jacobian       = currentElement%element_type%get_Jacobian()

                h     = 2*Jacobian
                sigma = 10*p**2/h

                call currentElement%element_type%get_elementQuadrature(quadPoints, quadWeights)

                call solution_dg_calculate_elementDoFs(this, k, elementDoFs)

                do j1 = 0, p
                    j = elementDoFs(j1+1)

                    ! You'd need to add something to RHS if there were non-zero BCs.

                    do i1 = 0, p
                        i = elementDoFs(i1+1)

                        u  = currentElement%element_type%basisLegendre(i1, 0, -1.0_dp)
                        v  = currentElement%element_type%basisLegendre(j1, 0, -1.0_dp)
                        u_ = currentElement%element_type%basisLegendre(i1, 1, -1.0_dp)/Jacobian
                        v_ = currentElement%element_type%basisLegendre(j1, 1, -1.0_dp)/Jacobian

                        b = -this%epsilon*(u_)*(-v) + theta*this%epsilon*(v_)*(-u) + sigma*(-u)*(-v)

                        stiffnessMatrix(i, j) = stiffnessMatrix(i, j) + b
                    end do
                end do

                deallocate(elementDoFs)

                deallocate(quadPoints)
                deallocate(quadWeights)
            ! Right face.
            else if (f == o) then
                k = f-1

                currentElement = this%solutionMesh%elements(k)
                p              = currentElement%element_type%polynomialDegree
                Jacobian       = currentElement%element_type%get_Jacobian()

                h     = 2*Jacobian
                sigma = 10*p**2/h

                call currentElement%element_type%get_elementQuadrature(quadPoints, quadWeights)

                call solution_dg_calculate_elementDoFs(this, k, elementDoFs)

                do j1 = 0, p
                    j = elementDoFs(j1+1)

                    ! You'd need to add something to RHS if there were non-zero BCs.

                    do i1 = 0, p
                        i = elementDoFs(i1+1)

                        u  = currentElement%element_type%basisLegendre(i1, 0, 1.0_dp)
                        v  = currentElement%element_type%basisLegendre(j1, 0, 1.0_dp)
                        u_ = currentElement%element_type%basisLegendre(i1, 1, 1.0_dp)/Jacobian
                        v_ = currentElement%element_type%basisLegendre(j1, 1, 1.0_dp)/Jacobian

                        b = -this%epsilon*(u_)*(v) + theta*this%epsilon*(v_)*(u) + sigma*(-u)*(-v)

                        stiffnessMatrix(i, j) = stiffnessMatrix(i, j) + b
                    end do
                end do

                deallocate(elementDoFs)

                deallocate(quadPoints)
                deallocate(quadWeights)
            ! Interior face.
            else  
                k1 = f-1
                k2 = f

                previousElement = this%solutionMesh%elements(k1)
                nextElement     = this%solutionMesh%elements(k2)
                p1              = previousElement%element_type%polynomialDegree
                p2              = nextElement    %element_type%polynomialDegree
                Jacobian1       = previousElement%element_type%get_Jacobian()
                Jacobian2       = nextElement    %element_type%get_Jacobian()

                h     = 2*(Jacobian1+Jacobian2)/2
                p     = (p1+p2)/2
                sigma = 10*p**2/h

                call solution_dg_calculate_elementDoFs(this, k1, previousElementDoFs)
                call solution_dg_calculate_elementDoFs(this, k2, nextElementDoFs)

                do j1 = 0, p1
                    j = previousElementDoFs(j1+1)

                    ! --
                    do i1 = 0, p1
                        i = previousElementDoFs(i1+1)

                        um  = previousElement%element_type%basisLegendre(i1, 0, 1.0_dp)
                        vm  = previousElement%element_type%basisLegendre(j1, 0, 1.0_dp)
                        um_ = previousElement%element_type%basisLegendre(i1, 1, 1.0_dp)/Jacobian1
                        vm_ = previousElement%element_type%basisLegendre(j1, 1, 1.0_dp)/Jacobian1

                        b = -this%epsilon*(um_)*(vm)/2 + theta*this%epsilon*(vm_)*(um)/2 + sigma*(um)*(vm)

                        stiffnessMatrix(i, j) = stiffnessMatrix(i, j) + b
                    end do

                    ! +-
                    do i1 = 0, p2
                        i = nextElementDoFs(i1+1)

                        up  = nextElement    %element_type%basisLegendre(i1, 0, -1.0_dp)
                        vm  = previousElement%element_type%basisLegendre(j1, 0,  1.0_dp)
                        up_ = nextElement    %element_type%basisLegendre(i1, 1, -1.0_dp)/Jacobian2
                        vm_ = previousElement%element_type%basisLegendre(j1, 1,  1.0_dp)/Jacobian1

                        b = -this%epsilon*(up_)*(vm)/2 + theta*this%epsilon*(vm_)*(-up)/2 + sigma*(-up)*(vm)

                        stiffnessMatrix(i, j) = stiffnessMatrix(i, j) + b
                    end do
                end do

                do j1 = 0, p2
                    j = nextElementDoFs(j1+1)

                    ! -+
                    do i1 = 0, p1
                        i = previousElementDoFs(i1+1)

                        um  = previousElement%element_type%basisLegendre(i1, 0,  1.0_dp)
                        vp  = nextElement    %element_type%basisLegendre(j1, 0, -1.0_dp)
                        um_ = previousElement%element_type%basisLegendre(i1, 1,  1.0_dp)/Jacobian1
                        vp_ = nextElement    %element_type%basisLegendre(j1, 1, -1.0_dp)/Jacobian2

                        b = -this%epsilon*(um_)*(-vp)/2 + theta*this%epsilon*(vp_)*(um)/2 + sigma*(um)*(-vp)

                        stiffnessMatrix(i, j) = stiffnessMatrix(i, j) + b
                    end do

                    ! ++
                    do i1 = 0, p2
                        i = nextElementDoFs(i1+1)

                        up  = nextElement%element_type%basisLegendre(i1, 0, -1.0_dp)
                        vp  = nextElement%element_type%basisLegendre(j1, 0, -1.0_dp)
                        up_ = nextElement%element_type%basisLegendre(i1, 1, -1.0_dp)/Jacobian2
                        vp_ = nextElement%element_type%basisLegendre(j1, 1, -1.0_dp)/Jacobian2

                        b = -this%epsilon*(up_)*(-vp)/2 + theta*this%epsilon*(vp_)*(-up)/2 + sigma*(-up)*(-vp)

                        stiffnessMatrix(i, j) = stiffnessMatrix(i, j) + b
                    end do
                end do

                deallocate(nextElementDoFs)
                deallocate(previousElementDoFs)
            end if
        end do

        call direct(stiffnessMatrix, loadVector, this%uh)

        print *, this%uh

        deallocate(loadVector)
        deallocate(stiffnessMatrix)
    end subroutine

    subroutine solution_dg_calculate_elementDoFs(this, a_k, a_elementDoFs)
        class(solution_dg)                  :: this
        integer                             :: a_k
        integer, dimension(:), allocatable  :: a_elementDoFs
        
        integer :: elementDoFsStart, elementDoFsEnd
        integer :: j
        integer :: p

        p = this%solutionMesh%elements(a_k)%element_type%polynomialDegree

        allocate(a_elementDoFs(p+1))

        elementDoFsStart = this%startDoFs(a_k)
        elementDoFsEnd   = this%startDoFs(a_k+1)
        do j = 1, elementDoFsEnd-elementDoFsStart
            a_elementDoFs(j) = elementDoFsStart + j-1
        end do
    end subroutine

    ! subroutine calculate_DoFs(this) ! Um, be careful of CG and DG?!??!
    !     class(solution_dg) :: this

    !     integer :: i, m, p

    !     m = this%solutionMesh%noElements

    !     allocate(this%startDoFs(m+1))

    !     this%startDoFs(1) = 1
    !     do i = 2, m+1
    !         p = this%solutionMesh%elements(i-1)%element_type%polynomialDegree
    !         this%startDoFs(i) = this%startDoFs(i-1) + 2 + (p-1)
    !     end do

    !     this%DoFs = this%startDoFs(m+1)-1
    ! end subroutine

    ! subroutine calculate_DoFs(this)
    !     class(solution_dg) :: this

    !     integer :: i, m, p

    !     m = this%solutionMesh%noElements

    !     allocate(this%startDoFs(m+1))

    !     this%startDoFs(1) = 1
    !     do i = 2, m+1
    !         p = this%solutionMesh%elements(i-1)%element_type%polynomialDegree
    !         this%startDoFs(i) = this%startDoFs(i-1) + 1 + (p-1)
    !     end do

    !     this%DoFs = this%startDoFs(m+1)
    ! end subroutine

    subroutine solution_dg_calculate_DoFs(this)
        class(solution_dg) :: this

        integer :: i
        integer :: m
        integer :: p

        m = this%solutionMesh%noElements

        allocate(this%startDoFs(m+1))

        this%startDoFs(1) = 1
        do i = 2, m+1
            p = this%solutionMesh%elements(i-1)%element_type%polynomialDegree

            this%startDoFs(i) = this%startDoFs(i-1) + p+1
        end do

        this%DoFs = this%startDoFs(m+1) - 1
    end subroutine

    function solution_dg_get_typeName(this)
        class(solution_dg) :: this
        character(len=32)  :: solution_dg_get_typeName

        solution_dg_get_typeName = 'solution_dg'
    end function

    function solution_dg_a(this, currentElement, basis1, basis2, basis1_, basis2_)
        class(solution_dg)                  :: this
        real(dp)                            :: solution_dg_a
        type(element)                       :: currentElement
        real(dp), dimension(:), allocatable :: basis1
        real(dp), dimension(:), allocatable :: basis2
        real(dp), dimension(:), allocatable :: basis1_
        real(dp), dimension(:), allocatable :: basis2_

        real(dp)                            :: J
        real(dp), dimension(:), allocatable :: coordinates
        real(dp), dimension(:), allocatable :: weights
        integer                             :: i, n
        real(dp)                            :: bValue, cValue
        real(dp)                            :: x

        call currentElement%element_type%get_elementQuadrature(coordinates, weights)
        n = size(coordinates, 1)

        J = currentElement%element_type%get_Jacobian()
        solution_dg_a = 0

        do i = 1, n
            bValue = basis1_(i) * basis2_(i)
            solution_dg_a = solution_dg_a + this%epsilon*bValue*weights(i)/J
        end do

        do i = 1, n
            x             = currentElement%element_type%map_localToGlobal(coordinates(i))
            bValue        = basis1(i) * basis2(i)
            cValue        = this%c(x)
            solution_dg_a = solution_dg_a + cValue*bValue*weights(i)*J
        end do

        deallocate(coordinates, weights)
    end function

    function solution_dg_l(this, currentElement, basis)
        class(solution_dg)                  :: this
        real(dp)                            :: solution_dg_l
        type(element)                       :: currentElement
        real(dp), dimension(:), allocatable :: basis

        real(dp)                            :: J
        real(dp), dimension(:), allocatable :: coordinates
        real(dp), dimension(:), allocatable :: weights
        integer                             :: i, n
        real(dp)                            :: bValue, fValue
        real(dp)                            :: xValue

        call currentElement%element_type%get_elementQuadrature(coordinates, weights)
        n = size(coordinates, 1)

        J = currentElement%element_type%get_Jacobian()
        solution_dg_l = 0

        do i = 1, n
            xValue        = currentElement%element_type%map_localToGlobal(coordinates(i))
            bValue        = basis(i)
            fValue        = this%f(xValue)
            solution_dg_l = solution_dg_l + bValue*fValue*weights(i)*J
        end do

        deallocate(coordinates, weights)
    end function

    function solution_dg_compute_uh_single(this, a_elementNo, a_deriv, a_xi)
        class(solution_dg) :: this
        real(dp)           :: solution_dg_compute_uh_single
        integer            :: a_elementNo
        integer            :: a_deriv
        real(dp)           :: a_xi

        real(dp), dimension(:), allocatable :: xi
        class(element), pointer             :: currentElement
        real(dp)                            :: Jacobian
        integer                             :: startDoF, endDoF
        integer                             :: i, i1, j
        integer                             :: n
        integer                             :: p
        real(dp), dimension(:), allocatable :: basis
        integer, dimension(:), allocatable  :: elementDoFs

        currentElement => this%solutionMesh%elements(a_elementNo)
        Jacobian = currentElement%element_type%get_Jacobian()
        p        = currentElement%element_type%polynomialDegree

        allocate(xi(1))
        xi(1) = a_xi

        solution_dg_compute_uh_single = 0

        call solution_dg_calculate_elementDoFs(this, a_elementNo, elementDoFs)

        do i1 = 0, p
            i = elementDoFs(i1+1)
            call currentElement%element_type%evaluateBasisLegendre(basis, i1, a_deriv, xi)

            solution_dg_compute_uh_single = solution_dg_compute_uh_single + this%uh(i)*basis(1)

            deallocate(basis)
        end do

        deallocate(elementDoFs)
        deallocate(xi)

        solution_dg_compute_uh_single = solution_dg_compute_uh_single / (Jacobian ** a_deriv)

    end function

    function solution_dg_compute_uh(this, a_elementNo, a_deriv, a_xi)
        class(solution_dg)                  :: this
        real(dp), dimension(:), allocatable :: solution_dg_compute_uh
        integer                             :: a_elementNo
        integer                             :: a_deriv
        real(dp), dimension(:), allocatable :: a_xi

        class(element), pointer             :: currentElement
        real(dp)                            :: Jacobian
        integer                             :: startDoF, endDoF
        integer                             :: i, i1, j
        integer                             :: n
        integer                             :: p
        real(dp), dimension(:), allocatable :: basis
        integer, dimension(:), allocatable  :: elementDoFs

        currentElement => this%solutionMesh%elements(a_elementNo)
        Jacobian = currentElement%element_type%get_Jacobian()
        p        = currentElement%element_type%polynomialDegree

        n = size(a_xi, 1)
        allocate(solution_dg_compute_uh(n))
        solution_dg_compute_uh = 0

        call solution_dg_calculate_elementDoFs(this, a_elementNo, elementDoFs)

        do i1 = 0, p
            i = elementDoFs(i1+1)
            call currentElement%element_type%evaluateBasisLegendre(basis, i1, a_deriv, a_xi)

            do j = 1, n
                solution_dg_compute_uh(j) = solution_dg_compute_uh(j) + this%uh(i)*basis(j)
            end do

            deallocate(basis)
        end do

        deallocate(elementDoFs)

        solution_dg_compute_uh = solution_dg_compute_uh / (Jacobian ** a_deriv)
    end function

    subroutine solution_dg_output_solution(this)
        class(solution_dg) :: this

        integer           :: fileNo
        character(len=32) :: fileName
        integer           :: i, j, n
        integer           :: noSamples
        real(dp)          :: h
        real(dp)          :: x, xi, uh
        integer           :: x1

        n = this%solutionMesh%noElements

        noSamples = 20
        h = 2.0_dp/noSamples

        fileNo = 1
        fileName = './data/solution.dat'

        open(fileNo, file=fileName, status='old')

        write(fileNo, *) 0.0_dp, this%compute_uh_single(1, 0, -1.0_dp)
        do i = 1, n
            do j = 1, noSamples
                xi = -1 + j*h
                x = this%solutionMesh%elements(i)%element_type%map_localToGlobal(xi)

                uh = this%compute_uh_single(i, 0, xi)

                write(fileNo, *) x, uh
            end do
        end do

        close(fileNo)
    end subroutine

    subroutine solution_dg_output_solution_u(this, a_function)
        class(solution_dg)       :: this
        procedure(double_double) :: a_function

        integer           :: fileNo
        character(len=32) :: fileName
        integer           :: i, j, n
        integer           :: noSamples
        real(dp)          :: h
        real(dp)          :: x, xi, uh
        integer           :: x1

        n = this%solutionMesh%noElements

        noSamples = 20
        h = 2.0_dp/noSamples

        fileNo = 1
        fileName = './data/solution.dat'

        open(fileNo, file=fileName, status='old')

        write(fileNo, *) 0.0_dp, this%compute_uh_single(1, 0, -1.0_dp), a_function(0.0_dp)
        do i = 1, n
            do j = 1, noSamples
                xi = -1 + j*h
                x = this%solutionMesh%elements(i)%element_type%map_localToGlobal(xi)

                uh = this%compute_uh_single(i, 0, xi)

                write(fileNo, *) x, uh, a_function(x)
            end do
        end do

        close(fileNo)
    end subroutine

end module