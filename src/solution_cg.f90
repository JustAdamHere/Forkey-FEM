module class_solution_cg
    use common
    use class_element
    use class_element_interval
    use class_mesh
    use class_solution
    use solvers_linear

    implicit none

    public

    type, extends(solution) :: solution_cg
        procedure(double_double), pointer, nopass :: f => null()
        real(dp)                                  :: epsilon
        procedure(double_double), pointer, nopass :: c => null()
        integer, dimension(:), allocatable        :: startHigherDoFs
    contains
        ! Constructors.
        procedure :: constructor => solution_cg_constructor
        procedure :: destructor  => solution_cg_destructor

        ! Solvers.
        procedure :: solve => solution_cg_solve

        ! Calculators.
        procedure :: calculate_DoFs        => solution_cg_calculate_DoFs
        procedure :: calculate_elementDoFs => solution_cg_calculate_elementDoFs

        ! Getters.
        procedure :: get_typeName => solution_cg_get_typeName

        ! Internal methods.
        procedure :: a => solution_cg_a
        procedure :: l => solution_cg_l

        ! Computers.
        procedure :: compute_uh        => solution_cg_compute_uh
        procedure :: compute_uh_single => solution_cg_compute_uh_single

        ! Outputters.
        procedure :: output_solution   => solution_cg_output_solution
        procedure :: output_solution_u => solution_cg_output_solution_u
    end type

    ! interface
    !     function integer_integer(x)
    !         use common

    !         integer :: x
    !         integer :: integer_integer
    !     end function
    ! end interface

    ! interface
    !     subroutine nothing()

    !     end subroutine
    ! end interface

    ! interface
    !     function interface_basis(this, a_degree, a_deriv, a_point)
    !         use common
    !         use class_element

    !         class(element)                   :: this
    !         integer                                :: a_degree
    !         integer                                :: a_deriv
    !         real(dp)                               :: a_point
    !     end function
    ! end interface

contains
    subroutine solution_cg_constructor(this, a_mesh, a_f, a_epsilon, a_c, a_u, a_u_1)
        class(solution_cg)       :: this
        class(mesh)              :: a_mesh
        procedure(double_double) :: a_f
        real(dp)                 :: a_epsilon
        procedure(double_double) :: a_c
        procedure(double_double) :: a_u
        procedure(double_double) :: a_u_1

        this%solutionMesh =  a_mesh ! Is this making a copy?!
        this%f            => a_f
        this%epsilon      =  a_epsilon
        this%c            => a_c
        this%u            => a_u
        this%u_1          => a_u_1

        call solution_cg_calculate_DoFs(this)

        allocate(this%uh(this%DoFs))
    end subroutine

    subroutine solution_cg_destructor(this)
        class(solution_cg) :: this

        deallocate(this%startHigherDoFs)
        deallocate(this%uh)
    end subroutine

    subroutine solution_cg_solve(this)
        class(solution_cg) :: this

        integer                                :: n, m
        real(dp), dimension(:, :), allocatable :: stiffnessMatrix 
        real(dp), dimension(:), allocatable    :: loadVector
        integer                                :: i1, j1
        integer                                :: i, j, k
        type(element)                          :: currentElement
        integer                                :: p
        real(dp), dimension(:), allocatable    :: basis, basis1, basis2, basis1_, basis2_
        real(dp), dimension(:), allocatable    :: quadPoints, quadWeights
        real(dp), dimension(:), allocatable    :: loadVector_u0, u0
        integer, dimension(:), allocatable     :: elementDoFs

        real(dp) :: stiffnessNorm
        real(dp) :: loadNorm

        n = this%DoFs
        m = this%solutionMesh%noElements

        allocate(stiffnessMatrix(n, n))
        allocate(loadVector(n))

        stiffnessMatrix = 0
        loadVector      = 0

        do k = 1, m
            currentElement = this%solutionMesh%elements(k)
            p              = currentElement%element_type%polynomialDegree

            call currentElement%element_type%get_elementQuadrature(quadPoints, quadWeights)

            call solution_cg_calculate_elementDoFs(this, k, elementDoFs)

            do j1 = 0, p
                j = elementDoFs(j1+1)
                call currentElement%element_type%evaluateBasisLobatto(basis, j1, 0, quadPoints)

                loadVector(j) = loadVector(j) + solution_cg_l(this, currentElement, basis)

                do i1 = 0, p
                    i = elementDoFs(i1+1)

                    call currentElement%element_type%evaluateBasisLobatto(basis1,  i1, 0, quadPoints)
                    call currentElement%element_type%evaluateBasisLobatto(basis2,  j1, 0, quadPoints)
                    call currentElement%element_type%evaluateBasisLobatto(basis1_, i1, 1, quadPoints)
                    call currentElement%element_type%evaluateBasisLobatto(basis2_, j1, 1, quadPoints)

                    stiffnessMatrix(i, j) = stiffnessMatrix(i, j) &
                        + solution_cg_a(this, currentElement, basis1, basis2, basis1_, basis2_)

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

        ! Sets boundary conditions.
        !!! Only works in 1D !!!
        loadVector(1)   = 0
        loadVector(m+1) = 0

        stiffnessMatrix(1,   :)   = 0
        stiffnessMatrix(:,   1)   = 0
        stiffnessMatrix(m+1, :)   = 0
        stiffnessMatrix(:,   m+1) = 0

        allocate(u0(n))
        u0      = 0
        u0(1)   = 0
        u0(m+1) = 0

        allocate(loadVector_u0(n))

        loadVector_u0 = matmul(stiffnessMatrix, u0)

        loadVector = loadVector - loadVector_u0

        deallocate(loadVector_u0)
        deallocate(u0)

        stiffnessMatrix(1,   1)   = 1
        stiffnessMatrix(m+1, m+1) = 1

        call direct(stiffnessMatrix, loadVector, this%uh)



        loadNorm = norm(loadVector)

        stiffnessNorm = 0
        do k = 1, n
            !print *, stiffnessMatrix(k, :)
            stiffnessNorm = stiffnessNorm + norm(stiffnessMatrix(k, :))
        end do
        stiffnessNorm = sqrt(stiffnessNorm)

        print *, stiffnessNorm
        print *, loadNorm
        print *, norm(this%uh)

        ! print *, ""
        ! print *, loadVector

        ! print *, ""
        ! print *, this%uh




        deallocate(loadVector)
        deallocate(stiffnessMatrix)
    end subroutine

    subroutine solution_cg_calculate_elementDoFs(this, a_k, a_elementDoFs)
        class(solution_cg)                  :: this
        integer                             :: a_k
        integer, dimension(:), allocatable  :: a_elementDoFs
        
        integer                             :: elementHigherDoFsStart, elementHigherDoFsEnd
        integer                             :: j
        integer                             :: p

        p = this%solutionMesh%elements(a_k)%element_type%polynomialDegree

        allocate(a_elementDoFs(p+1))
        a_elementDoFs(1) = a_k
        a_elementDoFs(2) = a_k+1

        elementHigherDoFsStart = this%startHigherDoFs(a_k)
        elementHigherDoFsEnd   = this%startHigherDoFs(a_k+1)
        do j = 1, elementHigherDoFsEnd-elementHigherDoFsStart
            a_elementDoFs(j+2) = elementHigherDoFsStart + j-1
        end do
    end subroutine

    ! subroutine calculate_DoFs(this) ! Um, be careful of CG and DG?!??!
    !     class(solution_cg) :: this

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
    !     class(solution_cg) :: this

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

     subroutine solution_cg_calculate_DoFs(this)
        class(solution_cg) :: this

        integer :: i, m, p

        m = this%solutionMesh%noElements

        allocate(this%startHigherDoFs(m+1))

        this%startHigherDoFs(1) = m+2
        do i = 2, m+1
            p = this%solutionMesh%elements(i-1)%element_type%polynomialDegree
            ! this%startDoFs(i) = this%startDoFs(i-1) + 1 + (p-1)
            this%startHigherDoFs(i) = this%startHigherDoFs(i-1) + (p-1)
        end do

        this%DoFs = this%startHigherDoFs(m+1) - 1
    end subroutine

    function solution_cg_get_typeName(this)
        class(solution_cg) :: this
        character(len=32)  :: solution_cg_get_typeName

        solution_cg_get_typeName = 'solution_cg'
    end function

    function solution_cg_a(this, currentElement, basis1, basis2, basis1_, basis2_)
        class(solution_cg)                  :: this
        real(dp)                            :: solution_cg_a
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
        solution_cg_a = 0

        do i = 1, n
            bValue = basis1_(i) * basis2_(i)
            solution_cg_a = solution_cg_a + this%epsilon*bValue*weights(i)/J
        end do

        do i = 1, n
            x             = currentElement%element_type%map_localToGlobal(coordinates(i))
            bValue        = basis1(i) * basis2(i)
            cValue        = this%c(x)
            solution_cg_a = solution_cg_a + cValue*bValue*weights(i)*J
        end do

        deallocate(coordinates, weights)
    end function

    function solution_cg_l(this, currentElement, basis)
        class(solution_cg)                  :: this
        real(dp)                            :: solution_cg_l
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
        solution_cg_l = 0

        do i = 1, n
            xValue        = currentElement%element_type%map_localToGlobal(coordinates(i))
            bValue        = basis(i)
            fValue        = this%f(xValue)
            solution_cg_l = solution_cg_l + bValue*fValue*weights(i)*J
        end do

        deallocate(coordinates, weights)
    end function

    function solution_cg_compute_uh_single(this, a_elementNo, a_deriv, a_xi)
        class(solution_cg) :: this
        real(dp)           :: solution_cg_compute_uh_single
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

        solution_cg_compute_uh_single = 0

        call solution_cg_calculate_elementDoFs(this, a_elementNo, elementDoFs)

        do i1 = 0, p
            i = elementDoFs(i1+1)
            call currentElement%element_type%evaluateBasisLobatto(basis, i1, a_deriv, xi)

            solution_cg_compute_uh_single = solution_cg_compute_uh_single + this%uh(i)*basis(1)

            deallocate(basis)
        end do

        deallocate(elementDoFs)
        deallocate(xi)

        solution_cg_compute_uh_single = solution_cg_compute_uh_single / (Jacobian ** a_deriv)

    end function

    function solution_cg_compute_uh(this, a_elementNo, a_deriv, a_xi)
        class(solution_cg)                  :: this
        real(dp), dimension(:), allocatable :: solution_cg_compute_uh
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
        allocate(solution_cg_compute_uh(n))
        solution_cg_compute_uh = 0

        call solution_cg_calculate_elementDoFs(this, a_elementNo, elementDoFs)

        do i1 = 0, p
            i = elementDoFs(i1+1)
            call currentElement%element_type%evaluateBasisLobatto(basis, i1, a_deriv, a_xi)

            do j = 1, n
                solution_cg_compute_uh(j) = solution_cg_compute_uh(j) + this%uh(i)*basis(j)
            end do

            deallocate(basis)
        end do

        deallocate(elementDoFs)

        solution_cg_compute_uh = solution_cg_compute_uh / (Jacobian ** a_deriv)

    end function

    subroutine solution_cg_output_solution(this)
        class(solution_cg) :: this

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

    subroutine solution_cg_output_solution_u(this)
        class(solution_cg)       :: this

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

        write(fileNo, *) 0.0_dp, this%compute_uh_single(1, 0, -1.0_dp), this%u(0.0_dp)
        do i = 1, n
            do j = 1, noSamples
                xi = -1 + j*h
                x = this%solutionMesh%elements(i)%element_type%map_localToGlobal(xi)

                uh = this%compute_uh_single(i, 0, xi)

                write(fileNo, *) x, uh, this%u(x)
            end do
        end do

        close(fileNo)
    end subroutine

end module