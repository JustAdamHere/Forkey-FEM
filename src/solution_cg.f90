module class_solution_cg
    use common
    use class_element
    use class_element_interval
    use class_mesh
    use class_solution

    implicit none

    public

    type, extends(solution) :: solution_cg
        procedure(double_double), pointer, nopass :: f => null()
        real(dp)                                  :: epsilon
        procedure(double_double), pointer, nopass :: c => null()
        real(dp), dimension(:), allocatable       :: u 
        type(mesh), pointer                       :: solutionMesh => null()
    contains
        ! Constructors.
        procedure :: constructor
        procedure :: destructor

        ! Solvers.
        procedure :: solve

        ! Calculators.
        procedure :: calculate_DoFs

        ! Internal methods.
        procedure :: a
        procedure :: l
    end type

    interface
        function double_double(x)
            use common

            real(dp) :: x
            real(dp) :: double_double
        end function
    end interface

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
    subroutine constructor(this, a_mesh, a_f, a_epsilon, a_c)
        class(solution_cg)       :: this
        class(mesh)              :: a_mesh
        procedure(double_double) :: a_f
        real(dp)                 :: a_epsilon
        procedure(double_double) :: a_c

        integer :: n

        this%solutionMesh =  a_mesh ! Is this making a copy?!
        this%f            => a_f
        this%epsilon      =  a_epsilon
        this%c            => a_c

        n = 5

        allocate(this%u(n))

        call calculate_DoFs(this)
    end subroutine

    subroutine destructor(this)
        class(solution_cg) :: this

        deallocate(this%startDoFs)
        deallocate(this%u)
    end subroutine

    subroutine solve(this)
        class(solution_cg) :: this

        integer                                :: n, m
        real(dp), dimension(:, :), allocatable :: stiffnessMatrix 
        real(dp), dimension(:), allocatable    :: loadVector
        integer                                :: a, b
        integer                                :: i, j, k
        type(element)                          :: currentElement
        integer                                :: p
        integer                                :: startDoF, endDoF
        procedure(interface_basis), pointer    :: basisType1
        procedure(interface_basis), pointer    :: basisType2
        real(dp), dimension(:), allocatable    :: basis, basis1, basis2, basis1_, basis2_
        real(dp), dimension(:), allocatable    :: quadPoints, quadWeights

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

            startDoF = this%startDoFs(k)
            endDoF   = this%startDoFs(k+1)
            do b = 0, endDoF-startDoF-1
                j     = startDoF + b
                !basisType2 => currentElement%element_type%basisLobatto

                call evaluateBasis(basis, basisType2, b, 0, quadPoints) ! <-- Write this into a routine in element.

                !loadVector(j) = loadVector(j) + this%l(currentElement, basis)

                do a = 0, endDoF-startDoF-1
                    i = startDoF + a
                    !basisType1 => currentElement%element_type%basisLobatto

                    call evaluateBasis(basis1,  basisType1, a, 0, quadPoints)
                    call evaluateBasis(basis2,  basisType2, b, 0, quadPoints)
                    call evaluateBasis(basis1_, basisType1, a, 1, quadPoints)
                    call evaluateBasis(basis2_, basisType2, b, 1, quadPoints)


                    deallocate(basis2_)
                    deallocate(basis1_)
                    deallocate(basis2)
                    deallocate(basis1)
                end do

                deallocate(basis)
            end do

            deallocate(quadPoints)
            deallocate(quadWeights)
        end do


        deallocate(loadVector)
        deallocate(stiffnessMatrix)
    end subroutine

    subroutine calculate_DoFs(this)
        class(solution_cg) :: this

        integer :: i, m, p

        m = this%solutionMesh%noElements

        allocate(this%startDoFs(m+1))

        this%startDoFs(1) = 1
        do i = 2, m+1
            p = this%solutionMesh%elements(i-1)%element_type%polynomialDegree
            this%startDoFs(i) = this%startDoFs(i-1) + 2 + (p-1)
        end do

        this%DoFs = this%startDoFs(m+1)-1
    end subroutine

    function a(this, currentElement, basis1, basis2, basis1_, basis2_)
        class(solution_cg)       :: this
        real(dp)                 :: a
        type(element)            :: currentElement
        procedure(double_double) :: basis1
        procedure(double_double) :: basis2
        procedure(double_double) :: basis1_
        procedure(double_double) :: basis2_

        real(dp)                            :: J
        real(dp), dimension(:), allocatable :: coordinates
        real(dp), dimension(:), allocatable :: weights
        integer                             :: i, n
        real(dp)                            :: b, c
        real(dp)                            :: x

        call currentElement%element_type%get_elementQuadrature(coordinates, weights)
        n = size(coordinates, 1)

        J = currentElement%element_type%get_Jacobian()
        a = 0

        do i = 1, n
            b = basis1_(coordinates(i)) * basis2_(coordinates(i))
            a = a + this%epsilon*b*weights(i)/J
        end do

        do i = 1, n
            x = currentElement%element_type%map_localToGlobal(coordinates(i))
            b = basis1(coordinates(i)) * basis2(coordinates(i))
            c = this%c(x)
            a = a + c*b*weights(i)*J
        end do

        deallocate(coordinates, weights)
    end function

    function l(this, currentElement, basis)
        class(solution_cg)       :: this
        real(dp)                 :: l
        type(element)            :: currentElement
        procedure(double_double) :: basis

        real(dp)                            :: J
        real(dp), dimension(:), allocatable :: coordinates
        real(dp), dimension(:), allocatable :: weights
        integer                             :: i, n
        real(dp)                            :: b, f
        real(dp)                            :: x

        call currentElement%element_type%get_elementQuadrature(coordinates, weights)
        n = size(coordinates, 1)

        J = currentElement%element_type%get_Jacobian()
        l = 0

        do i = 1, n
            x = currentElement%element_type%map_localToGlobal(coordinates(i))
            b = basis(coordinates(i))
            f = this%f(x)
            l = l + b*f*weights(i)*J
        end do

        deallocate(coordinates, weights)
    end function



end module