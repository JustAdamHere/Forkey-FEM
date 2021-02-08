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

        procedure :: testadam
    end type

    interface
        function double_double(x)
            use common

            real(dp) :: x
            real(dp) :: double_double
        end function
    end interface

contains
    subroutine testadam(this)
        class(solution_cg) :: this
    end subroutine

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

        call this%testadam()
        !call this%calculate_DoFs()

        this%DoFs = 4

        n = this%DoFs
        allocate(this%startDoFs(n))

        allocate(this%u(n))
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
        integer                                :: i, j, k
        type(element)                          :: currentElement
        integer                                :: p

        n = this%DoFs
        m = this%solutionMesh%noElements

        allocate(stiffnessMatrix(n, n))
        allocate(loadVector(n))

        stiffnessMatrix = 0
        loadVector      = 0

        do k = 1, m
            currentElement = this%solutionMesh%elements(k)
            p              = currentElement%element_type%polynomialDegree

            
        end do


        deallocate(loadVector)
        deallocate(stiffnessMatrix)
    end subroutine

    subroutine calculate_DoFs(this)
        class(solution_cg) :: this

        integer :: m

        !m = this%solutionMesh%noElements

        !allocate(this%startDoFs(m+1))

        !this%DoFs = 10
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