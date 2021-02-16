module class_solution
    use common
    use class_element
    use class_element_interval
    use class_mesh

    implicit none

    public

    type, abstract :: solution
        type(mesh), pointer                 :: solutionMesh => null()
        real(dp), dimension(:), allocatable :: uh
        integer                             :: DoFs
    contains
        procedure(interface_compute_uh), deferred        :: compute_uh
        procedure(interface_compute_uh_single), deferred :: compute_uh_single
        procedure(interface_output_solution), deferred   :: output_solution
        procedure(interface_get_typeName), deferred      :: get_typeName
    end type

    abstract interface
        function interface_compute_uh(this, a_elementNo, a_deriv, a_xi)
            use common
            import solution

            class(solution)                     :: this
            real(dp), dimension(:), allocatable :: interface_compute_uh
            integer                             :: a_elementNo
            integer                             :: a_deriv
            real(dp), dimension(:), allocatable :: a_xi
        end function
    end interface

    abstract interface
        function interface_compute_uh_single(this, a_elementNo, a_deriv, a_xi)
            use common
            import solution

            class(solution) :: this
            real(dp)        :: interface_compute_uh_single
            integer         :: a_elementNo
            integer         :: a_deriv
            real(dp)        :: a_xi
        end function
    end interface

    abstract interface
        subroutine interface_output_solution(this)
            import solution

            class(solution) :: this
        end subroutine
    end interface

    abstract interface
        function interface_get_typeName(this)
            import solution

            class(solution)   :: this
            character(len=32) :: interface_get_typeName
        end function
    end interface

contains
    ! subroutine output_solution_vtk(this)
    !     class(solution) :: this

    !     integer           :: fileNo
    !     character(len=32) :: fileName
    !     integer           :: i

    !     fileNo = 1
    !     fileName = './data/solution.vtk'

    !     open(fileNo, file=fileName, status='old')

    !     !!!!!!!!!!!!!
    !     ! FILE INFO !
    !     !!!!!!!!!!!!!
    !     write(fileNo, *) "# vtk DataFile Version 3.1"
    !     write(fileNo, *) "<<<<<< Mesh and Solution - VTK Format >>>>>>"
    !     write(fileNo, *) "ASCII"
    !     write(fileNo, *) "DATASET UNSTRUCTURED_GRID"

    !     !!!!!!!!!!!!!!!
    !     ! MESH POINTS !
    !     !!!!!!!!!!!!!!!
    !     write(fileNo, *) "POINTS", this%solutionMesh%noElements+1, "float"
    !     do i = 1, this%solutionMesh%noElements+1
    !         write(fileNo, *) this%solutionMesh%nodeCoordinates(i), 0.0_dp, 0.0_dp
    !     end do
    !     write(fileNo, *)

    !     !!!!!!!!!
    !     ! CELLS !
    !     !!!!!!!!!
    !     write(fileNo, *) "CELLS", this%solutionMesh%noElements, 64
    !     do i = 1, this%solutionMesh%noElements
    !         write(fileNo, *) 3, this%solutionMesh%elements(i)%element_type%nodeIndices(1), &
    !             this%solutionMesh%elements(i)%element_type%nodeIndices(2)
    !     end do
    !     write(fileNo, *)

    !     write(fileNo, *) "CELL_TYPES", this%solutionMesh%noElements
    !     do i = 1, this%solutionMesh%noElements
    !         write(fileNo, *) 3
    !     end do

    !     write(fileNo, *) "CELL_DATA", this%solutionMesh%noElements
    !     write(fileNo, *) "SCALARS ", "Element_Region_Id ", "integer 1"
    !     write(fileNo, *) "LOOKUP_TABLE ", "default"
    !     do i = 1, this%solutionMesh%noElements
    !         write(fileNo, *) 0
    !     end do
    !     write(fileNo, *)

    !     !!!!!!!!!!!!!!
    !     ! POINT DATA !
    !     !!!!!!!!!!!!!!
    !     write(fileNo, *) "POINT_DATA", this%solutionMesh%noElements+1
    !     write(fileNo, *) "SCALARS ", "u double 1 "
    !     write(fileNo, *) "LOOKUP_TABLE default"
    !     !write(fileNo, *) compute_uh_single(this, 1, 0, -1.0_dp)
    !     write(fileNo, *) this%compute_uh_single(1, 0, -1.0_dp)
    !     do i = 1, this%solutionMesh%noElements
    !         !write(fileNo, *) compute_uh_single(this, i, 0, 1.0_dp)
    !         write(fileNo, *) this%compute_uh_single(i, 0, 1.0_dp)
    !     end do

    !     close(fileNo)
    ! end subroutine
end module