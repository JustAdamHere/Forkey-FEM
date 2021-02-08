module class_solution
    use common
    use class_element
    use class_element_interval
    use class_mesh

    implicit none

    public

    type, abstract :: solution
        type(mesh), pointer                 :: mesh
        real(dp), dimension(:), allocatable :: solution
        integer                             :: DoFs
        integer, dimension(:), allocatable  :: startDoFs
    contains
        
    end type

contains
    
end module