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
        integer, dimension(:), allocatable  :: higherDoFs
    contains
        
    end type

contains
    
end module