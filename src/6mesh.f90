module class_mesh
    use common
    use class_element

    implicit none

    type mesh
        private
            integer :: noElements
            class(element), dimension(:), allocatable :: elements
    end type mesh

contains

end module class_mesh