module class_point
    use common

    implicit none

    type point
        integer :: dimension
        real(dp), dimension(:), allocatable :: coordinates
    contains
        procedure :: point_constructor
        procedure :: point_destructor
    end type point

contains
    subroutine point_constructor(this, dim)
        class(point) :: this
        integer      :: dim

        this%dimension = dim
        allocate(this%coordinates(dim))
    end subroutine

    subroutine point_destructor(this)
        class(point) :: this

        deallocate(this%coordinates)
    end subroutine

end module class_point