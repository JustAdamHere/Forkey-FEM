module class_point
    use common

    implicit none

    type point
        private
            integer :: dimension
            real(dp), dimension(:), allocatable :: coordinates
    end type point

contains

end module class_point