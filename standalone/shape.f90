module shape_mod
type shape
    integer :: color
    logical :: filled
    integer :: x
    integer :: y
contains
    procedure :: initialize => initShape
end type shape

type, extends(shape) :: rectangle
    integer :: length
    integer :: width
contains
    procedure :: initialize => initRectangle
end type rectangle

type, extends(rectangle) :: square
end type square

contains

    subroutine initShape(this, color, filled, x, y, length, width)
        ! initialize shape objects
        class(shape) :: this
        integer :: color
        logical :: filled
        integer :: x
        integer :: y
        integer, optional :: length  ! ingnored for shape
        integer, optional :: width   ! ignored for shape

        this%color = color
        this%filled = filled
        this%x = x
        this%y = y
    end subroutine initShape

    subroutine initRectangle(this, color, filled, x, y, length, width)
        ! initialize rectangle objects
        class(rectangle) :: this
        integer :: color
        logical :: filled
        integer :: x
        integer :: y
        integer, optional :: length  
        integer, optional :: width   

        this%color = color
        this%filled = filled
        this%x = x
        this%y = y

        if (present(length)) then
           this%length = length
        else
           this%length = 0
        endif
        if (present(width)) then 
            this%width = width + 1
        else
             this%width = 0
        endif

    end subroutine initRectangle
    
end module

program shape_test
  use shape_mod
  implicit none

  type(rectangle) :: r 
  type(shape) :: s

  !r = rectangle(0, .true., 5, 7, 1, 2)
  !s = shape(0, .false., 12, 15)

  print *, r%width

  call r%initialize(0, .true., 5, 7, 1, 6)

  print *, r%width
end program shape_test