module class_interval
    implicit none

    type interval
        integer :: testInt
    contains
        procedure :: init
        procedure :: value => print_value
    end type interval

    ! Do I need an interface here for the default constructor?
    ! When I uncomment the below, I get the error:
    !  FUNCTION attribute conflicts with SUBROUTINE attribute in ‘interval’

    ! interface interval
    !     subroutine init(this, value1, value2)
    !         class(interval) :: this
    !         integer, intent(in) :: value1
    !         integer, intent(in) :: value2
    !     end subroutine init
    ! end interface interval

contains
    subroutine init(this, value1, value2)
        implicit none

        class(interval), intent(inout) :: this
        integer, intent(in) :: value1
        integer, intent(in) :: value2

        this%testInt = value1 + value2
    end subroutine init

    subroutine print_value(this)
        implicit none

        class(interval) :: this

        print *, this%testInt
    end subroutine print_value
end module class_interval

program test
    use class_interval

    implicit none

    type(interval) :: interval_instance

    interval_instance = interval(4, 6)
    !interval_instance = interval(10)

    call interval_instance%value()
end program test