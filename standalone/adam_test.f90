module test
    implicit none

    type element
        integer :: elementNo
    contains
        procedure :: get_elementNo
    end type element

    interface element
        procedure construct
    end interface

contains
    function get_elementNo(this) result(value)
        class(element), intent(inout) :: this
        integer :: value

        value = this%elementNo
    end function get_elementNo

    function construct(elementNo)
        type(element) :: construct
        integer, intent(in) :: elementNo
        !integer, intent(in) :: value

        construct%elementNo = elementNo + 1
    end function
end module test

program adam_test
    use test

    implicit none

    type(element) :: ele
    ele = element(5)

    print *, ele%get_elementNo()
end program adam_test