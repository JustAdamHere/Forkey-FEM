module class_solution
    type solution
        integer :: number
    contains
        procedure :: set_number
    end type
contains
    subroutine set_number(this, a_number)
        class(solution) :: this
        integer         :: a_number

        this%number = a_number
    end subroutine
end module

program test
    use class_solution

    type(solution), target  :: mySolution
    type(solution), pointer :: myNewSolution

    call mySolution%set_number(5)
    print *, mySolution%number

    myNewSolution => mySolution

    print *, myNewSolution%number
end program