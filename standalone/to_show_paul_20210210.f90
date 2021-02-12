program adam_test
    print *, "Why, hello there"
    call nothing
    call something(nothing)
end program adam_test

subroutine nothing()
    print *, "Another hello"
end subroutine

subroutine something(a_procedure)
    interface
        subroutine interface_nothing()

        end subroutine
    end interface

    procedure(interface_nothing) :: a_procedure

    call a_procedure
end subroutine