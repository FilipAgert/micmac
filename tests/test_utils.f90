module test_utils
    use constants
    implicit none
    private
    public :: eq, eq_r

    real(r_kind), parameter :: epsilon = 1e-9




    contains

    logical function eq(n1,n2)
        integer, intent(in) :: n1,n2
        eq = n1 == n2
    end function

    logical function eq_r(r1,r2)
        real(r_kind), intent(in) :: r1, r2
        if(abs(r1-r2) < epsilon) then
            eq_r = .true.
        else
            eq_r = .false.
        endif
    end function


end module test_utils