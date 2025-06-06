module test_optimise
    use optimise
    use constants, only: r_kind, default_num_restarts
    use test_utils

    implicit none
    type, extends(func_1d) :: poly2
    contains
        procedure :: eval =>  poly2d
    end type

    type, extends(func_nd) :: redh
    contains
        procedure :: eval =>  red_herring_poly
    end type



    private
    public :: run_tests_optimise
    contains 
    subroutine run_tests_optimise()
        implicit none
        integer :: n_passed = 0, n_failed = 0

        print *, "Running minimiser tests..."
        call test_1d_minval(n_passed, n_failed)
        call test_1d_minloc(n_passed, n_failed)
        call test_1d_findglobal(n_passed, n_failed)
        print *, "mimimiser results: ", n_passed, " passed,", n_failed, " failed"
    end subroutine run_tests_optimise


    subroutine test_1d_minloc(n_passed, n_failed)
        integer, intent(inout) :: n_passed, n_failed
        real(r_kind) :: expected, actual, fmin
        logical :: pass, conv
        type(poly2) :: func

        ! Dummy test: expected == actual
        expected = 5.0_r_kind/6.0_r_kind
        call find_optimal_point(actual, fmin,conv, func,-100.0_r_kind,100.0_r_kind,10.0_r_kind)  ! Replace with actual call: e.g. call micmac_result(actual)
        pass = eq_r(expected,actual)
        !print*, "niter = ", niter
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_example - expected", expected, "but got", actual
        end if

         
    end subroutine

    subroutine test_1d_minval(n_passed,n_failed)
        integer, intent(inout) :: n_passed, n_failed
        real(r_kind) :: expected, actual, fmin
        logical :: pass
        logical :: conv
        type(poly2) :: func
        call find_optimal_point(actual, fmin,conv, func,-100.0_r_kind,100.0_r_kind,10.0_r_kind)  ! Replace with actual call: e.g. call micmac_result(actual)
        !print*, "niter = ", niter
        expected = -1.0_r_kind/12.0_r_kind
        actual = fmin
        pass = eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_example - expected", expected, "but got", actual
        end if

    end subroutine

    subroutine test_1d_findglobal(n_passed,n_failed)
        integer, intent(inout) :: n_passed, n_failed
        real(r_kind) :: expected, actual, fmin,xmin(1), lb(1), ub(1)
        logical :: pass
        type(redh) :: func
        lb = -10
        ub = 10
        call find_min(xmin, fmin, func,lb, ub,10, 1)  ! Replace with actual call: e.g. call micmac_result(actual)
        !print*, "niter = ", niter
        expected = -2.63891740462711_r_kind
        actual = xmin(1)
        pass = eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_findglobal - xvalue expected:", expected, "but got", actual
        end if

        expected = -14.3327254685187_r_kind
        actual = fmin
        pass = eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_findglobal - yvalue expected:", expected, "but got", actual
        end if

    end subroutine

    pure elemental function poly2d(self, x) result(f)
        ! A simple quadratic function for testing
        class(poly2), intent(in) :: self
        real(r_kind), intent(in) :: x
        real(r_kind) :: f

        f = 3.0_r_kind * x**2 - 5.0_r_kind * x + 2.0_r_kind

    end function

    pure function red_herring_poly(self, xs) result(f)
        ! A simple quadratic function for testing
        real(r_kind), intent(in) :: xs(:)
        class(redh), intent(in) :: self
        real(r_kind) :: f

        f = xs(1)**4 + 3*xs(1)**3 - 3*xs(1)**2 -5*xs(1)

    end function

end module test_optimise