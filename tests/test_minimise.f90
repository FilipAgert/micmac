module test_minimise
    use minimise
    use constants, only: r_kind, default_num_restarts
    use test_utils

    implicit none
    private
    public :: run_tests_minimise
    contains 
    subroutine run_tests_minimise()
        implicit none
        integer :: n_passed = 0, n_failed = 0

        print *, "Running minimiser tests..."
        call test_1d_minval(n_passed, n_failed)
        call test_1d_minloc(n_passed, n_failed)
        call test_1d_findglobal(n_passed, n_failed)
        print *, "mimimiser results: ", n_passed, " passed,", n_failed, " failed"
    end subroutine run_tests_minimise


    subroutine test_1d_minloc(n_passed, n_failed)
        integer, intent(inout) :: n_passed, n_failed
        real(r_kind) :: expected, actual, fmin
        logical :: pass, conv


        ! Dummy test: expected == actual
        expected = 5.0_r_kind/6.0_r_kind
        call find_local_min(actual, fmin,conv, poly2d,-100.0_r_kind,100.0_r_kind,10.0_r_kind)  ! Replace with actual call: e.g. call micmac_result(actual)
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
        call find_local_min(actual, fmin,conv, poly2d,-100.0_r_kind,100.0_r_kind,10.0_r_kind)  ! Replace with actual call: e.g. call micmac_result(actual)
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
        real(r_kind) :: expected, actual, fmin,xmin
        logical :: pass
        call find_min(xmin, fmin, red_herring_poly,-10.0_r_kind,10.0_r_kind,10)  ! Replace with actual call: e.g. call micmac_result(actual)
        !print*, "niter = ", niter
        expected = -2.63891740462711_r_kind
        actual = xmin
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

    function poly2d(x) result(f)
        ! A simple quadratic function for testing
        real(r_kind), intent(in) :: x
        real(r_kind) :: f

        f = 3.0_r_kind * x**2 - 5.0_r_kind * x + 2.0_r_kind

    end function

    function red_herring_poly(x) result(f)
                ! A simple quadratic function for testing
        real(r_kind), intent(in) :: x
        real(r_kind) :: f

        f = x**4 + 3*x**3 - 3*x**2 -5*x

    end function

end module test_minimise