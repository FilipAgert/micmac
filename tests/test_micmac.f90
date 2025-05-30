module test_micmac
    use micmac
    use constants
    use test_utils
    implicit none
    private
    public :: run_tests_micmac

contains

    subroutine run_tests_micmac()
        implicit none
        integer :: n_passed = 0, n_failed = 0

        print *, "Running micmac tests..."
        call test_staircase_low(n_passed, n_failed)
        call test_staircase_mid(n_passed, n_failed)
        print *, "micmac results: ", n_passed, " passed,", n_failed, " failed"
    end subroutine run_tests_micmac

    subroutine test_staircase_low(n_passed, n_failed)
        integer, intent(inout) :: n_passed, n_failed
        real(r_kind) :: expected, actual
        logical :: pass

        ! Dummy test: expected == actual
        expected = 3.0_r_kind/5.0_r_kind * (2.0_r_kind**(5.0_r_kind/3.0_r_kind)) / (2.0_r_kind)
        actual = staircase(2.0_r_kind, magic_num_N)   ! Replace with actual call: e.g. call micmac_result(actual)
        pass = eq_r(expected,actual)

        pass = pass .and. eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_example - expected", expected, "but got", actual
        end if
    end subroutine test_staircase_low

    subroutine test_staircase_mid(n_passed, n_failed)
        integer, intent(inout) :: n_passed, n_failed
        real(r_kind) :: expected, actual
        logical :: pass

        ! Dummy test: expected == actual
        expected = 3.0_r_kind/5.0_r_kind * (2.0_r_kind**(5.0_r_kind/3.0_r_kind)) / (2.0_r_kind)
        actual = staircase(1.0_r_kind, magic_num_N)   ! Replace with actual call: e.g. call micmac_result(actual)
        pass = eq_r(expected,actual)

        pass = pass .and. eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_example - expected", expected, "but got", actual
        end if
    end subroutine test_staircase_mid

    
end module test_micmac
