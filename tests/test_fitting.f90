module test_fitting
    use constants
    use test_utils
    implicit none
    private
    public :: run_tests_fitting

contains

    subroutine run_tests_fitting()
        implicit none
        integer :: n_passed = 0, n_failed = 0

        print *, "Running fitting tests..."
        print *, "fitting results: ", n_passed, " passed,", n_failed, " failed"
    end subroutine run_tests_fitting
    
end module test_fitting