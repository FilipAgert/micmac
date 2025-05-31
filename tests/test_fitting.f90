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
        call test_solve_linsys_2x2(n_passed, n_failed)
        call test_solve_4x4(n_passed, n_failed)
        print *, "fitting results: ", n_passed, " passed,", n_failed, " failed"
    end subroutine run_tests_fitting

    subroutine test_solve_linsys_2x2(n_passed, n_failed)
        use fitting, only: solve_linsys
        implicit none
        integer, intent(inout) :: n_passed, n_failed
        real(r_kind), dimension(2,2) :: A
        real(r_kind), dimension(2)   :: b, x, expected
        logical :: pass

        ! Set up system Ax = b
        A = reshape([2.0_r_kind, 1.0_r_kind, &
                    5.0_r_kind, 7.0_r_kind], shape(A), order=[2,1])
        b = [11.0_r_kind, 13.0_r_kind]
        expected = [7.1111111111_r_kind, -3.2222222222_r_kind]

        ! Call solver
        call solve_linsys(A, b, 2)

        ! Check results
        pass = eq_r(b(1), expected(1)) .and. eq_r(b(2), expected(2))

        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_solve_linsys"
            print *, "Expected:", expected
            print *, "Got     :", b
        end if
    end subroutine test_solve_linsys_2x2

    subroutine test_solve_4x4(n_passed, n_failed)
        use fitting, only: solve_linsys
        implicit none
        integer, intent(inout) :: n_passed, n_failed
        real(r_kind), dimension(4,4) :: A
        real(r_kind), dimension(4)   :: b, expected
        integer :: i
        logical :: pass

        ! Initialize matrix A and vector b
        A = reshape([ &
            4.0_r_kind, -1.0_r_kind,  0.0_r_kind,  0.0_r_kind, &
        -1.0_r_kind,  4.0_r_kind, -1.0_r_kind,  0.0_r_kind, &
            0.0_r_kind, -1.0_r_kind,  4.0_r_kind, -1.0_r_kind, &
            0.0_r_kind,  0.0_r_kind, -1.0_r_kind,  3.0_r_kind  &
        ], shape(A), order=[2,1])

        b = [15.0_r_kind, 10.0_r_kind, 10.0_r_kind, 10.0_r_kind]
        expected = [5.0_r_kind, 5.0_r_kind, 5.0_r_kind, 5.0_r_kind]

        ! Solve the system
        call solve_linsys(A, b, 4)

        ! Validate solution
        pass = all([(eq_r(b(i), expected(i)), i = 1, 4)])

        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_solve_4x4"
            print *, "Expected:", expected
            print *, "Got     :", b
        end if
    end subroutine test_solve_4x4



    
end module test_fitting