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
        call test_find_gs(n_passed, n_failed)
        call test_F(n_passed,n_failed)
        print *, "micmac results: ", n_passed, " passed,", n_failed, " failed"
    end subroutine run_tests_micmac

    subroutine test_staircase_low(n_passed, n_failed)
        integer, intent(inout) :: n_passed, n_failed
        real(kind) :: expected, actual
        logical :: pass

        ! Dummy test: expected == actual
        expected = 3.0_kind/5.0_kind * (2.0_kind**(5.0_kind/3.0_kind)) / (2.0_kind)
        actual = staircase(2.0_kind, magic_num_N)   ! Replace with actual call: e.g. call micmac_result(actual)
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
        real(kind) :: expected, actual
        logical :: pass

        ! Dummy test: expected == actual
        expected = 3.0_kind/5.0_kind * (2.0_kind**(5.0_kind/3.0_kind)) / (2.0_kind)
        actual = staircase(1.0_kind, magic_num_N)   ! Replace with actual call: e.g. call micmac_result(actual)
        pass = eq_r(expected,actual)

        pass = pass .and. eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_example - expected", expected, "but got", actual
        end if
    end subroutine test_staircase_mid

    subroutine test_find_gs(n_passed, n_failed)
        integer, intent(inout) :: n_passed, n_failed
        real(kind) :: expected, actual
        logical :: pass, found
        integer :: Z, A
        real(kind) :: BE, BE2, BE3, BE4
        
        type(deformation):: def 
        real(kind), dimension(num_params) :: params
        params = starting_params
        Z = 92
        A = 238
        call find_gs(BE, def, found,params, Z, A)   
        ! BE2 = binding_energy(params, Z, A)
        ! BE3 = binding_energy_def(params, Z, A, 0.0_kind)
        ! BE4 = binding_energy_def(params, Z, A, def)
        ! write(*,*) "BE fit with def: ", BE
        ! WRITE(*,*) "BE def func with def = 0: ", BE3
        ! WRITE(*,*) "Be def func with def = ", def, ", :", BE4
        ! write(*,*) "BE with no def: ", BE2
    end subroutine

    subroutine test_F(n_passed,n_failed)
        integer, intent(inout) :: n_passed, n_failed
        real(kind) :: expected, actual
        logical :: pass
        integer :: Z, A
        real(kind) :: BE,  BE2, BE3, BE4
        type(deformation):: def 
        real(kind), dimension(num_params) :: params
        expected = 0.0_kind
        actual = F(82,.true.)
        pass = eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_example - expected", expected, "but got", actual
        end if

    end subroutine
    
end module test_micmac
