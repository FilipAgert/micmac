module test_ho
    use constants
    use def_ho
    use optimise, only: conj_grad_method
    use test_utils
    implicit none
    private
    public :: run_tests_ho

contains

    subroutine run_tests_ho()
        implicit none
        integer :: n_passed = 0, n_failed = 0

        print *, "Running fitting tests..."
        call test_makestates(n_passed,n_failed)
        call test_find_mindist(n_passed, n_failed)
        print *, "fitting results: ", n_passed, " passed,", n_failed, " failed"
    end subroutine run_tests_ho

    subroutine test_makestates(n_passed, n_failed)
        integer :: n_passed, n_failed, i
        integer, parameter :: N = 3
        type(an_ho_state), dimension(getnumstates(N)) :: states
        states = get_ho_states(N)

        write(*,'(A)') "nz    nr    ml    2omega"
        do i = 1, getnumstates(N)
            write(*,'(4I4)') states(i)
        end do


    end subroutine

    subroutine test_find_mindist(n_passed, n_failed)
        integer :: n_passed, n_failed, i
        real(r_kind) :: r, theta, theta_min, dist, expected, actual
        logical ::converged, pass
        type(dist_min) :: distfunc
        type(betadef) :: def
        def = betadef(beta2 = 0.2, beta4=0.0)
        r = 0.5
        theta = pi/2.0_r_kind
        distfunc = dist_min(r=r, theta=theta, def=def)

        call conj_grad_method(theta_min,dist,converged,distfunc,0.0_r_kind,pi,theta,1e-5_r_kind)

        expected = pi/2
        actual = theta_min
        pass = eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_example - expected", expected, "but got", actual
        end if
        actual = sqrt(dist)
        expected = 0.4369216869
        pass = eq_r(expected, actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_example - expected", expected, "but got", actual
        end if


        theta = pi-0.1
        distfunc = dist_min(r=r, theta=theta, def=def)

        call conj_grad_method(theta_min,dist,converged,distfunc,0.0_r_kind,pi,theta,1e-5_r_kind)
        expected = 2.9714242267_r_kind
        actual = theta_min
        pass = eq_r(expected, actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_example - expected", expected, "but got", actual
        end if
        actual = sqrt(dist)
        expected = 0.6229471088
        pass = eq_r(expected, actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test_example - expected", expected, "but got", actual
        end if

    end subroutine
    
end module test_ho