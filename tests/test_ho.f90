module test_ho
    use constants
    use def_ho
    use optimise, only: conj_grad_method
    use Hamiltonian
    use test_utils
    implicit none
    private
    public :: run_tests_ho

contains

    subroutine run_tests_ho()
        implicit none
        integer :: n_passed = 0, n_failed = 0

        print *, "Running fitting tests..."
        !call test_makestates(n_passed,n_failed)
        call test_find_mindist(n_passed, n_failed)
        call test_coul(n_passed, n_failed)
        call testST(n_passed, n_failed)
        call test_deriv(n_passed,n_failed)
        call test_her(n_passed, n_failed)
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
        distfunc = dist_min(r=r, theta=theta, def=def, radius=1)

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
        distfunc = dist_min(r=r, theta=theta, def=def, radius=1)

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

    subroutine test_coul(n_passed, n_failed)
        integer :: n_passed, n_failed
        real(r_kind) :: r, theta, theta_min, dist, expected, actual
        logical ::converged, pass
        type(betadef) :: def
        type(VC_pot)::Vcoul
        integer :: Z = 50
        
        def = betadef(beta2 = 0.0, beta4=0.0)
        vcoul%def = def
        vcoul%radius = 1.0_r_kind
        call vcoul%set_charge_dens(Z)


        actual = vcoul%eval(0.0_r_kind, 0.0_r_kind)
        expected = el_pot(0.0_r_kind, vcoul%radius, vcoul%charge_dens)
       
        pass = eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: Vc r = 0", expected, "but got", actual
        end if
        actual = vcoul%eval(1.0_r_kind, 0.0_r_kind)
        expected = el_pot(1.0_r_kind, vcoul%radius, vcoul%charge_dens)
        
        pass = eq_r(expected, actual, 1e-1_r_kind)  !
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: Vc r = R", expected, "but got", actual
        end if

        actual = vcoul%eval(3.0_r_kind, 0.0_r_kind)
        expected = el_pot(3.0_r_kind, vcoul%radius, vcoul%charge_dens)
        pass = eq_r(expected, actual, 1e-4_r_kind)  !
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: Vc r = 3R", expected, "but got", actual
        end if

    end subroutine test_coul

    subroutine testST(n_passed, n_failed)
        integer :: n_passed, n_failed
        real(r_kind) :: r, theta, theta_min, dist, expected, actual, val, eta, val2
        logical ::converged, pass
        type(an_ho_state) :: s1, s2
        type(WS_pot) :: pot
        type(betadef) :: def
        integer :: Z = 50
        def = betadef(beta2 = 0, beta4=0)
        s1%nr = 1
        s2%nr = 1
        s1%ml = 2
        s2%ml = 3
        s1%ms = 1
        s2%ms = -1
        s1%nz = 2
        s2%nz = 1
        actual = S0(1.5_r_kind,s1,s2)
        expected = S0(1.5_r_kind,s2,s1)
       
        pass = eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: testst", expected, "but got", actual
        end if

        actual = Sp(1.5_r_kind,s1,s2)

        expected = Sm(1.5_r_kind,s2,s1)
        pass = eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: Splus Sminus symmetry", expected, "but got", actual
        end if
        pot%radius = 5
        pot%V_0= 1
        pot%def = def

        eta = 1.0_r_kind
        expected = Tplus(eta,1,s1,s2,pot,1.0_r_kind,1.0_r_kind)
        actual = Tminus(eta,1,s2,s1,pot,1.0_r_kind,1.0_r_kind)
        pass = eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: testst TpTm", expected, "but got", actual
        end if

        expected = VSO_off_diag_elem_v2(s1, s2, pot, 1.0_r_kind, 1.0_r_kind)
        actual =VSO_off_diag_elem_v2(s2, s1, pot, 1.0_r_kind, 1.0_r_kind)
        val = VSO_off_diag_elem_v2(s2,s1,pot, 1.0_r_kind, 1.0_r_kind)
        print*, expected
        print*, actual
        print*, val
        pass = eq_r(expected, actual)  ! Replace with actual call: e.g. call micmac_result(actual)
        print*, "testing...", val, val2
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: Test SO off diag", expected, "but got", actual
        end if
    end subroutine



    subroutine test_deriv(n_passed, n_failed)
        real(r_kind) :: expected, actual, x
        real(r_kind), parameter :: dx = 1e-6_r_kind
        integer :: n_passed, n_failed, n, l
        logical :: pass
        x = 0.5_r_kind
        expected = exp((x*x)/2.0) * (exp(-(x+dx)*(x+dx)/2.0) *hmn(x+dx,1) - exp(-((x*x)/2.0)) *hmn(x,1)) / dx
        actual = hmnp(x,1)
        pass = eq_r(expected, actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test deriv", expected, "but got", actual
        end if
    end subroutine

    subroutine test_her(n_passed, n_failed)
        real(r_kind) :: expected, actual, x
        real(r_kind), parameter :: dx = 1e-6_r_kind
        integer :: n, n_passed, n_failed, l
        logical :: pass
        x = 3.0_r_kind
        n = 4
        expected = sqrt(real(2*n))*Hmn(x,n-1)
        actual =( Hmn(x+dx,n)-Hmn(x,n))/dx
        pass = eq_r(expected, actual,1e-2_r_kind)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test hermite", expected, "but got", actual
        end if
        l = 3
        expected = ((2.0*n + l)/2.0* gnl(x,n,l) - sqrt(real(n*(n+l)))* gnl(x,n-1,l))/x
        actual =( gnl(x+dx,n,l)-gnl(x,n,l))/dx
        pass = eq_r(expected, actual)
        if (pass) then
            n_passed = n_passed + 1
        else
            n_failed = n_failed + 1
            print *, "FAIL: test mod lag", expected, "but got", actual
        end if
    end subroutine


    
end module test_ho