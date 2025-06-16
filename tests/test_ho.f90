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
        call test_makestates(n_passed,n_failed)
        call test_find_mindist(n_passed, n_failed)
        call test_commute(n_passed, n_failed)
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

    subroutine test_commute(n_passed, n_failed)
        implicit none
        integer :: n_passed, n_failed
        real(r_kind), dimension(:,:), allocatable :: sigmaz, sigmap, sigmam, R, Rp, S, Sp, zero
        integer :: n, shelldegen, idx, numstates
        type(an_ho_state), dimension(:), allocatable :: states
        integer, parameter :: max_n = 2
        logical :: pass
        integer :: i, j
        character(len=10), dimension(3) :: sigma_names = ['sigmaz', 'sigmap', 'sigmam']
        character(len=10), dimension(4) :: P_names = ['R ', 'Rp', 'S ', 'Sp']
        real(r_kind), dimension(:,:,:), allocatable :: sigmas
        real(r_kind), dimension(:,:,:), allocatable :: Ps
        real(r_kind), dimension(:,:), allocatable :: A, B, C
    
        ! Generate states
        numstates = getnumstatesupto(max_n)
        allocate(states(numstates))
        idx = 1
        do n = 0, max_n
            shelldegen = getnumstates(n)
            states(idx:idx+shelldegen - 1) = get_ho_states(n)
            idx = idx + shelldegen
        end do
    
        ! Allocate matrices
        allocate(zero(numstates, numstates), sigmaz(numstates,numstates),sigmap(numstates, numstates), &
                 sigmam(numstates,numstates), R(numstates,numstates), Rp(numstates,numstates), &
                 S(numstates, numstates), Sp(numstates,numstates))
    
        zero = 0.0_r_kind
        sigmaz = pauli_z(states)
        sigmap = pauli_p(states)
        sigmam = pauli_m(states)
        Rp = Rp_mat(states)
        R = R_mat(states)
        Sp = Sp_mat(states)
        S = S_mat(states)
    
        ! Pack into arrays for looped testing
        allocate(sigmas(numstates, numstates, 3))
        allocate(Ps(numstates, numstates, 4))
    
        sigmas(:,:,1) = sigmaz
        sigmas(:,:,2) = sigmap
        sigmas(:,:,3) = sigmam
    
        Ps(:,:,1) = R
        Ps(:,:,2) = Rp
        Ps(:,:,3) = S
        Ps(:,:,4) = Sp
    
        ! Loop over all combinations and test commutators
        do i = 1, 3
            A = sigmas(:,:,i)
            do j = 1, 4
                B = Ps(:,:,j)
                C = commutator(A, B)
                pass = eq_mat(zero, C)
                if (pass) then
                    print *, 'PASS: [', trim(sigma_names(i)), ',', trim(P_names(j)), '] = 0'
                    n_passed = n_passed + 1
                else
                    print *, 'FAIL: [', trim(sigma_names(i)), ',', trim(P_names(j)), '] ≠ 0'
                    print *, 'Commutator matrix:'
                    call print_matrix(C)
                    n_failed = n_failed + 1
                end if
            end do
        end do
    
        ! Clean up
        deallocate(sigmaz, sigmap, sigmam, R, Rp, S, Sp, zero, sigmas, Ps, states)
    
    end subroutine test_commute
    


    
end module test_ho