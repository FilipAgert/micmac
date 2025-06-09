program run_ho
    use constants
    use def_ho, only: getnumstatesupto, getnumstates, get_ho_states, an_ho_state, Ws_mat_elem, betadef
    use hamiltonian, only: diagonalize
    implicit none


    integer, parameter :: max_n = 5
    integer, parameter :: A = 48, Z=20
    type(an_ho_state), dimension(:), allocatable :: states
    type(an_ho_state) :: state1, state2
    type(betadef) :: def 
    integer :: n, shelldegen, idx, m, numstates, ii
    real(r_kind), dimension(:,:), allocatable :: Vws, Tkin, H, V
    real(r_kind), dimension(:), allocatable :: E
    real(r_kind), parameter :: V0 = 49.6
    real(r_kind), parameter ::r0 = 1.347
    real(r_kind), parameter :: kappa=0.86
    real(r_kind) :: hbaromega0, radius, hbaromegaperp, hbaromegaz
    real(r_kind) :: I, Vwsdepth

    def = betadef(beta2 = 0.4, beta4=0.0)
    hbaromega0 = 41.0_r_kind * A**(-1.0_r_kind/3.0_r_kind) !!MeV
    hbaromegaperp = def%omega_perp(hbaromega0) !! omega = Mev/hbar
    hbaromegaz = def%omega_z(hbaromega0)
    radius = r0 * A**(1.0_r_kind/3.0_r_kind)

    I = (A-2.0_r_kind*Z)/A
    Vwsdepth = V0 * (1.0_r_kind-kappa*I)
    write(*,'(A,F10.3)') "Well depth:", Vwsdepth
    write(*,'(A,F10.3, A)')"omega:", hbaromega0, " MeV/hbar"
    write(*,'(A,F10.3, A)')"omegaz:", hbaromegaz, " MeV/hbar"
    write(*,'(A,F10.3, A)')"omegaperp", hbaromegaperp, " MeV/hbar"
    numstates = getnumstatesupto(max_n)
    allocate(states(numstates))

    idx = 1
    do n = 0, max_n
        shelldegen = getnumstates(n)
        write(*,*) "N: ", N
        write(*,*) "Degen: ", shelldegen
        states(idx:idx+shelldegen - 1) = get_ho_states(n)
        idx = idx + shelldegen
    end do

    write(*,'(A)') states(1)%header()
    do n = 1, numstates

        write(*,*) states(n)%text()
    end do

    allocate(Vws(numstates,numstates), Tkin(numstates,numstates), H(numstates,numstates))
    

    do n = 1, numstates
        state1 = states(n)
        do m = 1, n
            state2 = states(m)
            Vws(n,m) = Ws_mat_elem(state1, state2, def,VwsDepth,radius,mass_n,hbaromegaz,hbaromegaperp)

        end do
        write(*,'(I5,A,I5,A)') n, " out of ", numstates, " rows completed"
    end do

    do n = 1, numstates !!Use the fact that matrix elements are symmetric.
        state1 = states(n)
        do m = numstates, n + 1, -1
            state2 = states(m)
            Vws(n,m) = Vws(m,n)

        end do
        ! write(*,'(40F7.3)') Vws(n,:)
    end do

    Tkin = 0.0_r_kind
    write(*,*)
    write(*,*) "Kinetic energy"
    do n = 1, numstates
        Tkin(n,n) = states(n)%kinetic_energy(hbaromegaz,hbaromegaperp)
        ! write(*,'(15F10.3)') Tkin(n,:)
    end do
    write(*,*)
    ! write(*,*) "Hamiltonian"
    H = Tkin + Vws


    allocate(V(numstates,numstates), E(numstates))
    call diagonalize(E, V, H)
    write(*,'(A)') "n     E_n (MeV)  E_ho (MeV)  |  000    001-    001+     100    002-    002+    010    101-    101+    200"
    do n = 1, min(numstates,10)
        write(*,'(I3, F10.1, F10.1)',advance='no')n, E(n) - E(1), states(n)%kinetic_energy(hbaromegaz,hbaromegaperp) * 2 - states(1)%kinetic_energy(hbaromegaz,hbaromegaperp) * 2

        write(*,'(A, 20F8.3)') "     ", V(n,1:20)
    end do

end program