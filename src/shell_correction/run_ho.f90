program run_ho
    use constants
    use def_ho, only: getnumstatesupto, getnumstates, get_ho_states, an_ho_state, betadef
    use hamiltonian, only: diagonalize, mat_elem_axsym, WS_pot, VC_pot
    implicit none


    integer, parameter :: max_n = 3
    integer, parameter :: A = 238, Z=92
    type(an_ho_state), dimension(:), allocatable :: states
    type(an_ho_state) :: state1, state2
    type(betadef) :: def 
    type(WS_pot) :: WSpot
    type(VC_pot) :: VCpot
    integer :: n, shelldegen, idx, m, numstates, ii
    real(r_kind), dimension(:,:), allocatable :: Vws, Tkin, H, V, VC
    real(r_kind), dimension(:), allocatable :: E
    real(r_kind), parameter :: V0 = 49.6
    real(r_kind), parameter ::r0 = 1.275
    real(r_kind), parameter :: kappa=0.86
    real(r_kind) :: hbaromega0, radius, hbaromegaperp, hbaromegaz
    real(r_kind) :: I, Vwsdepth, r, theta, Vcone, vcold

    def = betadef(beta2 = 0.2, beta4=0.1)
    hbaromega0 = 41.0_r_kind * A**(-1.0_r_kind/3.0_r_kind) !!MeV
    hbaromegaperp = def%omega_perp(hbaromega0) !! omega = Mev/hbar
    hbaromegaz = def%omega_z(hbaromega0)
    radius = r0 * A**(1.0_r_kind/3.0_r_kind)

    I = (A-2.0_r_kind*Z)/A
    Vwsdepth = V0 * (1.0_r_kind+kappa*I)
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
    write(*,*) states(1)%header()

    do n = 1, numstates
        write(*,*) states(n)%text()
    end do

    allocate(Vws(numstates,numstates), Tkin(numstates,numstates), H(numstates,numstates), Vc(numstates, numstates))
    
    WSpot%def = def
    WSpot%radius = radius
    Wspot%V_0 = Vwsdepth

    VCpot%def = def
    VCpot%radius = radius
    call VCpot%set_charge_dens(Z)

    do n = 1, numstates
        state1 = states(n)
        do m = 1, n
            state2 = states(m)
            Vws(n,m) = mat_elem_axsym(state1, state2, Wspot,mass_p,hbaromegaz,hbaromegaperp)
            Vc(n,m) = mat_elem_axsym(state1, state2, VCpot, mass_p, hbaromegaz, hbaromegaperp)
        end do
        ! write(*,'(I5,A,I5,A)') n, " out of ", numstates, " rows completed"
    end do

    do n = 1, numstates !!Use the fact that matrix elements are symmetric.
        state1 = states(n)
        do m = numstates, n + 1, -1
            state2 = states(m)
            Vws(n,m) = Vws(m,n)
            Vc(n,m) = Vc(m,n)

        end do
        ! write(*,'(40F7.3)') Vws(n,:)
    end do
    write(*,*) "Coul"
    do n = 1,numstates
        ! write(*,'(40F7.3)') Vc(n,:)
    end do

    Tkin = 0.0_r_kind
    write(*,*)
    write(*,*) "Kinetic energy"
    do n = 1, numstates
        Tkin(n,n) = states(n)%kinetic_energy(hbaromegaz,hbaromegaperp)
        ! write(*,'(15F10.3)', advance='no') Tkin(n,n)
    end do
    write(*,*)
    ! write(*,*) "Hamiltonian"
    open(3,file="data/out/H.dat")
    open(4,file="data/out/WS.dat")
    open(5,file="data/out/VC.dat")
    H = Tkin + Vws + Vc
    do n = 1,numstates
        do idx = 1,numstates
            write(3, '(F10.3)',advance='no') H(n,idx)
            write(4, '(F10.3)',advance='no') Vws(n,idx)
            write(5, '(F10.3)',advance='no') Vc(n,idx)
        end do
        write(3,*)
        write(4,*)
        write(5,*)
    end do
    close(3)
    close(4)
    close(5)

    allocate(V(numstates,numstates), E(numstates))
    call diagonalize(E, V, H)
    write(*,'(A)') "n     E_n (MeV)  E_ho (MeV)  |  000    001-    001+     100    002-    002+    010    101-    101+    200"
    do n = 1, min(numstates,50)
        write(*,'(I3, F10.1, F10.1)',advance='no')n, E(n), states(n)%kinetic_energy(hbaromegaz,hbaromegaperp) * 2 

        write(*,'(A, 10F8.3)') "     ", V(n,1:min(numstates,10))
    end do



    r = 0.0
    theta = 0.0
    vcold = Vcpot%eval(r,theta)
    write(*,'(A,e15.5,A)') "Pot: ", vcold, " at, 0 fm"

    r = 0.5
    theta = 0.0
    Vcone = Vcpot%eval(r,theta)
    write(*,'(A,e15.5,A)') "Pot: ", Vcone, " at, 0.5 fm"
    write(*,'(A, f10.3)') "ratio: ", vcold/Vcone



end program