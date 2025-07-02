program run_ho
    use constants
    use def_ho, only: an_ho_state, betadef, getnumstates, getnumstatesupto, kin_en, get_ho_states
    use hamiltonian
    implicit none


    integer, parameter :: A = 240, Z=94
    type(an_ho_state), dimension(:), allocatable :: states
    type(betadef) :: def 
    integer :: n, shelldegen, idx, numstates, j
    real(r_kind), dimension(:,:), allocatable :: Vws, Tkin, V, VC,Vso,Hn,Hp
    real(r_kind), dimension(:), allocatable :: E, E_n, E_p
    real(r_kind), parameter :: V0 = V0_ws
    real(r_kind), parameter ::r0 = r0_p
    real(r_kind), parameter :: kappa=kappa_ws
    real(r_kind) :: hbaromega0, hbaromegaperp, hbaromegaz, rad,Vwsdepth

    real(r_kind), dimension(:,:), allocatable :: matr

    def = betadef(beta2 = 0, beta4=0)
    hbaromega0 = 41.0_r_kind * A**(-1.0_r_kind/3.0_r_kind) !!MeV
    hbaromegaperp = def%omega_perp(hbaromega0) !! omega = Mev/hbar
    hbaromegaz = def%omega_z(hbaromega0)
    rad = r0_p*A**(1.0/3.0)
    Vwsdepth = V0_ws*(1.0+kappa*(1.0*A-Z*2.0)/A)
    write(*,'(A,F10.3, A)')"omega:", hbaromega0, " MeV/hbar"
    write(*,'(A,F10.3, A)')"omegaz:", hbaromegaz, " MeV/hbar"
    write(*,'(A,F10.3, A)')"omegaperp", hbaromegaperp, " MeV/hbar"
    numstates = getnumstatesupto(max_n)
    allocate(states(numstates))

    
    idx = 1
    do n = 0, max_n
        shelldegen = getnumstates(n)
        states(idx:idx+shelldegen - 1) = get_ho_states(n)
        idx = idx + shelldegen
    end do
    ! write(*,*) states(1)%header()

    ! do n = 1, numstates
    !     write(*,*) states(n)%text()
    ! end do

    !Protons
    allocate(Vws(numstates,numstates), Tkin(numstates,numstates), Hn(numstates,numstates),Hp(numstates,numstates), Vc(numstates, numstates), Vso(numstates,numstates))
    !Hp = H_protons(states, Z, A, def, hbaromegaz, hbaromegaperp)
    !Hn = H_neutrons(states, Z, A, def, hbaromegaz, hbaromegaperp)
    allocate(E_p(numstates), E_n(numstates))
    call get_levels(E_p, E_n,Z,A,def, max_n)
    print*, "Neutrons:", E_n(1:5)
    print*, "Protons:" ,E_p(1:5)
    !neutrons
    Vc = coul_mat(states, def, rad, Z,mass_p,hbaromegaz, hbaromegaperp)!coul_mat(states, def, rad, Z, mass_p, hbaromegaz, hbaromegaperp)
    matr = Vc

    allocate(E(numstates), V(numstates, numstates))
    !call diagonalize(E,V,Hn)
    print *, "Protons"
    call print_levels(E_p,V)

    print*, "neutrons"
    call print_levels(E_n,V)


    open(3,file="data/out/mat.dat")
    write(3,*)
    write(3,'(14x, A5)',advance='no') "l"
    do n = 1, numstates
        write(3,'(I6)',advance='no') states(n)%ml
    end do
    write(3,*)

    write(3,'(14x, A5)',advance='no')"nz   "
    do n = 1, numstates
        write(3,'(I6)',advance='no') states(n)%nz
        !write(,*) states(n)%nr
    end do
    write(3,*)
    write(3,'(14x, A5)',advance='no')"ms   "
    do n = 1, numstates
        write(3,'(I6)',advance='no') states(n)%ms
    end do

    write(3,*)
    write(3,*)

    do n = 1, numstates
        write(3,'(3I5, 5x)', advance='no') states(n)%ml, states(n)%nz,states(n)%ms
        do j = 1, numstates
            if (abs(matr(n,j)) < 1e-8) then
                write(3,'(6x)', advance='no')
            else
                write(3,'(F6.3)', advance='no') matr(n,j)
            endif
            
        end do
        write(3,'(A)')""
    end do
    close(3)

    
end program