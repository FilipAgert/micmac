program run_ho
    use constants
    use def_ho, only: an_ho_state, betadef, getnumstates, getnumstatesupto, kin_en, get_ho_states
    use hamiltonian
    implicit none


    integer, parameter :: max_n = 8
    integer, parameter :: A = 208, Z=82
    type(an_ho_state), dimension(:), allocatable :: states
    type(betadef) :: def 
    integer :: n, shelldegen, idx, numstates
    real(r_kind), dimension(:,:), allocatable :: Vws, Tkin, V, VC,Vso,Hn,Hp
    real(r_kind), dimension(:), allocatable :: E
    real(r_kind), parameter :: V0 = V0_ws
    real(r_kind), parameter ::r0 = r0_p
    real(r_kind), parameter :: kappa=kappa_ws
    real(r_kind) :: hbaromega0, radius, hbaromegaperp, hbaromegaz, radius_so
    real(r_kind) :: I, Vwsdepth
    real(r_kind), dimension(:,:), allocatable :: matr

    def = betadef(beta2 = 0, beta4=0)
    hbaromega0 = 41.0_r_kind * A**(-1.0_r_kind/3.0_r_kind) !!MeV
    hbaromegaperp = def%omega_perp(hbaromega0) !! omega = Mev/hbar
    hbaromegaz = def%omega_z(hbaromega0)
    radius = r0 * A**(1.0_r_kind/3.0_r_kind)
    radius_so = r0_so_n* A**(1.0_r_kind/3.0_r_kind)

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

    !Protons
    allocate(Vws(numstates,numstates), Tkin(numstates,numstates), Hn(numstates,numstates),Hp(numstates,numstates), Vc(numstates, numstates), Vso(numstates,numstates))
    Hp = H_protons(states, Z, A, def, hbaromegaz, hbaromegaperp)
    Hn = H_neutrons(states, Z, A, def, hbaromegaz, hbaromegaperp)

    !neutrons
    

    matr = Hp

    allocate(E(numstates), V(numstates, numstates))
    call diagonalize(E,V,Hp)
    call print_levels(E,V)


    open(3,file="data/out/mat.dat")
    write(3,'(19x, A5)',advance='no') "r    "
    do n = 1, numstates
        write(3,'(I5)',advance='no') states(n)%r
    end do
    write(3,*)
    write(3,'(19x, A5)',advance='no') "s   "
    do n = 1, numstates
        write(3,'(I5)',advance='no') states(n)%s
    end do
    write(3,*)

    write(3,'(19x, A5)',advance='no')"nz   "
    do n = 1, numstates
        write(3,'(I5)',advance='no') states(n)%nz
        !write(,*) states(n)%nr
    end do
    write(3,*)
    write(3,'(19x, A5)',advance='no')"ms   "
    do n = 1, numstates
        write(3,'(I5)',advance='no') states(n)%ms
    end do

    write(3,*)
    write(3,*)

    do n = 1, numstates
        write(3,'(4I5, 5x, 3000F5.1)') states(n)%r, states(n)%s,states(n)%nz,  states(n)%ms,matr(n,:)
    end do
    close(3)

    contains

    function H_protons(states, Z,A,def,hbaromegaz, hbaromegaperp) result(H)
        implicit none
        integer :: Z, A
        type(betadef) :: def
        real(r_kind) :: hbaromegaz, hbaromegaperp, Vwsdepth, I
        type(an_ho_state) :: states(:)
        real(r_kind), dimension(size(states),size(states)) :: H, VSO, VC, VWS, Tkin
        radius = r0_p * A**(1.0_r_kind/3.0_r_kind)
        radius_so = r0_so_p* A**(1.0_r_kind/3.0_r_kind)
        I = (A-2.0_r_kind*Z)/A
        Vwsdepth = V0 * (1.0_r_kind+kappa_ws*I)
        Vc = coul_mat(states, def, radius, Z, mass_p, hbaromegaz, hbaromegaperp)
        Vws = Vws_mat(states,def,radius, hbaromegaz,hbaromegaperp,mass_p,Vwsdepth)
        Vso = Vso_mat(states, def, radius_so, hbaromegaz,hbaromegaperp, mass_p, Vwsdepth, lambda_p)
        Tkin = kin_en(states, hbaromegaz, hbaromegaperp)
        H = Vc + Vws + Vso + Tkin

    end function

    function H_neutrons(states, Z,A,def,hbaromegaz, hbaromegaperp) result(H)
        implicit none
        integer :: Z, A
        type(betadef) :: def
        real(r_kind) :: hbaromegaz, hbaromegaperp, Vwsdepth, I
        type(an_ho_state) :: states(:)
        real(r_kind), dimension(size(states),size(states)) :: H, VSO, VC, VWS, Tkin
        radius = r0_n * A**(1.0_r_kind/3.0_r_kind)
        radius_so = r0_so_n* A**(1.0_r_kind/3.0_r_kind)
        I = (A-2.0_r_kind*Z)/A
        Vwsdepth = V0 * (1.0_r_kind-kappa_ws*I)
        Vws = Vws_mat(states,def,radius, hbaromegaz,hbaromegaperp,mass_p,Vwsdepth)
        Vso = Vso_mat(states, def, radius_so, hbaromegaz,hbaromegaperp, mass_p, Vwsdepth, lambda_n)
        Tkin = kin_en(states, hbaromegaz, hbaromegaperp)
        H = Vws + Vso + Tkin

    end function
end program