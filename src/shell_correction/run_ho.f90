program run_ho
    use constants
    use def_ho, only: an_ho_state, betadef, getnumstates, getnumstatesupto, kin_en, get_ho_states, get_ho_states_upto
    use hamiltonian
    implicit none


    integer, parameter :: A = 208, Z=82
    type(an_ho_state), dimension(:), allocatable :: states
    type(betadef) :: def 
    type(VC_pot) :: VC
    integer :: n, shelldegen, idx, numstates, j
    real(kind), dimension(:,:), allocatable :: Vws, Tkin, V, VC_MAT,Vso,Hn,Hp
    real(kind), dimension(:), allocatable :: E, E_n, E_p
    real(kind), parameter :: V0 = V0_ws
    real(kind), parameter ::r0 = r0_p
    real(kind), parameter :: kappa=kappa_ws
    real(kind) :: hbaromega0, hbaromegaperp, hbaromegaz, rad,Vwsdepth,matelem
    character :: p
    real(kind), dimension(:,:), allocatable :: matr
    real(kind) :: r, theta, vcpot, phi
    real(kind) :: normal_vec(3)

    def = betadef(beta2 = 0.0, beta4=0.0)
    hbaromega0 = 41.0_kind * A**(-1.0_kind/3.0_kind) !!MeV
    hbaromegaperp = def%omega_perp(hbaromega0) !! omega = Mev/hbar
    hbaromegaz = def%omega_z(hbaromega0)
    rad = r0_p*A**(1.0/3.0)
    Vwsdepth = V0_ws*(1.0+kappa*(1.0*A-Z*2.0)/A)

    

    ! write(*,*) states(1)%header()

    ! do n = 1, numstates
    !     write(*,*) states(n)%text()
    ! end do

    !Protons
    !Hp = H_protons(states, Z, A, def, hbaromegaz, hbaromegaperp)
    !Hn = H_neutrons(states, Z, A, def, hbaromegaz, hbaromegaperp)
    VC%radius = rad
    call VC%set_charge_dens(Z)
    VC%def = def

    r = 10.0
    theta = pi/2.0
    phi = 0
    normal_vec = normal_sph(def, theta, phi)
    write(*,'(a,f7.2)') "r = ", r
    write(*,'(a,f7.2)') "Theta = ", theta
    vcpot = VC%eval(r, theta)
    write(*,'(a,f7.2)') "Pot = ", vcpot
    write(*,'(a,3f7.2)')"normal: ", normal_vec

    call exit
    allocate(E_p(num_p_states), E_n(num_n_states))
    call print_shell_params(Z,A,def)
    write(*,'(A,F10.3, A)')"omega:", hbaromega0, " MeV/hbar"
    write(*,'(A,F10.3, A)')"omegaz:", hbaromegaz, " MeV/hbar"
    write(*,'(A,F10.3, A)')"omegaperp", hbaromegaperp, " MeV/hbar"
    call get_levels(E_p, E_n,Z,A,def, max_n)

    ! print*, "Neutrons:", E_n(1:5)
    ! print*, "Protons:" ,E_p(1:5)
    !neutrons

    open(4, file="data/out/levels.dat")
    write(4,'(I3,A,I3,A)') Z, ",", A, ", = Z,A"
    do n = 1, num_p_states
        write(4, '(F10.5,A,F10.5)') E_p(n), "," , E_n(n)
    end do
    do n = num_p_states + 1, num_n_states
        write(4, '(F10.5,A,F10.5)') 0.0_kind, "," , E_n(n)
    end do
    close(4)
    
    
end program