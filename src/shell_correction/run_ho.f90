program run_ho
    use constants
    use def_ho, only: an_ho_state, betadef, getnumstates, getnumstatesupto, kin_en, get_ho_states, get_ho_states_upto
    use hamiltonian
    implicit none


    integer, parameter :: A = 208, Z=82
    type(an_ho_state), dimension(:), allocatable :: states
    type(betadef) :: def 
    type(VC_pot) :: pot
    integer :: n, shelldegen, idx, numstates, j
    real(r_kind), dimension(:,:), allocatable :: Vws, Tkin, V, VC,Vso,Hn,Hp
    real(r_kind), dimension(:), allocatable :: E, E_n, E_p
    real(r_kind), parameter :: V0 = V0_ws
    real(r_kind), parameter ::r0 = r0_p
    real(r_kind), parameter :: kappa=kappa_ws
    real(r_kind) :: hbaromega0, hbaromegaperp, hbaromegaz, rad,Vwsdepth,matelem
    character :: p
    real(r_kind), dimension(:,:), allocatable :: matr

    def = betadef(beta2 = 0.0, beta4=0.0)
    hbaromega0 = 41.0_r_kind * A**(-1.0_r_kind/3.0_r_kind) !!MeV
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
        write(4, '(F10.5,A,F10.5)') 0.0_r_kind, "," , E_n(n)
    end do
    close(4)

    allocate(states(getnumstatesupto(max_n)))
    states = get_ho_states_upto(max_n)
    ! Vc = H_neutrons(states, Z, A, def, hbaromegaz, hbaromegaperp)!coul_mat(states, def, rad, Z,mass_p,hbaromegaz, hbaromegaperp)!coul_mat(states, def, rad, Z, mass_p, hbaromegaz, hbaromegaperp)
    allocate(matr(size(Vc,1),size(vc,2)))
    matr = Vc
    numstates = size(states)

    open(3,file="data/out/mat.dat")
    write(3,*)
    write(3,'(14x, A5)',advance='no') "parity"
    do n = 1, numstates
        if(states(n)%pi < 0) then
            p = '-'
        else
            p = '+'
        endif
        write(3,'(a7)',advance='no') p
    end do
    write(3,*)

    write(3,'(14x, A5)',advance='no')"ml   "
    do n = 1, numstates
        write(3,'(I7)',advance='no') states(n)%ml
        !write(,*) states(n)%nr
    end do
    write(3,*)
    write(3,'(14x, A5)',advance='no')"mj   "
    do n = 1, numstates
        write(3,'(I7)',advance='no') states(n)%mj
    end do

    write(3,*)
    write(3,*)

    do n = 1, numstates
        if(states(n)%pi < 0) then
            p = '-'
        else
            p = '+'
        endif
        write(3,'(a5,2I5, 5x)', advance='no') p, states(n)%ml,states(n)%mj
        do j = 1, numstates
            if (abs(matr(n,j)) < 1e-8) then
                write(3,'(7x)', advance='no')
            else
                write(3,'(F7.3)', advance='no') matr(n,j)
            endif
            
        end do
        write(3,'(A)')""
    end do
    close(3)

    
end program