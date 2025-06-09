program run_ho
    use constants
    use def_ho, only: getnumstatesupto, getnumstates, get_ho_states, an_ho_state, Ws_mat_elem, betadef

    implicit none


    integer, parameter :: max_n = 2
    integer, parameter :: A = 48
    type(an_ho_state), dimension(:), allocatable :: states
    type(an_ho_state) :: state1, state2
    type(betadef) :: def 
    integer :: n, shelldegen, idx, m, numstates
    real(r_kind), dimension(:,:), allocatable :: Vws, Tkin, H
    real(r_kind), parameter :: VwsDepth = 50
    real(r_kind), parameter ::r0 = 1.2
    real(r_kind) :: hbaromega0, radius, hbaromegaperp, hbaromegaz

    def = betadef(beta2 = 0.2, beta4=0.0)
    hbaromega0 = 41.0_r_kind * A**(-1.0_r_kind/3.0_r_kind) !!MeV
    hbaromegaperp = def%omega_perp(hbaromega0) !! omega = Mev/hbar
    hbaromegaz = def%omega_z(hbaromega0)
    radius = r0 * A**(1.0_r_kind/3.0_r_kind)

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
        do m = 1, numstates
            state2 = states(m)
            Vws(n,m) = Ws_mat_elem(state1, state2, def,VwsDepth,radius,mass_p,hbaromegaz,hbaromegaperp)
            if(Vws(n,m) < -1e6) then
                write(*,*) state1%header()
                write(*,*) state1%text()
            endif
        end do
        write(*,'(15F10.3)') Vws(n,:)
    end do

    Tkin = 0.0_r_kind
    write(*,*)
    write(*,*) "Kinetic energy"
    do n = 1, numstates
        Tkin(n,n) = states(n)%kinetic_energy(hbaromegaz,hbaromegaperp)
        write(*,'(15F10.3)') Tkin(n,:)
    end do
    write(*,*)
    write(*,*) "Hamiltonian"
    H = Tkin + Vws
    do n = 1, numstates
        write(*,'(15F10.3)') H(n,:)
    end do

end program