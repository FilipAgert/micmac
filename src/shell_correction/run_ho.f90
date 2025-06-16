program run_ho
    use constants
    use def_ho
    use hamiltonian
    implicit none


    integer, parameter :: max_n = 10
    integer, parameter :: A = 238, Z=92
    type(an_ho_state), dimension(:), allocatable :: states
    type(an_ho_state) :: state1, state2
    type(betadef) :: def 
    type(WS_pot) :: WSpot
    type(inv_rho) :: inv_rho_pot
    type(rhopot) :: rh
    type(VC_pot) :: VCpot
    type(unit) :: unitel
    integer :: n, shelldegen, idx, m, numstates, ii
    real(r_kind), dimension(:,:), allocatable :: Vws, Tkin, H, V, VC, inprod, invrho, rhopotmat
    real(r_kind), dimension(:), allocatable :: E
    real(r_kind), parameter :: V0 = V0_ws
    real(r_kind), parameter ::r0 = r0_p
    real(r_kind), parameter :: kappa=kappa_ws
    real(r_kind) :: hbaromega0, radius, hbaromegaperp, hbaromegaz
    real(r_kind) :: I, Vwsdepth, r, theta, Vcone, vcold
    integer, dimension(:,:), allocatable :: mat
    real(r_kind), dimension(:,:), allocatable :: matr

    def = betadef(beta2 = 0, beta4=0)
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

    allocate(rhopotmat(numstates,numstates),invrho(numstates, numstates), Vws(numstates,numstates), Tkin(numstates,numstates), H(numstates,numstates), Vc(numstates, numstates), inprod(numstates,numstates))
    
    WSpot%def = def
    WSpot%radius = radius
    Wspot%V_0 = Vwsdepth

    VCpot%def = def
    VCpot%radius = radius
    call VCpot%set_charge_dens(Z)

    ! do n = 1, numstates
    !     state1 = states(n)
    !     do m = 1, numstates
    !         state2 = states(m)
    !         inprod(n,m) = mat_elem_axsym(state1, state2, unitel,mass_p,hbaromegaz,hbaromegaperp)
    !         invrho(n,m) = mat_elem_axsym(state1, state2, inv_rho_pot,mass_p,hbaromegaz,hbaromegaperp)
    !         rhopotmat(n,m) =mat_elem_axsym(state1, state2, rh,mass_p,hbaromegaz,hbaromegaperp)
    !        ! Vc(n,m) = mat_elem_axsym(state1, state2, VCpot, mass_p, hbaromegaz, hbaromegaperp)
    !     end do
    !     ! write(*,'(I5,A,I5,A)') n, " out of ", numstates, " rows completed"
    ! end do

    ! do n = 1, numstates !!Use the fact that matrix elements are symmetric.
    !     state1 = states(n)
    !     do m = numstates, n + 1, -1
    !         state2 = states(m)
    !         Vws(n,m) = Vws(m,n)
    !         Vc(n,m) = Vc(m,n)

    !     end do
    !     ! write(*,'(40F7.3)') Vws(n,:)
    ! end do
    ! write(*,*) "Coul"
    ! do n = 1,numstates
    !     ! write(*,'(40F7.3)') Vc(n,:)
    ! end do

    ! Tkin = 0.0_r_kind
    ! write(*,*)
    ! write(*,*) "Kinetic energy"
    ! do n = 1, numstates
    !     Tkin(n,n) = states(n)%kinetic_energy(hbaromegaz,hbaromegaperp)
    !     ! write(*,'(15F10.3)', advance='no') Tkin(n,n)
    ! end do
    ! write(*,*)
    ! ! write(*,*) "Hamiltonian"
    ! open(3,file="data/out/H.dat")
    ! open(4,file="data/out/WS.dat")
    ! open(5,file="data/out/VC.dat")
    ! H = Tkin + Vws + Vc
    ! do n = 1,numstates
    !     do idx = 1,numstates
    !         write(3, '(F10.3)',advance='no') H(n,idx)
    !         write(4, '(F10.3)',advance='no') Vws(n,idx)
    !         write(5, '(F10.3)',advance='no') Vc(n,idx)
    !     end do
    !     write(3,*)
    !     write(4,*)
    !     write(5,*)
    ! end do
    ! close(3)
    ! close(4)
    ! close(5)

    ! allocate(V(numstates,numstates), E(numstates))
    ! call diagonalize(E, V, H)
    ! write(*,'(A)') "n     E_n (MeV)  E_ho (MeV)  |  000    001-    001+     100    002-    002+    010    101-    101+    200"
    ! do n = 1, min(numstates,50)
    !     write(*,'(I3, F10.1, F10.1)',advance='no')n, E(n), states(n)%kinetic_energy(hbaromegaz,hbaromegaperp) * 2 

    !     write(*,'(A, 10F8.3)') "     ", V(n,1:min(numstates,10))
    ! end do



    ! r = 0.0
    ! theta = 0.0
    ! vcold = Vcpot%eval(r,theta)
    ! write(*,'(A,e15.5,A)') "Pot: ", vcold, " at, 0 fm"

    ! r = 0.5
    ! theta = 0.0
    ! Vcone = Vcpot%eval(r,theta)
    ! write(*,'(A,e15.5,A)') "Pot: ", Vcone, " at, 0.5 fm"
    ! write(*,'(A, f10.3)') "ratio: ", vcold/Vcone


    allocate(mat(numstates, numstates), matr(numstates, numstates))
    ! mat = pauli_z(states)
    ! write(*,*) "Test step matrices"
    ! write(*,*)
    ! write(*,*) "pauliz"
    ! write(*,*)" +  -  +  +"
    ! do n = 1, numstates
    !     write(*,'(100I3)') mat(n,:)
    ! end do

    ! mat = pauli_p(states)
    ! write(*,*) "pauli plus"
    ! write(*,*)" +  -  +  +"
    ! do n = 1, numstates
    !     write(*,'(100I3)') mat(n,:)
    ! end do

    ! mat = pauli_m(states)
    ! write(*,*) "pauli minus"
    ! write(*,*)" +  -  +  +"
    ! do n = 1, numstates
    !     write(*,'(100I3)') mat(n,:)
    ! end do

    ! mat = R_mat(states)
    ! write(*,*)
    ! write(*,*) "R minus"
    ! write(*,*)" 0  1  1  0"
    ! write(*,*)
    ! do n = 1, numstates
    !     write(*,'(100I3)') mat(n,:)
    ! end do

    ! mat = Rp_mat(states)
    ! write(*,*)
    ! write(*,*) "R plus"
    ! write(*,*)" 0  1  1  0"
    ! write(*,*)
    ! do n = 1, numstates
    !     write(*,'(100I3)') mat(n,:)
    ! end do

    ! mat = Sp_mat(states)
    ! write(*,*)
    ! write(*,*) "S plus"
    ! write(*,*)" 0  0  0  0"
    ! write(*,*)
    ! do n = 1, numstates
    !     write(*,'(100I3)') mat(n,:)
    ! end do

    ! mat = S_mat(states)
    ! write(*,*)
    ! write(*,*) "S minus"
    ! write(*,*)" 0  0  0  0"
    ! write(*,*)
    ! do n = 1, numstates
    !     write(*,'(100I3)') mat(n,:)
    ! end do

    ! matr = cosphi_m(states)
    ! write(*,*)
    ! write(*,*) "cosphi"
    ! write(*,*)"ml 0   0   -1   -1   1    1    0    0"
    ! write(*,*)"nz 0   0    0    0   0    0    1    1"
    ! write(*,*)"s  -   +    -    +   -    +    -    +"
    ! write(*,*)
    ! do n = 1, numstates
    !     write(*,'(100F5.1)') matr(n,:)
    ! end do

    !hbaromegaz = hbaromegaperp
    matr = matmul(inv_rho_mat(states, hbaromegaz, hbaromegaperp, mass_p),matmul(matmul(pauli_m(states), Rp_mat(states)+S_mat(states)) -matmul(pauli_p(states), R_mat(states)+Sp_mat(states)) - matmul(lz(states),pauli_z(states)), az_mat(states)-azp_mat(states))) !kin_en(states, hbaromegaz, hbaromegaperp)! + transpose(expiphi(states))
    open(3,file="data/out/mat.dat")
    write(*,*)
    write(3,*) "[cosphi, R]"
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

    !matr = Vso_mat(states, A, def, radius, hbaromegaz, hbaromegaperp, mass_p, Vwsdepth)
    ! write(*,*)
    ! write(*,*) "SO"
    ! write(*,*)" 0  1  1  0"
    ! write(*,*)
    ! do n = 1, numstates
    !     write(*,'(100F5.1)') matr(n,:)
    ! end do

end program