module pairing
    use constants
    use iso_fortran_env, only:iostat_end
    implicit none
    private 
    public :: Ebcs, pairing_correction, shell_energy, fermi_level_sh, inverse_SVD, solve_BCS_eq, smoothEpair, EbcsD0, pairing_corr_x
    contains
    pure real(kind) function shell_energy(A, levels) result(energy)
        real(kind), dimension(:), intent(in) :: levels
        integer, intent(in) :: A !!Particle number
        integer :: ii, idx
        real(kind) :: e
        idx = 0
        energy = 0.0_kind
        do ii = 1, A
            idx = (ii + 1) / 2 !!1: 2 => 2 / 2 = 1,  2:    3/2 = 1    3:    4/2 = 2...
            e = levels(idx)
            energy = energy + e
        end do
    end function


    pure real(kind) function fermi_level_sh(A, levels) result(fermi)
        real(kind), dimension(:), intent(in) :: levels
        integer, intent(in) :: A !!Particle number
        fermi = levels((A+1)/2) 
    end function
    real(kind) function pairing_correction(levels_p, levels_n, Z, A, leveldens_efp, leveldens_efn)
        real(kind), intent(in) :: leveldens_efp, leveldens_efn
        real(kind), intent(in) :: levels_p(:), levels_n(:)
        integer, intent(in) :: Z, A
        real(kind) :: I, deltaP, deltaN, Gp, Gn, pc_p, pc_n

        I = (A-2.0_kind*Z) / A
        Gp = (g0p + g1p * I) / A
        Gn = (g0n + g1n * I) / A
        deltaN = 5.72/(A-Z)**(1.0_kind/3.0_kind) * exp(-0.119*I -7.89*I**2)
        deltaP = 5.72/Z**(1.0_kind/3.0_kind) * exp(0.119*I -7.89*I**2)
        write(*,*)
        write(*,'(A)')"Proton pairing:"
        pc_p = pairing_corr_x(levels_p, Z, leveldens_efp, Gp, deltaP)
        write(*,'(A,F7.3,A)')"Proton pairing correction:", pc_p, " MeV"
        write(*,*)
        write(*,'(A)')"Neutron pairing:"
        pc_n = pairing_corr_x(levels_n, A-Z, leveldens_efn, Gn, deltaN)
        write(*,'(A,F7.3,A)')"Neutron pairing correction:", pc_n, " MeV"
        pairing_correction = pc_p + pc_n
    end function

    real(kind) function pairing_corr_x(levels, numparts, leveldens, G, deltaavg)
        implicit none
        real(kind), intent(in) :: leveldens
        real(kind), intent(in) :: levels(:)
        integer, intent(in) :: numparts
        real(kind), intent(in) :: G, deltaavg
        real(kind) :: smoothe, E0, epair, efermi, delta
        real(kind), allocatable :: levs(:)
        integer :: blocked_level, effective_part_nbr


         
        if(mod(numparts, 2) == 1) then !!If odd particle nbr, then one level is blocked from pairing due to the extra particle
            blocked_level = (numparts-1)/2 + 1
            allocate(levs(size(levels)-1))
            levs(1:blocked_level-1) = levels(1:blocked_level-1)
            levs(blocked_level:size(levels)-1) = levels(blocked_level+1:size(levels))
            effective_part_nbr = numparts - 1
        else
            allocate(levs(size(levels)))
            levs = levels
            effective_part_nbr = numparts
        endif
        write(*,'(a,f10.3)')"Level density at fermi:", leveldens
        smoothe = smoothEpair(leveldens, deltaavg, effective_part_nbr)
        E0 = EbcsD0(levs, G, effective_part_nbr)
        call solve_BCS_eq(efermi, delta, levs, effective_part_nbr, G)
        epair = Ebcs(levs, efermi, delta, G)

        write(*,'(A,F10.3, A)') "Pairing energy: ", epair - E0, " MeV"
        ! write(*,'(A,F10.3, A)') "Smooth pairing energy: ", smoothe, " MeV"
        pairing_corr_x = epair - E0 - smoothe

    end function

    pure real(kind) function smoothEpair(leveldens, deltaavg, num_parts)
        real(kind), intent(in) :: leveldens, deltaavg
        
        integer, intent(in) :: num_parts
        integer :: num_pairs
        real(kind) :: x, G
        num_pairs = num_parts / 2
        x = leveldens * deltaavg / num_pairs
        G = 1.0_kind / (leveldens * log((sqrt(1.0_kind+x**2) + 1.0_kind)/ x))

        smoothEpair = atan(1.0_kind/x) * G*num_pairs * x / 2.0_kind  - (sqrt(1.0_kind+x**2) - 1.0_kind) * num_pairs**2*1.0_kind / leveldens

    end function

    pure real(kind) function EbcsD0(levels, G, numparts) !!bcs energy at zero pairing gap
        real(kind), dimension(:), intent(in) :: levels
        real(kind), intent(in) ::  G
        integer, intent(in) :: numparts
        if(mod(numparts,2) .eq. 1) then
            !Then one of the occ numbers = 1/2 => - G sum occ^2 => need to adjust by adding
            EbcsD0 = shell_energy(numparts, levels) - G*(numparts-1)*1.0_kind/2.0_kind
            EbcsD0 = EbcsD0 - G*0.25_kind
        else
            EbcsD0 = shell_energy(numparts, levels) - G*numparts*1.0_kind/2.0_kind
        endif
    end function

    real(kind) function Ebcs(levels, efermi, delta, G)
        real(kind), dimension(:), intent(in) :: levels
        real(kind), intent(in) :: efermi, delta, G
        real(kind) :: ev, occ
        integer :: vv
        Ebcs = - delta**2 / G
        do vv = 1, size(levels)
            ev = levels(vv)

            occ = occ_nbr(levels, vv, efermi, delta)
            Ebcs = Ebcs + 2.0_kind * ev * occ - G * occ**2

        end do

    end function
    pure real(kind) function part_num(levels, fermi, delta) !!particle number in BCS equations
        implicit none
        real(kind), dimension(:), intent(in) :: levels
        real(kind), intent(in) :: fermi, delta
        integer :: vv
        part_num = 0
        do vv = 1,size(levels)
            part_num = part_num + 2.0_kind*occ_nbr(levels, vv, fermi, delta)
        end do
    end function

    pure real(kind) function pairing_gap(levels, fermi, delta)
        implicit none
        real(kind), dimension(:), intent(in) :: levels
        real(kind), intent(in) :: fermi, delta
        real(kind) :: ev
        integer :: vv
        pairing_gap = 0
        do vv = 1,size(levels)
            ev = levels(vv)
            pairing_gap = pairing_gap + 1.0_kind/sqrt((ev-fermi)**2 + delta**2)
        end do
    end function

    pure real(kind) function occ_nbr(levels,idx, fermi, delta)
        implicit none
        real(kind), dimension(:), intent(in) :: levels
        real(kind), intent(in) :: fermi, delta
        integer, intent(in) :: idx
        occ_nbr = 0.5_kind * ( 1.0_kind- (levels(idx)-fermi) / sqrt((levels(idx)-fermi)**2 + delta**2))

    end function

    function BCS_jacobian(levels, fermi, delta) result(J) !!f1 = part numb     J = [dn/df, dn/dD ; dg/df, dg/dD]
        implicit none
        real(kind), dimension(:), intent(in) :: levels
        real(kind), intent(in) :: fermi, delta
        real(kind), dimension(2,2) :: J
        real(kind) :: energy_diff, ev, fact
        integer :: vv
        J = 0
        do vv = 1, size(levels)
            ev = levels(vv)
            energy_diff = ev-fermi
            fact =1.0_kind/  (energy_diff**2 + delta**2)**(3.0_kind/2.0_kind)

            J(1,1) = J(1,1) + delta**2 * fact! dn/def
            J(1,2) = J(1,2)  + delta * energy_diff *fact !! dn/d Delta

            J(2,1) = J(2,1) + energy_diff  *fact !! dg/def
            J(2,2) = J(2,2) - delta *fact !! dg/d Delta
        end do
    end function

    subroutine solve_BCS_eq(efermi, delta, levels, num_part, G) !!https://math.stackexchange.com/questions/268991/how-to-solve-simultaneous-equations-using-newton-raphsons-method
        implicit none
        real(kind), dimension(:), intent(in) :: levels
        integer, intent(in) :: num_part
        real(kind), intent(in) :: G !!pairing strength
        real(kind), intent(out) :: efermi !!fermi level
        real(kind), intent(out) :: delta 
        real(kind), dimension(2,2) :: J, invJ !!jacobian
        real(kind), dimension(2) :: x, f, px
        real(kind), parameter :: tol = 1e-9
        logical :: converged
        integer :: iter
        !!use multivariate newtons method
        
        delta = 1.2_kind !!initial guess for Delta
        efermi = fermi_level_sh(num_part, levels) !!initial guess for fermi level
        x = [efermi, delta]
        converged = .false.
        iter = 1
        ! write(*,*) "particle number: ", num_part
        ! write(*,*) "numlevels: ", size(levels)
        ! write(*,'(A,10F10.3)') "levels", levels(1:10)
        do while(.not. converged)
            J = BCS_jacobian(levels, efermi, delta)
            efermi = x(1)
            delta = x(2)
           

            f = [part_num(levels, efermi, delta) - num_part, pairing_gap(levels, efermi, delta) - 2.0_kind / G]
            invJ = J
            call inverse_SVD(invJ,2)

            px = x
            x = x - matmul(invJ,f)
            ! write(*,'(A,2F10.3)') "Fermi, delta ", efermi, delta
            ! write(*,'(A, F10.3)') "Eq (1): ", f(1)
            ! write(*,'(A, F10.3)') "Eq (2): ", f(2)
           ! pause
            if(dot_product(x-px,x-px) < tol) then
                converged = .true.
            endif
            iter = iter + 1
            !pause
        end do
        efermi = x(1)
        delta = x(2)
        write(*,'(A,F8.3, A, F8.3,A,I3,A)') "Gap equations solved with values: Fermi level = ", efermi, " MeV, Delta = ", delta, " MeV, after ", iter, " iterations"

    end subroutine

    !!subroutine uses singular value decomposition to invert a square matrix
    subroutine inverse_SVD(A,N)
        implicit none
        real(kind), intent(inout) :: A(N,N)
        integer, intent(in) :: N

        character(len=1) :: JOBU, JOBVT
        integer :: LDA, LDU, LDVT, INFO, LWORK
        integer :: i
        real(kind) :: S(N)
        real(kind) :: WORK(5*N), U(N,N), VT(N,N), SIGMAINV(N,N), V(N,N)
        external :: dgesvd
        JOBU = 'A'  ! Compute all left singular vectors
        JOBVT = 'A' ! Compute all right singular vectors
        LDA = N
        LDVT = N
        LDU = N
        LWORK = 5*N
        call dgesvd(JOBU, JOBVT, N, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)
        if(INFO /= 0) then
            WRITE(*,*) "Error in SVD computation. INFO: ", INFO
            stop
        end if
        SIGMAINV = 0.0_kind
        do i = 1,N
            SIGMAINV(i,i) = 1.0_kind/S(i)
        end do
        V = transpose(VT)

        A = MATMUL(V, MATMUL(SIGMAINV, TRANSPOSE(U)))
        ! WRITE(*,*) "SVD Inverse matrix A:"
        ! do i = 1, N
        !     WRITE(*,*) A(i,:)
        ! end do

    end subroutine inverse_SVD


end module pairing