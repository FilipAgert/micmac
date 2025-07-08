module pairing
    use constants
    use iso_fortran_env, only:iostat_end
    implicit none
    private 
    public :: Ebcs, pairing_correction, shell_energy, fermi_level_sh, inverse_SVD
    contains
    pure real(r_kind) function shell_energy(A, levels) result(energy)
        real(r_kind), dimension(:), intent(in) :: levels
        integer, intent(in) :: A !!Particle number
        integer :: ii, idx
        real(r_kind) :: e
        idx = 0
        energy = 0.0_r_kind
        do ii = 1, A
            idx = (ii + 1) / 2 !!1: 2 => 2 / 2 = 1,  2:    3/2 = 1    3:    4/2 = 2...
            e = levels(idx)
            energy = energy + e
        end do
    end function


    pure real(r_kind) function fermi_level_sh(A, levels) result(fermi)
        real(r_kind), dimension(:), intent(in) :: levels
        integer, intent(in) :: A !!Particle number
        fermi = levels((A+1)/2) 
    end function
    real(r_kind) function pairing_correction(levels_p, levels_n, Z, A, leveldens_efp, leveldens_efn)
        real(r_kind), intent(in) :: leveldens_efp, leveldens_efn
        real(r_kind), intent(in) :: levels_p(:), levels_n(:)
        integer, intent(in) :: Z, A
        real(r_kind) :: I, deltaP, deltaN, Gp, Gn, pc_p, pc_n

        I = (A-2.0_r_kind*Z) / A
        Gp = (g0p + g1p * I) / A
        Gn = (g0n + g1n * I) / A
        deltaN = 5.72/(A-Z)**(1.0_r_kind/3.0_r_kind) * exp(-0.119*I -7.89*I**2)
        deltaP = 5.72/Z**(1.0_r_kind/3.0_r_kind) * exp(0.119*I -7.89*I**2)

        pc_p = pairing_corr_x(levels_p, Z, leveldens_efp, Gp, deltaP)
        pc_n = pairing_corr_x(levels_n, A-Z, leveldens_efn, Gn, deltaN)
        pairing_correction = pc_p + pc_n
    end function

    real(r_kind) function pairing_corr_x(levels, numparts, leveldens, G, deltaavg)
        implicit none
        real(r_kind), intent(in) :: leveldens
        real(r_kind), intent(in) :: levels(:)
        integer, intent(in) :: numparts
        real(r_kind), intent(in) :: G, deltaavg
        real(r_kind) :: smoothe, E0, epair, efermi, delta
        smoothe = smoothEpair(leveldens, deltaavg, numparts)
        E0 = EbcsD0(levels, G, numparts)
        call solve_BCS_eq(efermi, delta, levels, numparts, G)
        epair = Ebcs(levels, efermi, delta, G)

        write(*,'(A,F10.3)') "Pairing energy: ", epair - E0
        write(*,'(A,F10.3)') "Smooth pairing energy: ", smoothe
        pairing_corr_x = epair - E0 - smoothe

    end function

    pure real(r_kind) function smoothEpair(leveldens, deltaavg, num_parts)
        real(r_kind), intent(in) :: leveldens, deltaavg
        
        integer, intent(in) :: num_parts
        integer :: num_pairs
        real(r_kind) :: x, G
        num_pairs = num_parts / 2
        x = leveldens * deltaavg / num_pairs
        G = 1.0_r_kind / (leveldens * log((sqrt(1.0_r_kind+x**2) + 1.0_r_kind)/ x))

        smoothEpair = atan(1.0_r_kind/x) * G*num_pairs * x / 2.0_r_kind  - (sqrt(1.0_r_kind+x**2) - 1.0_r_kind) * num_pairs**2*1.0_r_kind / leveldens

    end function

    pure real(r_kind) function EbcsD0(levels, G, numparts) !!bcs energy at zero pairing gap
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) ::  G
        integer, intent(in) :: numparts
        EbcsD0 = shell_energy(numparts, levels) + G*numparts*1.0_r_kind/2.0_r_kind
    end function

    pure real(r_kind) function Ebcs(levels, efermi, delta, G)
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: efermi, delta, G
        real(r_kind) :: ev, occ
        integer :: vv
        Ebcs = - delta**2 / G
        do vv = 1, size(levels)
            ev = levels(vv)

            occ = occ_nbr(levels, vv, efermi, delta)
            Ebcs = Ebcs + 2 * ev * occ + G * occ**2

        end do

    end function
    pure real(r_kind) function part_num(levels, fermi, delta) !!particle number in BCS equations
        implicit none
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: fermi, delta
        integer :: vv
        part_num = 0
        do vv = 1,size(levels)
            part_num = part_num + 2*occ_nbr(levels, vv, fermi, delta)
        end do
    end function

    pure real(r_kind) function pairing_gap(levels, fermi, delta)
        implicit none
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: fermi, delta
        real(r_kind) :: ev
        integer :: vv
        pairing_gap = 0
        do vv = 1,size(levels)
            ev = levels(vv)
            pairing_gap = pairing_gap + 1.0_r_kind/sqrt((ev-fermi)**2 + delta**2)
        end do
    end function

    pure real(r_kind) function occ_nbr(levels,idx, fermi, delta)
        implicit none
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: fermi, delta
        integer, intent(in) :: idx
        occ_nbr = 0.5_r_kind * ( 1.0_r_kind- (levels(idx)-fermi) / sqrt((levels(idx)-fermi)**2 + delta**2))

    end function

    function BCS_jacobian(levels, fermi, delta) result(J) !!f1 = part numb     J = [dn/df, dn/dD ; dg/df, dg/dD]
        implicit none
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: fermi, delta
        real(r_kind), dimension(2,2) :: J
        real(r_kind) :: energy_diff, ev, fact
        integer :: vv
        J = 0
        do vv = 1, size(levels)
            ev = levels(vv)
            energy_diff = ev-fermi
            fact =1.0_r_kind/  (energy_diff**2 + delta**2)**(3.0_r_kind/2.0_r_kind)

            J(1,1) = J(1,1) + delta**2 * fact! dn/def
            J(1,2) = J(1,2)  + delta * energy_diff *fact !! dn/d Delta

            J(2,1) = J(2,1) + energy_diff  *fact !! dg/def
            J(2,2) = J(2,2) - delta *fact !! dg/d Delta
        end do
    end function

    subroutine solve_BCS_eq(efermi, delta, levels, num_part, G) !!https://math.stackexchange.com/questions/268991/how-to-solve-simultaneous-equations-using-newton-raphsons-method
        implicit none
        real(r_kind), dimension(:), intent(in) :: levels
        integer, intent(in) :: num_part
        real(r_kind), intent(in) :: G !!pairing strength
        real(r_kind), intent(out) :: efermi !!fermi level
        real(r_kind), intent(out) :: delta 
        real(r_kind), dimension(2,2) :: J, invJ !!jacobian
        real(r_kind), dimension(2) :: x, f, px
        real(r_kind), parameter :: tol = 1e-6
        logical :: converged
        integer :: iter
        !!use multivariate newtons method
        
        delta = 1.0_r_kind !!initial guess for Delta
        efermi = fermi_level_sh(num_part, levels) !!initial guess for fermi level
        write(*,*)"Fermi level shell:", efermi
        x = [efermi, delta]
        converged = .false.
        iter = 1
        write(*,*) "G: ", G
        ! write(*,*) "particle number: ", num_part
        ! write(*,*) "numlevels: ", size(levels)
        ! write(*,'(A,10F10.3)') "levels", levels(1:10)
        do while(.not. converged)
            J = BCS_jacobian(levels, efermi, delta)
            efermi = x(1)
            delta = x(2)
           

            f = [part_num(levels, efermi, delta) - num_part, pairing_gap(levels, efermi, delta) - 2.0_r_kind / G]
            invJ = J
            call inverse_SVD(invJ,2)

            px = x
            x = x - matmul(invJ,f) / 2
            write(*,'(A,2F10.3)') "Fermi, delta ", efermi, delta
            write(*,'(A, F10.3)') "Eq (1): ", f(1)
            write(*,'(A, F10.3)') "Eq (2): ", f(2)
           ! pause
            if(dot_product(x-px,x-px) < tol) then
                converged = .true.
            endif
            iter = iter + 1
            !pause
        end do
        efermi = x(1)
        delta = x(2)
        write(*,'(A,F8.3, A, F8.3,A,I3,A)') "Converged with values e_fermi = ", efermi, " MeV, Delta = ", delta, " MeV, after ", iter, " iterations"

    end subroutine

    !!subroutine uses singular value decomposition to invert a square matrix
    subroutine inverse_SVD(A,N)
        implicit none
        real(r_kind), intent(inout) :: A(N,N)
        integer, intent(in) :: N

        character(len=1) :: JOBU, JOBVT
        integer :: LDA, LDU, LDVT, INFO, LWORK
        integer :: i
        real(r_kind) :: S(N)
        real(r_kind) :: WORK(5*N), U(N,N), VT(N,N), SIGMAINV(N,N), V(N,N)
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
        SIGMAINV = 0.0_r_kind
        do i = 1,N
            SIGMAINV(i,i) = 1.0_r_kind/S(i)
        end do
        V = transpose(VT)

        A = MATMUL(V, MATMUL(SIGMAINV, TRANSPOSE(U)))
        ! WRITE(*,*) "SVD Inverse matrix A:"
        ! do i = 1, N
        !     WRITE(*,*) A(i,:)
        ! end do

    end subroutine inverse_SVD


end module pairing