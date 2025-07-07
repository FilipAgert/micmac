program pairing
    use constants
    use strutinsky, only: fermi_level_sh, shell_energy
    use iso_fortran_env, only:iostat_end
    implicit none
    
    real(r_kind), allocatable :: levels_p(:), levels_n(:), temp_p(:), temp_n(:)
    integer :: ii, Z, A, iostat
    character(len=500) :: trash, line
    real(r_kind) :: ef, delta, G, I, E
    allocate(temp_n(5000), temp_p(5000))
    open(3, file='data/out/levels.dat', status='old', action='read')
    read(3, '(I3,A1,I3,A)', iostat = iostat) Z, trash, A, trash
    ii = 0
    do 
        read(3, '(A)', iostat = iostat) line
        Select case(iostat)
        case(0)
            ii = ii + 1

            read(line,'(F10.5,A,F10.5)')temp_p(ii), trash, temp_n(ii)
        case(iostat_end)
            exit
        case default
            write(*,*) "error in reading file"
            exit
        end select
    end do

    allocate(levels_p(Z), levels_n(A-Z))
    levels_n = temp_n(1:Z)
    levels_p = temp_p(1:(A-Z))
    deallocate(temp_n, temp_p)
    I = (A-2.0_r_kind*Z) / A
    G = (g0p + g1p * I) / A
    write(*,'(A,2I4)') "Pairing for nucleus:", Z, A
    write(*,'(A,F10.3,A,F10.3)') "G:", G, ", 2/G:", 2.0_r_kind / G
    call solve_BCS_eq(ef, delta, levels_p, Z, G)
    E = Ebcs(levels_p, ef, delta, G)
    write(*,'(A,F10.3)') "BCS energy: ", E


    contains

    pure real(r_kind) function EbcsD0(levels, G, numparts) !!bcs energy at zero pairing gap
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) ::  G
        integer, intent(in) :: numparts
        EbcsD0 = shell_energy(numparts, levels) - G*numparts/2
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
            Ebcs = Ebcs + 2 * ev * occ - G * occ**2

        end do

    end function
    pure real(r_kind) function part_num(levels, fermi, delta) !!particle number in BCS equations
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: fermi, delta
        real(r_kind) :: ev
        integer :: vv
        part_num = 0
        do vv = 1,size(levels)
            ev = levels(vv)
            part_num = part_num + 1.0_r_kind - (ev-fermi)/sqrt((ev-fermi)**2 + delta**2)
        end do
    end function

    pure real(r_kind) function pairing_gap(levels, fermi, delta)
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
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: fermi, delta
        integer, intent(in) :: idx
        occ_nbr = 0.5_r_kind * ( 1- (levels(idx)-fermi) / sqrt((levels(idx)-fermi)**2 + delta**2))

    end function

    function BCS_jacobian(levels, fermi, delta) result(J) !!f1 = part numb     J = [dn/df, dn/dD ; dg/df, dg/dD]
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
        real(r_kind), dimension(:), intent(in) :: levels
        integer, intent(in) :: num_part
        real(r_kind), intent(in) :: G !!pairing strength
        real(r_kind), intent(out) :: efermi !!fermi level
        real(r_kind), intent(out) :: delta 
        real(r_kind), dimension(2,2) :: J, invJ !!jacobian
        real(r_kind), dimension(2) :: x, f, px
        real(r_kind), parameter :: tol = 1e-6
        logical :: converged
        !!use multivariate newtons method
        delta = 1.3_r_kind !!initial guess for Delta
        efermi = fermi_level_sh(num_part, levels) !!initial guess for fermi level
        x = [efermi, delta]
        converged = .false.
        do while(.not. converged)
            J = BCS_jacobian(levels, efermi, delta)
            efermi = x(1)
            delta = x(2)
            f = [part_num(levels, efermi, delta) - num_part, pairing_gap(levels, efermi, delta) - 2.0_r_kind / G]
            invJ = J
            call inverse_SVD(invJ,2)

            px = x
            x = x - matmul(invJ,f)
            ! write(*,'(A,2F10.3)') "Fermi, delta ", efermi, delta
            ! write(*,'(A, F10.3)') "Eq (1): ", f(1)
            ! write(*,'(A, F10.3)') "Eq (2): ", f(2)
            ! write(*,'(A, 2F10.3)') "Occupation number at lowest level & fermi level", occ_nbr(levels, 1, efermi, delta), occ_nbr(levels, num_part/2, efermi, delta)
           ! pause
            if(dot_product(x-px,x-px) < tol) then
                converged = .true.
            endif
        end do
        efermi = x(1)
        delta = x(2)
        write(*,'(A,F8.3, A, F8.3)') "Converged with values e_f = ", efermi, ", Delta = ", delta
        write(*,'(A,F8.3, A, F8.3)') "Sum: Particle number:", part_num(levels, efermi, delta), ", pairing gap ", pairing_gap(levels, efermi, delta)
        write(*,'(A,I8, A, F8.3)') "RHS: Particle number:", num_part, ", pairing gap ", 2.0_r_kind/G

    end subroutine

    !!subroutine uses singular value decomposition to invert a square matrix
    subroutine inverse_SVD(A,N)
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


end program pairing