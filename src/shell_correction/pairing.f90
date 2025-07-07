module pairing
    use constants
    use strutinsky, only: fermi_level_sh
    implicit none
    private


    contains

    pure real(r_kind) function part_num(levels, fermi, delta) !!particle number in BCS equations
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: fermi, delta
        real(r_kind) :: ev
        integer :: vv
        part_num = 0
        do vv = 1,size(levels)
            ev = levels(vv)
            part_num = part_num + 1 - (ev-fermi)/sqrt((ev-fermi)**2 + delta**2)
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

    subroutine solve_BCS_eq(efermi, delta, levels, num_part, G) !!https://math.stackexchange.com/questions/268991/how-to-solve-simultaneous-equations-using-newton-raphsons-method
        real(r_kind), dimension(:), intent(in) :: levels
        integer, intent(in) :: num_part
        real(r_kind), intent(in) :: G !!pairing strength
        real(r_kind), intent(out) :: efermi !!fermi level
        real(r_kind), intent(out) :: delta 
        real(r_kind), dimension(2,2) :: J !!jacobian
        real(r_kind), dimension(2) :: x, f
        !!use multivariate newtons method
        delta = 1 !!initial guess
        efermi = fermi_level_sh(num_part, levels) !!initial guess
        x = [delta, efermi]
        





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


end module pairing