module fitting
    use constants
    use micmac
    use mass_table
    implicit none
    
    private
    public :: read_fit_exp_data, fit_param_cov, fit_iterative, fit_rms, write_table
    real(r_kind), dimension(:),public,allocatable  :: exp_be, exp_be_unc, exp_me, exp_me_unc
    integer, dimension(:),allocatable, public :: exp_Z, exp_A
    character(len=3), dimension(:), allocatable, public :: exp_elname
    integer, public :: num_fit_vals

    contains
    !!Reads relevant exp data to fit against
    subroutine read_fit_exp_data()
        integer :: idx, Z, A, idx2,N
        real(r_kind) :: unc
        num_fit_vals = 0

        do idx = 1, num_vals
            Z = AME_Z(idx)
            A = AME_A(idx)
            N = A-Z
            unc = AME_BE_unc(idx)
            if(Z >= Z_fit_minval .and. N >= N_fit_minval .and. unc < max_unc_mev ) then
                num_fit_vals = num_fit_vals + 1
            else
                cycle
            endif
        end do
        allocate(exp_Z(num_fit_vals), exp_A(num_fit_vals), exp_elname(num_fit_vals), exp_me(num_fit_vals))
        allocate(exp_me_unc(num_fit_vals), exp_be(num_fit_vals), exp_be_unc(num_fit_vals))
        idx2 = 1
        do idx = 1, num_vals
            Z = AME_Z(idx)
            A = AME_A(idx)
            N = A - Z
            unc = AME_BE_unc(idx)
            if(Z >= Z_fit_minval .and. N >= N_fit_minval .and. unc < max_unc_mev) then

            else
                cycle
            endif
            exp_Z(idx2) = AME_Z(idx)
            exp_A(idx2) = AME_A(idx)
            exp_elname(idx2) = AME_EL(idx)
            exp_me(idx2) = AME_ME(idx)
            exp_me_unc(idx2) = AME_ME_unc(idx)
            exp_be(idx2) = AME_BE(idx)
            exp_be_unc(idx2) = AME_BE_unc(idx)
            idx2 = idx2 + 1
        end do
    end subroutine


    !!Subroutien uses singular value decomposition to invert a square matrix
    subroutine inverse_SVD(A,N)
        real(r_kind), intent(inout) :: A(N,N) !!Matrix to invert A(N,N)
        integer, intent(in) :: N !!Size of matrix to invert

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
    end subroutine inverse_SVD

    subroutine solve_linsys(A,b,dim) !Solve Ax=b for x
        real(r_kind), dimension(dim,dim), intent(in) :: A
        real(r_kind), dimension(dim), intent(inout) :: b !!In: b. Out: Solution
        integer,intent(in) :: dim
        real(r_kind), dimension(dim,1) :: Bmat
        integer, dimension(dim) :: ipiv
        integer ::info
        external :: dgesv
        Bmat(:,1) = b

        ! Solve the linear system A * x = b, with 
        !   XTX (Np, Np).   P: (Np)    Xt(Np,Nv)     BE(Nv)
        !Call LAPACK routine dgesv to solve the linear system
        call dgesv(dim,1,A,dim, ipiv,Bmat,dim,info)
        if(info /= 0) then
            print*, "Error when solving linsys, info = ", info
        endif
        b = Bmat(:,1)
    end subroutine

    !!Covariances of fit parameters
    function fit_param_cov(params, num_nuclei, Zs, As, defs) result(cov)
        real(kind=r_kind), intent(in) :: params(num_params)
        type(deformation), intent(in) :: defs(num_nuclei) 
        real(kind=r_kind) :: cov(num_params,num_params)
        integer, dimension(num_nuclei), intent(in) :: Zs, As
        integer, intent(in) :: num_nuclei
        real(kind=r_kind) ::  J(num_nuclei,num_params), JT(num_params,num_nuclei), JTJ(num_params,num_params), rms, res(num_nuclei), sqsum
        J = J_mat(Zs,As,defs, num_nuclei,params)
        JT = transpose(J)
        JTJ = matmul(JT,J)

        !!Get RMS error
        rms = fit_rms(params)

        call inverse_SVD(JTJ,num_params)
        cov = JTJ * rms**2
    end function

    !!Returns covariances of BE predictions
    function be_cov(params, num_nuclei, Zs, As, defs) 
        real(kind=r_kind), intent(in) :: params(num_params)
        type(deformation), intent(in) :: defs(num_nuclei) 
        real(kind=r_kind) :: cov(num_params,num_params), be_cov(num_nuclei,num_nuclei)
        integer, dimension(num_nuclei), intent(in) :: Zs, As
        integer, intent(in) :: num_nuclei
        real(kind=r_kind) ::  J(num_nuclei,num_params), JT(num_params,num_nuclei), JTJ(num_params,num_params), rms, res(num_nuclei), sqsum
        J = J_mat(Zs,As,defs, num_nuclei,params)
        JT = transpose(J)
        cov = fit_param_cov(params, num_nuclei, zs, as, defs)

        be_cov = matmul(J, matmul(cov, JT))
    end function

    function fit_rms(params)
        real(r_kind), intent(in) :: params(num_params)
        real(r_kind) :: resid(num_fit_vals), sum_sq_err, fit_rms, BEs(num_fit_vals)
        type(deformation) :: defs(num_fit_vals) 
        call find_gs_multiple(BEs,defs, params, exp_Z, exp_A, num_fit_vals)
        resid = exp_be - BEs
        sum_sq_err = dot_product(resid, resid)
        fit_rms = sqrt(sum_sq_err/(num_fit_vals-num_params))
    end function

    subroutine fit_iterative(params, converged)
        real(kind=r_kind), intent(inout) :: params(num_params) !!In: Starting guess of parameters. Out: Converged solution
        !!Assume that we have already read exp data.
        real(kind=r_kind) :: Jmat(num_fit_vals,num_params), resid(num_fit_vals), BEs(num_fit_vals)
        type(deformation) :: defs(num_fit_vals) 
        real(kind=r_kind) :: sum_sq_err, JT(num_params, num_fit_vals), JTJ(num_params, num_params), RHS(num_params), rms, oldrms, bestrms, bestparams(num_params)
        integer :: num_nuclei, maxitr, ii, numitr
        logical, intent(out) :: converged
        converged = .false.
        maxitr =1000
        num_nuclei = num_fit_vals
        bestrms = huge(bestrms)
        WRITE(*,'(9F15.3)') params
        do ii = 1,maxitr
            ! write(*,*) 
            ! write(*,*) 
            ! write(*,*) "############################"
            ! write(*,'(A,I5)') "ITER: ", ii
            ! write(*,'(A,7F10.3)')"Params: ", params
            call find_gs_multiple(BEs, defs, params, exp_Z, exp_A, num_nuclei)
            Jmat =J_mat(exp_Z,exp_A,defs, num_nuclei,params) !!Setup matrix.
            
            resid = exp_be  - BEs
            ! call write_table(params, Jmat)
            JT = transpose(Jmat)
            JTJ = matmul(JT,Jmat)
            RHS = matmul(JTJ, params) + matmul(JT, resid) 
            !!We have to solve J^T J * P' = JTJ*P + JT*resid for P' 
            call solve_linsys(JTJ,RHS,num_params)
            params = params + (RHS-params) / 5.0_r_kind

            rms = fit_rms(params)
            if(rms < bestrms) then
                bestrms = rms
                bestparams = params
                numitr = 0
            else
                numitr = numitr + 1
                if(numitr > 20 .and. abs(rms-bestrms) < 1e-2) then
                    converged = .true.
                    write(*,'(A,I6,A)')"Fit converged after ", II, " iterations"
                    exit
                endif
            endif

            if(abs(oldrms-rms) < 1e-6) then
                write(*,'(A,I6,A)')"Fit converged after ", II, " iterations"
                converged = .true.
                exit
            endif

            if(mod(ii,100) == 0) then
                write(*,*)
                write(*,*)
                write(*,*) "Iteration: ", ii
                WRITE(*,'(8F15.3)') params
                WRITE(*,'(A,F15.7)') "Rms: ", rms

            endif
            oldrms = rms
        end do
        params = bestparams
    end subroutine

    subroutine write_table(params)
        real(r_kind), intent(in) :: params(num_params)
        integer :: idx
        real(r_kind) :: BE, BE_exp, BEA, BEA_exp, UNC, unc_per_a,def
        real(r_kind) :: BEs(num_fit_vals), BEcov(num_fit_vals,num_fit_vals)
        type(deformation) :: defs(num_fit_vals) 
        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,*) "Z    A           MicMac BE/A   UNC          EXP BE/A     DELTA         MicMac BE UNC      EXP BE   DELTA"
        WRITE(*,*) "Z    A           (Mev)                      (Mev)        (MEV)         (Mev)              (Mev)    (MeV)"
        call find_gs_multiple(BEs, defs,params,exp_Z,exp_A,num_fit_vals)
        BEcov = be_cov(params, num_fit_vals, exp_Z, exp_A, defs)
        do idx = 1, num_fit_vals
            def = defs(idx)%alphas(2)
            BE = BEs(idx)
            UNC = sqrt(BEcov(idx,idx))
            unc_per_a = UNC/(exp_A(idx)*1.0)
            BE_exp = exp_be(idx)
            BEA = BE/(exp_A(idx)*1.0)
            BEA_exp = BE_exp/(exp_A(idx)*1.0)
            WRITE(*,'(I4, I4, I4,1x,      A3, 1x,F10.3,1x,F8.4,A,F6.4, 2x,F12.4, 4x,F8.4, 4x,         F10.2, A, F4.2,2x F12.2, 2x,F6.2)') &
            exp_A(idx)-exp_Z(idx), exp_Z(idx), exp_A(idx), exp_elname(idx), def, BEA, ' ± ',   unc_per_a, BEA_exp, ABS(BEA-BEA_exp), BE,' ± ',        UNC,      BE_exp, abs(BE - BE_exp)
        end do
        write(*,*)
        WRITE(*,*) "N    Z    A         beta2   MicMac BE/A             EXP BE/A    DELTA        MicMac BE            EXP BE   DELTA"
        WRITE(*,*) "N    Z    A                 (Mev)                   (Mev)       (MEV)        (Mev)                (Mev)    (MeV)"

    end subroutine write_table
end module fitting