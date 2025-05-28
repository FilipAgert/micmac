program fitting
    use constants
    use micmac
    use mass_table
    implicit none
    
    
    real(r_kind), dimension(:),allocatable  :: FIT_BE, FIT_BE_unc, FIT_ME, FIT_ME_unc
    integer, dimension(:),allocatable  :: FIT_Z, FIT_A
    character(len=3), dimension(:), allocatable :: FIT_EL
    real(r_kind), dimension(:,:), allocatable :: X
    real(r_kind) :: params(num_params), val, cov(num_params,num_params), stdev, testinv(2,2)
    integer :: num_fit_vals, idx
    WRITE(*,*)
    WRITE(*,*) "########################################"
    WRITE(*,*) "MicMac - Fitting Binding Energy and Mass Excess"
    WRITE(*,*) "########################################"
    WRITE(*,*)
    WRITE(*,*) "Read experimental masses from AME 2020"
    WRITE(*,*) "########################################"
    call read_ame()
    WRITE(*,*) "Mass reading completed. ", num_vals, " values read."
    WRITE(*,*)
    call read_fit_exp_data()


    WRITE(*,*) "Experimental masses chosen for fit. ", num_fit_vals, " nuclei chosen."
    WRITE(*,*) "########################################"
    WRITE(*,*) "Fitting parameters to experimental data"
    stdev = 1.0

    call fit_linsys(params, cov)
    allocate(X(num_fit_vals, num_params))
    X = BE_mat(FIT_Z, FIT_A, num_fit_vals)
    

    val = RMS(params, BE_mat(FIT_Z, FIT_A, num_fit_vals), num_fit_vals)

    write(*,*)
    write(*,'(A,F10.5,A)') "RMS: ", val, " (MeV)"
    write(*,'(A,5F10.5)') "Parameters: ", params(1), params(2), params(3), params(4), params(5)

    write(*,*)
    write(*,*) "Covariance matrix:"
    do idx = 1, num_params
        write(*,'(5F20.10)') cov(idx,1), cov(idx,2), cov(idx,3), cov(idx,4), cov(idx,5)
    end do

    call write_table(params,cov,X)
    contains
    !!Reads relevant exp data to fit against
    subroutine read_fit_exp_data()
        integer :: idx, Z, A, idx2
        num_fit_vals = 0
        do idx = 1, num_vals
            Z = AME_Z(idx)
            A = AME_A(idx)
            if(Z >= Z_fit_minval .and. A >= A_fit_minval) then
                num_fit_vals = num_fit_vals + 1
            else
                cycle
            endif
        end do
        allocate(FIT_Z(num_fit_vals), FIT_A(num_fit_vals), FIT_EL(num_fit_vals), FIT_ME(num_fit_vals))
        allocate(FIT_ME_unc(num_fit_vals), FIT_BE(num_fit_vals), FIT_BE_unc(num_fit_vals))
        idx2 = 1
        do idx = 1, num_vals
            Z = AME_Z(idx)
            A = AME_A(idx)
            if(Z >= Z_fit_minval .and. A >= A_fit_minval) then

            else
                cycle
            endif
            FIT_Z(idx2) = AME_Z(idx)
            FIT_A(idx2) = AME_A(idx)
            FIT_EL(idx2) = AME_EL(idx)
            FIT_ME(idx2) = AME_ME(idx)
            FIT_ME_unc(idx2) = AME_ME_unc(idx)
            FIT_BE(idx2) = AME_BE(idx)
            FIT_BE_unc(idx2) = AME_BE_unc(idx)
            idx2 = idx2 + 1
        end do
    end subroutine

    real(r_kind) function RMS(parameters, X, num_nuclei) 
        ! This function calculates the loss function for the fitting process
        ! It compares the calculated binding energy with the experimental data
        ! and returns the sum of squared differences
        real(r_kind), intent(in) :: parameters(num_params)
        real(r_kind), intent(in), dimension(num_nuclei, num_params) :: X
        integer :: num_nuclei
        real(r_kind), dimension(num_nuclei) :: BEs, Diff, BEsold, DIFFTOOLD
        RMS = 0.0
        BEs = MATMUL(X,parameters)

        DIFF = BEs-FIT_BE * FIT_A
        DIFFTOOLD = BEsold- BEs

        RMS = dot_product(DIFF,DIFF)
        RMS = RMS/(num_nuclei-num_params)
        RMS = SQRT(RMS)
    end function

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
        WRITE(*,*) "SVD Inverse matrix A:"
        do i = 1, N
            WRITE(*,*) A(i,:)
        end do

    end subroutine inverse_SVD

    subroutine inverse(A,N)
        real(r_kind), intent(inout) :: A(N,N)
        integer, intent(in) :: N
        integer :: IPIV(N)
        integer :: LWORK, INFO
        real(r_kind) :: WORK(N)
        external :: dgetrf, dgetri
        LWORK = N
        WRITE(*,*) "Input matrix A:"
        do idx = 1, N
            WRITE(*,*) A(idx,:)
        end do
        call dgetrf(N,N,A,N,IPIV,INFO)
        if(INFO/=0) then
            WRITE(*,*) "Error in matrix inversion. INFO: ", INFO
            stop
        end if
        call dgetri(N,A,N,IPIV,WORK,LWORK,INFO)
        if(INFO/=0) then
            WRITE(*,*) "Error in matrix inversion. INFO: ", INFO
            stop
        end if
    end subroutine


    subroutine fit_linsys(fit_parameters, covariance)
        real(r_kind), dimension(num_params), intent(out) :: fit_parameters
        real(r_kind), dimension(num_params,num_params), intent(out) :: covariance
        real(r_kind) ::stdev
        ! This subroutine performs the fitting process
        ! It uses a simple gradient descent method to minimize the loss function
        real(r_kind), dimension(num_params,1) :: B
        real(r_kind), dimension(num_fit_vals,num_params) :: X
        real(r_kind), dimension(num_params, num_fit_vals) :: XT
        real(r_kind), dimension(num_params,num_params) ::XTX, XTX2
        real(r_kind), dimension(num_params) :: params2

        integer, dimension(num_params) :: ipiv
        integer ::info, idx
        external :: dgesv

        X = BE_mat(FIT_Z, FIT_A, num_fit_vals)
        XT = TRANSPOSE(X)
        XTX = MATMUL(XT, X)
        B(:,1) = MATMUL(XT, FIT_BE * FIT_A) 

        ! Solve the linear system XTX * parameters = XT * FIT_BE, with 
        !   XTX (Np, Np).   P: (Np)    Xt(Np,Nv)     BE(Nv)
        !Call LAPACK routine dgesv to solve the linear system
        call dgesv(num_params,1,XTX,num_params, ipiv,B,num_params,info)
        if(info /= 0) then
            WRITE(*,*) "Error solving linear system. Info: ", info
            !>          = 0:  successful exit
            !>          < 0:  if INFO = -i, the i-th argument had an illegal value
            !>          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
            !>                has been completed, but the factor U is exactly
            !>                singular, so the solution could not be computed.
            !> 
            stop
        end if

        fit_parameters = B(:,1)
        stdev = RMS(fit_parameters, X, num_fit_vals)
        XTX2=  XTX
        CALL inverse_SVD(XTX,num_params) !!invert XTX to get covariance matrix
        write(*,*) "Inverted XTX matrix:"
        do idx = 1, num_params
            WRITE(*,*) XTX(idx,:)
        end do



        params2 = MATMUL(XTX,MATMUL(XT,FIT_BE*FIT_A))
        write(*,*) "Fitted parameters inv: ", params2   
        write(*,*) "Prev fitted: ", fit_parameters
        covariance = XTX * stdev**2
    end subroutine fit_linsys

    subroutine write_table(parameters, covariance, X)
        real(r_kind), intent(in) :: parameters(num_params), covariance(num_params,num_params), X(num_fit_vals,num_params)
        integer :: idx
        real(r_kind) :: BE, BE_exp, BEA, BEA_exp, UNC
        real(r_kind) :: BEs(num_fit_vals), BEcov(num_fit_vals,num_fit_vals)
        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,*) "Z    A     MicMac BE/A    AME 2020 BE/A "
        WRITE(*,*) "Z    A     (Mev)          (kev)        "
        BEs = MATMUL(X, parameters)
        BEcov = MATMUL(X,MATMUL(covariance,TRANSPOSE(X)))


        do idx = 1, num_fit_vals
            BE = BEs(idx)
            write(*,*) BEcov(idx,idx)
            UNC = SQRT(BEcov(idx,idx))
            BE_exp = FIT_BE(idx) * FIT_A(idx)
            BEA = BE/(FIT_A(idx)*1.0)
            BEA_exp = BE_exp/(FIT_A(idx)*1.0)
            WRITE(*,'(I4, I4,1x, A3, F12.4,2x,F12.4, 2x,F12.4, 4x,F8.4, 2x, F14.2, 2x, F12.2,2x F6.2)') FIT_Z(idx), FIT_A(idx), FIT_EL(idx), BEA, UNC/(FIT_A(IDX)*1.0), BEA_exp,ABS(BEA-BEA_exp), BE, UNC, BE_exp, abs(BE - BE_exp)
        end do
        WRITE(*,*) "Z    A         MicMac BE/A      UNC       AME 2020 BE/A  DELTA   MicMac BE   UNC          AME 2020 BE   DELTA"
        WRITE(*,*) "Z    A         (Mev)                      (Mev)          (MEV)     (Mev)                  (Mev)       (MeV)"

    end subroutine write_table
end program fitting