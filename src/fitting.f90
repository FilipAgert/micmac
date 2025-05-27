program fitting
    use constants
    use micmac
    use mass_table
    implicit none
    
    
    real(r_kind), dimension(:),allocatable  :: FIT_BE, FIT_BE_unc, FIT_ME, FIT_ME_unc
    integer, dimension(:),allocatable  :: FIT_Z, FIT_A
    character(len=3), dimension(:), allocatable :: FIT_EL
    real(r_kind) :: params(num_params), val
    integer :: num_fit_vals
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
    params = fit_linsys()
    val = RMS(params, BE_mat(FIT_Z, FIT_A, num_fit_vals), num_fit_vals)
    call write_table(params)
    write(*,*)
    write(*,'(A,F10.5,A)') "RMS: ", val, " (MeV)"
    write(*,'(A,5F10.5)') "Parameters: ", params(1), params(2), params(3), params(4), params(5)
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

    subroutine inverse(A,N)
        real(r_kind), intent(inout) :: A(N,N)
        integer, intent(in) :: N
        integer :: IPIV(N)
        integer :: LWORK, INFO
        real(r_kind) :: WORK(N*N)
        external :: dgetrf, dgetri
        LWORK = N * N
        call dgetrf(N,N,A,N,IPIV,INFO)
        call dgetri(N,A,N,IPIV,WORK,LWORK,INFO)

        if(INFO/=0) then
            WRITE(*,*) "Error in matrix inversion. INFO: ", INFO
            stop
        end if
    end subroutine

    function fit_linsys()
        real(r_kind), dimension(num_params) :: fit_linsys
        ! This subroutine performs the fitting process
        ! It uses a simple gradient descent method to minimize the loss function
        real(r_kind), dimension(num_params,1) :: B
        real(r_kind), dimension(num_fit_vals,num_params) :: X
        real(r_kind), dimension(num_params, num_fit_vals) :: XT
        real(r_kind), dimension(num_params,num_params) ::XTX
        integer, dimension(num_params) :: ipiv
        integer ::info
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
        fit_linsys = B(:,1)
    end function fit_linsys

    subroutine write_table(parameters)
        real(r_kind), intent(in) :: parameters(num_params)
        integer :: idx
        real(r_kind) :: BE, BE_exp, BEA, BEA_exp
        real(r_kind) :: BEs(num_fit_vals)
        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,*) "Z    A     MicMac BE/A    AME 2020 BE/A "
        WRITE(*,*) "Z    A     (Mev)          (kev)        "
        BEs = binding_energies(parameters, FIT_Z,FIT_A,num_fit_vals)
        do idx = 1, num_fit_vals
            BE = BEs(idx)
            BE_exp = FIT_BE(idx) * FIT_A(idx)
            BEA = BE/(FIT_A(idx)*1.0)
            BEA_exp = BE_exp/(FIT_A(idx)*1.0)
            WRITE(*,'(I4, I4,1x, A3, F12.4, 2x,F12.4, 4x,F8.4, 2x, F14.2, 2x, F12.2,2x F6.2)') FIT_Z(idx), FIT_A(idx), FIT_EL(idx), BEA, BEA_exp,ABS(BEA-BEA_exp), BE, BE_exp, abs(BE - BE_exp)
        end do
        WRITE(*,*) "Z    A         MicMac BE/A    AME 2020 BE/A  DELTA   MicMac BE    AME 2020 BE   DELTA"
        WRITE(*,*) "Z    A         (Mev)          (Mev)          (MEV)     (Mev)          (Mev)       (MeV)"

    end subroutine write_table
end program fitting