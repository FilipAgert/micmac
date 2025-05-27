program fitting
    use constants
    use micmac
    use mass_table
    implicit none
    
    
    real(r_kind), dimension(4000)  :: FIT_BE, FIT_BE_unc, FIT_ME, FIT_ME_unc
    integer, dimension(4000)  :: FIT_Z, FIT_A
    character(len=3), dimension(4000) :: FIT_EL
    real(r_kind) :: params(num_params)
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
    params = fit_parameters()
    call write_table(params)
    contains
    !!Reads relevant exp data to fit against
    subroutine read_fit_exp_data()
        integer :: idx, Z, A
        num_fit_vals = 0
        do idx = 1, num_vals
            Z = AME_Z(idx)
            A = AME_A(idx)
            if(Z >= Z_fit_minval .and. A >= A_fit_minval) then
                num_fit_vals = num_fit_vals + 1
            else
                cycle
            endif
            
            FIT_Z(num_fit_vals) = AME_Z(idx)
            FIT_A(num_fit_vals) = AME_A(idx)
            FIT_EL(num_fit_vals) = AME_EL(idx)
            FIT_ME(num_fit_vals) = AME_ME(idx)
            FIT_ME_unc(num_fit_vals) = AME_ME_unc(idx)
            FIT_BE(num_fit_vals) = AME_BE(idx)
            FIT_BE_unc(num_fit_vals) = AME_BE_unc(idx)
        end do
    end subroutine

    real(r_kind) function RMS(parameters) 
        ! This function calculates the loss function for the fitting process
        ! It compares the calculated binding energy with the experimental data
        ! and returns the sum of squared differences
        real(r_kind), intent(in) :: parameters(num_params)
        integer :: idx
        real(r_kind) :: BE, BE_exp
        RMS = 0.0
        do idx = 1, num_fit_vals
            BE = binding_energy(parameters, FIT_Z(idx), FIT_A(idx))
            BE_exp = FIT_BE(idx) * FIT_A(idx)
            RMS = RMS + ((BE - BE_exp))**2
        end do
        RMS = RMS/num_fit_vals
        RMS = SQRT(RMS)

    end function

    function gradient(parameters) result(grad)
        ! This function calculates the gradient of the loss function
        ! with respect to the parameters
        real(r_kind), intent(in) :: parameters(num_params)
        real(r_kind), dimension(num_params) :: grad
        integer :: i
        real(r_kind) :: delta
        real(r_kind) :: fac(num_params)
        delta = 1.0e-6
        fac = 0
        do i = 1, num_params
            fac(i) = delta
            grad(i) = (RMS(parameters + fac) - RMS(parameters)) / delta
            fac(i) = 0
        end do
    end function

    function fit_parameters()
        real(r_kind), dimension(num_params) :: fit_parameters
        ! This subroutine performs the fitting process
        ! It uses a simple gradient descent method to minimize the loss function
        real(r_kind), dimension(num_params) :: parameters, grad, delta, prevParam, prevGrad, graddiff
        integer :: i,  max_iter
        real(r_kind) :: learning_rate, tolerance, val, stepsize

        ! Initialize parameters
        parameters = [-10.0,10.0, 0.1, 10.0,20.0]![-16.08773,  18.31489,   0.74291,  25.21314,  51.36916]!standard_values - 5 ! Initial guess for the parameters
        learning_rate = 0.000001
        tolerance = 1.0e-3
        max_iter = 1000000
        delta = [208.0, 35.0, 1135.0, 9.3, 0.02]

        WRITE(*,*) parameters
        val = RMS(parameters)

        do i = 1, max_iter

            grad = gradient(parameters)    
            
            if(i > 1) then
                graddiff = grad-prevgrad
                stepsize = abs(sum((parameters-prevParam)*(graddiff),1))/sum(graddiff**2,1)

            else
                stepsize = learning_rate 
            end if
            stepsize = min(stepsize, .05)
            
            
            prevParam = parameters 
            prevGrad = grad
            parameters = parameters - stepsize*grad
            
            VAL = RMS(parameters)
            if(MOD(i, 10000) == 0) then
                WRITE(*,*)
                WRITE(*,'(A,I8)') "Iteration: ", i
                WRITE(*,'(A,D15.5)') "stepsize: ", stepsize
                WRITE(*,'(A,5F10.5)') "PARAMS", parameters(1), parameters(2), parameters(3), parameters(4), parameters(5)
                WRITE(*,'(A,5F10.5)') "GRAD  ", grad(1), grad(2), grad(3), grad(4), grad(5)
                WRITE(*,'(A,1F15.4)') "RMS", VAL
            end if

        end do

        WRITE(*,*) "Standard parameters:"
        val = RMS(standard_values)
        WRITE(*,*) standard_values, ", RMS = ", val
        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,*) "Fitted Parameters:"
        val = RMS(parameters)
        WRITE(*,*) parameters, ", RMS = ", val
        fit_parameters = parameters
    end function

    subroutine write_table(parameters)
        real(r_kind), intent(in) :: parameters(num_params)
        integer :: idx
        real(r_kind) :: BE, BE_exp, val, BEA, BEA_exp
        WRITE(*,*)
        WRITE(*,*)
        WRITE(*,*) "Z    A     MicMac BE/A    AME 2020 BE/A "
        WRITE(*,*) "Z    A     (Mev)          (kev)        "
        do idx = 1, num_fit_vals
            BE = binding_energy(parameters, FIT_Z(idx), FIT_A(idx))
            BE_exp = FIT_BE(idx) * FIT_A(idx)
            BEA = BE/(FIT_A(idx)*1.0)
            BEA_exp = BE_exp/(FIT_A(idx)*1.0)
            WRITE(*,'(I4, I4,1x, A3, F12.4, 2x,F12.4, 4x,F8.4, 2x, F14.2, 2x, F12.2,2x F6.2)') FIT_Z(idx), FIT_A(idx), FIT_EL(idx), BEA, BEA_exp,ABS(BEA-BEA_exp), BE, BE_exp, abs(BE - BE_exp)
        end do
        WRITE(*,*) "Z    A         MicMac BE/A    AME 2020 BE/A  DELTA   MicMac BE    AME 2020 BE   DELTA"
        WRITE(*,*) "Z    A         (Mev)          (Mev)          (MEV)     (Mev)          (Mev)       (MeV)"
        val = RMS(parameters)
        WRITE(*,*)
        WRITE(*,*) parameters, ", RMS = ", val, " MeV"
    end subroutine write_table
end program fitting