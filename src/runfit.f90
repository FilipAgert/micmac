program runfit
    use constants
    use micmac
    use mass_table
    use fitting
    implicit none
    real(kind=r_kind) :: params(num_params), param_cov(num_params,num_params)
    real(r_kind), allocatable, dimension(:) :: BEs
    type(deformation), allocatable, dimension(:):: defs 
    integer :: i
    logical :: converged
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
    params = fitted_params
    call fit_iterative(params, converged)
    allocate(BEs(num_fit_vals))
    allocate(defs(num_fit_vals))

    call find_gs_multiple(BEs, defs, params, exp_Z, exp_A, num_fit_vals)
    param_cov = fit_param_cov(params, num_fit_vals, exp_Z,exp_A, defs)

    call write_table(params)

    write(*,*)
    if(.not. converged) then
        write(*,*) "ERROR: Fit did not converge"
    endif
    WRITE(*,'(A,F10.3, A)') "Rms: ", fit_rms(params), " (MeV)"

    write(*,*)
    write(*,*) "         a_vol     a_sur    k_v         r0        C         c         adef"
    write(*,*) "         (MeV)     (MeV)                (fm)      (MeV)               (fm)  "
    write(*,'(A, 7F10.3)') "val: ", params
    write(*,'(A, 7F10.3)') "unc: ", sqrt(param_cov(1,1)), sqrt(param_cov(2,2)), sqrt(param_cov(3,3)), sqrt(param_cov(4,4)), &
    sqrt(param_cov(5,5)), sqrt(param_cov(6,6)), sqrt(param_cov(7,7))
    write(*,'(A, 7F10.3)') "re;unc: ", sqrt(param_cov(1,1))/abs(params(1)), sqrt(param_cov(2,2))/params(2), sqrt(param_cov(3,3))/params(3), sqrt(param_cov(4,4))/params(4), &
    sqrt(param_cov(5,5))/params(5), sqrt(param_cov(6,6))/params(6), sqrt(param_cov(7,7))/params(7)

    write(*,*)
    write(*,*) "Parameter covariance matrix:"
    do i = 1,num_params
        write(*,'(8F10.3)') param_cov(i,:)
    end do

end program runfit