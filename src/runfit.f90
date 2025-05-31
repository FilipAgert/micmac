program runfit
    use constants
    use micmac
    use mass_table
    use table_writer
    use fitting
    implicit none
    real(kind=r_kind) :: params(num_params), param_cov(num_params,num_params)
    integer :: i
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
    params = starting_params
    call fit_iterative(params)
    param_cov = fit_param_cov(params, num_fit_vals, exp_Z,exp_A)
    write(*,*)
    write(*,*) "        a_vol      a_sur       k         r0        C         c"
    write(*,*) "       (MeV)      (MeV)                (fm)      (MeV)       "
    write(*,'(A, 6F10.3)') "val: ", params
    write(*,'(A, 6F10.3)') "unc: ", param_cov(1,1), param_cov(2,2), param_cov(3,3), param_cov(4,4), param_cov(5,5), param_cov(6,6)

end program runfit