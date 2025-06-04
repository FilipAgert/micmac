program run_fissbarr
    use  constants, only: fitted_params, r_kind, num_params
    use fiss_barr
    use micmac, only: alpha_to_beta

    implicit none

    integer, parameter :: Z = 92, A = 238
    real(r_kind) :: params(num_params)
    real(r_kind) :: defBf, Bf, betaBf
    params = fitted_params
    ! call find_fiss_barr(Bf, defBf, params, Z, A)


    write(*,'(A,I3,A,I3)') "Z = ", Z, ", A = ", A
    betaBf = alpha_to_beta(defBf)
    write(*,'(A,F10.3,A,F5.3)') "Bf = ", Bf, " (MeV), beta2Bf = ", betaBf


end program run_fissbarr