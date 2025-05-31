program main
    use constants
    use micmac
    use mass_table
    implicit none


    real(r_kind) :: BE, ME, ME_exp, BE_exp, BE_exp_unc, ME_exp_unc, params(num_params), BE2
    character(len=3) :: El
    integer :: idx
    integer :: Z, A, N
    Z = 20
    A = 48
    N = A - Z

    params = starting_params
    params(3) = 50000
    BE = binding_energy(params, Z,A)
    params(3)=0
    BE2 = binding_energy(params,Z,A)
    WRITE(*,*) "DIFF:"
    WRITE(*,*) BE-BE2
    ME = mass_excess(BE, Z,A)
    WRITE(*,*)
    WRITE(*,*) "########################################"
    WRITE(*,*) "Read experimental masses"
    WRITE(*,*) "########################################"
    WRITE(*,*)
    call read_ame()
    WRITE(*,*)"Z    A     Mass Excess    unc          BE/A      unc"
    WRITE(*,*)"Z    A     (kev)          (kev)        (kev)     (kev)"
    do idx = 1, num_vals

        !WRITE(*, '(I4, I4,1x, A3, F14.6,F12.6, F13.5, 1x, F10.5)') AME_Z(idx), AME_A(idx), AME_EL(idx), AME_ME(idx), AME_ME_unc(idx), AME_BE(idx), AME_BE_unc(idx)

    end do

    WRITE(*,*) "########################################"
    WRITE(*,*) "MicMac - Binding Energy and Mass Excess"
    WRITE(*,*) "########################################"
    call read_elem(BE_exp, BE_exp_unc, ME_exp, ME_exp_unc, El, Z, A)
    WRITE(*,'(A, I3, A, I3, A, A)') "Z = ", Z, ", A = ", A, ", ", El
    WRITE(*,*) "CALC         EXP      EXP UNC"
    WRITE(*,*) "(MeV)        (MeV)    (keV)"
    WRITE(*,'(F10.5,2x, F10.5, F10.5)')  BE/A, BE_exp/A, 1000*BE_exp_unc/A


    WRITE(*,'(A, F10.3, A)') "Binding Energy = ", BE, " MeV"
    WRITE(*,'(A, F8.3, A)') "Binding Energy per nucleon = ", BE/A, " MeV"
    WRITE(*,*)
    WRITE(*,'(A, F8.3, A)') "Mass Excess = ", ME, " MeV/c^2"
    WRITE(*,'(A, F8.3, A)') "Mass Excess = ", ME/dalton, " u"
    WRITE(*,'(A, F8.3, A)') "Mass Excess per nucleon = ", ME/A, " MeV"



end program main