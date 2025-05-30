module constants
    implicit none


    integer, parameter :: r_kind = 8
    integer, parameter :: num_params = 6
    real(r_kind), parameter :: standard_values(num_params) = [-15.68, 18.56, 28.0,28.0,0.717, 34.0]
    real(r_kind), parameter :: mass_p = 938.27208943
    real(r_kind), parameter :: mass_e = 0.5109989461
    real(r_kind), parameter :: mass_n = 939.56542194
    real(r_kind), parameter :: dalton = 931.49410372

    integer, parameter :: Z_fit_minval = 84
    integer, parameter :: N_fit_minval = 126
    real(r_kind), parameter :: max_unc_mev = 0.15_r_kind


    integer, parameter, dimension(9) :: magic_num_Z = [0,2,8,20,28,50,82,114,164]
    integer, parameter, dimension(10) :: magic_num_N = [0,2,8,20,28,50,82,126,184,228]
end module constants
