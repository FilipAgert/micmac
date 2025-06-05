module constants
    implicit none


    integer, parameter :: r_kind = 8
    
    integer, parameter :: num_params = 7
    !! volume, surf, coul, asym, pairing: -15.46541  15.94866   0.71158  24.03928  50.66491
    real(r_kind), parameter :: starting_params(num_params) = [-15.63695, 18.40844, 1.78978,  1.21103,3.69404,0.22215, 0.28959] ![-16.786,29.36,1.926,1.2,2.75,9.651] my minimisation
    real(r_kind), parameter :: fitted_params(num_params) = [-15.398,17.722, 1.774, 1.246, 6.643, 0.239, 0.431] !!Z .GE. 84. RMS = 1.037
    real(r_kind), parameter :: cauchois_params(num_params) = [-15.63695,18.40844,1.78709,1.21103,3.69404,0.2215,0.28959] !!Z .GE. 84. RMS = 1.037
    real(r_kind), parameter :: mass_p = 938.27208943
    real(r_kind), parameter :: mass_e = 0.5109989461
    real(r_kind), parameter :: mass_n = 939.56542194
    real(r_kind), parameter :: dalton = 931.49410372
    real(r_kind), parameter :: d =0.5461
    real(r_kind), parameter :: e_squared = 1.4399764_r_kind !! MeV fm
    real(r_kind), parameter :: pi = ACOS(-1.0_r_kind)
    real(r_kind), parameter :: pi2 = pi*pi

    integer, parameter :: Z_fit_minval = 84!!84
    integer, parameter :: N_fit_minval = 0!8 !!126
    real(r_kind), parameter :: max_unc_mev = 0.15_r_kind


    integer, parameter, dimension(9) :: magic_num_Z = [0,2,8,20,28,50,82,126,184]
    integer, parameter, dimension(10) :: magic_num_N = [0,2,8,20,28,50,82,126,184,258]



    !!Finding G.S. minimization parameters
    integer, parameter :: default_num_restarts = 50
    real(r_kind), parameter :: gamma_damp_fac = 0.5_r_kind
    real(r_kind), parameter :: default_def_bounds = 0.5_r_kind
end module constants
