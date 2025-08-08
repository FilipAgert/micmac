module constants
    implicit none


    integer, parameter :: kind = 8
    integer, parameter :: i_kind = 8
    
    integer, parameter :: num_params = 7
    !! volume, surf, coul, asym, pairing: -15.46541  15.94866   0.71158  24.03928  50.66491
    real(kind), parameter :: starting_params(num_params) = [-15.63695, 18.40844, 1.78978,  1.21103,3.69404,0.22215, 0.28959] ![-16.786,29.36,1.926,1.2,2.75,9.651] my minimisation
    real(kind), parameter :: fitted_params(num_params) = [-15.432, 17.886, 1.778, 1.244, 6.601, 0.246, 0.396] !!Z .GE. 84. RMS = 1.037
    real(kind), parameter :: cauchois_params(num_params) = [-15.63695,18.40844,1.78709,1.21103,3.69404,0.2215,0.28959] !!Z .GE. 84. RMS = 1.037
    real(kind), parameter :: mass_p = 938.27208943
    real(kind), parameter :: mass_e = 0.5109989461
    real(kind), parameter :: mass_n = 939.56542194
    real(kind), parameter :: dalton = 931.49410372
    real(kind), parameter :: d =0.5461 !! fm
    real(kind), parameter :: e_squared = 1.4399764_kind !! MeV fm
    real(kind), parameter :: pi = ACOS(-1.0_kind)
    real(kind), parameter :: pi2 = pi*pi
    real(kind), parameter :: epsilonzero = 0.0552634936 !! e^2 / (MeV * fm)

    integer, parameter :: Z_fit_minval = 84!!84
    integer, parameter :: N_fit_minval = 0!8 !!126

    logical, parameter :: use_ml_sym = .true.
    real(kind), parameter :: max_unc_mev = 0.15_kind


    integer, parameter, dimension(9) :: magic_num_Z = [0,2,8,20,28,50,82,126,184]
    integer, parameter, dimension(10) :: magic_num_N = [0,2,8,20,28,50,82,126,184,258]

    real(kind), parameter :: hbarc = 197.327! [MeV fm]


    !!Finding G.S. minimization parameters
    integer, parameter :: default_num_restarts = 50
    real(kind), parameter :: gamma_damp_fac = 0.5_kind
    real(kind), parameter :: default_def_bounds = 0.5_kind
    integer, parameter :: num_def_params = 3


    integer, parameter :: num_threads = 1


    !! Shell model parameters
    real(kind), parameter :: kappa_ws = 0.86_kind !
    real(kind), parameter :: V0_ws = 49.6_kind !MeV
    real(kind), parameter :: a_ws = 0.70_kind !fm


    real(kind), parameter :: r0_p = 1.275_kind !fm. for proton
    real(kind), parameter :: r0_n = 1.347_kind !fm. for neutron
    real(kind), parameter :: r0_so_p = 1.32_kind !fm. for proton
    real(kind), parameter :: r0_so_n = 1.31_kind !fm. for neutron
    real(kind), parameter :: lambda_p = 36.0_kind ! for proton
    real(kind), parameter :: lambda_n = 35.0_kind !for neutron


    integer, parameter :: N_max =19
    integer, parameter :: num_p_states = 450
    integer, parameter :: num_n_states = 550
    logical :: printflag = .false.

    integer, parameter :: nquad =64
    !!pairing parameters
    real(kind), parameter :: g0p = 13.40
    real(kind), parameter :: g1p = 44.89
    real(kind), parameter :: g0n = 17.67
    real(kind), parameter :: g1n = -13.11



end module constants
