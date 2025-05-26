module constants
    implicit none


    integer, parameter :: r_kind = 8
    integer, parameter :: num_params = 5
    real(r_kind), parameter :: standard_values(num_params) = [-15.68, 18.56, 0.717, 28.1, 34.0]
    real(r_kind), parameter :: mp = 938.27208943
    real(r_kind), parameter :: mn = 939.56542194
    real(r_kind), parameter :: dalton = 931.49410372

    integer, parameter :: Z_fit_minval = 82, A_fit_minval = 208
end module constants
