module micmac

    use constants
    implicit none


    private
    public :: binding_energy, mass_excess
    contains

pure function mass_excess(BE, Z, A) result(ME)
    real(r_kind), intent(in) :: BE
    integer, intent(in) :: Z, A
    real(r_kind) :: ME
    integer :: N
    N = A - Z

    ME = BE + (N*mn + Z*mp) - dalton*(N+Z)

end function mass_excess    


pure function binding_energy(parameters, Z, A) result(BE)
    ! This function calculates the binding energy of a nucleus
    real(r_kind), intent(in) :: parameters(num_params)
    integer, intent(in) :: Z, A
    real(r_kind) :: BE
    real(r_kind) :: av, as, ac, ai, ap
    integer :: N
    N = A - Z
    ! Unpack parameters
    av = parameters(1)
    as = parameters(2)
    ac = parameters(3)
    ai = parameters(4)
    ap = parameters(5)
    ! Calculate binding energy
    BE = volume_term(av, N, Z) + &
            surface_term(as, N, Z) + &
            coulomb_term(ac, N, Z) + &
            asymmetry_term(ai, N, Z) - &
            pairing_term(ap, N, Z)
end function binding_energy

pure real(r_kind) function volume_term(av, N, Z)
    real(r_kind), intent(in) :: av
    integer, intent(in) :: N, Z
    integer :: A
    A = N + Z
    volume_term = av * A
end function volume_term

pure real(r_kind) function surface_term(as, N, Z)
    real(r_kind), intent(in) :: as
    integer, intent(in) :: N, Z
    real(r_kind), parameter :: two_thirds = 2.0/3.0
    integer :: A
    A = N + Z
    surface_term = as * A**two_thirds
end function surface_term
pure real(r_kind) function coulomb_term(ac, N, Z)
    real(r_kind), intent(in) :: ac
    integer, intent(in) :: Z, N
    integer :: A
    real(r_kind), parameter :: one_third = 1.0/3.0
    A = N + Z
    coulomb_term = ac * Z*Z / A**one_third
end function coulomb_term

pure real(r_kind) function asymmetry_term(ai, N, Z)
    real(r_kind), intent(in) :: ai
    integer, intent(in) :: Z, N
    integer :: A
    A = N + Z
    asymmetry_term = ai * (N-Z)**2 / A
end function asymmetry_term

pure real(r_kind) function pairing_term(ap, N, Z)
    real(r_kind), intent(in) :: ap
    integer, intent(in) :: N, Z
    integer :: A
    integer :: even_ctr
    even_ctr = 0
    A = N + Z
    if(mod(N,2) == 0) then
        even_ctr = even_ctr + 1
    end if
    if(mod(Z,2) == 0) then
        even_ctr = even_ctr + 1
    end if
    if (even_ctr == 2) then
        pairing_term = ap * A ** (-3.0/4.0)
    else if (even_ctr == 0) then
        pairing_term = -ap * A ** (-3.0/4.0)
    else
        pairing_term = 0.0
    end if
end function pairing_term

end module micmac