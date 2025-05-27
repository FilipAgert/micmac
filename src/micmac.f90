module micmac

    use constants
    implicit none


    private
    public :: binding_energy, mass_excess, binding_energies, BE_mat
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
    real(r_kind), dimension(num_params) :: vec_format

    vec_format = BE_vec(Z,A)
    BE = dot_product(vec_format,parameters)
    
    ! Calculate binding energy
end function binding_energy


!!For a vector of Zs and As, compute binding energy.
pure function binding_energies(parameters, Zs, As, num_nuclei) result(BEs)
    real(r_kind), intent(in) :: parameters(num_params)
    integer, intent(in), dimension(num_nuclei) :: Zs, As
    integer, intent(in) :: num_nuclei
    real(r_kind), dimension(num_nuclei) :: BEs
    real(r_kind), dimension(num_nuclei,num_params) :: mat_format

    mat_format = BE_mat(Zs,As,num_nuclei)
    BEs = MATMUL(mat_format,parameters)
end function

!!Get the Binding energy matrix X. If multiplied by parameters y = Xb, y will be filled with predictions of binding energy
pure function BE_mat(Zs, As, num_nuclei) result(mat)
    integer, intent(in), dimension(num_nuclei) :: Zs, As
    integer, intent(in) :: num_nuclei
    real(r_kind), dimension(num_nuclei, num_params) :: mat
    integer :: i
    do i = 1, num_nuclei
        mat(i,:) = BE_vec(Zs(i), As(i))
    end do
end function

pure function BE_vec(Z,A) result(vec)
    integer, intent(in) :: Z,A
    real(r_kind), dimension(num_params) :: vec
    integer :: N
    N = A-Z
    vec = [volume_term(N,Z), surface_term(N,Z), coulomb_term(N,Z), asymmetry_term(N,Z), -pairing_term(N,Z)]
end function

pure elemental real(r_kind) function volume_term(N, Z)
    integer, intent(in) :: N, Z
    integer :: A
    A = N + Z
    volume_term = A
end function volume_term

pure elemental real(r_kind) function surface_term(N, Z)
    integer, intent(in) :: N, Z
    real(r_kind), parameter :: two_thirds = 2.0/3.0
    integer :: A
    A = N + Z
    surface_term = A**two_thirds
end function surface_term

pure elemental real(r_kind) function coulomb_term(N, Z)
    integer, intent(in) :: Z, N
    integer :: A
    real(r_kind), parameter :: one_third = 1.0/3.0
    A = N + Z
    coulomb_term = Z*Z / A**one_third
end function coulomb_term

pure elemental real(r_kind) function asymmetry_term(N, Z)
    integer, intent(in) :: Z, N
    integer :: A
    A = N + Z
    asymmetry_term = (N-Z)**2 / A
end function asymmetry_term

pure elemental real(r_kind) function pairing_term(N, Z)
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
        pairing_term = A ** (-3.0/4.0)
    else if (even_ctr == 0) then
        pairing_term = -A ** (-3.0/4.0)
    else
        pairing_term = 0.0
    end if
end function pairing_term

end module micmac