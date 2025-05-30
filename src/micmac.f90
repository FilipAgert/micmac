module micmac

    use constants
    implicit none


    private
    public :: binding_energy, mass_excess, binding_energies, BE_mat, shell_correction, mass_excesses, staircase
    contains

pure function mass_excess(BE, Z, A) result(ME)
    real(r_kind), intent(in) :: BE
    integer, intent(in) :: Z, A
    real(r_kind) :: ME
    integer :: N
    N = A - Z

    ME = BE + (N*mass_n + Z*(mass_p+mass_e)) - dalton*(N+Z)

end function mass_excess    

pure function mass_excesses(parameters, Zs, As) result(MEs)
    real(r_kind), intent(in) :: parameters(num_params)
    integer, intent(in), dimension(:) :: Zs, As
    integer :: num_nuclei
    real(r_kind), dimension(size(Zs)) :: MEs, BEs
    integer :: i

    num_nuclei = size(Zs)
    BEs = binding_energies(parameters, Zs, As, num_nuclei)
    MEs = 0.0_r_kind
    do i = 1, num_nuclei
        MEs(i) = mass_excess(BEs(i), Zs(i), As(i))
    end do
end function


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
    vec = [volume_term(N,Z), surface_term(N,Z), vol_as_term(N,Z), sur_as_term(N,Z),col_term(N,Z),pairing_term(N,Z)]
end function

pure elemental real(r_kind) function volume_term(N, Z)
    integer, intent(in) :: N, Z
    integer :: A
    A = N + Z
    volume_term = A
end function volume_term

pure elemental real(r_kind) function surface_term(N, Z)
    integer, intent(in) :: N, Z
    real(r_kind), parameter :: two_thirds = 2.0_r_kind/3.0_r_kind
    integer :: A
    A = N + Z
    surface_term = A**two_thirds
end function surface_term

pure elemental real(r_kind) function col_term(N, Z)
    integer, intent(in) :: Z, N
    integer :: A
    real(r_kind), parameter :: one_third = 1.0_r_kind/3.0_r_kind
    A = N + Z
    col_term = Z*(Z-1) / A**one_third
end function col_term

pure elemental real(r_kind) function vol_as_term(N,Z)
    integer, intent(in) :: Z, N
    integer :: A
    real(r_kind) :: I
    A = N + Z
    I = (N-Z)*1.0_r_kind/A
    vol_as_term = I*I * A
end function

pure elemental real(r_kind) function sur_as_term(N,Z)
    integer, intent(in) :: Z, N
    integer :: A
    real(r_kind) :: I
    A = N + Z
    I = (N-Z)*1.0_r_kind/A
    sur_as_term = I*I * A**(2.0_r_kind/3.0_r_kind)
end function


pure elemental real(r_kind) function pairing_term(N, Z)
    integer, intent(in) :: N, Z
    integer :: A
    real(r_kind), parameter :: three_fourths = 3.0_r_kind/4.0_r_kind
    real(r_kind) :: p
    A = N + Z
    p = ((-1)**(N) + (-1)**(Z))/2.0_r_kind
    pairing_term = p / sqrt(A*1.0_r_kind)
end function pairing_term



pure elemental real(r_kind) function shell_correction(N,Z)
    integer, intent(in) :: N, Z
    integer :: A
    shell_correction = 0.0_r_kind
end function

pure elemental real(r_kind) function F(N,proton)
    integer, intent(in) :: N !!Nucleon number
    logical,intent(in) :: proton !!proton or neutron

    F = intstaircase(N, proton) - intn23(N)

end function

pure real(r_kind) function intn23(N) !!Integrate from 0 to N n^(2/3)
    integer, intent(in) :: N
    intn23 = real(N,kind=r_kind)**(5.0_r_kind/3.0_r_kind) *3.0_r_kind/5.0_r_kind
end function

!!Integrate staircase function from 0 to n
pure real(r_kind) function intstaircase(N,proton)
    integer, intent(in) :: N !!Nucleon number
    logical,intent(in) :: proton !!proton or neutron
    integer ::ii, nbr, sz, delta
    integer, dimension(:), allocatable :: magics
    real(r_kind) :: height, rolling
    if(proton) then
        sz = size(magic_num_Z)
        allocate(magics(sz))
        magics = magic_num_Z
    else
        sz = size(magic_num_N)
        allocate(magics(sz))
        magics = magic_num_N
    endif
    ii = 2
    rolling=0.0_r_kind
    do while(magics(ii) < N)
        delta = magics(ii) - magics(ii-1)
        height = staircase(real(magics(ii-1),r_kind),magics)
        rolling = rolling + height*delta
        ii = ii + 1
    end do
    !!Okay so now the magic number is greater than N. 
    delta = N - magics(ii-1)
    height = staircase(real(N,kind=r_kind),magics)
    rolling = rolling + height*delta
end function

pure real(r_kind) function staircase(n, magics) !!Gets staircase function at value. If value exactly equals a magic number, use that as the lower bound
    real(r_kind), intent(in) :: n
    integer, intent(in), dimension(:) :: magics
    integer :: num_l, num_gr, ii

    do ii = 1, size(magics)
        if(magics(ii) >= n) then
            num_l = magics(ii-1)
            num_gr = magics(ii)
            exit
        endif
    end do
    staircase = 3.0_r_kind/5.0_r_kind * (num_gr**(5.0_r_kind/3.0_r_kind) - num_l**(5.0_r_kind/3.0_r_kind)) / (num_gr-num_l)
end function

end module micmac