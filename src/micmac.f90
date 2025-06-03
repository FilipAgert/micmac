module micmac

    use constants
    use minimise, only: find_min, one_dim, rand_d, find_min_brute_force
    implicit none


    private
    public :: binding_energy, mass_excess, binding_energies, shell_corr, mass_excesses, staircase, J_mat, find_gs, binding_energy_def, find_gs_multiple, alpha_to_beta
    contains

pure function mass_excess(BE, Z, A) result(ME)
    real(r_kind), intent(in) :: BE
    integer, intent(in) :: Z, A
    real(r_kind) :: ME
    integer :: N
    N = A - Z

    ME = BE + (N*mass_n + Z*(mass_p+mass_e)) - dalton*(N+Z)

end function mass_excess

function mass_excesses(parameters, Zs, As) result(MEs)
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

pure function binding_energy(params, Z, A) result(BE)
    ! This function calculates the binding energy of a nucleus
    real(r_kind), intent(in) :: params(num_params)
    integer, intent(in) :: Z, A
    real(r_kind) :: BE
    integer :: N
    real(r_kind) :: av, as, r0, k, smallC, C
    av = params(1)
    as = params(2)
    k = params(3)
    r0 = params(4)
    C = params(5)
    smallC = params(6)
    N = A - Z
    BE = volume_term(N,Z,av,k) + surface_term(N,Z,as,k) + col_term(N,Z,r0) + pairing_term(N,Z) + shell_corr(N,Z,C, smallC)  
    
    ! Calculate binding energy
end function binding_energy

!!Gets binding energy for a given deformtaion
function binding_energy_def(params, Z, A, def) result(BE)
    ! This function calculates the binding energy of a nucleus
    real(r_kind), intent(in) :: params(num_params), def
    integer, intent(in) :: Z, A
    real(r_kind) :: BE
    
    integer :: N
    real(r_kind) :: av, as, r0, k, smallC, C, adef, alph0sq
    av = params(1)
    as = params(2)
    k = params(3)
    r0 = params(4)
    C = params(5)
    smallC = params(6)
    adef = params(7)
    alph0sq = alpha0sq(adef, r0, A)
    
    N = A - Z
    BE = volume_term(N,Z,av,k) + surface_term(N,Z,as,k)*def_f(def) + &
        col_vol_term(N,Z,r0)*def_g(def) + col_corr_term(N,Z,r0) + &
        pairing_term(N,Z) + shell_corr(N,Z,C, smallC)*shell_damp_fac(N,Z,def,adef,r0)
    ! WRITE(*,*) "Def:", def, ", BE:", BE
    ! write(*,*) volume_term(N,Z,av,k), surface_term(N,Z,as,k)*def_f(def), col_vol_term(N,Z,r0)*def_g(def), col_corr_term(N,Z,r0),pairing_term(N,Z), shell_corr(N,Z,C, smallC)*exp(-(def**2)/alph0sq)
    ! Calculate binding energy
end function binding_energy_def




!!For a vector of Zs and As, compute binding energy.
pure function binding_energies(parameters, Zs, As, num_nuclei) result(BEs)
    real(r_kind), intent(in) :: parameters(num_params)
    integer, intent(in), dimension(num_nuclei) :: Zs, As
    integer, intent(in) :: num_nuclei
    real(r_kind), dimension(num_nuclei) :: BEs
    real(r_kind), dimension(num_nuclei,num_params) :: mat_format
    integer :: i
    do i = 1, num_nuclei
        BEs(i) = binding_energy(parameters, Zs(i), As(i))
    end do
end function

!!Computes ground state energy and deformation for given nuclei.
subroutine find_gs_multiple(BEs, defs, params, Zs, As, num_nuclei)
    real(r_kind), intent(out) :: BEs(num_nuclei), defs(num_nuclei)
    real(r_kind), intent(in) :: params(num_params)
    integer, intent(in), dimension(num_nuclei) :: Zs, As
    integer, intent(in) :: num_nuclei
    real(r_kind) :: BE, def
    integer :: ii, Z , A

    do ii = 1, num_nuclei
        Z = Zs(ii)
        A = As(ii)
        call find_gs(BE, def, params, Z, A)
        BEs(ii) = BE
        defs(ii) = def
        ! WRITE(*,*) "Z,   A"
        ! write(*,'(2I4)') Z, A
        ! write(*,'(A)') "BE (MeV),         def"
        ! write(*,'(E15.3,2x, F10.3)') BE, def

        ! if(ii > num_nuclei - 1) then
        !     call exit
        ! endif
    end do
end subroutine

!!Minimizes BE wrt deformation and returns G.S. energy and G.S. deformation
subroutine find_gs(BE, def, params, Z, A)
    real(r_kind), intent(out) :: BE, def
    real(r_kind), intent(in) :: params(num_params)
    integer, intent(in) :: Z,A
    procedure(one_dim), pointer  :: func_def
    func_def => be_def_func(params, Z, A) !!Gets binding energy as a function of deformation

    !call find_min_brute_force(def, BE, func_def, -default_def_bounds, default_def_bounds)
    call find_min(def, BE, func_def, -default_def_bounds, default_def_bounds,default_num_restarts) !!Minimize binding energy as a function of deformation
    ! write(*,'(A,I5, A)') "Found ground state after ", Niters, " iterations"
    ! write(*,'(A, f10.3, A, f10.3, A)') "Deformation: ", def, ", binding energy: ", BE, " MeV"
end subroutine

!!Function takes the binding energy formula and creates a function on the form A + Bx^2 + Cx^3 + Dexp(-x^2/a^2)
!!Where x is the deformation. This function is then minimzed to find the deformation x and G.S. energy.
function be_def_func(parameters, Z, A) result(func)
    procedure(one_dim), pointer :: func
    real(r_kind), intent(in) :: parameters(:)
    integer, intent(in) :: Z, A

    ! Local variables
    real(r_kind) :: av, as, k, r0, C, smallC, adef
    real(r_kind) :: constant, quadratic, cubic, exponent, expconst
    real(r_kind) :: surf, colvol
    integer :: N

    ! Internal function pointer
    procedure(one_dim), pointer :: f_ptr

    ! Assign parameter values
    av = parameters(1)
    as = parameters(2)
    k = parameters(3)
    r0 = parameters(4)
    C = parameters(5)
    smallC = parameters(6)
    adef = parameters(7)
    N = A - Z

    ! Compute components of the function
    surf     = surface_term(N, Z, as, k)
    colvol   = col_vol_term(N, Z, r0)
    constant = volume_term(N,Z,av,k) + surf + col_term(N,Z,r0) + pairing_term(N,Z)
    quadratic = 2.0_r_kind / 5.0_r_kind * surf - colvol / 5.0_r_kind
    cubic     = -(surf + colvol) * 4.0_r_kind / 105.0_r_kind
    expconst  = shell_corr(N, Z, C, smallC)
    exponent  = -1.0_r_kind / (alpha0sq(adef, r0, A))

    ! write(*,*) "const", constant
    ! write(*,*) "quad", quadratic
    ! write(*,*) "cubic", cubic
    ! write(*,*) "exp", exponent
    ! write(*,*) "shellcorr", expconst
    f_ptr => make_func(constant, quadratic, cubic, exponent, expconst)
    func => f_ptr
end function be_def_func


function make_func(constant, quadratic, cubic, exponent, expconst) result(f_ptr)
    procedure(one_dim), pointer :: f_ptr
    real(r_kind), intent(in) :: constant, quadratic, cubic, exponent, expconst

    procedure(one_dim), pointer :: f
    f => internal_f
    f_ptr => f
    

    contains

    function internal_f(x) result(res)
        real(r_kind), intent(in) :: x
        real(r_kind) :: res
        res = constant + x**2 * quadratic + x**3 * cubic + expconst * exp(exponent * x**2)
    end function internal_f

end function make_func



pure elemental real(r_kind) function volume_term(N, Z,av,k)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: av,k
    real(r_kind), parameter :: two_thirds = 2.0_r_kind/3.0_r_kind
    real(r_kind) :: I
    integer :: A
    A = N + Z
    I = (N-Z)*1.0_r_kind/real(A,kind=r_kind)
    volume_term = av*(1-k*I*I)*A
end function volume_term

pure elemental real(r_kind) function surface_term(N, Z,as,k)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: as,k
    real(r_kind), parameter :: two_thirds = 2.0_r_kind/3.0_r_kind
    real(r_kind) :: I
    integer :: A
    A = N + Z
    I = (N-Z)*1.0_r_kind/real(A,kind=r_kind)
    surface_term = as*(1-k*I*I)*A**two_thirds
end function surface_term

pure elemental real(r_kind) function col_term(N, Z, r0)
    implicit none
    integer, intent(in) :: Z, N
    real(r_kind), intent(in) :: r0
    integer :: A
    real(r_kind), parameter :: one_third = 1.0_r_kind/3.0_r_kind
    A = N + Z
    col_term = col_vol_term(N,Z,r0) + col_corr_term(N,Z,r0)
end function col_term

pure elemental real(r_kind) function col_vol_term(N,Z,r0)
    implicit none
    integer, intent(in) :: Z, N
    real(r_kind), intent(in) :: r0
    integer :: A
    real(r_kind), parameter :: one_third = 1.0_r_kind/3.0_r_kind
    A = N + Z
    col_vol_term = Z*Z *1.0_r_kind / (A**(one_third)) * 3.0_r_kind * e_squared / (5.0_r_kind * r0)


end function

pure elemental real(r_kind) function col_corr_term(N,Z,r0)
    implicit none
    integer, intent(in) :: Z, N
    real(r_kind), intent(in) :: r0
    integer :: A
    real(r_kind), parameter :: one_third = 1.0_r_kind/3.0_r_kind
    A = N + Z
    col_corr_term = - (Z * Z * 1.0_r_kind / A) * pi2/2 * (d/r0)**2 * e_squared/r0 

end function

pure elemental real(r_kind) function pairing_term(N, Z)
    implicit none
    integer, intent(in) :: N, Z
    integer :: A
    real(r_kind), parameter :: P = 11
    A = N + Z
    if(mod(A,2) == 1) then
        pairing_term = 0.0_r_kind
    else if(mod(N,2) == 1 .and. mod(Z,2) == 1) then
        pairing_term = P
    else
        pairing_term = -P
    end if
    pairing_term = pairing_term / sqrt(A*1.0_r_kind)
end function pairing_term

pure elemental real(r_kind) function shell_damp_fac(N,Z, def, adef, r0) result(fac)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: def, adef, r0
    real(r_kind) :: alph0sq
    integer :: A
    A = N + Z
    alph0sq = alpha0sq(adef, r0, A)
    fac = exp(-def**2/alph0sq)
end function


pure function def_f(alpha)
    real(r_kind), intent(in) :: alpha
    real(r_kind) :: def_f
    def_f = 1+(alpha**2)*2.0_r_kind/5.0_r_kind - alpha**3 * 4.0_r_kind/105.0_r_kind

end function

pure function def_g(alpha)
    real(r_kind), intent(in) :: alpha
    real(r_kind) :: def_g
    def_g = 1-(alpha**2)/5.0_r_kind - alpha**3 * 4.0_r_kind/105.0_r_kind

end function


pure function alpha0sq(adef, r0, A)
    real(r_kind), intent(in) :: adef, r0
    integer, intent(in) :: A
    real(r_kind) :: alpha0sq

    alpha0sq = 5.0_r_kind * (adef/r0)**2 * A**(-2.0_r_kind/3.0_r_kind)

end function

 function J_mat(Zs, As, defs, num_nuclei, params)
    implicit none
    integer, intent(in), dimension(num_nuclei) :: As, Zs
    integer,intent(in) :: num_nuclei
    real(r_kind), intent(in) :: params(num_params), defs(num_nuclei)
    integer :: i, Z,A
    real(r_kind) :: J_mat(num_nuclei,num_params), def
    do i = 1,num_nuclei
        Z = Zs(i)
        A = As(i)
        def = defs(i)
        J_mat(i,:) = J_vec(Z,A,params, def)
    end do
end function

 function J_vec(Z,A,params, def)
    integer, intent(in) :: Z,A
    integer :: N
    real(r_kind),intent(in) :: params(num_params), def
    real(r_kind) :: J_vec(num_params)
    real(r_kind) :: av, as, r0, k, smallC, C, adef
    av = params(1)
    as = params(2)
    k = params(3)
    r0 = params(4)
    C = params(5)
    smallC = params(6)
    adef = params(7)
    N = A - Z
    J_vec = [dBe_daV(N,Z,k),dBe_daS(N,Z,k, def),dBe_dk(N,Z,av,as, def),dBe_dr0(N,Z,r0, def, C, smallC, adef),&
            dBe_dC(N,Z,smallC, def, adef, r0),dBe_dsc(N,Z,C, def,adef,r0), &
            dBe_dadef(N,Z,C,smallC,def,adef, r0)]
    if(isnan(J_vec(1)) .or. J_vec(1) > 1e6) then
        call exit
    endif
    ! write(*,*)
    ! write(*,'(A,2I4)') "Z, A: ", Z, A
    ! write(*,'(A, F10.3)') 'Deformation: ', def
    ! write(*,'(7F10.3)') J_vec
end function


pure real(r_kind) function dBe_daV(N,Z,k)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: k
    real(r_kind) :: I
    integer :: A
    A = N + Z
    I = (N-Z)*1.0_r_kind/A
    dBe_daV = 1.0_r_kind-k * I*I
    dBe_daV = dBe_daV * A
end function

pure real(r_kind) function dBe_daS(N,Z,k, def)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: k, def
    real(r_kind) :: I
    integer :: A
    A = N + Z
    I = (N-Z)*1.0_r_kind/A
    dBe_daS = 1.0_r_kind-k * I*I
    dBe_daS = dBe_daS * A**(2.0_r_kind/3.0_r_kind) * def_f(def)
end function

pure real(r_kind) function dBe_dk(N,Z,av,as, def)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: av,as, def
    real(r_kind) :: I
    integer :: A
    A = N + Z
    I = (N-Z)*1.0_r_kind/A
    dBe_dk = -I*I*(av*A+def_f(def)*as*A**(2.0_r_kind/3.0_r_kind))
end function

pure real(r_kind) function dBe_dr0(N,Z,r0, def, C, smallC, adef)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: r0, def, c, smallC, adef
    integer :: A
    A = N + Z
    !!Coulomb terms where one term has macroscopic deformation 
    dBe_dr0 = 3.0_r_kind*e_squared*Z*Z*(5.0_r_kind*pi2*d*d-2*A**(2.0_r_kind/3.0_r_kind)*r0*r0*def_g(def))/(10.0_r_kind*A*r0**4)


    !!Shell correction term as alpha0 has r0
    dBe_dr0 = dBe_dr0 - 2.0_r_kind/5.0_r_kind * r0*(def/adef)**2 * A**(2.0_r_kind/3.0_r_kind) * shell_corr(N,Z,C,smallC)*shell_damp_fac(N,Z,def,adef,r0) 
end function dBe_dr0
!!Partial derivate with respect to small C
pure elemental real(r_kind) function dBe_dsc(N,Z,C, def, adef, r0)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: C, def, adef, r0
    integer :: A
    A = Z + N
    dBe_dsc = -C* A**(1.0_r_kind/3.0_r_kind) * shell_damp_fac(N,Z,def,adef,r0)
end function dBe_dsc

!!Partial derivate of shell correction with respect to C
real(r_kind) function dBe_dC(N,Z,smallC, def, adef, r0)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: smallC, def, adef, r0
    integer :: A
    A = Z + N
    dBe_dC = (F(Z,.true.) + F(N,.false.)) / ((A*1.0_r_kind/2.0_r_kind)**(2.0_r_kind/3.0_r_kind)) - smallC * A**(1.0_r_kind/3.0_r_kind)
    dBe_dC = dBe_dC * shell_damp_fac(N,Z,def,adef,r0)

end function dBe_dC

!!Partial derivate of shell correction with respect to C
pure real(r_kind) elemental function dBe_dadef(N,Z,C, smallC, def, adef, r0)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: smallC, def, adef, r0, C
    integer :: A
    A = Z + N
    dBe_dadef = 2.0_r_kind/5.0_r_kind * A**(2.0_r_kind/3.0_r_kind) * def**2 * r0**2/ (adef**3) * shell_corr(N,Z,C,smallC) * shell_damp_fac(N,Z,def,adef,r0)
end function dBe_dadef

!! ########################## Shell correction Part
pure elemental real(r_kind) function shell_corr(N,Z, C, smallC)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: C, smallC
    integer :: A
    A = N + Z
    shell_corr = (F(Z,.true.) + F(N,.false.)) / ((A*1.0_r_kind/2.0_r_kind)**(2.0_r_kind/3.0_r_kind)) - smallC * A**(1.0_r_kind/3.0_r_kind)
    shell_corr = shell_corr * C
end function




pure real(r_kind) function F(N,proton)
    integer, intent(in) :: N !!Nucleon number
    logical,intent(in) :: proton !!proton or neutron
    real(r_kind) ::fold
    integer, allocatable :: magics(:)
    integer ::sz, ii, magic
    
    if(proton) then
        sz = size(magic_num_Z)
        allocate(magics(sz))
        magics = magic_num_Z
    else
        sz = size(magic_num_N)
        allocate(magics(sz))
        magics = magic_num_N
    endif
    do ii = 1,sz
        magic = magics(ii)
        if(magic > N) then
            magic = magics(ii-1)
            exit
        endif
    end do
    F = staircase(real(N,r_kind), magics)*(N-magic) - 3.0_r_kind/5.0_r_kind * (N**(5.0_r_kind/3.0_r_kind) - magic**(5.0_r_kind/3.0_r_kind))
end function

pure real(r_kind) function staircase(n, magics) !!Gets staircase function at value. If value exactly equals a magic number, use that as the lower bound
    real(r_kind), intent(in) :: n
    integer, intent(in), dimension(:) :: magics
    integer :: num_l, num_gr, ii

    do ii = 2, size(magics)
        if(magics(ii) >= n) then
            num_l = magics(ii-1)
            num_gr = magics(ii)
            exit
        endif
    end do
    staircase = 3.0_r_kind/5.0_r_kind * (num_gr**(5.0_r_kind/3.0_r_kind) - num_l**(5.0_r_kind/3.0_r_kind)) / (num_gr-num_l)
end function


pure real(r_kind) function alpha_to_beta(alpha) result(beta)
    real(r_kind), intent(in) :: alpha
    beta = sqrt(5/(4*pi)) * alpha

end function

end module micmac