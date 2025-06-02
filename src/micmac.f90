module micmac

    use constants
    use minimise, only: find_min, one_dim
    implicit none


    private
    public :: binding_energy, mass_excess, binding_energies, shell_corr, mass_excesses, staircase, J_mat, find_gs, binding_energy_def
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
    adef = aparam
    alph0sq = alpha0sq(adef, r0, A)
    ! write(*,*) "alpha0: ", sqrt(alph0sq)
    
    N = A - Z
    BE = volume_term(N,Z,av,k) + surface_term(N,Z,as,k)*def_f(def) + col_vol_term(N,Z,r0)*def_g(def) + col_corr_term(N,Z,r0) + pairing_term(N,Z) + shell_corr(N,Z,C, smallC)*exp(-(def**2)/alph0sq)
    ! WRITE(*,*) "Def:", def, ", BE:", BE
    ! write(*,*) volume_term(N,Z,av,k), surface_term(N,Z,as,k)*def_f(def), col_vol_term(N,Z,r0)*def_g(def), col_corr_term(N,Z,r0),pairing_term(N,Z), shell_corr(N,Z,C, smallC)*exp(-(def**2)/alph0sq)
    ! Calculate binding energy
end function binding_energy_def

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

subroutine find_gs(BE, def, params, Z, A)
    real(r_kind), intent(out) :: BE, def
    real(r_kind), intent(in) :: params(num_params)
    integer, intent(in) :: Z,A
    integer :: Niters
    procedure(one_dim), pointer  :: func_def
    func_def => be_def_func(params, Z, A)

    call find_min(def, BE, Niters, func_def, -1.0_r_kind, 1.0_r_kind,0.2_r_kind)
    ! write(*,'(A,I5, A)') "Found ground state after ", Niters, " iterations"
    ! write(*,'(A, f10.3, A, f10.3, A)') "Deformation: ", def, ", binding energy: ", BE, " MeV"
end subroutine

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
    adef = aparam
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


pure function alpha0sq(adef, r0, A)
    real(r_kind), intent(in) :: adef, r0
    integer, intent(in) :: A
    real(r_kind) :: alpha0sq

    alpha0sq = 5.0_r_kind * (adef/r0)**2 * A**(-2.0_r_kind/3.0_r_kind)

end function

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
    surface_term = as*(1-k*I*I/(1.0_r_kind))*A**two_thirds
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



pure function J_mat(Zs, As, num_nuclei, params)
    implicit none
    integer, intent(in), dimension(num_nuclei) :: As, Zs
    integer,intent(in) :: num_nuclei
    real(r_kind), intent(in) :: params(num_params)
    integer :: i, Z,A
    real(r_kind) :: J_mat(num_nuclei,num_params)
    do i = 1,num_nuclei
        Z = Zs(i)
        A = As(i)
        J_mat(i,:) = J_vec(Z,A,params)
    end do
end function

pure function J_vec(Z,A,params)
    integer, intent(in) :: Z,A
    integer :: N
    real(r_kind),intent(in) :: params(num_params)
    real(r_kind) :: J_vec(num_params)
    real(r_kind) :: av, as, r0, k, smallC, C
    av = params(1)
    as = params(2)
    k = params(3)
    r0 = params(4)
    C = params(5)
    smallC = params(6)
    N = A - Z
    J_vec = [dBe_daV(N,Z,k),dBe_daS(N,Z,k),dBe_dk(N,Z,av,as),dBe_dr0(N,Z,r0),dBe_dC(N,Z,smallC),dBe_dsc(N,Z,C)]
end function

function J_shell_part(As,Zs,params, num_nuclei)
    integer, intent(in), dimension(num_nuclei) :: As, Zs
    real(r_kind), intent(in) :: params(num_params)
    integer,intent(in) :: num_nuclei
    integer :: i, Ns(num_nuclei)
    real(r_kind) :: J_shell_part(num_nuclei,2), C, smallC
    Ns = As-Zs
    C = params(num_params-1)
    smallC = params(num_params)
    J_shell_part = 0.0_r_kind
    do i = 1, num_nuclei
        J_shell_part(i,1) = dBe_dC(Ns(i), Zs(i), smallC)
        J_shell_part(i,2) = dBe_dsc(Ns(i), Zs(i), C)
    end do
    
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

pure real(r_kind) function dBe_daS(N,Z,k)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: k
    real(r_kind) :: I
    integer :: A
    A = N + Z
    I = (N-Z)*1.0_r_kind/A
    dBe_daS = 1.0_r_kind-k * I*I
    dBe_daS = dBe_daS * A**(2.0_r_kind/3.0_r_kind)
end function

pure real(r_kind) function dBe_dk(N,Z,av,as)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: av,as
    real(r_kind) :: I
    integer :: A
    A = N + Z
    I = (N-Z)*1.0_r_kind/A
    dBe_dk = -I*I*(av*A+as*A**(2.0_r_kind/3.0_r_kind))
end function

pure real(r_kind) function dBe_dr0(N,Z,r0)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: r0
    integer :: A
    A = N + Z
    dBe_dr0 = 3.0_r_kind*e_squared*Z*Z*(5.0_r_kind*pi2*d*d-2*A**(2.0_r_kind/3.0_r_kind)*r0*r0)/(10.0_r_kind*A*r0**4)
end function dBe_dr0
!!Partial derivate with respect to small C
pure elemental real(r_kind) function dBe_dsc(N,Z,C)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: C
    integer :: A
    A = Z + N
    dBe_dsc = -C* A**(1.0_r_kind/3.0_r_kind)
end function dBe_dsc

!!Partial derivate of shell correction with respect to C
pure real(r_kind) elemental function dBe_dC(N,Z,smallC)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: smallC
    integer :: A
    A = Z + N
    dBe_dC = (F(Z,.true.) + F(N,.false.)) / ((A*1.0_r_kind/2.0_r_kind)**(2.0_r_kind/3.0_r_kind)) - smallC * A**(1.0_r_kind/3.0_r_kind)
end function dBe_dC
!! ########################## Shell correction Part
pure elemental real(r_kind) function shell_corr(N,Z, C, smallC)
    integer, intent(in) :: N, Z
    real(r_kind), intent(in) :: C, smallC
    shell_corr = dBe_dC(N,Z,smallC) * C
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

end module micmac