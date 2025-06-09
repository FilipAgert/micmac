module def_ho
    use constants, only: r_kind, hbarc, pi
    use micmac, only: deformation
    use optimise, only:func_1d, conj_grad_method
    use quadrule, only:laguerre_ss_compute, hermite_ek_compute
    implicit none
    private
    public :: an_ho_state, get_ho_states, getnumstates, betadef, dist_min, getnumstatesupto, Ws_mat_elem

    integer, parameter :: gauss_order =64
    real(r_kind), dimension(gauss_order) :: her_x, her_w, lag_x, lag_w
    logical :: precomputed = .false.

    real(r_kind), parameter :: aws = 0.70_r_kind !! (fm)
    type :: betadef
        real(r_kind) :: beta2, beta4

        contains
        procedure :: omega_perp => compute_omega_perp
        procedure :: omega_z => compute_omega_Z
        procedure :: omega_def => compute_omega_def

    end type

    type :: ax_deformed_ho !!
        real(r_kind) :: freq(3) 
        type(betadef) :: def

    end type


    type, extends(func_1d) :: dist_min
        real(r_kind) :: r, theta, radius! radius normalized as r/R0
        type(betadef) :: def
    contains
        procedure :: eval => surfdisteval

    end type


    type :: an_ho_state !!mj is half integer. so mj = 3 in reality means 3/2
        integer :: nz, nr, ml, mj
        contains
        procedure :: text => aniho_text
        procedure :: header => aniho_header
        procedure :: kinetic_energy => ho_kinen
    end type
    !! Computes single particle energies in a shell model potential.
    contains 


pure function compute_omega_def(self, omega0) result(omegadef)
    class(betadef), intent(in) :: self
    real(r_kind), intent(in) :: omega0
    real(r_kind) :: omegadef
    real(r_kind) :: delta
    delta = self%beta2 / 1.057
    omegadef = omega0 * (1.0_r_kind + 2.0_r_kind/3.0_r_kind * delta**2)
end function

pure function compute_omega_perp(self, omega0) result(omegaperp)
    class(betadef), intent(in) :: self
    real(r_kind), intent(in) :: omega0
    real(r_kind) :: omegaperp
    real(r_kind) :: delta
    delta = self%beta2 / 1.057
    omegaperp = self%omega_def(omega0)*sqrt((1.0_r_kind + 2.0_r_kind/3.0_r_kind * delta))
end function

pure function compute_omega_Z(self, omega0) result(omegaz)
    class(betadef), intent(in) :: self
    real(r_kind), intent(in) :: omega0
    real(r_kind) :: omegaz
    real(r_kind) :: delta
    delta = self%beta2 / 1.057
    omegaz = self%omega_def(omega0) * sqrt( (1.0_r_kind - 4.0_r_kind/3.0_r_kind * delta))
end function

pure function ho_kinen(self, hbaromega_z, hbaromega_xy) result(en)
    class(an_ho_state), intent(in) :: self
    real(r_kind) :: en
    real(r_kind), intent(in) :: hbaromega_z, hbaromega_xy

    en = ani_ho_en(self, hbaromega_z, hbaromega_xy) / 2.0_r_kind !as per the virial theorem

end function

function aniho_text(self) result(text)
    class(an_ho_state) :: self
    character(len=100) :: text
    integer :: N
    N = self%ml + self%nr * 2 + self%nz
    write(text, '(5I5,A)') N, self%nz, self%nr, self%ml, self%mj, "/2"
end function

function aniho_header(self) result(hdr)
    class(an_ho_state) :: self
    character(len=100) :: hdr
    write(hdr, '(A)') "    N    nz   nr   ml    mj "

end function
subroutine precompute_weights()
    if(precomputed) then
        return
    endif
    precomputed = .true.

    !!64 point gaussian quadrature.
    call hermite_ek_compute(gauss_order,her_x,her_w)
    !!Gets weights and locations for where to evaluate 64 point integral. 
    !!Integral must be of form 
    !! -infty < x < infty  dx f(x) * exp(-x^2)
    !! where x = alpha_z*z => z = x/alpha_z


    call laguerre_ss_compute(gauss_order, lag_x, lag_w)
    !!Gets weights and locations for where to evaluate 64 point integral. 
    !!Integral must be of form 
    !! 0 < x < infty  dx f(x) * exp(-x)
    !! where x = alpha^2 * rho^2


end subroutine


pure elemental real(r_kind) function Y20(theta)
    real(r_kind), intent(in) :: theta
    Y20 = sqrt(5.0_r_kind / (16 * pi)) * (3.0_r_kind * cos(theta)**2 - 1.0_r_kind)
end function

pure elemental real(r_kind) function Y40(theta)
    real(r_kind) ,intent(in) :: theta
    Y40 = sqrt(9.0_r_kind / (256 * pi)) * (35.0_r_kind * cos(theta)**4 - 30.0_r_kind * cos(theta)**2 + 3)

end function

pure elemental real(r_kind) function surfdisteval(self, x) !!Gets distance to surface at angle x.
    class(dist_min), intent(in) :: self
    real(r_kind), intent(in) :: x !!Angle theta in this case.

    real(r_kind) :: radius
    radius = surfRadius(x, self%def, self%radius)

    surfdisteval = (self%r*cos(self%theta) - radius*cos(x))**2 &
                 + (self%r*sin(self%theta) - radius*sin(x))**2
end function

pure elemental real(r_kind) function surfRadius(theta, def, R0) !!computes distance to surface
    real(r_kind), intent(in) :: theta, R0
    type(betadef), intent(in) :: def
    real(r_kind) :: dist
    surfRadius = R0 * (1.0_r_kind + def%beta2 * Y20(theta) + def%beta4*Y40(theta))
end function


real(r_kind) function Ws_mat_elem(state1, state2, def, Ws_depth, radius, mass, omegaz, omegaperp) result(elem) 
    !!Calculates matrix element for central part of WS potential.
    !!uses 64 points gaussian quadrature
    type(an_ho_state), intent(in) :: state1, state2
    type(betadef) :: def
    real(r_kind), intent(in) :: Ws_depth, radius
    real(r_kind), intent(in) :: mass !!In units MeV/c^2
    real(r_kind), intent(in) :: omegaperp, omegaz !!In units MeV/hbar
    integer :: il, ih
    real(r_kind) :: rolling, z, rho, alphaz, alpharho, u, x, factz, factrho, xpart, r, ang
    !The phi integral turns out to be a kronecker delta * 2pi
    if (state1%ml /= state2%ml) then 
        elem = 0
        return
    else if(state1%mj /= state2%mj) then !Spin.
        elem = 0
        return
    endif
    elem = 2*pi !!from phi part.

    call precompute_weights()

    
    alphaz = sqrt(mass*omegaz/(hbarc*hbarc))!Replace with sqrt(m * omega/hbar)
    alpharho = sqrt(mass*omegaperp/(hbarc*hbarc))

    factz = 1.0_r_kind/sqrt(pi * 2**(state1%nz+state2%nz) * fac(state1%nz)*fac(state2%nz))
    factrho = 1.0_r_kind/(pi * sqrt(real(fac(state1%nr) * fac(state2%nr),r_kind) / real(fac(state1%nr+state1%ml) * fac(state2%nr+state2%ml),r_kind)))


    rolling = 0.0_r_kind
    do il = 1, gauss_order
        x = lag_x(il)
        rho = sqrt(x)/alpharho
        xpart = x ** ((state1%ml + state2%ml)/2) * lna(x,state1%nr,state1%ml) * lna(x,state2%nr,state2%ml)
        do ih = 1, gauss_order
            u = her_x(ih)
            z = u/alphaz
            r = rad_cyl(rho, z)
            ang = theta_cyl(rho,z)

            rolling = rolling + lag_w(il) * her_w(ih) * Hn(u,state1%nr) * Hn(u, state2%nr) & !Z part
                    * xpart * Ws(r,ang,def, Ws_depth, radius)!!Xpart and mixed part.
            
            ! if(rolling < -1e6) then
            !     write(*,*) rolling, lag_w(il), her_w(ih), u, z, r, ang, xpart
            !     pause 'press enter to continue'
            ! endif
        end do
    end do

    
    elem = elem * factz * factrho * rolling

end function

real(r_kind) pure elemental function theta_cyl(rho,z)  !!Coordinate transform from cylindrical coordinates to spherical theta
    real(r_kind), intent(in):: rho,z
    theta_cyl = ATAN(z/rho)
end function

real(r_kind) pure elemental function rad_cyl(rho,z)  !!Coordinate transform from cylindrical coordinates to spherical radius
    real(r_kind), intent(in):: rho,z
    rad_cyl = sqrt(rho*rho+z*z)
end function

real(r_kind) function Ws(r, theta, def, ws_depth, radius) !!axially symmetric. function of z and rho
    real(r_kind), intent(in) :: r, theta, ws_depth, radius
    type(betadef), intent(in) :: def
    real(r_kind) ::dist
    dist = sqrt(surfdist(r,theta,def, radius)) 

    Ws = -ws_depth/ (1 + exp(dist/aws))
end function

real(r_kind) function surfdist(r,theta,def, radius) result(dist) !!Finds shortest distance to surface of nucleus.
    real(r_kind), intent(in) :: r, theta, radius
    type(betadef), intent(in) :: def
    type(dist_min) :: distfunc
    real(r_kind) :: ang
    logical :: conv

    distfunc = dist_min(r=r, theta=theta, def=def, radius=radius)
    call conj_grad_method(ang, dist, conv, distfunc, 0.0_r_kind, pi,theta,1e-5_r_kind)
end function

pure elemental real(r_kind) function phiVcphi(state1, state2, def)
    type(betadef), intent(in) :: def
    type(an_ho_state), intent(in) :: state1, state2 !!Compute a matrix element for coulomb repulsion
    
    !Computes triple integral 1/(abs)

end function

pure elemental integer function getnumstates(N) result(num) !!get number of states in one major shell of h.o
    integer, intent(in) :: N
    num = (N+1)*(N+2)/2 !!incl spin degeneracy
end function

pure elemental integer function getnumstatesupto(N) result(num) !!get number of states in all major shells up to N
    integer, intent(in) :: N
    num = (N+1)*(N+2)*(N+3)/6 !!incl spin degeneracy
end function

pure function get_ho_states(N) result(states)
    integer, intent(in) :: N !!ho quantum number
    type(an_ho_state), dimension(getnumstates(N)) :: states
    type(an_ho_state) :: state
    integer :: nz, nr, ml, omega, qrem, idx
    idx = 0
    do nz = 0,N
        !we need nz + 2nr + ml = N.
        qrem = N - nz
        do nr = 0, qrem/2 !!integer division. if remaining is 3, we want n = 0, 1
            ml = qrem - 2*nr
            do omega = abs(2*ml - 1), 2*ml + 1, 2 
                idx = idx + 1
                state = an_ho_state(nz = nz, nr = nr, ml=ml, mj=omega)
                states(idx) = state
            end do
        end do
    end do
end function



pure function pl(x,l)
    real(r_kind), intent(in) :: x
    integer, intent(in) :: l
    real(r_kind) :: pl

    if (l == 0) then
        pl = 1.0_r_kind
    else if (l == 1) then
        pl = x
    else if (l == 2) then
        pl = 0.5_r_kind * (3.0_r_kind * x**2 - 1.0_r_kind)
    else if (l == 3) then
        pl = 0.5_r_kind * (5.0_r_kind * x**3 - 3.0_r_kind * x)
    else if(l == 4) then
        pl = 0.125_r_kind * (35.0_r_kind * x**4 - 30.0_r_kind * x**2 + 3.0_r_kind)
    else
        error stop "Higher order Legendre polynomials not implemented"
    end if

end function

pure real(r_kind) elemental recursive function Hn(x,n) result(res)
    real(r_kind), intent(in) :: x
    integer, intent(in) :: n
    if(n == 0) then
        res = 1.0_r_kind
    else if(n==1) then
        res = 2.0_r_kind * x
    else
        res = 2.0_r_kind * x * Hn(x,n-1) - 2.0_r_kind * (n-1) * Hn(x,n-2)
    endif
end function

pure real(r_kind) elemental function ani_ho_en(state,hbaromega_z,hbaromega_xy) result(energy)
    real(r_kind), intent(in) :: hbaromega_z, hbaromega_xy
    class(an_ho_state), intent(in) :: state
    energy = hbaromega_xy*(2*state%nr+state%ml + 1) + hbaromega_z*(state%nz + 1.0_r_kind/2.0_r_kind)
end function

pure real(r_kind) elemental function phin(x,n, omega, mass) !!Eigenfunctions of 1d harmonic oscillator
    real(r_kind), intent(in) :: x !! [fm]
    real(r_kind), intent(in) :: omega !!In units MeV/hbar
    real(r_kind), intent(in) :: mass  !!In units MeV/c^2
    integer, intent(in) :: n
    real(r_kind) :: alpha, xi !!coordinate transform
    alpha = sqrt(mass*omega/(hbarc*hbarc)) !!Alpha is sqrt(mass*omega/hbar). We can use better units to get sqrt(mass/c^2 * omega/hbar / hbar)
    xi = alpha * x

    phin = exp(-xi**2 / 2.0_r_kind) * Hn(xi,n) *hofact(n,alpha)

end function

pure real(r_kind) elemental function hofact(n,alpha)
    integer, intent(in) :: n
    real(r_kind), intent(in) :: alpha
    hofact = sqrt(alpha)*pi**(-1.0_r_kind/4.0_r_kind) / sqrt(real(2**n * fac(n),r_kind))

end function

pure real(r_kind) elemental function phi2drad(r,nr,ml, mass, omega) result(val) !!Eigenfunction radial part of 2d harmonic oscllator
    real(r_kind), intent(in) :: r !!In units fm
    real(r_kind), intent(in) :: omega !!In units MeV/hbar
    real(r_kind), intent(in) :: mass  !!In units MeV/c^2
    integer, intent(in) :: nr, ml 
    real(r_kind) :: alpha, ar2

    alpha = sqrt(mass*omega/(hbarc*hbarc))
    ar2 = (alpha*r)**2
    val = r ** ml * exp(-ar2/2) * lna(ar2,nr,ml)
    val = val * ho_fact2d(alpha,ml,nr)!!normalization constant
end function

pure real(r_kind) elemental function ho_fact2d(alpha,ml,nr)
    real(r_kind), intent(in) :: alpha
    integer, intent(in) :: ml, nr
    ho_fact2d =  alpha**(ml+1) * sqrt(real(2*fac(nr),r_kind) / real(pi * fac(nr+ml),r_kind))
end function

pure elemental function phi2radang(ml, ang) result(phase) !!Eigenfunction ang part of 2d harmonic oscillator.
    integer, intent(in) :: ml
    real(r_kind), intent(in) :: ang
    complex(kind=r_kind) :: phase
    phase = complex(cos(ang),sin(ang))

end function

pure real(r_kind) recursive elemental function lna(x,n,a) result(val)
    real(r_kind), intent(in) :: x
    integer, intent(in) :: n, a


    if(n == 0) then
        val = 1.0_r_kind
    elseif(n==1) then
        val = 1.0_r_kind + a - x
    else
        val = (2*n-1+a-x)*lna(x,n-1,a) - (n-1+a)*lna(x,n-2,a)
        val = val/real(n,r_kind)
    endif

end function

pure integer recursive function fac(n) result(res)
    integer, intent(in) :: n
    if(n == 0 .or. n == 1) then
        res = 1
    else
        res = n * fac(n-1)
    endif
end function




end module def_ho