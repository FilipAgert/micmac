module def_ho
    use constants, only: r_kind, hbarc, pi, i_kind
    use micmac, only: deformation
    use optimise, only:func_1d, conj_grad_method
    implicit none
    private
    public :: an_ho_state, get_ho_states, getnumstates, betadef, getnumstatesupto, fac, Hn, lna



    type :: betadef
        real(r_kind) :: beta2, beta4

        contains
        procedure :: omega_perp => compute_omega_perp
        procedure :: omega_z => compute_omega_Z
        procedure :: omega_def => compute_omega_def

    end type



    type :: an_ho_state !!mj is half integer. so mj = 3 in reality means 3/2
        integer :: nz, nr, ml, mj, r, s, pi, N, nperp, ms
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
    character(len=1) :: par
    if(self%pi == -1) then
        par = "-"
    else
        par = "+"
    endif
    write(text, '(A,5I5,A,I5,A, 3I5)')par, self%N, self%nz, self%nr, self%ml, self%mj, "/2", self%ms, "/2", self%r, self%s, self%nperp
end function

function aniho_header(self) result(hdr)
    class(an_ho_state) :: self
    character(len=100) :: hdr
    write(hdr, '(A)') "p    N    nz   nr   ml    mj   ms    r    s   n_p"

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
    integer :: nz, nr, ml, omega, qrem, idx,r , s, p, nperp, sign, msigned, ms
    idx = 0
    do nz = 0,N
        !we need nz + 2nr + ml = N.
        qrem = N - nz
        do nr = 0, qrem/2 !!integer division. if remaining is 3, we want n = 0, 1
            ml = qrem - 2*nr
           
            do omega = abs(2*ml - 1), 2*ml + 1, 2 
                ms = omega-2*ml
                idx = idx + 1
                nperp = 2*nr + abs(ml)
                r = (nperp + ml)/2
                s = (nperp - ml)/2


                p = (-1)**N
                state = an_ho_state(nz = nz, nr = nr, ml=ml, mj=omega, r= r, s=s, N=N, pi=p, nperp = nperp, ms=ms)
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


pure real(r_kind) elemental function Hn(x,n) result(res)
    real(r_kind), intent(in) :: x
    integer, intent(in) :: n
    real(r_kind), dimension(0:n) :: Hns
    integer :: nn
    Hns(0) = 1.0_r_kind
    if(n > 0) then
        Hns(1) = 2.0_r_kind * x
        do nn = 2,n
            Hns(nn) = 2.0_r_kind*x*Hns(nn-1) - 2.0_r_kind*(nn-1)*Hns(nn-2)
        end do
    endif
    res = Hns(n)

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




pure real(r_kind) elemental function lna(x,n,a) result(res) !gen laguerre polynomial iterative version
    real(r_kind), intent(in) :: x
    integer, intent(in) :: n, a
    real(r_kind), dimension(0:n) :: lnas
    integer :: nn
    lnas(0) = 1.0_r_kind

    if(n > 0) then
        lnas(1) = 1.0_r_kind + real(a,r_kind) - x

        do nn = 2, n
            lnas(nn) = (real(2*nn-1+a,r_kind)-x)*lnas(nn-1) - (real(nn-1+a,r_kind))*lnas(nn-2)
            lnas(nn) = lnas(nn)/real(nn,r_kind)
        end do
    endif
    res = lnas(n)

end function


pure recursive function fac(n) result(res)
    integer(kind=i_kind) :: res
    integer, intent(in) :: n
    if(n == 0 .or. n == 1) then
        res = 1
    else
        res = n* fac(n-1)
    endif
end function


pure function Rp_mat(states)!!Matrix for the R^+ matrix
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: Rp_mat
    integer :: row, col
    type(an_ho_state) :: sr, sc

    do row = 1, size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            Rp_mat(row, col) = kronecker(sr%ms,sc%ms)*kronecker(sr%nz,sc%nz) * kronecker(sr%s, sc%s) * kronecker(sr%r, sc%r + 1) * sqrt(real(sc%r + 1,r_kind))
        end do
    end do

end function

pure function R_mat(states) !!Matrix for the R matrix
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: R_mat
    R_mat = transpose(Rp_mat(states))
end function

pure function Sp_mat(states)
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: Sp_mat
    integer :: row, col
    type(an_ho_state) :: sr, sc

    do row = 1, size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            Sp_mat(row, col) = kronecker(sr%ms,sc%ms)*kronecker(sr%nz,sc%nz) * kronecker(sr%r, sc%r) * kronecker(sr%s, sc%s + 1) * sqrt(real(sc%s + 1,r_kind))
        end do
    end do

end function
pure function S_mat(states) !!Matrix for the R matrix
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: S_mat
    S_mat = transpose(Sp_mat(states))
end function

pure function azp_mat(states) !!Constructor operator for z quanta
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: azp_mat
    integer :: row, col
    type(an_ho_state) :: sr, sc

    do row = 1, size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            azp_mat(row, col) = kronecker(sr%ms,sc%ms)*kronecker(sr%nz,sc%nz+1) * kronecker(sr%r, sc%r) * kronecker(sr%s, sc%s) * sqrt(real(sc%nz + 1,r_kind))
        end do
    end do
end function

pure function az_mat(states)!!Destructor operator for z quanta
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: az_mat
    az_mat = transpose(azp_mat(states))
end function

pure integer function kronecker(a,b)
    integer, intent(in) :: a,b
    if(a == b) then
        kronecker = 1
    else
        kronecker = 0
    endif

end function

!!momentum operators




!! pauli matrices
pure function pauli_z(states)
    type(an_ho_state), intent(in) :: states(:)
    integer, dimension(size(states), size(states)) :: pauli_z
    integer :: n
    pauli_z = 0
    do n = 1,size(states)
        pauli_z(n,n) = states(n)%ms
    end do
end function

pure function pauli_x(states)
    type(an_ho_state), intent(in) :: states(:)
    integer, dimension(size(states), size(states)) :: pauli_x
    pauli_x = pauli_p(states) + pauli_m(states)
end function

pure function pauli_yim(states) !!This is the imaginary part of the pauli_y matrix.
    type(an_ho_state), intent(in) :: states(:)
    integer, dimension(size(states), size(states)) :: pauli_yim
    pauli_yim =  pauli_m(states) - pauli_p(states)
end function

pure function pauli_p(states)
    type(an_ho_state), intent(in) :: states(:)
    integer, dimension(size(states), size(states)) :: pauli_p
    integer :: row, col
    type(an_ho_state) :: sr, sc
    pauli_p = 0
    do row = 1,size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            pauli_p(row, col) = kronecker(sr%ms,sc%ms+2)*kronecker(sr%nz,sc%nz) * kronecker(sr%r, sc%r) * kronecker(sr%s, sc%s)
        end do
    end do
end function

pure function pauli_m(states) !!Destructor operator for pauli matrix
    type(an_ho_state), intent(in) :: states(:)
    integer, dimension(size(states), size(states)) :: pauli_m
    pauli_m = transpose(pauli_p(states))
end function

end module def_ho