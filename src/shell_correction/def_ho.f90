module def_ho
    use constants, only: r_kind, hbarc, pi, i_kind, use_ml_sym
    use micmac, only: deformation
    use optimise, only:func_1d, conj_grad_method
    implicit none
    private
    public :: an_ho_state, get_ho_states, getnumstates, betadef, getnumstatesupto, fac, Hn, lna, pauli_p, pauli_m, pauli_z, R_mat, Rp_mat, S_mat, Sp_mat, alpha, cosphi_m, azp_mat, az_mat, isinphi_m
    public :: kronecker, lz, expiphi, kin_en, gnl, gnlp, hmn, hmnp



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
    num = (N+1)*(N+2)!!incl spin degeneracy
    if(use_ml_sym) num = num / 2
end function

pure elemental integer function getnumstatesupto(N) result(num) !!get number of states in all major shells up to N
    integer, intent(in) :: N
    num = (N+1)*(N+2)*(N+3)/3 !!incl spin degeneracy
    if(use_ml_sym) num = num / 2
end function

function get_ho_states(N) result(states)
    integer, intent(in) :: N !!ho quantum number
    type(an_ho_state), dimension(getnumstates(N)) :: states
    type(an_ho_state) :: state
    integer :: nz, nr, ml, omega, qrem, idx,r , s, p, nperp, ms, lb
    idx = 0
    do nz = 0,N
        !we need nz + 2nr + ml = N.
        qrem = N - nz
        do r = 0, qrem
            s = qrem - r
            ml = r - s
            if(use_ml_sym .and. ml < 0) continue
                
            
            nperp = r + s

            if(use_ml_sym) then
                lb = abs(2*ml-1)
            else
                lb = 2*ml-1
            endif
            do omega = lb, 2*ml + 1, 2 !abs(2*mlp - 1), 2*mlp + 1, 2  !2*mlp - 1, 2*mlp + 1, 2 
                ! write(*,*) "  N     nz    ml    mj"
                ! write(*,'(5I6)') N, nz, mlp, omega
                ms = omega-2*ml
                idx = idx + 1
                nr = (nperp-abs(ml))/2
                if(r-s /= ml) then
                    print*, "ERROR! R-S is not ml"
                elseif(r+s /= nperp) then
                    print*, "ERROR! R+S is not nperp"
                endif

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

pure real(r_kind) elemental function alpha(mass, omega)
    real(r_kind), intent(in) :: omega !!In units MeV/hbar
    real(r_kind), intent(in) :: mass  !!In units MeV/c^2
    alpha = sqrt(mass*omega/(hbarc*hbarc))
end function

pure real(r_kind) elemental function phin(x,n, omega, mass) !!Eigenfunctions of 1d harmonic oscillator
    real(r_kind), intent(in) :: x !! [fm]
    real(r_kind), intent(in) :: omega !!In units MeV/hbar
    real(r_kind), intent(in) :: mass  !!In units MeV/c^2
    integer, intent(in) :: n
    real(r_kind) :: alph, xi !!coordinate transform
    alph =  alpha(mass, omega) !!Alpha is sqrt(mass*omega/hbar). We can use better units to get sqrt(mass/c^2 * omega/hbar / hbar)
    xi = alph * x

    phin = exp(-xi**2 / 2.0_r_kind) * Hn(xi,n) *hofact(n,alph)

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
    real(r_kind) :: alph, ar2

    alph =  alpha(mass, omega)
    ar2 = (alph*r)**2
    val = r ** ml * exp(-ar2/2) * lna(ar2,nr,ml)
    val = val * ho_fact2d(alph,ml,nr)!!normalization constant
end function

pure real(r_kind) elemental function ho_fact2d(alpha,ml,nr)
    real(r_kind), intent(in) :: alpha
    integer, intent(in) :: ml, nr
    ho_fact2d =  alpha**(ml+1) * sqrt(real(2*fac(nr),r_kind) / real(pi * fac(nr+ml),r_kind))
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


pure elemental recursive function fac(n) result(res)
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
    integer :: row, col
    type(an_ho_state) :: sr, sc

    do row = 1, size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            R_mat(row, col) = (kronecker(sr%ms,sc%ms)*kronecker(sr%nz,sc%nz) * kronecker(sr%s, sc%s) * kronecker(sr%r, sc%r - 1)) * sqrt(real(sc%r,r_kind))
        end do
    end do

end function

pure function Sp_mat(states)
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: Sp_mat
    integer :: row, col, kron
    type(an_ho_state) :: sr, sc

    do row = 1, size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            kron = kronecker(sr%ms,sc%ms)*kronecker(sr%nz,sc%nz) * kronecker(sr%r, sc%r) * kronecker(sr%s, sc%s + 1)
            Sp_mat(row, col) = kron * sqrt(real(sc%s + 1,r_kind))
        end do
    end do

end function
pure function S_mat(states) !!Matrix for the R matrix
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: S_mat
    integer :: row, col, kron
    type(an_ho_state) :: sr, sc

    do row = 1, size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            kron = kronecker(sr%ms,sc%ms)*kronecker(sr%nz,sc%nz) * kronecker(sr%r, sc%r) * kronecker(sr%s, sc%s - 1)
            S_mat(row, col) =  kron * sqrt(real(sc%s,r_kind))
        end do
    end do
end function

pure function azp_mat(states) !!Constructor operator for z quanta
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: azp_mat
    integer :: row, col, kron
    type(an_ho_state) :: sr, sc

    do row = 1, size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            kron = kronecker(sr%ms,sc%ms)*kronecker(sr%nz,sc%nz+1) * kronecker(sr%r, sc%r) * kronecker(sr%s, sc%s)
            azp_mat(row, col) = kron * sqrt(real(sc%nz + 1,r_kind))
        end do
    end do
end function

pure function az_mat(states)!!Destructor operator for z quanta
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: az_mat
    integer :: row, col, kron
    type(an_ho_state) :: sr, sc

    do row = 1, size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            kron = kronecker(sr%ms,sc%ms)*kronecker(sr%nz,sc%nz-1) * kronecker(sr%r, sc%r) * kronecker(sr%s, sc%s)
            az_mat(row, col) = kron*  sqrt(real(sc%nz,r_kind))
        end do
    end do
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
    integer :: row, col
    type(an_ho_state) :: sr, sc
    pauli_m = 0
    do row = 1,size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            pauli_m(row, col) = kronecker(sr%ms,sc%ms-2)*kronecker(sr%nz,sc%nz) * kronecker(sr%r, sc%r) * kronecker(sr%s, sc%s)
        end do
    end do
end function

pure function cosphi_m(states) !! matrix w. elements <m_l' | cosphi | m_l > = 1/2 if |m_l-m_l'| = 1,
                               !! 0, otherwise
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: cosphi_m
    integer :: row, col
    type(an_ho_state) :: sr, sc
    do row = 1,size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            cosphi_m(row, col) = 0.5_r_kind*(kronecker(sr%ml,sc%ml+1) + kronecker(sr%ml,sc%ml-1))*   kronecker(sr%nz,sc%nz) * kronecker(sr%nperp, sc%nperp)* kronecker(sr%ms, sc%ms)
        end do
    end do
end function

pure function expiphi(states) !! matrix w. elements <m_l' | cosphi | m_l > = 1/2 if |m_l-m_l'| = 1,
                               !! 0, otherwise
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: expiphi
    integer :: row, col
    type(an_ho_state) :: sr, sc
    do row = 1,size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            expiphi(row, col) = (kronecker(sr%ml,sc%ml+1))*   kronecker(sr%nz,sc%nz) * kronecker(sr%nperp, sc%nperp)* kronecker(sr%ms, sc%ms)
        end do
    end do
end function

pure function isinphi_m(states) !! matrix w. elements <m_l' | isinphi | m_l > = \pm 1/2 if m_l-m_l' = \mp 1,
                               !! 0, otherwise
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: isinphi_m
    integer :: row, col
    type(an_ho_state) :: sr, sc
    do row = 1,size(states)
        sr = states(row)
        do col = 1,size(states)
            sc = states(col)
            isinphi_m(row, col) = 0.5_r_kind*(kronecker(sr%ml,sc%ml+1) - kronecker(sr%ml,sc%ml-1))*   kronecker(sr%nz,sc%nz) * kronecker(sr%nperp, sc%nperp) * kronecker(sr%ms, sc%ms)
        end do
    end do
end function


pure function lz(states)
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: lz
    integer :: row
    type(an_ho_state) :: sr
    lz = 0
    do row = 1,size(states)
        sr = states(row)
        lz(row,row) = sr%ml
    end do
end function

pure function kin_en(states, hbaromega_z, hbaromega_perp) !Nuclear Physics A A method for solving the independent-particle Schrödinger equation with a deformed average field
                                                            !Volume 135, Issue 2, 2 October 1969, Pages 432-444
    type(an_ho_state), intent(in) :: states(:)
    real(r_kind), dimension(size(states), size(states)) :: kin_en
    integer :: row, col, nzr, nzc, nrr, nrc, mlr, mlc, msr, msc, npc
    real(r_kind), intent(in) :: hbaromega_z, hbaromega_perp
    type(an_ho_state) :: sr,sc
    kin_en = 0
    do row = 1,size(states)
        sr = states(row)
        nzr = sr%nz
        nrr = sr%nr
        mlr = sr%ml
        msr = sr%ms
        do col = row, size(states)
            sc = states(col)
            nzc = sc%nz
            nrc = sc%nr
            mlc = sc%ml
            msc = sc%ms
            npc = sc%nperp
            !kin_en(row, col) = kronecker(nrr, nrc) * kronecker(nzr, nzc) * kronecker(mlr, mlc) * &
            !0.5_r_kind * (hbaromega_perp * (npc+ 1) + hbaromega_z * (nzc + 1))
            if(msr /= msc .or. mlr /= mlc) continue

            kin_en(row, col) = kin_en(row,col) + kronecker(nzr, nzc) * kronecker(nrc, nrr) * 0.5_r_kind * (hbaromega_perp * (npc + 1) + hbaromega_z*(nzc + 1))!diag

            kin_en(row, col) = kin_en(row,col) + kronecker(nzr, nzc) * kronecker(nrc-1, nrr) * 0.5_r_kind*hbaromega_perp*sqrt(real(nrr*(nrr+mlr)))!r -1

            kin_en(row, col) = kin_en(row,col) - kronecker(nzr, nzc-2) * kronecker(nrc, nrr) * 0.25_r_kind*hbaromega_z*sqrt(real(nzr*(nzr-1)))!z - 2

            kin_en(row, col) = kin_en(row,col) * kronecker(msr, msc) * kronecker(mlr, mlc)
        end do
    end do

    do row = 1,size(states)
        do col = 1, row
            kin_en(row,col) = kin_en(col,row)
        end do
    end do

end function




pure real(r_kind) elemental function gnl(x,n,l) !!modified laguerre polynomial
    integer, intent(in) :: n, l
    real(r_kind), intent(in) :: x
    gnl = x**(l*1.0_r_kind/2.0) * sqrt(1.0_r_kind*fac(n)/(fac(n+l))) * Lna(x,n,l)
end function

pure real(r_kind) elemental function gnlp(x,n,l) !!modified laguerre polynomial derivative
    integer, intent(in) :: n, l
    real(r_kind), intent(in) :: x
    if(n == 0) then
        gnlp = 0
    else
        gnlp = (gnl(x,n,l) * (2*n+l-x) - gnl(x,n-1,l) * 2* sqrt(real(n*(n+l),r_kind)))/sqrt(x)
    endif
end function


pure real(r_kind) elemental function Hmn(x,n) !!modified hermite polynomial
    integer, intent(in) :: n
    real(r_kind), intent(in) :: x
    Hmn = Hn(x,n)/sqrt(real(2**n * fac(n) * sqrt(pi),r_kind))
end function

pure real(r_kind) elemental function Hmnp(x,n) !!mod hermite polynomial derivative
    integer, intent(in) :: n
    real(r_kind), intent(in) :: x
    if(n == 0) then
        Hmnp = 0
    else
        Hmnp = Hmn(x,n-1)*sqrt(0.5_r_kind*n) - Hmn(x,n) *sqrt(0.5_r_kind*(n+1))
    endif
end function


end module def_ho