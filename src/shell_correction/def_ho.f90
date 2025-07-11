module def_ho
    use constants, only: r_kind, hbarc, pi, i_kind, use_ml_sym
    use micmac, only: deformation
    use optimise, only:func_1d, conj_grad_method
    implicit none
    private
    public :: an_ho_state, get_ho_states, getnumstates, betadef, getnumstatesupto, fac, Hn, lna, alpha, get_lowest_ho_states, get_ho_states_upto
    public :: kronecker, kin_en, gnl, gnlp, hmn, hmnp, spherical_def



    type :: betadef
        real(r_kind) :: beta2, beta4

        contains
        procedure :: eq => def_equals
        procedure :: omega_perp => compute_omega_perp
        procedure :: omega_z => compute_omega_Z
        procedure :: omega_def => compute_omega_def

    end type

    type(betadef), parameter :: spherical_def = betadef(beta2=0, beta4=0)



    type :: an_ho_state !!mj is half integer. so mj = 3 in reality means 3/2
        integer :: nz, nr, ml, mj, r, s, pi, N, nperp, ms
        contains
        procedure :: text => aniho_text
        procedure :: header => aniho_header
        procedure :: energy => ho_en
    end type
    !! Computes single particle energies in a shell model potential.
    contains 

pure logical function def_equals(self, other) result(equality)
    class(betadef), intent(in) :: self,other
    equality = self%beta2==other%beta2 .and. self%beta4 == other%beta4 

end function
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

pure function ho_en(self, hbaromega_z, hbaromega_perp) result(en)
    class(an_ho_state), intent(in) :: self
    real(r_kind) :: en
    real(r_kind), intent(in) :: hbaromega_z, hbaromega_perp

    en = hbaromega_z * (self%nz + 0.5_r_kind) + hbaromega_perp * (self%nperp + 1.0_r_kind)
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

function get_ho_states_upto(N) result(states)
    integer, intent(in) :: N !!ho quantum number
    type(an_ho_state), dimension(getnumstatesupto(N)) :: states
    integer :: nn, idx, shelldegen
    idx = 1
    do nn = 0, N
        shelldegen = getnumstates(nn)
        states(idx:idx+shelldegen - 1) = get_ho_states(nn)
        idx = idx + shelldegen
    end do
end function

function get_lowest_ho_states(N_max, n_states, hbaromega_z, hbaromega_perp) result(states)
    integer, intent(in) :: N_max, n_states !!ho quantum number
    type(an_ho_state), dimension(n_states) :: states
    type(an_ho_state), dimension(getnumstatesupto(N_max)) :: states_all
    type(an_ho_state) :: state
    real(r_kind), intent(in) :: hbaromega_z, hbaromega_perp
    real(r_kind) :: E
    integer :: idx, jj
    states_all = get_ho_states_upto(N_max)
    do idx = 2, size(states_all)
        state = states_all(idx)
        E = state%energy(hbaromega_z, hbaromega_perp)
        jj = idx - 1

        do while(jj .ge. 1 .and. states_all(jj)%energy(hbaromega_z, hbaromega_perp) .ge. E)
            states_all(jj + 1) = states_all(jj)
            jj = jj - 1
        end do
        states_all(jj+1) = state

    end do
    states = states_all(1:n_states)

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

                p = (-1)**(2*nr+ml+nz)
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

real(r_kind) function phin(x,n, omega, mass) !!Eigenfunctions of 1d harmonic oscillator
    real(r_kind), intent(in) :: x !! [fm]
    real(r_kind), intent(in) :: omega !!In units MeV/hbar
    real(r_kind), intent(in) :: mass  !!In units MeV/c^2
    integer, intent(in) :: n
    real(r_kind) :: alph, xi !!coordinate transform
    alph =  alpha(mass, omega) !!Alpha is sqrt(mass*omega/hbar). We can use better units to get sqrt(mass/c^2 * omega/hbar / hbar)
    xi = alph * x

    phin = exp(-xi**2 / 2.0_r_kind) * Hn(xi,n) *hofact(n,alph)

end function

real(r_kind) function hofact(n,alpha)
    integer, intent(in) :: n
    real(r_kind), intent(in) :: alpha
    hofact = sqrt(alpha)*pi**(-1.0_r_kind/4.0_r_kind) / sqrt(real(2**n * fac(n),r_kind))

end function

real(r_kind) function phi2drad(r,nr,ml, mass, omega) result(val) !!Eigenfunction radial part of 2d harmonic oscllator
    real(r_kind), intent(in) :: r !!In units fm
    real(r_kind), intent(in) :: omega !!In units MeV/hbar
    real(r_kind), intent(in) :: mass  !!In units MeV/c^2
    integer, intent(in) :: nr, ml 
    real(r_kind) :: alph, ar2

    alph =  alpha(mass, omega)
    ar2 = (alph*r)**2
    val = r ** ml * exp(-ar2/2) * lna(ar2,nr,real(ml,r_kind))
    val = val * ho_fact2d(alph,ml,nr)!!normalization constant
end function

real(r_kind) function ho_fact2d(alpha,ml,nr)
    real(r_kind), intent(in) :: alpha
    integer, intent(in) :: ml, nr
    ho_fact2d =  alpha**(ml+1) * sqrt(real(2*fac(nr),r_kind) / real(pi * fac(nr+ml),r_kind))
end function

pure real(r_kind) elemental function lna(x,n,a) result(res) !gen laguerre polynomial iterative version
    real(r_kind), intent(in) :: x
    integer, intent(in) :: n
    real(r_kind), intent(in) :: a
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


real(r_kind) function fac(n) result(f) !!factorial function which saves result. Real as integer*8 overflows at 21!
    integer, intent(in) :: n
    real(r_kind), save :: fc(0:169)
    logical, save :: first = .true.
    integer :: ii
    if(n < 0) THEN
        write(*,*) "error: negative argument for factorial"
        stop
    endif
    IF (first) THEN
        fc(0) = 1.0
        fc(1) = 1.0
        DO ii = 2, 169
        fc(ii) = fc(ii-1)*ii
        END DO
        first = .FALSE.
    END IF

    f = fc(n)
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
        do col = 1, size(states)
            sc = states(col)
            nzc = sc%nz
            nrc = sc%nr
            mlc = sc%ml
            msc = sc%ms
            npc = sc%nperp
            !kin_en(row, col) = kronecker(nrr, nrc) * kronecker(nzr, nzc) * kronecker(mlr, mlc) * &
            !0.5_r_kind * (hbaromega_perp * (npc+ 1) + hbaromega_z * (nzc + 1))
            ! if(msr /= msc .or. mlr /= mlc) then! .or. mlr /= mlc) then
            !     cycle

            ! endif
           
            ! if(nzr == nzc .and. nrc == nrr) then 
            !     kin_en(row, col) = 0.5_r_kind * (hbaromega_perp * (npc + 1) + hbaromega_z*(nzc + 0.5))!diag
            ! else if(nzr == nzc .and.  nrr == nrc + 1) then
            !     kin_en(row, col) = 0.5_r_kind*hbaromega_perp*sqrt(real(nrr*(nrr+mlr)))
            ! else if(nzr == nzc + 2 .and. nrc == nrr) then
            !     kin_en(row, col) = - 0.25_r_kind*hbaromega_z*sqrt(real(nzr*(nzr-1)))
            ! endif
            kin_en(row, col) = kin_en(row,col) + kronecker(nzr, nzc) * kronecker(nrc, nrr) * 0.5_r_kind * (hbaromega_perp * (npc + 1) + hbaromega_z*(nzc + 0.5))!diag

            kin_en(row, col) = kin_en(row,col) + kronecker(nzr, nzc) * kronecker(nrc+1, nrr) * 0.5_r_kind*hbaromega_perp*sqrt(real(nrr*(nrr+mlr)))!r -1

            kin_en(row, col) = kin_en(row,col) - kronecker(nzr, nzc+2) * kronecker(nrc, nrr) * 0.25_r_kind*hbaromega_z*sqrt(real(nzr*(nzr-1)))!z - 2

            kin_en(row, col) = kin_en(row,col) * kronecker(msr, msc) * kronecker(mlr, mlc)
        end do
    end do
    

    do row = 1,size(states)
        do col = row, size(states)
            kin_en(row,col) = kin_en(col,row)
        end do
    end do

end function




real(r_kind) function gnl(x,n,l) !!modified laguerre polynomial
    integer, intent(in) :: n, l
    real(r_kind), intent(in) :: x
    if(n < 0)then
        gnl = 0
    else
        gnl = x**(l*1.0_r_kind/2.0) * sqrt(1.0_r_kind*fac(n)/(fac(n+l)*1.0_r_kind)) * Lna(x,n,real(l,r_kind))
    endif
end function

real(r_kind) function gnlp(x,n,l) !!modified laguerre polynomial derivative
    integer, intent(in) :: n, l
    real(r_kind), intent(in) :: x

    gnlp = (gnl(x,n,l) * (2*n+l-x) &
            - gnl(x,n-1,l) * 2.0* sqrt(real(n*(n+l),r_kind)))/sqrt(x)

end function


real(r_kind) function Hmn(x,n) !!modified hermite polynomial
    integer, intent(in) :: n
    real(r_kind), intent(in) :: x
    if(n < 0) then
        Hmn = 0
    else
        Hmn = Hn(x,n)/(sqrt(real(2_i_kind**n * sqrt(pi))) * sqrt(real(fac(n),r_kind)))
    endif
end function

real(r_kind) function Hmnp(x,n) !!mod hermite polynomial derivative
    integer, intent(in) :: n
    real(r_kind), intent(in) :: x

    Hmnp = Hmn(x,n-1)*sqrt(0.5_r_kind*n) - Hmn(x,n+1) *sqrt(0.5_r_kind*(n+1))

end function


end module def_ho