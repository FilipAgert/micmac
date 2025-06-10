module Hamiltonian
    use constants
    use def_ho, only: betadef, an_ho_state, fac, Hn, lna
    use quadrule, only:legendre_dr_compute, hermite_ek_compute, laguerre_ss_compute
    use optimise, only:func_1d, conj_grad_method
    implicit none


    private
    public :: diagonalize, Ws_mat_elem, nucleus

    integer, parameter :: gauss_order =99
    real(r_kind), dimension(gauss_order) :: her_x, her_w, lag_x, lag_w, leg_x, leg_w
    real(r_kind), dimension(gauss_order, gauss_order) :: V_Coul_precomputed!!Cached values of the coulomb potential.
    logical :: precomputed = .false.


    real(r_kind), parameter :: aws = 0.70_r_kind !! (fm) WS potential diffusiveness



    type, extends(func_1d) :: dist_min
        real(r_kind) :: r, theta, radius! radius normalized as r/R0
        type(betadef) :: def
    contains
        procedure :: eval => surfdisteval

    end type

    type nucleus
        integer :: ZZ
        type(betadef) :: def
        real(r_kind) :: radius


    contains
        procedure :: charge_density => coul_charge_density
        procedure :: Vc => getVc
    end type
    contains

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

        
        call legendre_dr_compute(gauss_order, leg_x, leg_w)


    end subroutine

        pure real(r_kind) function coul_charge_density(self) result(rho)
        class(nucleus), intent(in) :: self

        rho = 1.0_r_kind

    end function

    subroutine diagonalize(E, V, H) !!Diagnoalize H and return V and E as eigenvectors sorted with lowest energy first.
        real(r_kind), intent(in) :: H(:,:)
        real(r_kind), intent(out) :: E(size(H,1)), V(size(H,1),size(H,1))
        external :: dsyev !!lapack routine for diagnoalizing symmetric matrix

        !!lapack variables!!
        character(len=1), parameter :: jobz = 'V', uplo = 'U'
        integer :: N
        integer :: lda
        real(r_kind) :: w(size(H,1))
        real(r_kind), allocatable :: work(:)
        integer :: lwork
        integer :: info

        N = size(H,1)
        write(*,*) N
        lda = N
        lwork = 10*N
        allocate(work(lwork))
        call dsyev(jobz, uplo, N, H, lda, w, work, lwork, info)

        if (info /= 0) then
            write(*,'(A,I10)') "Error using dsyev, info:", info
            call exit
        endif
        E = W
        V = H
    end subroutine


    pure real(r_kind) elemental function spher_dist(r1,r2,theta1,theta2,phi2)
        real(r_kind), intent(in) ::r1,r2,theta1,theta2,phi2
        spher_dist = sqrt(r1**2+r2**2-2*r1*r2*(cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi2)))
    end function



    real(r_kind) function getVc(self, r, theta) result(pot)
        class(nucleus), intent(in) :: self
        real(r_kind), intent(in) :: r, theta
        !Do the integral int_V     dr^3/|r-r'| at the point r.
        integer :: iu, it, ix
        real(r_kind) :: u, t, x, int_theta, int_radius_ub, rolling
        call precompute_weights()

        rolling = 0.0_r_kind
        do iu = 1, gauss_order
            u = leg_x(iu)
            int_theta = acos(u)
            int_radius_ub = surfRadius(int_theta, self%def, self%radius)
            do it = 1, gauss_order
                t = leg_x(it)
                do ix = 1,gauss_order
                    x = leg_x(ix)

                    rolling = rolling + leg_w(iu)*leg_w(it)*leg_w(ix)*pi*int_radius_ub**3*(x+1.0_r_kind)**2 / (8.0_r_kind) &
                            / sqrt(r**2 + int_radius_ub**2 * (x+1.0_r_kind)**2 / 4.0_r_kind &
                                    + r*int_radius_ub*(x+1.0_r_kind)* (sin(theta)*sin(acos(u))*cos(pi+pi*t) + u*cos(theta)))
                end do

            end do

        end do

        pot = rolling
    end function




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
        real(r_kind) :: rolling, z, rho, alphaz, alpharho, u, x, factz, factrho, xpart, r, ang, dist
        !The phi integral turns out to be a kronecker delta * 2pi
        if (state1%ml /= state2%ml) then 
            elem = 0
            return
        else if(state1%mj /= state2%mj) then !Spin orthogonality
            elem = 0
            return
        endif
        elem = 1.0_r_kind !!from phi part.

        call precompute_weights()

        
        alphaz = sqrt(mass*omegaz/(hbarc*hbarc))
        alpharho = sqrt(mass*omegaperp/(hbarc*hbarc))
        factz = 1.0_r_kind/sqrt(pi * 2**(state1%nz+state2%nz) * fac(state1%nz)*fac(state2%nz))
        factrho = sqrt(real(fac(state1%nr) * fac(state2%nr),r_kind) / real(fac(state1%nr+state1%ml) * fac(state2%nr+state2%ml),r_kind))

        rolling = 0.0_r_kind
        do il = 1, gauss_order
            x = lag_x(il)
            rho = sqrt(x)/alpharho
            ! write(*,'(A,F10.3, A, F10.3)')"x:", x, " w:", lag_w(il)
            ! pause
            xpart = x ** ((state1%ml + state2%ml * 1.0_r_kind)/2.0_r_kind) * lna(x,state1%nr,state1%ml) * lna(x,state2%nr,state2%ml)
            do ih = 1, gauss_order
                u = her_x(ih)
                z = u/alphaz
                r = rad_cyl(rho, z)
                ang = theta_cyl(rho,z)
                ! write(*,'(A,F10.3,A,F10.3,A,F10.3)')"z:", z, ", r:", r, " ang(deg):", ang*180/pi
                rolling = rolling + lag_w(il) * her_w(ih) * Hn(u,state1%nz) * Hn(u, state2%nz) & !Z part
                        * xpart * Ws(r,ang,def, Ws_depth, radius)!!Xpart and mixed part.
                ! write(*,'(A,F10.3)')"Pot:", Ws(r,ang,def, Ws_depth, radius)
                ! write(*,*)
                ! if(rolling < -1e6) then
                !     write(*,*) rolling, lag_w(il), her_w(ih), u, z, r, ang, xpart
                !     pause 'press enter to continue'
                ! endif
            end do
            ! pause
        end do

        
        elem = elem * factz * factrho * rolling

    end function

    real(r_kind) pure elemental function theta_cyl(rho,z)  !!Coordinate transform from cylindrical coordinates to spherical theta
        real(r_kind), intent(in):: rho,z
        if (z < 0) then
            theta_cyl = ATAN(rho/z) + pi
        else
            theta_cyl = ATAN(rho/z)
        endif
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
        if(r < radius) then
            dist = - dist
        endif
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

end module