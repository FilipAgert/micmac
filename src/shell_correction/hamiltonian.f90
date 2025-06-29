module Hamiltonian
    use constants
    use def_ho
    use quadrule, only:legendre_dr_compute, hermite_ek_compute, laguerre_ss_compute
    use optimise, only:func_1d, conj_grad_method
    implicit none


    private
    public :: diagonalize, mat_elem_axsym, WS_pot, VC_pot, dist_min, Vso_mat, commutator, inv_rho, unit, rhopot, ddrhomat, rhomat, inv_rho_mat

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



    type, abstract :: potential
        type(betadef) :: def
        real(r_kind) :: radius
        real(r_kind), dimension(gauss_order, gauss_order) :: precomputed_pot
        logical, dimension(gauss_order, gauss_order) :: is_computed = .false.
    contains
        procedure(eval), deferred :: eval
        procedure :: set => set_precomp
        procedure :: get => get_precomp
        procedure :: eval_pre => eval_precomp
    end type
    !!!!!!!!
    type, extends(potential) :: WS_pot
        real(r_kind) :: V_0
    contains
        procedure :: eval => eval_WS
    end type

    type,extends(potential) :: inv_rho

    contains
        procedure :: eval => eval_inv_rho
    end type

    type,extends(potential) :: rhopot

    contains
        procedure :: eval => eval_rho
    end type

    type,extends(potential) :: unit

    contains
        procedure :: eval => eval_unit
    end type
    !!!!!!!!!!!!!!!!!!!!!
    type, abstract, extends(WS_pot) :: WS_derivative
    contains
        procedure :: pot => eval_WS_wrapper
    end type
    !!!!!!!!!!!!!!!!!!!!!
    type, extends(WS_derivative) :: dWs_dz
    contains
        procedure :: eval => eval_dWs_dz
    end type

    !!!!!!!!!!!!!!
    type, extends(WS_derivative) :: dWs_drho

    contains
        procedure :: eval => eval_dWs_drho
    end type


    type, extends(potential) :: VC_pot
        real(r_kind) :: charge_dens
    contains
        procedure :: eval => eval_VC
        procedure :: set_charge_dens => set_charge_VC
    end type

    abstract interface
        real(r_kind) function eval(self, r, theta)
            import:: potential, r_kind
            class(potential), intent(in) :: self
            real(r_kind), intent(in) :: r, theta
        end function
    end interface


    contains
    real(r_kind) function eval_precomp(self, ir, it, r, theta)
        class(potential) :: self
        integer, intent(in) :: ir, it
        real(r_kind), intent(in) :: r, theta

        if(self%is_computed(ir, it)) then
            eval_precomp = self%precomputed_pot(ir,it)
        else
            eval_precomp = self%eval(r, theta)
            self%precomputed_pot(ir, it) = eval_precomp
            self%is_computed(ir,it) = .true.
        endif

    end function
    pure real(r_kind) function get_precomp(self, ir, it)
        class(potential), intent(in) :: self
        integer, intent(in) :: ir, it
        get_precomp = self%precomputed_pot(ir, it)
    end function

    subroutine set_precomp(self, ir, it, val)
        class(potential) :: self
        integer, intent(in) :: ir, it
        real(r_kind), intent(in) :: val
        self%precomputed_pot(ir, it) = val
    end subroutine



    subroutine precompute()
        if(precomputed) then
            return
        endif


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

        precomputed = .true.




    end subroutine

    subroutine set_charge_VC(self, ZZ)
        class(VC_pot) :: self
        integer, intent(in) :: ZZ !!proton number
        self%charge_dens = coul_charge_density(ZZ, self%radius)
    end subroutine

    pure real(r_kind) function coul_charge_density(ZZ, radius) result(rho)
        integer, intent(in) :: ZZ !!proton number
        real(r_kind), intent(in) :: radius
        rho = ZZ *3.0_r_kind/(4.0_r_kind * pi*4*pi * epsilonzero * radius**3)  !! Charge /( Volume * 4piepsilonzero)

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


    real(r_kind) function eval_vc(self,r, theta) !!axially symmetric. function of z and rho
        class(VC_pot), intent(in) :: self
        real(r_kind), intent(in) :: r, theta
        !Do the integral int_V     dr^3/|r-r'| at the point r.
        integer :: iu, it, ix
        real(r_kind) :: u, t, x, int_theta, int_radius_ub, rolling
        call precompute()

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
        eval_vc = rolling * self%charge_dens
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


    real(r_kind) function mat_elem_axsym(state1, state2, pot, mass, omegaz, omegaperp) result(elem) 
        implicit none
        !!Calculates matrix element for central part of WS potential.
        !!uses 64 points gaussian quadrature
        type(an_ho_state), intent(in) :: state1, state2
        class(potential)  :: pot
        real(r_kind), intent(in) :: mass !!In units MeV/c^2
        real(r_kind), intent(in) :: omegaperp, omegaz !!In units MeV/hbar
        integer :: il, ih
        real(r_kind) :: rolling, z, rho, alphaz, alpharho, u, x, factz, factrho, xpart, r, ang
        !The phi integral turns out to be a kronecker delta * 2pi
        if (state1%ml /= state2%ml) then 
            elem = 0.0_r_kind
            return
        endif
        
        if(state1%ms /= state2%ms) then !Spin orthogonality
            elem = 0.0_r_kind
            return
        endif
        elem = 1.0_r_kind !!from phi part.

        call precompute()

        
        alphaz = alpha(mass, omegaz)
        alpharho = alpha(mass,omegaperp)

        !in factors, split up the square roots to avoid integer overflow
        factz = 1.0_r_kind/(sqrt(pi * 2_i_kind**(state1%nz+state2%nz))*sqrt(real(fac(state1%nz),r_kind))*sqrt(real(fac(state2%nz),r_kind))) 
        factrho = sqrt(real(fac(state1%nr),r_kind))*sqrt(real( fac(state2%nr),r_kind)) / (sqrt(real(fac(state1%nr+abs(state1%ml)),r_kind)) * sqrt(real(fac(state2%nr+abs(state2%ml)),r_kind)))

        rolling = 0.0_r_kind
        do il = 1, gauss_order
            x = lag_x(il)
            rho = sqrt(x)/alpharho
            ! write(*,'(A,F10.3, A, F10.3)')"x:", x, " w:", lag_w(il)
            ! pause
            xpart = x ** ((abs(state1%ml)*1.0_r_kind + abs(state2%ml) * 1.0_r_kind)/2.0_r_kind) * lna(x,state1%nr,abs(state1%ml)) * lna(x,state2%nr,abs(state2%ml))
            do ih = 1, gauss_order
                u = her_x(ih)
                z = u/alphaz
                r = rad_cyl(rho, z)
                ang = theta_cyl(rho,z)
                ! write(*,'(A,F10.3,A,F10.3,A,F10.3)')"z:", z, ", r:", r, " ang(deg):", ang*180/pi
                rolling = rolling + lag_w(il) * her_w(ih) * Hn(u,state1%nz) * Hn(u, state2%nz) & !Z part
                        * xpart * pot%eval_pre(il, ih, r,ang)   !!Xpart and mixed part.
                if(isnan(xpart) .or. isnan(rolling)) then
                    write(*,*) "u, x", u, x
                    write(*,*) lag_w(il),her_w(ih), Hn(u,state1%nz), Hn(u, state2%nz) & !Z part
                    ,xpart , pot%eval_pre(il, ih, r,ang)

                    pause
                endif
                if(rho < 1e-6_r_kind) then
                    write(*,*) rho
                    write(*,*) lag_w(il) * her_w(ih) * Hn(u,state1%nz) * Hn(u, state2%nz) & !Z part
                    * xpart

                endif

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
        if(isnan(elem)) then
            write(*,*) "elem, z fact, rho fact, sum", elem, factz, factrho, rolling
            write(*,*) "nr1, ml1, nr2, ml2", state1%nr, state1%ml, state2%nr, state2%ml
            write(*,*) fac(state1%nr),fac(state2%nr),fac(state1%nr+state1%ml), fac(state2%nr+state2%ml)
            write(*,*) real(fac(state1%nr)*fac(state2%nr),r_kind),(sqrt(real(fac(state1%nr+state1%ml),r_kind)) * sqrt(real(fac(state2%nr+state2%ml),r_kind)))
            pause
        endif

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

    real(r_kind) function eval_ws(self,r, theta) !!axially symmetric. function of z and rho
        real(r_kind), intent(in) :: r, theta
        class(WS_pot), intent(in) :: self
        real(r_kind) ::dist
        dist = sqrt(surfdist(r,theta,self%def, self%radius)) 
        if(r < self%radius) then
            dist = - dist
        endif
        eval_ws = -self%V_0/ (1 + exp(dist/aws))
    end function

    real(r_kind) function eval_inv_rho(self,r, theta) !!axially symmetric. function of z and rho
        real(r_kind), intent(in) :: r, theta
        class(inv_rho), intent(in) :: self
        real(r_kind) :: rho
        rho = r*sin(theta)
        if (rho == 0.0_r_kind) rho = 1e-12_r_kind
        eval_inv_rho = 1.0_r_kind/rho
    end function

    real(r_kind) function eval_rho(self,r, theta) !!axially symmetric. function of z and rho
        real(r_kind), intent(in) :: r, theta
        class(rhopot), intent(in) :: self
        real(r_kind) :: rho
        rho = r*sin(theta)
        eval_rho = rho
    end function


    real(r_kind) function eval_ws_wrapper(self,r, theta) !!Calls the 
        real(r_kind), intent(in) :: r, theta
        class(WS_derivative), intent(in) :: self
        eval_ws_wrapper = eval_ws(self, r, theta)
    end function


    real(r_kind) function eval_dWs_dz(self, r, theta)
        real(r_kind), intent(in) :: r, theta
        class(dWs_dz), intent(in) :: self
        real(r_kind) :: rprime, thetaprime
        real(r_kind), parameter :: dz = 1e-6
        
        rprime = r + cos(theta)*dz
        thetaprime = theta - sin(theta) * dz /r

        eval_dWs_dz = (self%pot(rprime, thetaprime) - self%pot(r, theta)) / dz
    end function

    real(r_kind) function eval_dWs_drho(self, r, theta)
        real(r_kind), intent(in) :: r, theta
        class(dWs_drho), intent(in) :: self
        real(r_kind) :: rprime, thetaprime
        real(r_kind), parameter :: drho = 1e-6
        
        rprime = r + sin(theta)*drho !! new r after we have rho' = rho+drho
        thetaprime = theta + cos(theta) * drho / r !!new theta after rho' = rho+drho

        eval_dWs_drho = (self%pot(rprime, thetaprime) - self%pot(r, theta)) / drho
    end function

    real(r_kind) function eval_unit(self, r, theta)
        real(r_kind), intent(in) :: r, theta
        class(unit), intent(in) :: self
        eval_unit = 1.0_r_kind
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


    function inv_rho_mat(states,omegaz, omegaperp, mass) !!matrix for 1/rho
        type(an_ho_state), intent(in) :: states(:)
        real(r_kind), intent(in) :: omegaz, omegaperp !!Frequencies in units MeV/hbar
        real(r_kind), intent(in) :: mass!!Frequencies in units MeV/c^2
        real(r_kind), dimension(size(states), size(states)) :: inv_rho_mat
        integer :: row, col
        type(an_ho_state) :: sr, sc
        type(inv_rho) :: invrhopot
        inv_rho_mat = 0
        do row = 1,size(states)
            sr = states(row)
            do col = 1,size(states)
                sc = states(col)
                inv_rho_mat(row,col) = mat_elem_axsym(sr, sc, invrhopot, mass ,omegaz, omegaperp)
            end do
        end do
    end function

    function rhomat(states,omegaz, omegaperp, mass) !!rho
        type(an_ho_state), intent(in) :: states(:)
        real(r_kind), intent(in) :: omegaz, omegaperp !!Frequencies in units MeV/hbar
        real(r_kind), intent(in) :: mass!!Frequencies in units MeV/c^2
        real(r_kind), dimension(size(states), size(states)) :: rhomat
        integer :: row, col
        type(an_ho_state) :: sr, sc
        type(rhopot) :: rhop
        rhomat = 0
        do row = 1,size(states)
            sr = states(row)
            do col = 1,size(states)
                sc = states(col)
                rhomat(row,col) = mat_elem_axsym(sr, sc, rhop, mass ,omegaz, omegaperp)
            end do
        end do
    end function

    function ddrhomat(states, omegaperp, mass,rmat, invrhomat) !!matrix for 1/rho
        type(an_ho_state), intent(in) :: states(:)
        real(r_kind), intent(in) :: omegaperp !!Frequencies in units MeV/hbar
        real(r_kind), intent(in) :: mass!!Frequencies in units MeV/c^2
        real(r_kind), dimension(size(states), size(states)) :: ddrhomat, rmat, invrhomat
        integer :: row, col, np, ml
        type(an_ho_state) :: sr, sc
        
        do row = 1,size(states)
            sr = states(row)
            do col = 1,size(states)
                sc = states(col)
                np = sc%nr
                ml = sc%ml
                ddrhomat(row,col) = sc%ml * invrhomat(row, col)  &
            - 2.0_r_kind*  kronecker(sr%nr, sc%nr - 1) * kronecker(sr%ml, sc%ml + 1) * sqrt(np*(np+ml+1.0_r_kind)*(np+ml+2.0_r_kind))
            end do
        end do
        ddrhomat = ddrhomat -alpha(mass, omegaperp)**2 * rmat
    end function

    function VWS_mat(states, def, R0, omegaz, omegaperp, mass, WS_depth, lambda)
        implicit none
        type(an_ho_state), intent(in) :: states(:)
        type(betadef), intent(in) :: def !!deformation of body
        real(r_kind), intent(in) :: R0 !!radius as r0*A**(1/3). Where r0 is r0_so from constants.f90
        real(r_kind), intent(in) :: WS_depth !!WS potental well depth
        real(r_kind), dimension(size(states),size(states)) :: VWS_mat
        real(r_kind), intent(in) :: omegaz, omegaperp !!Frequencies in units MeV/hbar
        real(r_kind), intent(in) :: mass!!Frequencies in units MeV/c^2
        real(r_kind), intent(in) ::lambda !!multiplicative factor
        real(r_kind) :: alpha_z, alpha_perp
        type(WS_pot) :: VWS
        integer ::numstates, row, col
        type(an_ho_state) :: s1, s2
        call precompute()

        !setup potential
        VWS%def = def
        VWS%radius = R0
        VWS%V_0 = WS_depth

        alpha_z = alpha(mass, omegaz)
        alpha_perp = alpha(mass, omegaperp)

        numstates = size(states)
        VWS_mat = 0
        do row = 1, numstates
            s1 = states(row)
            do col = 1,row
                s2 = states(col)
                VWS_mat(row,col) = VWS_elem(s1,s2,VWS, alpha_z, alpha_perp)
            end do
        end do

        do row = 1,numstates
            do col = row +1, numstates
                VWS_mat(row, col) = VWS_mat(col, row)
            end do
        end do
    end function

    function Vso_mat(states, def, R0, omegaz, omegaperp, mass, WS_depth, lambda)
        implicit none
        type(an_ho_state), intent(in) :: states(:)
        type(betadef), intent(in) :: def !!deformation of body
        real(r_kind), intent(in) :: R0 !!radius as r0*A**(1/3). Where r0 is r0_so from constants.f90
        real(r_kind), intent(in) :: WS_depth !!WS potental well depth
        real(r_kind), dimension(size(states),size(states)) :: Vso_mat
        real(r_kind), intent(in) :: omegaz, omegaperp !!Frequencies in units MeV/hbar
        real(r_kind), intent(in) :: mass!!Frequencies in units MeV/c^2
        real(r_kind), intent(in) ::lambda !!multiplicative factor
        real(r_kind) :: alpha_z, alpha_perp, kappa
        type(WS_pot) :: so_ws
        integer ::numstates, row, col, l1, l2, ms1, ms2
        type(an_ho_state) :: s1, s2
        call precompute()

        kappa = -lambda*(hbarc/(2*mass))**2
        !setup potential
        so_ws%def = def
        so_ws%radius = R0
        so_ws%V_0 = WS_depth

        alpha_z = alpha(mass, omegaz)
        alpha_perp = alpha(mass, omegaperp)

        numstates = size(states)
        Vso_mat = 0
        do row = 1, numstates
            s1 = states(row)
            l1 = s1%ml
            ms1 = s1%ms
            do col = 1,row
                s2 = states(col)
                l2 = s2%ml
                ms2 = s2%ms
                print*, l1, ms1, l2, ms2
                
                if(ms1 == ms2 .and. l1 == l2) then
                    Vso_mat(row,col) = VSO_diag_elem(s1,s2,so_ws,alpha_z,alpha_perp)
                    print*, "diag"
                elseif((ms1 == ms2 + 1 .and. l1 == l2 - 1) .or. (ms1 == ms2 - 1 .and. l1 == l2 + 1)) then
                    Vso_mat(row,col) = VSO_off_diag_elem(s1,s2,so_ws,alpha_z,alpha_perp)
                endif        
            end do
        end do

        do row = 1,numstates
            do col = row +1, numstates
                Vso_mat(row, col) = Vso_mat(col, row)
            end do
        end do
        Vso_mat = Vso_mat*kappa
    end function


    pure function commutator(A1, A2)
        real(r_kind), dimension(:,:), intent(in) ::  A1, A2
        real(r_kind), dimension(size(A1,1),size(A1,2)) :: commutator
        commutator = matmul(A1,A2)-matmul(A2,A1)
    end function

    pure real(r_kind) elemental function S0(eta,s1,s2)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: eta
        S0 = (gnl(eta, s1%nr, s1%ml) * gnlp(eta, s2%nr, s2%ml) +gnl(eta, s2%nr, s2%ml) * gnlp(eta, s1%nr, s1%ml) ) / sqrt(eta)
    end function

    pure real(r_kind) elemental function Sp(eta,s1,s2)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: eta
        real(r_kind) :: K
        K = s2%ml + s2%ms*0.5_r_kind
        Sp = - gnl(eta,s1%nr,s1%ml)*gnl(eta,s2%nr,s2%ml) *(K+0.5_r_kind) / sqrt(eta) - gnlp(eta,s1%nr, s1%ml) * gnl(eta,s2%nr,s2%ml)
    end function

    pure real(r_kind) elemental function Sm(eta,s1,s2)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: eta
        Sm = Sp(eta, s2, s1)
    end function

    real(r_kind) function T0(eta,s1,s2, SO_WS,alpha_z, alpha_perp)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: eta, alpha_z, alpha_perp
        type(WS_pot), intent(in) :: SO_WS
        integer :: i
        real(r_kind) :: zeta, z, rho,r ,theta
        T0 = 0
        call precompute()
        rho = eta_to_rho(eta, alpha_perp)
        do i = 1,gauss_order
            zeta = her_x(i)
            z = zeta_to_z(zeta, alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho, z)


            T0 = T0 + her_w(i)*SO_WS%eval(r,theta) * Hmn(zeta, s1%nz) * Hmn(zeta, s2%nz)
        end do
    end function

    real(r_kind) function Tminus(eta,s1,s2, SO_WS,alpha_z, alpha_perp)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: eta, alpha_z, alpha_perp
        type(WS_pot), intent(in) :: SO_WS
        integer :: i
        real(r_kind) :: zeta, z, rho,r ,theta
        Tminus = 0
        call precompute()
        rho = eta_to_rho(eta, alpha_perp)
        do i = 1,gauss_order
            zeta = her_x(i)
            z = zeta_to_z(zeta, alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho, z)


            Tminus = Tminus + her_w(i)*SO_WS%eval(r,theta) * Hmn(zeta, s1%nz) * Hmnp(zeta, s2%nz)
        end do
    end function

    real(r_kind) function Tplus(eta,s1,s2, SO_WS,alpha_z, alpha_perp)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: eta, alpha_z, alpha_perp
        type(WS_pot), intent(in) :: SO_WS
        Tplus = Tminus(eta, s2, s1, SO_WS, alpha_z, alpha_perp) !!change of indices·
    end function



    pure elemental real(r_kind) function zeta_to_z(zeta, alpha)
        real(r_kind), intent(in) :: zeta, alpha
        zeta_to_z = zeta / alpha
    end function

    pure elemental real(r_kind) function eta_to_rho(eta, alpha)
        real(r_kind), intent(in) :: eta, alpha
        eta_to_rho = sqrt(eta) / alpha
    end function

    

    real(r_kind) function VWS_elem(s1,s2,V_WS,alpha_z,alpha_perp) result(matelem)
            type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: alpha_z, alpha_perp
        type(WS_pot), intent(in) :: V_WS
        integer :: ii
        real(r_kind) :: eta
        call precompute()
        matelem = 0.0_r_kind
        if(s1%ml /= s2%ml .or. s1%ms /= s2%ms) then
            matelem = 0.0_r_kind
            return
        endif

        do ii = 1, gauss_order
            eta = lag_x(ii)
            matelem = matelem + lag_w(ii) * Z_matelem(eta, s1,s2,V_WS,alpha_z,alpha_perp) * gnl(eta, s1%nr,s1%ml) * gnl(eta, s2%nr, s2%ml) 
        end do

    end function

    real(r_kind) function Z_matelem(eta, s1,s2,V_WS,alpha_z, alpha_perp)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: alpha_z, alpha_perp, eta
        type(WS_pot), intent(in) :: V_WS
        integer :: ii
        real(r_kind) :: zeta, z, rho, r, theta
        Z_matelem = 0.0_r_kind
        rho = eta_to_rho(eta, alpha_perp)
        do ii = 1, gauss_order
            zeta = her_x(ii)
            z = zeta_to_z(zeta,alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho ,z)
            Z_matelem = Z_matelem + her_w(ii) * Hmn(zeta, s1%nz) * Hmn(zeta, s2%nz) * V_WS%eval(r, theta) 
        end do
    end function


    real(r_kind)  function VSO_diag_elem(s1,s2,SO_WS,alpha_z,alpha_perp) result(matelem)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: alpha_z, alpha_perp
        type(WS_pot), intent(in) :: SO_WS
        integer :: ii
        real(r_kind) :: eta
        matelem = 0.0_r_kind
        if(s1%ml /= s2%ml .and. s1%ms /= s2%ms) then
            write(*,*) "Error: states are non diagonal in diagonal element calculation"
            call exit
        endif

        do ii = 1, gauss_order
            eta = lag_x(ii)
            matelem = matelem + S0(eta, s1, s2) * T0(eta, s1, s2, SO_WS, alpha_z, alpha_perp) * lag_w(ii)*alpha_perp* s1%ml * s1%ms !== 2*lambda*sigma where sigma = +- 1/2
        end do
    end function

    real(r_kind)  function VSO_off_diag_elem(s1,s2,SO_WS,alpha_z,alpha_perp) result(matelem)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: alpha_z, alpha_perp
        type(WS_pot), intent(in) :: SO_WS
        integer :: ii
        real(r_kind) :: eta
        matelem = 0

        if(s1%ml == s2%ml + 1 .and. s1%ms == s2%ms - 1) then

        elseif( s1%ml == s2%ml - 1 .and. s1%ms == s2%ms + 1)then

        else
            write(*,*) "Error: states dont have the correct quantum numbers in ms and ml"
            call exit
        endif

        do ii = 1, gauss_order
            eta = lag_x(ii)
            matelem = matelem + lag_w(ii) * (Sp(eta, s1, s2) * Tplus(eta, s1, s2, SO_WS, alpha_z, alpha_perp) + &
                                            Sm(eta, s1, s2) * Tminus(eta, s1, s2, SO_WS, alpha_z, alpha_perp))
        end do
        matelem = matelem * alpha_perp*alpha_z
    end function
end module