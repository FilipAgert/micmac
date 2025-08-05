module Hamiltonian
    use constants
    use def_ho
    use quad
    !!use quadrule, only:legendre_dr_compute, hermite_ek_compute, laguerre_ss_compute
    use optimise, only:func_1d, conj_grad_method
    implicit none


    private
    public :: diagonalize, WS_pot, VC_pot, dist_min, Vso_mat, commutator, VWS_mat, coul_mat, print_levels, H_protons, H_neutrons, get_levels, surfRadius
    public :: S0, T0, Tplus, Tminus, VSO_off_diag_elem_v2, el_pot, print_shell_params, write_result, time_diag, time_ws, time_Vc, time_Vso, surface_elem
    integer, parameter :: nquad =8
    real(kind), dimension(nquad) :: her_x, her_w, lag_x, lag_w, leg_x, leg_w
    logical :: precomputed_quad = .false.

    real(kind) :: time_diag=0, time_ws=0, time_Vc=0, time_Vso=0






    type, extends(func_1d) :: dist_min
        real(kind) :: r, theta, radius! radius normalized as r/R0
        type(betadef) :: def
    contains
        procedure :: eval => surfdisteval

    end type



    type, abstract :: potential
        type(betadef) :: def
        real(kind) :: radius
        real(kind), dimension(nquad, nquad) :: precomputed_pot
        logical, dimension(nquad, nquad) :: is_computed = .false.
    contains
        procedure(eval), deferred :: eval
        procedure :: set => set_precomp
        procedure :: get => get_precomp
        procedure :: eval_pre => eval_precomp
    end type
    !!!!!!!!
    type, extends(potential) :: WS_pot
        real(kind) :: V_0
    contains
        procedure :: eval => eval_WS
    end type


    type, extends(potential) :: VC_pot
        real(kind) :: charge_dens
    contains
        procedure :: eval => eval_VC
        procedure :: set_charge_dens => set_charge_VC
    end type

    abstract interface
        real(kind) function eval(self, r, theta)
            import:: potential, kind
            class(potential), intent(in) :: self
            real(kind), intent(in) :: r, theta
        end function
    end interface


    contains
    real(kind) function eval_precomp(self, ir, it, r, theta)
        class(potential) :: self
        integer, intent(in) :: ir, it
        real(kind), intent(in) :: r, theta

        if(self%is_computed(ir, it)) then
            eval_precomp = self%precomputed_pot(ir,it)
        else
            eval_precomp = self%eval(r, theta)
            self%precomputed_pot(ir, it) = eval_precomp
            self%is_computed(ir,it) = .true.
        endif

    end function
    pure real(kind) function get_precomp(self, ir, it)
        class(potential), intent(in) :: self
        integer, intent(in) :: ir, it
        get_precomp = self%precomputed_pot(ir, it)
    end function

    subroutine set_precomp(self, ir, it, val)
        class(potential) :: self
        integer, intent(in) :: ir, it
        real(kind), intent(in) :: val
        self%precomputed_pot(ir, it) = val
    end subroutine



    subroutine precompute_quad()
        if(precomputed_quad) then
            return
        endif


        !!64 point gaussian quadrature.
        call herquad(her_x, her_w)
        !call hermite_ek_compute(gauss_order,her_x,her_w)
        !!Gets weights and locations for where to evaluate 64 point integral. 
        !!Integral must be of form 
        !! -infty < x < infty  dx f(x) * exp(-x^2)
        !! where x = alpha_z*z => z = x/alpha_z

        call LAGQUAD(lag_x, lag_w)
        !!call laguerre_ss_compute(gauss_order, lag_x, lag_w)
        !!Gets weights and locations for where to evaluate 64 point integral. 
        !!Integral must be of form 
        !! 0 < x < infty  dx f(x) * exp(-x)
        !! where x = alpha^2 * rho^2

        call LEGQUAD(leg_x, leg_w)
        !!call legendre_dr_compute(gauss_order, leg_x, leg_w)

        precomputed_quad = .true.




    end subroutine

    subroutine set_charge_VC(self, ZZ)
        class(VC_pot) :: self
        integer, intent(in) :: ZZ !!proton number
        if(self%radius < 1e-6)then
            print*, "error: radius not set correctly for coulomb potential"
            call exit
        endif
        self%charge_dens = coul_charge_density(ZZ, self%radius)
    end subroutine

    pure real(kind) function coul_charge_density(ZZ, radius) result(rho)
        integer, intent(in) :: ZZ !!proton number
        real(kind), intent(in) :: radius
        rho = (ZZ-1) *3.0_kind/(4.0_kind * pi*4*pi * epsilonzero * radius**3)  !! Charge /( Volume * 4piepsilonzero)

    end function

    subroutine diagonalize(E, V, H) !!Diagnoalize H and return V and E as eigenvectors sorted with lowest energy first.
        real(kind), intent(in) :: H(:,:)
        real(kind), intent(out) :: E(size(H,1)), V(size(H,1),size(H,1))
        external :: dsyev !!lapack routine for diagnoalizing symmetric matrix

        !!lapack variables!!
        character(len=1), parameter :: jobz = 'V', uplo = 'U'
        integer :: N
        integer :: lda
        real(kind) :: w(size(H,1))
        real(kind), allocatable :: work(:)
        integer, allocatable :: iwork(:)
        integer :: liwork
        integer :: lwork
        integer :: info
        real(kind) :: time0, time1
        call cpu_time(time0)
        N = size(H,1)
        ! write(*,*) N
        lda = N
        lwork = (1+6*N+2*N*N)*2
        liwork = 3+5*N
        allocate(work(lwork), iwork(liwork))
        !call dsyev(jobz, uplo, N, H, lda, w, work, lwork, info)
        call dsyevd(jobz,uplo, N, H, lda, w, work, lwork, iwork, liwork,info)

        if (info /= 0) then
            write(*,'(A,I10)') "Error using dsyev in diagonalisation, info:", info
            call exit
        endif
        if(W(1) > W(size(E))) then !Reverse order if wrong order
            E = W(size(W):1:-1) 
            V = H(:,size(W):1:-1)
        else
            E = w
            V = H
        endif
        deallocate(work, iwork)
        call cpu_time(time1)
        time_diag = time_diag + time1-time0
    end subroutine

    pure real(kind) function el_pot(r,rad,charge_dens) !!electric potential of a uniformly charged sphere
        real(kind),intent(in) :: r, rad, charge_dens
        if(r .le. rad) then !!inside nucleus
            el_pot = 4.0_kind*pi/2.0_kind * (rad**2 - r**2/3.0_kind)!4(3.0-(r/rad)**2)/(2.0*rad)
        else!outside nucleus
            el_pot = 4.0_kind*pi/3.0_kind * rad**3/r
        endif
        el_pot = el_pot * charge_dens
    end function

    function surface_elem(def, theta, phi, r0) result(out)
        !!Non-normalized Normal vector of surface of nucleus at theta,phi in cartesian coordinates
        !!is Jacobian as well. so if v is this vector, dS = n dA = v dtheta dphi
        type(betadef), intent(in) :: def
        real(kind), intent(in) :: phi, theta
        real(kind), dimension(3) :: out
        real(kind), dimension(3) :: r, rt, rp
        real(kind), intent(in) :: r0
        real(kind) :: st, ct, sp, cp
        real(kind) :: rad, draddt
        st = sin(theta)
        ct = cos(theta)
        sp = sin(phi)
        cp = cos(phi)
        r = [st*cp, st*sp, ct]
        rp = [-st*sp, st*cp,0.0_kind]
        rt = [ct*cp, ct*sp, -st]
        rad = surfRadius(theta, def, r0)
        draddt = dsurfRadiusdtheta(theta, def, r0)

        out = rad*draddt * cross(r,rp) + rad*rad * cross(rt,rp)
    end function

    pure function cross(v1,v2)
        !!computes cross product in cartesian coordinates
        real(kind), dimension(3), intent(in) :: v1,v2
        real(kind), dimension(3):: cross
        cross(1) = v1(2)*v2(3) - v1(3)*v2(2)
        cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
        cross(3) = v1(1)*v2(2) - v1(2)*v2(1)
    end function

    real(kind) function eval_vc(self,r, theta) !!axially symmetric. function of z and rho
        class(VC_pot), intent(in) :: self
        real(kind), intent(in) :: r, theta
        !Do the integral int_V     dr^3/|r-r'| at the point r.
        integer :: iu, it, ix
        real(kind) :: u, t, x, int_theta, int_radius_ub, rolling, sintheta, sinacosu, cospipit, costheta, vinc(3), vint(3), diff(3), dist, nj(3)
        real(kind) :: phi, thetap, sinp, cosp, sint, cost, radp
        logical, save :: first_time = .true.
        real(kind), save, dimension(nquad) :: sinthetas, costhetas, sinphis, cosphis, surfrads
        real(kind), save, dimension(3,nquad,nquad) :: normJac

        if(self%def%eq(spherical_def) )then !!if spherical, use analytical formula
            eval_vc = el_pot(r, self%radius, self%charge_dens)
           return
        endif
        call precompute_quad()
        if(first_time) then

            do iu = 1, nquad
                u = leg_x(iu)
                thetap = (u + 1)*pi/2.0_kind
                sinthetas(iu) = sin(thetap)
                costhetas(iu) = cos(thetap)
                surfrads(iu) = surfRadius(thetap, self%def, self%radius)

                phi = (u+1)*pi
                sinphis(iu) = sin(phi)
                cosphis(iu) = cos(phi)

                do it = 1,nquad
                    t = leg_x(it)
                    phi = (t+1)*pi
                    normJac(:,it,iu) = surface_elem(self%def, thetap, phi, self%radius)
                end do
            end do

            first_time = .false.
        endif
        call precompute_quad()

        sintheta = sin(theta)
        costheta = cos(theta)
        vinc = r*[sin(theta), 0.0_kind, cos(theta)] ! input vector in cartesian coordinates
        rolling = 0.0_kind
        call omp_set_num_threads(num_threads)
        !//$omp parallel do private (iu, it, ix, u, int_theta, int_radius_ub, t, x) &
        !//$omp& reduction(+:rolling)
        !https://en.wikipedia.org/wiki/Green%27s_identities

        do iu = 1, nquad
            u = leg_x(iu)
            thetap = (u + 1)*pi/2.0_kind
            sint = sinthetas(iu)
            cost = costhetas(iu)
            radp = surfrads(iu)!surfRadius(thetap, self%def, self%radius)!!r for r'
            do it = 1,nquad
                t = leg_x(it)
                phi = (t+1)*pi
                sinp = sinphis(it)
                cosp = cosphis(it)
                vint = radp*[sint*cosp,sint*sinp,cost] 
                diff = vinc - vint
                dist = sqrt(dot_product(diff,diff))
                nj = normJac(:,it,iu)!surface_elem(self%def, thetap, phi, self%radius) !!normal of surface at integration point
                rolling = rolling - dot_product(diff, nj)/dist * leg_w(iu)*leg_w(it)
            end do
        end do
        eval_vc = rolling * self%charge_dens*(pi/2.0_kind * pi)/2 !from variable transform. last /2 unknown origin but makes result good.
        ! write(*,*) eval_vc
        ! pause
        ! write(*,'(a15,f12.6)') "new method:", eval_vc
        ! rolling = 0
        ! do iu = 1, nquad
        !     u = leg_x(iu)
        !     int_theta = acos(u)
        !     int_radius_ub = surfRadius(int_theta, self%def, self%radius)
        !     sinacosu = sin(acos(u))
        !     do it = 1, nquad
        !         t = leg_x(it)
        !         cospipit = cos(pi + pi*t)

        !         do ix = 1,nquad
        !             x = leg_x(ix)

        !             rolling = rolling + leg_w(iu)*leg_w(it)*leg_w(ix)*pi*int_radius_ub**3*(x+1.0_kind)**2 / (8.0_kind) &
        !                     / sqrt(r**2 + int_radius_ub**2 * (x+1.0_kind)**2 / 4.0_kind &
        !                             + r*int_radius_ub*(x+1.0_kind)* (sintheta*sinacosu*cospipit + u*costheta))
                    
        !             ! if(isnan(rolling)) then
        !             !     print*, "isnan in coul: iu, it, ix:"
        !             !     print*, iu, it, ix
        !             ! endif
        !         end do
        !     end do
        ! end do
        ! !//$omp end parallel do
        ! eval_vc = rolling * self%charge_dens
        ! write(*,'(a15,f12.6)') "old method:", eval_vc
    end function





    pure elemental real(kind) function Y20(theta)
        real(kind), intent(in) :: theta
        Y20 = sqrt(5.0_kind / (16 * pi)) * (3.0_kind * cos(theta)**2 - 1.0_kind)
    end function

    pure elemental real(kind) function Y40(theta)
        real(kind) ,intent(in) :: theta
        Y40 = sqrt(9.0_kind / (256 * pi)) * (35.0_kind * cos(theta)**4 - 30.0_kind * cos(theta)**2 + 3)
    end function

    pure elemental real(kind) function dY20dtheta(theta) result(dydt)
        real(kind), intent(in) :: theta
        dydt = -1.5_kind * sqrt(5.0_kind/pi) * sin(theta)*cos(theta)
    end function

    pure elemental real(kind) function dY40dtheta(theta) result(dydt)
        real(kind), intent(in) :: theta
        dydt = 15*sin(theta)*cos(theta)*(3.0_kind-7.0_kind*cos(theta)*cos(theta))/(4.0_kind*sqrt(pi))
    end function

    pure elemental real(kind) function surfdisteval(self, x) !!Gets distance to surface at angle x.
        class(dist_min), intent(in) :: self
        real(kind), intent(in) :: x !!Angle theta in this case.

        real(kind) :: radius
        radius = surfRadius(x, self%def, self%radius)

        surfdisteval = (self%r*cos(self%theta) - radius*cos(x))**2 &
                    + (self%r*sin(self%theta) - radius*sin(x))**2
    end function

    pure elemental real(kind) function surfRadius(theta, def, R0) !!computes distance to surface
        real(kind), intent(in) :: theta, R0
        type(betadef), intent(in) :: def
        surfRadius = R0 * (1.0_kind + def%beta2 * Y20(theta) + def%beta4*Y40(theta))
    end function

    pure elemental real(kind) function dsurfRadiusdtheta(theta, def, R0)result(dRdt) !!partial derivative of radius wrt theta.
        real(kind), intent(in) :: theta, R0
        type(betadef), intent(in) :: def
        drdt = R0 * (def%beta2 * dY20dtheta(theta) + def%beta4*dY40dtheta(theta))
    end function


    real(kind) pure elemental function theta_cyl(rho,z)  !!Coordinate transform from cylindrical coordinates to spherical theta
        real(kind), intent(in):: rho,z
        if (z < 0) then
            theta_cyl = ATAN(rho/z) + pi
        else
            theta_cyl = ATAN(rho/z)
        endif
    end function

    real(kind) pure elemental function rad_cyl(rho,z)  !!Coordinate transform from cylindrical coordinates to spherical radius
        real(kind), intent(in):: rho,z
        rad_cyl = sqrt(rho*rho+z*z)
    end function

    real(kind) function eval_ws(self,r, theta) !!axially symmetric. function of z and rho
        real(kind), intent(in) :: r, theta
        class(WS_pot), intent(in) :: self
        real(kind) ::dist
        dist = sqrt(surfdist(r,theta,self%def, self%radius)) 
        if(r < surfRadius(theta, self%def,self%radius)) then
            dist = - dist
        endif
        eval_ws = -self%V_0/ (1 + exp(dist/a_ws))
    end function


    real(kind) function surfdist(r,theta,def, radius) result(dist) !!Finds shortest distance to surface of nucleus.
        real(kind), intent(in) :: r, theta, radius
        type(betadef), intent(in) :: def
        type(dist_min) :: distfunc
        real(kind) :: ang
        logical :: conv

        distfunc = dist_min(r=r, theta=theta, def=def, radius=radius)
        call conj_grad_method(ang, dist, conv, distfunc, 0.0_kind, pi,theta,1e-5_kind)
    end function

    function coul_mat(states, def, R0,ZZ, mass, omegaz, omegaperp)
        implicit none
        type(an_ho_state), intent(in) :: states(:)
        type(betadef), intent(in) :: def !!deformation of body
        real(kind), intent(in) :: R0 !!radius as r0*A**(1/3). Where r0 is r0_so from constants.f90
        real(kind), dimension(size(states),size(states)) :: coul_mat
        integer, intent(in) :: ZZ!!proton number
        real(kind) :: alpha_z, alpha_perp
        type(VC_pot) :: cp
        real(kind), intent(in) :: omegaz, omegaperp !!Frequencies in units MeV/hbar
        real(kind), intent(in) :: mass!!Frequencies in units MeV/c^2
        real(kind), allocatable, dimension(:,:,:) :: Z_mat
        integer ::numstates, row, col
        type(an_ho_state) :: s1, s2
        real(kind) :: time0, time1
        allocate(Z_mat(nquad,0:max_n,0:max_n))
        call precompute_quad()
        call cpu_time(time0)

        !setup potential

        cp%def = def
        cp%radius = R0
        call cp%set_charge_dens(ZZ)
        
        Z_mat = 0
        ! printflag = .true.
        alpha_z = alpha(mass, omegaz)
        alpha_perp = alpha(mass, omegaperp)

        call comp_zmat(Z_mat, cp, alpha_z, alpha_perp)

        numstates = size(states)
        coul_mat = 0
        do row = 1, numstates
            s1 = states(row)
            do col = 1,row
                s2 = states(col)
                coul_mat(row,col) = pot_elem_zmatpre(s1,s2,Z_mat)!

            end do
        end do

        do row = 1,numstates
            do col = row +1, numstates
                coul_mat(row, col) = coul_mat(col, row)
            end do
        end do

        if(hasnan2d(coul_mat)) then
            write(*,*) "Error: Coul Matrix has a NaN"
        endif
        call cpu_time(time1)
        time_Vc = time_Vc + time1-time0
    end function

    subroutine comp_zmat(Zmat, pot,alpha_z, alpha_perp) !!use the reccurence relations to compute the Zmat completely from the main diagonals
        real(kind), dimension(nquad,0:max_n,0:max_n), intent(out) :: Zmat
        logical, dimension(0:max_n, 0:max_n) :: computed
        class(potential) :: pot
        real(kind), intent(in) :: alpha_z, alpha_perp
        real(kind) :: eta
        integer :: nz, nzp, ii, diag

        call precompute_quad()
        Zmat = 0
        computed = .false.

        do nz = 0, max_n !!compute main diagonals
            do nzp = nz, min(max_n, nz+1)
                do ii = 1, nquad
                    eta = lag_x(ii)
                    Zmat(ii,nzp, nz) = Z_matelem(eta, ii,nzp,nz,pot,alpha_z,alpha_perp)

                end do
                computed(nzp, nz) = .true.
            end do
        end do

        do nzp = 1,max_n !reflect off diagonal.
            nz = nzp -1
            Zmat(:,nzp, nz) = Zmat(:,nz, nzp)
            computed(nzp, nz) = .true.
        end do
        ! write(str, '(A,I2,A)')'(',max_n+1,'F10.3)'
        ! write(*,*) "After main diagonals:"
        ! do nz = 0, max_n
        !     write(*, str) Zmat(nz,:,1)
        ! end do

        !!compute off diagonals using recurrence relation:
        do diag = 2, max_n
            do nz = diag, max_n
                nzp = nz - diag

                if(computed(nzp, nz)) cycle


                if(within_bounds(nzp + 1,0, max_n) .and. within_bounds(nz - 1, 0, max_n)) then
                    if(.not. computed(nzp + 1, nz - 1)) then
                        do ii = 1, nquad
                            eta = lag_x(ii)
                            Zmat(ii,nzp + 1, nz - 1) = Z_matelem(eta, ii,nzp + 1,nz - 1,pot,alpha_z,alpha_perp)
                            
                        end do
                        computed(nzp + 1, nz - 1) = .true.
                    endif


                    Zmat(:,nzp, nz) = Zmat(:,nzp, nz) + Zmat(:,nzp+1, nz-1) * sqrt((nzp+1.0_kind)/nz)
                    if(nz==0) then
                        write(*,*) "ERR: divide by zero nz, nzp: ", nz, nzp
                    endif
                endif

                if(within_bounds(nzp - 1,0, max_n) .and. within_bounds(nz - 1, 0, max_n)) then
                    if(.not. computed(nzp - 1, nz - 1)) then
                        do ii = 1, nquad
                            eta = lag_x(ii)
                            Zmat(ii,nzp - 1, nz - 1) = Z_matelem(eta, ii,nzp - 1,nz - 1,pot,alpha_z,alpha_perp)
        
                        end do
                        computed(nzp - 1, nz - 1) = .true.
                    endif


                    Zmat(:,nzp, nz) = Zmat(:,nzp, nz) + Zmat(:,nzp-1, nz-1) * sqrt(real(nzp,kind)/nz)
                    if(nz==0) then
                        write(*,*) "ERR: divide by zero nz, nzp: ", nz, nzp
                    endif
                endif

                if(within_bounds(nzp,0, max_n) .and. within_bounds(nz - 2, 0, max_n)) then
                    if(.not. computed(nzp, nz - 2)) then
                        do ii = 1, nquad
                            eta = lag_x(ii)
                            Zmat(ii,nzp, nz - 2) = Z_matelem(eta, ii,nzp,nz - 2,pot,alpha_z,alpha_perp)
        
                        end do
                        computed(nzp, nz - 2) = .true.
                    endif


                    Zmat(:,nzp, nz) = Zmat(:,nzp, nz) - Zmat(:,nzp, nz-2) * sqrt((nz-1.0_kind)/nz)
                    if(nz==0) then
                        write(*,*) "ERR: divide by zero: nz, nzp: ", nz, nzp
                    endif
                endif
                computed(nzp, nz) = .true.
            end do

        end do
        ! write(*,*) "after recurrence:"
        ! do nz = 0, max_n
        !     write(*, str) Zmat(nz,:,1)
        ! end do
        do nz = 0, max_n !reflect to lower triangular matrix.
            do nzp = nz, max_n
                Zmat(:,nzp, nz)= Zmat(:,nz, nzp)
            end do
        end do
        ! write(*,*) "after reflection:"
        ! do nz = 0, max_n
        !     write(*, str) Zmat(nz,:,1)
        ! end do

        if(hasnan3d(Zmat)) then
            write(*,*) "Error: Zmat Matrix has a NaN"
        endif
    end subroutine

    pure logical function within_bounds(val, lb, ub)
        integer, intent(in) :: val, lb, ub
        within_bounds = val .ge. lb .and. val .le. ub
    end function

    function VWS_mat(states, def, R0, omegaz, omegaperp, mass, WS_depth)
        implicit none
        type(an_ho_state), intent(in) :: states(:)
        type(betadef), intent(in) :: def !!deformation of body
        real(kind), intent(in) :: R0 !!radius as r0*A**(1/3). Where r0 is r0_so from constants.f90
        real(kind), intent(in) :: WS_depth !!WS potental well depth
        real(kind), dimension(size(states),size(states)) :: VWS_mat
        real(kind), intent(in) :: omegaz, omegaperp !!Frequencies in units MeV/hbar
        real(kind), intent(in) :: mass!!Frequencies in units MeV/c^2
        real(kind) :: alpha_z, alpha_perp
        real(kind), allocatable, dimension(:,:,:) :: Z_mat
        type(WS_pot) :: VWS
        integer ::numstates, row, col
        type(an_ho_state) :: s1, s2
        real(kind) :: time0, time1


        call precompute_quad()
        call cpu_time(time0)
        allocate(Z_MAT(nquad,0:max_n,0:max_n))
        !setup potential
        VWS%def = def
        VWS%radius = R0
        VWS%V_0 = WS_depth

        alpha_z = alpha(mass, omegaz)
        alpha_perp = alpha(mass, omegaperp)
        numstates = size(states)
        VWS_mat = 0
        call omp_set_num_threads(num_threads)
        call comp_zmat(Z_mat, VWS, alpha_z, alpha_perp)

        !//$omp parallel shared(states,VWS,alpha_z,alpha_perp, VWS_mat) private(s1,s2)
        !//$omp do schedule(static)
        do row = 1, numstates
            s1 = states(row)
            do col = 1,row
                s2 = states(col)
                VWS_mat(row,col) = pot_elem_zmatpre(s1,s2,Z_mat)!
                ! if(isnan(VWS_mat(row,col))) then
                !     write(*,*) "error: ", row, col, " is nan"
                ! endif
            end do
        end do
        !//$omp end do
        !//$omp end parallel
        do row = 1,numstates
            do col = row +1, numstates
                VWS_mat(row, col) = VWS_mat(col, row)
                ! if(isnan(VWS_mat(row,col))) then
                !     write(*,*) "error: ", row, col, " is nan"
                ! endif
            end do
        end do
        call cpu_time(time1)
        time_ws = time_ws + time1-time0
        !call write_mat_to_file(VWS_mat)
    end function

    function Vso_mat(states, def, R0, omegaz, omegaperp, mass, WS_depth, lambda)
        implicit none
        type(an_ho_state), intent(in) :: states(:)
        type(betadef), intent(in) :: def !!deformation of body
        real(kind), intent(in) :: R0 !!radius as r0*A**(1/3). Where r0 is r0_so from constants.f90
        real(kind), intent(in) :: WS_depth !!WS potental well depth
        real(kind), dimension(size(states),size(states)) :: Vso_mat
        real(kind), intent(in) :: omegaz, omegaperp !!Frequencies in units MeV/hbar
        real(kind), intent(in) :: mass!!Frequencies in units MeV/c^2
        real(kind), intent(in) ::lambda !!multiplicative factor
        real(kind) :: alpha_z, alpha_perp, kappa
        real(kind), allocatable, dimension(:,:,:) :: T0, W0
        type(WS_pot) :: so_ws
        real(kind) :: time0, time1
        integer ::numstates, row, col, l1, l2, ms1, ms2
        type(an_ho_state) :: s1, s2
        allocate(T0(nquad,0:max_n, 0:max_n), W0(nquad,0:max_n, 0:max_n))
        call precompute_quad()
        call cpu_time(time0)
        kappa = lambda*(hbarc/(2*mass))**2
        ! print*, "kappa:", kappa
        !setup potential
        so_ws%def = def
        so_ws%radius = R0
        so_ws%V_0 = WS_depth

        alpha_z = alpha(mass, omegaz)
        alpha_perp = alpha(mass, omegaperp)

        numstates = size(states)
        Vso_mat = 0
        call omp_set_num_threads(num_threads)
        T0 = T0_mat(So_ws, alpha_z, alpha_perp)
        W0 = W0_mat(so_ws, alpha_z, alpha_perp)
        !//$omp parallel shared(states,Vso_mat,alpha_z,alpha_perp,so_ws) private(s1,l1,ms1,s2,l2,ms2)
        !//$omp do schedule(static)
        do row = 1, numstates
            s1 = states(row)
            l1 = s1%ml
            ms1 = s1%ms
            do col = 1,row
                s2 = states(col)
                l2 = s2%ml
                ms2 = s2%ms
                
                if(ms1 == ms2 .and. l1 == l2) then
                    Vso_mat(row,col) = VSO_diag_elem(s1,s2,T0,alpha_perp)
                    ! if(isnan(Vso_mat(row,col))) then
                    !     write(*,*) "diagonal VSO element is not a number:"
                    !     call exit 
                    ! endif
                elseif((ms1 == ms2 + 2 .and. l1 == l2 - 1) .or. (ms1 == ms2 - 2 .and. l1 == l2 + 1)) then
                    Vso_mat(row,col) = VSO_off_diag_elem_v2(s1,s2,W0,alpha_z,alpha_perp)
                    ! if(isnan(Vso_mat(row,col))) then
                    !     write(*,*) "Off diagonal VSO element is not a number:"
                    !     call exit 
                    ! endif
                endif

            end do
        end do
        !//$omp end do
        !//$omp end parallel
        do row = 1,numstates
            do col = row +1, numstates
                Vso_mat(row, col) = Vso_mat(col, row)
            end do
        end do
        Vso_mat = Vso_mat*kappa
        call cpu_time(time1)
        time_Vso = time_Vso + time1-time0
        !call write_mat_to_file(Vso_mat)
    end function


    pure function commutator(A1, A2)
        real(kind), dimension(:,:), intent(in) ::  A1, A2
        real(kind), dimension(size(A1,1),size(A1,2)) :: commutator
        commutator = matmul(A1,A2)-matmul(A2,A1)
    end function


    real(kind) function S0(ii,s1,s2)
        type(an_ho_state), intent(in) :: s1, s2
        integer, intent(in) :: ii
        real(kind), save :: invsqrteta(1:nquad)
        logical, save :: first_time = .true.
        integer :: nn
        real(kind) :: eta

        if(first_time) then
            call precompute_quad()
            do nn = 1, nquad
                eta = lag_x(nn)
                invsqrteta(nn) = 1.0_kind/sqrt(eta)
            end do
            first_time = .false.
        endif
        S0 = (get_quad_gnl(ii, s1%nr, s1%ml) * get_quad_gnlp(ii, s2%nr, s2%ml) +get_quad_gnl(ii, s2%nr, s2%ml) * get_quad_gnlp(ii, s1%nr, s1%ml) ) *invsqrteta(ii)
    end function

    function T0_mat(SO_WS, alpha_z, alpha_perp)
        real(kind), dimension(nquad,0:max_n, 0:max_n) :: T0_mat
        type(WS_pot), intent(in) :: SO_WS
        real(kind), intent(in) :: alpha_z, alpha_perp

        call comp_zmat(T0_mat, SO_WS, alpha_z, alpha_perp)
    end function

    function W0_mat(SO_WS, alpha_z, alpha_perp) !!int tilde(h) h S
        !!TOdo: use recurence relation to not have to compute entire matrix
        real(kind), dimension(nquad,0:max_n, 0:max_n) :: W0_mat
        type(WS_pot), intent(in) :: SO_WS
        real(kind), intent(in) :: alpha_z, alpha_perp
        real(kind) :: eta
        integer :: nz, nzp, ii
        W0_mat = 0
        call precompute_quad()
        do nz = 0, max_n
            do nzp = 0, max_n
                do ii = 1, nquad
                    eta = lag_x(ii)
                    W0_mat(ii,nzp, nz) = W_matelem(eta, ii, nzp, nz, SO_WS, alpha_z, alpha_perp)
                end do
            end do
        end do


    end function

    real(kind) function T0(eta,eta_ii,nz1,nz2, SO_WS,alpha_z, alpha_perp)
        integer, intent(in) :: nz1, nz2
        real(kind), intent(in) :: eta, alpha_z, alpha_perp
        type(WS_pot), intent(in) :: SO_WS
        integer :: i
        integer, intent(in) :: eta_ii
        real(kind) :: zeta, z, rho,r ,theta
        T0 = 0
        call precompute_quad()
        rho = eta_to_rho(eta, alpha_perp)
        do i = 1,nquad
            zeta = her_x(i)
            z = zeta_to_z(zeta, alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho, z)


            T0 = T0 + her_w(i)*SO_WS%eval_pre(eta_ii,i,r,theta) * Hmn(zeta, nz1) * Hmn(zeta, nz2)
        end do
    end function

    real(kind) function Tminus(eta,eta_ii,s1,s2, SO_WS,alpha_z, alpha_perp)
        type(an_ho_state), intent(in) :: s1, s2
        real(kind), intent(in) :: eta, alpha_z, alpha_perp
        type(WS_pot), intent(in) :: SO_WS
        integer, intent(in) :: eta_ii
        integer :: i
        real(kind) :: zeta, z, rho,r ,theta
        Tminus = 0
        call precompute_quad()
        rho = eta_to_rho(eta, alpha_perp)
        do i = 1,nquad
            zeta = her_x(i)
            z = zeta_to_z(zeta, alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho, z)


            Tminus = Tminus + her_w(i)*SO_WS%eval_pre(eta_ii,i,r,theta) * Hmn(zeta, s1%nz) * Hmnp(zeta, s2%nz)
        end do
    end function

    real(kind) function Tplus(eta,eta_ii,s1,s2, SO_WS,alpha_z, alpha_perp)
        type(an_ho_state), intent(in) :: s1, s2
        real(kind), intent(in) :: eta, alpha_z, alpha_perp
        type(WS_pot), intent(in) :: SO_WS
        integer, intent(in) :: eta_ii
        integer :: i
        real(kind) :: zeta, z, rho,r ,theta
        Tplus = 0
        call precompute_quad()
        rho = eta_to_rho(eta, alpha_perp)
        do i = 1,nquad
            zeta = her_x(i)
            z = zeta_to_z(zeta, alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho, z)


            Tplus = Tplus + her_w(i)*SO_WS%eval_pre(eta_ii,i,r,theta) * Hmn(zeta, s2%nz) * Hmnp(zeta, s1%nz)
        end do
    end function

    real(kind) function I1(s1,s2, W0,alpha_z, alpha_perp)
        type(an_ho_state), intent(in) :: s1, s2
        real(kind), intent(in) ::  alpha_z, alpha_perp
        real(kind), dimension(nquad,0:max_n, 0:max_n), intent(in) :: W0
        integer :: nn
        real(kind) ::rho, eta, weta, term
        real(kind), save :: sqrteta(1:nquad)
        logical, save :: first_time = .true.
        call precompute_quad()
        if(first_time) then
            do nn = 1, nquad
                eta = lag_x(nn)
                sqrteta(nn) = sqrt(eta)
            end do
            first_time = .false.
        endif



        I1 = 0
        do nn = 1, nquad
            eta = lag_x(nn)
            weta = lag_w(nn)
            term = 0
            rho = eta_to_rho(eta, alpha_perp)

            term = W0(nn,s1%nz, s2%nz) * s2%ml + W0(nn,s2%nz, s1%nz) * s1%ml
            term = term * weta * get_quad_gnl(nn,s1%nr, s1%ml) * get_quad_gnl(nn, s2%nr,s2%ml) / sqrteta(nn)
            I1 = I1 + term
        end do
        I1 = I1 * alpha_perp*alpha_z
        ! if(isnan(I1)) then
        !     write(*,*) "ERR: I1 is not a number"
        !     write(*,*) "W0:", W0(s1%nz, s2%nz, 1), W0(s2%nz, s1%nz, 1)
        !     write(*,*) "gnl: ",gnl(eta,s1%nr, s1%ml), gnl(eta, s2%nr,s2%ml)
        ! endif
    end function

    real(kind) function I2(s1,s2, W0,alpha_z, alpha_perp)
        type(an_ho_state), intent(in) :: s1, s2
        real(kind), intent(in) ::  alpha_z, alpha_perp
        integer :: nn
        real(kind), dimension(nquad,0:max_n, 0:max_n), intent(in) :: W0
        real(kind) :: rho, eta, weta, term1, term2,etafac
        real(kind), save :: invsqrteta(1:nquad)
        logical, save :: first_time = .true.
        call precompute_quad()
        if(first_time)then
            do nn = 1, nquad
                eta = lag_x(nn)
                invsqrteta(nn) = 1.0_kind/sqrt(eta)
            end do
            first_time = .false.
        endif

        I2 = 0
        do nn = 1, nquad
            eta = lag_x(nn)
            weta = lag_w(nn)
            term1 = 0
            term2 = 0
            rho = eta_to_rho(eta, alpha_perp)

            term1 = W0(nn,s1%nz, s2%nz)
            term2 = W0(nn,s2%nz, s1%nz)

            etafac = weta *invsqrteta(nn)
            term1 = term1 * etafac * get_quad_gnl(nn, s1%nr, s1%ml) * get_quad_gnlp(nn, s2%nr, s2%ml)
            term2 = term2 * etafac * get_quad_gnlp(nn, s1%nr, s1%ml) * get_quad_gnl(nn, s2%nr, s2%ml)

            I2 = I2 + term1 - term2
        end do
        I2 = I2 * alpha_perp*alpha_z

    end function


    pure elemental real(kind) function zeta_to_z(zeta, alpha)
        real(kind), intent(in) :: zeta, alpha
        zeta_to_z = zeta / alpha
    end function

    pure elemental real(kind) function eta_to_rho(eta, alpha)
        real(kind), intent(in) :: eta, alpha
        eta_to_rho = sqrt(eta) / alpha
    end function

    

    real(kind) function pot_elem_zmatpre(s1,s2,Zmat) result(matelem)
        type(an_ho_state), intent(in) :: s1, s2
        real(kind), intent(in), dimension(nquad,0:max_N, 0:max_N) :: Zmat
        integer :: ii
        real(kind) :: eta
        call precompute_quad()
        matelem = 0.0_kind
        if(s1%ml /= s2%ml .or. s1%ms /= s2%ms) then
            return
        endif

        do ii = 1, nquad
            eta = lag_x(ii)
            matelem = matelem + lag_w(ii) * Zmat(ii,s1%nz, s2%nz) * get_quad_gnl(ii, s1%nr,s1%ml) * get_quad_gnl(ii, s2%nr, s2%ml) 
        end do

    end function

    logical function is_reflection_sym(def) result(sym)
        type(betadef), intent(in) :: def
        sym = .true.

    end function

    real(kind) function Z_matelem(eta, eta_ii,nz1,nz2,pot,alpha_z, alpha_perp)
        integer, intent(in) :: nz1, nz2
        real(kind), intent(in) :: alpha_z, alpha_perp, eta
        class(potential), intent(in) :: pot
        integer :: ii, UB
        integer, intent(in) ::eta_ii
        real(kind) :: zeta, z, rho, r, theta
        Z_matelem = 0.0_kind
        rho = eta_to_rho(eta, alpha_perp)
        if(is_reflection_sym(pot%def)) then
            UB = nquad/2
            if(mod(nz1+nz2,2)==1)then
                 Z_matelem = 0
                return
            endif
        else
            UB = nquad
        endif

        do ii = 1, UB
            zeta = her_x(ii)
            z = zeta_to_z(zeta,alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho ,z)
     
            ! if(printflag) then
            !     write(*,'(A, 2F10.3)') 'z, rho', z, rho
            !     write(*,'(A, 2F10.3)') 'eta, zeta', eta, zeta
            !     write(*,'(A, 2F10.3)') 'alpha_z, alpha_perp', alpha_z, alpha_perp
            !     write(*,'(A, 2F10.3)') 'R, theta', r, theta
            !     write(*,'(A,F10.3,F10.3)'),"Eval vs preeval: " ,(pot%eval(r,theta)), pot%eval_pre(eta_ii, ii, r, theta)
            ! endif

            Z_matelem = Z_matelem + her_w(ii) * get_quad_Hmn(ii, nz1) * get_quad_Hmn(ii, nz2)*pot%eval_pre(eta_ii,ii,r, theta) 
            ! if(printflag) then
            !     write(*,'(F10.3)'),(pot%eval(0.0_kind,0.0_kind))
            ! endif

            
        end do
        write(*,*) "HALFWAY:" , Z_matelem, " dub:",  Z_matelem*2
        do ii = UB+1,nquad
            zeta = her_x(ii)
            z = zeta_to_z(zeta,alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho ,z)
     
            ! if(printflag) then
            !     write(*,'(A, 2F10.3)') 'z, rho', z, rho
            !     write(*,'(A, 2F10.3)') 'eta, zeta', eta, zeta
            !     write(*,'(A, 2F10.3)') 'alpha_z, alpha_perp', alpha_z, alpha_perp
            !     write(*,'(A, 2F10.3)') 'R, theta', r, theta
            !     write(*,'(A,F10.3,F10.3)'),"Eval vs preeval: " ,(pot%eval(r,theta)), pot%eval_pre(eta_ii, ii, r, theta)
            ! endif

            Z_matelem = Z_matelem + her_w(ii) * get_quad_Hmn(ii, nz1) * get_quad_Hmn(ii, nz2)*pot%eval_pre(eta_ii,ii,r, theta) 
            ! if(printflag) then
            !     write(*,'(F10.3)'),(pot%eval(0.0_kind,0.0_kind))
            ! endif

            
        end do
        write(*,*) "FIN:" , Z_matelem

        pause
        if(is_reflection_sym(pot%def)) then
            Z_matelem = Z_matelem*2
            return
            ! print*, "here"
            if(mod(nquad,2) == 1) then
                ii = UB + 1
                zeta = her_x(ii) !!should be zero
                if(zeta /= 0.0_kind) error stop "not zero"
                z = zeta_to_z(zeta,alpha_z)
                theta = theta_cyl(rho, z)
                r = rad_cyl(rho ,z)
                Z_matelem = Z_matelem + her_w(ii) * get_quad_Hmn(ii, nz1) * get_quad_Hmn(ii, nz2)*pot%eval_pre(eta_ii,ii,r, theta) 
            endif

        endif
    end function

    real(kind) function W_matelem(eta, eta_ii,nz1,nz2,pot,alpha_z, alpha_perp) !! W(eta) = int zeta    h(zeta) h'(zeta) S(eta, zeta)
        integer, intent(in) :: nz1, nz2
        real(kind), intent(in) :: alpha_z, alpha_perp, eta
        class(potential), intent(in) :: pot
        integer :: ii
        integer, intent(in) ::eta_ii
        real(kind) :: zeta, z, rho, r, theta
        W_matelem = 0.0_kind
        rho = eta_to_rho(eta, alpha_perp)
        do ii = 1, nquad
            zeta = her_x(ii)
            z = zeta_to_z(zeta,alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho ,z)
            W_matelem = W_matelem + her_w(ii) * get_quad_Hmnp(ii, nz1) *get_quad_Hmn(ii,nz2)*pot%eval_pre(eta_ii,ii,r, theta) 
            ! if(isnan(W_matelem)) then
            !     print*, "Err Wmatelem is Nan"
            !     print*, "hmnp:", Hmnp(zeta, nz1)
            !     print*,  "Hmn:", Hmn(zeta, nz2)
            !     print*, "pot:", pot%eval_pre(eta_ii,ii,r, theta) 
            !     call exit 
            ! endif
        end do

    end function


    real(kind)  function VSO_diag_elem(s1,s2,T0_mat,alpha_perp) result(matelem)
        type(an_ho_state), intent(in) :: s1, s2
        real(kind), intent(in) :: alpha_perp
        real(kind), dimension(nquad,0:max_n, 0:max_n) :: T0_mat
        integer :: ii
        real(kind) :: eta
        matelem = 0.0_kind
        if(s1%ml /= s2%ml .and. s1%ms /= s2%ms) then
            write(*,*) "Error: states are non diagonal in diagonal element calculation"
            call exit
        endif

        do ii = 1, nquad
            eta = lag_x(ii)
            matelem = matelem + S0(ii, s1, s2) * T0_mat(ii,s1%nz, s2%nz)* lag_w(ii) !== 2*lambda*sigma where sigma = +- 1/2
            
        end do
        matelem = matelem * alpha_perp**2 * s1%ml * s1%ms
    end function


    real(kind)  function VSO_off_diag_elem_v2(s1,s2,W0,alpha_z,alpha_perp) result(matelem)
        type(an_ho_state), intent(in) :: s1, s2
        real(kind), intent(in) :: alpha_z, alpha_perp
        real(kind), intent(in) :: W0(nquad,0:max_n, 0:max_n)
        matelem = 0

        if(s1%ml == s2%ml + 1 .and. s1%ms == s2%ms - 2) then
            matelem = I1(s1,s2,W0,alpha_z,alpha_perp) - I2(s1,s2,W0,alpha_z,alpha_perp)
        elseif( s1%ml == s2%ml - 1 .and. s1%ms == s2%ms + 2)then
            matelem = I1(s1,s2,W0,alpha_z,alpha_perp) + I2(s1,s2,W0,alpha_z,alpha_perp)
        else
            write(*,*) "Error: states dont have the correct quantum numbers in ms and ml"
            call exit
        endif
        matelem = matelem * (-0.5)
    end function


    subroutine print_levels(E, V)
        real(kind), intent(in) :: E(:), V(:,:)
        integer :: sz, i, numstates
        sz = size(E)
        numstates = 0

        ! write(*, '(A,A)')"Energy (MeV)"
        do i = sz, 1, -1
            !if(E(i) > 0) then!skip unbound states
            !    cycle
            !else
            numstates = numstates + 1
            write(*, '(1F11.5,A)', advance='no')E(i) ,","
            !endif
        end do
        write(*,*)
        write(*,'(A,I5)') "Number of bound particles:", numstates*2


    end subroutine


    subroutine get_levels(E_P, E_N,Z,A,def, max_N)
        real(kind),allocatable, intent(out) :: E_P(:), E_N(:)
        integer, intent(in) :: Z, A, max_N
        type(betadef), intent(in) :: def
        real(kind) :: Vn(num_n_states,num_n_states) ,Hp(num_p_states,num_p_states), Hn(num_n_states,num_n_states), Vp(num_p_states,num_p_states)
        type(an_ho_state) :: states_p(num_p_states), states_n(num_n_states)
        real(kind) :: hbaromega0, hbaromegaz, hbaromegaperp
        allocate(E_p(num_p_states),E_n(num_n_states))
        Vn = 0
        Vp = 0
        E_p = 0
        E_n =0
        hbaromega0 = 41.0_kind * A**(-1.0_kind/3.0_kind) !!MeV
        hbaromegaperp = def%omega_perp(hbaromega0) !! omega = Mev/hbar
        hbaromegaz = def%omega_z(hbaromega0)
        states_p = get_lowest_ho_states(max_n, num_p_states, hbaromegaz, hbaromegaperp)
        states_n = get_lowest_ho_states(max_n, num_n_states, hbaromegaz, hbaromegaperp)
        write(*,'(A)') 
        write(*,'(A)') "Calculating proton hamiltonian..."
        Hp = H_protons(states_p, Z, A, def, hbaromegaz, hbaromegaperp)
        write(*,'(A)') "Diagonalizing..."
        call diagonalize(E_p, Vp, Hp)
        write(*,'(A)') 
        write(*,'(A)') "Calculating neutron hamiltonian..."
        Hn = H_neutrons(states_n, Z, A, def, hbaromegaz, hbaromegaperp)
        write(*,'(A)') "Diagonalizing..."
        call diagonalize(E_n, Vn, Hn)
        write(*,*)
        write(*,'(A26,f5.2,a10)') "Time spent in V_ws:", time_ws, " seconds"
        write(*,'(A26,f5.2,a10)') "Time spent in V_so:", time_Vso, " seconds"
        write(*,'(A26,f5.2,a10)') "Time spent in V_C:", time_Vc, " seconds"
        write(*,'(A26,f5.2,a10)') "Time spent in diagonalise:", time_diag, " seconds"
        write(*,*)
        call write_result(Z, A, E_n, E_p, states_n, states_p, Vn, Vp)
    end subroutine

    function H_protons(states, Z,A,def,hbaromegaz, hbaromegaperp) result(H)
        implicit none
        integer :: Z, A
        type(betadef) :: def
        real(kind) :: hbaromegaz, hbaromegaperp, Vwsdepth, I,radius,radius_so
        type(an_ho_state) :: states(:)
        real(kind), dimension(size(states),size(states)) :: H, VSO, VC, VWS, Tkin
        radius = r0_p * A**(1.0_kind/3.0_kind)
        radius_so = r0_so_p* A**(1.0_kind/3.0_kind)
        I = (A-2.0_kind*Z)/A
        Vwsdepth = V0_ws * (1.0_kind+kappa_ws*I)
        write(*,'(A)') "Wood-Saxon term..."
        Vws = Vws_mat(states,def,radius, hbaromegaz,hbaromegaperp,mass_p,Vwsdepth)
        write(*,'(A)') "Spin-orbit term..."
        Vso = Vso_mat(states, def, radius_so, hbaromegaz,hbaromegaperp, mass_p, Vwsdepth, lambda_p)
        write(*,'(A)') "Kinetic-Energy term..."  
        Tkin = kin_en(states, hbaromegaz, hbaromegaperp)
        write(*,'(A)') "Coulomb term..."
        Vc = coul_mat(states, def, radius, Z, mass_p, hbaromegaz, hbaromegaperp)


        H = Vc + Vws + Vso + Tkin

    end function

    function H_neutrons(states, Z,A,def,hbaromegaz, hbaromegaperp) result(H)
        implicit none
        integer :: Z, A
        type(betadef) :: def
        real(kind) :: hbaromegaz, hbaromegaperp, Vwsdepth, I, radius,radius_so
        type(an_ho_state) :: states(:)
        real(kind), dimension(size(states),size(states)) :: H, VSO, VWS, Tkin
        radius = r0_n * A**(1.0_kind/3.0_kind)
        radius_so = r0_so_n* A**(1.0_kind/3.0_kind)
        I = (A-2.0_kind*Z)/A
        Vwsdepth = V0_ws * (1.0_kind-kappa_ws*I)
        write(*,'(A)') "Wood-Saxon term..."
        Vws = Vws_mat(states,def,radius, hbaromegaz,hbaromegaperp,mass_n,Vwsdepth)
        write(*,'(A)') "Spin-orbit term..."
        Vso = Vso_mat(states, def, radius_so, hbaromegaz,hbaromegaperp, mass_n, Vwsdepth, lambda_n)
        write(*,'(A)') "Kinetic-Energy term..."  
        Tkin = kin_en(states, hbaromegaz, hbaromegaperp)
        H = Vws + Vso + Tkin

    end function

    subroutine write_mat_to_file(mat)
        real(kind), dimension(:,:) :: mat
        integer :: ii
        character(len=100) :: str

        write(str, '(A,I5,A)') '(', size(mat,2), 'F10.3)'
        open(5, file='data/out/matout.dat')
        do ii = 1, size(mat,1)
            write(5, str) mat(ii,:)
        end do

    end subroutine


    logical function hasnan2d(mat)
        real(kind), dimension(:,:) :: mat

        integer :: ii, jj
        hasnan2d = .false.
        do ii = 1, size(mat,1)
            do jj = 1, size(mat,2)
                ! if (isnan(mat(ii,jj)))then
                !     hasnan2d = .true.
                !     return
                ! endif
            end do
        end do

    end function

    logical function hasnan3d(mat)
    real(kind), dimension(:,:,:) :: mat

    integer :: ii, jj
    hasnan3d = .false.
    do ii = 1, size(mat,1)
        if(hasnan2d(mat(ii,:,:))) then
            hasnan3d = .true.
            return
        endif
    end do

end function


subroutine print_shell_params(Z, A, def)
    integer :: Z, A
    type(betadef) :: def
    character(len=90) :: comm
    character :: symb
    character(len=87) :: line
    character(len=67) :: text
    integer :: ii
    symb = '+'
    do ii = 1, 90
        comm(ii:ii) = symb
    end do

    write(*,*) comm
    write(line,'(A)') "Shell model version 0.1"
    call print_line(symb, line)
    write(line,'(A)') "F. Agert"
    call print_line(symb, line)
    write(line,'(A)') "2025"
    call print_line(symb, line)
    write(line,'(A)') "Axially deformed harmonic oscillator basis"
    call print_line(symb, line)
    write(line,'(A)') "Wood-Saxon term, Coulomb term and spin-orbit term included"
    call print_line(symb, line)


    write(*,*) comm
    write(line,'(A)') "Input parameters for this calculation:"
    call print_line(symb, line)
!!!!!!!!!!!!
    text = ": proton number"
    write(line,'(A10,I10, A67)') "Z = ", Z, adjustl(text)
    call print_one_line(symb, line)
    text = ": neutron number"
    write(line,'(A10,I10, A67)') "N = ", A-Z, adjustl(text)
    call print_one_line(symb, line)
    text = ": mass number"
    write(line,'(A10,I10, A67)') "A = ", A, adjustl(text)
    call print_one_line(symb, line)
    line = ""
    call print_one_line(symb, line)
    text = ": beta_2 deformation value"
    write(line,'(A10,F10.3, A67)') "beta_2 = ", def%beta2, adjustl(text)
    call print_one_line(symb, line)

    text = ": beta_4 deformation value"
    write(line,'(A10,F10.3, A67)') "beta_4 = ", def%beta4, adjustl(text)
    call print_one_line(symb, line)

    line = ""
    call print_one_line(symb, line)

    text = ": Wood saxon potential depth [MeV]"
    write(line,'(A10,F10.1, A67)') "Vws = ", V0_ws, adjustl(text)
    call print_one_line(symb, line)
    text = ": Wood saxon isospin strength factor"
    write(line,'(A10,F10.2, A67)') "K = ", kappa_ws, adjustl(text)
    call print_one_line(symb, line)
    text = ": Wood saxon surface diffuseness [fm]"
    write(line,'(A10,F10.2, A67)') "a_ws = ", a_ws, adjustl(text)
    call print_one_line(symb, line)
    line = ""
    call print_one_line(symb, line)

    
    text = ": Wood saxon proton radius parameter [fm]"
    write(line,'(A10,F10.2, A67)') "r0_p = ", r0_p, adjustl(text)
    call print_one_line(symb, line)
    text = ": Wood saxon neutron radius parameter [fm]"
    write(line,'(A10,F10.2, A67)') "r0_n = ", r0_n, adjustl(text)
    call print_one_line(symb, line)
    text = ": Wood saxon spin orbit proton radius parameter [fm]"
    write(line,'(A10,F10.2, A67)') "r0_so_p = ", r0_so_p, adjustl(text)
    call print_one_line(symb, line)
    text = ": Wood saxon spin orbit neutron radius parameter [fm]"
    write(line,'(A10,F10.2, A67)') "r0_so_n = ", r0_so_n, adjustl(text)
    call print_one_line(symb, line)
    line = ""
    call print_one_line(symb, line)
    text = ": Spin orbit strength for protons"
    write(line,'(A10,F10.1, A67)') "lambda_p= ", lambda_p, adjustl(text)
    call print_one_line(symb, line)
    text = ": Spin orbit strength for neutrons"
    write(line,'(A10,F10.1, A67)') "lambda_n= ", lambda_n, adjustl(text)
    call print_one_line(symb, line)
    line = ""
    call print_one_line(symb, line)
    line = ""
    call print_one_line(symb, line)
    text = ": Maximum harmonic oscillator shell"
    write(line,'(A10,I10, A67)') "N_max = ", max_n, adjustl(text)
    call print_one_line(symb, line)
    text = ": Number of proton levels in diagonalisation"
    write(line,'(A10,I10, A67)') "N_p = ", num_p_states, adjustl(text)
    call print_one_line(symb, line)
    text = ": Number of proton levels in diagonalisation"
    write(line,'(A10,I10, A67)') "N_n = ", num_n_states, adjustl(text)
    call print_one_line(symb, line)
    line = ""
    call print_one_line(symb, line)
    line = ""
    call print_one_line(symb, line)
    text = ": Number of integration points used in the Gaussian quadrature"
    write(line,'(A10,I10, A67)') "N_GQ = ", nquad, adjustl(text)
    call print_one_line(symb, line)


    write(*,*) comm

    write(*,*)


end subroutine      
subroutine print_line(symb, line)
    character :: symb
    character(len=87) :: line
    write(*,'(1x,a,88x,a)') symb,symb
    write(*,'(1x,a,1x,A87,a)') symb,line,symb
    write(*,'(1x,a,88x,a)') symb,symb
end subroutine

subroutine print_one_line(symb, line)
    character :: symb
    character(len=87) :: line
    write(*,'(1x,a,1x,A87,a)') symb,line,symb
end subroutine

subroutine write_result(Z,A,E_n,E_p, states_n, states_p, V_n, V_p)
    real(kind) :: E_n(:), E_p(:)
    integer :: Z, A, ii, idx_p, idx_n, fermi_i_p, fermi_i_n
    type(betadef) :: def
    type(an_ho_state), dimension(:) :: states_n, states_p
    type(an_ho_state) :: state_p, state_n
    real(kind), dimension(:,:) :: V_n, V_p
    character(len=90) :: comm
    character :: symb
    character(len=87) :: line
    character(len=67) :: text
    character :: pp, pn
    character(len=8) :: fermi_p, fermi_n
    
    symb = '+'
    do ii = 1, 90
        comm(ii:ii) = symb
    end do

    write(*,*) comm

    write(line,'(A)') "Shell model calculation completed"
    call print_line(symb, line)


    line = " N   |m_j| p      E (MeV)         |m_j|  p    E (MeV)"
    call print_one_line(symb, line)
    if(mod(Z,2) == 0) then
        fermi_i_p = Z/2
    else
        fermi_i_p = Z/2 + 1
    endif

    if(mod(A-Z,2) == 0) then
        fermi_i_n = (A-Z)/2
    else
        fermi_i_n = (A-Z)/2 + 1
    endif
    do ii = 1, (A-Z)/2+2
        idx_p = maxloc(abs(V_p(ii,:)),1)
        state_p = states_p(idx_p)
        idx_n = maxloc(abs(V_n(ii,:)),1)
        state_n = states_n(idx_n)
        if(state_p%pi < 0) then
            pp = '-'
        else
            pp = '+'
        endif

        if(state_n%pi < 0) then
            pn = '-'
        else
            pn = '+'
        endif

        if(ii == fermi_i_p) then
            fermi_p = " <- e_f "
        else
            fermi_p = ""
        endif


        if(ii == fermi_i_n) then
            fermi_n = " <- e_f "
        else
            fermi_n = ""
        endif


        write(line,'(I3,1x, I3,A3,A1, 2x,F10.3, a8,1x, I3,A3,A1,2x, F10.3,a8)')ii, state_p%mj, "/2 ", pp, E_p(ii), fermi_p,state_n%mj, "/2 ", pn, E_n(ii), fermi_n
        call print_one_line(symb, line)
    end do


    write(*,*) comm
end subroutine

real(kind) function get_quad_Hmn(ii,n)
    !!Gets Hmn at quadrature integration point ii
    !!Precomputes and saves
    use def_ho, only:Hmn
    integer, intent(in) :: n !!order of modified polynomial
    integer, intent(in) :: ii !!index of quadrature root
    integer :: jj
    integer :: mm
    real(kind), save :: H(1:nquad, -1:max_n+1)
    logical, save :: first_time = .true.

    if(first_time) then
        call precompute_quad()
        H(:,-1) = 0
        do mm = 0,max_n+1
            do jj = 1,nquad
                H(jj,mm) = Hmn(her_x(jj), mm)
            end do
        end do
        first_time = .false.
    endif
    get_quad_Hmn = H(ii,n)
end function
real(kind) function get_quad_Hmnp(ii,n)
    !!Gets Hmnp at quadrature integration point ii
    !!Precomputes and saves
    integer, intent(in) :: n !!order of modified polynomial
    integer, intent(in) :: ii !!index of quadrature root
    integer :: jj
    integer :: mm
    real(kind), save :: Hp(1:nquad, -1:max_n)
    logical, save :: first_time = .true.

    if(first_time) then
        call precompute_quad()
        do mm = 0,max_n
            do jj = 1,nquad
                Hp(jj,mm) =get_quad_Hmn(jj,mm-1)*sqrt(0.5_kind*mm) - get_quad_Hmn(jj,mm+1)*sqrt(0.5_kind*(mm+1))
            end do
        end do
        first_time = .false.
    endif
    get_quad_Hmnp = Hp(ii,n)
end function

real(kind) function get_quad_gnl(ii,n,l)
    !!Gets Hmn at quadrature integration point ii
    !!Precomputes and saves
    use def_ho, only:gnl
    integer, intent(in) :: n,l !!order of modified polynomial
    integer, intent(in) :: ii !!index of quadrature root
    integer :: ll
    integer :: nn, jj
    real(kind), save :: g(1:nquad, 0:(max_n+1), 0:(max_n+1))
    logical, save :: first_time = .true.
    if(l < 0) print*, "ERR: l < 0"
    if(first_time) then
        call precompute_quad()
        do nn = 0,max_n+1
            do ll = 0, max_n+1
                do jj = 1,nquad
                    g(jj,nn,ll) = gnl(lag_x(jj),nn,ll)
                end do
            end do
        end do
        first_time = .false.
    endif
    if(n < 0) then
        get_quad_gnl = 0
    else
        get_quad_gnl = g(ii,n,l)
    endif
end function

real(kind) function get_quad_gnlp(ii,n,l)
    !!Gets Hmn at quadrature integration point ii
    !!Precomputes and saves
    integer, intent(in) :: n,l !!order of modified polynomial
    integer, intent(in) :: ii !!index of quadrature root
    integer :: ll
    integer :: nn, jj
    real(kind), save ::  gp(1:nquad, 0:(max_n), 0:(max_n))
    logical, save :: first_time = .true.

    if(first_time) then
        call precompute_quad()
        do nn = 0,max_n
            do ll = 0,max_n
                do jj = 1,nquad
                    gp(jj,nn,ll) = (get_quad_gnl(jj,nn,ll) * (2*nn+ll-lag_x(jj)) &
                    - get_quad_gnl(jj,nn-1,ll) * 2.0* sqrt(real(nn*(nn+ll),kind)))/sqrt(lag_x(jj))
        
                end do
            end do
        end do
        first_time = .false.
    endif
    if(n < 0) then
        get_quad_gnlp = 0
    else
        get_quad_gnlp = gp(ii,n,l)
    endif

end function
end module