module Hamiltonian
    use constants
    use def_ho
    use quadrule, only:legendre_dr_compute, hermite_ek_compute, laguerre_ss_compute
    use optimise, only:func_1d, conj_grad_method
    implicit none


    private
    public :: diagonalize, WS_pot, VC_pot, dist_min, Vso_mat, commutator, VWS_mat, coul_mat, print_levels, H_protons, H_neutrons, get_levels
    public :: S0, Sp, Sm, T0, Tplus, gnlp, gnl, Tminus, VSO_off_diag_elem_v2, VSO_off_diag_elem, el_pot
    integer, parameter :: gauss_order =50
    real(r_kind), dimension(gauss_order) :: her_x, her_w, lag_x, lag_w, leg_x, leg_w
    logical :: precomputed_quad = .false.



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



    subroutine precompute_quad()
        if(precomputed_quad) then
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

    pure real(r_kind) function coul_charge_density(ZZ, radius) result(rho)
        integer, intent(in) :: ZZ !!proton number
        real(r_kind), intent(in) :: radius
        rho = (ZZ-1) *3.0_r_kind/(4.0_r_kind * pi*4*pi * epsilonzero * radius**3)  !! Charge /( Volume * 4piepsilonzero)

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
        ! write(*,*) N
        lda = N
        lwork = 20*N
        allocate(work(lwork))
        call dsyev(jobz, uplo, N, H, lda, w, work, lwork, info)

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
        deallocate(work)

    end subroutine

    pure real(r_kind) function el_pot(r,rad,charge_dens) !!electric potential of a uniformly charged sphere
        real(r_kind),intent(in) :: r, rad, charge_dens
        if(r .le. rad) then !!inside nucleus
            el_pot = 4.0_r_kind*pi/2.0_r_kind * (rad**2 - r**2/3.0_r_kind)!4(3.0-(r/rad)**2)/(2.0*rad)
        else!outside nucleus
            el_pot = 4.0_r_kind*pi/3.0_r_kind * rad**3/r
        endif
        el_pot = el_pot * charge_dens
    end function

    real(r_kind) function eval_vc(self,r, theta) !!axially symmetric. function of z and rho
        class(VC_pot), intent(in) :: self
        real(r_kind), intent(in) :: r, theta
        !Do the integral int_V     dr^3/|r-r'| at the point r.
        integer :: iu, it, ix
        real(r_kind) :: u, t, x, int_theta, int_radius_ub, rolling, sintheta, sinacosu, cospipit, costheta

        if(self%def%eq(spherical_def) )then !!if spherical, use analytical formula
            eval_vc = el_pot(r, self%radius, self%charge_dens)
            return
        endif


        call precompute_quad()

        sintheta = sin(theta)
        costheta = cos(theta)

        rolling = 0.0_r_kind
        call omp_set_num_threads(num_threads)
        !//$omp parallel do private (iu, it, ix, u, int_theta, int_radius_ub, t, x) &
        !//$omp& reduction(+:rolling)

        do iu = 1, gauss_order
            u = leg_x(iu)
            int_theta = acos(u)
            int_radius_ub = surfRadius(int_theta, self%def, self%radius)
            sinacosu = sin(acos(u))
            do it = 1, gauss_order
                t = leg_x(it)
                cospipit = cos(pi + pi*t)

                do ix = 1,gauss_order
                    x = leg_x(ix)

                    rolling = rolling + leg_w(iu)*leg_w(it)*leg_w(ix)*pi*int_radius_ub**3*(x+1.0_r_kind)**2 / (8.0_r_kind) &
                            / sqrt(r**2 + int_radius_ub**2 * (x+1.0_r_kind)**2 / 4.0_r_kind &
                                    + r*int_radius_ub*(x+1.0_r_kind)* (sintheta*sinacosu*cospipit + u*costheta))
                    
                    ! if(isnan(rolling)) then
                    !     print*, "isnan in coul: iu, it, ix:"
                    !     print*, iu, it, ix
                    ! endif
                end do
            end do
        end do
        !//$omp end parallel do
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
        surfRadius = R0 * (1.0_r_kind + def%beta2 * Y20(theta) + def%beta4*Y40(theta))
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
        if(r < surfRadius(theta, self%def,self%radius)) then
            dist = - dist
        endif
        eval_ws = -self%V_0/ (1 + exp(dist/aws))
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

    function coul_mat(states, def, R0,ZZ, mass, omegaz, omegaperp)
        implicit none
        type(an_ho_state), intent(in) :: states(:)
        type(betadef), intent(in) :: def !!deformation of body
        real(r_kind), intent(in) :: R0 !!radius as r0*A**(1/3). Where r0 is r0_so from constants.f90
        real(r_kind), dimension(size(states),size(states)) :: coul_mat
        integer, intent(in) :: ZZ!!proton number
        real(r_kind) :: alpha_z, alpha_perp
        type(VC_pot) :: cp
        real(r_kind), intent(in) :: omegaz, omegaperp !!Frequencies in units MeV/hbar
        real(r_kind), intent(in) :: mass!!Frequencies in units MeV/c^2
        real(r_kind), dimension(0:max_n,0:max_n, gauss_order) :: Z_mat
        integer ::numstates, row, col
        type(an_ho_state) :: s1, s2
        call precompute_quad()

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

    end function

    subroutine comp_zmat(Zmat, pot,alpha_z, alpha_perp) !!use the reccurence relations to compute the Zmat completely from the main diagonals
        real(r_kind), dimension(0:max_n,0:max_n, gauss_order), intent(out) :: Zmat
        logical, dimension(0:max_n, 0:max_n) :: computed
        class(potential) :: pot
        real(r_kind), intent(in) :: alpha_z, alpha_perp
        real(r_kind) :: eta
        integer :: nz, nzp, ii, diag
        character(len=50) :: str

        call precompute_quad()
        Zmat = 0
        computed = .false.

        do nz = 0, max_n !!compute main diagonals
            do nzp = nz, min(max_n, nz+1)
                do ii = 1, gauss_order
                    eta = lag_x(ii)
                    Zmat(nzp, nz, ii) = Z_matelem(eta, ii,nzp,nz,pot,alpha_z,alpha_perp)

                end do
                computed(nzp, nz) = .true.
            end do
        end do

        do nzp = 1,max_n !reflect off diagonal.
            nz = nzp -1
            Zmat(nzp, nz, :) = Zmat(nz, nzp, :)
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
                        do ii = 1, gauss_order
                            eta = lag_x(ii)
                            Zmat(nzp + 1, nz - 1, ii) = Z_matelem(eta, ii,nzp + 1,nz - 1,pot,alpha_z,alpha_perp)
                            
                        end do
                        computed(nzp + 1, nz - 1) = .true.
                    endif


                    Zmat(nzp, nz, :) = Zmat(nzp, nz, :) + Zmat(nzp+1, nz-1,:) * sqrt((nzp+1.0_r_kind)/nz)
                    if(nz==0) then
                        write(*,*) "ERR: divide by zero nz, nzp: ", nz, nzp
                    endif
                endif

                if(within_bounds(nzp - 1,0, max_n) .and. within_bounds(nz - 1, 0, max_n)) then
                    if(.not. computed(nzp - 1, nz - 1)) then
                        do ii = 1, gauss_order
                            eta = lag_x(ii)
                            Zmat(nzp - 1, nz - 1, ii) = Z_matelem(eta, ii,nzp - 1,nz - 1,pot,alpha_z,alpha_perp)
        
                        end do
                        computed(nzp - 1, nz - 1) = .true.
                    endif


                    Zmat(nzp, nz, :) = Zmat(nzp, nz, :) + Zmat(nzp-1, nz-1,:) * sqrt(real(nzp,r_kind)/nz)
                    if(nz==0) then
                        write(*,*) "ERR: divide by zero nz, nzp: ", nz, nzp
                    endif
                endif

                if(within_bounds(nzp,0, max_n) .and. within_bounds(nz - 2, 0, max_n)) then
                    if(.not. computed(nzp, nz - 2)) then
                        do ii = 1, gauss_order
                            eta = lag_x(ii)
                            Zmat(nzp, nz - 2, ii) = Z_matelem(eta, ii,nzp,nz - 2,pot,alpha_z,alpha_perp)
        
                        end do
                        computed(nzp, nz - 2) = .true.
                    endif


                    Zmat(nzp, nz, :) = Zmat(nzp, nz, :) - Zmat(nzp, nz-2,:) * sqrt((nz-1.0_r_kind)/nz)
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
                Zmat(nzp, nz,:) = Zmat(nz, nzp,:)
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
        real(r_kind), intent(in) :: R0 !!radius as r0*A**(1/3). Where r0 is r0_so from constants.f90
        real(r_kind), intent(in) :: WS_depth !!WS potental well depth
        real(r_kind), dimension(size(states),size(states)) :: VWS_mat
        real(r_kind), intent(in) :: omegaz, omegaperp !!Frequencies in units MeV/hbar
        real(r_kind), intent(in) :: mass!!Frequencies in units MeV/c^2
        real(r_kind) :: alpha_z, alpha_perp
        real(r_kind), dimension(0:max_n,0:max_n, gauss_order) :: Z_mat
        type(WS_pot) :: VWS
        integer ::numstates, row, col
        logical :: hasnan
        type(an_ho_state) :: s1, s2
        call precompute_quad()

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
        call write_mat_to_file(VWS_mat)
        call nanloc(hasnan, row, col, VWS_mat)
        if(hasnan) then
            write(*,*) "Error: VWS Matrix has a NaN at row, col:", row, col
        endif
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
        real(r_kind), dimension(0:max_n, 0:max_n, gauss_order) :: T0, W0
        type(WS_pot) :: so_ws
        logical :: hasnan
        integer ::numstates, row, col, l1, l2, ms1, ms2
        type(an_ho_state) :: s1, s2
        call precompute_quad()

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
        !call write_mat_to_file(Vso_mat)
        call nanloc(hasnan, row, col, Vso_mat)
        if(hasnan) then
            write(*,*) "Error: VSO Matrix has a NaN at row, col:", row, col
        endif
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

     real(r_kind)  function Sp(eta,s1,s2)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: eta
        real(r_kind) :: K
        
        K = s2%ml + s2%ms*0.5_r_kind
        print*, "K+0.5: ", K+0.5, " l1, l2:", s1%ml , s2%ml
        print*, "first term: ",- gnl(eta,s1%nr,s1%ml)*gnl(eta,s2%nr,s2%ml) *(-K+0.5_r_kind) /sqrt(eta)
        print*, "second term:", - gnlp(eta,s1%nr, s1%ml) * gnl(eta,s2%nr,s2%ml)
        Sp = - gnl(eta,s1%nr,s1%ml)*gnl(eta,s2%nr,s2%ml) *(K+0.5_r_kind)/sqrt(eta) - gnlp(eta,s1%nr, s1%ml) * gnl(eta,s2%nr,s2%ml)
        
    end function

     real(r_kind)  function Sm(eta,s1,s2)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: eta
        real(r_kind) :: K
        K = s2%ml + s2%ms*0.5_r_kind
        print*, "K-0.5 ", K-0.5, " l1, l2:", s1%ml , s2%ml
        print*, "first term: ",-gnl(eta,s1%nr, s1%ml)*gnl(eta,s2%nr,s2%ml)*(K-0.5_r_kind)/sqrt(eta)
        print*, "second term:", +gnl(eta,s1%nr,s1%ml)*gnlp(eta,s2%nr,s2%ml)
        Sm = Sp(eta,s2,s1)!-gnl(eta,s1%nr, s1%ml)*gnl(eta,s2%nr,s2%ml)*(K-0.5_r_kind)/sqrt(eta)+gnl(eta,s1%nr,s1%ml)*gnlp(eta,s2%nr,s2%ml)

    end function

    function T0_mat(SO_WS, alpha_z, alpha_perp)
        real(r_kind), dimension(0:max_n, 0:max_n, gauss_order) :: T0_mat
        type(WS_pot), intent(in) :: SO_WS
        real(r_kind), intent(in) :: alpha_z, alpha_perp

        call comp_zmat(T0_mat, SO_WS, alpha_z, alpha_perp)
    end function

    function W0_mat(SO_WS, alpha_z, alpha_perp) !!int tilde(h) h S
        !!TOdo: use recurence relation to not have to compute entire matrix
        real(r_kind), dimension(0:max_n, 0:max_n, gauss_order) :: W0_mat
        type(WS_pot), intent(in) :: SO_WS
        real(r_kind), intent(in) :: alpha_z, alpha_perp
        real(r_kind) :: eta
        integer :: nz, nzp, ii
        W0_mat = 0
        call precompute_quad()
        do nz = 0, max_n
            do nzp = 0, max_n
                do ii = 1, gauss_order
                    eta = lag_x(ii)
                    W0_mat(nzp, nz, ii) = W_matelem(eta, ii, nzp, nz, SO_WS, alpha_z, alpha_perp)
                end do
            end do
        end do


    end function

    real(r_kind) function T0(eta,eta_ii,nz1,nz2, SO_WS,alpha_z, alpha_perp)
        integer, intent(in) :: nz1, nz2
        real(r_kind), intent(in) :: eta, alpha_z, alpha_perp
        type(WS_pot), intent(in) :: SO_WS
        integer :: i
        integer, intent(in) :: eta_ii
        real(r_kind) :: zeta, z, rho,r ,theta
        T0 = 0
        call precompute_quad()
        rho = eta_to_rho(eta, alpha_perp)
        do i = 1,gauss_order
            zeta = her_x(i)
            z = zeta_to_z(zeta, alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho, z)


            T0 = T0 + her_w(i)*SO_WS%eval_pre(eta_ii,i,r,theta) * Hmn(zeta, nz1) * Hmn(zeta, nz2)
        end do
    end function

    real(r_kind) function Tminus(eta,eta_ii,s1,s2, SO_WS,alpha_z, alpha_perp)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: eta, alpha_z, alpha_perp
        type(WS_pot), intent(in) :: SO_WS
        integer, intent(in) :: eta_ii
        integer :: i
        real(r_kind) :: zeta, z, rho,r ,theta
        Tminus = 0
        call precompute_quad()
        rho = eta_to_rho(eta, alpha_perp)
        do i = 1,gauss_order
            zeta = her_x(i)
            z = zeta_to_z(zeta, alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho, z)


            Tminus = Tminus + her_w(i)*SO_WS%eval_pre(eta_ii,i,r,theta) * Hmn(zeta, s1%nz) * Hmnp(zeta, s2%nz)
        end do
    end function

    real(r_kind) function Tplus(eta,eta_ii,s1,s2, SO_WS,alpha_z, alpha_perp)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: eta, alpha_z, alpha_perp
        type(WS_pot), intent(in) :: SO_WS
        integer, intent(in) :: eta_ii
        integer :: i
        real(r_kind) :: zeta, z, rho,r ,theta
        Tplus = 0
        call precompute_quad()
        rho = eta_to_rho(eta, alpha_perp)
        do i = 1,gauss_order
            zeta = her_x(i)
            z = zeta_to_z(zeta, alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho, z)


            Tplus = Tplus + her_w(i)*SO_WS%eval_pre(eta_ii,i,r,theta) * Hmn(zeta, s2%nz) * Hmnp(zeta, s1%nz)
        end do
    end function

    real(r_kind) function I1(s1,s2, W0,alpha_z, alpha_perp)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) ::  alpha_z, alpha_perp
        real(r_kind), dimension(0:max_n, 0:max_n, gauss_order), intent(in) :: W0
        integer :: nn, zz
        real(r_kind) :: zeta, z, rho,r ,theta, eta, weta, wzeta, term


        call precompute_quad()
        I1 = 0
        do nn = 1, gauss_order
            eta = lag_x(nn)
            weta = lag_w(nn)
            term = 0
            rho = eta_to_rho(eta, alpha_perp)

            term = W0(s1%nz, s2%nz, nn) * s2%ml + W0(s2%nz, s1%nz, nn) * s1%ml
            term = term * weta * gnl(eta,s1%nr, s1%ml) * gnl(eta, s2%nr,s2%ml) / sqrt(eta)
            I1 = I1 + term
        end do
        I1 = I1 * alpha_perp*alpha_z
        ! if(isnan(I1)) then
        !     write(*,*) "ERR: I1 is not a number"
        !     write(*,*) "W0:", W0(s1%nz, s2%nz, 1), W0(s2%nz, s1%nz, 1)
        !     write(*,*) "gnl: ",gnl(eta,s1%nr, s1%ml), gnl(eta, s2%nr,s2%ml)
        ! endif
    end function

    real(r_kind) function I2(s1,s2, W0,alpha_z, alpha_perp)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) ::  alpha_z, alpha_perp
        integer :: nn, zz
        real(r_kind), dimension(0:max_n, 0:max_n, gauss_order), intent(in) :: W0
        real(r_kind) :: zeta, z, rho,r ,theta, eta, weta, wzeta, term1, term2, pot, etafac


        call precompute_quad()
        I2 = 0
        do nn = 1, gauss_order
            eta = lag_x(nn)
            weta = lag_w(nn)
            term1 = 0
            term2 = 0
            rho = eta_to_rho(eta, alpha_perp)

            term1 = W0(s1%nz, s2%nz, nn)
            term2 = W0(s2%nz, s1%nz, nn)

            etafac = weta / sqrt(eta)
            term1 = term1 * etafac * gnl(eta, s1%nr, s1%ml) * gnlp(eta, s2%nr, s2%ml)
            term2 = term2 * etafac * gnlp(eta, s1%nr, s1%ml) * gnl(eta, s2%nr, s2%ml)

            I2 = I2 + term1 - term2
        end do
        I2 = I2 * alpha_perp*alpha_z

    end function


    pure elemental real(r_kind) function zeta_to_z(zeta, alpha)
        real(r_kind), intent(in) :: zeta, alpha
        zeta_to_z = zeta / alpha
    end function

    pure elemental real(r_kind) function eta_to_rho(eta, alpha)
        real(r_kind), intent(in) :: eta, alpha
        eta_to_rho = sqrt(eta) / alpha
    end function

    

    real(r_kind) function pot_elem_zmatpre(s1,s2,Zmat) result(matelem)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in), dimension(0:max_N, 0:max_N, gauss_order) :: Zmat
        integer :: ii
        real(r_kind) :: eta
        call precompute_quad()
        matelem = 0.0_r_kind
        if(s1%ml /= s2%ml .or. s1%ms /= s2%ms) then
            return
        endif

        do ii = 1, gauss_order
            eta = lag_x(ii)
            matelem = matelem + lag_w(ii) * Zmat(s1%nz, s2%nz, ii) * gnl(eta, s1%nr,s1%ml) * gnl(eta, s2%nr, s2%ml) 
        end do

    end function

    real(r_kind) function Z_matelem(eta, eta_ii,nz1,nz2,pot,alpha_z, alpha_perp)
        integer, intent(in) :: nz1, nz2
        real(r_kind), intent(in) :: alpha_z, alpha_perp, eta
        class(potential), intent(in) :: pot
        integer :: ii
        integer, intent(in) ::eta_ii
        real(r_kind) :: zeta, z, rho, r, theta
        Z_matelem = 0.0_r_kind
        rho = eta_to_rho(eta, alpha_perp)

        do ii = 1, gauss_order
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

            Z_matelem = Z_matelem + her_w(ii) * Hmn(zeta, nz1) * Hmn(zeta, nz2)*pot%eval_pre(eta_ii,ii,r, theta) 
            ! if(printflag) then
            !     write(*,'(F10.3)'),(pot%eval(0.0_r_kind,0.0_r_kind))
            ! endif

        end do
    end function

    real(r_kind) function W_matelem(eta, eta_ii,nz1,nz2,pot,alpha_z, alpha_perp) !! W(eta) = int zeta    h(zeta) h'(zeta) S(eta, zeta)
        integer, intent(in) :: nz1, nz2
        real(r_kind), intent(in) :: alpha_z, alpha_perp, eta
        class(potential), intent(in) :: pot
        integer :: ii
        integer, intent(in) ::eta_ii
        real(r_kind) :: zeta, z, rho, r, theta
        W_matelem = 0.0_r_kind
        rho = eta_to_rho(eta, alpha_perp)
        do ii = 1, gauss_order
            zeta = her_x(ii)
            z = zeta_to_z(zeta,alpha_z)
            theta = theta_cyl(rho, z)
            r = rad_cyl(rho ,z)
            W_matelem = W_matelem + her_w(ii) * Hmnp(zeta, nz1) * Hmn(zeta, nz2)*pot%eval_pre(eta_ii,ii,r, theta) 
            ! if(isnan(W_matelem)) then
            !     print*, "Err Wmatelem is Nan"
            !     print*, "hmnp:", Hmnp(zeta, nz1)
            !     print*,  "Hmn:", Hmn(zeta, nz2)
            !     print*, "pot:", pot%eval_pre(eta_ii,ii,r, theta) 
            !     call exit 
            ! endif
        end do

    end function


    real(r_kind)  function VSO_diag_elem(s1,s2,T0_mat,alpha_perp) result(matelem)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: alpha_perp
        real(r_kind), dimension(0:max_n, 0:max_n, gauss_order) :: T0_mat
        integer :: ii
        real(r_kind) :: eta
        matelem = 0.0_r_kind
        if(s1%ml /= s2%ml .and. s1%ms /= s2%ms) then
            write(*,*) "Error: states are non diagonal in diagonal element calculation"
            call exit
        endif

        do ii = 1, gauss_order
            eta = lag_x(ii)
            matelem = matelem + S0(eta, s1, s2) * T0_mat(s1%nz, s2%nz, ii)* lag_w(ii) !== 2*lambda*sigma where sigma = +- 1/2
            
        end do
        matelem = matelem * alpha_perp**2 * s1%ml * s1%ms
    end function

    real(r_kind)  function VSO_off_diag_elem(s1,s2,SO_WS,alpha_z,alpha_perp) result(matelem)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: alpha_z, alpha_perp
        type(WS_pot), intent(in) :: SO_WS
        integer :: ii
        real(r_kind) :: eta
        matelem = 0

        if(s1%ml == s2%ml + 1 .and. s1%ms == s2%ms - 2) then

        elseif( s1%ml == s2%ml - 1 .and. s1%ms == s2%ms + 2)then

        else
            write(*,*) "Error: states dont have the correct quantum numbers in ms and ml"
            call exit
        endif

        do ii = 1, gauss_order
            eta = lag_x(ii)
            matelem = matelem + lag_w(ii) * (Sm(eta, s1, s2) * Tplus(eta, ii,s1, s2, SO_WS, alpha_z, alpha_perp) + &
                                            Sp(eta, s1, s2) * Tminus(eta, ii,s1, s2, SO_WS, alpha_z, alpha_perp))
        end do
        matelem = matelem * alpha_perp*alpha_z
    end function

    real(r_kind)  function VSO_off_diag_elem_v2(s1,s2,W0,alpha_z,alpha_perp) result(matelem)
        type(an_ho_state), intent(in) :: s1, s2
        real(r_kind), intent(in) :: alpha_z, alpha_perp
        integer :: ii
        real(r_kind) :: eta
        real(r_kind), intent(in) :: W0(0:max_n, 0:max_n, gauss_order)
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
        real(r_kind), intent(in) :: E(:), V(:,:)
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
        real(r_kind),allocatable, intent(out) :: E_P(:), E_N(:)
        integer, intent(in) :: Z, A, max_N
        type(betadef), intent(in) :: def
        integer :: numstates, idx, shelldegen, n,ii
        real(r_kind), allocatable :: Vn(:,:), H(:,:), Vws(:,:), Tkin(:,:), Vso(:,:), Vc(:,:), gs(:), Hp(:,:), Hn(:,:), Vp(:,:)
        type(an_ho_state), allocatable :: states(:)
        real(r_kind) :: hbaromega0, hbaromegaz, hbaromegaperp, Egs
        numstates = getnumstatesupto(max_n)
        allocate(E_p(numstates),E_n(numstates), Vn(numstates,numstates), Vp(numstates,numstates),states(numstates), Hn(numstates,numstates), gs(numstates), Hp(numstates,numstates))
        Vn = 0
        Vp = 0
        E_p = 0
        E_n =0
        hbaromega0 = 41.0_r_kind * A**(-1.0_r_kind/3.0_r_kind) !!MeV
        hbaromegaperp = def%omega_perp(hbaromega0) !! omega = Mev/hbar
        hbaromegaz = def%omega_z(hbaromega0)
        numstates = getnumstatesupto(max_n)
        idx = 1
        do n = 0, max_n
            shelldegen = getnumstates(n)
            states(idx:idx+shelldegen - 1) = get_ho_states(n)
            idx = idx + shelldegen
        end do

        
        print*, "Calculating proton hamiltonian..."
        Hp = H_protons(states, Z, A, def, hbaromegaz, hbaromegaperp)
        call diagonalize(E_p, Vp, Hp)

        Hn = H_neutrons(states, Z, A, def, hbaromegaz, hbaromegaperp)
        call diagonalize(E_n, Vn, Hn)
            
    end subroutine

    function H_protons(states, Z,A,def,hbaromegaz, hbaromegaperp) result(H)
        implicit none
        integer :: Z, A
        type(betadef) :: def
        real(r_kind) :: hbaromegaz, hbaromegaperp, Vwsdepth, I,radius,radius_so
        type(an_ho_state) :: states(:)
        real(r_kind), dimension(size(states),size(states)) :: H, VSO, VC, VWS, Tkin
        radius = r0_p * A**(1.0_r_kind/3.0_r_kind)
        radius_so = r0_so_p* A**(1.0_r_kind/3.0_r_kind)
        I = (A-2.0_r_kind*Z)/A
        Vwsdepth = V0_ws * (1.0_r_kind+kappa_ws*I)
        write(*,'(A,F10.3,A)') "Wood-Saxon hamiltonian. V0 = ", Vwsdepth, " MeV"
        Vws = Vws_mat(states,def,radius, hbaromegaz,hbaromegaperp,mass_p,Vwsdepth)
        write(*,*) "Spin-orbit hamiltonian..."
        Vso = Vso_mat(states, def, radius_so, hbaromegaz,hbaromegaperp, mass_p, Vwsdepth, lambda_p)
        write(*,*) "Kinetic-Energy..."  
        Tkin = kin_en(states, hbaromegaz, hbaromegaperp)
        write(*,*) "Coulomb-potential..."
        Vc = coul_mat(states, def, radius, Z, mass_p, hbaromegaz, hbaromegaperp)


        H = Vc + Vws + Vso + Tkin

    end function

    function H_neutrons(states, Z,A,def,hbaromegaz, hbaromegaperp) result(H)
        implicit none
        integer :: Z, A
        type(betadef) :: def
        real(r_kind) :: hbaromegaz, hbaromegaperp, Vwsdepth, I, radius,radius_so
        type(an_ho_state) :: states(:)
        real(r_kind), dimension(size(states),size(states)) :: H, VSO, VWS, Tkin
        radius = r0_n * A**(1.0_r_kind/3.0_r_kind)
        radius_so = r0_so_n* A**(1.0_r_kind/3.0_r_kind)
        I = (A-2.0_r_kind*Z)/A
        Vwsdepth = V0_ws * (1.0_r_kind-kappa_ws*I)
        write(*,'(A,F10.3,A)') "Wood-Saxon hamiltonian. V0 = ", Vwsdepth, " MeV"
        Vws = Vws_mat(states,def,radius, hbaromegaz,hbaromegaperp,mass_n,Vwsdepth)
        write(*,*) "Spin-orbit hamiltonian..."
        Vso = Vso_mat(states, def, radius_so, hbaromegaz,hbaromegaperp, mass_n, Vwsdepth, lambda_n)
        write(*,*) "Kinetic-Energy hamiltonian..."  
        Tkin = kin_en(states, hbaromegaz, hbaromegaperp)
        H = Vws + Vso + Tkin

    end function

    subroutine write_mat_to_file(mat)
        real(r_kind), dimension(:,:) :: mat
        integer :: ii
        character(len=100) :: str

        write(str, '(A,I5,A)') '(', size(mat,2), 'F10.3)'
        open(5, file='data/out/matout.dat')
        do ii = 1, size(mat,1)
            write(5, str) mat(ii,:)
        end do

    end subroutine

    subroutine nanloc(nan, row,col,mat)

        real(r_kind), dimension(:,:) :: mat
        integer :: row, col
        logical :: nan
        integer :: ii, jj
        nan = .false.
        do ii = 1, size(mat,1)
            do jj = 1, size(mat,2)
                ! if (isnan(mat(ii,jj)))then
                !     row = ii
                !     col = jj
                !     nan = .true.
                !     return
                ! endif
            end do
        end do

    end subroutine

    logical function hasnan2d(mat)
        real(r_kind), dimension(:,:) :: mat

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
    real(r_kind), dimension(:,:,:) :: mat

    integer :: ii, jj
    hasnan3d = .false.
    do ii = 1, size(mat,1)
        if(hasnan2d(mat(ii,:,:))) then
            hasnan3d = .true.
            return
        endif
    end do

end function
end module