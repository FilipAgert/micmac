module strutinsky
    use constants, only: r_kind, max_n, pi
    use Hamiltonian, only: get_levels
    use def_ho, only: betadef, getnumstatesupto, lna
    use quadrule, only: laguerre_ss_compute, legendre_dr_compute
    use brent, only: zero_rc, zero_rc2
    use pairing


    implicit none

    private

    public :: shell_correction, fermi_level_osc, fermi_level_sh, Peven, get_shell, shell_energy, microscopic_corrections, level_dens
    logical :: precomputed_quad = .false.
    integer, parameter :: gauss_order = 64, M_ord = 6
    real(r_kind), dimension(gauss_order) :: lag_x, lag_w, leg_x, leg_w
    contains

     subroutine precompute_quad()
        if(precomputed_quad) then
            return
        endif

        call laguerre_ss_compute(gauss_order, lag_x, lag_w)
        !!Gets weights and locations for where to evaluate 64 point integral. 
        !!Integral must be of form F
        !! 0 < x < infty  dx f(x) * exp(-x)
        !! where x = alpha^2 * rho^2

        
        call legendre_dr_compute(gauss_order, leg_x, leg_w)

        precomputed_quad = .true.
    end subroutine

    real(r_kind) function microscopic_corrections(Z, A, def, fix_gamma) result(Emic_corr)
        integer, intent(in) :: Z, A
        real(r_kind), dimension(:), allocatable :: E_n, E_p
        type(betadef), intent(in) :: def
        logical, intent(in) :: fix_gamma
        real(r_kind) :: E_sh_corr_p, E_sh_corr_n, hbaromega0
        real(r_kind) :: gamma_p, efp, gamma_n, efn, epair_corr, dens_p, dens_n
        
        write(*,'(A)') "Calculating single particle energies..."
        call get_levels(E_p, E_n, Z, A, def, max_n)
        hbaromega0 = 41.0_r_kind * A**(-1.0_r_kind/3.0_r_kind) !!MeV
        call purge_levels(E_p, 3.2*hbaromega0, Z)
        call purge_levels(E_n, 3.2*hbaromega0, A-Z)
        write(*,*)
        write(*,'(A)') "Calculating shell correction for protons..."
        E_sh_corr_p = get_shell(gamma_p, efp, Z, E_p, hbaromega0, fix_gamma)

        write(*,*)
        write(*,'(A)') "Calculating shell correction for neutrons..."
        E_sh_corr_n = get_shell(gamma_n, efn,A-Z, E_n, hbaromega0, fix_gamma)
        dens_p = level_dens(efp, E_p, gamma_p) / 2
        dens_n = level_dens(efn, E_n, gamma_n) / 2
        write(*,'(A,F10.3,A)') "Shell correction energy: ", E_sh_corr_n + E_sh_corr_p, " MeV" 

        write(*,'(A)') "" 
        write(*,'(A)') "" 
        write(*,'(A)') "Calculating pairing correction with BCS method" 
        !write(*,'(A,2F10.3)') "Proton/neutron level density at fermi level:", dens_p * 2, dens_n * 2
        epair_corr = pairing_correction(E_p(1:Z), e_n(1:(A-Z)), Z, A, dens_p, dens_n)

        write(*,'(A,F10.3, A)') "Pairing correction: ", epair_corr, " MeV"
        Emic_corr = E_sh_corr_p + E_sh_corr_n + epair_corr
    end function

    real(r_kind) function shell_correction(Z, A, def, fix_gamma) result(Eshellcorr)
        integer, intent(in) :: Z, A
        real(r_kind), dimension(:), allocatable :: E_n, E_p
        type(betadef), intent(in) :: def
        logical, intent(in) :: fix_gamma
        real(r_kind) :: E_sh_corr_p, E_sh_corr_n, hbaromega0, gamma_p, efp, gamma_n, efn
        
        write(*,'(A)') "Calculating single particle energies..."
        call get_levels(E_p, E_n, Z, A, def, max_n)
        hbaromega0 = 41.0_r_kind * A**(-1.0_r_kind/3.0_r_kind) !!MeV
        call purge_levels(E_p, 3.2*hbaromega0, Z) !cutoff of levels g.r.t. 3.2 hbaromega above fermi level
        call purge_levels(E_n, 3.2*hbaromega0, A-Z)


        write(*,'(A)') "Calculating shell correction for protons..."
        E_sh_corr_p = get_shell(gamma_p, efp, Z, E_p, hbaromega0, fix_gamma)

        write(*,*)
        write(*,'(A)') "Calculating shell correction for neutrons..."
        E_sh_corr_n = get_shell(gamma_n, efn,A-Z, E_n, hbaromega0, fix_gamma)
        Eshellcorr = E_sh_corr_p + E_sh_corr_n
    end function

    subroutine purge_levels(levels, lim, N) !!purge levels lim above fermi level
        real(r_kind), dimension(:), allocatable, intent(inout) :: levels
        real(r_kind), dimension(size(levels)) :: temp
        real(r_kind), intent(in) :: lim
        integer, intent(in) :: N !!particle number
        real(r_kind) :: fermi
        integer :: ii
        fermi = fermi_level_sh(N, levels)
        
        ii = maxloc(levels, 1, levels .le. fermi+lim) !!find index of largest energy which is less than fermi level + limit

        temp(1:ii) = levels(1:ii)
        deallocate(levels)
        allocate(levels(ii))
        levels(1:ii) = temp(1:ii)
    end subroutine

    function Jac_strut(gamma, efermi, levels, n_parts)
        integer, intent(in)::n_parts
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: gamma, efermi
        real(r_kind) :: Jac_strut(2,2)
        integer :: ll
        real(r_kind) :: el, tl
        Jac_strut = 0
        do ll = 1, size(levels)
            el = levels(ll)
            tl = (efermi - el) / gamma
            Jac_strut(1,1) = Jac_strut(1,1) + f_m(tl, M_ord) / gamma !dN/def
            Jac_strut(1,2) = Jac_strut(1,2) - f_m(tl,M_ord) * tl / gamma !dn/dgamma

            Jac_strut(2,1) = Jac_strut(2,1) + f_m(tl, M_ord) * tl / gamma
            Jac_strut(2,2) = Jac_strut(2,2) - f_m(tl,M_ord) * tl**2 / gamma



        end do



    end function
    real(r_kind) function get_shell(gamma, efermi, n_parts, levels, hbaromega0, fix_gamma) result(E_shell) !gets shell correction energy
        integer, intent(in)::n_parts
        real(r_kind), intent(in) :: hbaromega0
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind) :: E_sh, fermi_sh,  J(2,2), func(2), x(2), invJ(2,2), px(2)
        logical, intent(in) :: fix_gamma
        logical ::converged
        integer :: iter
        real(r_kind), intent(out) :: gamma, efermi
        E_sh = shell_energy(n_parts, levels) !!gets full shell model shell energy
        gamma = hbaromega0*1.2
        fermi_sh = fermi_level_sh(n_parts, levels)
        converged = .false.
        if(fix_gamma) then !!fix gamma to 1.2 hbaromega0
            efermi = fermi_level_osc(n_parts, levels, gamma)
            E_shell = E_sh - smooth_e(levels, gamma, efermi)
            write(*,'(A,F10.3, A)') "Converged gamma: ", gamma, " MeV"
            write(*,'(A,F10.3, A)') "Fermi level: ", efermi, " MeV"
            write(*,'(A,F10.3, A)') "Converged shell correction energy: ", E_shell, " MeV"
            return
        endif
        x(1) = fermi_sh !!starting guess
        x(2) = gamma !!starting guess
        iter = 0
        do while(.not. converged)
            iter = iter + 1
            efermi = x(1)
            gamma = x(2)
            !write(*,'(A,2F10.3)') "Fermi level/ gamma:", efermi, gamma
            J = Jac_strut(gamma, efermi,levels, n_parts)
            !write(*,'(2F10.3)') J(1,:)
            !write(*,'(2F10.3)') J(2,:)
            invJ = J
            call inverse_SVD(invJ,2)
            func = [sum(get_occ_numbers(levels, gamma, efermi),1) - n_parts, F(levels, gamma, efermi)/gamma]
            px = x
            x = x - matmul(invJ, func) / 2
            if(dot_product(x-px,x-px) < 1e-6 ) then
                converged = .true.
                exit
            endif

            if(iter > 20) then
                gamma = 1.2 * hbaromega0
                efermi = fermi_level_osc(n_parts, levels, gamma)
                E_shell = E_sh - smooth_e(levels, gamma, efermi)
                write(*,'(A)') "Did not converge. Set gamma = 1.2 hbaromega"
                write(*,'(A,F10.3, A)') "Shell energy: ", E_shell, " MeV"
                return
            endif

        end do
        efermi = x(1)
        gamma = x(2)
        E_shell = E_sh - smooth_e(levels, gamma, efermi)
        Efermi = efermi
        write(*,'(A,I4, A)') "Gamma and fermi level converged after ", iter, " iterations"
        write(*,'(A,F10.3, A)') "Converged gamma: ", gamma, " MeV"
        write(*,'(A,F10.3, A)') "Fermi level: ", efermi, " MeV"
        write(*,'(A,F10.3, A)') "Converged shell correction energy: ", E_shell, " MeV"


    end function

    function get_occ_numbers(levels, gamma, fermi) result(occ_numbers)
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: gamma, fermi
        real(r_kind) :: u, E, t
        integer :: ii, numlevels, ll
        real(r_kind), dimension(size(levels)) :: occ_numbers

        call precompute_quad()
        numlevels = size(levels)
        occ_numbers = 0.0_r_kind
        do ll = 1, numlevels
            E = levels(ll)
            t = (fermi - E)/gamma

            
            do ii = 1, gauss_order
                u = lag_x(ii)
                occ_numbers(ll) = occ_numbers(ll) + lag_w(ii) * exp(u) * f_m((-u+t),M_ord)
            end do
        end do
        occ_numbers = occ_numbers * 2
    end function

    function F(levels, gamma, fermi) 
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: gamma, fermi
        real(r_kind) :: u, E, t
        integer :: ii, numlevels, ll
        real(r_kind) :: F

        call precompute_quad()
        numlevels = size(levels)
        F = 0.0_r_kind
        do ll = 1, numlevels
            E = levels(ll)
            !if(E > 0) cycle
            t = (fermi - E)/gamma
            do ii = 1, gauss_order
                u = lag_x(ii)
                F = F + lag_w(ii) * exp(u) * f_m((-u+t),M_ord)*(-u + t)
            end do
        end do
        F = F * gamma * 2
    end function

    pure elemental real(r_kind) function Peven(x,M)
        real(r_kind), intent(in) :: x
        integer, intent(in) :: M
        !!Even polynomial in the sturtinsky weighting integral
        Peven = lna(x*x,M,0.5_r_kind)
    end function

    pure elemental real(r_kind) function f_m(x,M)
        real(r_kind), intent(in) :: x
        integer, intent(in) :: M
        !!Even polynomial in the sturtinsky weighting integral
        f_m = Peven(x,M) * exp(-x*x)/sqrt(pi)
    end function




    real(r_kind) function smooth_e(levels, gamma, fermi) !!Energy of the smooth shell correction
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: gamma, fermi
        real(r_kind) :: occ_nbr(size(levels))

        call precompute_quad()
        occ_nbr = get_occ_numbers(levels, gamma, fermi)
        !write(*,'(A,200F10.3)') "Occupation numbers:", (occ_nbr)
        smooth_e = dot_product(occ_nbr, levels) + F(levels, gamma, fermi)
    end function

    real(r_kind) function fermi_level_osc(num_parts, levels, gamma) result(fermi_osc)
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: gamma
        integer, intent(in) :: num_parts
        real(r_kind) :: Ef_sh, lb, ub, arg, val, flb, fub
        integer :: status
        real(r_kind), parameter :: tol = 1e-3


        Ef_sh = fermi_level_sh(num_parts, levels)
        lb = Ef_sh - 1.2*gamma
        ub = Ef_sh + 2*gamma
        flb = 1
        lb = Ef_sh
        do while(flb > 0)
            lb = lb - gamma
            flb = sum(get_occ_numbers(levels, gamma, lb),1) - num_parts 
        end do
        fub = -1
        ub = Ef_sh
        do while(fub < 0)
            ub = ub + gamma
            fub = sum(get_occ_numbers(levels, gamma, ub),1) - num_parts 
        end do


        if( flb * fub > 0.0_r_kind) then
            write(*,'(A)') "Error: Same signs of lb and ub"
            write(*,'(A,2F10.3)') "F(lb), F(ub) : ", flb, fub
            call exit
        endif

        status = 0
        call zero_rc(lb,ub,tol,arg,status,val) !Find occupation number by varying fermi level. Finds zero of sum(occ_nbr) - num_particles 
        ! Using the brent method to find the zero of the function
        do while(status /= 0)

            val = sum(get_occ_numbers(levels, gamma, arg),1) - num_parts
            call zero_rc(lb, ub, tol, arg, status, val)
            ! print*, "val: ", val, "fermi: ", arg
        end do
        fermi_osc = arg
        !write(*,*) "Fermi level smooth: ", fermi_osc
        ! write(*,*) "Fermi level shell: ", Ef_sh
        ! write(*,*) "Occupation number:", sum(get_occ_numbers(levels, gamma, fermi_osc),1)
        ! write(*,*) "Num particles:", num_parts
    end function
    
    pure real(r_kind) function level_dens(E, levels, gamma)
        real(r_kind), intent(in) :: E, gamma, levels(:)
        integer :: ii
        real(r_kind) :: Eii
        level_dens = 0
        do ii = 1, size(levels)
            Eii = levels(ii)
            level_dens = level_dens + f_m((Eii - E)/gamma,M_ord)

        end do
        level_dens = level_dens / gamma
    end function

end module