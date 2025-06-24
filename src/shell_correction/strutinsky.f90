module strutinsky
    use constants, only: r_kind, max_n, pi
    use Hamiltonian, only: get_levels
    use def_ho, only: betadef, getnumstatesupto, lna
    use quadrule, only: laguerre_ss_compute, legendre_dr_compute
    use brent, only: zero_rc


    implicit none

    private

    public :: shell_correction, fermi_level_osc, fermi_level_sh, Peven, get_shell
    logical :: precomputed_quad = .false.
    integer, parameter :: gauss_order = 99
    real(r_kind), dimension(gauss_order) :: lag_x, lag_w, leg_x, leg_w
    contains

     subroutine precompute_quad()
        if(precomputed_quad) then
            return
        endif

        call laguerre_ss_compute(gauss_order, lag_x, lag_w)
        !!Gets weights and locations for where to evaluate 64 point integral. 
        !!Integral must be of form 
        !! 0 < x < infty  dx f(x) * exp(-x)
        !! where x = alpha^2 * rho^2

        
        call legendre_dr_compute(gauss_order, leg_x, leg_w)

        precomputed_quad = .true.
    end subroutine

    real(r_kind) function shell_correction(Z, A, def) result(Eshellcorr)
        integer, intent(in) :: Z, A
        real(r_kind), dimension(:), allocatable :: E_n, E_p
        type(betadef), intent(in) :: def
        real(r_kind) :: E_sh_corr_p, E_sh_corr_n, hbaromega0
        

        call get_levels(E_p, E_n, Z, A, def, max_n)
        hbaromega0 = 41.0_r_kind * A**(-1.0_r_kind/3.0_r_kind) !!MeV
        E_sh_corr_p = get_shell(Z, E_p, hbaromega0)
        E_sh_corr_n = get_shell(A-Z, E_n, hbaromega0)
        Eshellcorr = E_sh_corr_p + E_sh_corr_n
    end function

    real(r_kind) function get_shell(n_parts, levels, hbaromega0) result(E_shell)
        integer, intent(in)::n_parts
        real(r_kind), intent(in) :: hbaromega0
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind) :: E_sh, gamma, fermi_sh, fermi_smooth, E_smooth, E_smoothp, E_smoothm, df, d2f, prev
        real(r_kind), parameter :: dg =1e-5
        logical ::converged
        E_sh = shell_energy(n_parts, levels) !!gets full shell model shell energy
        gamma = hbaromega0
        fermi_sh = fermi_level_sh(n_parts, levels)
        converged = .false.
        prev = 0.0_r_kind
        do while(.not. converged)
            fermi_smooth = fermi_level_osc(n_parts, levels, gamma)
            
            E_smooth = smooth_e(levels, gamma, fermi_smooth)
            E_smoothp = smooth_e(levels, gamma+dg, fermi_smooth)
            E_smoothm = smooth_e(levels, gamma-dg, fermi_smooth)
            df = (E_smoothp - E_smoothm)/(2.0_r_kind*dg)
            d2f = (E_smoothp - 2.0_r_kind*E_smooth + E_smoothm)/(dg*dg)
            
            E_shell = E_sh - E_smooth
            print*, "Shell correction: ", E_shell   
            if(d2f .eq. 0.0_r_kind .or. abs(df/d2f) <1e-2 .or. abs(prev - E_shell) < 1e-3) then
                converged = .true.
    
            end if
            gamma = gamma - df/d2f
            print*, "gamma: ", gamma    
            prev = E_shell
        end do
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
                occ_numbers(ll) = occ_numbers(ll) + lag_w(ii) * exp(u) * f_m((-u+t),3)
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
            t = (fermi - E)/gamma
            do ii = 1, gauss_order
                u = lag_x(ii)
                F = F + lag_w(ii) * exp(u) * f_m((-u+t),3)*(-u + t)
            end do
        end do
        F = F * gamma**2 * 2
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

    pure real(r_kind) function avg_lev_dens(e, levels, gamma) result(g)
        real(r_kind), intent(in) :: e, gamma
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind) :: level_e
        integer :: i, numlevels
        integer, parameter :: degen = 2 !due to time reversal symmetry
        integer, parameter:: M = 3
        numlevels = size(levels)
        g = 0.0_r_kind
        do i = 1, numlevels
            level_e = levels(i)
            g = g +  f_m((level_e- e) / gamma, M)
        end do
        g = degen *g/gamma
    end function

    pure real(r_kind) function shell_energy(A, levels) result(energy)
        real(r_kind), dimension(:), intent(in) :: levels
        integer, intent(in) :: A !!Particle number
        integer :: ii, idx
        real(r_kind) :: e
        idx = 0
        energy = 0.0_r_kind
        do ii = 1, A
            idx = (ii + 1) / 2 !!1, 2 => 2 / 2 = 1, 3/2 => 1
            e = levels(idx)
            energy = energy + e
        end do
    end function

    pure real(r_kind) function fermi_level_sh(A, levels) result(fermi)
        real(r_kind), dimension(:), intent(in) :: levels
        integer, intent(in) :: A !!Particle number
        fermi = levels((A+1)/2) 
    end function

    real(r_kind) function smooth_e(levels, gamma, fermi) !!Energy of the smooth shell correction
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: gamma, fermi
        real(r_kind) :: u, occ_nbr(size(levels))
        integer :: ii

        call precompute_quad()
        occ_nbr = get_occ_numbers(levels, gamma, fermi)
        smooth_e = dot_product(occ_nbr, levels) + F(levels, gamma, fermi)
    end function

    real(r_kind) function fermi_level_osc(num_parts, levels, gamma) result(fermi_osc)
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: gamma
        integer, intent(in) :: num_parts
        real(r_kind) :: u, Ef_sh, lb, ub, arg, val
        integer :: ii, status
        real(r_kind), parameter :: tol = 1e-3


        Ef_sh = fermi_level_sh(num_parts, levels)
        lb = Ef_sh - gamma
        ub = Ef_sh + gamma

        status = 0
        call zero_rc(lb,ub,tol,arg,status,val) !Find occupation number by varying fermi level. Finds zero of sum(occ_nbr) - num_particles 
        ! Using the brent method to find the zero of the function
        do while(status /= 0)

            val = sum(get_occ_numbers(levels, gamma, arg),1) - num_parts
            call zero_rc(lb, ub, tol, arg, status, val)
            ! print*, "val: ", val, "fermi: ", arg
        end do
        fermi_osc = arg
        ! write(*,*) "Fermi level smooth: ", fermi_osc
        ! write(*,*) "Fermi level shell: ", Ef_sh
        ! write(*,*) "Occupation number:", sum(get_occ_numbers(num_parts, levels, gamma, fermi_osc),1)
        ! write(*,*) "Num particles:", num_parts
    end function
    


end module