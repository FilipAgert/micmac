module strutinsky
    use constants, only: r_kind, max_n, pi
    use Hamiltonian, only: get_levels
    use def_ho, only: betadef, getnumstatesupto, lna
    use quadrule, only: laguerre_ss_compute, legendre_dr_compute
    use brent, only: zero_rc, zero_rc2


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

    real(r_kind) function shell_correction(Z, A, def, fix_gamma) result(Eshellcorr)
        integer, intent(in) :: Z, A
        real(r_kind), dimension(:), allocatable :: E_n, E_p
        type(betadef), intent(in) :: def
        logical, intent(in) :: fix_gamma
        real(r_kind) :: E_sh_corr_p, E_sh_corr_n, hbaromega0
        
        write(*,'(A)') "Calculating single particle energies..."
        call get_levels(E_p, E_n, Z, A, def, max_n)
        hbaromega0 = 41.0_r_kind * A**(-1.0_r_kind/3.0_r_kind) !!MeV
        write(*,'(A)') "Calculating shell correction..."
        E_sh_corr_p = get_shell(Z, E_p, hbaromega0, fix_gamma)
        E_sh_corr_n = get_shell(A-Z, E_n, hbaromega0, fix_gamma)
        Eshellcorr = E_sh_corr_p + E_sh_corr_n
    end function

    real(r_kind) function get_shell(n_parts, levels, hbaromega0, fix_gamma) result(E_shell)
        integer, intent(in)::n_parts
        real(r_kind), intent(in) :: hbaromega0
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind) :: E_sh, gamma, fermi_sh, fermi_smooth,df, d2f, prev
        real(r_kind), parameter :: dg =1e-5
        logical, intent(in) :: fix_gamma
        logical ::converged
        E_sh = shell_energy(n_parts, levels) !!gets full shell model shell energy
        gamma = hbaromega0
        fermi_sh = fermi_level_sh(n_parts, levels)
        converged = .false.
        prev = 0.0_r_kind
        write(*,'(A,F10.3)') "start gamma: ", gamma
        write(*,'(A,F10.3)') "1.2 hbaromega0", 1.2*hbaromega0

        do while(.not. converged)
            fermi_smooth = fermi_level_osc(n_parts, levels, gamma)
            if(fix_gamma) then
                gamma = 1.2 * hbaromega0
                converged = .true.
                exit
            end if
            df = F(levels, gamma, fermi_smooth)/gamma
            d2f = (df - F(levels, gamma-dg, fermi_smooth)/(gamma-dg))/dg

            if(d2f .eq. 0.0_r_kind .or. abs(df/d2f) <1e-5 ) then
                converged = .true.

            end if
            gamma = gamma - df/d2f
            E_shell = E_sh - smooth_e(levels, gamma, fermi_smooth)
        end do
        E_shell = E_sh - smooth_e(levels, gamma, fermi_smooth)
        write(*,'(A,F10.3)') "Converged gamma: ", gamma
        write(*,'(A,F10.3)') "Converged shell energy: ", E_shell

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
        real(r_kind) :: occ_nbr(size(levels))

        call precompute_quad()
        occ_nbr = get_occ_numbers(levels, gamma, fermi)
        smooth_e = dot_product(occ_nbr, levels) + F(levels, gamma, fermi)
    end function

    real(r_kind) function fermi_level_osc(num_parts, levels, gamma) result(fermi_osc)
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: gamma
        integer, intent(in) :: num_parts
        real(r_kind) :: Ef_sh, lb, ub, arg, val
        integer :: status
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