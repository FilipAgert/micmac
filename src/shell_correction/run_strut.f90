program run_strut
    use strutinsky
    use constants
    implicit none

    real(r_kind),allocatable :: levels(:)
    real(r_kind) :: omega_0, gamma, ef, ef_sh, e_shell_corr
    integer :: A, Z

    gamma = 1
    A = 70
    omega_0 = 6!41.0 * A**(-1.0/3.0)
    levels = getHOlevels(14, omega_0)

    ef = fermi_level_osc(A, levels, gamma)
    ef_sh = fermi_level_sh(A, levels)
    e_shell_corr= get_shell(A, levels, omega_0)
    
    contains

    function getHOlevels(N, omega_0)
        integer :: N
        real(r_kind), allocatable :: getHOlevels(:)
        real(r_kind) :: omega_0

        integer :: ii, jj
        integer :: particles
        integer :: numlevels
        
        numlevels = (N+1)*(N+2)*(N+3)/6
        allocate(getHOlevels(numlevels))
        particles = 0
        do ii =0, N
            do jj = 1, (ii+1)*(ii+2)/2
                particles = particles + 1
                getHOlevels(particles) = omega_0 * (ii + 1.5_r_kind)
            end do
        end do
    end function


end program