program run_strut
    use strutinsky
    use constants
    use def_ho, only: betadef
    implicit none

    real(r_kind):: levels(680)
    real(r_kind) :: omega_0, gamma, e_shell_corr
    integer :: A, Z

    gamma = 1
    A = 208
    Z = 82
    omega_0 = 41.0 * A**(-1.0/3.0)
    write(*,'(A,2I4)') "Running shell correction calculation for Z,A=",Z,A
    !e_shell_corr= get_shell(A-Z, levels, omega_0, .false.)
    e_shell_corr = shell_correction(Z, A, betadef(beta2=0.0_r_kind, beta4=0.0_r_kind), .false.)
    write(*,'(A,F10.3,A)') "Shell correction energy:", e_shell_corr, " MeV"
    
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