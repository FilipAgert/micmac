program run_strut
    use strutinsky
    use constants
    use def_ho, only: betadef
    use Hamiltonian, only: print_shell_params
    implicit none

    real(kind) :: e_shell_corr
    integer :: A, Z
    type(betadef) :: def

    A = 208
    Z = 82
    def = betadef(beta2 = 0.0_kind, beta4=0.0_kind)
    call print_shell_params(Z, A, def)

    !e_shell_corr= get_shell(A-Z, levels, omega_0, .false.)
    e_shell_corr = microscopic_corrections(Z, A, def, .true.)
    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)
    write(*,'(A)') "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    write(*,'(a,1x,A30,F10.3,A,a44)')"+", "Microscopic correction energy:", e_shell_corr, " MeV",'+'
    write(*,'(A)') "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
    
    contains

    function getHOlevels(N, omega_0)
        integer :: N
        real(kind), allocatable :: getHOlevels(:)
        real(kind) :: omega_0

        integer :: ii, jj
        integer :: particles
        integer :: numlevels
        
        numlevels = (N+1)*(N+2)*(N+3)/6
        allocate(getHOlevels(numlevels))
        particles = 0
        do ii =0, N
            do jj = 1, (ii+1)*(ii+2)/2
                particles = particles + 1
                getHOlevels(particles) = omega_0 * (ii + 1.5_kind)
            end do
        end do
    end function


end program