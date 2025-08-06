program run_strut
    use strutinsky
    use constants
    use def_ho, only: betadef
    use Hamiltonian, only: print_shell_params
    implicit none

    real(kind) :: e_shell_corr
    integer :: A, Z,n
    type(betadef) :: def
    real(kind) :: beta2, beta4
    character(len=50) :: arg1, arg2
    beta2 = 0
    beta4 = 0

    A = 208
    Z = 82

    ! Get the number of command-line arguments
    n = command_argument_count()
    if(n > 0) then
        if(n < 2) error stop "must input 0,2, 3 or 4 arguments: <Z> <A> <beta2> <beta4>" 
        call get_command_argument(1,arg1)
        call get_command_argument(2,arg2)
        read(arg1,'(I20)')Z
        read(arg2,'(I20)')A
        if(n>2) then
            call get_command_argument(3,arg1)
            read(arg1,*)beta2
        endif
        if(n>3) then
            call get_command_argument(4,arg1)
            read(arg1,*)beta4
        endif
    endif
    def = betadef(beta2 = beta2, beta4=beta4)
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