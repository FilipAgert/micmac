program run_ho
    use constants
    use def_ho, only: an_ho_state, betadef, getnumstates, getnumstatesupto, kin_en, get_ho_states, get_ho_states_upto
    use hamiltonian
    implicit none


    integer :: A = 208, Z=82
    type(betadef) :: def 
    real(kind), dimension(:), allocatable :: E_n, E_p
    integer, dimension(:,:), allocatable :: qnp, qnn
    integer :: n, status, ierr
    character(len=50) :: folder_name
    character(len=50) :: filename, arg1,arg2
    real(kind) :: beta2, beta4
    logical :: exists
    beta2 = 0
    beta4 = 0

    ! Get the number of command-line arguments
    n = command_argument_count()
    if(n > 0) then
        if(n < 2) error stop "must input 0,2, 3 or 4 arguments: <Z> <A> <beta2> <beta4>" 
        call get_command_argument(1,arg1, status)
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

    print *, "Length of folder_name: ", len(folder_name)
    write(folder_name,'(a,i0,a1,i0)')"data/out/", Z,'_',A

    def = betadef(beta2 =beta2, beta4=beta4)
    allocate(E_p(num_p_states), E_n(num_n_states))
    call print_shell_params(Z,A,def)
    call get_levels(E_p,qnp, E_n,qnn,Z,A,def, N_max, .true.)

    write(*,*) Z, A
    inquire (file=folder_name,exist=exists)

    if (.not. exists) then
        call execute_command_line('mkdir -p ' // trim(folder_name), exitstat=ierr)
        if (ierr == 0) then
            print *, 'Directory created.'
        else
            print *, 'Failed to create directory.'
        end if
    end if





    write(filename,'(a,a)')trim(folder_name), "/levels_p.dat"
    open(4, file=filename)
    write(4,'(a3,1x,a3,1x,a10)') 'm_j', 'p', 'E (MeV)'
    do n = 1, num_p_states
        write(4, '(2(i3,1x),f10.5)') qnp(1,n), qnp(2,n), E_p(n)
    end do
    close(4)
    write(filename,'(a,a)')trim(folder_name), "/levels_n.dat"
    open(4, file=filename)
    write(4,'(a3,1x,a3,1x,a10)') 'm_j', 'p', 'E (MeV)'
    do n = 1, num_n_states
        write(4, '(2(i3,1x),f10.5)') qnn(1,n), qnn(2,n), E_n(n)
    end do
    close(4)
    
    
end program