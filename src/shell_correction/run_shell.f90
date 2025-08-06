program run_ho
    use constants
    use def_ho, only: an_ho_state, betadef, getnumstates, getnumstatesupto, kin_en, get_ho_states, get_ho_states_upto
    use hamiltonian
    implicit none


    integer :: A = 208, Z=82
    type(betadef) :: def 
    real(kind), dimension(:), allocatable :: E_n, E_p
    integer :: n, status
    character(len=50) :: filename, arg1,arg2
    real(kind) :: beta2, beta4
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



    write(filename,'(a,i0,a1,i0,a)')"data/out/levels_", Z,'_',A,'.dat'


    def = betadef(beta2 =beta2, beta4=beta4)
    allocate(E_p(num_p_states), E_n(num_n_states))
    call print_shell_params(Z,A,def)
    call get_levels(E_p, E_n,Z,A,def, max_n, .true.)


    open(4, file=filename)
    write(4,'(I3,A,I3,A)') Z, ",", A, ", = Z,A"
    do n = 1, num_p_states
        write(4, '(F10.5,A,F10.5)') E_p(n), "," , E_n(n)
    end do
    do n = num_p_states + 1, num_n_states
        write(4, '(F10.5,A,F10.5)') 0.0_kind, "," , E_n(n)
    end do
    close(4)
    
    
end program