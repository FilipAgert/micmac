program run_pairing
    use constants
    use pairing
    use iso_fortran_env, only:iostat_end
    implicit none
    
    real(r_kind), allocatable :: levels_p(:), levels_n(:), temp_p(:), temp_n(:)
    integer :: ii, Z, A, iostat
    character(len=500) :: trash, line
    real(r_kind) :: ef, delta, G, I, E
    allocate(temp_n(5000), temp_p(5000))
    open(3, file='data/out/levels.dat', status='old', action='read')
    read(3, '(I3,A1,I3,A)', iostat = iostat) Z, trash, A, trash
    ii = 0
    do 
        read(3, '(A)', iostat = iostat) line
        Select case(iostat)
        case(0)
            ii = ii + 1

            read(line,'(F10.5,A,F10.5)')temp_p(ii), trash, temp_n(ii)
        case(iostat_end)
            exit
        case default
            write(*,*) "error in reading file"
            exit
        end select
    end do

    allocate(levels_p(Z), levels_n(A-Z))
    levels_n = temp_n(1:Z)
    levels_p = temp_p(1:(A-Z))
    deallocate(temp_n, temp_p)
    I = (A-2.0_r_kind*Z) / A
    G = (g0p + g1p * I) / A
    write(*,'(A,2I4)') "Pairing for nucleus:", Z, A
    write(*,'(A,F10.3,A,F10.3)') "G:", G, ", 2/G:", 2.0_r_kind / G
    call solve_BCS_eq(ef, delta, levels_p, Z, G)
    E = Ebcs(levels_p, ef, delta, G)
    write(*,'(A,F10.3)') "BCS energy: ", E

end program