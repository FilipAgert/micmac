program run_pairing
    use constants
    use pairing
    use iso_fortran_env, only:iostat_end
    implicit none
    
    real(r_kind), allocatable :: levels_p(:), levels_n(:), temp_p(:), temp_n(:)
    integer :: ii, Z, A, iostat, leastneg, numvals
    character(len=500) :: trash, line
    real(r_kind) :: ef, delta, G, I, E, smoothe, E0, epair
    allocate(temp_n(5000), temp_p(5000))
    open(3, file='data/out/levels.dat', status='old', action='read')
    read(3, '(I3,A1,I3,A)', iostat = iostat) Z, trash, A, trash
    ii = 0
    do 
        read(3, '(A)', iostat = iostat) line
        Select case(iostat)
        case(0)
            ii = ii + 1

            read(line,'(F10.5,A2,F10.5)')temp_p(ii), trash, temp_n(ii)
        case(iostat_end)
            exit
        case default
            write(*,*) "error in reading file"
            exit
        end select
    end do
    leastneg = maxloc(temp_p, 1, temp_p < 0)
    numvals = Z
    allocate(levels_p(numvals), levels_n(A-Z))

    levels_n = temp_n(1:A-Z)

    levels_p = temp_p(1:numvals)
    deallocate(temp_n, temp_p)
    I = (A-2.0_r_kind*Z) / A
    G = (g0p + g1p * I) / A
    write(*,'(A,2I4)') "Pairing for nucleus:", Z, A
    write(*,'(A,F10.3,A,F10.3)') "G:", G, ", 2/G:", 2.0_r_kind / G
    !call solve_BCS_eq(ef, delta, levels_p, Z, G)
    delta = 5.72/Z**(1.0_r_kind/3.0_r_kind) * exp(0.119*I -7.89*I**2)
    E = pairing_corr_x(levels_p, Z, 2.0_r_kind, G, delta) !pairin(levels_p, ef, delta, G)
    write(*,'(A,F10.3, A)') "Pairing correction: ", E, " MeV"
    delta = 5.72/(A-Z)**(1.0_r_kind/3.0_r_kind) * exp(-0.119*I -7.89*I**2)
    G = (g0n + g1n * I) / A
    E = pairing_corr_x( levels_n, A-Z, 4.0_r_kind, G, delta) !pairin(levels_p, ef, delta, G)
    write(*,'(A,F10.3, A)') "Pairing correction: ", E, " MeV"



end program