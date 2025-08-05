program table_writer
    use constants
    use micmac, only:find_gs, mass_excess, Eshell, deformation
    
    implicit none
    integer, parameter :: startZ = 92, startA = 220
    integer, parameter :: maxN = 258 !!Largest known magic number
    integer, parameter :: maxZ = 126
    integer, parameter :: numA = 64
    integer :: num_nuclei
    integer, dimension(2000) :: Zs, As
    integer :: Z, A, idx, Alow, Ahigh
    real(kind) :: params(num_params)
    logical :: found

    idx = 0
    Alow = startA
    Ahigh = startA+numA
    do Z = startZ, maxZ
        do A = Alow, Ahigh
            idx = idx + 1
            Zs(idx) = Z
            As(idx) = A
        end do
        Alow = Alow + 2
        Ahigh = Ahigh + 1
    end do
    
    num_nuclei = idx
    params = fitted_params

    call write_mass_table(params, Zs, As, .true.)
    contains
    
    subroutine write_mass_table(params, Zs, As, write_to_file)
        ! Write the table to a file
        real(kind), intent(in) :: params(num_params)
        logical, intent(in) :: write_to_file
        integer, intent(in), dimension(:) :: Zs, As 
        integer :: iunit, i
        real(kind) :: BE, ME, Esh
        type(deformation) :: def
        iunit = 20  ! Output unit number
        if(write_to_file) then 
            open(unit=iunit, file='data/out/table.dat', status='replace')
        else
            iunit = 6
        endif

        write(iunit, '(A)') "//N   Z   A     Mass defect (MeV)   Esh (MeV)   alpha2   alpha3   alpha4"
        do i = 1, num_nuclei
            ! Calculate binding energy and mass excess
            call find_gs(BE, def, found, params, Zs(i), As(i))
            if(.not. found) then
                write(*,*) "did not find G.S. for Z,A: ", Zs(i), As(i)
                cycle
            endif
            ME = mass_excess(BE, Zs(i), As(i))
            Esh = Eshell(As(i), Zs(i), params(5), params(6), params(4),params(7), def)
            ! Write the data for each nucleus
            write(iunit, '(I3, 1x, I3, 1x,I3,3x, F10.3,5x, F10.3,6x, F6.3, 4x, F6.3,4x, F6.3)') &
                As(i)-Zs(i), Zs(i), As(i), ME, Esh, def%alphas(2),def%alphas(3),def%alphas(4)
            end do

        if(write_to_file) then
            close(iunit)
        endif
    end subroutine write_mass_table

end program