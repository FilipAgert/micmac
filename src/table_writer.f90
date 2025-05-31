module table_writer
    use constants
    use micmac
    implicit none


    private
    public :: write_mass_table


    contains
    
    subroutine write_mass_table(params, Zs, As, write_to_file)
        ! Write the table to a file
        real(r_kind), intent(in) :: params(num_params)
        logical, intent(in) :: write_to_file
        integer, intent(in), dimension(:) :: Zs, As 
        integer :: iunit, i
        real(r_kind) :: BE, ME, Esh
        iunit = 20  ! Output unit number
        if(write_to_file) then 
            open(unit=iunit, file='data/out/table.dat', status='replace')
        else
            iunit = 6
        endif

        write(iunit, '(A)') "//Z   A     Mass defect (MeV)   Esh (MeV)   beta2    beta4    beta6"
        do i = 1, size(Zs)
            ! Calculate binding energy and mass excess

            BE = binding_energy(params, Zs(i), As(i))
            ME = mass_excess(BE, Zs(i), As(i))
            Esh = shell_corr(Zs(i), As(i), params)

            ! Write the data for each nucleus
            write(iunit, '(I3, 3x,I3,3x, F10.3,5x, F10.3,6x, F5.2, 4x, F5.2,4x, F5.2)') &
                Zs(i), As(i), ME, Esh, 0.0,0.0,0.0
            end do

        if(write_to_file) then
            close(iunit)
        endif
    end subroutine write_mass_table

end module table_writer