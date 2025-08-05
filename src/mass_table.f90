module mass_table
    use constants
    implicit none
    public :: read_ame, read_elem

    character(len=100), parameter :: ame_file = "data/in/ame2020.dat"  ! Path to the AME data file
    real(kind), dimension(4000)  :: AME_BE, AME_BE_unc, AME_ME, AME_ME_unc
    integer, dimension(4000)  :: AME_Z, AME_A
    character(len=3), dimension(4000) :: AME_EL
    integer :: num_vals

    contains

    subroutine read_elem(BE, BE_unc, ME, ME_unc, El, Z, A)
        ! This subroutine reads the AME data file and extracts the binding energy and mass excess
        implicit none
        real(kind), intent(out) :: BE, BE_unc, ME, ME_unc
        character(len=3), intent(out) :: El
        integer, intent(in) :: Z, A
        integer :: idx
        BE = 0.0
        BE_unc = 0.0
        ME = 0.0
        ME_unc = 0.0
        do idx = 1, num_vals
            if (AME_Z(idx) == Z .and. AME_A(idx) == A) then
                BE = AME_BE(idx)
                BE_unc = AME_BE_unc(idx)
                ME = AME_ME(idx)
                ME_unc = AME_ME_unc(idx)
                El = AME_EL(idx)
                return
            end if
        end do

    end subroutine read_elem

    subroutine read_ame()
        ! This subroutine reads the AME data file and extracts the binding energy and mass excess
        implicit none
        real(kind) :: BE, BE_unc, ME, ME_unc
        character(len=3) :: El
        integer :: N, Z, A
        integer :: iunit, ios
        integer :: idx
        logical :: success
        character(len=200) :: line

        iunit = 10
        ! Open the file for reading
        open(unit=iunit, file=ame_file, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening file: ", ame_file
            return
        end if

        ! Read the file line by line until we find the desired Z and A
        idx = 1
        do while (.true.)
            read(iunit, '(A)', iostat=ios) line
            if (ios == 0) then
                if (line(1:1) .eq. '#') then
                    ! Skip comment lines
                    cycle
                end if
                call read_line(line, N,Z,A,El,ME,ME_unc,BE,BE_unc, success)
                if (.not. success) cycle
                AME_EL(idx) = El
                AME_Z(idx) = Z
                AME_A(idx) = A
                AME_BE(idx) = -BE * A / 1000 ! divide by 1000 to convert to MeV
                AME_BE_unc(idx) = BE_unc * A / 1000
                AME_ME(idx) = -ME * A/ 1000
                AME_ME_unc(idx) = ME_unc * A / 1000
                idx = idx + 1
            else
                exit
            end if
            
        end do

        ! Close the file
        close(iunit)
        num_vals = idx - 1

    end subroutine read_ame

    subroutine read_line(line, N, Z, A, str_El, ME, ME_unc, BE, BE_unc, success)
        character(len=200) :: line
        character(len=5)   :: str_N, str_Z, str_A 
        character(len=3)   :: str_El
        character(len=14)  :: str_ME
        character(len=12)  :: str_ME_unc
        character(len=13)  :: str_BE
        character(len=10)  :: str_BE_unc
        integer, intent(out) :: Z, A
        integer :: N, ios
        logical :: success
        real(kind), intent(out) :: ME, ME_unc, BE, BE_unc

        ! slice substrings according to column positions
        str_N          = line(5:9)
        str_Z          = line(10:14)
        str_A          = line(15:19)
        str_El         = line(21:23)
        str_ME         = line(29:42)
        str_ME_unc     = line(43:54)
        str_BE         = line(55:67)
        str_BE_unc     = line(69:79)

        ! Convert fields, checking for '*'
        read(str_N, *) N
        read(str_Z, *) Z
        read(str_A, *) A
        read(str_ME, *, iostat=ios) ME
        success = .true.
        if (ios /= 0) then
            ME = -999.0  ! or IEEE_VALUE(ME, IEEE_QUIET_NAN)
            success = .false.
        endif
        ! Repeat for each numeric field:
        read(str_ME_unc, *, iostat=ios) ME_unc
        if (ios /= 0) then 
            ME_unc = -999.0
            success = .false.
        endif

        read(str_BE, *, iostat=ios) BE
        if (ios /= 0) then 
            BE = -999.0
            success = .false.
        endif

        read(str_BE_unc, *, iostat=ios) BE_unc
        if (ios /= 0)then
            BE_unc = -999.0
            success = .false.
        endif

        

end subroutine read_line



end module mass_table