module pairing
    use constants
    use strutinsky, only: fermi_level_sh
    implicit none
    private


    contains

    pure real(r_kind) function part_num(levels, fermi, delta) !!particle number in BCS equations
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: fermi, delta
        real(r_kind) :: ev
        integer :: vv
        part_num = 0
        do vv = 1,size(levels)
            ev = levels(vv)
            part_num = part_num + 1 - (ev-fermi)/sqrt((ev-fermi)**2 + delta**2)
        end do
    end function

    pure real(r_kind) function pairing_gap(levels, fermi, delta)
        real(r_kind), dimension(:), intent(in) :: levels
        real(r_kind), intent(in) :: fermi, delta
        real(r_kind) :: ev
        integer :: vv
        pairing_gap = 0
        do vv = 1,size(levels)
            ev = levels(vv)
            pairing_gap = pairing_gap + 1.0_r_kind/sqrt((ev-fermi)**2 + delta**2)
        end do
    end function

    subroutine solve_BCS_eq(efermi, delta, levels, num_part, G) !!https://math.stackexchange.com/questions/268991/how-to-solve-simultaneous-equations-using-newton-raphsons-method
        real(r_kind), dimension(:), intent(in) :: levels
        integer, intent(in) :: num_part
        real(r_kind), intent(in) :: G !!pairing strength
        real(r_kind), intent(out) :: efermi !!fermi level
        real(r_kind), intent(out) :: delta 
        real(r_kind), dimension(2,2) :: J !!jacobian
        real(r_kind), dimension(2) :: x, f
        !!use multivariate newtons method
        delta = 1 !!initial guess
        efermi = fermi_level_sh(num_part, levels) !!initial guess
        x = [delta, efermi]
        





    end subroutine

    subroutine inv_J(J) !!take inverse of jacobian
        real(r_kind), intent(inout) :: J(2,2)

        
    end subroutine


end module pairing