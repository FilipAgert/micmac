module fiss_barr
    use constants
    use micmac, only: def_func, find_gs
    use optimise, only: one_dim, find_optimal_point
    !!Module for finding fission barrier and evaluating its uncertainty wrt parameter uncertainty.
    implicit none
    private
    public :: find_fiss_barr



    contains



    subroutine find_fiss_barr(Bf, def, params, Z, A)
        real(r_kind), intent(out) :: Bf, def
        real(r_kind), intent(in) :: params(num_params)
        integer, intent(in) :: Z, A
        type(def_func) :: func !!Binding energy as a function of deformation
        real(r_kind), parameter :: search_lim = 5.0
        real(r_kind), parameter :: coarse_step = 0.05
        real(r_kind) :: curr, x, pot, gs, defgs, max
        logical :: conv

        call func%setup(params, Z, A)
        x = 0.0_r_kind
        print*, "test"
        max = func%eval(x)
        print*, "test"
        x = x + coarse_step
        max = -huge(max)

        do while(x < search_lim)
            pot = func%eval(x)
            if(pot > max) then
                max = pot
                def = x
            endif
        end do

        call find_optimal_point(def, pot, conv, func, def -coarse_step, def+coarse_step, def - coarse_step/2.0_r_kind) !!Finer search for fission barrier.
        call find_gs(gs, defgs, params, Z, A)
        Bf = curr - gs    
    end subroutine





end module