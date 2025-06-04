module fiss_barr
    use constants
    use micmac, only: be_def_func, alpha_to_beta, find_gs
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
        procedure(one_dim), pointer :: def_func !!Binding energy as a function of deformation
        real(r_kind), parameter :: search_lim = 5.0
        real(r_kind), parameter :: coarse_step = 0.05
        integer :: ii
        real(r_kind) :: prev, curr, x, currslope, prevslope, pot, gs, defgs
        logical :: conv

        !!Step from zero deformation. Save each point where slope changes (0 derivative.) The last point where slope changes
        !! and slope is negative is the fission barrier

        ! def_func => be_def_func(params, Z, A)
        ! if (.not. associated(def_func)) then
        !     print *, "ERROR: def_func is not associated!"
        !     stop
        ! endif
        
        x = 0.0_r_kind
        print*, "test"
        prev = def_func(x)
        print*, "test"
        prevslope = 0.0_r_kind
        x = x + coarse_step
        do while(x < search_lim)
            curr = def_func(x)
            currslope = curr - prev

            if(currslope * prevslope < 0.0_r_kind)  then !!Change of slope. Derivative was 0
                def = x
                pot = curr
            end if
            prev = curr
            prevslope = currslope
            x = x + coarse_step

        end do


        if(currslope .ge. 0.0_r_kind) then
            !!If slope positive, we have not reached a barrier
            Bf = 0.0_r_kind
            print*, "Error: did not find fission barrier"
        else
            call find_optimal_point(def, pot, conv, def_func, def -coarse_step, def+coarse_step, def - coarse_step/2.0_r_kind) !!Finer search for fission barrier.
            call find_gs(gs, defgs, params, Z, A)
            Bf = curr - gs
        endif
    end subroutine





end module