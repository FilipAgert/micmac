module fiss_barr
    use constants
    use micmac, only: def_func, find_gs, J_mat, Eshell, alpha0sq
    use optimise, only: find_optimal_point
    use fitting, only: fit_param_cov
    !!Module for finding fission barrier and evaluating its uncertainty wrt parameter uncertainty.
    implicit none
    private
    public :: find_fiss_barr, cov_barr_gs



    contains



    subroutine find_fiss_barr(Bf, defbf, gs, defgs, params, Z, A)
        real(r_kind), intent(out) :: Bf, defbf, gs, defgs
        real(r_kind), intent(in) :: params(num_params)
        integer, intent(in) :: Z, A
        type(def_func) :: func !!Binding energy as a function of deformation
        real(r_kind), parameter :: search_lim = 15.0
        real(r_kind), parameter :: coarse_step = 0.02
        real(r_kind) :: x, pot, max
        logical :: conv

        call func%setup(params, Z, A)
        call find_gs(gs, defgs, params, Z, A)
        x = 0.0_r_kind
 
        x = defgs
        max = -huge(max)

        do while(x < search_lim)
            pot = func%eval(x)
            if(pot > max) then
                max = pot
                defbf = x
            endif
            x = x + coarse_step
        end do
        write(*,'(A,F10.3,A,F10.3)')"find optimal point. Def: ", defbf, ", with pot: ", max 
        call find_optimal_point(defbf, pot, conv, func, defbf -coarse_step*5, defbf+coarse_step*5, defbf - coarse_step/2.0_r_kind) !!Finer search for fission barrier.
        write(*,'(A,F10.3,A,F10.3)')"find optimal point. Def: ", defbf, ", with pot: ", pot 
        call find_gs(gs, defgs, params, Z, A)
        write(*,'(A,F10.3,A,F10.3)')"find gs point. Def: ", defgs, ", with pot: ", gs 
        Bf = pot - gs
    end subroutine

    !!Produce covariance matrix of fission barriers and ground state energies and the covariance between these for a number of nuclei
    subroutine cov_barr_gs(Bfs, BEs, defgss, defbfs, covBf, covBe, covAll, Zs, As, params, param_cov)
        implicit none
        integer, dimension(:), intent(in) :: Zs, As
        integer, parameter :: macparams = num_params - 3
        real(r_kind), intent(in) :: params(num_params),param_cov(num_params,num_params)
        real(r_kind), intent(out), dimension(size(Zs), size(Zs)) :: covBf, covBe !!Covariance of fission barrier and binding energy
        real(r_kind), intent(out), dimension(2*size(Zs), 2*size(Zs)) :: covAll  !!Covariance of fission barrier and be in the same matrix. 
        real(r_kind), dimension(size(Zs), size(Zs)) :: covBeSp
        !!Order is first entries are cov for binding energies, last entries are cov for barriers. in same order as nuclei in Z, A
        !!So for N nuclei, covAll(1, N+1) is the covariance between the barrier and g.s. for this nucleus.
        real(r_kind), intent(out), dimension(size(Zs)) :: Bfs, BEs, defgss, defbfs
        real(r_kind), dimension(size(Zs), num_params) :: Jmat_Sp, Jmat_Be
        real(r_kind), dimension(num_params, size(Zs)) :: Jmat_SpT, Jmat_BeT
        real(r_kind), dimension(size(Zs), size(Zs)) :: covSp
        integer :: ii , Z, A, N, jj
        real(r_kind) :: Be, defgs, Bf, defbf, alpha0
        N = size(Zs)

        do ii = 1, N
            Z = Zs(ii)
            A = As(ii)
            call find_fiss_barr(Bf, defbf, Be, defgs, params, Z, A)
            Bfs(ii) = Bf
            defbfs(ii) = defbf
            BEs(ii) = Be
            defgss(ii) = defgs      
            alpha0 = sqrt(alpha0sq(params(7), params(4),A))
            write(*,*) "frac gs:", defgs/alpha0
            write(*,*) "frac bf:", defbf/alpha0
        end do

        !When deformations are found, create J_matrices
        Jmat_Be = J_mat(Zs, As, defgss, N, params)
        Jmat_Sp = J_mat(Zs, As, defbfs, N, params)
        write(*,*) "JvecBe:"
        write(*,'(8F10.3)') Jmat_Be
        write(*,*) "JvecSp:"
        write(*,'(8F10.3)') Jmat_Sp
        Jmat_SpT = transpose(Jmat_Sp)
        Jmat_BeT = transpose(Jmat_Be)

        covBe = matmul(Jmat_Be, matmul(param_cov, Jmat_BeT))
        covSp = matmul(Jmat_Sp, matmul(param_cov, Jmat_SpT))

        covAll(1:N,1:N) = covBe
        
        covBeSp = matmul(Jmat_Be, matmul(param_cov, Jmat_SpT))

        covBf = covSp + covBe - covBeSp - transpose(covBeSp) !!Covariance of fission barriers.
        write(*,*) "covSp:"
        write(*,*) covSp
        write(*,*) "covBe:"
        write(*,*) covBe
        write(*,*) "covBeSp:"
        write(*,*) covBeSp
        covAll(N+1:2*N, N+1:2*N) = covBf
        covAll(1:N, N+1:2*N) = covBeSp - covBe
        covAll(N+1:2*N, 1:N) = transpose(covAll(1:N, N+1:2*N))  
    end subroutine


    





end module