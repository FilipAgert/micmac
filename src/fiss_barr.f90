module fiss_barr
    use constants
    use micmac, only: def_func_ho, find_gs, J_mat, Eshell, alpha0sqho, deformation
    use optimise, only: find_optimal_point
    use fitting, only: fit_param_cov
    use stack_module
    !!Module for finding fission barrier and evaluating its uncertainty wrt parameter uncertainty.
    implicit none
    private
    public :: find_fiss_barr
    integer, parameter :: max_num_neighbours = num_def_params**3 - 1

    contains

    subroutine find_fiss_barr(Bf, defbf, gs, defgs, params, Z, A)
        real(r_kind), intent(out) :: Bf, gs
        real(r_kind), intent(in) :: params(num_params)
        integer, intent(in) :: Z, A
        type(def_func_ho) :: func !!Binding energy as a function of deformation
        type(deformation), intent(out) :: defgs, defbf
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

    subroutine fill_tub(Bf, bfidx, potsurf, gsidx, fissidx)
        real(r_kind), intent(out) :: Bf
        real(r_kind), dimension(:,:,:) :: potsurf
        integer, dimension(num_def_params), intent(in) :: gsidx, fissidx
        integer, dimension(num_def_params), intent(out) :: bfidx
        integer, dimension(num_def_params) :: c
        integer, dimension(max_num_neighbours, num_def_params) :: neigh
        integer :: num, dimsize(num_def_params), n(num_def_params), ii
        logical, dimension(size(potsurf,1), size(potsurf,2), size(potsurf,3)) :: filled, visited
        real(r_kind), parameter :: delta_e = 0.1, maxe = 50
        real(r_kind) :: Egs, E, sp
        type(Stack) :: idx_stack, next_iter
        dimsize = shape(potsurf)
        idx_stack%top_index = 0
        next_iter%top_index = 0
        filled = .false.

        Egs = potsurf(gsidx(1),gsidx(2),gsidx(3))
        E = Egs

        filled(gsidx(1),gsidx(2),gsidx(3)) = .true.
        visited(gsidx(1),gsidx(2),gsidx(3)) = .true.
        call idx_stack%push(gsidx)

        eloop: do while(E < maxE)
            visited = .false.
            do while(.not. idx_stack%is_empty())
                call idx_stack%pop(c)

                if(coord_eq(c,fissidx)) then !!Reached fission state can exit loop. Found energy of barrier.
                    sp = E
                    exit eloop
                end if

                call get_neighbours(neigh, num, c, dimsize)
                call quicksort_coords(neigh, num, fissidx)

                do ii = num, 1, -1 !!Put neighborus in order highest dist to lowest dist onto stack.
                    n = neigh(ii,:)
                    

                    if (.not. filled(n(1),n(2),n(3))) then !!if not filled, check if its visited or not.
                        if(.not. visited(n(1),n(2),n(3))) then
                            visited(n(1),n(2),n(3)) = .true.
                            if(potsurf(n(1),n(2),n(3)) > E) then
                                call next_iter%push(n)
                            else !!Here, for each energy iteration, also push energies into a stack or queue. When we exit loop, we can use this stack to find position of fission barrier.
                                filled(n(1),n(2),n(3)) = .true.
                                call idx_stack%push(n)
                            endif
                        endif
                    endif
                end do
            end do

            idx_stack = next_iter
            E = E + delta_e
        end do eloop
        Bf = sp - Egs
    end subroutine

    pure logical function coord_eq(c1,c2)
        integer, dimension(num_def_params), intent(in) :: c1,c2
        integer :: ii
        coord_eq = .true.
        do ii = 1,num_def_params
            coord_eq = coord_eq .and. c1(ii) .eq. c2(ii)
        end do  
    end function

    !!Gets all neighbouring points given a coordinate.
    subroutine get_neighbours(neighbours, num, coord, dimsize)
        integer, dimension(max_num_neighbours,3), intent(out) :: neighbours
        integer, intent(out) :: num
        integer, intent(in), dimension(num_def_params) :: coord, dimsize
        integer :: i,j,k, idx
        idx = 0

        do i = MAX(1, coord(1)-1), MIN(dimsize(1), coord(1)+1)
            do j = MAX(1, coord(2)-1), MIN(dimsize(1), coord(2)+1)
                do k = MAX(1, coord(3)-1), MIN(dimsize(1), coord(3)+1)
                    if(coord(1) == i .and. coord(2) == j .and. coord(3) == k) then
                    else
                        idx = idx + 1
                        neighbours(idx,:) = [i,j,k]
                    endif

                end do
            end do
        end do
        num = idx
    end subroutine

    !!Sorts coordinates in list according to shortest distance to end_coord (DFS)
    subroutine quicksort_coords(coords, N_elems, end_coord)
        implicit none
        integer, intent(inout), dimension(max_num_neighbours, num_def_params) :: coords !!Coordinates to sort.
        integer,intent(in) :: N_elems !!Number of elements in list to sort.
        integer, intent(in), dimension(num_def_params) :: end_coord !!Coordinate to sort against.
        real(kind=r_kind), dimension(max_num_neighbours) :: dist
        integer :: i

        do i = 1,N_elems !Finds distance to end to sort against.
            dist(i) = norm2(real(coords(:,i),r_kind)-real(end_coord,r_kind))
        end do
        call quicksort(coords, dist, 1, N_elems) !Sorts coords according to distances in dist
    end subroutine quicksort_coords

    recursive subroutine quicksort(coord, a, first, last) !Sorts coord and a according to a.
        implicit none
        real(kind=r_kind)  a(max_num_neighbours), x, t
        integer, dimension(max_num_neighbours,num_def_params), intent(inout) :: coord
        integer, dimension(num_def_params) :: temp_coord
        integer first, last
        integer i, j

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
            do while (a(i) < x)
                i=i+1
            end do
            do while (x < a(j))
                j=j-1
            end do
            if (i >= j) exit
            temp_coord = coord(i,:)
            coord(i,:) = coord(j,:)
            coord(j,:) = temp_coord
            t = a(i);  a(i) = a(j);  a(j) = t
            i=i+1
            j=j-1
        end do
        if (first < i-1) call quicksort(coord, a, first, i-1)
        if (j+1 < last)  call quicksort(coord, a, j+1, last)
    end subroutine quicksort

    !!Produce covariance matrix of fission barriers and ground state energies and the covariance between these for a number of nuclei
    ! subroutine cov_barr_gs(Bfs, BEs, defgss, defbfs, covBf, covBe, covAll, Zs, As, params, param_cov)
    !     implicit none
    !     integer, dimension(:), intent(in) :: Zs, As
    !     integer, parameter :: macparams = num_params - 3
    !     real(r_kind), intent(in) :: params(num_params),param_cov(num_params,num_params)
    !     real(r_kind), intent(out), dimension(size(Zs), size(Zs)) :: covBf, covBe !!Covariance of fission barrier and binding energy
    !     real(r_kind), intent(out), dimension(2*size(Zs), 2*size(Zs)) :: covAll  !!Covariance of fission barrier and be in the same matrix. 
    !     real(r_kind), dimension(size(Zs), size(Zs)) :: covBeSp
    !     !!Order is first entries are cov for binding energies, last entries are cov for barriers. in same order as nuclei in Z, A
    !     !!So for N nuclei, covAll(1, N+1) is the covariance between the barrier and g.s. for this nucleus.
    !     real(r_kind), intent(out), dimension(size(Zs)) :: Bfs, BEs, defgss, defbfs
    !     real(r_kind), dimension(size(Zs), num_params) :: Jmat_Sp, Jmat_Be
    !     real(r_kind), dimension(num_params, size(Zs)) :: Jmat_SpT, Jmat_BeT
    !     real(r_kind), dimension(size(Zs), size(Zs)) :: covSp
    !     integer :: ii , Z, A, N, jj
    !     real(r_kind) :: Be, defgs, Bf, defbf, alpha0
    !     N = size(Zs)

    !     do ii = 1, N
    !         Z = Zs(ii)
    !         A = As(ii)
    !         call find_fiss_barr(Bf, defbf, Be, defgs, params, Z, A)
    !         Bfs(ii) = Bf
    !         defbfs(ii) = defbf
    !         BEs(ii) = Be
    !         defgss(ii) = defgs      
    !         alpha0 = sqrt(alpha0sq(params(7), params(4),A))
    !         write(*,*) "frac gs:", defgs/alpha0
    !         write(*,*) "frac bf:", defbf/alpha0
    !     end do

    !     !When deformations are found, create J_matrices
    !     Jmat_Be = J_mat(Zs, As, defgss, N, params)
    !     Jmat_Sp = J_mat(Zs, As, defbfs, N, params)
    !     write(*,*) "JvecBe:"
    !     write(*,'(8F10.3)') Jmat_Be
    !     write(*,*) "JvecSp:"
    !     write(*,'(8F10.3)') Jmat_Sp
    !     Jmat_SpT = transpose(Jmat_Sp)
    !     Jmat_BeT = transpose(Jmat_Be)

    !     covBe = matmul(Jmat_Be, matmul(param_cov, Jmat_BeT))
    !     covSp = matmul(Jmat_Sp, matmul(param_cov, Jmat_SpT))

    !     covAll(1:N,1:N) = covBe
        
    !     covBeSp = matmul(Jmat_Be, matmul(param_cov, Jmat_SpT))

    !     covBf = covSp + covBe - covBeSp - transpose(covBeSp) !!Covariance of fission barriers.
    !     write(*,*) "covSp:"
    !     write(*,*) covSp
    !     write(*,*) "covBe:"
    !     write(*,*) covBe
    !     write(*,*) "covBeSp:"
    !     write(*,*) covBeSp
    !     covAll(N+1:2*N, N+1:2*N) = covBf
    !     covAll(1:N, N+1:2*N) = covBeSp - covBe
    !     covAll(N+1:2*N, 1:N) = transpose(covAll(1:N, N+1:2*N))  
    ! end subroutine


    





end module

module stack_module
    use constants, only: r_kind, num_def_params
    implicit none

    private
    public :: Stack
    integer, parameter :: max_size = 1e5
    type :: Stack
        integer :: elements(max_size, num_def_params)
        integer, public :: top_index = 0
    contains
        procedure, public :: pop => stack_pop
        procedure, public :: push => stack_push
        procedure, public :: top => stack_top
        procedure, public :: is_empty => stack_is_empty
    end type

    contains
    subroutine stack_push(self, value)
        class(Stack), intent(inout) :: self
        integer, intent(in) :: value(num_def_params)
        if(self%top_index < max_size) then
            self%top_index = self%top_index + 1
            self%elements(self%top_index, :) = value
        else
            write(*,'(A)') "Stack overflow: Cannot push, stack is full"
        end if
    end subroutine


    subroutine stack_pop(self, value)
        class(Stack), intent(inout) :: self
        integer, intent(out) :: value(num_def_params)
        if(self%top_index > 0) then
            value = self%elements(self%top_index,:)
            self%top_index = self%top_index - 1
        else
            write(*,'(A)') "Stack underflow: Cannot pop, stack is empty"
            value = -1
        end if
    end subroutine

    function stack_top(self) result(value)
        class(Stack), intent(in) :: self
        integer :: value(num_def_params)
        if(self%top_index > 0) then
            value = self%elements(self%top_index,:)
        else
            write(*,'(A)') "Stack underflow: Cannot top, stack is empty"
            value = -1
        end if
    end function

    pure logical function stack_is_empty(self)
        class(Stack), intent(in) :: self
        stack_is_empty = self%top_index == 0
    end function

end module