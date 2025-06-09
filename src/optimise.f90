module optimise
    use constants, only:r_kind, gamma_damp_fac
    implicit none
    !!Module for finding ground state of a nucleus by minimising BE through deformation
    private
    public :: find_min,  rand_d, find_optimal_point, func_1d, func_nd, conj_grad_method

    type, abstract :: func_1d
    contains
        procedure(f_int), deferred :: eval 
    end type

    type, abstract :: func_nd
    contains
        procedure(f_mul), deferred :: eval

    end type

    abstract interface
        pure elemental real(r_kind) function f_int(self, x)
            import:: func_1d, r_kind
            class(func_1d), intent(in) :: self
            real(r_kind), intent(in) :: x
        end function

        pure real(r_kind) function f_mul(self, xs)
            import :: func_nd, r_kind
            class(func_nd), intent(in) :: self
            real(r_kind), intent(in) :: xs(:)
        end function
        
    end interface
    contains

    !!Finds m
    recursive subroutine find_min(x_min, f_min, foundmin, func, x0, x1, num_restarts, dim)    
        real(r_kind), intent(out) :: x_min(dim), f_min
        logical ,intent(out) :: foundmin !!True if minimum exists.
        real(r_kind), intent(in), dimension(dim) :: x0, x1 !!Bounds on the minimisation
        integer, intent(in) :: num_restarts !!how many times to minimize? Ensures to find the best local min in the bounds
        class(func_nd), intent(in) :: func
        integer :: ii
        real(r_kind) :: init_pos(num_restarts,dim), xstart(dim), hess(dim,dim)
        real(r_kind) :: x(dim), y
        logical :: converged
        integer, intent(in) :: dim

        do ii = 1, dim

            init_pos(:,ii) = rand_ds(x0(ii),x1(ii),num_restarts)
        end do
        init_pos(1,:) = 0.0_r_kind

        foundmin = .false.
        f_min = huge(f_min)
        x_min = 0.0_r_kind
        do ii = 1, num_restarts !!Compute multiple times with a random starting position to ensure the minimum of the local minima is picked as the g.s.
            ! write(*,*)
            ! write(*,*) "running..."
            xstart = init_pos(ii,:)
            call conj_grad_method_nd(x, y, converged, func, x0, x1, xstart,dim)
            hess = hessian(func,x,dim)
            ! write(*,*) "Converged: ", converged
            ! write(*,*) "xval: ", x
            ! write(*,*) "With yval:", y
            ! write(*,*) "Init pos:", init_pos(ii)
            ! write(*,*)
            if(y < f_min .and. converged) then !!Hess must be positive for it to be a local min
                f_min = y
                x_min = x
                foundmin = .true.
            endif
        end do

        if(.not. foundmin) then!!If not converged, restart.
            call find_min(x_min, f_min, foundmin, func, x0,x1, num_restarts - 1, dim)
        endif
        ! write(*,*) "Best: "
        ! write(*,*) "xval: ", x_min
        ! write(*,*) "With yval:", f_min
    end subroutine

    subroutine conj_grad_method(x_min, f_min, converged, func, x0, x1, x_start, tol) !!https://en.wikipedia.org/wiki/Conjugate_gradient_method
        real(r_kind), intent(out) :: x_min, f_min
        real(r_kind), intent(in) :: x0, x1 !!Bounds on the minimisation
        real(r_kind), intent(in) :: x_start !!Initial guess
        real(r_kind), intent(in):: tol !Tolerance for convergence
        class(func_1d), intent(in) :: func
        real(r_kind) :: df, d2f, fpl, fm, alpha, A, g
        real(r_kind), parameter :: dx = 1.0e-6_r_kind !!step size for numerical derivatives
        integer, parameter :: max_iter = 100
        real(r_kind) :: x
        logical, intent(out) :: converged
        integer :: iter
        x = x_start
        converged = .false.

        do iter = 1, max_iter
            if(x < x0 .or. x > x1) then
                exit
            endif
            fpl = func%eval(x + dx)
            fm = func%eval(x - dx)
            df = (fpl - fm) / (2.0_r_kind * dx)
            d2f = (fpl - 2.0_r_kind * func%eval(x) + fm) / (dx * dx)
            A = d2f

            g = df
            ! write(*,*) g, A, x
            if(abs(g) < tol) then
                converged = .true.
                exit
            endif
            alpha = g*g/(g*A*g)
            x = x - alpha * g
        end do

        x_min = x
        f_min = func%eval(x)
        if (iter == max_iter) then
            ! print *, "Warning: Maximum iterations reached in find_min"
            converged = .false.
        end if
        ! Ensure the minimum is within bounds
        if (x_min < x0 .or. x_min > x1) then
            converged = .false.
            ! print*, "Warning: X_min outside bounds with value: ", x_min
        endif
    end subroutine conj_grad_method


    subroutine conj_grad_method_nd(x_min, f_min, converged, func, lb, ub, x_start, dim) !!https://en.wikipedia.org/wiki/Conjugate_gradient_method
        real(r_kind), intent(out) :: x_min(dim), f_min
        integer, intent(in) :: dim !!Dimension of system f(x1,x2,...,xn)
        real(r_kind), intent(in) :: lb(dim), ub(dim) !!Bounds on the minimisation
        real(r_kind), intent(in) :: x_start(dim) !!Initial guess
        real(r_kind), parameter :: tol = 1e-9_r_kind !!Tolerance for convergence
        class(func_nd), intent(in) :: func
        real(r_kind), parameter :: dx = 1.0e-6_r_kind !!step size for numerical derivatives
        integer, parameter :: max_iter = 100
        real(r_kind) :: x(dim)
        logical, intent(out) :: converged
        integer :: iter
        real(r_kind) :: dfdx(dim), dfdx2(dim,dim),A(dim,dim), g(dim), alpha
        x = x_start
        converged = .false.


        do iter = 1, max_iter

            if(oob(x, lb, ub, dim)) exit !!break iteration if outside bounds

            dfdx = grad(func, x, dim)
            dfdx2 = hessian(func, x, dim)
            g = dfdx
            A = dfdx2

            if(sum(abs(g)) < tol) then
                converged = .true.
                exit
            endif
            
            alpha = dot_product(g,g)/(dot_product(g,matmul(A,g)))

            x = x - alpha * g
        end do

        x_min = x
        f_min = func%eval(x)
        if (iter == max_iter) then
            ! print *, "Warning: Maximum iterations reached in find_min"
            converged = .false.
        end if
        ! Ensure the minimum is within bounds
        if (oob(x,lb,ub,dim)) then
            converged = .false.
            ! print*, "Warning: X_min outside bounds with value: ", x_min
        endif
    end subroutine

    pure logical function oob(x, lb, ub, dim) !!Check if coordinate in N-dims is outside bounds
        real(r_kind), dimension(dim), intent(in) :: x, lb, ub
        integer, intent(in) :: dim
        integer :: ii
        oob = .false.
        do ii = 1, dim
            if(x(ii) .lt. lb(ii) .or. x(ii) .gt. ub(ii)) then
                oob = .true.
                exit
            endif
        end do
    end function

    subroutine find_optimal_point(x_min, f_min, converged, func, x0, x1, x_start) !!Finds nearest point with derivative 0
        real(r_kind), intent(out) :: x_min, f_min
        real(r_kind), intent(in) :: x0, x1 !!Bounds on the minimisation
        real(r_kind), intent(in) :: x_start !!Initial guess
        real(r_kind), parameter :: tol = 1e-6 !!Tolerance for convergence
        class(func_1d), intent(in) :: func
        real(r_kind) :: df, d2f, fpl, fm
        real(r_kind), parameter :: dx = 1e-6 !!step size for numerical derivatives
        integer, parameter :: max_iter = 100
        real(r_kind) :: x, xprev
        logical, intent(out) :: converged
        integer :: iter
        x = x_start
        converged = .false.
        do iter = 1, max_iter
   
            if (x < x0) x = x0
            if (x > x1) x = x1
            fpl = func%eval(x + dx)
            fm = func%eval(x-dx)
            ! Calculate the function value and its derivatives
            df = (fpl - fm) / (2.0_r_kind * dx)
            d2f = (fpl - 2.0_r_kind * func%eval(x) + fm) / (dx * dx)

            ! Update the position using Newton's method
            xprev = x
            x = x - gamma_damp_fac * df / d2f

            ! Check for convergence
            if (abs(x - xprev) < tol) then
                converged = .true.
                exit
            endif

        end do
        x_min = x
        f_min = func%eval(x)
        if (iter == max_iter) then
            print *, "Warning: Maximum iterations reached in find_min"
        end if
        ! Ensure the minimum is within bounds
        if (x_min < x0) then
            print*, "Warning: x_min is below x0, adjusting to x0"
            x_min = x0
            f_min = func%eval(x0)
        else if (x_min > x1) then
            print*, "Warning: x_min is above x1, adjusting to x1"
            x_min = x1
            f_min = func%eval(x1)
        end if
    end subroutine find_optimal_point

    

    function grad(f,x, dim)
        implicit none
        ! Calculate the gradient of a function at a point x
        real(r_kind), intent(in), dimension(dim) :: x
        integer, intent(in) :: dim
        class(func_nd), intent(in) :: f
        real(r_kind) :: grad(dim)
        real(r_kind), parameter :: dx = 1e-6
        real(r_kind) :: dxv(dim)
        integer :: i
        dxv = 0.0_r_kind
        do i =1,dim
            dxv(i) = dx
            grad(i) = (f%eval(x + dxv) - f%eval(x - dxv)) / (2.0_r_kind * dx)
            dxv(i) = 0.0_r_kind  ! Reset the perturbation for the next dimension
        end do
    end function grad

    function hessian(f, x, dim)
        implicit none
        ! Calculate the Hessian matrix of a function at a point x
        real(r_kind), intent(in), dimension(dim) :: x
        integer, intent(in) :: dim
        class(func_nd), intent(in) :: f
        real(r_kind), dimension(dim,dim) :: hessian
        real(r_kind), parameter :: dx = 1e-6
        real(r_kind) :: dxv(dim), dxv2(dim)
        integer :: i, j

        dxv = 0.0_r_kind
        hessian = 0.0_r_kind
        do i = 1, dim

            dxv(i) = dx
            dxv2 = 0.0_r_kind  
            do j = 1, dim
                if (i == j) then
                    hessian(i,j) = (f%eval(x + dxv) - 2.0_r_kind * f%eval(x) + f%eval(x - dxv)) / (dx * dx)
                else
                    dxv2(j) = dx
                    hessian(i,j) = (f%eval(x + dxv + dxv2) - f%eval(x+dxv-dxv2) - f%eval(x-dxv + dxv2) + f%eval(x-dxv-dxv2)) / (4.0_r_kind * dx*dx)
                    dxv2(j) = 0.0_r_kind 
                end if
            end do
            dxv(i) = 0.0_r_kind
        end do
        
    end function hessian


    function rand_d(a, b)
        real(r_kind), intent(in) :: a, b !!Bounds
        real(r_kind) :: rand_d, r

        call random_number(r)

        rand_d = a + r * (b-a)
    end function

    function rand_ds(a,b,num)
        real(r_kind), intent(in) :: a, b !!Bounds
        real(r_kind) :: rand_ds(num)
        integer, intent(in) :: num
        integer :: i
        do i = 1,num
            rand_ds(i) = rand_d(a,b)
        end do

    end function

    function linspace(a,b,N) result(xs)
        real(r_kind), intent(in) :: a,b
        integer, intent(in) :: N
        real(r_kind) :: step
        real(r_kind) :: xs(N)
        integer :: ii
        step = (b-a)/(N-1)
        
        do ii = 1, N
            xs(ii) = a + step*(ii-1)
        end do
    end function


end module optimise