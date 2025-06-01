module minimise
    use constants, only:r_kind
    implicit none
    !!Module for finding ground state of a nucleus by minimising BE through deformation
    private
    public :: find_min
    interface
        function one_dim(x)
            import r_kind
            real(r_kind), intent(in) :: x
            real(r_kind) :: one_dim
        end function one_dim

        function n_dim(x)
            import r_kind
            real(r_kind), intent(in) :: x(:)
            real(r_kind) :: n_dim
        end function n_dim

    end interface
    contains

    subroutine find_min(x_min, f_min, iter, f, x0, x1, x_start)
        real(r_kind), intent(out) :: x_min, f_min
        real(r_kind), intent(in) :: x0, x1 !!Bounds on the minimisation
        real(r_kind), intent(in) :: x_start !!Initial guess
        real(r_kind), parameter :: tol = 1e-6 !!Tolerance for convergence
        procedure(one_dim) :: f
        real(r_kind) :: df, d2f
        real(r_kind), parameter :: dx = 1e-6 !!step size for numerical derivatives
        real(r_kind), parameter :: gamma = 0.5_r_kind !!Damping factor for Newton's method
        integer, parameter :: max_iter = 1000
        real(r_kind) :: x, xprev
        integer, intent(out) :: iter
        x = x_start
        do iter = 1, max_iter
   
            if (x < x0) x = x0
            if (x > x1) x = x1



            ! Calculate the function value and its derivatives
            df = (f(x + dx) - f(x - dx)) / (2.0_r_kind * dx)
            d2f = (f(x + dx) - 2.0_r_kind * f(x) + f(x - dx)) / (dx * dx)

            ! Update the position using Newton's method
            xprev = x
            x = x - gamma * df / d2f

            ! Check for convergence
            if (abs(x - xprev) < tol) exit

        end do
        x_min = x
        f_min = f(x)
        if (iter == max_iter) then
            print *, "Warning: Maximum iterations reached in find_min"
        end if
        ! Ensure the minimum is within bounds
        if (x_min < x0) then
            print*, "Warning: x_min is below x0, adjusting to x0"
            x_min = x0
            f_min = f(x0)
        else if (x_min > x1) then
            print*, "Warning: x_min is above x1, adjusting to x1"
            x_min = x1
            f_min = f(x1)
        end if
    end subroutine find_min

    function grad(f,x, dim)
        implicit none
        ! Calculate the gradient of a function at a point x
        real(r_kind), intent(in), dimension(dim) :: x
        integer, intent(in) :: dim
        procedure(n_dim) :: f
        real(r_kind) :: grad(dim)
        real(r_kind), parameter :: dx = 1e-6
        real(r_kind) :: dxv(dim)
        integer :: i
        dxv = 0.0_r_kind
        do i =1,dim
            dxv(i) = dx
            grad(i) = (f(x + dxv) - f(x - dxv)) / (2.0_r_kind * dx)
            dxv(i) = 0.0_r_kind  ! Reset the perturbation for the next dimension
        end do
    end function grad

    function hessian(f, x, dim)
        implicit none
        ! Calculate the Hessian matrix of a function at a point x
        real(r_kind), intent(in), dimension(dim) :: x
        integer, intent(in) :: dim
        procedure(n_dim) :: f
        real(r_kind), dimension(dim,dim) :: hessian, gradient(dim)
        real(r_kind), parameter :: dx = 1e-6
        real(r_kind) :: dxv(dim), dxv2(dim)
        integer :: i, j

        dxv = 0.0_r_kind
        hessian = 0.0_r_kind
        do i = 1, dim

            dxv(i) = dx
            dxv2 = 0.0_r_kind  ! Copy dxv for the second derivative calculationd
            do j = 1, dim
                if (i == j) then
                    hessian(i,j) = (f(x + dxv) - 2.0_r_kind * f(x) + f(x - dxv)) / (dx * dx)
                else
                    dxv2(j) = dx
                    hessian(i,j) = (f(x + dxv + dxv2) - f(x+dxv-dxv2) - f(x-dxv + dxv2) + f(x-dxv-dxv2)) / (4.0_r_kind * dx*dx)
                    dxv2(j) = 0.0_r_kind  ! Reset the perturbation for the next dimension
                end if
            end do
            dxv(i) = 0.0_r_kind
        end do
        
    end function hessian


end module minimise