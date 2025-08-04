module quad
    !!Source https://github.com/FilipAgert/fquad
    use constants, only: kind => r_kind
    implicit none
    private
    public :: legquad, lagquad, herquad
    type, abstract :: func
        contains
        procedure(eval_interface), deferred :: eval
        procedure(eval_interface), deferred :: evaldfdx
    end type func

    type, abstract :: diffeq
        contains
        procedure(diffeq_int), deferred :: evaldydx
    end type


    integer :: NQUAD

    abstract interface
        function eval_interface(self, x) result(y)
            import :: func, kind
            class(func), intent(in) :: self
            real(kind=kind), intent(in) :: x
            real(kind=kind) :: y
        end function eval_interface

        pure  function diffeq_int(self, x, y)  result(dydx)
            import::diffeq, kind
            real(kind=kind), intent(in) :: x,y
            real(kind=kind) :: dydx
            class(diffeq), intent(in) :: self
        end function
    end interface


    type , extends(func) :: legendre
        integer :: n
        contains
        procedure ::eval => legendre_eval
        procedure :: evaldfdx => legendre_evaldfdx
    end type

    type , extends(func) :: laguerre
        integer :: n
        contains
        procedure ::eval => laguerre_eval
        procedure :: evaldfdx => laguerre_evaldfdx
    end type


    type , extends(func) :: hermite
        integer :: n
        contains
        procedure ::eval => hermite_eval
        procedure :: evaldfdx => hermite_evaldfdx
    end type

    type, extends(diffeq) :: laguerre_diff
        integer :: n
        contains
        procedure :: evaldydx => diffeqlag
    end type

    type, extends(diffeq) :: hermite_diff
        integer :: n
        contains
        procedure :: evaldydx => diffeqher
    end type

    interface
        function f2d(x,y)
            import kind
            real(kind=kind), intent(in) :: x,y
            real(kind=kind) :: f2d
        end function
    end interface


    contains


    !! > Subroutine for computing legendre quadrature weights and absiccas 
    !! > 
    !! > Integrate: I =  int f(x) dx with bounds [-1, 1]
    !! > is calculated as:
    !! > I = sum_i=1,n   leg_w(i) * f(leg_x(i))
    !! > 
    !! > n weights can integrate a polynomial of order 2n + 1 exactly.
    !! > If a 2n+1 polynomial can approximate f(x) well, legendre quadrature of order n is a good choice of quadrature points
    !! > 
    !! > Author: Filip Agert (2025)
    subroutine LEGQUAD(leg_x,leg_w)
        real(kind=kind), dimension(:), intent(out) :: leg_x !!integration absiccas
        real(kind=kind), dimension(:), intent(out) :: leg_w!!integration weights
        real(kind=kind) :: guess, x
        type(legendre) ::leg
        integer :: ub, k, lb, n
        n = size(leg_x)
        if(n /= size(leg_w)) then
            WRITE(*,*) "Error: Node and weight vector must be the same size. leg_x, leg_w: ", leg_x, leg_w
        endif
        NQUAD = n
        leg = legendre(n=n)
        leg_x =0
        if(mod(n,2) ==0) then
            ub = n/2
            lb = n/2 +1
        else
            ub = (n-1)/2
            lb = ub+2
        endif
        do k =1,ub
            guess =leg_root_approx(n,k-1) !Approximate root
            leg_x(ub-k+1) = find_root(leg, -guess) !Refine root computation
        end do
        
        do k = 1, ub
            leg_x(n-k + 1) = -leg_x(k) !!use mirror symmetry in the roots

        end do


        do k = 1, n !!Calculates weights
            x = leg_x(k)
            leg_w(k) = 2.0_kind/((1.0_kind-leg_x(k)**2)* Pnd(leg_x(k),n)**2)
        end do
    end subroutine

    !! > Subroutine for computing laguerre quadrature weights and absiccas 
    !! > 
    !! > 
    !! > Integrate: I = int f(x) e^-x dx with bounds [0, +infty]
    !! > is calculated as:
    !! > I = sum_i=1,n   lag_w(i) * f(lag_x(i))
    !! > 
    !! > n weights can integrate a polynomial of order 2n + 1 exactly.
    !! > If a 2n+1 polynomial can approximate f(x) well, laguerre quadrature of order n is a good choice of quadrature points
    !! > 
    !! > Author: Filip Agert (2025)
    subroutine LAGQUAD(lag_x,lag_w)
        real(kind=kind), dimension(:), intent(out) :: lag_x !!integration absiccas
        real(kind=kind), dimension(:), intent(out) :: lag_w! !integration weights
        real(kind=kind) :: guess, x, init_cond_y, init_cond_x
        real(kind=kind), parameter :: pi =acos(-1.0_kind)
        type(laguerre) ::lag
        type(laguerre_diff) :: lagdiff
        integer :: k,n
        n = size(lag_x)
        if(n /= size(lag_w)) then
            WRITE(*,*) "Error: Node and weight vector must be the same size. leg_x, leg_w: ", lag_x, lag_w
        endif
        NQUAD = n
        lag = laguerre(n=n)
        lagdiff = laguerre_diff(n=n)
        lag_x = 0
        guess = lag_root_approx(n,1) !First root
        lag_x(1) = find_root(lag, guess) !Refine root computation


        do k =2,n
            init_cond_y = lag_x(k-1)
            init_cond_x = pi /2.0_kind
            call solveRK(guess, -pi/2.0_kind, init_cond_x, init_cond_y, lagdiff)
            !write(*,'(A,F10.3,A,F10.3)', advance='no') "Guess: ", guess
            lag_x(k) = find_root(lag, guess) !Refine root computation
            !write(*,'(A,F10.3,A,F10.3)')", found root:", lag_x(k)
            if (k>1) then 
                if(lag_x(k) .le. lag_x(k-1)) then

                    write(*,*) "wrong zero found: ", lag_x(k-1:k)
                    stop
                endif
            endif
        end do
        
        do k = 1, n !!Calculates weights
            x = lag_x(k)
            lag_w(k) = x / ((n+1) * Ln(x,n+1) )/ ( (n+1) * Ln(x,n+1))
        end do
    end subroutine
    !! > Subroutine for computing hermite quadrature weights and absiccas 
    !! > UNFINISHED. Does not give high precision calculation
    !! > 
    !! > Integrate: I = int f(x) e^-x^2 dx with bounds [-infty, +infty]
    !! > is calculated as:
    !! > I = sum_i=1,n   her_w(i) * f(her_x(i))
    !! > 
    !! > n weights can integrate a polynomial of order 2n + 1 exactly.
    !! > If a 2n+1 polynomial can approximate f(x) well, hermite quadrature of order n is a good choice of quadrature points
    !! > 
    !! > Author: Filip Agert (2025)
    subroutine HERQUAD(her_x,her_w)
        real(kind=kind), dimension(:), intent(out) :: her_x !!integration absiccas
        real(kind=kind), dimension(:), intent(out) :: her_w! !integration weights
        real(kind=kind) :: guess, x, init_cond_y, init_cond_x
        real(kind=kind), parameter :: pi =acos(-1.0_kind)
        type(hermite) ::her
        type(hermite_diff) :: herdiff
        integer :: k, first,n
        n = size(her_x)
        if(n /= size(her_w)) then
            WRITE(*,*) "Error: Node and weight vector must be the same size. leg_x, leg_w: ", her_x, her_w
        endif
        NQUAD = n
        her = hermite(n=n)
        herdiff = hermite_diff(n=n)
        her_x = 0


        if((mod(n,2)) == 0) then
            first = n/2 + 1
        else
            first = n/2 + 2
            her_x(first - 1) = 0
        endif
        guess = her_root_approx(n,1) !First root
        her_x(n) = find_root(her, guess) !Refine root computation
        her_x(1)=-her_x(n)
        !write(*,*) "First root:", her_x(n), ", guess:", guess

        do k =n-1,first,-1
            init_cond_y = her_x(k+1)
            init_cond_x = -pi /2.0_kind
            call solveRK(guess, pi/2.0_kind, init_cond_x, init_cond_y, herdiff)
            !write(*,'(A,F10.3,A,F10.3)', advance='no') "Guess: ", guess
            her_x(k) = find_root(her, guess) !Refine root computation
            !write(*,'(A,E30.10,A,E20.10)')", found root:", her_x(k)!, " h_n(x)", Hn(her_x(k),n)
            her_x(n-k+1) = -her_x(k)

            if (k<n) then 
                if(her_x(k) .ge. her_x(k+1)) then

                    write(*,*) "wrong zero found: ", her_x(k-1:k)
                    stop
                endif
            endif
        end do
        
        do k = 1, n !!Calculates weights
            x = her_x(k)
            her_w(k) = sqrt(pi) * 2.0_kind**(n-1) *fac(n-1) / (n*Hn(x,n-1)**2)
        end do
        !WRITE(*,*) "WARNING: Hermite quadrature does not work yet."
    end subroutine

    real(kind=kind) function fac(n)
        logical, save :: firsttime = .true.
        integer :: n
        real(kind=kind), save :: f(0:1750)
        integer :: i
        if(kind < 16 .and. n > 165) write(*,*) "ERR: need quad precision for this N"
        f(0) = 1
        if(firsttime) then
            do i = 1,1750
                f(i) = f(i-1) * i
            end do
            firsttime = .false.
        endif
        fac = f(n)
    end function

    real(kind=kind) function find_root(f, x0) result(x)
        class(func), intent(in) :: f !!function to evaluate
        real(kind=kind), intent(in) :: x0 !!starting guess for root
        real(kind=kind) :: fx, dfdx, xprev

        integer :: ii
        integer, parameter :: max_iter = 10000
        real(kind=kind), parameter :: tol = 16.0_kind*epsilon(1.0_kind)
        real(kind) :: factor
        factor = 1.0
        x = x0
        do ii = 1, max_iter
            
            fx = f%eval(x)
            dfdx = f%evaldfdx(x)
            xprev = x
            if(dfdx .eq. 0.0_kind) dfdx = 1000
            x = x  -fx/ dfdx*0.25

            if(abs(x) > 1.0_kind) factor = abs(x)
            if(abs(x-xprev) < tol*factor) exit
        end do

        if (ii .ge. max_iter) then
            write(*,*) "Error: Max iterations reached in finding root. Stopping"
            STOP
        endif

    end function




    real(kind=kind) function hmn(x,n)
        ! Computes the value of the n-th Hermite polynomial at x
        integer, intent(in) :: n
        real(kind=kind), intent(in) :: x
        real(kind=kind), parameter :: pi = ACOS(-1.0_kind), rootrootpi = sqrt(sqrt(pi))
        real(kind=kind), dimension(-1:n) ::h
        integer ::i
        logical, save :: first_time = .true.
        real(kind=kind), save, allocatable ,dimension(:) :: c1, c2
        if(first_time) then
            allocate(c1(0:nquad-1), c2(0:NQUAD-1))
            do i = 0,NQUAD-1
                c1(i) = sqrt(2.0_kind/(i+1.0_kind))
                c2(i) = sqrt(real(i,kind)/(i+1.0_kind))
            end do
            first_time = .false.
        end if


        h(-1) = 0.0_kind
        h(0) = exp(-x*x / 2.0_kind) / rootrootpi
        if(n>0) then
            do i = 0, n-1
                h(i + 1) = c1(i) * x * h(i) - c2(i) * h(i-1)
            end do
        endif
        hmn = h(n)
    end function

    pure real(kind=kind) elemental function hn(x,n)
        ! Computes the value of the n-th Hermite polynomial at x
        integer, intent(in) :: n
        real(kind=kind), intent(in) :: x
        real(kind=kind), dimension(0:n) ::h
        integer ::i

  
        h(0) = 1 
        
        if(n>0) then
            h(1) = 2*x
            do i = 1, n-1
                h(i + 1) = 2 *x * h(i) - 2*i * h(i-1)
            end do
        endif
        hn = h(n)
    end function

    pure real(kind=kind) elemental function hnd(x,n)
        ! Computes the value of the derivative of the n-th Hermite polynomial at x
        integer, intent(in) :: n
        real(kind=kind), intent(in) :: x
        hnd = 2*n*hn(x,n-1)
    end function


    real(kind=kind) function hmnd(x,n)
        ! Computes the value of the derivative n-th modified Hermite polynomial at x
        integer, intent(in) :: n
        real(kind=kind), intent(in) :: x

        hmnd = 2*n*hmn(x,n-1) - x*hmn(x,n)
    end function

    function hermite_eval(self, x) result(y)
        class(hermite), intent(in) :: self
        real(kind=kind), intent(in) :: x
        real(kind=kind) :: y

        y = hmn(x,self%n)
    end function hermite_eval

    function hermite_evaldfdx(self, x) result(y)
        class(hermite), intent(in) :: self
        real(kind=kind), intent(in) :: x
        real(kind=kind) :: y

        y = hmnd(x,self%n)
    end function hermite_evaldfdx

    pure real(kind=kind) elemental function Lmn(x,n)
        ! Computes the value of the n-th modified laguerre polynomial at x
        integer, intent(in) :: n
        real(kind=kind), intent(in) :: x
        real(kind=kind), parameter :: pi = ACOS(-1.0_kind), rootrootpi = sqrt(sqrt(pi))
        real(kind=kind), dimension(-1:n) ::L
        integer ::i

        L(-1) = 0.0_kind
        L(0) = exp(-x/ 2.0_kind)
        if(n>0) then
            do i = 0, n-1
                L(i + 1) = ((2*i + 1.0_kind  - x)*L(i) - real(i,8)* L(i-1))/(i+1.0_kind)
            end do
        endif
        Lmn = L(n)
    end function

    pure real(kind=kind) elemental function Lmnd(x,n)
        ! Computes the derivative of the n-th modified laguerre polynomial at x
        integer, intent(in) :: n
        real(kind=kind), intent(in) :: x
        !!d/dx L_n = (nL_n - nL_n-1)    /x
        Lmnd = n*(Lmn(x,n) - Lmn(x,n-1))/x - 0.5 *Lmn(x,n) 
    end function

    pure real(kind=kind) elemental function Ln(x,n)
        ! Computes the value of the n-th laguerre polynomial at x
        integer, intent(in) :: n
        real(kind=kind), intent(in) :: x
        real(kind=kind), dimension(0:n) ::L
        integer ::i

        L(0) = 1
 
        if(n>0) then
            L(1) = 1-x
            do i = 1, n-1
                L(i + 1) = ((2*i + 1.0_kind  - x)*L(i) - real(i,8)* L(i-1))/(i+1.0_kind)
            end do
        endif
        Ln = L(n)
    end function

    pure function laguerre_eval(self, x) result(y)
        class(laguerre), intent(in) :: self
        real(kind=kind), intent(in) :: x
        real(kind=kind) :: y

        y = Lmn(x,self%n)
    end function laguerre_eval

        pure function laguerre_evaldfdx(self, x) result(y)
        class(laguerre), intent(in) :: self
        real(kind=kind), intent(in) :: x
        real(kind=kind) :: y

        y = Lmnd(x,self%n)
    end function laguerre_evaldfdx

    pure real(kind=kind) elemental function lag_root_approx(n, k) result(root)
        ! Approximate the k-th root of the n-th laguerre polynomial
        integer, intent(in) :: n, k
        real(kind=kind) :: invn, rho,bz, phi
        real(kind=kind), parameter :: pi =ACOS(-1.0_kind)
        if(k == 1) then
            invn = 1.0_kind/n
            root = invn + invn**2 * (n-1.0_kind)/2 - invn**3*(n**2 +3*n - 4.0_kind)/12 + invn**4*(7.0_kind*n**3+6*n**2+23*n-36)/144 ! + O(1/n^5)
            ! src: https://en.wikipedia.org/wiki/Laguerre_polynomials
        else
            rho = n + 0.5_kind
            bz = bessel_zero_approx(k)
            phi = bz/rho
            root = phi + (-0.25_kind)*(1.0_kind-phi*cotan(phi))/(2*phi) / rho**2 
        endif
    end function

    pure real(kind=kind) elemental function her_root_approx(n, k) result(root)
        ! Approximate the 1st root of the n-th hermite polynomial
        integer, intent(in) :: n, k
        real(kind=kind) :: airy_firstroot
        root = 0
        if(k == 1) then
            airy_firstroot = -2.33810741045977_kind
            root = sqrt(n*2+1.0_kind) + 2.0_kind**(-1.0_kind/3.0_kind) * airy_firstroot*(2.0*n+1)**(-1.0_kind/6.0_kind)
            ! https://en.wikipedia.org/wiki/Hermite_polynomials#Zeroes
        endif
    end function

    pure real(kind=kind) elemental function bessel_zero_approx(k) result(root)
        ! Approximate the k-th root of J_0 bessel function
        ! https://dlmf.nist.gov/10.21#vi
        integer, intent(in) :: k !!kth zero of bessel
        real(kind=kind), parameter :: pi =ACOS(-1.0_kind)
        real(kind=kind) :: a, inv8a
        a = (k - 0.25_kind) * pi
        inv8a = 1.0_kind/(8*a)
        root = a + inv8a - inv8a**3* 4.0_kind*31.0_kind / (3) + inv8a**5*32*(3779.0_kind)/(15) - inv8a**7 *(64*6277237.0_kind)/(105)

    end function

    !! ######################### LEGENDRE ########################
    pure real(kind=kind) elemental function Pn(x,n)
        ! Computes the value of the n-th legendre polynomial at x
        integer, intent(in) :: n
        real(kind=kind), intent(in) :: x
        real(kind=kind), dimension(0:n) ::P
        integer ::i

        P(0) = 1.0_kind

        if(n>0) then
            P(1) = x
            do i = 1, n-1
                P(i + 1) = ((2*i + 1)*x * P(i) - i * P(i-1))/(i+1)
            end do
        endif
        Pn = P(n)
    end function

    pure real(kind=kind) elemental function Pnd(x,n)
        ! Computes the value of derivative of the n-th legendre polynomial at x
        integer, intent(in) :: n
        real(kind=kind), intent(in) :: x

        Pnd = (n)*(x*Pn(x,n)-Pn(x,n-1)) / (x**2 - 1.0_kind)
    end function

    pure function legendre_eval(self, x) result(y)
        class(legendre), intent(in) :: self
        real(kind=kind), intent(in) :: x
        real(kind=kind) :: y

        y = Pn(x,self%n)
    end function legendre_eval

    pure function legendre_evaldfdx(self, x) result(y)
        class(legendre), intent(in) :: self
        real(kind=kind), intent(in) :: x
        real(kind=kind) :: y

        y = Pnd(x,self%n)
    end function legendre_evaldfdx

    pure real(kind=kind) elemental function leg_root_approx(n, k)
        ! Approximate the k-th root of the n-th Legendre polynomial via an asymptotic formula with error on order (1/n^5)
        integer, intent(in) :: n, k
        real(kind=kind), parameter :: pi = ACOS(-1.0_kind)
        real(kind=kind) :: theta
        integer :: intpart
        intpart = n/2 - k
        theta = pi * (intpart * 4.0_kind - 1.0_kind ) / (4*n+2)
        leg_root_approx = cos(theta) * (1 - (n-1.0_kind) / (8*n**3) - (39.0_kind - 28.0_kind/sin(theta)**2)/(384.0_kind*n**4))
    end function

    pure function diffeqlag(self, x, y) result(val)
        class(laguerre_diff), intent(in) :: self
        real(kind=kind), intent(in) :: x,y
        real(kind=kind) :: val

        val = -1.0_kind/(sqrt(r(self,y)/p(y)) + (-1.0_kind + 2*q(y)) / (2*p(y)) * sin(2.0_kind*x)/2 )
        contains

        pure real(kind=kind) function r(self, x)
            class(laguerre_diff), intent(in) :: self
            real(kind=kind), intent(in) :: x
            r = self%n
        end function
        pure real(kind=kind) function q(x)
            real(kind=kind), intent(in) :: x
            q = 1-x
        end function
        pure real(kind=kind) function p(x)
            real(kind=kind),intent(in) :: x
            p=x
        end function
    end function

        pure function diffeqher(self, x, y) result(val)
        class(hermite_diff), intent(in) :: self
        real(kind=kind), intent(in) :: x,y
        real(kind=kind) :: val

        val = -1.0_kind/(sqrt(r(self,y)/p(y)) + (q(y)) / (p(y)) * sin(2.0_kind*x)/2 )
        contains

        pure real(kind=kind) function r(self, x)
            class(hermite_diff), intent(in) :: self
            real(kind=kind), intent(in) :: x
            r = self%n*2
        end function
        pure real(kind=kind) function q(x)
            real(kind=kind), intent(in) :: x
            q = -2*x
        end function
        pure real(kind=kind) function p(x)
            real(kind=kind),intent(in) :: x
            p = 1
        end function
    end function

    !! Solve for y given y' = f(x,y) and y(x0)=y0
    subroutine solveRK(yval, xval,x0,y0,eq)
        real(kind=kind), intent(in) :: xval, x0, y0
        real(kind=kind), intent(out) :: yval
        class(diffeq), intent(in) :: eq
        integer, parameter :: steps = 100
        real(kind=kind) :: L, h, k(0:steps), y(0:steps), x
        integer ::i


        L = xval-x0
        h = L/steps
        y(0) = y0
        k(0) = eq%evaldydx(x0, y0)
        x = x0
        do i = 0,steps-1
            x = x + h
            k(i+1) = h * eq%evaldydx(x, y(i)+k(i))
            y(i+1) = y(i) + 0.5 * (k(i)+k(i+1))
        end do
        yval = y(steps)
    end subroutine
end module quad