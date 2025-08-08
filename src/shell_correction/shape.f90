module shape
    use constants
    use quad, only: LEGQUAD
    
    private
    public ::betadef, spherical_def, vol_fac, def_equals, surfRadius, Y20,Y40,dY20dtheta,d2Y20dtheta2,d3Y20dtheta3,dY40dtheta,d2Y40dtheta2,d3Y40dtheta3
    type :: betadef
        real(kind) :: beta2, beta4

        contains
        procedure :: eq => def_equals

    end type
    real(kind), private :: leg_x(nquad), leg_w(nquad)
    type(betadef), parameter :: spherical_def = betadef(beta2=0, beta4=0)

    contains

    subroutine compute_quad()
        logical, save :: first_time = .true.
        if(first_time) then
            call legquad(leg_x, leg_w)
            first_time = .false.
        endif
    end subroutine
    real(kind) function vol_fac(def)
        !!Gets V(def)/V(sph)
        type(betadef), intent(in) ::def
        integer :: ii
        real(kind) :: V, Vsph, x, theta
        call compute_quad()

        Vsph = 4.0*pi/3.0
        V = 0
        do ii =1,nquad !dt/dx = pi/2 => dt = dx pi/2
            x = leg_x(ii)
            theta = x*pi/2 + pi/2
            V = V + leg_w(ii)*surfRadius(theta,def,1.0_kind)**3 * sin(theta)
        end do
        V = V*pi*pi/3.0_kind 
        vol_fac = (Vsph/V)**(1.0_kind/3.0_kind)
        !Multiply 
    end function


    pure logical function def_equals(self, other) result(equality)
        class(betadef), intent(in) :: self,other
        equality = self%beta2==other%beta2 .and. self%beta4 == other%beta4 

    end function

    pure elemental real(kind) function surfRadius(theta, def, R0) !!computes distance to surface
        real(kind), intent(in) :: theta, R0
        type(betadef), intent(in) :: def
        surfRadius = R0 * (1.0_kind + def%beta2 * Y20(theta) + def%beta4*Y40(theta))
    end function


    pure elemental real(kind) function Y20(theta)
        real(kind), intent(in) :: theta
        Y20 = sqrt(5.0_kind / (16 * pi)) * (3.0_kind * cos(theta)**2 - 1.0_kind)
    end function

    pure elemental real(kind) function Y40(theta)
        real(kind) ,intent(in) :: theta
        Y40 = sqrt(9.0_kind / (256 * pi)) * (35.0_kind * cos(theta)**4 - 30.0_kind * cos(theta)**2 + 3)
    end function

    pure elemental real(kind) function dY20dtheta(theta) result(dydt)
        real(kind), intent(in) :: theta
        dydt = -0.75_kind * sqrt(5.0_kind/pi) * sin(2*theta)
    end function

    pure elemental real(kind) function d2Y20dtheta2(theta) result(dy2dt2)!!second
        real(kind), intent(in) :: theta
        dy2dt2 = -1.5_kind * sqrt(5.0_kind/pi) * cos(2*theta)
    end function

    pure elemental real(kind) function d3Y20dtheta3(theta) result(dy3dt3)
        real(kind), intent(in) :: theta
        dy3dt3 = 6* sqrt(5.0_kind/pi) * sin(theta)*cos(theta)
    end function

    pure elemental real(kind) function dY40dtheta(theta) result(dydt)
        real(kind), intent(in) :: theta
        dydt = 15*sin(theta)*cos(theta)*(3.0_kind-7.0_kind*cos(theta)*cos(theta))/(4.0_kind*sqrt(pi))
    end function

    pure elemental real(kind) function d2Y40dtheta2(theta) result(dy2dt2)!!second 
        real(kind), intent(in) :: theta
        real(kind)::ct2,st2
        ct2 = cos(theta)**2
        st2 = sin(theta)**2
        dy2dt2 = 15.0_kind*(3.0*st2+7.0*ct2*ct2 - 3.0*(7.0*st2+1)*ct2)/(4.0_kind*sqrt(pi))
    end function

    pure elemental real(kind) function d3Y40dtheta3(theta) result(dy3dt3)
        real(kind), intent(in) :: theta
        dy3dt3 = 15.0*(sin(2.0*theta)+14.0*sin(4.0*theta))/(4.0*sqrt(pi))
    end function


end module