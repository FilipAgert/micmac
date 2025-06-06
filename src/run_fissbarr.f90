program run_fissbarr
    use  constants, only: fitted_params, r_kind, num_params, cauchois_params
    !use fiss_barr
    use micmac, only: find_gs_multiple, def_func_ho
    use fitting, only: read_fit_exp_data, fit_param_cov, num_fit_vals, exp_A, exp_Z
    use mass_table, only: read_ame

    implicit none

    integer, parameter :: Z = 92, A = 238
    integer, parameter :: nn = 6
    real(r_kind) :: params(num_params), param_cov(num_params,num_params), temp(num_params)
    real(r_kind) :: defBf, Bf, betaBf
    real(r_kind), dimension(:), allocatable :: BE_fits, def_fits
    integer, dimension(0:nn) :: ZCns, ACNs
    real(r_kind), dimension(nn+1) :: BEs, BFs, defGss, defBfs
    real(r_kind), dimension(nn+1, nn+1) :: covBE, covBF
    real(r_kind), dimension(2*(nn+1), 2*(nn+1)) :: covAll
    integer :: n, i
    character(len=100) :: fmt
    type(def_func_ho) :: func
    params = fitted_params
    call func%setup(params, Z, A)
    call func%print()
!     call read_ame()
!     call read_fit_exp_data()


!     allocate(BE_fits(num_fit_vals))
!     allocate(def_fits(num_fit_vals))

!     call find_gs_multiple(BE_fits, def_fits, params, exp_Z, exp_A, num_fit_vals)
!     param_cov = fit_param_cov(params, num_fit_vals, exp_Z,exp_A, def_fits)
!     ! param_cov = reshape([ &
! !     0.0012450119d0, 0.0044954844d0, 0.0001273946d0, -0.0001482834d0, 0.0007801377d0, 0.0002802700d0, -0.0002276940d0, &
! !     0.0044954844d0, 0.0170887151d0, 0.0004767682d0, -0.0005115828d0, 0.0032982450d0, 0.0017919812d0, -0.0014759566d0, &
! !     0.0001273946d0, 0.0004767682d0, 0.0000244730d0, -0.0000129862d0, 0.0001483606d0, 0.0000459051d0, -0.0000394312d0, &
! !    -0.0001482834d0, -0.0005115828d0, -0.0000129862d0, 0.0000188338d0, -0.0000756432d0, -0.0000157301d0, 0.0000145733d0, &
! !     0.0007801377d0, 0.0032982450d0, 0.0001483606d0, -0.0000756432d0, 0.0206094480d0, 0.0011127571d0, -0.0020027538d0, &
! !     0.0002802700d0, 0.0017919812d0, 0.0000459051d0, -0.0000157301d0, 0.0011127571d0, 0.0008791023d0, -0.0007569549d0, &
! !    -0.0002276940d0, -0.0014759566d0, -0.0000394312d0, 0.0000145733d0, -0.0020027538d0, -0.0007569549d0, 0.0012938000d0 ], &
! !     shape(param_cov))

    
    
!     ZCns = Z
!     ACns = A
!     do n = 1,nn
!         ACns(n) = ACns(n) - n
!     end do


    

    
!     call cov_barr_gs(BFs, BEs, defgss, defBfs, covBf, covBe, covAll, ZCNs, ACns, params, param_cov)

!     write(*,'(A,I3,A,I3)') "Z = ", Z, ", A = ", A

!     write(*,*)
!     write(*,*) "BE/ Bf covariance matrix:"
!     write(*,'(A)', advance='no') "n   "
!     do n = 0,nn
!         write(*,'(5x, I1,4x)',advance="no") n
!     end do
!     do n = 0,nn
!         write(*,'(5x,I1,4x)',advance="no") n
!     end do
!     write(*,*)

!     write(fmt,'(A,I5,A)') "(", 2*(nn+1),"F10.3)"
!     do n = 1, (2*nn+2)
!         write(*,fmt) covAll(n,:)
!     end do

!     write(*,*)
!     write(*,*)
!     write(*,*) "Deformation and barriers"
!     write(*,*) "A    BE (MeV)      Bf (MeV)  defGs     defBf"
!     do n = 0,nn
!         write(*,'(I3,2x,4F10.3)') ACns(n),   BEs(n+1),     BFs(n+1),      defgss(n+1),       defBfs(n+1)
!     end do
!     write(*,*)

!     ! write(*,*) "Parameter covariance matrix:"
!     ! do i = 1,num_params
!     !     write(*,'(8F10.6)') param_cov(i,:)
!     ! end do


!     write(*,*)


end program run_fissbarr