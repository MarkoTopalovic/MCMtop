SUBROUTINE mcm_corr_strain
!************************************************************************
!
!    Purpose: Calculate the strain rates for every particle
!
!  Called by: strain
!
!       Date: 09-08-2002
!
!     Errors: 
!
!      Notes: Uses a mixed correction approach
!             In this routine  x and v are held at the n-1/2 time level.
!             However h and rho are only known at the n-1 time level, as they
!             are only advance to the n time level later in the calculation
!             once tracerod is known. For this routine then, both rho and h are
!             kept at the n-1 time level.
!
!************************************************************************
!
use mcm_database
use mcm_math_func
!
IMPLICIT NONE
!
INTEGER :: i, j, n, m, k, l
REAL(kind=real_acc)    :: Volj, dwdx(mcm_ndim),havg,dwdr
real(kind=real_acc), dimension(3,3) :: grad_v, f
real(kind=real_acc), dimension(3) :: vi
!
real(kind=real_acc) :: massj, rhoj, hj, holdj
real(kind=real_acc), dimension(3) :: xj,vj
!
real(kind=real_acc) :: sum_w, sum_gradw(3), corr_gradw(3), B(3,3), invB(3,3), w
! 
!                                              
DO i=mcm_ssp,mcm_esp 
 !
 if(i.eq.10) then
  continue
 endif
 !
 grad_v = 0.0_d
 par(i)%rod = 0.0_d
 par(i)%spin = 0.0_d
 vi = par(i)%v
 !
 !-------------------------------------------------------------------------------------------------
 ! Loop over all neighbours to calculate the sum_w and sum_gradw terms for this neighbourhood
 !
 sum_w = 0.0_d
 sum_gradw = 0.0_d
 do k=1,par(i)%nnbr + par(i)%g_nnbr
  !
  call mcm_get_j_strain_info(i,k,xj,vj,massj,rhoj,hj,holdj)
  !
  Volj = massj/rhoj
  havg = 0.5_d*(par(i)%h+hj)
  !
  ! Note: kernel centred at i particle
  call mcm_kernel(w,par(i)%x,xj,havg)
  call mcm_gradw(dwdx,dwdr,par(i)%x,xj,havg)
  !
  sum_w = sum_w + Volj*w
  do l=1,mcm_ndim
   sum_gradw(l) = sum_gradw(l) + Volj*dwdx(l)
  enddo
  !
 enddo
 !
 !-------------------------------------------------------------------------------------------------
 ! Now calculate B tensor
 !
 B = 0.0_d
 !
 do k=1,par(i)%nnbr + par(i)%g_nnbr
  !
  call mcm_get_j_strain_info(i,k,xj,vj,massj,rhoj,hj,holdj)
  !
  Volj = massj/rhoj
  havg = 0.5_d*(par(i)%h+hj)
  !
  ! Note: kernel centred at i particle
  call mcm_kernel(w,par(i)%x,xj,havg)
  call mcm_gradw(dwdx,dwdr,par(i)%x,xj,havg)
  !
  do m=1,mcm_ndim
   corr_gradw(m) = (dwdx(m)*sum_w - w*sum_gradw(m)) / sum_w**2
  enddo
  !
  do m=1,mcm_ndim
   do n = 1,mcm_ndim
    B(n,m) = B(n,m) + Volj * (xj(n) - par(i)%x(n)) * corr_gradw(m)
   enddo
  enddo
  !
 enddo
 !
 ! Calculate inverse of B tensor
 !
 invB = inverse(B)
 !
 !-------------------------------------------------------------------------------------------------
 !
 ! Calculate strain rate tensor
 !
 f = 0.0_d
 !
 DO k=1,par(i)%nnbr + par(i)%g_nnbr
  !
  call mcm_get_j_strain_info(i,k,xj,vj,massj,rhoj,hj,holdj)
  !
  Volj = massj/rhoj
  havg = 0.5_d*(par(i)%h+hj)
  !
  !
  ! Note: kernel centred at i particle
  call mcm_kernel(w,par(i)%x,xj,havg)
  call mcm_gradw(dwdx,dwdr,par(i)%x,xj,havg)
  !
  !
  do m=1,mcm_ndim
   corr_gradw(m) = (dwdx(m)*sum_w - w*sum_gradw(m)) / sum_w**2
  enddo
  !
  corr_gradw=matmul(invB,corr_gradw)
  !
  DO n = 1,mcm_ndim
   DO m = 1,mcm_ndim
    !f(n,m) = f(n,m) + Volj*(vj(n)-par(i)%v(n))*corr_gradw(m)
    f(n,m) = f(n,m) + Volj*vj(n)*corr_gradw(m)
   ENDDO
  ENDDO
 ENDDO
 !
 grad_v = f
 !
 do n=1,mcm_ndim
  do m=1,mcm_ndim
   par(i)%rod(n,m)  = 0.5_d * (grad_v(n,m) + grad_v(m,n))
   par(i)%spin(n,m) = 0.5_d * (grad_v(n,m) - grad_v(m,n))
  enddo
 enddo
 !
 par(i)%tracerod=par(i)%rod(1,1)+par(i)%rod(2,2)+par(i)%rod(3,3)	!rate of deformation tensor trace
ENDDO
! 
!
END SUBROUTINE mcm_corr_strain
