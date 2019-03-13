SUBROUTINE mcm_oldstrain
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
!      Notes: Uses the 'standard' sph approach.
!             In this routine  x and v are held at the n-1/2 time level.
!             However h and rho are only known at the n-1 time level, as they
!             are only advance to the n time level later in the calculation
!             once tracerod is known. For this routine then, both rho and h are
!             kept at the n-1 time level.
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i, j, n, m, k
REAL(kind=real_acc)    :: Volj, dwdx(mcm_ndim),havg,dwdr
real(kind=real_acc), dimension(3,3) :: grad_v
real(kind=real_acc), dimension(3) :: vi
!
real(kind=real_acc) :: massj, rhoj, hj, holdj
real(kind=real_acc), dimension(3) :: xj,vj
! 
!                                              
DO i=mcm_ssp,mcm_esp 
 grad_v = 0.0_d
 par(i)%rod = 0.0_d
 par(i)%spin = 0.0_d
 vi = par(i)%v
 !
 DO k=1,par(i)%nnbr + par(i)%g_nnbr
  !
  call mcm_get_j_strain_info(i,k,xj,vj,massj,rhoj,hj,holdj)
  !
  Volj = massj/rhoj
  havg = 0.5_d*(par(i)%h+hj)
  !
  call mcm_gradw(dwdx,dwdr,par(i)%x,xj,havg)
  !
  DO n = 1,mcm_ndim
   DO m = 1,mcm_ndim
    grad_v(n,m) = grad_v(n,m) + Volj*(vj(n)-vi(n))*dwdx(m)
   ENDDO
  ENDDO
 ENDDO
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
END SUBROUTINE mcm_oldstrain
