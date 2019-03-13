SUBROUTINE mcm_rstrn(i)
!************************************************************************
!
!    Purpose: Rotate Backstress 
!
!  Called by: Constitutive
!
!       Date: 29-01-02
!
! Written by: Tom De Vuyst
!
!     Errors: 
!
!      Notes: Translated from DYNA3D to ensure correct rotation of back stress
!             mat52 - NOT YET FINISHED SHOULD USE ALFA(1:3,1:3,NP) instead of epx1..6 - TOM
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i, iorder=3
!
!     rotate back stresses
!
!
!  used as working space only
!
!
real(kind=real_acc) :: straux(3,3),qs(3,3),qqinc(3,3),qqaux(3,3)
real(kind=real_acc) :: wxx,wyy,wzz,dt2,q1,q2,q3,q4,q5,q6,q1s23,q2s13,q3s12
real(kind=real_acc) :: aa,bb,cc,dd,qm,qm2,wxx2,wyy2,wzz2,etol,bwxx,bwyy,bwzz
real(kind=real_acc) :: z1,z2,z3,z4
REAL(kind=real_acc) :: s1, s2, s3, s4, s5, s6
!
data etol/1.0e-3/
data z1/0.33333333333333333333/  ! 1/3
data z2/0.13333333333333333333/  ! 2/15
data z3/0.053968253968253971  /  ! 17/315
data z4/0.021869488536155203  /  ! 62/2835
!
s1 = par(i)%alfa(1,1)
s2 = par(i)%alfa(2,2)
s3 = par(i)%alfa(3,3)
s4 = par(i)%alfa(1,2)
s5 = par(i)%alfa(2,3)
s6 = par(i)%alfa(3,1)
!
!	calculate spin increments
!
dt2 = 0.5*mcm_dt     
!      wzz=dt2*(dyx(i)-dxy(i))		! the convention used in DYNA
!      wyy=dt2*(dxz(i)-dzx(i))
!      wxx=dt2*(dzy(i)-dyz(i))
wzz =-dt2*par(i)%spin(1,2)
wyy = dt2*par(i)%spin(1,3)
wxx =-dt2*par(i)%spin(2,3)
!
if (iorder.le.2) then
   !
   ! first order terms
   !
   q1 = 2.*s4*wzz
   q2 = 2.*s6*wyy
   q3 = 2.*s5*wxx
   !
   par(i)%alfa(1,1) = s1 - q1 + q2
   par(i)%alfa(2,2) = s2 + q1 - q3
   par(i)%alfa(1,2) = s4 + wzz*(s1-s2) + wyy*s5 - wxx*s6
   par(i)%alfa(2,3) = s5 + wxx*(s2-s3) + wzz*s6 - wyy*s4
   par(i)%alfa(3,1) = s6 + wyy*(s3-s1) + wxx*s4 - wzz*s5
   !
   ! second order terms  (must do first order terms as well)
   !
   q1 = wxx*wxx
   q2 = wyy*wyy
   q3 = wzz*wzz
   q4 = wxx*wyy
   q5 = wyy*wzz
   q6 = wzz*wxx
   !
   q1s23 = q1*(s2-s3)
   q2s13 = q2*(s1-s3)
   q3s12 = q3*(s1-s2)
   !
   par(i)%alfa(1,1) = par(i)%alfa(1,1) - q3s12 - q2s13 + q4*s4 - 2.0*q5*s5 + q6*s6
   par(i)%alfa(2,2) = par(i)%alfa(2,2) + q3s12 - q1s23 + q4*s4 + q5*s5 - 2.0*q6*s6
   par(i)%alfa(1,2) = par(i)%alfa(1,2) + q4*(0.5*(s1+s2) - s3) - 0.5*s4*(q1 + q2 + 4.0*q3) + 1.5*(q6*s5 + q5*s6)
   par(i)%alfa(2,3) = par(i)%alfa(2,3) + q5*(0.5*(s2+s3) - s1) - 0.5*s5*(4.0*q1 + q2 + q3) + 1.5*(q6*s4 + q4*s6)
60 par(i)%alfa(3,1) = par(i)%alfa(3,1) + q6*(0.5*(s1+s3) - s2) - 0.5*s6*(q1 + 4.0*q2 + q3) + 1.5*(q5*s4 + q4*s5)
   !
   !
elseif (iorder.eq.3) then
   !
   ! exponential map
   !
   wxx2 = wxx*wxx
   wyy2 = wyy*wyy
   wzz2 = wzz*wzz
   !
   qm2 = wxx2 + wyy2 + wzz2
   !
   qm = sqrt(qm2)
   !
   dd = tan(qm*0.5)
   !
   if (qm.lt.etol) then
      aa = qm2*0.25        !dummy variable
      cc = 0.5*(1. + aa*(z1 + aa*(z2 + aa*(z3 + aa*z4))))
   else
      cc = dd/qm
   endif
   !
   bb = 2.0*cc/(1.0 + dd*dd)  
   !
   aa = bb*cc  
   !
   ! form the rotation matrix
   !
   qqinc(1,1) = 1.0 - aa*(wyy2+wzz2)
   qqinc(2,2) = 1.0 - aa*(wxx2+wzz2)
   qqinc(3,3) = 1.0 - aa*(wxx2+wyy2)
   qqinc(1,2) = aa*wxx*wyy
   qqinc(1,3) = aa*wxx*wzz
   qqinc(2,3) = aa*wyy*wzz
   !
   bwxx = bb*wxx
   bwyy = bb*wyy
   bwzz = bb*wzz
   qqinc(2,1) = qqinc(1,2)
   qqinc(3,1) = qqinc(1,3)
   qqinc(3,2) = qqinc(2,3)
   !
   qqinc(1,2) = qqinc(1,2) - bwzz
   qqinc(2,1) = qqinc(2,1) + bwzz
   qqinc(1,3) = qqinc(1,3) + bwyy
   qqinc(3,1) = qqinc(3,1) - bwyy
   qqinc(2,3) = qqinc(2,3) - bwxx
   qqinc(3,2) = qqinc(3,2) + bwxx
   !
   !  do t=q.s.q^t  as  qs=q.s then  t=qs.q^t
   !
   qs(1,1) = qqinc(1,1)*s1 + qqinc(1,2)*s4 + qqinc(1,3)*s6
   qs(1,2) = qqinc(1,1)*s4 + qqinc(1,2)*s2 + qqinc(1,3)*s5
   qs(1,3) = qqinc(1,1)*s6 + qqinc(1,2)*s5 + qqinc(1,3)*s3
   qs(2,1) = qqinc(2,1)*s1 + qqinc(2,2)*s4 + qqinc(2,3)*s6
   qs(2,2) = qqinc(2,1)*s4 + qqinc(2,2)*s2 + qqinc(2,3)*s5
   qs(2,3) = qqinc(2,1)*s6 + qqinc(2,2)*s5 + qqinc(2,3)*s3
   qs(3,1) = qqinc(3,1)*s1 + qqinc(3,2)*s4 + qqinc(3,3)*s6
   qs(3,2) = qqinc(3,1)*s4 + qqinc(3,2)*s2 + qqinc(3,3)*s5
   qs(3,3) = qqinc(3,1)*s6 + qqinc(3,2)*s5 + qqinc(3,3)*s3
   !
   par(i)%alfa(1,1) = qs(1,1)*qqinc(1,1) + qs(1,2)*qqinc(1,2) + qs(1,3)*qqinc(1,3)
   par(i)%alfa(2,2) = qs(2,1)*qqinc(2,1) + qs(2,2)*qqinc(2,2) + qs(2,3)*qqinc(2,3)
   par(i)%alfa(1,2) = qs(1,1)*qqinc(2,1) + qs(1,2)*qqinc(2,2) + qs(1,3)*qqinc(2,3)
   par(i)%alfa(2,3) = qs(2,1)*qqinc(3,1) + qs(2,2)*qqinc(3,2) + qs(2,3)*qqinc(3,3)
   par(i)%alfa(3,1) = qs(1,1)*qqinc(3,1) + qs(1,2)*qqinc(3,2) + qs(1,3)*qqinc(3,3)
   par(i)%alfa(3,3) =-par(i)%alfa(1,1) - par(i)%alfa(2,2)
   par(i)%alfa(2,1) = par(i)%alfa(1,2)
   par(i)%alfa(3,2) = par(i)%alfa(2,3)
   par(i)%alfa(1,3) = par(i)%alfa(3,1)
   !
endif      
!
RETURN
!
ENDSUBROUTINE mcm_rstrn