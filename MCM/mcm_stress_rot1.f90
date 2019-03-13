      subroutine mcm_stress_rot1(i)
!************************************************************************
!
!    Purpose: Rotate stresses into the current configuration
!     
!  Called by: constitutive
!
!       Date: 15-12-98
!
!     Errors: 
!
!      Notes: 
!             
!             
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer:: i,iorder=3
!
real(kind=real_acc) :: straux(3,3),qs(3,3),qqinc(3,3),qqaux(3,3)
real(kind=real_acc) :: wxx,wyy,wzz,dt2,q1,q2,q3,q4,q5,q6,q1s23,q2s13,q3s12
real(kind=real_acc) :: aa,bb,cc,dd,qm,qm2,wxx2,wyy2,wzz2,etol,bwxx,bwyy,bwzz
real(kind=real_acc) :: z1,z2,z3,z4
!
data etol/1.0e-3/
!data z1/0.33333333333333333333/  ! 1/3
!data z2/0.13333333333333333333/  ! 2/15
!data z3/0.053968253968253971  /  ! 17/315
!data z4/0.021869488536155203  /  ! 62/2835
!
z1 =  1.0_d /    3.0_d
z2 =  2.0_d /   15.0_d
z3 = 17.0_d /  315.0_d
z4 = 62.0_d / 2835.0_d
!
straux(1,1) = par(i)%sigma(1,1)
straux(2,2) = par(i)%sigma(2,2)
straux(3,3) = par(i)%sigma(3,3)
straux(1,2) = par(i)%sigma(1,2)
straux(2,3) = par(i)%sigma(2,3)
straux(3,1) = par(i)%sigma(3,1)
straux(2,1) = par(i)%sigma(2,1)
straux(3,2) = par(i)%sigma(3,2)
straux(1,3) = par(i)%sigma(1,3)
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
if(iorder.le.2) then
   !
   ! first order terms
   !
   !        q1=2.*s4(i)*wzz
   !        q2=2.*s6(i)*wyy
   !        q3=2.*s5(i)*wxx
   q1 = 2.0_d*straux(1,2)*wzz
   q2 = 2.0_d*straux(3,1)*wyy
   q3 = 2.0_d*straux(2,3)*wxx
   !
   !	sigma=sigma+spin^t*sigma+spin*sigma+(spin^t*sigma+spin*sigma)*spin
   !
   par(i)%sigma(1,1) = straux(1,1) - q1 + q2
   par(i)%sigma(2,2) = straux(2,2) + q1 - q3
   par(i)%sigma(3,3) = straux(3,3) - q2 + q3
   par(i)%sigma(1,2) = straux(1,2) + wzz*(straux(1,1) - straux(2,2)) &
   				                   + wyy*straux(2,3) - wxx*straux(3,1)
   par(i)%sigma(2,3) = straux(2,3) + wxx*(straux(2,2) - straux(3,3)) &
				                   + wzz*straux(3,1) - wyy*straux(1,2)
   par(i)%sigma(3,1) = straux(3,1) + wyy*(straux(3,3) - straux(1,1)) &
				                   + wxx*straux(1,2) - wzz*straux(2,3)
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
   q1s23 = q1*(straux(2,2)-straux(3,3))
   q2s13 = q2*(straux(1,1)-straux(3,3))
   q3s12 = q3*(straux(1,1)-straux(2,2))
   !
   par(i)%sigma(1,1) = par(i)%sigma(1,1) - q3s12 - q2s13 &
                        + q4*straux(1,2) - 2.0_d*q5*straux(2,3) + q6*straux(3,1)
   par(i)%sigma(2,2) = par(i)%sigma(2,2) + q3s12-q1s23 &
                        + q4*straux(1,2) + q5*straux(2,3) - 2.0_d*q6*straux(3,1)
   par(i)%sigma(3,3) = par(i)%sigma(3,3) + q2s13+q1s23 &
                        - 2.0_d*q4*straux(1,2) + q5*straux(2,3) + q6*straux(3,1) 
   par(i)%sigma(1,2) = par(i)%sigma(1,2) + q4*(0.5_d*(straux(1,1) + straux(2,2)) &
  					    - straux(3,3)) - 0.5_d*straux(1,2)*(q1+q2+4.0_d*q3) &
                        + 1.5_d*(q6*straux(2,3) + q5*straux(3,1))
   par(i)%sigma(2,3) = par(i)%sigma(2,3) + q5*(0.5_d*(straux(2,2) + straux(3,3)) &
					    - straux(1,1)) - 0.5_d*straux(2,3)*(4.0_d*q1+q2+q3) &
                        + 1.5_d*(q6*straux(1,2) + q4*straux(3,1))
   par(i)%sigma(3,1) = par(i)%sigma(3,1) + q6*(0.5_d*(straux(1,1) + straux(3,3)) &
					    - straux(2,2)) - 0.5_d*straux(3,1)*(q1+4.0_d*q2+q3) &
                        + 1.5_d*(q5*straux(1,2) + q4*straux(2,3))
   par(i)%sigma(2,1) = par(i)%sigma(1,2)
   par(i)%sigma(3,2) = par(i)%sigma(2,3)
   par(i)%sigma(1,3) = par(i)%sigma(3,1)
   !
elseif (iorder.eq.3) then
   !
   ! exponential map
   !
   wxx2 = wxx*wxx
   wyy2 = wyy*wyy
   wzz2 = wzz*wzz
   qm2  = wxx2+wyy2+wzz2
   qm   = sqrt(qm2)
   dd   = tan(qm*0.5_d)
   !
   if (qm.lt.etol) then
      !
	  aa = qm2*0.25_d        !dummy variable
      !        c=0.5*(1.+a*(z1+a*(z2+a*(z3+a*z4))))
      cc = 0.5_d*(1.0_d + aa*(z1 + aa*(z2 + aa*(z3 + aa*z4))))	! approximation for tg(aa)
   else
      !
	  cc = dd/qm	
   endif
   !
   bb = 2.0_d*cc / (1.0_d + dd*dd)  
   aa = bb*cc  
   !
   ! form the rotation matrix
   !
   qqinc(1,1) = 1.0_d - aa*(wyy2+wzz2)
   qqinc(2,2) = 1.0_d - aa*(wxx2+wzz2)
   qqinc(3,3) = 1.0_d - aa*(wxx2+wyy2)
   qqinc(1,2) = aa*wxx*wyy
   qqinc(1,3) = aa*wxx*wzz
   qqinc(2,3) = aa*wyy*wzz
   qqinc(2,1) = qqinc(1,2)
   qqinc(3,1) = qqinc(1,3)
   qqinc(3,2) = qqinc(2,3)
   !
   bwxx = bb*wxx	
   bwyy = bb*wyy
   bwzz = bb*wzz
   !
   qqinc(1,2) = qqinc(1,2) - bwzz
   qqinc(2,1) = qqinc(2,1) + bwzz
   qqinc(1,3) = qqinc(1,3) + bwyy
   qqinc(3,1) = qqinc(3,1) - bwyy
   qqinc(2,3) = qqinc(2,3) - bwxx
   qqinc(3,2) = qqinc(3,2) + bwxx
   !
   !  calculate t=q.s.q^t  as  qs=q.s then  t=qs.q^t
   !
   qs(1,1) = qqinc(1,1)*straux(1,1) + qqinc(1,2)*straux(1,2) + qqinc(1,3)*straux(3,1)
   qs(1,2) = qqinc(1,1)*straux(1,2) + qqinc(1,2)*straux(2,2) + qqinc(1,3)*straux(2,3)
   qs(1,3) = qqinc(1,1)*straux(3,1) + qqinc(1,2)*straux(2,3) + qqinc(1,3)*straux(3,3)
   qs(2,1) = qqinc(2,1)*straux(1,1) + qqinc(2,2)*straux(1,2) + qqinc(2,3)*straux(3,1)
   qs(2,2) = qqinc(2,1)*straux(1,2) + qqinc(2,2)*straux(2,2) + qqinc(2,3)*straux(2,3)
   qs(2,3) = qqinc(2,1)*straux(3,1) + qqinc(2,2)*straux(2,3) + qqinc(2,3)*straux(3,3)
   qs(3,1) = qqinc(3,1)*straux(1,1) + qqinc(3,2)*straux(1,2) + qqinc(3,3)*straux(3,1)
   qs(3,2) = qqinc(3,1)*straux(1,2) + qqinc(3,2)*straux(2,2) + qqinc(3,3)*straux(2,3)
   qs(3,3) = qqinc(3,1)*straux(3,1) + qqinc(3,2)*straux(2,3) + qqinc(3,3)*straux(3,3)
   !
   par(i)%sigma(1,1) = qs(1,1)*qqinc(1,1) + qs(1,2)*qqinc(1,2) + qs(1,3)*qqinc(1,3)
   par(i)%sigma(2,2) = qs(2,1)*qqinc(2,1) + qs(2,2)*qqinc(2,2) + qs(2,3)*qqinc(2,3)
   par(i)%sigma(3,3) = qs(3,1)*qqinc(3,1) + qs(3,2)*qqinc(3,2) + qs(3,3)*qqinc(3,3)
   par(i)%sigma(1,2) = qs(1,1)*qqinc(2,1) + qs(1,2)*qqinc(2,2) + qs(1,3)*qqinc(2,3)
   par(i)%sigma(2,3) = qs(2,1)*qqinc(3,1) + qs(2,2)*qqinc(3,2) + qs(2,3)*qqinc(3,3)
   par(i)%sigma(3,1) = qs(1,1)*qqinc(3,1) + qs(1,2)*qqinc(3,2) + qs(1,3)*qqinc(3,3)
   par(i)%sigma(2,1) = par(i)%sigma(1,2)
   par(i)%sigma(3,2) = par(i)%sigma(2,3)
   par(i)%sigma(1,3) = par(i)%sigma(3,1)
   !
endif      
!
!	calculate global rotation matrix
!
!		qqaux(1:3,1:3)=qq(i,1:3,1:3)
qqaux(1,1) = par(i)%qq(1,1)
qqaux(2,2) = par(i)%qq(2,2)
qqaux(3,3) = par(i)%qq(3,3)
qqaux(1,2) = par(i)%qq(1,2)
qqaux(1,3) = par(i)%qq(1,3)
qqaux(2,3) = par(i)%qq(2,3)
qqaux(3,2) = par(i)%qq(3,2)
qqaux(3,1) = par(i)%qq(3,1)
qqaux(2,1) = par(i)%qq(2,1)
!
par(i)%qq(1:3,1:3) = matmul(qqaux,qqinc)
!
return
!
end subroutine mcm_stress_rot1
