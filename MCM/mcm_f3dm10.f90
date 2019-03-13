SUBROUTINE mcm_f3dm10(cm,i)
!************************************************************************
!
!    Purpose: Elastic plastic hydrodynamic material model
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
integer:: i,ispall
real(kind=real_acc) :: cm(*)
real(kind=real_acc) :: deltap,gg,g2,qs,qh,a1,a2,epf,aj2,sj2,ak,akt,seps
real(kind=real_acc) :: davg,davg1,defps,ywh,scle,thrg
real(kind=real_acc) :: p_incr, othird
!
othird = 1.0_d / 3.0_d
gg     = cm(1)				!shear modulus
qs     = cm(2)				!initial yield stress
qh     = cm(3)				!plastic modulus
!	pcut=cm(4)				!initial cutoff pressure
a1     = cm(5)				!linear pressure hardening coefficient
a2     = cm(6)				!quadratic pressure hardening coefficient
ispall = cm(7)			!spall type
epf    = cm(8)				!effective plastic strain at failure
g2     = 2.*mcm_dt*gg				!2*G*dt
!
!	davg(i)=-third*(d1(i)+d2(i)+d3(i))
davg   = othird*par(i)%tracerod
par(i)%p   =-othird*(par(i)%sigma(1,1)+par(i)%sigma(2,2)+par(i)%sigma(3,3))
!
akt      = qh
par(i)%p_cut = par(i)%pcut   
ak       = qs + qh*par(i)%efps + (a1 + a2*par(i)%p)*max(0.0_d,par(i)%p)	!old yield stress
!
!trial deviatoric stress values
!
par(i)%s(1,1) = par(i)%sigma(1,1) + par(i)%p + g2*(par(i)%rod(1,1) + davg)
par(i)%s(2,2) = par(i)%sigma(2,2) + par(i)%p + g2*(par(i)%rod(2,2) + davg)
par(i)%s(3,3) = par(i)%sigma(3,3) + par(i)%p + g2*(par(i)%rod(3,3) + davg)
par(i)%s(1,2) = par(i)%sigma(1,2) + 0.5*g2*par(i)%rod(1,2)
par(i)%s(1,3) = par(i)%sigma(1,3) + 0.5*g2*par(i)%rod(1,3)
par(i)%s(2,3) = par(i)%sigma(2,3) + 0.5*g2*par(i)%rod(2,3)
!
if (ispall.ne.3) go to 100
!
deltap   = par(i)%p - par(i)%p_cut
!
scle     = .50*(1. + sign(1.0,deltap))
!
par(i)%pcut  = scle*par(i)%p_cut
par(i)%s(1,1) = scle*par(i)%s(1,1)
par(i)%s(2,2) = scle*par(i)%s(2,2)
par(i)%s(3,3) = scle*par(i)%s(3,3)
par(i)%s(1,2) = scle*par(i)%s(1,2)
par(i)%s(1,3) = scle*par(i)%s(1,3)
par(i)%s(2,3) = scle*par(i)%s(2,3)
par(i)%s(2,1) = par(i)%s(1,2) ! tom
par(i)%s(3,1) = par(i)%s(1,3)
par(i)%s(3,2) = par(i)%s(2,3)
!
par(i)%p_cut = par(i)%pcut
!
100 continue
aj2 = .5*(par(i)%s(1,1)*par(i)%s(1,1) + par(i)%s(2,2)*par(i)%s(2,2) + par(i)%s(3,3)*par(i)%s(3,3)) + par(i)%s(1,2)*par(i)%s(1,2)&
    + par(i)%s(1,3)*par(i)%s(1,3) + par(i)%s(2,3)*par(i)%s(2,3)		! aj2=(Sij*Sij/2)
sj2 = sqrt(3.*aj2)					! effective trial stress sj2=sqrt(3*Sij*Sij/2)
!
!	  calculate effective plastic strain increment and efective plastic strain
!		efps=(sj2-current yield stress)/(3*G+Ep)
if (par(i)%fail.eq.1.0) then
   thrg  = othird/gg
   defps = dim(sj2,ak)/(1.e-30+dim(1.0,-thrg*qh))
   defps = thrg*  min(sj2,defps)
   par(i)%efps = par(i)%efps + defps
endif
!
ak   =  max(ak+qh*defps,0.0)				!new yield stress - isotropic hardening
seps = 1.e-30*gg
scle = (ak+seps) / (max(ak,sj2)+seps)   ! scale factor
!
!	radial return stresses to the yield surface
par(i)%s(1,1) = scle*par(i)%s(1,1)
par(i)%s(2,2) = scle*par(i)%s(2,2)
par(i)%s(3,3) = scle*par(i)%s(3,3)
par(i)%s(1,2) = scle*par(i)%s(1,2)
par(i)%s(1,3) = scle*par(i)%s(1,3)
par(i)%s(2,3) = scle*par(i)%s(2,3)
par(i)%s(2,1) = par(i)%s(1,2)
par(i)%s(3,1) = par(i)%s(1,3)
par(i)%s(3,2) = par(i)%s(2,3)
!
!	set stresses to zero for failed particles, failure effective plastic strain based
if (epf.ne.0.0) then
   if (par(i)%efps.gt.epf) then
      !
	  par(i)%fail=0.0
	  par(i)%s(1,1)=0.
	  par(i)%s(2,2)=0.
	  par(i)%s(3,3)=0.
	  par(i)%s(1,2)=0.
	  par(i)%s(1,3)=0.
	  par(i)%s(2,3)=0.
	  par(i)%s(2,1)=0.
	  par(i)%s(3,1)=0.
	  par(i)%s(3,2)=0.
      !
   endif
endif
!
if(ispall.ne.2) goto 300
!
aj2   = .5*(par(i)%s(1,1)*par(i)%s(1,1) + par(i)%s(2,2)*par(i)%s(2,2) + par(i)%s(3,3)*par(i)%s(3,3))&
	  + par(i)%s(1,2)*par(i)%s(1,2) + par(i)%s(1,3)*par(i)%s(1,3) + par(i)%s(2,3)*par(i)%s(2,3)&
	  + 1.e-12
!      sj2(i)=sign1(i)*sign5(i)**2+sign2(i)*sign6(i)**2+sign3(i)*sign4(i)
!     *      **2-sign1(i)*sign2(i)*sign3(i)-2.*sign4(i)*sign5(i)*sign6(i)
sj2   = par(i)%s(1,1)*par(i)%s(2,3)*par(i)%s(2,3) + par(i)%s(2,2)*par(i)%s(1,3)*par(i)%s(1,3)&
      + par(i)%s(3,3)*par(i)%s(1,2)*par(i)%s(1,2) - par(i)%s(1,1)*par(i)%s(2,2)*par(i)%s(3,3)&
	  - 2.0*par(i)%s(1,2)*par(i)%s(1,3)*par(i)%s(2,3)
!
akt   =-sqrt(27.0/aj2)*sj2*0.5/aj2
akt   = sign(min(abs(akt),1.0),akt)
ywh   = acos(akt)*othird
!
sj2   = 2.*sqrt(aj2*othird)*cos(ywh)
!
davg1 = par(i)%p - sj2 - par(i)%p_cut	! maximum principal tensile stress > cotoff pressure
!
scle  = .50*(1. + sign(1.0,davg1))
!
par(i)%pcut  = scle*par(i)%p_cut
!
par(i)%s(1,1) = scle*par(i)%s(1,1)
par(i)%s(2,2) = scle*par(i)%s(2,2)
par(i)%s(3,3) = scle*par(i)%s(3,3)
par(i)%s(1,2) = scle*par(i)%s(1,2)
par(i)%s(1,3) = scle*par(i)%s(1,3)
par(i)%s(2,3) = scle*par(i)%s(2,3)
par(i)%s(2,1) = par(i)%s(1,2)
par(i)%s(3,1) = par(i)%s(1,3)
par(i)%s(3,2) = par(i)%s(2,3)
!
par(i)%p_cut = par(i)%pcut
!
300 continue
!
return
!
end
