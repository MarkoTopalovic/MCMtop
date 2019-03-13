subroutine mcm_f3dm3 (cm,i)
!************************************************************************
!
!    Purpose: Elastic-plastic material with isotropic and 
!			  kinematic hardening
!
!  Called by: constitutive
!
!       Date: 25-01-99
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
integer:: i, efps_flag
real(kind=real_acc) :: cm(5,*)
!
real(kind=real_acc) ::defps,des,depn,othird,press
real(kind=real_acc) :: davg,ym,pr,gg,gdt,gd2,blk,qh,qb,qs,qbqh,k,k2,j2,j1
real(kind=real_acc) :: efps_f
real(kind=real_acc) :: p_incr,q1,q2,q3,q4,q5,q6,scle,fac1,fac2,fac3
!
othird = 1.0_d / 3.0_d
qh     = cm(2,4)			! hardening modulus
qb   = cm(1,5)			! beta
qs   = cm(1,3)			! initial yield stress
ym   = cm(1,1)			! Young's modulus
pr   = cm(1,2)			! Poison's ratio
efps_f = cm(2,6)            ! Effective plastic strain at failure
gg     = 2.0_d*cm(5,4)		! 2.0*shear modulus
gdt  = mcm_dt*gg
gd2  = .5*gdt
qbqh = qb*qh			! beta*hardening modulus
!
efps_flag = nint(cm(1,6))
!
!	calculate pressure increment and pressure 
!
blk    = mcm_dt*ym / ((1.-2.*pr))		! Bulk modulus *3 *dt
davg   = othird*par(i)%tracerod
p_incr = blk*davg				! Presure increment
par(i)%p   = par(i)%p - p_incr				! Pressure
!
!	particle internal energy increment befor the stress update 
!
par(i)%einc = (par(i)%rod(1,1)*(par(i)%sigma(1,1)+par(i)%q(1,1)) + par(i)%rod(2,2)*(par(i)%sigma(2,2)+par(i)%q(2,2))            &
             + par(i)%rod(3,3)*(par(i)%sigma(3,3)+par(i)%q(3,3)) + 2.0*(par(i)%rod(1,2)*(par(i)%sigma(1,2)+par(i)%q(1,2))       &
	         + par(i)%rod(1,3)*(par(i)%sigma(1,3)+par(i)%q(1,3)) + par(i)%rod(2,3)*(par(i)%sigma(2,3)+par(i)%q(2,3)))) * mcm_dt
!
!   compute trial stress
!
par(i)%sigma(1,1) = par(i)%sigma(1,1) + p_incr + gdt*(par(i)%rod(1,1)-davg)
par(i)%sigma(2,2) = par(i)%sigma(2,2) + p_incr + gdt*(par(i)%rod(2,2)-davg)
par(i)%sigma(3,3) = par(i)%sigma(3,3) + p_incr + gdt*(par(i)%rod(3,3)-davg)
!      sigma(1,2,i)=sigma(1,2,i)+gd2*rod(1,2,i) 		!+0.5*dt*gg*rod(1,2,i)*2.0
!      sigma(1,3,i)=sigma(1,3,i)+gd2*rod(1,3,i)			!+0.5*dt*gg*rod(1,3,i)*2.0
!      sigma(2,3,i)=sigma(2,3,i)+gd2*rod(2,3,i)			!+0.5*dt*gg*rod(2,3,i)*2.0
par(i)%sigma(1,2) = par(i)%sigma(1,2) + gdt*par(i)%rod(1,2) 			
par(i)%sigma(1,3) = par(i)%sigma(1,3) + gdt*par(i)%rod(1,3)			
par(i)%sigma(2,3) = par(i)%sigma(2,3) + gdt*par(i)%rod(2,3)			
par(i)%sigma(2,1) = par(i)%sigma(1,2)
par(i)%sigma(3,1) = par(i)%sigma(1,3)
par(i)%sigma(3,2) = par(i)%sigma(2,3)
!
k=qs+qbqh*par(i)%efps	! yield stress update
!
!	deviatoric trial stress components
!
par(i)%s(1,1) = par(i)%p + par(i)%sigma(1,1) - par(i)%alfa(1,1)
par(i)%s(2,2) = par(i)%p + par(i)%sigma(2,2) - par(i)%alfa(2,2)
par(i)%s(3,3) = par(i)%p + par(i)%sigma(3,3) - par(i)%alfa(3,3)
par(i)%s(1,2) = par(i)%sigma(1,2) - par(i)%alfa(1,2)
par(i)%s(1,3) = par(i)%sigma(1,3) - par(i)%alfa(1,3)
par(i)%s(2,3) = par(i)%sigma(2,3) - par(i)%alfa(2,3)
par(i)%s(2,1) = par(i)%s(1,2)
par(i)%s(3,1) = par(i)%s(1,3)
par(i)%s(3,2) = par(i)%s(2,3)
!
!	J2 stress invariant
!
j2   = par(i)%s(2,3)*par(i)%s(2,3) + par(i)%s(1,3)*par(i)%s(1,3) + par(i)%s(1,2)*par(i)%s(1,2) &
     - par(i)%s(1,1)*par(i)%s(2,2) - par(i)%s(2,2)*par(i)%s(3,3) - par(i)%s(1,1)*par(i)%s(3,3)
!
!	k2  von Mises yield function
!
k2   = 3.0*j2-k*k	
!
!	calculate scale factors	  	
!
scle = 0.5*(1.0+sign(1.0,k2))
fac1 = 1.0/(1.5*gg+qh)
fac2 = 1.5*gg
fac3 = (1.-qb)*qh
!
!	j1  effective stress
!
j1   = sqrt(3.0*abs(j2))+1.0-scle
!
!	calculate effective plastic strain increment and efective plastic strain
!
defps   = scle*fac1*(j1-k)
par(i)%efps = par(i)%efps+defps
!	
des     = scle*fac2*defps/j1			!stress correction factor
depn    = scle*fac3*defps/j1			!yield surface translation
!
par(i)%sigma(1,1) = par(i)%sigma(1,1) - des*par(i)%s(1,1)
par(i)%sigma(2,2) = par(i)%sigma(2,2) - des*par(i)%s(2,2)
par(i)%sigma(3,3) = par(i)%sigma(3,3) - des*par(i)%s(3,3)
par(i)%sigma(1,2) = par(i)%sigma(1,2) - des*par(i)%s(1,2)
par(i)%sigma(1,3) = par(i)%sigma(1,3) - des*par(i)%s(1,3)
par(i)%sigma(2,3) = par(i)%sigma(2,3) - des*par(i)%s(2,3)
par(i)%sigma(2,1) = par(i)%sigma(1,2)
par(i)%sigma(3,1) = par(i)%sigma(1,3)
par(i)%sigma(3,2) = par(i)%sigma(2,3)
!
par(i)%alfa(1,1) = par(i)%alfa(1,1) + depn*par(i)%s(1,1)
par(i)%alfa(2,2) = par(i)%alfa(2,2) + depn*par(i)%s(2,2)
par(i)%alfa(3,3) = par(i)%alfa(3,3) + depn*par(i)%s(3,3)
par(i)%alfa(1,2) = par(i)%alfa(1,2) + depn*par(i)%s(1,2)
par(i)%alfa(1,3) = par(i)%alfa(1,3) + depn*par(i)%s(1,3)
par(i)%alfa(2,3) = par(i)%alfa(2,3) + depn*par(i)%s(2,3)
par(i)%alfa(2,1) = par(i)%alfa(1,2)
par(i)%alfa(3,1) = par(i)%alfa(1,3)
par(i)%alfa(3,2) = par(i)%alfa(2,3)
!
! Test implementation of efps filure (cf dyna model 13)
!
if(efps_flag.eq.1) then
 if(par(i)%efps.gt.efps_f) then
  press = othird*(par(i)%sigma(1,1)+par(i)%sigma(2,2)+par(i)%sigma(3,3))
  par(i)%sigma(1,1) = press
  par(i)%sigma(1,2) = 0.0_d
  par(i)%sigma(1,3) = 0.0_d
  par(i)%sigma(2,2) = press
  par(i)%sigma(2,1) = 0.0_d
  par(i)%sigma(2,3) = 0.0_d
  par(i)%sigma(3,3) = press
  par(i)%sigma(3,1) = 0.0_d
  par(i)%sigma(3,2) = 0.0_d
 endif
endif
!
!	total particle internal energy increment after the stress update
!
par(i)%einc = 0.5*((par(i)%rod(1,1)*par(i)%sigma(1,1) + par(i)%rod(2,2)*par(i)%sigma(2,2) + par(i)%rod(3,3)* &
		            par(i)%sigma(3,3) + 2.0*(par(i)%rod(1,2)*par(i)%sigma(1,2) + par(i)%rod(1,3)* &
		            par(i)%sigma(1,3) + par(i)%rod(2,3)*par(i)%sigma(2,3)))*mcm_dt + par(i)%einc)
!
return
!
end
