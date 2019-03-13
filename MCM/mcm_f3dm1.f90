subroutine mcm_f3dm1 (cm,i)
!************************************************************************
!
!    Purpose: Isotropic elastic material
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
integer:: i
real(kind=real_acc) :: cm(5,*)
real(kind=real_acc) :: davg,ym,pr,gg,gdt,gd2,blk
real(kind=real_acc) :: p_incr, othird
!
othird = 1.0_d/3.0_d
ym=cm(1,1)		!Young's modulus
pr=cm(1,2)		!Poison's ratio
! gg=ym/(1.+pr)
gg=cm(5,4)		!shear modulus
gdt= mcm_dt*gg
! gd2=0.5*gdt
gd2=2.0_d*gdt
!
!	calculate volumetric strain increment and pressure 
!
! xm(i)=1./volo(i)
! vlrho(i)=crho*volo(i)
blk=mcm_dt*ym/((1.0_d-2.0_d*pr))		!Bulk modulus *3 *dt
davg=othird*par(i)%tracerod
p_incr=blk*davg				!Presure increment=volume strain *bulk modulus
!
!	particle internal energy increment befor the stress update 
!
par(i)%einc = (par(i)%rod(1,1)*(par(i)%sigma(1,1)+par(i)%q(1,1)) + par(i)%rod(2,2)*(par(i)%sigma(2,2)+par(i)%q(2,2))   &
          + par(i)%rod(3,3)*(par(i)%sigma(3,3)+par(i)%q(3,3)) + 2.0_d*(par(i)%rod(1,2)*(par(i)%sigma(1,2)+par(i)%q(1,2)) &
	      + par(i)%rod(1,3)*(par(i)%sigma(1,3)+par(i)%q(1,3)) + par(i)%rod(2,3)*(par(i)%sigma(2,3)+par(i)%q(2,3)))) * mcm_dt
!
!	update stress
!
par(i)%sigma(1,1) = par(i)%sigma(1,1) + p_incr + gd2*(par(i)%rod(1,1)-davg)
par(i)%sigma(2,2) = par(i)%sigma(2,2) + p_incr + gd2*(par(i)%rod(2,2)-davg)
par(i)%sigma(3,3) = par(i)%sigma(3,3) + p_incr + gd2*(par(i)%rod(3,3)-davg)
par(i)%sigma(1,2) = par(i)%sigma(1,2) + gd2*par(i)%rod(1,2) 			!+0.5*dt*gg*rod(1,2,i)*2.0
par(i)%sigma(1,3) = par(i)%sigma(1,3) + gd2*par(i)%rod(1,3)			!+0.5*dt*gg*rod(1,3,i)*2.0
par(i)%sigma(2,3) = par(i)%sigma(2,3) + gd2*par(i)%rod(2,3)			!+0.5*dt*gg*rod(2,3,i)*2.0
par(i)%sigma(2,1) = par(i)%sigma(1,2)
par(i)%sigma(3,1) = par(i)%sigma(1,3)
par(i)%sigma(3,2) = par(i)%sigma(2,3)
!
!	total particle internal energy increment after the stress update
!
par(i)%einc = 0.5_d*((par(i)%rod(1,1)*par(i)%sigma(1,1) + par(i)%rod(2,2)*par(i)%sigma(2,2) + par(i)%rod(3,3)* &
		            par(i)%sigma(3,3) + 2.0_d*(par(i)%rod(1,2)*par(i)%sigma(1,2) + par(i)%rod(1,3)* &
		            par(i)%sigma(1,3) + par(i)%rod(2,3)*par(i)%sigma(2,3)))*mcm_dt + par(i)%einc)
!
! Calculate pressure so it is available for post-processing
!
!p(i)=p(i)-p_incr
par(i)%p = -othird * (par(i)%sigma(1,1) + par(i)%sigma(2,2) + par(i)%sigma(3,3))
!
return
!
end subroutine mcm_f3dm1
