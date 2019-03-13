subroutine mcm_f3dm9(cm,i)
!************************************************************************
!
!    Purpose: Fluid material model
!
!  Called by: constitutive
!
!       Date: 07-01-98
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
real(kind=real_acc) :: davg,pc,xmu,othird
!
othird = 1.0_d / 3.0_d 
pc=cm(1,par(i)%mat)		!cutof pressure
xmu=cm(2,par(i)%mat)		!dynamic viscosity
!
if (xmu.eq.0.0) then
   !
   par(i)%s(1:3,1:3)=0.0	!deviatoric part of the stress tensor
   !
else
   !
   davg=-othird*par(i)%tracerod
   !
   par(i)%s(1,1)=xmu*(par(i)%rod(1,1)+davg)	!deviatoric part of the stress tensor
   par(i)%s(2,2)=xmu*(par(i)%rod(2,2)+davg)
   par(i)%s(3,3)=xmu*(par(i)%rod(3,3)+davg)
   par(i)%s(1,3)=xmu*par(i)%rod(1,2)
   par(i)%s(1,3)=xmu*par(i)%rod(1,3)
   par(i)%s(2,3)=xmu*par(i)%rod(2,3)
   par(i)%s(2,1)=par(i)%s(1,2)
   par(i)%s(3,1)=par(i)%s(1,3)
   par(i)%s(3,2)=par(i)%s(2,3)
   !
endif
!
return
!
end
