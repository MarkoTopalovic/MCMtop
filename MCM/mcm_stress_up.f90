subroutine mcm_stress_up(cm,i)
!************************************************************************
!
!    Purpose: Stress update for the material models which use EOS
!
!  Called by: Constitutive
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
real(kind=real_acc) :: davg,pc,xmu
!
!      pc=cm(1,pmat(i))		!cutof pressure
!
par(i)%sigma(1,1) = par(i)%s(1,1) - par(i)%p	!stress tensor
par(i)%sigma(2,2) = par(i)%s(2,2) - par(i)%p
par(i)%sigma(3,3) = par(i)%s(3,3) - par(i)%p
par(i)%sigma(1,2) = par(i)%s(1,2)
par(i)%sigma(1,3) = par(i)%s(1,3)
par(i)%sigma(2,3) = par(i)%s(2,3)
par(i)%sigma(2,1) = par(i)%sigma(1,2)
par(i)%sigma(3,1) = par(i)%sigma(1,3)
par(i)%sigma(3,2) = par(i)%sigma(2,3)
!
return
!
end
