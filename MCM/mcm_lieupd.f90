subroutine mcm_lieupd (nn)
!************************************************************************
!
!    Purpose: Calculates particle and material internal energy
!
!  Called by: CONSTITUTIVE
!
!       Date: 06-01-99
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
integer:: nn,nn1,davg
real(kind=real_acc) :: qp,qpm
real(kind=real_acc) :: sss(3,3),ss(9),ppp(3,3),pp(9) 
!
sss = par(nn)%q(1:3,1:3)
ppp = par(nn)%rod(1:3,1:3)
pp  = pack(ppp,.true.)
ss  = reshape(source=sss,shape=shape(ss))
qpm = dot_product(ss,pp)
!	qp=sum(qpm)
!
par(nn)%e = par(nn)%e + (par(nn)%einc - mcm_dt*qpm) * par(nn)%mass/par(nn)%rho
!
end
