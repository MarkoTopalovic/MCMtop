SUBROUTINE mcm_eos4p(i)
!************************************************************************
!
!		Purpose: Calculate pressure, dp/drho, dp/dE, c FOR GRUNEISEN EOS
!
!	  Called by: Eoscalc
!
!	     Author: Tom De Vuyst
!
!          Date: 02-08-2002
!
! Last Modified: 02-08-2002 by J. Campbell
!
!        Errors: 
!
!         Notes:
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i,j
REAL(kind=real_acc) :: mu,s1,s2,s3,c1,gamma,a1
REAL(kind=real_acc) :: mu11, mu21, mu22, mu32, mu33
REAL(kind=real_acc) :: termone, termtwo, eosnum, eosden, termthree, termfour
REAL(kind=real_acc) :: dpdrho, dpde, e1
!
! ASSIGN EOS PARAMETERS PER MATERIAL
!
c1    = mcm_mat(par(i)%mat)%eosinput(1)
s1    = mcm_mat(par(i)%mat)%eosinput(2)
s2    = mcm_mat(par(i)%mat)%eosinput(3)
s3    = mcm_mat(par(i)%mat)%eosinput(4)
gamma = mcm_mat(par(i)%mat)%eosinput(5)
a1    = mcm_mat(par(i)%mat)%eosinput(6)
! this EoS uses internal energy per unit initial volume
e1    = par(i)%e * par(i)%rho0/par(i)%mass
mu    = par(i)%rho / par(i)%rho0 - 1.0_d
!
! CHECK WHETHER MATERIAL IS IN TENSION/COMPRESSION
! 
IF(mu.gt.0.0_d) THEN
   !
   ! CALCULATE PRESSURE IN TENSION
   !
   termone = mu*(1.0_d-(gamma/2.0_d))
   termtwo = (a1/2.0_d)*mu*mu
   eosnum  = mu * par(i)%rho0 * (c1**2) * (1.0_d+termone-termtwo)
   mu11 = mu/(mu+1.0_d)
   mu21 = (mu*mu)/(mu+1.0_d)
   mu32 = (mu*mu*mu)/((mu+1.0_d)*(mu+1.0_d))
   eosden  = 1.0_d-(s1-1.0_d)*mu-s2*mu21-s3*mu32
   par(i)%p = eosnum/eosden**2 + (gamma+a1*mu)*e1
ELSE
   !
   ! CALCULATE p(i) IN COMPRESSION
   !
   par(i)%p = par(i)%rho0*c1**2*mu + (gamma+a1*mu)*e1
ENDIF
!
END SUBROUTINE mcm_eos4p