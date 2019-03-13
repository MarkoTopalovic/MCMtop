SUBROUTINE mcm_eos4c(i)
!************************************************************************
!
!		Purpose: Calculate pressure, dp/drho, dp/dE, c FOR GRUNEISEN EOS
!
!	  Called by: Eoscalc(i)
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
INTEGER :: i
REAL(kind=real_acc) :: mu,s1,s2,s3,c1,gamma,a1
REAL(kind=real_acc) :: mu11, mu21, mu22, mu32, mu33
REAL(kind=real_acc) :: termone, termtwo, eosnum, eosden, termthree, termfour
REAL(kind=real_acc) :: dpdrho, dpde, e1
REAL(kind=real_acc) :: c_sq
!
! ASSIGN EOS PARAMETERS PER MATERIAL
!
c1    = mcm_mat(par(i)%mat)%eosinput(1)
s1    = mcm_mat(par(i)%mat)%eosinput(2)
s2    = mcm_mat(par(i)%mat)%eosinput(3)
s3    = mcm_mat(par(i)%mat)%eosinput(4)
gamma = mcm_mat(par(i)%mat)%eosinput(5)
a1    = mcm_mat(par(i)%mat)%eosinput(6)
!
! CONVERT TO SPECIFIC INTERNAL ENERGY
!
! this EoS uses internal energy per unit initial volume
e1    = par(i)%e * par(i)%rho0/par(i)%mass
mu    = par(i)%rho / par(i)%rho0 - 1.0_d
!
! CHECK WHETHER MATERIAL IS IN TENSION/COMPRESSION
! 
IF(mu.gt.0.0_d) THEN
   !
   ! CALCULATE dpdrho IN TENSION
   !
   termone = mu*(1.0_d-(gamma/2.0_d))
   termtwo = (a1/2.0_d)*mu*mu
   eosnum  = mu*par(i)%rho0*(c1**2)*(1.0_d+termone-termtwo)
   mu11 = mu/(mu+1.0_d)
   mu21 = (mu*mu)/(mu+1.0_d)
   mu22 = (mu*mu)/((mu+1.0_d)*(mu+1.0_d))
   mu32 = (mu*mu*mu)/((mu+1.0_d)*(mu+1.0_d))
   mu33 = (mu*mu*mu)/((mu+1.0_d)*(mu+1.0_d)*(mu+1.0_d))
   eosden  = 1.0_d-(s1-1.0_d)*mu-s2*mu21-s3*mu32
   termthree = par(i)%rho0*(c1**2)*(1.0_d+2.0_d*termone-3.0*termtwo)
   termfour  = 1.0_d-s1-2.0_d*s2*mu11+s2*mu22-3.0_d*s3*mu22+2.0_d*s3*mu33
   dpdrho  = (termthree/eosden**2-2.0_d*eosnum*termfour/eosden**3)/par(i)%rho0 &
            & + a1*e1/par(i)%rho0
ELSE
   dpdrho  = c1**2+a1*e1/par(i)%rho0
ENDIF
!
! CALCULATE dpde
!
dpde   = gamma+a1*mu
!
! CALCULATE c(i)
! **********************************
c_sq = ( (4.0_d/3.0_d)*mcm_mat(par(i)%mat)%g/par(i)%rho0 + dpdrho + &
         par(i)%p/par(i)%rho0*(par(i)%mass/par(i)%rho)**2*dpde)
if(c_sq.gt.0.0_d) then
 par(i)%c = sqrt(c_sq)
else
 write(*,1000)
 write(13,1000)
 write(*,1100) i,c_sq,mcm_mat(par(i)%mat)%g,dpdrho,par(i)%p,par(i)%rho,dpde
 write(13,1100) i,c_sq,mcm_mat(par(i)%mat)%g,dpdrho,par(i)%p,par(i)%rho,dpde
!
! call mcm_shutdown(2)
!
 par(i)%c = sqrt(abs(c_sq))		! force calculation to continue
!
endif
!
!
!
RETURN
!
1000 format(//5x,'Error in subroutine eos4c',/ &
              5x,'  Complex speed of sound, calculation forced to continue.')
1100 format(/5x,' Particle: ',i6,/5x,'c_sq: 'e14.6,/5x,'g(pmat(i): ',e14.6,&
            /5x,'dpdrho: ',e14.6,/5x,'p(i): ',e14.6,/5x,'rho(i): ',e14.6,&
			/5x,'dpde: ',e14.6)
!
END SUBROUTINE mcm_eos4c
