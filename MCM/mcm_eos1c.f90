SUBROUTINE mcm_eos1c(i)
!************************************************************************
!
!		Purpose: Calculate pressure, dp/drho, dp/dE, c FOR Polynomial EOS
!
!	  Called by: Eoscalc(i)
!
!	     Author: J.Reveles
!
!          Date: 16-08-2005
!
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
REAL(kind=real_acc) :: c0,c1,c2,c3,c4,c5,c6
REAL(kind=real_acc) :: dpdrho, dpde, e1,dmudrho
REAL(kind=real_acc) :: termone, termtwo
REAL(kind=real_acc) :: c_sq
real(kind=real_acc) :: mu,mu_t
!
! ASSIGN EOS PARAMETERS PER MATERIAL
!
c0    = mcm_mat(par(i)%mat)%eosinput(1)
c1    = mcm_mat(par(i)%mat)%eosinput(2)
c2    = mcm_mat(par(i)%mat)%eosinput(3)
c3    = mcm_mat(par(i)%mat)%eosinput(4)
c4    = mcm_mat(par(i)%mat)%eosinput(5)
c5    = mcm_mat(par(i)%mat)%eosinput(6)
c6    = mcm_mat(par(i)%mat)%eosinput(7)
!
mu    = par(i)%rho / par(i)%rho0 - 1.0_d
mu_t= max(mu,0.0_d)
!
! CONVERT TO SPECIFIC INTERNAL ENERGY
!
! this EoS uses internal energy per unit initial volume
e1    = par(i)%e * par(i)%rho0/par(i)%mass
!
e1=0.0_d  !just for ab simulation
!
! CALCULATE dpdrho

dmudrho=1/par(i)%rho0
!
!
termone=c1*dmudrho +(c2*(dmudrho*dmudrho))+(3.0_d*c3*dmudrho)*mu*mu
!
termtwo=(c5*dmudrho+c6*dmudrho*dmudrho)*e1
!
   dpdrho  = termone+termtwo
!
! CALCULATE dpde
!
dpde   = c4+c5*mu+c6*(mu_t*mu_t)
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
1000 format(//5x,'Error in subroutine eos1c',/ &
              5x,'  Complex speed of sound, calculation forced to continue.')
1100 format(/5x,' Particle: ',i6,/5x,'c_sq: 'e14.6,/5x,'g(pmat(i): ',e14.6,&
            /5x,'dpdrho: ',e14.6,/5x,'p(i): ',e14.6,/5x,'rho(i): ',e14.6,&
			/5x,'dpde: ',e14.6)
!
END SUBROUTINE mcm_eos1c

