SUBROUTINE mcm_eqos1(i)
!************************************************************************
!
!    Purpose: Calculate pressure and internal energy based on trial value 
!             of internal energy 
!             Uses Polynomial Eos
!
!  Called by: epupdate
!
!       Date: 16-08-2005
!
!     Errors: 
!
!      Notes:
!
!************************************************************************
!
USE mcm_database
!
IMPLICIT NONE
!
REAL(kind=real_acc) :: mu11, mu21, mu22, mu32, mu33
REAL(kind=real_acc) :: termone, termtwo, eosnum, eosden
REAL(kind=real_acc) :: dpdrho, dpde, e1try, dvol, vol0
INTEGER :: i,j
REAL(kind=real_acc) :: mu,c0,c1,c2,c3,c4,c5,c6,mu_t
REAL(kind=real_acc) :: e1,dmudrho
REAL(kind=real_acc) :: c_sq
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
! specific 'trial' value of int. energy
!
dvol  = par(i)%mass/par(i)%rho - par(i)%mass/par(i)%rhoold
vol0  = par(i)%mass/par(i)%rho0 
e1try  = par(i)%etry/vol0
!
   !
   ! CALCULATE p(i) 
   !
termone=c0+(c1+c3*mu**2)*mu+c2*mu_t**2
termtwo=c4+c5*mu+c6*mu_t**2

      par(i)%p = (termone+termtwo*e1try)/(1+termtwo*(dvol/vol0))
!
! Pressure cutoff
!
if(par(i)%p.lt.par(i)%pcut) then
   !
   par(i)%p   = par(i)%pcut
   par(i)%rho = par(i)%rho / ( 1. - par(i)%tracerod * mcm_dt )
   !
else
   !
   par(i)%e = par(i)%etry - 0.5*dvol*par(i)%p
   !
endif
!
!if(fail1(i).eq.0.0) p(i) = max(p(i),0.0)
!
! CALCULATE e(i)
!
RETURN

END SUBROUTINE mcm_eqos1
