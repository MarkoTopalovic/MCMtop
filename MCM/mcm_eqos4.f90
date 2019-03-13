SUBROUTINE mcm_eqos4(i)
!************************************************************************
!
!    Purpose: Calculate pressure and internal energy based on trial value 
!             of internal energy 
!             Uses Gruneisen Eos, identical to that used in  DYNA3D
!
!  Called by: epupdate
!
!       Date: 6-1-99
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
INTEGER :: i
REAL(kind=real_acc) :: mu,s1,s2,s3,c1,gamma,a1
REAL(kind=real_acc) :: mu11, mu21, mu22, mu32, mu33
REAL(kind=real_acc) :: termone, termtwo, eosnum, eosden, termthree, termfour
REAL(kind=real_acc) :: dpdrho, dpde, e1try, dvol, vol0
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
! specific 'trial' value of int. energy
!
mu    = par(i)%rho/par(i)%rho0 - 1.0
dvol  = par(i)%mass/par(i)%rho - par(i)%mass/par(i)%rhoold
vol0  = par(i)%mass/par(i)%rho0 
e1try  = par(i)%etry/vol0
!
! CHECK WHETHER MATERIAL IS IN TENSION/COMPRESSION
! 
IF(mu.gt.0.0) THEN
   !
   ! CALCULATE p(i) IN COMPRESSION
   !
   termone = mu*(1.0-(gamma/2.0))
   termtwo = (a1/2.0)*mu*mu
   eosnum  = mu*par(i)%rho0*(c1*c1)*(1.0+termone-termtwo)
   mu11    = mu/(mu+1.0)
   mu21    = (mu*mu)/(mu+1.0)
   mu32    = (mu*mu*mu)/((mu+1.0)*(mu+1.0))
   eosden  = 1.0-(s1-1.0)*mu-s2*mu21-s3*mu32
   !
   !	DYNA current gama is = gama+a1*mu
   par(i)%p    = (eosnum/(eosden*eosden) + (gamma+a1*mu)*e1try) / &
             & (1.0 + 0.5*(gamma+a1*mu)*dvol/vol0)
   !
ELSE
   !
   ! CALCULATE p(i) IN TENSION
   !
   par(i)%p    = (par(i)%rho0*(c1*c1)*mu + (gamma+a1*mu)*e1try) / &
             & (1.0 + 0.5*(gamma+a1*mu)*dvol/vol0)
   !
ENDIF
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

END SUBROUTINE mcm_eqos4