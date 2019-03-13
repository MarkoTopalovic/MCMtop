SUBROUTINE mcm_eos13c(i)
!************************************************************************
!
!		Purpose: CALCULATE SPEED OF SOUND FOR IDEAL GAS EOS
!
!	  Called by: eoscalc
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
REAL(kind=real_acc) :: gamma, e1
!
! ASSIGN EOS PARAMETERS PER MATERIAL 
!
gamma  = mcm_mat(par(i)%mat)%eosinput(1)
! EoS requires internal energy per unit mass
e1     = par(i)%e/par(i)%mass
!
! CALCULATE c
!

!par(i)%c = sqrt( gamma*par(i)%p/par(i)%rho )
par(i)%c = sqrt( (gamma-1.0_d)*gamma*e1 )
!
END SUBROUTINE mcm_eos13c