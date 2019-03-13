SUBROUTINE mcm_eos13p(i)
!************************************************************************
!
!		Purpose: CALCULATE p(i) AND c(i) FOR IDEAL GAS EOS
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
! CALCULATE p(i)
!
par(i)%p   = par(i)%rho * (gamma-1.0_d) * e1
!
END SUBROUTINE mcm_eos13p