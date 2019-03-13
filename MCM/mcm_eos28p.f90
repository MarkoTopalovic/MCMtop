SUBROUTINE mcm_eos28p(i)
!************************************************************************
!
!		Purpose: CALCULATE p(i) AND c(i) FOR MURNAGHAN EOS
!
!	  Called by: eoscalc
!
!	     Author: Tom De Vuyst
!
!          Date: 13-5-02
!
! Last Modified: 13-5-02 by Tom De Vuyst
!
!        Errors: 
!
!         Notes:
!
!************************************************************************
USE mcm_database
!
IMPLICIT NONE
!
INTEGER :: i
REAL(kind=real_acc) :: gamma, e1, b
!
! ASSIGN EOS PARAMETERS PER MATERIAL 
!
B      = mcm_mat(par(i)%mat)%eosinput(1)
gamma  = mcm_mat(par(i)%mat)%eosinput(2)
!
! CALCULATE p(i)
!
par(i)%p   = B*((par(i)%rho/par(i)%rho0)**gamma-1.0)
  

RETURN

END SUBROUTINE mcm_eos28p