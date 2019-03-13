SUBROUTINE mcm_eos28c(i)
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
REAL(kind=real_acc) :: gamma, e0, dpdrho, dpde, e1, patmos, b
!
! ASSIGN EOS PARAMETERS PER MATERIAL 
!
b      = mcm_mat(par(i)%mat)%eosinput(1)
gamma  = mcm_mat(par(i)%mat)%eosinput(2)
!
! CALCULATE c(i) 
!
par(i)%c   = sqrt(b*gamma/par(i)%rho0)    

RETURN

END SUBROUTINE mcm_eos28c