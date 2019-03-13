subroutine mcm_initrod
!************************************************************************
!
!		Purpose: Initialise Rate of Deformation Tensor
!
!	  Called by: Initial
!
!	     Author: Tom De Vuyst
!
!          Date: 05-08-2002
!
! Last Modified: 05-08-2002 by J. Campbell
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
INTEGER :: i, j, k
!
DO i = mcm_ssp,mcm_esp
 DO j = 1,3
  DO k = 1,3
   par(i)%rod(k,j) = 0.0_d
   par(i)%spin(k,j) = 0.0_d
  ENDDO
 ENDDO
ENDDO
!
END SUBROUTINE mcm_initrod   
