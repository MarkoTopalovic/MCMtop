SUBROUTINE mcm_initsigma
!************************************************************************
!
!		Purpose: Initialise the stress tensor
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
INTEGER :: i, j , k
!
DO i = mcm_ssp,mcm_esp
 !
 par(i)%epx(1) = 0.0
 par(i)%epx(2) = 0.0
 par(i)%epx(3) = 0.0
 par(i)%epx(4) = 0.0
 par(i)%epx(5) = 0.0
 par(i)%epx(6) = 0.0
 !
 DO j = 1,3
  DO k = 1,3
   IF (j.eq.k) THEN
    par(i)%sigma(k,j) = -par(i)%p
   ELSE
    par(i)%sigma(k,j) = 0.0_d
   ENDIF
  ENDDO
 ENDDO
 !
 ! Initialise deviatoric stress and back stress
 !
 par(i)%s = 0.0_d
 par(i)%alfa = 0.0_d
ENDDO
!
END SUBROUTINE mcm_initsigma