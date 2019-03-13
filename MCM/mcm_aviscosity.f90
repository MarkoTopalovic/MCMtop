SUBROUTINE mcm_aviscosity(i)
!************************************************************************
!
!		Purpose: Calculate artificial viscosity tensor 
!
!	  Called by: Initial
!
!	     Author: Tom De Vuyst
!
!          Date: 05-08-2002
!
! Last Modified: 05-08-2002
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
!
! SELECT TYPE OF ARTIFICIAL VISCOSITY
!
SELECT CASE (mcm_mat(par(i)%mat)%visc_type)
   !
   ! NODAL ARTIFICIAL VISCOSITY (G. JOHNSON 96c p2728)
   !
   CASE (1)
      CALL mcm_nodalav(i)
      !
   CASE DEFAULT
      WRITE(*,1000)
	  WRITE(*,1010) mcm_mat(par(i)%mat)%visc_type
      !
END SELECT

1000 FORMAT(//5X, 'Error in subroutine aviscosisty')
1010 FORMAT(/5X,'Artificial viscosity number ',I3,' does not exist.')
!
END SUBROUTINE mcm_aviscosity
