SUBROUTINE mcm_initold
!************************************************************************
!
!    Purpose: Initialise rhoold and qold
!
!  Called by: Initial
!
!       Date: 05-08-2002
!
!     Errors: 
!
!      Notes:
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i,j,k
!
DO i = mcm_ssp,mcm_esp
 par(i)%rhoold = par(i)%rho
 par(i)%qold(1:3,1:3) = par(i)%q(1:3,1:3)
 !
ENDDO
!
END SUBROUTINE mcm_initold
