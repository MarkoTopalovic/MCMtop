SUBROUTINE mcm_sueos(i)
!************************************************************************
!
!    Purpose: Calculate speed of sound for particle i
!
!  Called by: 
!
!       Date: 09-08-2002
!
!       
!	  Modified : 16-08-05, J. Reveles to support EOS 1
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
INTEGER :: i
!
SELECT CASE (mcm_mat(par(i)%mat)%eos)

CASE (1)
 ! Polynomial
 CALL mcm_eos1c(i)
 !
 CASE (4)
  ! Dyna Gruneisen
  CALL mcm_eos4c(i)
  !
 CASE (13)
  ! Perfect gas
  CALL mcm_eos13c(i)
  !
 CASE (28)
  ! Murnaghan quasi-incompressible fluid
  CALL mcm_eos28c(i)
  !
 CASE (41)
  ! Mie-Gruneisen
  CALL mcm_eos41c(i)
  !
END SELECT
!
END SUBROUTINE mcm_sueos
      