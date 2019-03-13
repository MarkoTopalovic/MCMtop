SUBROUTINE mcm_eqos(i)
!************************************************************************
!
!    Purpose: Update pressure and internal energy
!
!  Called by: 
!
!       Date: 6-1-99
!
!	  Modified : 16-08-05, J. Reveles to support EOS 1
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
!
SELECT CASE (mcm_mat(par(i)%mat)%eos)

CASE (1) !Polynomial
!
   ! CALCULATE PRESSURE AND INTERNAL ENERGY
    CALL mcm_eqos1(i)
!
   CASE (4) ! GRUNEISEN
	  !
	  ! CALCULATE PRESSURE AND INTERNAL ENERGY
	  !
	  CALL mcm_eqos4(i)
	!
   CASE (13) ! IDEAL GAS
	  !
	  ! CALCULATE PRESSURE AND INTERNAL ENERGY
	  !
	  CALL mcm_eqos13(i)
	!
   CASE (28) ! MURNAGHAN (QUASI-INCOMPRESSIBLE)
	  !
	  ! CALCULATE PRESSURE AND INTERNAL ENERGY
	  !
	  CALL mcm_eqos28(i)
	  !
   CASE (41) ! MIE - GRUNEISEN
	  !
	  ! CALCULATE PRESSURE, INTERNAL ENERGY speed of sound temperature
	  !
	  CALL mcm_eqos41(i)
		!
END SELECT

RETURN

END SUBROUTINE mcm_eqos