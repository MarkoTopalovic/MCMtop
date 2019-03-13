SUBROUTINE mcm_bodyforce 
!************************************************************************
!
!    Purpose: 
!
!  Called by: 
!
!       Date: 
!
!     Errors: 
!
!      Notes: Currently this routine is not an option accessable through
!             the input file. Call is commented out. This routine 
!             hardwires the acceleration due to gravity to be 9.81 in the
!             negative y direction.
!
!************************************************************************
!
USE mcm_database
!
IMPLICIT NONE
!
INTEGER :: i, mgi
REAL    :: gravity
! hardwired gravity in m/s^2
gravity = -9.81_d !E-11
!
!F (numdp.ne.0) THEN
   !
!  IF (nmigs.ne.0) THEN
      !
!  DO i = 1,np
	     !
!            DO mgi=1,nmigs
            !
!           IF (mtigs(mgi).eq.pmat(i)) THEN
	           !
!              a(idirgv,i) = a(idirgv,i) + gravity
            !
!	        ENDIF
            !
!         ENDDO
		 !
!	  ENDDO
      !
!   ELSE
      !
      DO i = 1,mcm_np
	     !
         par(i)%a(2) = par(i)%a(2) + gravity
		 !
	  ENDDO
      !
!   ENDIF
   !
!ENDIF
!
RETURN
!
END SUBROUTINE mcm_bodyforce