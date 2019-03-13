subroutine mcm_calcrhopc
!************************************************************************
!
!		Purpose: Calculate rho based on initial relative volume if 
!                applicable and assign rho to every particle
!                Calculate pressure and speed of sound
!
!	  Called by: initial
!
!	     Author: Tom De Vuyst
!
!          Date: 02-08-2002
!
! Last Modified: 02-08-2002 by J. Campbell
!
!        Errors: Wrong material model no. or EOS no.
!
!         Notes:
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i,j
REAL(kind=real_acc) :: pmax
!
! LOOP OVER ALL PARTICLES
!
DO i = mcm_ssp, mcm_esp
   !
   ! CALCULATE rho0 
   !
   par(i)%rho0 = mcm_mat(par(i)%mat)%rho
   !
   ! CALCULATE rho DEPENDING INITIAL REL. VOLUME IN EOS
   ! CASE SELECTS MATERIAL TYPES
   !
   SELECT CASE (mcm_mat(par(i)%mat)%model)
      !
      CASE (1,3)
	     !
		 ! Elastic and elastic-plastic
         !
		 if(mcm_init_rhoe.ne.1) then
		 par(i)%rho  = par(i)%rho0
		 par(i)%e    = 0.0_d
         endif
		 par(i)%p    = 0.0_d
		 par(i)%c    = sqrt(mcm_mat(par(i)%mat)%strinput(2)/par(i)%rho)
		 !par(i)%c    = sqrt((7*3000)/par(i)%rho)
		 !par(i)%c    = 10*sqrt(2*0.00981*mcm_mat(par(i)%mat)%strinput(5))
         !
         CASE (4)
	     !
		 ! Elastic and elastic-plastic
         !
		 if(mcm_init_rhoe.ne.1) then
		 par(i)%rho  = par(i)%rho0
		 par(i)%e    = 0.0_d
         endif
		 par(i)%p    = 0.0_d
		 par(i)%c    = sqrt(mcm_mat(par(i)%mat)%strinput(10)/par(i)%rho)
		  !par(i)%c    = sqrt((7*3000)/par(i)%rho)
		 !par(i)%c    = 10*sqrt(2*0.00981*mcm_mat(par(i)%mat)%strinput(5))
         !
      CASE (9)
         !
         ! fluid
         !
         CALL mcm_eoscalc(i)
		 !
		 ! PRESSURE CUT-OFF
		 !
		 pmax = mcm_mat(par(i)%mat)%strinput(1)
		 IF (par(i)%p.lt.pmax) THEN
		    par(i)%p = pmax
		 ENDIF
		 !
      CASE (10)
         !
         ! Elastic-plastic-hydrodynamic 
         !
         CALL mcm_eoscalc(i)
         !
		 !
      CASE DEFAULT
	     !
         WRITE(*,1000)
         WRITE(*,1010) par(i)%mat
         !
   END SELECT
   !
ENDDO
!
RETURN
!
1000 FORMAT(//5X, 'Error in subroutine calcrho')
1010 FORMAT(/5X,' Material number ',I3,' does not exist.')
!
END SUBROUTINE mcm_calcrhopc