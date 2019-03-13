 SUBROUTINE mcm_calcg
!************************************************************************
!
!		Purpose: Calculate or assign shear modulus for each material
!
!	  Called by: initial
!
!	     Author: Tom De Vuyst
!
!          Date: 02-08-2002
!
! Last Modified: 02-08-2002 by J. Campbell
!
!        Errors: specified material model number does not exist 
!
!         Notes: the shear modulus is needed to calculate the speed of sound
!             for the time step calculation. If G is not specified in the
!             input (J-C), or mu/2 (fluid) it is calculated using
!                               E
!                         G= -------       where v=poissons ratio
!                            2(1+v)
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i
!
DO i = 1,mcm_nummat
 !
 SELECT CASE (mcm_mat(i)%model)
   !
   CASE (1,3)  ! ELASTIC, ELASTO-PLASTIC
      !
      mcm_mat(i)%g = mcm_mat(i)%strinput(1) / (2.0_d* (1.0_d+mcm_mat(i)%strinput(6)) )
      !
      CASE (4)  ! ELASTIC, ELASTO-PLASTIC
      !
      mcm_mat(i)%g = mcm_mat(i)%strinput(6)
      !
   CASE (9)    ! FLUID
      !
      mcm_mat(i)%g = 0.0_d
      !
   CASE(10) ! HYDRO-DYNAMIC
      !
      mcm_mat(i)%g = mcm_mat(i)%strinput(1)
	  !
      !
   CASE DEFAULT
      !
      WRITE(*,1000)
      WRITE(*,1010) mcm_mat(i)%model
      !
 END SELECT
 !
ENDDO
!
RETURN
!
1000 FORMAT(//5X, 'Error in subroutine initial')
1010 FORMAT(/5X,' Material number ',I3,' does not exist.')
!
END SUBROUTINE mcm_calcg