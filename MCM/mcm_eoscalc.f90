SUBROUTINE mcm_eoscalc(i)
!************************************************************************
!
!		Purpose: Select EOS model to calculate pressure, dp/drho, dp/dE
!
!	  Called by: Initial
!
!	     Author: Tom De Vuyst
!
!          Date: 24-11-98
!
! Last Modified: 24-11-98 by Tom De Vuyst
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
integer :: pmat
!
pmat = par(i)%mat
!
SELECT CASE (mcm_mat(pmat)%eos)
 !
 ! CASE FOR EOS NUMBERS
 !
Case(1)  ! polynomial EOS
   if(mcm_init_rhoe.ne.1) then
    par(i)%rho = par(i)%rho0 / mcm_mat(pmat)%eosinput(9)
    par(i)%e   = mcm_mat(pmat)%eosinput(8) * par(i)%mass / par(i)%rho
   endif

!CALCULATES p AND c 
   !
   CALL mcm_eos1p(i)
   CALL mcm_eos1c(i)
!

 CASE (4) ! GRUNEISEN
   !
   ! CALCULATES DENSITY
   !
   if(mcm_init_rhoe.ne.1) then
   par(i)%rho = par(i)%rho0 / mcm_mat(pmat)%eosinput(8)
   par(i)%e   = mcm_mat(pmat)%eosinput(7) * par(i)%mass / par(i)%rho
   endif
   !
   ! CALCULATES p AND c 
   !
   CALL mcm_eos4p(i)
   CALL mcm_eos4c(i)
   !
 CASE(13) ! IDEAL GAS
   !
   ! CALCULATES DENSITY
   !
   if(mcm_init_rhoe.ne.1) then
   par(i)%rho = par(i)%rho0 / mcm_mat(pmat)%eosinput(4)
   par(i)%e   = mcm_mat(pmat)%eosinput(3) * par(i)%mass / par(i)%rho
   endif
   !
   ! CALCULATES p AND c 
   !
   CALL mcm_eos13p(i)
   CALL mcm_eos13c(i)
   !
   !
 CASE(28) ! MURNAGHAN INCOMPRESSIBLE FLUID
   !
   ! CALCULATES DENSITY rho
   !
   if(mcm_init_rhoe.ne.1) then
   par(i)%rho = par(i)%rho0/mcm_mat(pmat)%eosinput(4)
   par(i)%e   = mcm_mat(pmat)%eosinput(3)*par(i)%mass
   endif
   !
   ! CALCULATES p AND c 
   !
   CALL mcm_eos28p(i)
   CALL mcm_eos28c(i)
   !
 CASE (41) ! MIE - GRUNEISEN
   !
   ! CALCULATES DENSITY
   !
   if(mcm_init_rhoe.ne.1) then
   par(i)%rho = par(i)%rho0 / mcm_mat(pmat)%eosinput(6)
   par(i)%e   = mcm_mat(pmat)%eosinput(3)*mcm_mat(pmat)%eosinput(4) * & 
                  par(i)%mass / par(i)%rho
   endif
   !
   ! CALCULATES p AND c 
   !
   CALL mcm_eos41p(i)
   CALL mcm_eos41c(i)
   !
 CASE DEFAULT
   !
   WRITE(*,1000)
   WRITE(*,1010) mcm_mat(pmat)%eos
   !
   WRITE(13,1000)
   WRITE(13,1010) mcm_mat(pmat)%eos
   !
   call mcm_shutdown(2)
   !
END SELECT
!
1000 FORMAT(//5X,'Error in subroutine eoscalc')
1010 FORMAT(/5X,'EOS number ',I3,' does not exist.')
!
END SUBROUTINE mcm_eoscalc
