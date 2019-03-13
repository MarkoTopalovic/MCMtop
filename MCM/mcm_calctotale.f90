SUBROUTINE mcm_calctotale
!************************************************************************
!
!		Purpose: Calculate total energy in the system
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
!         Notes: einternal = total internal energy in the problem
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i,j,mat
REAL(kind=real_acc)    :: einternal
!
do i=1,mcm_nummat
 mcm_mat(i)%ie = 0.0_d
 mcm_mat(i)%ke = 0.0_d
 mcm_mat(i)%mom = 0.0_d
enddo

do i = mcm_svp,mcm_evp
 ! kinetic energy
 if(par(i)%active) then
 mat = par(i)%mat
 mcm_mat(mat)%ke = mcm_mat(mat)%ke + 0.5_d * par(i)%mass * par(i)%vabs**2
  !
  do j=1,3
   mcm_mat(mat)%mom(j) = mcm_mat(mat)%mom(j) + par(i)%mass * abs(par(i)%v(j))
  enddo
 endif
enddo
!
do i = mcm_ssp, mcm_esp
 ! internal energy
 if(par(i)%active) then
 mat = par(i)%mat
 mcm_mat(mat)%ie = mcm_mat(mat)%ie + par(i)%e
 endif
enddo
!
mcm_kinetice = 0.0_d
mcm_internale = 0.0_d
do i=1,mcm_nummat
 mcm_kinetice  = mcm_kinetice  + mcm_mat(i)%ke
 mcm_internale = mcm_internale + mcm_mat(i)%ie
enddo
!
mcm_totale = mcm_kinetice + mcm_internale
!
end subroutine mcm_calctotale