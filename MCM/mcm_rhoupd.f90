SUBROUTINE mcm_rhoupd
!************************************************************************
!
!    Purpose: Update the density for Eulerian kernels
!
!  Called by: Constitutive
!
!       Date: 7-1-99
!
! Last Modified: 03-03-2005 by J. Campbell
!
!     Errors: 
!
!      Notes:
!
!************************************************************************
USE mcm_database
!
IMPLICIT NONE

INTEGER :: i
REAL(kind=real_acc) :: rhomax, rhomin, othird
!
do i=mcm_ssp,mcm_esp
 par(i)%rhoold=par(i)%rho
 !
 rhomax = +1.0e+04_d
 rhomin = -1.0e+04_d
 othird = 1.0_d / 3.0_d
 !
 par(i)%rho = par(i)%rho * ( 1. - par(i)%tracerod * mcm_dt )
 !
 !	James introduced this to limit the density
 !
 if(par(i)%rho.gt.mcm_mat(par(i)%mat)%rho_max) then
   !
   par(i)%rho = par(i)%rho / ( 1. - par(i)%tracerod * mcm_dt )
   par(i)%rod(1,1) = par(i)%rod(1,1) - othird*par(i)%tracerod
   par(i)%rod(2,2) = par(i)%rod(2,2) - othird*par(i)%tracerod
   par(i)%rod(3,3) = par(i)%rod(3,3) - othird*par(i)%tracerod
   par(i)%tracerod = 0.0_d
   par(i)%h = mcm_mat(par(i)%mat)%h * ( mcm_mat(par(i)%mat)%rho / par(i)%rho ) ** ( 1.0_d/real(mcm_ndim))
   !
 endif
 !
 if(par(i)%rho.gt.mcm_mat(par(i)%mat)%rho_max) then
   !
   par(i)%rho = par(i)%rho / ( 1. - par(i)%tracerod * mcm_dt )
   par(i)%rod(1,1) = par(i)%rod(1,1) - othird*par(i)%tracerod
   par(i)%rod(2,2) = par(i)%rod(2,2) - othird*par(i)%tracerod
   par(i)%rod(3,3) = par(i)%rod(3,3) - othird*par(i)%tracerod
   par(i)%tracerod = 0.0_d
   par(i)%h = par(i)%hold
   par(i)%h = mcm_mat(par(i)%mat)%h * ( mcm_mat(par(i)%mat)%rho / par(i)%rho ) ** ( 1.0_d/real(mcm_ndim))
   !
 endif
 !
 IF ( par(i)%rho .lt. rhomin ) THEN
  WRITE(*,1000)
  WRITE(*,1010) rhomin, i
  WRITE(13,1000)
  WRITE(13,1010) rhomin, i
 ENDIF
 IF ( par(i)%rho .gt. rhomax ) THEN
  WRITE(*,1000)
  WRITE(*,1020) rhomax, i
  WRITE(13,1000)
  WRITE(13,1020) rhomax, i
 ENDIF
 !
enddo
!
1000 format(//5x,'WARNING in subroutine rhoupd')
1010 format(/5x,'rhomin = ',E10.4,' for particle ',I7)
1020 format(/5x,'rhomax = ',E10.4,' for particle ',I7)
!
END SUBROUTINE mcm_rhoupd