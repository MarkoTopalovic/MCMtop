SUBROUTINE mcm_eos1p(i)
!************************************************************************
!
!		Purpose: Calculate pressure, dp/drho, dp/dE, c FOR Polynomial EOS
!
!	  Called by: Eoscalc
!
!	     Author: J.Reveles
!
!          Date: 12-08-2005
!
! Last Modified: 
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
INTEGER :: i,j
REAL(kind=real_acc) :: mu,c0,c1,c2,c3,c4,c5,c6,mu_t,e1

!
! ASSIGN EOS PARAMETERS PER MATERIAL
!
c0    = mcm_mat(par(i)%mat)%eosinput(1)
c1    = mcm_mat(par(i)%mat)%eosinput(2)
c2    = mcm_mat(par(i)%mat)%eosinput(3)
c3    = mcm_mat(par(i)%mat)%eosinput(4)
c4    = mcm_mat(par(i)%mat)%eosinput(5)
c5    = mcm_mat(par(i)%mat)%eosinput(6)
c6    = mcm_mat(par(i)%mat)%eosinput(7)
!
mu    = par(i)%rho / par(i)%rho0 - 1.0_d
mu_t= max(mu,0.0_d)
e1    = par(i)%e * par(i)%rho0/par(i)%mass

!
   ! CALCULATE PRESSURE 
 
e1=0.0_d  !just to simulate ab problem
!
   par(i)%p = c0+(c1+c3*mu**2)*mu+c2*mu_t**2+(c4+c5*mu+c6*mu_t**2)*e1
!
END SUBROUTINE mcm_eos1p
