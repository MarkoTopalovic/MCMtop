SUBROUTINE mcm_eqos28(i)
!************************************************************************
!
!		Purpose: CALCULATE p(i) AND c(i) FOR IDEAL GAS EOS
!
!	  Called by: epupdate
!
!	     Author: Tom De Vuyst
!
!          Date: 7-1-99
!
! Last Modified: 7-1-99 by Tom De Vuyst
!
!        Errors: 
!
!         Notes: *** Ideal Gas equation of state uses:
!                    INTERNAL ENERGY PER UNIT MASS
!
!************************************************************************
!
USE mcm_database
!
IMPLICIT NONE
!
INTEGER :: i
REAL(kind=real_acc) :: gamma, e1try, dvol, vol0,volold,volnew,vavg, b
!
! ASSIGN EOS PARAMETERS PER MATERIAL 
!
b        = mcm_mat(par(i)%mat)%eosinput(1)
gamma    = mcm_mat(par(i)%mat)%eosinput(2)
!
! specific 'trial' value for internal energy
! 
e1try    = par(i)%etry/par(i)%mass
dvol     = par(i)%mass/par(i)%rho - par(i)%mass/par(i)%rhoold
vol0     = par(i)%mass/par(i)%rho0
!
! CALCULATE p(i)
!
par(i)%p = b*((par(i)%rho/par(i)%rho0)**gamma - 1.0)
!
! CALCULATE e(i)
!
par(i)%e = par(i)%etry - 0.5*dvol*par(i)%p
!
RETURN
!
END SUBROUTINE mcm_eqos28