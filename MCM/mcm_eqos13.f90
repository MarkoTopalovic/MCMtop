SUBROUTINE mcm_eqos13(i)
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
REAL(kind=real_acc) :: gamma_m1, e1try, dvol, vol0,volold,volnew,vavg
!
! ASSIGN EOS PARAMETERS PER MATERIAL 
!
gamma_m1  = mcm_mat(par(i)%mat)%eosinput(1)-1.0_d
!
! specific 'trial' value for internal energy
! 
e1try  = par(i)%etry/par(i)%mass
dvol   = par(i)%mass/par(i)%rho - par(i)%mass/par(i)%rhoold
vol0   = par(i)%mass/par(i)%rho0
!
! CALCULATE p(i)
!         (according to J ANDERSON: Modern Compressible Flow
!           James' notes)
!
par(i)%p = (par(i)%rho*gamma_m1*e1try) / (1.0_d+0.5_d*par(i)%rho*gamma_m1*dvol/par(i)%mass)
!
! CALCULATE e(i)
!
par(i)%e = par(i)%etry - 0.5_d*dvol*par(i)%p
!
RETURN
!
END SUBROUTINE mcm_eqos13