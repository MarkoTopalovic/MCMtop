SUBROUTINE mcm_nodalav(i)
!************************************************************************
!
!		Purpose: Calculate  Nodal Art. Visc. (Pressure viscosity)
!
!	  Called by: aviscosity
!
!	     Author: Tom De Vuyst
!
!          Date: 05-08-2002
!
! Last Modified: 
!
!        Errors: 
!
!         Notes: Based on eq (14) in G. Johnsons 'Normalised Smoothing 
!                Functions for SPH Impact Computations' Int. J. Num. 
!                Meth. Engng. Vol. 39, 2725-2741 (1996)
!
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i, j, k, l, mat
real(kind=real_acc) :: q
!
! CALCULATE  ARTIFICIAL VISCOSITY
!
mat = par(i)%mat
!
DO j = 1,3
 DO k = 1,3
  par(i)%q(j,k) = 0.0_d
 ENDDO
ENDDO
!
if(par(i)%tracerod.lt.0.0_d) then
 q = mcm_mat(mat)%av_l * par(i)%rho * par(i)%c * par(i)%h * abs(par(i)%tracerod) + & ! linear
     mcm_mat(mat)%av_q * par(i)%rho * (par(i)%h * par(i)%tracerod)**2                ! quadratic
 !
 do j=1,3
  par(i)%q(j,j) = q
 enddo
endif
!
END SUBROUTINE mcm_nodalav