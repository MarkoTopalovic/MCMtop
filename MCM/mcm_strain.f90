subroutine mcm_strain
!************************************************************************
!
!    Purpose: control routine for calculation of rate of deformation tensor
!
!  Called by: solution
!
!       Date: 09-08-2002
!
! Last Modified: 03-03-2005 by J. Campbell
!
!     Errors: 
!
!      Notes: calls correct subroutine for interpolation method used
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,j
!
! 
! Half step positions back in time to that poitions and velocities are held at the same time
!
DO i=1,mcm_np
 DO j = 1,mcm_ndim
  par(i)%x(j) = par(i)%x(j) - 0.5_d * par(i)%v(j) * mcm_dt 
 ENDDO 
ENDDO   
!              
!
select case (mcm_disctype)
 !
 case(0)
 ! 
 ! Basic SPH
 !
 select case (mcm_axopt)
  case(4)
   ! 2D axisymmetric
   !call mcm_axistrain
  case default
   ! 1D, 2D, 3D Cartesian
   call mcm_oldstrain
  end select
  !
 case(1)
 !
 ! Mixed correction SPH
 call mcm_corr_strain
 !
end select
!
DO i=1,mcm_np
 DO j = 1,mcm_ndim
  par(i)%x(j) = par(i)%x(j) + 0.5_d * par(i)%v(j) * mcm_dt 
 ENDDO 
ENDDO   
!
end subroutine mcm_strain