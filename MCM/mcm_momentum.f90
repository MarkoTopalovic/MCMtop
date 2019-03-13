subroutine mcm_momentum
!************************************************************************
!
!    Purpose: control routine for calculation of acceleration
!
!  Called by: solution
!
!       Date: 09-08-2002
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
! If ghost particles are active update stress
!
if(mcm_boundary) call mcm_update_ghost_stress
!
select case (mcm_disctype)
 case(0)
 ! 
 ! Basic SPH
 !
 select case (mcm_axopt)
  case(4)
   ! 2D axisymmetric
   !call mcm_axiaccel
  case default
   ! 1D, 2D, 3D Cartesian
   call mcm_oldacceleration
  end select
  !
 case(1)
 !
 ! Mixed correction SPH
 call mcm_corr_acceleration
 !
end select
!
end subroutine mcm_momentum