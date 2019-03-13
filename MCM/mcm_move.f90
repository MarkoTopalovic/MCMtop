subroutine mcm_move
!************************************************************************
!
!    Purpose: update particle velocity and position
!
!  Called by: solution
!
!       Date: 09-08-2002
!
!     Errors: 
!
!      Notes: Originally this routine updated both the velocity and position
!             The velocity update has been moved to a separate routine
!             so that the velocity can be modified/used before the position
!             update without requiring the call to be placed in this routine.
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,j
!
do i=1,mcm_ndim
 mcm_coord_maxmin(1,i) =  1.0e+20_d
 mcm_coord_maxmin(2,i) = -1.0e+20_d
enddo
!
! Update positions and store max and min values for linked list.
!
do i=1,mcm_np
 do j=1,mcm_ndim
  par(i)%x(j) = par(i)%x(j) + par(i)%v(j)*mcm_dt
  !
  if(par(i)%x(j).lt.mcm_coord_maxmin(1,j)) mcm_coord_maxmin(1,j) = par(i)%x(j)
  if(par(i)%x(j).gt.mcm_coord_maxmin(2,j)) mcm_coord_maxmin(2,j) = par(i)%x(j)
  !
 enddo
enddo
!
if(mcm_axopt.eq.4) then
 ! ensure that particles do not pass through axis
 do i=1,mcm_np
  if(par(i)%x(1).le.0.2_d*par(i)%h) then
   par(i)%x(1) = par(i)%x(1) - par(i)%v(1)*mcm_dt
   par(i)%v(1) = 0.0_d
  endif
 enddo
endif
!
! If symmetry planes are active then elastically reflect any particle that has
!  passed through an active symmetry plane
!
if(mcm_boundary) call mcm_check_sym_pen
!
end subroutine mcm_move