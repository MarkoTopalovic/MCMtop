subroutine mcm_check_sym_pen
!************************************************************************
!
!    Purpose: Prevent particles penetrating active symmetry or periodic planes
!
!  Called by: move
!
!       Date: 06-02-2006
!
!     Errors: 
!
!      Notes: 
!             
!
!************************************************************************
!
use mcm_database
!
implicit none
!  
integer :: i
real(kind=real_acc) :: dx
!
! Check if there has been penetration through an active plane
!
! xmin
if(mcm_boundary_code(1,1).eq.1) then
 ! plane active
 if(mcm_coord_maxmin(1,1).lt.mcm_boundary_x(1,1)) then
  ! there is penetration, search for particle and correct
  do i=mcm_svp,mcm_evp
   if(par(i)%x(1).lt.mcm_boundary_x(1,1)) then
    par(i)%x(1) = 2.0_d*mcm_boundary_x(1,1) - par(i)%x(1)
	par(i)%v(1) = -par(i)%v(1)
   endif
  enddo
 endif
endif
!
! xmax
if(mcm_boundary_code(2,1).eq.1) then
 ! plane active
 if(mcm_coord_maxmin(2,1).gt.mcm_boundary_x(2,1)) then
  ! there is penetration, search for particle and correct
  do i=mcm_svp,mcm_evp
   if(par(i)%x(1).gt.mcm_boundary_x(2,1)) then
    par(i)%x(1) = 2.0_d*mcm_boundary_x(2,1) - par(i)%x(1)
	par(i)%v(1) = -par(i)%v(1)
   endif
  enddo
 endif
endif
!
if(mcm_boundary_code(1,1).eq.2) then
 !
 dx = mcm_boundary_x(2,1) - mcm_boundary_x(1,1)
 ! periodic plane active
 if(mcm_coord_maxmin(1,1).lt.mcm_boundary_x(1,1)) then
  ! there is penetration, search for particle and correct
  do i=mcm_svp,mcm_evp
   if(par(i)%x(1).lt.mcm_boundary_x(1,1)) then
    par(i)%x(1) =  par(i)%x(1) + dx
   endif
  enddo
 endif
 !
 if(mcm_coord_maxmin(2,1).gt.mcm_boundary_x(2,1)) then
  ! there is penetration, search for particle and correct
  do i=mcm_svp,mcm_evp
   if(par(i)%x(1).gt.mcm_boundary_x(2,1)) then
    par(i)%x(1) =  par(i)%x(1) - dx
   endif
  enddo
 endif
endif
!
!-------------------------------------------------------------
! ymin
if(mcm_boundary_code(1,2).eq.1) then
 ! plane active
 if(mcm_coord_maxmin(1,2).lt.mcm_boundary_x(1,2)) then
  ! there is penetration, search for particle and correct
  do i=mcm_svp,mcm_evp
   if(par(i)%x(2).lt.mcm_boundary_x(1,2)) then
    par(i)%x(2) = 2.0_d*mcm_boundary_x(1,2) - par(i)%x(2)
	par(i)%v(2) = -par(i)%v(2)
   endif
  enddo
 endif
endif
!
! ymax
if(mcm_boundary_code(2,2).eq.1) then
 ! plane active
 if(mcm_coord_maxmin(2,2).gt.mcm_boundary_x(2,2)) then
  ! there is penetration, search for particle and correct
  do i=mcm_svp,mcm_evp
   if(par(i)%x(2).gt.mcm_boundary_x(2,2)) then
    par(i)%x(2) = 2.0_d*mcm_boundary_x(2,2) - par(i)%x(2)
	par(i)%v(2) = -par(i)%v(2)
   endif
  enddo
 endif
endif
!
if(mcm_boundary_code(1,2).eq.2) then
 !
 dx = mcm_boundary_x(2,2) - mcm_boundary_x(1,2)
 ! periodic plane active
 if(mcm_coord_maxmin(1,2).lt.mcm_boundary_x(1,2)) then
  ! there is penetration, search for particle and correct
  do i=mcm_svp,mcm_evp
   if(par(i)%x(2).lt.mcm_boundary_x(1,2)) then
    par(i)%x(2) =  par(i)%x(2) + dx
   endif
  enddo
 endif
 !
 if(mcm_coord_maxmin(2,2).gt.mcm_boundary_x(2,2)) then
  ! there is penetration, search for particle and correct
  do i=mcm_svp,mcm_evp
   if(par(i)%x(2).gt.mcm_boundary_x(2,2)) then
    par(i)%x(2) =  par(i)%x(2) - dx
   endif
  enddo
 endif
endif
!
!-------------------------------------------------------------
! zmin
if(mcm_boundary_code(1,3).eq.1) then
 ! plane active
 if(mcm_coord_maxmin(1,3).lt.mcm_boundary_x(1,3)) then
  ! there is penetration, search for particle and correct
  do i=mcm_svp,mcm_evp
   if(par(i)%x(3).lt.mcm_boundary_x(1,3)) then
    par(i)%x(3) = 2.0_d*mcm_boundary_x(1,3) - par(i)%x(3)
	par(i)%v(3) = -par(i)%v(3)
   endif
  enddo
 endif
endif
!
!
if(mcm_boundary_code(2,3).eq.1) then
 ! plane active
 if(mcm_coord_maxmin(2,3).gt.mcm_boundary_x(2,3)) then
  ! there is penetration, search for particle and correct
  do i=mcm_svp,mcm_evp
   if(par(i)%x(3).gt.mcm_boundary_x(2,3)) then
    par(i)%x(3) = 2.0_d*mcm_boundary_x(2,3) - par(i)%x(3)
	par(i)%v(3) = -par(i)%v(3)
   endif
  enddo
 endif
endif
!
if(mcm_boundary_code(1,3).eq.2) then
 !
 dx = mcm_boundary_x(2,3) - mcm_boundary_x(1,3)
 ! periodic plane active
 if(mcm_coord_maxmin(1,3).lt.mcm_boundary_x(1,3)) then
  ! there is penetration, search for particle and correct
  do i=mcm_svp,mcm_evp
   if(par(i)%x(3).lt.mcm_boundary_x(1,3)) then
    par(i)%x(3) =  par(i)%x(3) + dx
   endif
  enddo
 endif
 !
 if(mcm_coord_maxmin(2,3).gt.mcm_boundary_x(2,3)) then
  ! there is penetration, search for particle and correct
  do i=mcm_svp,mcm_evp
   if(par(i)%x(3).gt.mcm_boundary_x(2,3)) then
    par(i)%x(3) =  par(i)%x(3) - dx
   endif
  enddo
 endif
endif
!
end subroutine mcm_check_sym_pen
