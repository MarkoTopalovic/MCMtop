subroutine mcm_setuplist
!************************************************************************
!
!    Purpose: set up linked list
!
!  Called by: neighbours
!
!       Date: 05-08-2002
!
!     Errors: Insufficient memory for linked list arrays
!
!      Notes: This subroutine calculates the linked list arrays.
!             It uses current coordinate max and min values to set up
!             the grid.
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,j,k, error
integer :: boxcrd(3)
logical :: in
real(kind=real_acc) :: factor
!
! setup underlying grid
! cell size is 2 times the maximum smoothing length
!
factor = 2.0_d
!if(axopt.eq.4) factor=3.0
do i=1,3
 mcm_gridmin(i)=0.0_d
 mcm_gridlim(i)=1
 mcm_gridsize(i)=1.0_d/(factor*mcm_hmax)
enddo
!
! Work out box coordinates
!
do i=1,mcm_ndim
 ! Initial setting of minimum coord
 select case (mcm_boundary_type(i))
  case(1)
   ! Both free 
   mcm_gridmin(i) = mcm_coord_maxmin(1,i)-0.01_d*mcm_hmax
   mcm_gridlim(i) = int(mcm_gridsize(i)*(mcm_coord_maxmin(2,i)-mcm_gridmin(i))) + 1
  case(2)
   ! Min fixed, max free
   mcm_gridmin(i) = mcm_boundary_x(1,i)
   mcm_gridlim(i) = int(mcm_gridsize(i)*(mcm_coord_maxmin(2,i)-mcm_gridmin(i))) + 1
  case(3)
   ! Min free, max fixed
   mcm_gridlim(i) = int(mcm_gridsize(i)*(mcm_boundary_x(2,i)-mcm_coord_maxmin(1,i))) + 1
   mcm_gridmin(i) = mcm_boundary_x(2,i) - mcm_gridlim(i)*factor*mcm_hmax
  case(4)
   ! Both fixed
   mcm_gridmin(i) = mcm_boundary_x(1,i)
   ! Adjust gridsize to ensure an integer number of cells
   mcm_gridlim(i) = int(mcm_gridsize(i)*(mcm_boundary_x(2,i)-mcm_boundary_x(1,i)))
   mcm_gridlim(i) = max(mcm_gridlim(i),1)
   mcm_gridsize(i) = real(mcm_gridlim(i))/(mcm_boundary_x(2,i)-mcm_boundary_x(1,i))
 end select
enddo
!
! Adjust box coordinates to allow for symmetry and periodic planes
!
do i=1,mcm_ndim
 ! Min coord edge
 select case(mcm_boundary_code(1,i))
  case(1,2)
   ! periodic or symmetry bc
   mcm_gridlim(i) = mcm_gridlim(i) + 1
   mcm_gridmin(i) = mcm_gridmin(i) - (1.0_d/mcm_gridsize(i))
 end select
 ! Max coord edge
 select case(mcm_boundary_code(2,i))
  case(1,2)
   ! periodic or symmetry bc
   mcm_gridlim(i) = mcm_gridlim(i) + 1
 end select
enddo
!
! Correct box coordinates have now been set so 
! allocate memory for linked list array
allocate (mcm_llgrid(mcm_gridlim(1),mcm_gridlim(2),mcm_gridlim(3)),STAT=error)  
if(error.ne.0) then
  write(*,1000)
  write(13,1000)
  write(*,1100) mcm_timestep
  write(13,1100) mcm_timestep
  call mcm_shutdown(2)
endif 
!
! initialise linked list
!
do k=1,mcm_gridlim(3)
 do j=1,mcm_gridlim(2)
  do i=1,mcm_gridlim(1)
   mcm_llgrid(i,j,k)=0
  enddo
 enddo
enddo
!
!  now set up list
!
do i=1,3
 boxcrd(i) = 1
enddo
!
do i=1,mcm_np
 if(par(i)%active) then  ! only add to list if particle is active
  in = .true.
 do j=1,mcm_ndim
   boxcrd(j)=int(mcm_gridsize(j)*(par(i)%x(j)-mcm_gridmin(j))) + 1
   !
   ! check that particle is in the sort domain
   !
   if(boxcrd(j).lt.1) in=.false.
   if(boxcrd(j).gt.mcm_gridlim(j)) in=.false.
 enddo
 !
  if(in) then
   !
 ! move pointer to first particle in box to the particle pointer
 !
 par(i)%llpointer=mcm_llgrid(boxcrd(1),boxcrd(2),boxcrd(3))
 !
 ! set the pointer to the first particle in the box
 !
 mcm_llgrid(boxcrd(1),boxcrd(2),boxcrd(3))=i
  else
   par(i)%llpointer = -1
  endif
 endif
enddo
!
return
!
1000 format(//5x,'Error in subroutine setuplist')
1100 format(/5x,'Insufficient memory for linked list array. Timestep: ',i6)
!
end subroutine mcm_setuplist