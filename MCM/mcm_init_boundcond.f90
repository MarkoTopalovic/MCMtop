subroutine mcm_init_boundcond
!************************************************************************
!
!    Purpose: Apply boundary conditions to initial velocities
!
!  Called by: Initialise
!
!       Date: 01-08-2002
!
!     Errors: 
!
!      Notes: Makes sure that initial velocity boundary conditions
!             do not conflict with nodal constraints. The nodal constraints
!             will override the initial velocities.
!
!************************************************************************
!
use mcm_database
!
implicit none
integer :: i
!
! Apply nodal displacement boundary conditions
!
do i=mcm_svp,mcm_evp
 select case (par(i)%dispbc)
  case(1)
   ! constrained x displacement
   par(i)%v(1) = 0.0
  case(2)
   ! constrained y displacement
   par(i)%v(2) = 0.0
  case(3)
   ! constrained z displacement
   par(i)%v(3) = 0.0
  case(4)
   ! constrained x and y displacement
   par(i)%v(1) = 0.0
   par(i)%v(2) = 0.0
  case(5)
   ! constrained y and z displacement
   par(i)%v(2) = 0.0
   par(i)%v(3) = 0.0
  case(6)
   ! constrained x and z displacement
   par(i)%v(1) = 0.0
   par(i)%v(3) = 0.0
  case(7)
   ! constrained x,y and z displacement
   par(i)%v(1) = 0.0
   par(i)%v(2) = 0.0
   par(i)%v(3) = 0.0
 end select
enddo
!
end subroutine mcm_init_boundcond