subroutine mcm_boundcond
!************************************************************************
!
!    Purpose: Apply boundary conditions to particle accelerations
!
!  Called by: Solution
!
!       Date: 09-08-2002
!
!     Errors: 
!
!      Notes:
!
!************************************************************************
!
use mcm_database
!
implicit none
integer :: i,j
!
! Apply base accelerations
!
if(mcm_baseaccel) then
 do i=mcm_svp,mcm_evp
  do j=1,mcm_ndim
   par(i)%a(j) = par(i)%a(j) + mcm_base_a(j)
  enddo
 enddo
endif
!
! Apply nodal displacement boundary conditions
!
do i=mcm_svp,mcm_evp
 select case (par(i)%dispbc)
  case(1)
   ! constrained x displacement
   par(i)%a(1) = 0.0_d
  case(2)
   ! constrained y displacement
   par(i)%a(2) = 0.0_d
  case(3)
   ! constrained z displacement
   par(i)%a(3) = 0.0_d
  case(4)
   ! constrained x and y displacement
   par(i)%a(1) = 0.0_d
   par(i)%a(2) = 0.0_d
  case(5)
   ! constrained y and z displacement
   par(i)%a(2) = 0.0_d
   par(i)%a(3) = 0.0_d
  case(6)
   ! constrained x and z displacement
   par(i)%a(1) = 0.0_d
   par(i)%a(3) = 0.0_d
  case(7)
   ! constrained x,y and z displacement
   par(i)%a(1) = 0.0_d
   par(i)%a(2) = 0.0_d
   par(i)%a(3) = 0.0_d
 end select
enddo
!
end subroutine mcm_boundcond