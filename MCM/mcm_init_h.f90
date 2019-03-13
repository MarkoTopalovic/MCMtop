subroutine mcm_init_h
!************************************************************************
!
!    Purpose: Assign initial h to all particles
!
!  Called by: initial
!
!       Date: 05-08-2002
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
!
integer :: i
!
mcm_hmax = 0.0_d
select case(mcm_init_h_opt)
 case(0)
  ! h defined for all particles of a given material
  do i=1,mcm_np
   par(i)%h = mcm_mat(par(i)%mat)%h
  enddo
  do i=1,mcm_nummat
   mcm_hmax = max(mcm_hmax,mcm_mat(par(i)%mat)%h)
  enddo
 case(1)
  ! h defined for each particle in input file
  do i=1,mcm_np
   mcm_hmax = max(mcm_hmax,par(i)%h)
  enddo
 case default
   write(*,1000) mcm_init_h_opt
   write(13,1000) mcm_init_h_opt
   call mcm_shutdown(2)
end select
!
! Set history variables
!
do i=1,mcm_np
 par(i)%h0 = par(i)%h
 par(i)%hold = par(i)%h
enddo
!
1000 format(5x,'Error in init_h, h initialisation option ',i1,' not recognised.')
!
end subroutine mcm_init_h