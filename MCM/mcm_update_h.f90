subroutine mcm_update_h
!************************************************************************
!
!    Purpose: calculate new particle smoothing lengths
!
!  Called by: solution
!
!       Date: 09-08-2002
!
!     Errors: 
!
!      Notes: Only called if mcm_h_opt > 0
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i
integer :: j,k,l
real(kind=real_acc) :: power
real(kind=real_acc) :: min_dist, dist, factor
!
select case (mcm_h_opt)
 !
 case(1)
  ! dh/dt term from Benz's paper
  mcm_hmax = 0.0_d
  do i=mcm_ssp,mcm_esp
   !
   if (par(i)%p.ge.par(i)%pcut) then
    par(i)%hold = par(i)%h
    par(i)%h = par(i)%hold * (1.0_d + par(i)%tracerod*mcm_dt/real(mcm_ndim))
    mcm_hmax = max(mcm_hmax,par(i)%h)
   endif
  enddo
  !
 case(2)
  ! alternative dh/dt described in Benz's paper
  mcm_hmax = 0.0_d
  power = 1.0_d/real(mcm_ndim)
  do i=mcm_ssp,mcm_esp
   if(par(i)%p.ge.par(i)%pcut) then
    par(i)%hold = par(i)%h
    !
    ! THIS rho0 SHOULD BE rho AT t=to TOM 14-4-99 
    !
    par(i)%h = par(i)%h0 * ( (par(i)%rho0/par(i)%rho)**power )
    mcm_hmax = max(mcm_hmax,par(i)%h)
   endif
  enddo
  !
end select
!
end subroutine mcm_update_h