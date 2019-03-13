subroutine mcm_baseaccel_in
!************************************************************************
!
!    Purpose: Read base acceleration values or set to zero
!
!  Called by: getinput
!
!       Date: 9/5/2005
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
integer lcount
character (len=80):: mssg,textline
!
if(mcm_baseaccel) then
 mssg='Error reading base accelerations card'
 call mcm_gttxsg(textline,lcount)
 !
 read(Unit=textline,fmt=100,err=500) mcm_base_a(1), mcm_base_a(2), mcm_base_a(3)
 !
endif
!
if(mcm_nthpx.ne.1) mcm_base_a(1) = 0.0_d
if(mcm_nthpy.ne.1) mcm_base_a(2) = 0.0_d
if(mcm_nthpz.ne.1) mcm_base_a(3) = 0.0_d
!
return
!
100 format(3e20.0)
500 call mcm_termin(textline,mssg,lcount,1)
!
end subroutine mcm_baseaccel_in
