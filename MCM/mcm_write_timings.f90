subroutine mcm_write_timings
!************************************************************************
!
!    Purpose: Write status to a file
!
!  Called by: solution
!
!       Date: 09-08-2002
!
!     Errors: 
!
!      Notes: This routine uses some compiler/hardware specific routines
!
!************************************************************************
!
use mcm_database
!
implicit none
integer :: time
real(4) :: etime
character(len=24) :: ctime
external :: time, ctime, etime
!
integer :: n
real(4) :: i, ta(2)
logical :: fexist
character(len=24) :: realtime
character(len=25) :: filename
!
! filename
!
filename=mcm_fileout(1:mcm_filelen(2))//'.status'
!
if(mcm_timestep.eq.1) then
 !
 ! if old file exists then delete it
 !
 inquire(file=filename,exist=fexist)
 if(fexist) then
  open(unit=2,file=filename,status='old',form='formatted')
  close(unit=2,status='delete')
 endif
endif
!
! open new time history file
!
open(unit=2,file=filename,form='formatted',position='append')
!
! For first timesep write column headers
!
if(mcm_timestep.eq.1) then
 write(2,1000)
endif
!
n=time()
realtime=CTIME(n)
i = etime(ta)
!
write(2,2000) mcm_timestep,mcm_ptime,mcm_dt,i,mcm_maxnbr,realtime
!
close(unit=2)
return
!
1000 format('Timestep    Problem time            dt        CPU Time   Maxnbr   Real Time')
2000 format(i7,5x,e12.5,2x,e12.5,2x,e14.7,2x,i7,2x,a24)
!
end subroutine mcm_write_timings