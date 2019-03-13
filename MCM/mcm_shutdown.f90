subroutine mcm_shutdown(class)
!************************************************************************
!
!    Purpose: Stop program tidily (close all files etc)
!
!  Called by: Lots
!
!       Date: 25-07-2002
!
!     Errors: None
!
!      Notes: class 1 = normal termination
!             class 2 = error termination
!
!************************************************************************
use mcm_database
!
implicit none
integer :: time
real(4) :: etime
character(len=24) :: ctime
external :: time, ctime, etime
!
integer :: class,n,i
real(4) :: cputime, ta(2)
character(len=24) :: realtime
logical :: test
!
!  write date and time of end of run to log file
!
if(class.eq.1) then
 write(*,1000)
 if (mcm_openfiles(13).eq.1) write(13,1000)
else
 write(*,1001)
 if (mcm_openfiles(13).eq.1) write(13,1001)
endif
!
 n=time()
 realtime=CTIME(n)
 cputime = etime(ta)
 if(mcm_openfiles(13).eq.1) write(13,1100) realtime
!
write(*,1200) mcm_timestep, mcm_ptime
write(*,1300) cputime
if (mcm_openfiles(13).eq.1) write(13,1200) mcm_timestep, mcm_ptime
 if(mcm_openfiles(13).eq.1) write(13,1300) cputime
!
! If writing to d3plot files, flush buffer to disk
!
if(mcm_state_opt.eq.3) then
 inquire(unit=22,OPENED= test)
 if(test) then
  call mcm_write_buffer
  close(unit=2)
 endif
endif
!
!  close open files
!
do i=1,60
 if(mcm_openfiles(i).eq.1) then
  close(unit=mcm_openfiles(i))
 endif
enddo
!
write(*,2000)
read'(A12)',realtime
stop
!
1000 format(//5x,'N o r m a l   T e r m i n a t i o n')
1001 format(//5x,'E R R O R   T E R M I N A T I O N')
1100 format(/5x,'Analysis stopped: ',a24)
1200 format(/5x,'At timestep: ',i6,'      Problem time: ',e10.4)
1300 format(/5x,'Total CPU time: ',e12.5)
2000 format(//5x,'Press RETURN to continue')
end