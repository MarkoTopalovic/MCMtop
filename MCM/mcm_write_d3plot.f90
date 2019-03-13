subroutine mcm_write_d3plot(w,nw)
!************************************************************************
!
!    Purpose: Routine to control actual write of data to d3plot files
!
!  Called by: d3plot routines
!
!       Date: 7-4-2005
!
!     Errors: 
!
!      Notes: This routine writes LS-DYNA format d3plot files
!             Writes files in 512 word blocks
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i, j, limit, remaining, nw
real(4) :: w(nw)
!
! Copy words into buffer
!
if(mcm_d3buffcnt+nw.lt.514) then
 ! won't fill buffer, just copy data into buffer
 do i=1,nw
  mcm_d3plot_buffer(mcm_d3buffcnt) = w(i)
  mcm_d3buffcnt = mcm_d3buffcnt + 1
 enddo
 !
 if(mcm_d3buffcnt.eq.513) call mcm_write_buffer
 !
else
 !
 ! will fill buffer, so have to fill buffer and write as many times as necessary
 !
 i=1
 remaining = nw
 do
  limit = 513 - mcm_d3buffcnt
  do j=1,limit
   mcm_d3plot_buffer(mcm_d3buffcnt) = w(i)
   mcm_d3buffcnt = mcm_d3buffcnt + 1
   i=i+1
  enddo
  remaining = remaining - limit
  !
  call mcm_write_buffer
  !
  if(remaining.lt.512) then
   !
   ! Have reached last section of the copy
   !
   do j=1,remaining
    mcm_d3plot_buffer(mcm_d3buffcnt) = w(i)
    mcm_d3buffcnt = mcm_d3buffcnt + 1
    i=i+1
   enddo
   !
   exit
  endif
  !
  ! Still more to do so go round again
  !
 enddo
endif
!
end subroutine mcm_write_d3plot
!
!======================================================
!
subroutine mcm_write_buffer
!
use mcm_database
!
implicit none
!
write(22,rec=mcm_d3count) mcm_d3plot_buffer
mcm_d3count = mcm_d3count + 1
mcm_d3buffcnt = 1
!
end subroutine mcm_write_buffer
