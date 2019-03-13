subroutine mcm_print_version(opt)
!************************************************************************
!
!    Purpose: Write program and version information to log file for a new
!             problem
!
!  Called by: startup
!
!       Date: 27-07-2002
!
!     Errors: 
!
!      Notes:
!
!************************************************************************
use mcm_database
!
implicit none
integer :: time
character(len=24) :: ctime
external :: time, ctime
!
character(len=16) :: filename
character(len=24) :: realtime
integer :: n,opt
!
select case (opt)
 case(1)
  ! write main information to screen and log file
  write(*,1000)
  write(*,1010)
  write(*,1020)
  !
  ! write filenames
  !
  filename=mcm_filein(1:mcm_filelen(1))//'.mcm'
  write(*,1100) filename
  filename=mcm_fileout(1:mcm_filelen(2))//'.log'
  write(*,1110) filename 
  !
 case(2)
  ! write main information to log file
  write(13,1000)
  write(13,1010)
  write(13,1020)
  !
  ! write filenames
  !
  filename=mcm_filein(1:mcm_filelen(1))//'.mcm'
  write(13,1100) filename
  !
  ! write date and time of start to log file
  ! 
  n=time()
  realtime=CTIME(n)
  write(13,1200) realtime
  !
end select
!
return
!
1000 format(/5x,'                  MCM',&
           //5x,'  A 1/2/3D Meshless Continuum Mechanics Program',&
            /5x,'            with material strength')
1010 format(//5x,'         Version: 2.07 W')
1020 format(5x,'       Code Date: 11-05-2007')
!
1100 format(/5x,'  Input filename: ',a16)
1110 format(5x,'    Log filename: ',a16)
1120 format(5x,' Main State File: ',a16)
!
1200 format(/5x,'Analysis started: ',a24) 
!
end subroutine mcm_print_version