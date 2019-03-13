subroutine mcm_startup(newproblem)
!************************************************************************
!
!    Purpose: Get number of command line arguments, parse command line,
!             if no arguments are given ask analyst for them. Check that
!             input file does exist, output files do not exist.
!
!  Called by: MAIN
!
!       Date: 25-07-2002
!
!     Errors: too many command line arguments
!             input file does not exist
!             output files already exist for new problem
!
!      Notes:
!
!************************************************************************
use mcm_database
!
implicit none
integer :: iargc
external iargc
!
character(len=1) :: question
character(len=16) :: filetest
logical :: newproblem, fexist
integer :: numarg
integer :: i
!
! Calculate constants
!
pi = acos(-1.0_d)
!
!  get number of command line arguments
!
numarg=IARGC()
!
!  error termination if more than two arguments
!
if(numarg.gt.2) then
 write(*,1000)
 write(*,1010)
 call mcm_shutdown(2)
endif
!
!  if no command line arguments are given then ask for them
!
if(numarg.eq.0) then
 write(*,1050)
! ask new problem or restart
 write(*,1060)
 read '(A1)',question
 if(question.eq."n") then
!new problem
  newproblem=.TRUE.
  write(*,1070)
  read '(A12)',mcm_filein
  write(*,1080)
  read '(A12)',mcm_fileout
  if(mcm_fileout.eq."") mcm_fileout=mcm_filein
 else
! restart
  newproblem=.FALSE.
  write(*,1090)
  read '(A16)',mcm_filerestart
  write(*,1100)
  read '(A12)',mcm_filein
 endif
else
! if command line arguments are given then parse command line
 call mcm_parse_cl(numarg,newproblem)
endif
!
! for a new problem
!
if(newproblem) then
 ! get filename lengths
 mcm_filelen=12
 do i=1,12
  if(mcm_filelen(1).eq.12) then
   if(mcm_filein(i:i).eq.' ') then
    mcm_filelen(1)=i-1
   endif
  endif
  if(mcm_filelen(2).eq.12) then
   if(mcm_fileout(i:i).eq.' ') then
    mcm_filelen(2)=i-1
   endif
  endif
 enddo
 !initialise list of open files
 do i=1,60
  mcm_openfiles(i)=0
 enddo
 ! check that input file exists, if it does not then error termination
 filetest=mcm_filein(1:mcm_filelen(1))//'.mcm'
 inquire(file=filetest,exist=fexist)
 if(.not.fexist) then
  write(*,1000)
  write(*,1110) filetest
  call mcm_shutdown(2)
 endif
 !
 ! write program and version information to log file
 !
 call mcm_print_version(1)
else
 !
 ! for a restart
 !
 ! get filename lengths
 mcm_filelen=12
 do i=1,12
  if(mcm_filelen(1).eq.12) then
   if(mcm_filein(i:i).eq.' ') then
    mcm_filelen(1)=i-1
   endif
  endif
  if(mcm_filelen(3).eq.12) then
   if(mcm_filerestart(i:i).eq.' ') then
    mcm_filelen(3)=i-1
   endif
  endif
 enddo
! check that restart file exists, if it does not then error termination
 inquire(file=mcm_filerestart,exist=fexist)
 if(.not.fexist) then
  write(*,1000)
  write(*,1140) mcm_filerestart
  call mcm_shutdown(2)
 endif
 ! check that input file exists (if specified), if it does not then error termination
 if(mcm_filelen(1).gt.0) then
  inquire(file=mcm_filein,exist=fexist)
  if(.not.fexist) then
   write(*,1000)
   write(*,1150) mcm_filein
   call mcm_shutdown(2)
  endif
 endif
 ! note that the log file will be opened in the restart subroutine once the old
 ! restart filename has been loaded
endif
return
!
1000 format(//5x,'Error in subroutine startup')
1010 format(/5x,'Too many command line arguments, maximum is 2')
!
1050 format(//5x,'No command line arguments')
1060 format(/5x,' New problem or restart? (n/r): ',$)
1070 format(/5x,' Input file: ',$)
1080 format(/5x,'Output file: ',$)
1090 format(/5x,'Restart file (include correct file extension): ',$)
1100 format(/5x,'Restart input file (leave blank if none): ',$)
!
1110 format(/5x,'Input file ',a12,' does not exist.')
1120 format(/5x,'State file ',a16,' already exists.')
1130 format(/5x,'Log file ',a16,' already exists.')
1140 format(/5x,'Restart file ',a16,' does not exist.')
1150 format(/5x,'Restart input file ',a12,' does not exist.')
!
end subroutine mcm_startup