subroutine mcm_parse_cl(numarg,newproblem)
!************************************************************************
!
!    Purpose: get filenames from command line
!
!  Called by: startup
!
!       Date: 25-07-2002
!
!     Errors: various syntax errors in command line cause error terminations
!
!      Notes: format of comand line arguments is
!               -i <input file name> -o <output file name>
!          or:  -i <input file name>
!          or:  -r <restart file name> -i <restart input filename>
!
!************************************************************************
use mcm_database
!
implicit none
external getarg
!
integer :: numarg,inputloc,outloc,restloc,argfound
integer :: i
! integer (2) :: i
logical :: newproblem
character(len=18) :: arguments(2)
!
! load in command line arguments
!
do i=1,numarg
 call getarg(i,arguments(i))
enddo
!
!  identify file names
!  check that an input or restartfile is specified, check that two input/output/restart
!  filenames have not been specified.
!
inputloc=0
outloc=0
restloc=0
do i=1,numarg
 argfound=0
 if(arguments(i)(1:1).eq.'i') then
  ! found an input file
  if(inputloc.eq.0) then
   ! no input file already given
   inputloc=i
   argfound=1
  else
   ! two input files given error termination
   write(*,1000)
   write(*,1010)
   call mcm_shutdown(2)
  endif
 endif
 if(arguments(i)(1:1).eq.'o') then
  ! found an output file
  if(outloc.eq.0) then
   ! no output file already specified
   outloc=i
   argfound=1
  else
   ! two output files given, error termination
   write(*,1000)
   write(*,1020)
   call mcm_shutdown(2)
  endif
 endif
 if(arguments(i)(1:1).eq.'r') then
  ! found a restart file
  if(restloc.eq.0) then
   ! no output file already specified
   restloc=i
   argfound=1
  else
   ! two restart files given, error termination
   write(*,1000)
   write(*,1030)
   call mcm_shutdown(2)
  endif
 endif
 if(argfound.eq.0) then
  ! incorrect command line argument
  write(*,1000)
  write(*,1040) i
  call mcm_shutdown(2)
 endif
enddo
!
! if no input or restart file given then error termination
!
if(inputloc.eq.0.and.restloc.eq.0) then
 write(*,1000)
 write(*,1050)
 call mcm_shutdown(2)
endif
!
! if both restart and output files given warn that output file name will be ignored
!
if(restloc.ne.0.and.outloc.ne.0) then
 write(*,1070)
endif
!
!  check that 2nd character of argument is =
!
do i=1,numarg
 if(arguments(i)(2:2).ne.'=') then
  write(*,1080) i
  call mcm_shutdown(2)
 endif
enddo
!
!  get input and output names or get restart name
!
if(inputloc.ne.0.and.restloc.eq.0) then
 newproblem=.TRUE.
 mcm_filein=arguments(inputloc)(3:)
 if(outloc.eq.0) then
  mcm_fileout=arguments(inputloc)(3:)
 else
  mcm_fileout=arguments(outloc)(3:)
 endif
else
 newproblem=.FALSE.
 mcm_filerestart=arguments(restloc)(3:)
 if(inputloc.ne.0) then
  mcm_filein=arguments(inputloc)(3:)
 endif
endif
!
return
1000 format(//5x,'Error in subroutine parse_cl')
1010 format(/5x,'Two input filenames specified')
1020 format(/5x,'Two output filenames specified')
1030 format(/5x,'Two restart filenames specified')
1040 format(/5x,'Command line argument ',i1,' is unknown.',&
            /5x,'Only i=<input filename>, o=<output filename> and r=<restart filename>',&
            /5x,'are recognised command line arguments.')
1050 format(/5x,'No input or restart file specified')
1070 format(/5x,'WARNING: both restart and output filenames given',&
           &/5x,'         Output filename will be ignored and output filename of',&
		    /5x,'         original run used')
1080 format(/5x,'Error in command line argument ',i1,&
            /5x,'Second character of argument must be "="')
end subroutine mcm_parse_cl