subroutine mcm_stateplot
!************************************************************************
!
!    Purpose: main subroutine for writing state output file
!
!  Called by: initial, solution 
!
!       Date: 07-08-2002
!
!     Errors: Warning if 999 state plit files have been written as no
!             more will be written.
!
!      Notes: called whenever a state plot is required, increments a 
!             counter of number of state plots (max=999)
!
!************************************************************************
!
use mcm_database
!
implicit none
!
logical fexist
character(len=5) :: extension
character(len=16) :: stfilename
!
!  only write ouput if less than 999 files have been written
!
if(mcm_istate.lt.1000) then
 !
 !  work out filename
 !
 if(mcm_istate.lt.10) then
  write(unit=extension,fmt=1000) mcm_istate
 else if(mcm_istate.lt.100) then
  write(unit=extension,fmt=1010) mcm_istate
 else
  write(unit=extension,fmt=1020) mcm_istate
 endif
 stfilename=mcm_fileout(1:mcm_filelen(2))//extension
 !
 ! if old state plot file exists then delete it
 !
 inquire(file=stfilename,exist=fexist)
 if(fexist) then
  open(unit=2,file=stfilename,status='old',form='formatted')
  close(unit=2,status='delete')
 endif
 !
 ! open new state plot file
 !
 open(unit=2,file=stfilename,status='new',form='formatted')
 !
 ! write data. MCM is only for 1D and 2D simulations
 !
 select case (mcm_disctype)
  case(0,1)
   ! colocated
   select case (mcm_ndim)
    case(1)
     call mcm_write1dasc
    case(2)
     call mcm_write2dasc
   end select
 end select
 !
 ! close state plot file
 !
 close(unit=2)
 write(13,1050) mcm_istate,mcm_ptime
 !
 ! write root state plot file (*.s000)
 !
 call mcm_writeszero
 !
 ! increment state number and check whether greater than 999
 !
 mcm_istate=mcm_istate+1
 if(mcm_istate.gt.999) then
  write(*,1100) mcm_ptime
  write(13,1100) mcm_ptime
 endif
endif
!
return
!
1000 format('.s00',i1)
1010 format('.s0',i2)
1020 format('.s',i3)
!
1050 format(/5x,'State plot ',i3,' written at ',e10.3)
!
1100 format(/5x,'* * W A R N I N G * *', &
            /5x,'State file 999 written at time ',e10.3, &
            /5x,'No further state files will be written.')
end subroutine mcm_stateplot
