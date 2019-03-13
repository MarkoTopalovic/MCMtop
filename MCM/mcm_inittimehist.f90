subroutine mcm_inittimehist
!************************************************************************
!
!    Purpose: create timehistory files and write header information to them
!
!  Called by: initial
!
!       Date: 08-08-2002
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
integer :: i, id, size, j
character(len=14) :: extension
character(len=17) :: endfmt
character(len=25) :: filename
real(kind=real_acc) :: columns(18)
logical :: fexist
!
! Start global time history file
!
filename=mcm_fileout(1:mcm_filelen(2))//'_thist.csv'
!
! if old time history file exists then delete it
!
inquire(file=filename,exist=fexist)
if(fexist) then
 open(unit=2,file=filename,status='old',form='formatted')
 close(unit=2,status='delete')
endif
!
! open new time history file
!
open(unit=2,file=filename,status='new',form='formatted')
!
! write data to file
!
write(2,1000) mcm_title
write(2,500)
!
! close file
!
close(unit=2)
!
! if no time histories have been requested then exit
!
if(mcm_thnode.gt.0) then
 !
 ! loop over specified nodes and write data to disk
 !
 do i=1,mcm_thnode
  id=mcm_thnodeid(i)
  !
  ! determine filename
  !
   size = int(log10(real(id)))+1
   write(unit=endfmt, fmt=4000) size
   write(unit=extension, fmt=endfmt) id
   filename=mcm_fileout(1:mcm_filelen(2))//extension
   !
   ! if old time history file exists then delete it
   !
   inquire(file=filename,exist=fexist)
   if(fexist) then
    open(unit=2,file=filename,status='old',form='formatted')
    close(unit=2,status='delete')
   endif
  !
  ! open new time history file
  !
  open(unit=2,file=filename,status='new',form='formatted')
  !
  ! write data to file
  !
  write(2,1000) mcm_title
  write(2,1010) id
  select case (mcm_disctype)
   case(0,1)
    select case (mcm_axopt)
     case(1)
      write(2,1100)
     case(2)
      write(2,1200)
     case(3)
      write(2,1300)
     case(4)
      write(2,1400)
    end select
  end select
  !
  ! close file
  !
  close(unit=2)
 enddo
endif
!
! if transducers have been specifed
!
if(mcm_num_transducer.gt.0) then
 !
 ! Start pressure transducer time history file
 !
 filename=mcm_fileout(1:mcm_filelen(2))//'_ptrn.csv'
 !
 ! if old time history file exists then delete it
 !
 inquire(file=filename,exist=fexist)
 if(fexist) then
  open(unit=2,file=filename,status='old',form='formatted')
  close(unit=2,status='delete')
 endif
 !
 ! open new time history file
 !
 open(unit=2,file=filename,status='new',form='formatted')
 !
 ! write data to file
 !
 write(2,1000) mcm_title
 write(2,600)
 !
 ! close file
 !
 close(unit=2)
 !
 ! Start velocity transducer time history file
 !
 filename=mcm_fileout(1:mcm_filelen(2))//'_vtrn.csv'
 !
 ! if old time history file exists then delete it
 !
 inquire(file=filename,exist=fexist)
 if(fexist) then
  open(unit=2,file=filename,status='old',form='formatted')
  close(unit=2,status='delete')
 endif
 !
 ! open new time history file
 !
 open(unit=2,file=filename,status='new',form='formatted')
 !
 ! write data to file
 !
 write(2,1000) mcm_title
 write(2,610)
 !
 ! close file
 !
 close(unit=2)
endif
!
return
!
500  format('Global time history data')
600  format('Time,Pressure transducer time histories')
610  format('Time,Velocity transducer time histories')
1000 format(a78)
1010 format('Time history data for particle ',i6)
1100 format(',Time,Coordinate,Velocity,Acceleration,Density,Pressure,Stress,Eff. Plastic Strain,', &
            'Total Internal Energy,Speed of Sound,Critical timestep,Temperature')
1200 format(',Time,x coordinate,y coordinate,x velocity,y velocity,x acceleration,y acceleration,', &
          'Density,Pressure,11 Stress,22 Stress,12 Stress, Eff. Plastic Strain,', &
		  'Total Internal Energy,Speed of Sound,Critical timestep,Temperature')
1300 format(',Time,x coordinate,y coordinate,z coordinate,x velocity,y velocity,z velocity,', &
          'Density,Pressure,11 Stress,22 Stress,33 Stress,12 Stress,13 Stress,23 Stress,', &
		  ' Eff. Plastic Strain,Total Internal Energy,Speed of Sound,Critical timestep,Temperature')
1400 format(',Time,x coordinate,y coordinate,x velocity,y velocity,x acceleration,y acceleration,', &
          'Density,Pressure,11 Stress,22 Stress,12 Stress,Hoop Stress, Eff. Plastic Strain,', &
		  'Total Internal Energy,Speed of Sound,Critical timestep,Temperature,Volumetric strain rate,rod(1,1),rod(2,2), Hoop strain rate')
!
2100 format(',Time,Coordinate,Velocity,Density,Pressure,Stress,Eff. Plastic Strain,', &
            'Total Internal Energy,Speed of Sound,Critical timestep,Temperature')
2110 format(',Time,Coordinate,Velocity,Acceleration,Density,Pressure')
2200 format(',Time,x coordinate,y coordinate,x velocity,y velocity,', &
          'Density,Pressure,11 Stress,22 Stress,12 Stress, Eff. Plastic Strain,', &
		  'Total Internal Energy,Speed of Sound,Critical timestep,Temperature')
2210 format(',Time,x coordinate,y coordinate,x velocity,y velocity,x acceleration,y acceleration,', &
          'Density')

!
4000 format('(''_th'',i',i1,',''.csv'')')
end subroutine mcm_inittimehist