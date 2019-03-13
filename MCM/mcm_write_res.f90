subroutine mcm_write_res
!*****************************************************************
!
!    Purpose : This module writes out the .res file
!
!  Called by : state_output
!
!       Date : 07-08-2002
!
!     Errors : 
!
!      Notes : 
!
!*****************************************************************
!
use mcm_database
!
IMPLICIT NONE 
!
integer :: rest,l,r
REAL(kind=real_acc), DIMENSION(6) :: timenum
CHARACTER(LEN=25) :: file_write
!
file_write=mcm_fileout(1:mcm_filelen(2))//".res"
open(unit=2,file=file_write,status='unknown',form='formatted')
!
! Write header
! these three numbers mean
! # of scalars  -----  in 3D add sigma33, sigma13 and sigma23
! # of of vectors
! geom_chang_flag
!
select case(mcm_ndim)
 case(2)
  if (mcm_axopt.eq.4) then
	write(2,*) "16 3 0" ! Hoop stress must be added
  else 
	write(2,*) "15 3 0"
  end if
 case(3)
  write(2,*) "18 3 0"
end select
!
write(2,fmt='(I3)') mcm_istate-1   ! # of time steps
!
! Read the time for all plots and write them in the .res file
!
open(unit=3,file=mcm_file_time,status='unknown',form='formatted',position='rewind')
if (mcm_istate>6) then
 do l=1,6*((mcm_istate-1)/6)-5,6
  do r=1,6
   read(3,fmt='(E13.6)') timenum(r)
  end do
  write(2,fmt='(6(E10.3))') timenum(1),timenum(2),timenum(3),&
						    timenum(4),timenum(5),timenum(6)
 end do
end if
!
rest=mcm_istate-(6*((mcm_istate-1)/6))-1
if (rest/=0) then
 do l=1,rest
  read(3,fmt='(1X,E13.6)') timenum(l)
 end do
 select case(rest)
   case(1)
	write(2,fmt='(E10.3)') timenum(1)
   case(2)
	write(2,fmt='(2(E10.3))') timenum(1),timenum(2)
   case(3)
	write(2,fmt='(3(E10.3))') timenum(1),timenum(2),timenum(3)
   case(4)
	write(2,fmt='(4(E10.3))') timenum(1),timenum(2),timenum(3),timenum(4)
   case(5)
	write(2,fmt='(5(E10.3))') timenum(1),timenum(2),timenum(3),timenum(4),timenum(5)
  end select
end if
!
write(2,*) "1 1"  ! start_file_# and skip_by_value
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_density***"," density"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_pressure***"," pressure"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_material***"," material"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_mass***"," mass"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_sigma11***"," sigma11"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_sigma22***"," sigma22"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_sigma12***"," sigma12"
if (mcm_ndim.eq.3) then					! in the 3D case some variables must be added
	write(2,*) mcm_fileout(1:mcm_filelen(2))//"_sigma33***"," sigma33"
	write(2,*) mcm_fileout(1:mcm_filelen(2))//"_sigma13***"," sigma13"
	write(2,*) mcm_fileout(1:mcm_filelen(2))//"_sigma23***"," sigma23"
end if
if (mcm_axopt.eq.4) write(2,*) mcm_fileout(1:mcm_filelen(2))//"_hoop***"," hoop_stress"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_energy***"," Internal_Energy"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_strain***"," Eff_Plastic_Strain"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_C***"," Sound Speed"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_PCT***"," Part_Crit_Timestep"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_Temp***"," Temperature"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_Bound***"," Boundary_Particles"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_Smooth***"," Smoothing_length"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_Neigh***"," Nber_of_Neigbours"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_displace***"," displacement"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_velocity***"," velocity"
write(2,*) mcm_fileout(1:mcm_filelen(2))//"_accelera***"," acceleration"
!
close(unit=3)
!
1000 format("00",I1)
1010 format("0",I2)
1020 format("",I3)
!
END SUBROUTINE mcm_write_res 

