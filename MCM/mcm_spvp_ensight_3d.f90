SUBROUTINE mcm_spvp_ensight_3d
!*****************************************************************
!
!    Purpose : This module writes out all the files for all variables              
!              
!
!  Called by : state_output
!
!       Date : 08-08-2002
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
real(kind=real_acc),dimension(6) :: p,rho,mass,efps,energy
real(kind=real_acc),dimension(6) :: c,ptimestep,temper,bound
real(kind=real_acc),dimension(6) :: h,neighbours,mat
real(kind=real_acc),dimension(3,6) :: disp, v, a
real(kind=real_acc),dimension(3,3,6) :: sigma
INTEGER r,l,rest,i,begin,j,k,m
real(kind=real_acc) :: marker(6)
character(len=25) :: extension
! 
! write the time in the time.txt file
!
open(unit=3,file=mcm_file_time,status='unknown',position='append',form='formatted')
write(3,fmt='(E13.6)') mcm_ptime
close(unit=3)

if (mcm_istate.lt.10) then
	write(unit=extension,fmt=1000) mcm_istate
else if (mcm_istate.lt.100) then
	write(unit=extension,fmt=1010) mcm_istate
else
	write(unit=extension,fmt=1020) mcm_istate
end if
!
! Open all the files and write the title
!
mcm_filename_mass       =mcm_fileout(1:mcm_filelen(2))//"_mass"//extension
mcm_filename_dens       =mcm_fileout(1:mcm_filelen(2))//"_density"//extension
mcm_filename_sig11      =mcm_fileout(1:mcm_filelen(2))//"_sigma11"//extension
mcm_filename_sig22      =mcm_fileout(1:mcm_filelen(2))//"_sigma22"//extension
mcm_filename_sig33      =mcm_fileout(1:mcm_filelen(2))//"_sigma33"//extension
mcm_filename_sig12      =mcm_fileout(1:mcm_filelen(2))//"_sigma12"//extension
mcm_filename_sig13      =mcm_fileout(1:mcm_filelen(2))//"_sigma13"//extension
mcm_filename_sig23      =mcm_fileout(1:mcm_filelen(2))//"_sigma23"//extension
mcm_filename_vel        =mcm_fileout(1:mcm_filelen(2))//"_velocity"//extension
mcm_filename_dis        =mcm_fileout(1:mcm_filelen(2))//"_displace"//extension
mcm_filename_pres       =mcm_fileout(1:mcm_filelen(2))//"_pressure"//extension
mcm_filename_strain     =mcm_fileout(1:mcm_filelen(2))//"_strain"//extension
mcm_filename_energy     =mcm_fileout(1:mcm_filelen(2))//"_energy"//extension
mcm_filename_Bound      =mcm_fileout(1:mcm_filelen(2))//"_Bound"//extension
mcm_filename_Material   =mcm_fileout(1:mcm_filelen(2))//"_material"//extension
!
open(unit=4,file=mcm_filename_pres,status='unknown',form='formatted')
open(unit=2,file=mcm_filename_dis,status='unknown',form='formatted')
open(unit=5,file=mcm_filename_vel,status='unknown',form='formatted')
open(unit=7,file=mcm_filename_dens,status='unknown',form='formatted')
open(unit=8,file=mcm_filename_mass,status='unknown',form='formatted')
open(unit=9,file=mcm_filename_sig11,status='unknown',form='formatted')
open(unit=10,file=mcm_filename_sig22,status='unknown',form='formatted')
open(unit=11,file=mcm_filename_sig33,status='unknown',form='formatted')
open(unit=12,file=mcm_filename_sig12,status='unknown',form='formatted')
open(unit=25,file=mcm_filename_sig13,status='unknown',form='formatted')
open(unit=14,file=mcm_filename_sig23,status='unknown',form='formatted')
open(unit=15,file=mcm_filename_strain,status='unknown',form='formatted')
open(unit=16,file=mcm_filename_energy,status='unknown',form='formatted')
open(unit=17,file=mcm_filename_Bound,status='unknown',form='formatted')
open(unit=21,file=mcm_filename_Material,status='unknown',form='formatted')

!
write(2,*)"displacement"
write(4,*) "pressure"
write(5,*)"velocity"
write(7,*)"density"
write(8,*)"mass"
write(9,*)"sigma11"
write(10,*)"sigma22"
write(11,*)"sigma33"
write(12,*)"sigma12"
write(25,*)"sigma13"
write(14,*)"sigma23"
write(15,*)"strain"
write(16,*)"energy"
write(17,*)"Boundary"
write(21,*)"Material"
!
do i=1,6*(mcm_np/6)-5,6
    !
	! set up arrays to be printed
	do j=1,6
	 k=i+j-1
	 !
	 p(j)=par(k)%p
     rho(j)=par(k)%rho
	 mass(j)=par(k)%mass
	 efps(j)=par(k)%efps
     energy(j)=par(k)%e/par(k)%mass
	 c(j)=par(k)%c
	 ptimestep(j)=par(k)%critts
	 temper(j)=par(k)%temper
	 bound(j)=real(par(k)%boundary)
	 h(j)=par(k)%h
	 neighbours(j)=real(par(k)%nnbr)
	 mat(j)=real(par(k)%mat)
	 !
	 do l=1,3
	  disp(l,j) = par(k)%x(l) - par(k)%xzero(l)
	  v(l,j) = par(k)%v(l)
	  a(l,j)= par(k)%a(l)
	  do m = 1,3
	   sigma(l,m,j)=par(k)%sigma(l,m)
	  enddo
	 enddo
	enddo
	!
	call mcm_write_vector(2,disp,6,i)
	call mcm_write_scalar(4,p,6,i)
	call mcm_write_vector(5,v,6,i)
	call mcm_write_scalar(7,rho,6,i)
	call mcm_write_scalar(8,mass,6,i)
	call mcm_write_tensor(9,sigma,1,1,6,i)
	call mcm_write_tensor(10,sigma,2,2,6,i)
	call mcm_write_tensor(11,sigma,3,3,6,i)
	call mcm_write_tensor(12,sigma,1,2,6,i)
	call mcm_write_tensor(25,sigma,1,3,6,i)
	call mcm_write_tensor(14,sigma,2,3,6,i)
	call mcm_write_scalar(15,efps,6,i)
	call mcm_write_scalar(16,energy,6,i)
	call mcm_write_scalar(17,bound,6,i)
	call mcm_write_scalar(21,mat,6,i)
	!
end do
!

begin=6*(mcm_np/6)
rest=mcm_np-(6*((mcm_np)/6))
if(rest/=0) then
    !
	! set up arrays to be printed
	do j=1,rest
	 k=i+j-1
	 !
	 p(j)=par(k)%p
     rho(j)=par(k)%rho
	 mass(j)=par(k)%mass
	 efps(j)=par(k)%efps
     energy(j)=par(k)%e/par(k)%mass
	 c(j)=par(k)%c
	 ptimestep(j)=par(k)%critts
	 temper(j)=par(k)%temper
	 bound(j)=real(par(k)%boundary)
	 h(j)=par(k)%h
	 neighbours(j)=real(par(k)%nnbr)
	 mat(j)=real(par(k)%mat)
	 !
	 do l=1,3
	  disp(l,j) = par(k)%x(l) - par(k)%xzero(l)
	  v(l,j) = par(k)%v(l)
	  a(l,j)= par(k)%a(l)
	  do m = 1,3
	   sigma(l,m,j)=par(k)%sigma(l,m)
	  enddo
	 enddo
	enddo
	!
	call mcm_write_vector(2,disp,rest,begin)
	call mcm_write_scalar(4,p,rest,begin)
	call mcm_write_vector(5,v,rest,begin)
	call mcm_write_scalar(7,rho,rest,begin)
	call mcm_write_scalar(8,mass,rest,begin)
	call mcm_write_tensor(9,sigma,1,1,rest,begin)
	call mcm_write_tensor(10,sigma,2,2,rest,begin)
	call mcm_write_tensor(11,sigma,3,3,rest,begin)
	call mcm_write_tensor(12,sigma,1,2,rest,begin)
	call mcm_write_tensor(25,sigma,1,3,rest,begin)
	call mcm_write_tensor(14,sigma,2,3,rest,begin)
	call mcm_write_scalar(15,efps,rest,begin)
	call mcm_write_scalar(16,energy,rest,begin)
	call mcm_write_scalar(17,bound,rest,begin)
	call mcm_write_scalar(21,mat,rest,begin)
endif
!
close(unit=2)
close(unit=4)
close(unit=5)
close(unit=7)
close(unit=8)
close(unit=9)
close(unit=10)
close(unit=11)
close(unit=12)
close(unit=25)
close(unit=14)
close(unit=15)
close(unit=16)
close(unit=17)
close(unit=21)
1000 format("00",I1)
1010 format("0",I2)
1020 format("",I3)
2000 format(22(1X,E13.6))
2010 format(28(1X,E13.6))
3000 format(3(1X,E13.6))
!
end SUBROUTINE mcm_spvp_ensight_3d