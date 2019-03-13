subroutine mcm_pressure_transducer
!************************************************************************
!
!    Purpose: Interpolate  for pressures at specified coordinates
!
!  Called by: Time history
!
!       Date: 12_8_2005
!
!     Errors: 
!
!      Notes:
!
!************************************************************************
use mcm_database
!
implicit none
!
integer :: i,j,k,l,size
real(kind=real_acc) :: sum, Vj, W,interp_pres_old
real(kind=real_acc) :: interp_pres(mcm_num_transducer), interp_v(3,mcm_num_transducer)
character(len=14) :: extension
character(len=17) :: endfmt
character(len=22) :: thformat
character(len=25) :: filename
!
! interpolate for pressures at transducer points
!
interp_pres = 0.0_d
interp_v = 0.0_d
!
do i=1,mcm_num_transducer
 !
 sum = 0.0_d
 if(mcm_num_tr_nbr(i).gt.1) then
   do k=1,mcm_num_tr_nbr(i)
    j=mcm_tr_nbr(k,i)
    Vj = par(j)%mass/par(j)%rho
    call mcm_kernel(W,mcm_trans_x(1,i),par(j)%x,par(j)%h)
    interp_pres(i) = interp_pres(i) + Vj*par(j)%p*W
	do l=1,3
	 interp_v(l,i) = interp_v(l,i) + Vj*par(j)%v(l)*W
	enddo
   enddo
 endif
enddo
!
!==================================================================
! write out pressures to file
!
filename=mcm_fileout(1:mcm_filelen(2))//'_ptrn.csv'
!
! open file, with file positioned at its end
!
open(unit=2,file=filename,form='formatted',position='append')
!
! Workout format for output
!
if (mcm_num_transducer.lt.9) then
 write(unit=thformat,fmt=500) mcm_num_transducer+1 
else
 write(unit=thformat,fmt=600) mcm_num_transducer+1
endif
!
! Write pressure to file
!
do j=1,mcm_num_transducer
 if(abs(interp_pres(j)).lt.1.0e-40_d) interp_pres(j) = 0.0
enddo
write(unit=2,fmt=thformat) mcm_ptime,(interp_pres(j), j=1,mcm_num_transducer)
!
! close file
!
close(unit=2)
!
!==================================================================
! write out velocities to file
!
filename=mcm_fileout(1:mcm_filelen(2))//'_vtrn.csv'
!
! open file, with file positioned at its end
!
open(unit=2,file=filename,form='formatted',position='append')
!
! Workout format for output
!
if (mcm_num_transducer.lt.3) then
 write(unit=thformat,fmt=500) (3*mcm_num_transducer)+1 
else
 write(unit=thformat,fmt=600) (3*mcm_num_transducer)+1
endif
!
! Write velocity to file
!
do j=1,mcm_num_transducer
 do k=1,3
  if(abs(interp_v(k,j)).lt.1.0e-40_d) interp_v(k,j) = 0.0
 enddo
enddo
write(unit=2,fmt=thformat) mcm_ptime,((interp_v(k,j),k=1,3), j=1,mcm_num_transducer)
!
! close file
!
close(unit=2)
!====================================================================================
return
!
500 format('(',i1,'(e14.6,'',''))')
600 format('(',i2,'(e14.6,'',''))')
!
end subroutine mcm_pressure_transducer

