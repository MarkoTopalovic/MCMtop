subroutine mcm_timehist
!************************************************************************
!
!    Purpose: Write time history data for specified nodes to disk
!
!  Called by: initial, solution
!
!       Date: 08-08-2002
!
!     Errors: 
!
!      Notes: Writes data in csv format
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
character(len=22) :: thformat
character(len=25) :: filename
real(kind=real_acc) :: columns(30)
!
! Write global time history file
!
filename=mcm_fileout(1:mcm_filelen(2))//'_thist.csv'
!
! open file, with file positioned at its end
!
open(unit=2,file=filename,form='formatted',position='append')
!
! Workout format for output
!
if (mcm_nummat.lt.2) then
 write(unit=thformat,fmt=500) mcm_nummat*5+1 
else
 write(unit=thformat,fmt=600) mcm_nummat*5+1
endif
!
! Write ke and ie by material to file
!
write(unit=2,fmt=thformat) mcm_ptime,(mcm_mat(j)%ke,mcm_mat(j)%ie, &
         & mcm_mat(j)%xcf, mcm_mat(j)%ycf, mcm_mat(j)%zcf, j=1,mcm_nummat)
!
! close file
!
close(unit=2)
!
!
! Transducers
!
if(mcm_num_transducer.gt.0) call mcm_pressure_transducer
!
! if no time histories have been requested then exit
!
if(mcm_thnode.eq.0) return
!
! loop over specified nodes and write data to disk
!
do i=1,mcm_thnode
 id=mcm_thnodeid(i)
 !
 ! determine filename
 !
 size = int(log10(real(id)))+1
 write(unit=endfmt, fmt=2000) size
 write(unit=extension, fmt=endfmt) id
 filename=mcm_fileout(1:mcm_filelen(2))//extension
 !
 ! open file, with file positioned at its end
 !
 open(unit=2,file=filename,form='formatted',position='append')
 !
 ! write data to file
 !
 select case (mcm_disctype)
  case(0,1)
  ! colocated
   select case (mcm_axopt)
    case(1)
     columns(1)=par(id)%x(1)
     columns(2)=par(id)%v(1)
     columns(3)=par(id)%a(1)
     columns(4)=par(id)%rho
     columns(5)=par(id)%p
     columns(6)=par(id)%sigma(1,1)
     columns(7)=par(id)%efps
     columns(8)=par(id)%e/par(id)%mass
     columns(9)=par(id)%c
     columns(10)=par(id)%critts
     columns(11)=par(id)%temper
     !
     do j=1,11
      if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
     enddo
	 !
     write(2,1100) mcm_ptime,(columns(j), j=1,11)
    case(2)
     columns(1)=par(id)%x(1)
     columns(2)=par(id)%x(2)
     columns(3)=par(id)%v(1)
     columns(4)=par(id)%v(2)
     columns(5)=par(id)%a(1)
     columns(6)=par(id)%a(2)
     columns(7)=par(id)%rho
     columns(8)=par(id)%p
     columns(9)=par(id)%sigma(1,1)
     columns(10)=par(id)%sigma(2,2)
     columns(11)=par(id)%sigma(1,2)
     columns(12)=par(id)%efps
     columns(13)=par(id)%e/par(id)%mass
     columns(14)=par(id)%c
     columns(15)=par(id)%critts
     columns(16)=par(id)%temper
     !
     do j=1,16
      if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
     enddo
	 !
     write(2,1100) mcm_ptime,(columns(j), j=1,16)
	case(3)
	 columns(1)=par(id)%x(1)
     columns(2)=par(id)%x(2)
     columns(3)=par(id)%x(3)
     columns(4)=par(id)%v(1)
     columns(5)=par(id)%v(2)
     columns(6)=par(id)%v(3)
     columns(7)=par(id)%rho
     columns(8)=par(id)%p
     columns(9)=par(id)%sigma(1,1)
     columns(10)=par(id)%sigma(2,2)
     columns(11)=par(id)%sigma(3,3)
     columns(12)=par(id)%sigma(1,2)
     columns(13)=par(id)%sigma(1,3)
     columns(14)=par(id)%sigma(2,3)
     columns(15)=par(id)%efps
     columns(16)=par(id)%e/par(id)%mass
     columns(17)=par(id)%c
     columns(18)=par(id)%critts
     columns(19)=par(id)%temper
     !
     do j=1,19
      if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
     enddo
	 !
     write(2,1200) mcm_ptime,(columns(j), j=1,19)
	case(4)
     columns(1)=par(id)%x(1)
     columns(2)=par(id)%x(2)
     columns(3)=par(id)%v(1)
     columns(4)=par(id)%v(2)
     columns(5)=par(id)%a(1)
     columns(6)=par(id)%a(2)
     columns(7)=par(id)%rho
     columns(8)=par(id)%p
     columns(9)=par(id)%sigma(1,1)
     columns(10)=par(id)%sigma(2,2)
     columns(11)=par(id)%sigma(1,2)
	 columns(12)=par(id)%sigma(3,3)
     columns(13)=par(id)%efps
     columns(14)=par(id)%e/par(id)%mass
     columns(15)=par(id)%c
     columns(16)=par(id)%critts
     columns(17)=par(id)%temper
	 columns(18)=par(id)%tracerod
	 columns(19)=par(id)%rod(1,1)
	 columns(20)=par(id)%rod(2,2)
	 columns(21)=par(id)%rod(3,3)
     !
     do j=1,21
      if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
     enddo
	 !
     write(2,1200) mcm_ptime,(columns(j), j=1,21)
   end select
  case(2)
   ! non-colocated
   select case (mcm_axopt)
    case(1)
	 if(id.le.mcm_esp) then
      columns(1)=par(id)%x(1)
      columns(2)=par(id)%v(1)
      columns(3)=par(id)%rho
      columns(4)=par(id)%p
      columns(5)=par(id)%sigma(1,1)
      columns(6)=par(id)%efps
      columns(7)=par(id)%e/par(id)%mass
      columns(8)=par(id)%c
      columns(9)=par(id)%critts
      columns(10)=par(id)%temper
      !
      do j=1,10
       if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
      enddo
	  !
      write(2,1100) mcm_ptime,(columns(j), j=1,10)
     else 
	  columns(1)=par(id)%x(1)
      columns(2)=par(id)%v(1)
	  columns(3)=par(id)%a(1)
      columns(4)=par(id)%rho
      !
      do j=1,4
       if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
      enddo
	  !
      write(2,1100) mcm_ptime,(columns(j), j=1,4)
     endif
    case(2)
	 if(id.le.mcm_esp) then
      columns(1)=par(id)%x(1)
      columns(2)=par(id)%x(2)
      columns(3)=par(id)%v(1)
      columns(4)=par(id)%v(2)
      columns(5)=par(id)%rho
      columns(6)=par(id)%p
      columns(7)=par(id)%sigma(1,1)
      columns(8)=par(id)%sigma(2,2)
      columns(9)=par(id)%sigma(1,2)
      columns(10)=par(id)%efps
      columns(11)=par(id)%e/par(id)%mass
      columns(12)=par(id)%c
      columns(13)=par(id)%critts
      columns(14)=par(id)%temper
      !
      do j=1,14
       if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
      enddo
	  !
      write(2,1100) mcm_ptime,(columns(j), j=1,14)
     else
      columns(1)=par(id)%x(1)
      columns(2)=par(id)%x(2)
      columns(3)=par(id)%v(1)
      columns(4)=par(id)%v(2)
      columns(5)=par(id)%a(1)
      columns(6)=par(id)%a(2)
      columns(7)=par(id)%rho
      !
      do j=1,7
       if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
      enddo
	  !
      write(2,1100) mcm_ptime,(columns(j), j=1,7)
	 endif
   end select
 end select
 !
 ! close file
 !
! IF ((ABS(mcm_ptime-48.0).LE.1.0).OR.(ABS(mcm_ptime-300.0).LE.1.0)) THEN ! hardwired - tom
! DO j = 45,mcm_np,100
!    WRITE(2,'(2E12.5)') xforce,yforce
! ENDDO
! ENDIF
 close(unit=2)
!
enddo
!
return
!
500 format('(',i1,'('','',e14.6))')
600 format('(',i2,'('','',e14.6))')
!
1100 format(','16(e14.6,','),e14.6)
1200 format(','21(e14.6,','),e14.6)
2000 format('(''_th'',i',i1,',''.csv'')')
!
end subroutine mcm_timehist
