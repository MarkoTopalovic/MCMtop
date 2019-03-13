subroutine mcm_write2dspvp
!************************************************************************
!
!    Purpose: write 2Dascii state plot to opened file for non colocated 
!             particle case
!
!  Called by: stateplot
!
!       Date: 20-1-99
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
integer :: i,j
real(kind=real_acc) :: columns(23)
!
! write time and timestep
!
write(2,1000) mcm_ptime, mcm_timestep
!
! write variables and check for minimum and maximum values
!
do i=mcm_svp,mcm_evp
 !velocity points
  columns(1)=par(i)%x(1)
  columns(2)=par(i)%x(2)
  columns(3)=real(par(i)%mat)
  columns(4)=par(i)%v(1)
  columns(5)=par(i)%v(2)
  columns(6)=par(i)%a(1)
  columns(7)=par(i)%a(2)
  columns(8)=par(i)%rho
 !
 ! delphi crashes when numbers e-100 are written to the output file
 ! this fix prevents that from happening
 do j=1,8
  if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
 enddo
 !
 write(2,1100) (columns(j), j=1,8)
 !
 ! Max and Min
 !
 do j=1,8
  if(columns(j).lt.mcm_maxmin(1,j)) mcm_maxmin(1,j)=columns(j)
  if(columns(j).gt.mcm_maxmin(2,j)) mcm_maxmin(2,j)=columns(j)
 enddo
enddo

do i=mcm_ssp,mcm_esp
 !stress points
  columns(1)=par(i)%x(1)
  columns(2)=par(i)%x(2)
  columns(9)=real(par(i)%mat)
  columns(10)=par(i)%v(1)
  columns(11)=par(i)%v(2)
  columns(12)=par(i)%a(1)
  columns(13)=par(i)%a(2)
  columns(14)=par(i)%rho
  columns(15)=par(i)%p
  columns(16)=par(i)%sigma(1,1)
  columns(17)=par(i)%sigma(2,2)
  columns(18)=par(i)%sigma(1,2)
  columns(19)=par(i)%efps
  columns(20)=par(i)%e/par(i)%mass
  columns(21)=par(i)%c
  columns(22)=par(i)%critts
  columns(23)=par(i)%temper
 !
 ! delphi crashes when numbers e-100 are written to the output file
 ! this fix prevents that from happening
 do j=1,2
  if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
 enddo
 do j=9,23
  if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
 enddo
 !
 write(2,1100) columns(1),columns(2),(columns(j), j=9,23)
 !
 ! Max and Min
 !
 do j=1,2
  if(columns(j).lt.mcm_maxmin(1,j)) mcm_maxmin(1,j)=columns(j)
  if(columns(j).gt.mcm_maxmin(2,j)) mcm_maxmin(2,j)=columns(j)
 enddo
 do j=9,23
  if(columns(j).lt.mcm_maxmin(1,j)) mcm_maxmin(1,j)=columns(j)
  if(columns(j).gt.mcm_maxmin(2,j)) mcm_maxmin(2,j)=columns(j)
 enddo
enddo
!
1000 format(e14.6,i8)
1100 format(20e14.6)
!
end subroutine mcm_write2dspvp