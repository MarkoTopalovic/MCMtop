subroutine mcm_write1dspvp
!************************************************************************
!
!    Purpose: write 1Dascii state plot to opened file for non-colocated 
!             particle case
!
!  Called by: stateplot
!
!       Date: 07-08-2002
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
real(kind=real_acc) :: columns(18)
!
! write time and timestep
!
write(2,1000) mcm_ptime, mcm_timestep
!
! write variables and check for minimum and maximum values
!
do i=mcm_svp,mcm_evp
 ! velocity particles
  columns(1)=par(i)%x(1)
  columns(2)=real(par(i)%mat)
  columns(3)=par(i)%v(1)
  columns(4)=par(i)%a(1)
  columns(5)=par(i)%rho
 !
 ! delphi crashes when numbers e-100 are written to the output file
 ! this fix prevents that from happening
 do j=1,5
  if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
 enddo
 !
 write(2,1100) (columns(j), j=1,5)
 !
 ! Max and Min
 !
 do j=1,5
  if(columns(j).lt.mcm_maxmin(1,j)) mcm_maxmin(1,j)=columns(j)
  if(columns(j).gt.mcm_maxmin(2,j)) mcm_maxmin(2,j)=columns(j)
 enddo
enddo
!
do i=mcm_ssp,mcm_esp
 ! stress particles
  columns(1)=par(i)%x(1)
  columns(6)=real(par(i)%mat)
  columns(7)=par(i)%v(1)
  columns(8)=par(i)%a(1)
  columns(9)=par(i)%rho
  columns(10)=par(i)%p
  columns(11)=par(i)%sigma(1,1)
  columns(12)=par(i)%efps
  columns(13)=par(i)%e/par(i)%mass
  columns(14)=par(i)%c
  columns(15)=par(i)%critts
  columns(16)=par(i)%temper
 !
 ! delphi crashes when numbers e-100 are written to the output file
 ! this fix prevents that from happening
 if(abs(columns(1)).lt.1.0e-50_d) columns(1) = 0.0_d
 do j=6,16
  if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
 enddo
 !
 write(2,1100) columns(1),(columns(j), j=6,16)
 !
 ! Max and Min
 !
 if(columns(1).lt.mcm_maxmin(1,1)) mcm_maxmin(1,1)=columns(1)
 if(columns(1).gt.mcm_maxmin(2,1)) mcm_maxmin(2,1)=columns(1)
 do j=6,16
  if(columns(j).lt.mcm_maxmin(1,j)) mcm_maxmin(1,j)=columns(j)
  if(columns(j).gt.mcm_maxmin(2,j)) mcm_maxmin(2,j)=columns(j)
 enddo
enddo
!
1000 format(e14.6,i8)
1100 format(18e14.6)
!
end subroutine mcm_write1dspvp