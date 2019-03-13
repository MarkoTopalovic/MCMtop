subroutine mcm_write1dasc
!************************************************************************
!
!    Purpose: write 1Dascii state plot to opened file for colocated 
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
real(kind=real_acc) :: columns(20)
!
! write time and timestep
!
write(2,1000) mcm_ptime, mcm_timestep
!
! write variables and check for minimum and maximum values
!
do i=1,mcm_np
 !
  columns(1)=par(i)%x(1)
  columns(2)=real(par(i)%mat)
  columns(3)=par(i)%v(1)
  columns(4)=par(i)%a(1)
  columns(5)=par(i)%rho
  columns(6)=par(i)%mass
  columns(7)=par(i)%p
  columns(8)=par(i)%sigma(1,1)
  columns(9)=par(i)%efps
  columns(10)=par(i)%e/par(i)%mass
  columns(11)=par(i)%c
  columns(12)=par(i)%critts
  columns(13)=par(i)%temper
  columns(14)=real(par(i)%nnbr)
  columns(15)=par(i)%h
 !
 ! delphi crashes when numbers e-100 are written to the output file
 ! this fix prevents that from happening
 !
 do j=1,mcm_out_cols(1)+mcm_ndim
  if(abs(columns(j)).lt.1.0e-50_d) columns(j) = 0.0_d
 enddo
 !
 write(2,1100) (columns(j), j=1,mcm_out_cols(1)+mcm_ndim)
 !
 ! Max and Min
 !
 do j=1,mcm_out_cols(1)+mcm_ndim
  if(columns(j).lt.mcm_maxmin(1,j)) mcm_maxmin(1,j)=columns(j)
  if(columns(j).gt.mcm_maxmin(2,j)) mcm_maxmin(2,j)=columns(j)
 enddo
enddo
!
1000 format(e14.6,i8)
1100 format(18e14.6)
!
end subroutine mcm_write1dasc