subroutine mcm_write2dasc
!************************************************************************
!
!    Purpose: write 2Dascii state plot to opened file for colocated 
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
real(kind=real_acc) :: columns(30)
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
  columns(2)=par(i)%x(2)
  columns(3)=real(par(i)%mat)
  columns(4)=par(i)%v(1)
  columns(5)=par(i)%v(2)
  columns(6)=par(i)%a(1)
  columns(7)=par(i)%a(2)
  columns(8)=par(i)%rho
  columns(9)=par(i)%mass
  columns(10)=par(i)%p
  columns(11)=par(i)%sigma(1,1)
  columns(12)=par(i)%sigma(2,2)
  columns(13)=par(i)%sigma(1,2)
  select case(mcm_axopt)
   case(2)
    !Plane strain
    columns(14)=par(i)%efps
    columns(15)=par(i)%e/par(i)%mass
    columns(16)=par(i)%c
    columns(17)=par(i)%critts
    columns(18)=par(i)%temper
    columns(19)=real(par(i)%boundary)
    columns(20)=par(i)%h
    columns(21)=real(par(i)%nnbr)
   case(4)
    ! Axsymmetric
	columns(14)=par(i)%sigma(3,3)
    columns(15)=par(i)%efps
    columns(16)=par(i)%e/par(i)%mass
    columns(17)=par(i)%c
    columns(18)=par(i)%critts
    columns(19)=par(i)%temper
    columns(20)=real(par(i)%boundary)
    columns(22)=par(i)%h
    columns(23)=real(par(i)%nnbr)
  end select
 !
 ! delphi crashes when numbers e-100 are written to the output file
 ! this fix prevents that from happening
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
1100 format(22e14.6)
!
end subroutine mcm_write2dasc