subroutine mcm_d3plot_state
!************************************************************************
!
!    Purpose: Set up and write current state to d3plot file
!
!  Called by: mcm_state_output
!
!       Date: 7-4-2005
!
!     Errors: 
!
!      Notes: This routine writes LS-DYNA format d3plot files
!             From LSDYNA database manual - default max size for
!             plot file is 7x512x512 words (ie 3584 calls of write_buffer)
!             This routine is set up to write at least one state to a file,
!             and not to split states across files.
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,j
integer :: req_space, num_blocks

character(len=7) :: extension
character(len=19) :: stfilename
!
! file is written in 32 bit words

real(4) :: data_type
real(4), dimension(3) :: x
real(4), dimension(19) :: d3data
!
! Check if need to open a new file
!
if(mcm_d3count.gt.1) then  ! write at least one state into file
 req_space = 7 + 6*mcm_nummat + 9*mcm_max_np + 19*mcm_max_np  ! no. words for one state
 num_blocks = req_space/512
 if(mcm_d3count+num_blocks.gt.3584) then
  ! will overrun file, close and open next one
  if(mcm_d3file.ge.999) return
  !
  call mcm_write_buffer  ! flush buffer
  close(unit=22)
  !
  ! Get new filename
  !
  mcm_d3file = mcm_d3file + 1
  if(mcm_d3file.lt.10) then
   write(unit=extension,fmt=1000) mcm_d3file
  else if(mcm_d3file.lt.100) then
   write(unit=extension,fmt=1010) mcm_d3file
  else
   write(unit=extension,fmt=1020) mcm_d3file
  endif 
  stfilename=mcm_fileout(1:mcm_filelen(2))//extension
  !
  open(unit=22,file=stfilename,form='unformatted',access='direct', convert='little_endian',recl=512)
  mcm_d3count = 1
 endif
endif
!
! First write header info
!
d3data(1) = mcm_ptime        ! Problem time
d3data(2) = mcm_kinetice     ! Total KE
d3data(3) = mcm_internale    ! Total IE
d3data(4) = mcm_totale       ! Total energy
d3data(5) = 0.0_d            ! X vel (currently don't calculate)
d3data(6) = 0.0_d            ! Y vel (currently don't calculate)
d3data(7) = 0.0_d            ! Z vel (currently don't calculate)
!
call mcm_write_d3plot(d3data,7)
!
! Write material IE
!
do i=1,mcm_nummat
 d3data(1) = mcm_mat(i)%ie
 call mcm_write_d3plot(d3data,1)
enddo
!
! Write material KE
!
do i=1,mcm_nummat
 d3data(1) = mcm_mat(i)%ke
 call mcm_write_d3plot(d3data,1)
enddo
!
! Write material x vel
!
do i=1,mcm_nummat
 d3data(1) = mcm_mat(i)%mom(1)
 call mcm_write_d3plot(d3data,1)
enddo
!
! Write material y vel
!
do i=1,mcm_nummat
 d3data(1) = mcm_mat(i)%mom(2)
 call mcm_write_d3plot(d3data,1)
enddo
!
! Write material z vel
!
do i=1,mcm_nummat
 d3data(1) = mcm_mat(i)%mom(3)
 call mcm_write_d3plot(d3data,1)
enddo
!
! Write material mass
!
do i=1,mcm_nummat
 d3data(1) = mcm_mat(i)%mass
 call mcm_write_d3plot(d3data,1)
enddo
!
! Write unknowm data - not actually needed. may be my mistake but when looking at
!  the d3plot files produced by LSDYNA there seemed to be an extra piece of data.
!
!do i=1,mcm_nummat
! d3data(1) = 0.0_d
! call mcm_write_d3plot(d3data,1)
!enddo
!
! NODAL VALUES
!
! Write particle coordinates
!
do i=1,mcm_max_np
 ! convert to 32 bit
 do j=1,3
  x(j) = par(i)%x(j)
 enddo
 call mcm_write_d3plot(x,3)
enddo
!
! Write particle velocities
!
do i=1,mcm_max_np
 ! convert to 32 bit
 do j=1,3
  x(j) = par(i)%v(j)
 enddo
 call mcm_write_d3plot(x,3)
enddo
!
! Write particle accelerations
!
do i=1,mcm_max_np
 ! convert to 32 bit
 do j=1,3
  x(j) = par(i)%a(j)
 enddo
 call mcm_write_d3plot(x,3)
enddo
!
! ELEMENT VALUES
!
do i=1,mcm_max_np
 d3data(1) = par(i)%mat              ! material
 d3data(2) = par(i)%h                ! Smoothing length
 d3data(3) = par(i)%p                ! Pressure
 d3data(4) = par(i)%sigma(1,1)       ! Stress 1
 d3data(5) = par(i)%sigma(2,2)       ! Stress 2
 d3data(6) = par(i)%sigma(3,3)       ! Stress 3
 d3data(7) = par(i)%sigma(1,2)       ! Stress 4
 d3data(8) = par(i)%sigma(1,3)       ! Stress 5
 d3data(9) = par(i)%sigma(2,3)       ! Stress 6
 d3data(10) = par(i)%efps            ! Plastic strain
 d3data(11) = par(i)%rho             ! Density
 d3data(12) = par(i)%e               ! Internal energy
 d3data(13) = par(i)%nnbr            ! Number of neighbours
 d3data(14) = 0.0_d                  ! Strain 1 - not currently calculated
 d3data(15) = 0.0_d                  ! Strain 2 - not currently calculated
 d3data(16) = 0.0_d                  ! Strain 3 - not currently calculated
 d3data(17) = 0.0_d                  ! Strain 4 - not currently calculated
 d3data(18) = 0.0_d                  ! Strain 5 - not currently calculated
 d3data(19) = 0.0_d                  ! Strain 6 - not currently calculated
 !
 call mcm_write_d3plot(d3data,19)
enddo
!
! write end of file marker
!
!call mcm_write_d3plot(-999999.0,1)
!
!call mcm_shutdown(1)
!
1000 format('.ptf0',i1)
1010 format('.ptf',i2)
1020 format('.ptf',i3)
!
end subroutine mcm_d3plot_state
