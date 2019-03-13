subroutine mcm_init_d3plot
!************************************************************************
!
!    Purpose: Set up and write first section of d3plot file
!
!  Called by: init_state_output
!
!       Date: 6-4-2005
!
!     Errors: 
!
!      Notes: This routine writes LS-DYNA format d3plot files
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,j,k, nglbv, narbs
logical :: fexist
real(4) :: data_type
character(len=40) :: title
character(len=16) :: stfilename
!
! file is written in 32 bit words
real(4), dimension(64) :: control_words
real(4), dimension(9) :: isphfg
real(4), dimension(3) :: x
!
mcm_d3buffcnt = 1
mcm_d3file = 0
mcm_d3count = 1
!
! Open first file
!
!
stfilename=mcm_fileout(1:mcm_filelen(2))//'.ptf'
!
! delete old file if it exists then open new state plot file
!
!inquire(file=stfilename,exist=fexist)
!if(fexist) then
! open(unit=22,file=stfilename)
! close(unit=22,status='delete')
!endif
open(unit=22,file=stfilename,form='unformatted',access='direct', convert='little_endian',recl=512)
!
!
! Place variables in contol words
!
! Zero array
do i=1,64
 control_words(i) = 0.0_d
enddo
! Following is for big endian, get the title correct
!do i=1,10
! do j=1,4
!  k=(i-1)*4
!  title(k+j:k+j) = mcm_title(k+5-j:k+5-j)
! enddo
!enddo
!
!do i=1,10
! k=(i-1)*4+1
! control_words(i) = transfer(title(k:k+3),data_type)
!enddo
!
title = mcm_title(1:40)
nglbv = 6+6*mcm_nummat
do i=1,10
 k=(i-1)*4+1
 control_words(i) = transfer(title(k:k+3),data_type)
enddo
!
control_words(15) = 970.0                            ! version
control_words(16) = transfer(4,data_type)            ! NDIM
control_words(17) = transfer(mcm_max_np,data_type)   ! Number of nodes
control_words(18) = transfer(6,data_type)            ! ICODE (pretend its LSDYNA)
control_words(19) = transfer(nglbv,data_type)        ! NGLBV
control_words(20) = transfer(0,data_type)            ! IT
control_words(21) = transfer(1,data_type)            ! IU
control_words(22) = transfer(1,data_type)            ! IV
control_words(23) = transfer(1,data_type)            ! IA
control_words(24) = transfer(0,data_type)            ! NEL8 - no. hex elements
control_words(25) = transfer(0,data_type)            ! NUMMAT8 - no. materials used by hex elements
control_words(26) = transfer(0,data_type)            !  blank
control_words(27) = transfer(0,data_type)            !  blank
control_words(28) = transfer(7,data_type)            ! NV3D - no. vbls. per hex element
control_words(29) = transfer(0,data_type)            ! NEL2 - no. 1d elements
control_words(30) = transfer(0,data_type)            ! NUMMAT2
control_words(31) = transfer(6,data_type)            ! NV2D
control_words(32) = transfer(0,data_type)            ! NEL4 - no. 4 node 2d elements
control_words(33) = transfer(0,data_type)            ! NUMMAT4
control_words(34) = transfer(33,data_type)           ! NV2D
control_words(35) = transfer(0,data_type)            ! NEIPH
control_words(36) = transfer(0,data_type)            ! NEIPS
control_words(37) = transfer(-10003,data_type)       ! Maxint
control_words(38) = transfer(mcm_max_np,data_type)   ! Number of SPH nodes
control_words(39) = transfer(mcm_nummat,data_type)   ! Number of SPH materials
! narbs = 13 + mcm_np ! not sure if this is needed - works OK with 0
control_words(40) = transfer(0,data_type)            ! NARBS - use sequential numbering
control_words(41) = transfer(0,data_type)            ! NELT - no tshell elements
control_words(42) = transfer(0,data_type)            ! NUMMATT
control_words(43) = transfer(21,data_type)           ! NV3DT
control_words(44) = transfer(1000,data_type)         ! IOSHL(1)
control_words(45) = transfer(1000,data_type)         ! IOSHL(2)
control_words(46) = transfer(1000,data_type)         ! IOSHL(3)
control_words(47) = transfer(1000,data_type)         ! IOSHL(4)
! All further values are 0
!
! Write values to file
!
call mcm_write_d3plot(control_words,64)
!
! SPH data flags
!
isphfg(1) = transfer(9,data_type)      ! length of flags array
isphfg(2) = transfer(1,data_type)      ! Radius of influence
isphfg(3) = transfer(1,data_type)      ! Pressure
isphfg(4) = transfer(6,data_type)      ! 6 true stress components
isphfg(5) = transfer(1,data_type)      ! plastic strain
isphfg(6) = transfer(1,data_type)      ! density
isphfg(7) = transfer(1,data_type)      ! internal energy
isphfg(8) = transfer(1,data_type)      ! number of neighbours
isphfg(9) = transfer(6,data_type)      ! 6 true strain components
!
call mcm_write_d3plot(isphfg,9)
!
! Now write coordinates
!
do i=1,mcm_max_np
 ! convert to 32 bit
 do j=1,3
  x(j) = par(i)%x(j)
 enddo
 call mcm_write_d3plot(x,3)
enddo
!
! Write SPH node and its material number
!
do i=1,mcm_max_np
 x(1) = transfer(i,data_type)
 x(2) = transfer(par(i)%mat,data_type)
 !
 call mcm_write_d3plot(x,2)
enddo
!
!call mcm_shutdown(1)
end subroutine mcm_init_d3plot
