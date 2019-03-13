      subroutine mcm_getinput
!************************************************************************
!
!    Purpose: Control routine for reading test input file
!
!  Called by: MAIN
!
!       Date: 30-07-2002
!
!     Errors: Wrong input file version number
!
!      Notes: Only the title and the input file version number are
!             read in this routine, all other data is read in subroutines
!             called by this routine.
!             Supports versions 1 and 2 of input file.
!
!************************************************************************
!
use mcm_database
!
implicit none       
!
character(len=16) :: filename
character (len=80):: mssg,txts
integer :: vern,lcount, i, j
!
! Set default boundary plane values
!
mcm_boundary_type = 1  ! array of dimension (3)
mcm_boundary_code = 0  ! array of dimension (2,3)
!
! open input file
!
filename=mcm_filein(1:mcm_filelen(1))//'.mcm'
open(unit=101,file=filename,status='old')
mcm_openfiles(1)=1
!
!     read and write  title		problem title, 
!					  vern		input version
!
call mcm_gttxsg (txts,lcount)
read (unit=txts,fmt=160) mcm_title,vern
!
! choose version number of input file, currently only version 1
! is defined, but this adds the support for future changes
!
select case (vern)
 !
 !=============================================================
 ! Version 1. Original input file version
 case(1) 
  !
  !  read control card data
  !
  call mcm_read_cont_v1
  !
  !  Allocate memory for particle arrays
  !
  call mcm_allocate_memory
  call mcm_initialise_database
  !
  !  read and write material properties
  !
  call mcm_matin 
  !
  !  read and write particle data
  !  routine nodein will read colocated and non-colocated data
  !
  call mcm_nodein
  !
  !  read and write particles for time histories
  !
  if(mcm_thnode.gt.0) call mcm_nodehist
  !
  ! real and write coordinates for pressure transducers
  !
  if(mcm_num_transducer.gt.0) call mcm_transducers_in
  !
  !  read and write initial velocity conditions
  !
  call mcm_inputvel
  !
  ! read symmetry planes
  !
  if(mcm_boundary) call mcm_inp_symmetry_v1
  !
  !=================================================================
  ! Format 2. Revised input format. Differs from version 1 in 
  !           main control cards and symmetry planes.
 case(2) 
  !
  !  read control card data
  !
  call mcm_control_v2
  !
  !  Allocate memory for particle arrays
  !
  call mcm_allocate_memory
  call mcm_initialise_database
  !
  !  read and write material properties
  !
  call mcm_matin 
  !
  !  read and write particle data
  !  routine nodein will read colocated and non-colocated data
  !
  call mcm_nodein
  !
  !  read and write particles for time histories
  !
  if(mcm_thnode.gt.0) call mcm_nodehist
  !
  ! real and write coordinates for pressure transducers
  !
  if(mcm_num_transducer.gt.0) call mcm_transducers_in
  !
  !  read and write initial velocity conditions
  !
  call mcm_inputvel
  !
  ! Boundary planes
  !
  if(mcm_boundary) call mcm_inp_symmetry_v2
  !
  if(mcm_contacttype.gt.0) call mcm_contact_in
  !
  if(mcm_baseaccel) call mcm_baseaccel_in
  !
 case default
  !
  ! error termination if input file version is not supported
  write(*,1000)
  write(*,1100) vern
  call mcm_shutdown(2)
end select
!
! close the input file 
!
close(unit=101,status='keep')
mcm_openfiles(1)=0
!
return
!
160 format(a78,i2)
221 format( &
     4x,'viscosity reset option.........................',i7/ &
     10x,'eq.0:  default viscosities set by dyna3d',/ &
     10x,'eq.1:  default viscosities read in',// &
     4x,'time step scale factor.........................',1pe10.2/)
229 format(&
     4x,'initial min time step for mass augmentation ...',1pe10.3// &
     4x,'sustained min time step for mass augmentaion ..',1pe10.3//)
230 format(&
     4x,'ratio of current to init energy for abort .....',1pe10.3/)
320 format(3i5)
!
!     error termination due to badly formatted data
!
400 call termin (txts,mssg,lcount,1)
!
1000 format(/,'Error in subroutine getinput')
1100 format('Input file version ',i2,' is not supported in this version of the code')
!
end subroutine mcm_getinput
