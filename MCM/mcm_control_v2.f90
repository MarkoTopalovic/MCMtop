subroutine mcm_control_v2
!************************************************************************
!
!    Purpose: reads and writes version 2 control card data
!
!  Called by: getinput
!
!       Date: 01-08-2002
!
!     Errors: Improperly formatted input data
!
!      Notes: For input file version 2
!
!************************************************************************
!
use mcm_database
!
implicit none       
!
character (len=80):: mssg,txts
integer ::vern,lcount,mirr_opt, keep_files
integer:: i,j
logical :: freeflag
character(len=16) :: filetest
logical :: fexist
!
!     set default values in case value are not specified on input
!
mcm_disctype = 0              ! default discretization type =0: basic SPH
!
mcm_max_np = 0                ! initialise max number of particles
!
mcm_itss = 0.0_d              ! default initial time step size
!
mcm_state_opt = 0             ! MCMGUI is default output
!
mcm_thnode = 0                ! default no time history particles
!
mcm_num_transducer = 0        ! default no pressure transducers
!
mcm_status_interval = 100     ! default 100 steps between status reports to screen
!
mcm_restart_interval = 0      ! default no restart fie written
!
mcm_run_restart = 0           ! default no running restart file
!
keep_files = 0                ! default do no overwrite output files
!
mcm_veloc_opt = 0             ! default zero initial velocity
!
mcm_massopt = 0               ! default mass given for material
!
mcm_init_h_opt = 0            ! default smoothing length given for material
!
mcm_init_rhoe = 0             ! default initial density and energy from material cards
!	
mcm_contacttype = 0           ! default contact type = 0 kernel contact
!
mcm_tcrit_opt = 0             ! default crtitical timestep calculation
!
mcm_veloc_opt = 0             ! default particles move with own velocity
!
mcm_boundary = .false.        ! no boundary planes
!
mirr_opt = 0
!
mcm_nlcur = 0                 ! number of  load curves
!
mcm_h_opt = 0                 ! default constant smoothing length
!
mcm_krtype = 1                ! default kernel type =1 b-spline
!
mcm_maxnbr = 40               ! default maximum number of neighbours
!
mcm_nthpx = 0                 ! no base acceleration in X-direction
!
mcm_nthpy = 0                 ! no base acceleration in Y-direction 
!
mcm_nthpz = 0                 ! no base acceleration in Z-direction  
!
mcm_drelax = .false.          ! Default standard analysis
!
!========================================================================================
!
! Read control cards
!
! Control card 1: problem definition
mssg =' error reading 1st control card'
call mcm_gttxsg (txts,lcount)
read (unit=txts,fmt=100,err=400) mcm_axopt,mcm_disctype,mcm_nummat,mcm_np, mcm_max_np, &
                                 mcm_nstressp, mcm_nvelocp
!
! Control card 2: time control
mssg =' error reading 2nd control card'
call mcm_gttxsg (txts,lcount)
read (unit=txts,fmt=110,err=400) mcm_endtime,mcm_tssf,mcm_itss,mcm_drelax_scale
!
! Control card 3: output file control
mssg =' error reading 3rd control card'
call mcm_gttxsg (txts,lcount)
read (unit=txts,fmt=120,err=400) mcm_stpltime,mcm_state_opt,mcm_thpltime,mcm_thnode,mcm_num_transducer, &
                                 mcm_status_interval, mcm_restart_interval, mcm_run_restart
!
! Control card 4: input and initialization options
mssg =' error reading 4th control card'
call mcm_gttxsg (txts,lcount)
read (unit=txts,fmt=130,err=400) keep_files,mcm_init_v_opt,mcm_massopt,mcm_init_h_opt,mcm_init_rhoe
!
! Control card 5: analysis options
mssg =' error reading 5th control card'
call mcm_gttxsg (txts,lcount)
read (unit=txts,fmt=140,err=400) mcm_contacttype, mcm_tcrit_opt, mcm_veloc_opt, mirr_opt, mcm_nlcur, &
                                 mcm_nthpx, mcm_nthpy, mcm_nthpz 
!
! Control card 6: Interpolation options
mssg =' error reading 6th control card'
call mcm_gttxsg (txts,lcount)
read (unit=txts,fmt=150,err=400) mcm_h_opt, mcm_krtype,mcm_maxnbr
!
! Control card 7: Options
mssg =' error reading 7th control card'
call mcm_gttxsg (txts,lcount)
! blank at this time
!read (unit=txts,fmt=160,err=400) 
!
!========================================================================================
!
if(keep_files.ne.1) then
 !
 ! Stop if output files exist
 !
 if(mcm_state_opt.eq.1) then
  ! check that ensight geometry file does not exist
  filetest=mcm_fileout(1:mcm_filelen(2))//'.geo'
  inquire(file=filetest,exist=fexist)
  if(fexist) then
   write(*,1000)
   write(*,171) filetest
   call mcm_shutdown(2)
  endif
  !
 elseif (mcm_state_opt.eq.2) then
  ! check that CASE file does not exist
  filetest=mcm_fileout(1:mcm_filelen(2))//'.case'
  inquire(file=filetest,exist=fexist)
  if(fexist) then
   write(*,1000)
   write(*,171) filetest
   call mcm_shutdown(2)
  endif
  !
 elseif (mcm_state_opt.eq.3) then
  ! check that d3plot file does not exist
  filetest=mcm_fileout(1:mcm_filelen(2))//'.ptf'
  inquire(file=filetest,exist=fexist)
  if(fexist) then
   write(*,1000)
   write(*,170) filetest
   call mcm_shutdown(2)
  endif
  !
 else
  ! check that main state plot file does not exist
  filetest=mcm_fileout(1:mcm_filelen(2))//'.s000'
  inquire(file=filetest,exist=fexist)
  if(fexist) then
   write(*,1000)
   write(*,170) filetest
   call mcm_shutdown(2)
  endif
 endif
endif
!
! open log file
filetest=mcm_fileout(1:mcm_filelen(2))//'.log'
inquire(file=filetest,exist=fexist)
if(fexist) then
 ! delete old log file
 open(unit=13,file=filetest,status='old',form='formatted')
 close(unit=13,status='delete')
endif
open(unit=13,file=filetest,status='new',form='formatted')
mcm_openfiles(13)=1  !record that file 13 is open
!
! Write start of log file
call mcm_print_version(2)
write(13,165) mcm_title,2
!
!==========================================================================================
!
! set variables to default values if required
!
if(mcm_tssf.eq.0.0_d) mcm_tssf=0.8_d
if(mcm_krtype.eq.0) mcm_krtype=1
if(mcm_maxnbr.eq.0) mcm_maxnbr=40
if(mcm_status_interval.lt.1) mcm_status_interval = 100
if(mirr_opt.eq.1) mcm_boundary = .true.
mcm_maxcont = mcm_maxnbr
mcm_g_maxnbr = mcm_maxnbr/2
mcm_g_maxcont = mcm_g_maxnbr
!
!
! Check for maximum number of particles
!
if(mcm_max_np.lt.mcm_np) mcm_max_np = mcm_np
!
! Process axis options
!
select case (mcm_axopt)
 case(1,2,3)
  ! 1,2,3D Cartesian
  mcm_ndim=mcm_axopt
 case(4)
  ! 2D Axisymmetric
  write(*,1200) 1
  write(13,1200) 1
  call mcm_shutdown(2)
end select
!
!
! Dynamic realxation
!
if(mcm_drelax_scale.gt.0.1_d.and.mcm_drelax_scale.lt.1.0_d) then
 mcm_drelax = .true.
endif
!
! Set up loop variables for discretisation options
!
select case (mcm_disctype)
 case(0,1)
  ! colocated
  mcm_ssp=1
  mcm_svp=1
  mcm_esp=mcm_np
  mcm_evp=mcm_np
  !
end select
!
! Check for base accelerations
!
if(mcm_nthpx.eq.1.or.mcm_nthpy.eq.1.or.mcm_nthpz.eq.1) then
 mcm_baseaccel = .true.
endif
!
! write out control card data
!
write(13,200)
write(13,210) 1
write(13,300) mcm_axopt,mcm_disctype,mcm_nummat,mcm_np,mcm_max_np,mcm_ndim
!
write(13,210) 2
write(13,310) mcm_endtime,mcm_tssf,mcm_itss
!
if(mcm_drelax) then
 write(13,311) mcm_drelax_scale
else
 write(13,312)
endif
!
write(13,210) 3
write(13,320) mcm_stpltime,mcm_state_opt,mcm_thpltime,mcm_thnode,mcm_num_transducer, &
                                 mcm_status_interval, mcm_restart_interval, mcm_run_restart
!
write(13,210) 4
write(13,330) keep_files,mcm_init_v_opt,mcm_massopt,mcm_init_h_opt,mcm_init_rhoe
!
write(13,210) 5
write(13,340) mcm_contacttype, mcm_tcrit_opt, mcm_veloc_opt, mirr_opt, mcm_nlcur, &
                                  mcm_nthpx, mcm_nthpy, mcm_nthpz
! 
write(13,210) 6
write(13,360) mcm_h_opt, mcm_krtype,mcm_maxnbr
!
write(13,210) 7
write(13,370) 
!
return
!
100 format(2i5,5i10)
110 format(4e10.0)
120 format(e10.0,i5,e10.0,i10,i5,i10,2i5)
130 format(5i5)
140 format(10i5)
150 format(3i5)
160 format(i5,i5,e10.0)
!
165 format(//1x,70('*'),//1x,a78,//1x,70('*'),//10x,'input format version: ',i5,'')
!
170 format(/5x,'State file ',a16,' already exists.')
171 format(/5x,'Geometry file ',a16,' already exists.')
!
200 format(//' c o n t r o l   i n f o r m a t i o n')
210 format(//, &
      4x,'*******************************************************',/ &
      4x,'*                                                     *',/ &
      4x,'*                   Control Card ',i2,19x,'*',/ &
      4x,'*                                                     *',/ &
      4x,'*******************************************************',/)
!
300 format( &
      4x,'axis option.................................    ',i7// &
      4x,'discretisation type.........................    ',i7// &
      4x,'number of materials.........................',i11// &
      4x,'number of nodes.............................',i11// &
	  4x,'maximum number of nodes.....................',i11// &
      4x,'number of dimensions........................    ',i7//)

310 format( &
      4x,'termination time............................ ',1pe10.2// &
      4x,'time step scale factor...................... ',1pe10.2// &
      4x,'initial time step size...................... ',1pe10.2/ &
     10x,'eq.0.0,  program picks initial step size        ',  //)
311 format( &
      4x,'Dynamic relaxation active with scale factor  ',1pe10.2//)
312 format(4x,'Dynamic relaxation not active.'//)
320 format( &
      4x,'time interval between state plots........... ',1pe10.2//, &
	  4x,'output file fomat........................... ',i10// &
      4x,'time interval between time history plots.... ',1pe10.2// &
      4x,'number of time history nodes................ ',i10// &
	  4x,'Number of pressure transducers..............      ',i5// &
	  4x,'Number of steps between status reports...... ',i10// &
	  4x,'Number of steps between restart files....... ',i10// &
	  4x,'Number of steps between running restart files',i10//)
330 format( &
      4x,'Preserve output files flag.................. ',i10// &
	 10x,'eq. 1, output files will be overwritten      ',//  &
	  4x,'Prescribed initial velocity flag............ ',i10// &
	  4x,'Mass initialisation option.................. ',i10// &
	  4x,'Smoothing length initialisation option...... ',i10// &
	  4x,'Particle density and energy initialisation   ',i10// )
340 format( &
      4x,'contact type................................    ',i7// &
	  4x,'critical timestep calculation option........    ',i7// &
	  4x,'Velocity smoothing option...................    ',i7// &
      4x,'Symmetry planes flag........................    ',i7// &
	  4x,'Number of load curves.......................    ',i7// &
      4x,'X-dir base acceleration.....................    ',i7// &
     10x,'eq.0,  no                                       ',  // &
     10x,'eq.1,  yes                                      ',  // &
      4x,'Y-dir base acceleration.....................    ',i7// &
     10x,'eq.0,  no                                       ',  // &
     10x,'eq.1,  yes                                      ',  // &
      4x,'Z-dir base acceleration.....................    ',i7// &
     10x,'eq.0,  no                                       ',  // &
     10x,'eq.1,  yes                                      ',  //)
350 format( &
      4x,'number of stress points.....................',i11// &
	  4x,'number of velocity points...................',i11//)
360 format( &
      4x,'Variable smoothing length option............    ',i7// &
      4x,'kernel type.................................',i11// &
	  4x,'Maximum number of neighbours................    ',i7// )
370 format(4x,'This card is not used in this version.')
!
400 call mcm_termin(txts,mssg,lcount,1)
!
1000 format(/,'Error in subroutine read_cont')
1100 format(/,'Number of velocity points plus number of stress points does not equal',/,&
              'total number of particles.')
1200 format(/,'Error: 2D axisymmetry is not supported in version 2.')
1300 format(/,'Error: Discretisation type ',i1,' is not supported in version 2.')
!
end subroutine mcm_control_v2
