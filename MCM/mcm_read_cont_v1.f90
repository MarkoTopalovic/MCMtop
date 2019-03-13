subroutine mcm_read_cont_v1
!************************************************************************
!
!    Purpose: reads and writes version 1 control card data
!
!  Called by: getinput
!
!       Date: 11-12-98
!
!     Errors: Improperly formatted input data
!
!      Notes: For input file version 1
!
!************************************************************************
!
use mcm_database
!
implicit none       
!
character (len=80):: mssg,txts
integer ::vern,lcount,mirr_opt
integer:: i,j
logical :: freeflag
character(len=16) :: filetest
logical :: fexist
!
!     set default values in case value are not specified on input
!
!	default contact type =0 kernel contact	
mcm_contacttype=0
!	default discretization type =0 colocated
mcm_disctype=0
!	default kernel type =1 b-spline
mcm_krtype=1
!	default initial time step size and time step scale factor
mcm_itss = 0.0_d  
! no symmetry planes
mcm_boundary = .false.
! MCMGUI is default output
mcm_state_opt = 0  
!
mcm_boundary_type = 1         ! Default to free-free boundaries
!
mcm_status_interval = 100
mcm_init_v_opt = 1
!
mcm_repulsive_force = .false. ! Switch to select whether repulse is included in simulations
!
mirr_opt = 0
!
mssg =' error reading 1st control card'
call mcm_gttxsg (txts,lcount)
read (unit=txts,fmt=100,err=400) mcm_axopt,mcm_disctype,mcm_nummat,mcm_np, mcm_nstressp, mcm_nvelocp
!
mssg =' error reading 2nd control card'
call mcm_gttxsg (txts,lcount)
read (unit=txts,fmt=110,err=400) mcm_endtime,mcm_tssf,mcm_itss
!
mssg =' error reading 3rd control card'
call mcm_gttxsg (txts,lcount)
read (unit=txts,fmt=120,err=400) mcm_stpltime,mcm_thpltime,mcm_thnode,mcm_num_transducer,mcm_state_opt
!
mssg =' error reading 4th control card'
call mcm_gttxsg (txts,lcount)
read (unit=txts,fmt=130,err=400) mcm_krtype,mcm_init_h_opt,mcm_maxnbr
!
mssg =' error reading 5th control card'
call mcm_gttxsg (txts,lcount)
read (unit=txts,fmt=140,err=400) mcm_contacttype, mcm_massopt, mcm_tcrit_opt, mcm_veloc_opt, &
                                 mcm_h_opt, mirr_opt
if(mirr_opt.eq.1) mcm_boundary = .true.
!
! Check if output files exist
!
! check that main state plot file does not exist
filetest=mcm_fileout(1:mcm_filelen(2))//'.s000'
inquire(file=filetest,exist=fexist)
if(fexist) then
 write(*,1000)
 write(*,170) filetest
 call mcm_shutdown(2)
endif
!check that log file does not exist, if not then open it
filetest=mcm_fileout(1:mcm_filelen(2))//'.log'
inquire(file=filetest,exist=fexist)
if(fexist) then
 write(*,1000)
 write(*,171) filetest
 call mcm_shutdown(2)
else
 open(unit=13,file=filetest,status='new',form='formatted')
 mcm_openfiles(13)=1  !record that file 13 is open
endif
!
call mcm_print_version(2)
write(13,165) mcm_title,1
!
! set variables to default values if required
!
if(mcm_tssf.eq.0.0) mcm_tssf=0.8_d
if(mcm_krtype.eq.0) mcm_krtype=1
if(mcm_maxnbr.eq.0) mcm_maxnbr=30
mcm_maxcont = mcm_maxnbr
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
! Set up loop variables for discretisation options
!
select case (mcm_disctype)
 case(0,1)
  ! colocated
  mcm_ssp=1
  mcm_svp=1
  mcm_esp=mcm_np
  mcm_evp=mcm_np
end select
!
! correct from old format to new
if(mcm_disctype.eq.5) then
 mcm_disctype = 1
else if(mcm_disctype.eq.3) then
 mcm_disctype = 2
else if(mcm_disctype.eq.1) then
 write(*,1300) 1
 write(13,1300) 1
 call mcm_shutdown(2)
else if(mcm_disctype.eq.2) then
 write(*,1300) 1
 write(13,1300) 1
 call mcm_shutdown(2)
endif
!
! write out control card data
!
write(13,200)
write(13,210) 1
write(13,300) mcm_axopt,mcm_disctype,mcm_nummat,mcm_np,mcm_ndim
!
write(13,210) 2
write(13,310) mcm_endtime,mcm_tssf,mcm_itss
!
write(13,210) 3
write(13,320) mcm_stpltime,mcm_thpltime,mcm_thnode,mcm_num_transducer
!
write(13,210) 4
write(13,330) mcm_krtype,mcm_init_h_opt,mcm_maxnbr
!
write(13,210) 5
write(13,340) mcm_contacttype, mcm_massopt, mcm_tcrit_opt, mcm_veloc_opt, mcm_h_opt
if(mcm_disctype.eq.2) write(13,350) mcm_nstressp,mcm_nvelocp
!
return
!
100 format(2i5,4i10)
110 format(3e10.0)
120 format(2e10.0,i10,i5,i5)
130 format(i10,2i5)
140 format(7i5)
!
165 format(//1x,70('*'),//1x,a78,//1x,70('*'),//10x,'input format version: ',i5,'')
!
170 format(/5x,'State file ',a16,' already exists.')
171 format(/5x,'Log file ',a16,' already exists.')
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
      4x,'number of dimensions........................    ',i7//)

310 format( &
      4x,'termination time............................ ',1pe10.2// &
      4x,'time step scale factor...................... ',1pe10.2// &
      4x,'initial time step size...................... ',1pe10.2/ &
     10x,'eq.0.0,  program picks initial step size        ',  //)
320 format( &
      4x,'time interval between state plots........... ',1pe10.2//, &
      4x,'time interval between time history plots.... ',1pe10.2// &
      4x,'number of time history nodes................ ',i10// &
	  4x,'Number of pressure transducers..............      ',i5//)
330 format( &
      4x,'kernel type.................................',i11// &
      4x,'smoothing length input option...............    ',i7// &
	  4x,'Maximum number of neighbours................    ',i7// )
340 format( &
      4x,'contact type................................    ',i7// &
      4x,'initial mass calculation option.............    ',i7// &
	  4x,'critical timestep calculation option........    ',i7// &
	  4x,'Velocity smoothing option...................    ',i7// &
	  4x,'Variable smoothing length option............    ',i7// )
350 format( &
      4x,'number of stress points.....................',i11// &
	  4x,'number of velocity points...................',i11//)
!
400 call mcm_termin (txts,mssg,lcount,1)
!
1000 format(/,'Error in subroutine read_cont')
1100 format(/,'Number of velocity points plus number of stress points does not equal',/,&
              'total number of particles.')
1200 format(/,'Error: 2D axisymmetry is not supported in version 2.')
1300 format(/,'Error: Discretisation type ',i1,' is not supported in version 2.')
!
end subroutine mcm_read_cont_v1
