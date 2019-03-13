subroutine mcm_initial
!************************************************************************
!
!        Purpose: Initialise all data
!
!      Called by: main
!
!         Author: Tom De Vuyst
!
!           Date: 01-08-2002
!
! Last Modified : 01-08-2002
!
!         Errors: 
!
!          Notes: Initialises stress, ....
!		    	  Calculates initial timestep, particle masses
!
!************************************************************************
!
use mcm_database
!
implicit none
integer :: i, l
!
! Write status report to screen and log file
!
write(*,100)
write(13,100)
!
! Apply displacement boundary conditions to initial velocity
!
call mcm_init_boundcond
!
! CALCULATE PARTICLE MASSES
!
CALL mcm_calcmass
!
! INITIALISE DATABASE AND HISTORY VARIABLES
!
CALL mcm_initvbls
!
! CALCULATE SHEAR MODULUS
!
CALL mcm_calcg
!
! CALCULATE DENSITY, PRESSURE AND SPEED OF SOUND 
!
CALL mcm_calcrhopc
!
! CALCULATE STRESS TENSOR, INITIALISE RATE OF DEFORMATION TENSOR 
!
CALL mcm_initsigma
CALL mcm_initrod
!
! ASSIGN INITIAL SMOOTHING LENGTHS
!
call mcm_init_h
!
! INITIAL TIMESTEP CALCULATION AND ART. VISCOSITY CALCULATION
!
mcm_timestep = 1
do i=mcm_ssp,mcm_esp
 par(i)%tracerod = 0.0_d
 do l = 1,mcm_ndim
   par(i)%tracerod = par(i)%tracerod + par(i)%rod(l,l)
 enddo
 CALL mcm_aviscosity(i)
enddo
!
!
! Get neighbours and boundary particles
!
CALL mcm_init_neighb
!
! Boundary particles are not currently used, but retain routines in case
!Call mcm_id_boundary    
!
! CALCULATE CRITICAL TIMESTEP
!
CALL mcm_calccritts
!
! SET FIRST TIMESTEP SIZE dt
!
if (mcm_itss .ne. 0.0_d) THEN
    mcm_dt = min(mcm_itss, mcm_critts*mcm_tssf)
	if (mcm_itss .gt. mcm_critts*mcm_tssf) THEN
	    write(*,2000)  mcm_itss, mcm_critts, mcm_tssf
		write(13,2000) mcm_itss, mcm_critts, mcm_tssf 
	endif
else
	mcm_dt=mcm_critts*mcm_tssf
endif
!
! CALCULATE TOTAL ENERGY
!
CALL mcm_calctotale
!
! INITIALISE rhoold(np) and qold(3,3,np)  
!
CALL mcm_initold
!
!
mcm_next_restart = mcm_restart_interval
mcm_next_run_restart = mcm_run_restart
!
! Initialise state output variables and plot first state output file
!
call mcm_init_state_output
call mcm_state_output
!
! Initialise time history output files and write first data line
!
call mcm_inittimehist
call mcm_timehist
!
! write status report to screen and log file
!
write(*,110)
write(13,110)
!
return
!
100 FORMAT(//5x,'Problem Initialisation Begun')
110 FORMAT(//5x,'Problem Initialisation Complete')
!
1000 FORMAT(//5X, 'Error in subroutine initial')
1010 FORMAT(/5X,' Material number ',I3,' does not exist.')
2000 format(//5X, 'WARNING', //5X, 'Initial timestep size ', E10.4, ' ignored. Greater &
       &than ', /5X,'Crit. Timestep Size x Factor of Safety' &
	   &/5X, '( ',E10.4, ' x ', E10.4, ' )')

end subroutine mcm_initial