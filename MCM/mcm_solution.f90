subroutine mcm_solution
!************************************************************************
!
!    Purpose: main solution loop. Contains the time integration loop
!
!  Called by: MAIN
!
!       Date: 09-08-2002
!
!     Errors: 
!
!      Notes: The objective is to keep this routine as concise as possible
!             by moving as much computation as possible into subroutines.
!
!             By this point in the code all variables which require calculation
!             before the time integration begins must have been calculated.
!
!             This subroutine uses the same central difference 
!             time integration algorithm used in DYNA3D
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i
!
write(*,900)
write(*,1000) mcm_timestep, mcm_ptime, mcm_dt
write(13,1000) mcm_timestep, mcm_ptime, mcm_dt
!
mcm_init_ts = mcm_dt
!
! Start of time integration loop
!
do
 !
 ! calculate new neighbour set and idemtify boundary particles if required
 !
 call mcm_neighbours
 !
 ! calculate rate of deformation tensor
 !
 call mcm_strain
 !
 ! update smoothing length (if necessary)
 !
 if(mcm_h_opt.gt.0) call mcm_update_h
 !
 ! material modelling
 !
 call mcm_constitutive
 !
 !!
 ! update time and set new time step
 !
 mcm_dtold = mcm_dt
 mcm_dt = mcm_tssf*mcm_critts
 mcm_dt = min(mcm_dt,1.1_d*mcm_dtold)
 !mcm_dt = max(mcm_dt,0.05_d*mcm_init_ts)
 mcm_ptime = mcm_ptime + mcm_dt
 !
 if (mcm_contacttype.gt.0) call mcm_rep_cont
 !
 ! solve momentum equation to calculate new acceleration
 !
 ! *** IMPORTANT *** : the particle acceleration vector is set to zero at the start of the momentum routines
 call mcm_momentum
 !
 ! impose acceleration boundary conditions
 !
 call mcm_boundcond
 !
 ! Calculate ke and ie by material > available for plot files
 !
 call mcm_calctotale
 !
 ! check whether output files should be written
 !
 ! state plot file
 if(mcm_ptime.ge.mcm_nextsttime.or.mcm_ptime.ge.mcm_endtime) then
  call mcm_state_output
  write(*,1200) mcm_ptime
  write(13,1200) mcm_ptime
  mcm_nextsttime = mcm_nextsttime + mcm_stpltime
 endif
 ! time history files
 if(mcm_ptime.ge.mcm_nextthtime.or.mcm_ptime.ge.mcm_endtime) then
  call mcm_timehist
  mcm_nextthtime = mcm_nextthtime + mcm_thpltime
 endif
 ! write problem status to log file
 if(mod(mcm_timestep,mcm_status_interval).eq.0) then
  write(*,1000) mcm_timestep, mcm_ptime, mcm_dt
  write(13,1000) mcm_timestep, mcm_ptime, mcm_dt
 endif
 !
 ! check for run termination conditions
 !
 if (mcm_ptime-mcm_dt.ge.mcm_endtime) exit
 !
 if(mcm_ctrlc) then
  call mcm_state_output
  write(*,1200) mcm_ptime
  write(13,1200) mcm_ptime
  exit
 endif
 !
 ! update velocity of particles
 !
 call mcm_update_velocity
 !
 ! update position of particles
 !
 call mcm_move
 !
 ! do any book-keeping required at end of timestep
 !
 ! the following routine writes CPU time and memory data to a status file
 ! the code used may not be portable to non Windows NT systems
 ! the code is contained only in the following subroutine and subroutine shutdown, 
 ! this subroutine can be commented out without affecting the rest of the program
 !call mcm_write_timings
 ! increment timestep counter
 mcm_timestep = mcm_timestep +1
 !
 ! Finally check whether to write a restart file
 ! Do this last so that can immediately start at the beginning
 ! of this loop with the information contained in the file.
 !
  !if(mcm_ptime.ge.mcm_endtime) then
  !  call mcm_write_restart
  !endif
!
! end of time integration loop
!
enddo
!
return
!
 900 format(/5x,'Start of solution')
1000 format(/5x,'Problem status for time-step ',i6, &
            /5x,'Time: ',e10.4,5x,'dt: ',e10.4)
1200 format(/5x,'State plot written at time ',e10.4)
!
end subroutine mcm_solution