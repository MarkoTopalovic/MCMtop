subroutine mcm_ctrlc_action
!
! This subroutine is called when control c is pressed
!
! For now all this does is set mcm_ctrlc to true.  This is checked in the main loop
!  and allows mcm to be terminated neatly.
!
use mcm_database
!
implicit none
!
mcm_ctrlc = .true.
!
write(*,1000)
write(13,1000)
!
1000 format(' Control-C interrupt - mcm will write a plot file and terminate.')
!
end subroutine mcm_ctrlc_action