subroutine mcm_ctrlc_handle
!
! Catch control-c interrupt to ensure clean termination
!
! See dyna routines ctrlco.f and enablc.f to see how dyna handles ctrl-c
!
use mcm_database
use ifport
!
implicit none
!
external mcm_ctrlc_action
!
integer :: inum, signum
!
signum = 2  ! ctrl-c interrupt
!
inum = signal(signum,mcm_ctrlc_action,-1)
!
! Set logical value that is checked in main routines
!
mcm_ctrlc = .false.
!
end subroutine mcm_ctrlc_handle