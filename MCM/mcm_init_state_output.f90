subroutine mcm_init_state_output
!************************************************************************
!
!    Purpose: Control routine for initialising state plots 
!
!  Called by: initial, solution
!
!       Date: 05-08-2002
!
!     Errors: 
!
!      Notes: Routine allows coice of writing MCMGUI files
!             or Ensight files.
!			  Modified with option eq.2 which writes output in Ensight 
!             CASE format.  This reduces the amount of output files.  29-10-02 - tom
!
!************************************************************************
!
use mcm_database
!
implicit none
!
select case (mcm_state_opt)
 case(1)
  ! Ensight 5.0
  mcm_istate = 1 ! bug - 0 instead of 1?? - tom - 31-10-02
  !
  ! calculate time for second state plot
  !
  mcm_nextsttime = mcm_ptime + mcm_stpltime
  !
  call mcm_geometry
  !
 case(2)
  ! Ensight CASE
  mcm_istate = 0
  !
  ! calculate time for second state plot
  !
  mcm_nextsttime = mcm_ptime + mcm_stpltime
  !
  call mcm_write_geo
  !
 case(3)
  ! LSDYNA d3plot format
  !
  call mcm_init_d3plot
  !
  ! calculate time for second state plot
  !
  mcm_nextsttime = mcm_ptime + mcm_stpltime
  !
  !
  case(4)
  !
  ! Topalovic VTK
  !
  call WriteVtk
    call WriteVtkmat(1)
  call WriteVtkmat(2)
  ! calculate time for second state plot
  !
  mcm_nextsttime = mcm_ptime + mcm_stpltime
  !
  
 case default
  ! write MCMGUI output files
  call mcm_initstateplot
  !
end select
!
end subroutine mcm_init_state_output