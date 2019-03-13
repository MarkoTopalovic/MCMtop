program mcm
!***************************************************************************
!                                                                          *
!     mcm                                                                  *
!     A 1/2/3D Meshless Continuum Mechanics Code with Material Strength    *
!                                                                          *
!     Authors: J. Campbell                                                 *
!              T. De Vuyst                                                 *
!              R. Vignjevic                                                *
!              J. Reveles                                                  *
!                                                                          *
!     Version: 2.07                                                        *
!                                                                          *
!     Date of last modification: 11-05-2007                                *
!                                                                          *
!     School of Engineering                                                *
!     Cranfield University                                                 *
!     Cranfield                                                            *
!     Beds. MK43 0AL                                                       *
!                                                                          *
!     Copyright 1998-2005 Cranfield University                             *
!     All rights reserved.                                                 *
!                                                                          *
!     This program is intended as a development code for research on       *
!     meshless methods within the School of Engineering.                   *
!                                                                          *
!***************************************************************************
!                                                                          *
!     Version History                                                      *
!                                                                          *
!       1.00: Original MCM                                                 *
!                                                                          *
!       2.00: Rewrite to use Fortran90 data type for particle data, to     *
!             rename variables to prevent clashes with DYNA3D variable     *
!             names, to make symmetry planes easier to use and to tidy     *
!             up some routines.                                            *
!                                                                          *
!       2.01: Minor bugfixes and addition of contact.                      *                                               *
!                                                                          *
!       2.03: Addition of support for LSDYNA d3plot output format          *
!                                                                          *
!       2.04: Added base accelerations.                                    *
!                                                                          *
!       2.05: Added a basic dynamic relaxation capability.                 *
!                                                                          *
!       2.07: Rewrite of symmetry plane routines to allow symmetry and     *
!             periodic boundaries plus problem limit planes.               *
!                                                                          *
!                                                                          *
!***************************************************************************
!
! the following line must be used in ALL subroutines
!
implicit none
!
logical newproblem
!
! Set up trap for control c
!
call mcm_ctrlc_handle
!
! get input file names
!
call mcm_startup(newproblem)
!
if(newproblem) then
 !
 ! read input file
 !
 call mcm_getinput
 !
 ! problem initialisation
 !
 call mcm_initial
 !
else
 !
 ! read restart file
 !
 !call mcm_restart
endif
!
! main solution loop
!
call mcm_solution
!
! Write dynamic relaxation data
!
call mcm_write_drelax
!
! shutdown program
!
call mcm_shutdown(1)
!
end