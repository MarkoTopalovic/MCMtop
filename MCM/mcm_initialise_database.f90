subroutine mcm_initialise_database
!
!************************************************************************
!
!		Purpose: Initialise main database values
!
!	  Called by: mcm_allocate_memory
!
!	     Author: J. Campbell
!
!          Date: 29-01-2006
!
!        Errors: 
!
!         Notes: This routine has been added to explicitly initilise the 
!                particle database values (generally to zero). This is because
!                compilers are not consistent in the way they initialise
!                variables when they are forst defined.  It is safest to 
!                explicitly intitialise variables as part of the code, so that
!                there is no confusion.
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i
!
do i=1,mcm_max_np
 par(i)%active = .false.                            ! Particle deletion flag
 par(i)%mat = 0                                     ! Material model id
 par(i)%dispbc = 0                                  ! Displacement boundary condition
 par(i)%llpointer = 0                               ! linked list pointer
 par(i)%delpointer = 0                              ! inactive particle linked list pointer
 par(i)%nnbr = 0                                    ! number of neighbour particles
 par(i)%ncont = 0                                   ! number of neighbour contact particles
 !
 par(i)%g_nnbr = 0                                  ! number of ghost neighbours
 par(i)%g_ncont = 0                                 ! number of ghost contact neighbours
 par(i)%boundary = 0                                ! flag for boundary particles
 par(i)%nsym = 0                                    ! flag for symmetry planes
! par(i)%n_symnbr = 0                                ! number of symmetry neighbours
! par(i)%n_symcont = 0                               ! number of symmetry contact neighbours
 !
 par(i)%mass = 0.0_d                                ! Mass
 par(i)%h = 0.0_d                                   ! Smoothing length
 par(i)%hold = 0.0_d                                ! Smoothing length for previous timestep
 par(i)%h0 = 0.0_d                                  ! Initial smoothing length (fixed during initialisation)
 par(i)%rho = 0.0_d                                 ! Density
 par(i)%rho0 = 0.0_d                                ! Initial density (fixed during initialisation)
 par(i)%rhoold = 0.0_d
 par(i)%p = 0.0_d                                   ! pressure
 par(i)%e = 0.0_d                                   ! TOTAL particle internal energy
 par(i)%etry = 0.0_d                                ! trial particle internal energy
 par(i)%einc = 0.0_d                                ! particle internal energy increment
 par(i)%c = 0.0_d                                   ! speed of sound
 par(i)%vabs = 0.0_d                                ! Absolute value of velocity
 par(i)%tracerod = 0.0_d                            ! Trace of the rate-of-deformation tensor
 par(i)%pcut = 0.0_d                                ! Cutoff pressure
 par(i)%p_cut = 0.0_d                               ! Cutoff pressure at different time
 par(i)%efps = 0.0_d                                ! Effective plastic strain
 par(i)%fail = 0.0_d                                ! Failure, varies between undamaged (1.0) and failed (0.0)
 par(i)%temper = 0.0_d                              ! Temperature, unused at this time
 par(i)%mindist = 0.0_d                             ! Distance to nearest neighbour
 par(i)%critts = 0.0_d                              ! Particle critical timestep
 par(i)%epx = 0.0_d                                 ! material history variables
 par(i)%alfa = 0.0_d				             	! back stress (kinematic hardening)
 ! vectors
 par(i)%x = 0.0_d          		             	    ! Current particle coordinates
 par(i)%xzero = 0.0_d     		             	    ! Initial particle coordinates
 par(i)%v = 0.0_d         		             	    ! Current particle velocity
 par(i)%smooth_v = 0.0_d  		             	    ! Interpolated particle velocity, currently used only for XSPH
 par(i)%a = 0.0_d         		             	    ! Particle acceleration
 par(i)%bndnorm  = 0.0_d  		             	    ! Unit normal vector to surface for boundary particles
 par(i)%repulsion = 0.0_d		             		! contact acceleration
 ! tensors
 par(i)%rod = 0.0_d      		              	    ! Rate-of-deformation tensor
 par(i)%spin = 0.0_d    		              	    ! Spin tensor
 par(i)%sigma = 0.0_d   		             	    ! stress tensor
 par(i)%s = 0.0_d	    		             	    ! deviatoric stress tensor
 par(i)%q = 0.0_d       		             	    ! artificial viscosity stress tensor
 par(i)%qold = 0.0_d   		             	        ! artificial viscosity stress tensor from previous timestep
 par(i)%qq = 0.0_d       		             	    ! global to current configuration rotation matrix 
 par(i)%qr = 0.0_d       		             	    ! material axes rotation matrix for orthotropic materials
!
enddo
!
end subroutine mcm_initialise_database