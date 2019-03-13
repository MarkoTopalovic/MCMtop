subroutine mcm_initvbls
!************************************************************************
!
!    Purpose: initialise required variables
!
!  Called by: initial
!
!       Date: 02-08-2002
!
!     Errors: 
!
!      Notes:
!
!************************************************************************
!
use mcm_database
!
implicit none
INTEGER :: i,j
!
! Source/sink/adaptivity
!
do i=1,mcm_np
 par(i)%active = .true.    ! particle is active
 par(i)%delpointer = 0     ! initialise deleted particles linked list
enddo
!
! Store initial particle coordinates
!
do i=1,mcm_np
 do j=1,mcm_ndim
  par(i)%xzero(j) = par(i)%x(j)
 enddo
enddo
!
do i=1,mcm_np
  par(i)%a = 0.0_d
enddo
!
! SET EFFECTIVE PLASTIC STRAIN, TEMPERATURE AND THERMAL ENERGY TO ZERO
! FOR STRESS PARTICLES, CALCULATE TOTAL KINETIC ENERGY
!
mcm_thermale = 0.0_d
!
! Initialise stress point variables
! ssp = svp, esp = evp FOR THE CASE OF CO-LOCATED PARTICLES
!
DO i = mcm_ssp,mcm_esp
   par(i)%efps = 0.0_d
   par(i)%fail = 1.0_d
   par(i)%temper = 0.0_d  ! Not currently used
   !
   !	SET CUT OFF PRESURE
   !
   select case (mcm_mat(par(i)%mat)%model)
    case(1,3,4)
	 par(i)%pcut = -9.90E+20_d
	case(9)
     par(i)%pcut = mcm_mat(par(i)%mat)%strinput(1)
	case(10)
     par(i)%pcut = mcm_mat(par(i)%mat)%strinput(4)
   end select
ENDDO
!
! Initialise velocity point variables
!
do i = mcm_svp,mcm_evp
 ! absolute velocity, required for time step calculation options 1,2
 par(i)%vabs = 0.0_d
 do j=1,mcm_ndim
  par(i)%vabs = par(i)%vabs + par(i)%v(j)*par(i)%v(j)
 enddo
 !
 par(i)%vabs = sqrt(par(i)%vabs)
enddo
!
! Find maximum and minimum values of the coordinates. These are used for the
! linked list neighbour searching to set up the underlying grid.
! All particles are looped over as not concerned what type the particle is.
!
DO i=1,3
 mcm_coord_maxmin(1,i) =  1.0e+20_d
 mcm_coord_maxmin(2,i) = -1.0e+20_d
ENDDO
!
DO i=1,mcm_np
 DO j=1,mcm_ndim
  if(par(i)%x(j).lt.mcm_coord_maxmin(1,j)) mcm_coord_maxmin(1,j) = par(i)%x(j)
  if(par(i)%x(j).gt.mcm_coord_maxmin(2,j)) mcm_coord_maxmin(2,j) = par(i)%x(j)
 ENDDO
ENDDO
!
!
! Initialise unused particles and set mat = 0 so they are not plotted
!
if(mcm_max_np.gt.mcm_np) then
 do i=mcm_np+1,mcm_max_np
  par(i)%mat = 0
  par(i)%active = .false.
  par(i)%delpointer = i+1
  par(i)%x = par(1)%x
  par(i)%h = par(1)%h
 enddo
 par(mcm_max_np)%delpointer = 0
endif
end subroutine mcm_initvbls