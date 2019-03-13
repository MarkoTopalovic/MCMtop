SUBROUTINE mcm_calccritts
!************************************************************************
!
!		Purpose: calculate critical timestep
!
!	  Called by: Initial, COnstitutive
!
!	     Author: Tom De Vuyst
!
!          Date: 05-08-2002
!
! Last Modified: 05-08-2002 by J. Campbell
!
!        Errors: 
!
!         Notes: veae : Ve/Ae (Element Volume/Max. Element Area)
!                       as in DYNA3D's timestep calculation
!                qu   : same as Q in DYNA3D solid element timestep calc.
!                ae   : Ae; Max. Element Area
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i, j
REAL(kind=real_acc) :: veae, qu, trrod, ae
!
! INITIALISE TIMESTEP
!
mcm_critts = 1.0e20_d
!
! CALCULATE CRITICAL TIMESTEP PER PARCTICLE
!
DO i = mcm_ssp,mcm_esp
 !
 ! use specified form of critical timestep calculation
 !
 select case (mcm_tcrit_opt)
  case(0)
   !
   ! DYNA3D formula
   !
   ! CALCULATE Ve/Ae DEPENDING ON DIMENSION OF PROBLEM
   !
   SELECT CASE (mcm_ndim)

      CASE (1)
         veae = par(i)%mass/par(i)%rho
	     ae   = 1
      CASE (2)
         veae = SQRT(par(i)%mass/par(i)%rho)
	     ae   = SQRT(par(i)%mass/par(i)%rho)
      CASE (3)
         veae = (par(i)%mass/par(i)%rho)**(1.0_d/3.0_d)
	     ae   = (par(i)%mass/par(i)%rho)**(2.0_d/3.0_d)

   END SELECT
   !
   ! CALCULATE qu & particle timestep
   !
   qu = (mcm_mat(par(i)%mat)%av_l*par(i)%c +   mcm_mat(par(i)%mat)%av_q*par(i)%mass/par(i)%rho*par(i)%tracerod) / ae
   par(i)%critts =  veae/(qu + sqrt(qu**2 + par(i)%c**2))
  !
  case(1)
   !
   ! use h as length
   !
   par(i)%critts = par(i)%h/(par(i)%c + par(i)%vabs)
  !
  case(2)
   ! 
   ! use minimum interparticle distance as length
   !
   par(i)%critts = par(i)%mindist/(par(i)%c + par(i)%vabs)
  !
  end select
  !
  ! CHECK WHETHER PARTICLE i TIMESTEP IS LESS THAN CURRENT PROVISIONAL critts
  !
  if(par(i)%active) then  ! inactive particles can not control timestep
  IF (par(i)%critts .LT. mcm_critts) THEN
    mcm_critts = par(i)%critts
  ENDIF
  endif
ENDDO
!
END SUBROUTINE mcm_calccritts