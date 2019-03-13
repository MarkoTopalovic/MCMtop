subroutine mcm_constitutive 
!************************************************************************
!
!    Purpose: Constitutive and eos evaluation, bulk vis.
!
!  Called by: MAIN
!
!       Date: 09-08-2002
!
! Last Modified: 03-03-2005 by J. Campbell
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
!
integer:: i,nn,nn1,id
!
! Update density
!
  call mcm_rhoupd
! loop over all particles where the stress is known
!
do nn=mcm_ssp,mcm_esp
 !
 id = par(nn)%mat
 !
 select case (mcm_mat(id)%model)
  !
  case(1)
   ! Elastic
   call mcm_stress_rot1(nn)
   call mcm_f3dm1(mcm_mat(id)%strinput,nn)
   call mcm_aviscosity(nn)
   call mcm_lieupd(nn)
   !
  case(3)
   ! Elastic-plastic work hardening
   call mcm_stress_rot1(nn)
   call mcm_f3dm3(mcm_mat(id)%strinput,nn)
   call mcm_aviscosity(nn)
   call mcm_lieupd(nn)
   !
   case(4)
   ! Generalized cap model for granular materials
   call mcm_stress_rot1(nn)
   call mcm_f3dm25(mcm_mat(id)%strinput,nn)
   call mcm_aviscosity(nn)
   call mcm_lieupd(nn)
   !
  case(9)
   ! Fluid Hydrodynamic
   call mcm_stress_rot1(nn)
   call mcm_f3dm9(mcm_mat(id)%strinput,nn)
   call mcm_sueos(nn)
   call mcm_aviscosity(nn)
   call mcm_hieupd(nn)
   call mcm_eqos(nn)
   call mcm_stress_up(mcm_mat(id)%strinput,nn)
   !
  case(10)
   ! Elastic-plastic hydrodynamic
   call mcm_stress_rot1(nn)
   call mcm_f3dm10(mcm_mat(id)%strinput,nn)
   call mcm_sueos(nn)
   call mcm_aviscosity(nn)
   call mcm_hieupd(nn)
   call mcm_eqos(nn)
   call mcm_stress_up(mcm_mat(id)%strinput,nn)
   !
   !
  case default
   nn1 = mcm_mat(id)%model
   write( *,1001)nn1,nn
   write(13,1001)nn1,nn
   call mcm_shutdown(2)
 end select
enddo
!
! CALCULATE CRITICAL TIMESTEP
!
CALL mcm_calccritts
!
!
return
!
1001  format(//,'material type: ',i5,5x,//,&
      '*** illegal material  for particle ***',i8,&
             /5x,'     execution aborted ')
end subroutine mcm_constitutive
