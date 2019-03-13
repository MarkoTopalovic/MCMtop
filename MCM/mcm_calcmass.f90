subroutine mcm_calcmass
!************************************************************************
!
!    Purpose: call correct subroutine for calculation of initial mass
!
!  Called by: Initial
!
!       Date: 01-08-2002
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
select case (mcm_massopt)
 !
 case(0)
  ! mass defined as total material mass.
  call mcm_calcmzero
  call mcm_printmass
  !
 case(1)
  ! mass of each particle given in input file
  ! just print total material masses
  call mcm_printmass
  !
end select
!
end subroutine mcm_calcmass
!
!
subroutine mcm_printmass
!************************************************************************
!
!    Purpose: Print material and total masses to log file
!
!  Called by: calcmass
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
!
integer :: i,j
integer :: count
real(kind=real_acc) :: mass
!
write(13,1000)
write(13,1100)
!
do i=1,mcm_nummat
 ! loop over particles and count mass of material i
 count = 0
 mass = 0.0_d
 !
 do j=mcm_svp,mcm_evp
  if(par(j)%mat.eq.i) then
   count = count + 1
   mass = mass + par(j)%mass
  endif
 enddo
 !
 ! Print results to log file
 !
 write(13,1200) i,mass,count
 !
enddo
!
return
!
1000 format(//5x,'PARTICLE MASS CALCULATION')
1010 format(/5x,'The list above is for the velocity particles.',&
            /5x,'The list below is for the stress particles.') 
1100 format(/5x,'Material     Total Mass  No. particles')
1200 format(5x,i8,3x,e12.5,7x,i8)
!
end subroutine mcm_printmass
!
!
SUBROUTINE mcm_calcmzero
!************************************************************************
!
!		Purpose: calculate particle masses by dividing up a supplied material
!                mass equally
!
!	  Called by: calcmass
!
!	     Author: Tom De Vuyst
!
!          Date: 01-08-2002
!
! Last Modified: 01-08-2002
!
!        Errors: 
!
!         Notes: nppermat : number of particles per material number
!                this routine is mass initialisation option 0
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i, nppermat(mcm_nummat)
real(kind=real_acc) :: r, real_mass, axi_c
!
! INITIALISE nppermat
!
DO i = 1,mcm_nummat
   nppermat(i) = 0
ENDDO
!
! CALCULATE NUMBER OF PARTICLES PER MATERIAL NUMBER
! for velocity points, where mass must be known
!
DO i = mcm_svp,mcm_evp
   nppermat(par(i)%mat) = nppermat(par(i)%mat) + 1
ENDDO
!
! CALCULATE MASS OF EVERY PARTICLE
!
DO i = mcm_svp,mcm_evp
   par(i)%mass = mcm_mat(par(i)%mat)%mass / nppermat(par(i)%mat)
   !
   ! MULTIPLY BY 2*pi*r FOR AXISYMMETRY
   !
   IF (mcm_axopt.EQ.4) THEN
      if(par(i)%x(1).eq.0.0_d) then
       par(i)%rho = mcm_mat(par(i)%mat)%rho
	   ! calculate the real mass of an axis particle 
       r = sqrt(par(i)%mass/(pi*par(i)%rho))
	   real_mass = (4.0_d/3.0_d)*pi*r*r*r*par(i)%rho
	   axi_c = real_mass/par(i)%mass
	  else
       axi_c = 2*pi*par(i)%x(1)
      endif
	  par(i)%mass = par(i)%mass*axi_c
   ENDIF
ENDDO
!
END SUBROUTINE mcm_calcmzero