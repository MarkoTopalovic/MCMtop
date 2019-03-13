SUBROUTINE mcm_hieupd(i)
!************************************************************************
!
!    Purpose: Calculate trial value of e
!
!  Called by: epupdate(i)
!
!       Date: 14-1-99
!
!     Errors: 
!
!      Notes: WARNING: Assumes that the viscosity tensor is symmetric
!                      According to Benson (1992) this may not always be the case
!
!************************************************************************
!
USE mcm_database
!
IMPLICIT NONE
!
REAL(kind=real_acc) :: de(3,3)
REAL(kind=real_acc) :: pold, volold, volnew, vavg, dvol, eincr
REAL(kind=real_acc) :: trace_q, mean_incr,temp
INTEGER				:: i, j, l, k
REAL(kind=real_acc) :: othird
!
othird = 1.0_d / 3.0_d
pold = -othird * (par(i)%sigma(1,1) + par(i)%sigma(2,2) + par(i)%sigma(3,3))
trace_q = othird*(par(i)%q(1,1) + par(i)%q(2,2) + par(i)%q(3,3))
!
volold = par(i)%mass/par(i)%rhoold
volnew = par(i)%mass/par(i)%rho
!
! calculate energy increment due to deviatoric terms
! 
do l=1,3
 do k=1,3
  if(k.eq.l) then
   de(k,l) = par(i)%rod(k,l)*(0.5*(par(i)%sigma(k,l)+pold  + par(i)%s(k,l))  - par(i)%q(k,l)+trace_q)
  else
   de(k,l) = par(i)%rod(k,l)*(0.5*(par(i)%sigma(k,l) + par(i)%s(k,l)) - par(i)%q(k,l))
  endif
 enddo
enddo
!
vavg = 0.5*(volnew + volold)
eincr = 0.0
do l=1,3
 do k=1,3
  eincr = eincr + de(k,l)
 enddo
enddo
eincr = vavg*eincr
!
! calculate energy increment due to hydrostatic (mean) terms
!
dvol = volnew - volold
mean_incr = dvol * (0.5*pold + trace_q)
!
! calculate trial value of internal energy
!
par(i)%etry = par(i)%e + mcm_dt*eincr - mean_incr
!
par(i)%qold(1:3,1:3) = par(i)%q(1:3,1:3)
!
RETURN

END SUBROUTINE mcm_hieupd