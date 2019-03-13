SUBROUTINE mcm_eos41p(i)
!************************************************************************
!
!    Purpose: Calculate pressure  
!            
!             Uses Mie - Gruneisen Eos
!
!  Called by: epupdate
!
!       Date: 02-08-2002
!
!     Errors: 
!
!      Notes: MAGI - subroutine eosmg
!			  
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i
REAL(kind=real_acc) :: mu, sl, cv, t0, cc, gamma, v0
REAL(kind=real_acc) :: f, ao, bo, co, a1, b1, c1, ph, dpdrho, didrho, c2, dpdv
REAL(kind=real_acc) :: eosnum, eosden, pn, dpdi, cs, temp
REAL(kind=real_acc) :: dpde, dvol, vol, e1
!
! ASSIGN EOS PARAMETERS PER MATERIAL
!
cc    = mcm_mat(par(i)%mat)%eosinput(1)	!velocity curve intercept
sl    = mcm_mat(par(i)%mat)%eosinput(2)	!slope coefficient
cv    = mcm_mat(par(i)%mat)%eosinput(3)	!specific heat at constant volume
t0    = mcm_mat(par(i)%mat)%eosinput(4)	!initial temperature
gamma = mcm_mat(par(i)%mat)%eosinput(5)	!initial Grunesian gamma
v0    = mcm_mat(par(i)%mat)%eosinput(6)	!initial relative volume
!
! specific int. energy = e1
!
! This EoS uses internal energy per unit mass
e1=par(i)%e/par(i)%mass
vol  = 1.0_d / par(i)%rho
ao   = par(i)%rho0 * cc ** 2
mu   = par(i)%rho / par(i)%rho0 - 1.0_d
!     
bo   = ao * (1.0_d + 2.0_d * (sl - 1.0_d))
co   = ao * (2.0_d * (sl - 1.0_d) + 3.0_d * (sl - 1.0_d) ** 2)
gamma   = gamma * par(i)%rho0 / par(i)%rho
!
f    = gamma * mu
a1   = ao * (1.0_d -         f)
b1   = bo * (2.0_d - 1.5_d * f)
c1   = co * (3.0_d - 2.0_d * f)
!
if (mu .gt. 0.0_d) then
   ph     = ao * mu + bo * mu ** 2 + co * mu ** 3
   dpdrho = (a1 + b1 * mu + c1 * mu ** 2) / par(i)%rho0 + gamma * e1
else
   ph     = ao * mu
   dpdrho = ao * (1.0_d - f) / par(i)%rho0 + gamma * e1
endif
!     
pn = ph * (1.0_d - 0.5_d * f) + par(i)%rho * gamma * e1
!     
if (pn .lt. par(i)%pcut) then
    pn  = par(i)%pcut
	par(i)%rho=par(i)%rho/(1.0_d-par(i)%tracerod*mcm_dt)
else
    dpdi   = gamma * par(i)%rho
    didrho = pn / par(i)%rho ** 2
    dpdrho = dpdrho + dpdi * didrho
endif
!     
dpdv   =  - dpdrho * par(i)%rho ** 2
!
par(i)%p =  pn						!pressure
!
end subroutine mcm_eos41p
