SUBROUTINE mcm_eos41c(i)
!************************************************************************
!
!    Purpose: Calculate dp/drho, dp/dE, c FOR MIE GRUNEISEN EOS
!             Uses Mie - Gruneisen Eos
!
!  Called by: epupdate
!
!       Date: 05-08-2002
!
!     Errors: 
!
!      Notes: speed of sound c(i)
!			  MAGI - subroutine eosmg
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i
REAL(kind=real_acc) :: mu, sl, cv, t0, cc, gamma,gamma0, v0
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
gamma0= mcm_mat(par(i)%mat)%eosinput(5)	!initial Grunesian gamma
v0    = mcm_mat(par(i)%mat)%eosinput(6)	!initial relative volume
! specific 'trial' value of int. energy = etry
!
e1=par(i)%e/par(i)%mass
vol  = 1.0_d / par(i)%rho
ao   = par(i)%rho0 * cc ** 2
mu   = par(i)%rho / par(i)%rho0 - 1.0_d
dvol = par(i)%mass/par(i)%rho-par(i)%mass/par(i)%rhoold
!
!	Correction for gama according to Segletes
!gama=gama0/(1.0+beta*mu)
!	Larry's correction for gama  
gamma   = gamma0 * par(i)%rho0 / par(i)%rho
!          
bo   = ao * (1.0_d + 2.0_d * (sl - 1.0_d))
co   = ao * (2.0_d * (sl - 1.0_d) + 3.0_d * (sl - 1.0_d) ** 2)
!
f    = gamma * mu
a1   = ao * (1.0_d -         f)
b1   = bo * (2.0_d - 1.5_d * f)
c1   = co * (3.0_d - 2.0_d * f)
!
if (mu .gt. 0.0_d) then
   dpdrho = (a1 + b1 * mu + c1 * mu ** 2) / par(i)%rho0 + gamma * e1
else
   dpdrho = ao * (1.0_d - f) / par(i)%rho0 + gamma * e1
endif
!     
if (par(i)%p .eq. par(i)%pcut) then
    c2  = cc * cc
else
    dpdi   = gamma * par(i)%rho
    didrho = par(i)%p / par(i)%rho ** 2
    dpdrho = dpdrho + dpdi * didrho
    c2     = dpdrho
endif
!     
if (c2 .lt. 0.0_d) then
	write(*,1000)
	write(13,1000)
	 write(*,1100) i,c2,mcm_mat(par(i)%mat)%g,dpdrho,par(i)%p,par(i)%rho
	write(13,1100) i,c2,mcm_mat(par(i)%mat)%g,dpdrho,par(i)%p,par(i)%rho
	call mcm_shutdown(2)
endif
!   
1000 format(//5x,'Error in subroutine eos41c')
1100 format(/5x,' Particle: ',i6,/5x,'c2: 'e14.6,/5x,'g(pmat(i): ',e14.6,&
            /5x,'dpdrho: ',e14.6,/5x,'p(i): ',e14.6,/5x,'rho(i): ',e14.6)
!
par(i)%c      =  sqrt(c2)				!speed of sound
!
end subroutine mcm_eos41c
