SUBROUTINE mcm_eqos41(i)
!************************************************************************
!
!    Purpose: Calculate pressure and internal energy based on trial value 
!             of internal energy 
!            
!             Uses Mie - Gruneisen Eos
!
!  Called by: epupdate
!
!       Date: 6-1-99
!
!     Errors: 
!
!      Notes: MAGI - subroutine eosmg
!			  
!
!************************************************************************
!
USE mcm_database
!
IMPLICIT NONE
!
INTEGER :: i
REAL(kind=real_acc) :: mu, sl, cv, t0, cc, gama,gama0,beta,v0
REAL(kind=real_acc) :: f, ao, bo, co, a1, b1, c1, ph, dpdrho, didrho, c2, dpdv
REAL(kind=real_acc) :: eosnum, eosden, pn, dpdi, cs, temp
REAL(kind=real_acc) :: dpde, dvol, vol,vol0, etry1
!
! ASSIGN EOS PARAMETERS PER MATERIAL
!
cc    = mcm_mat(par(i)%mat)%eosinput(1)	!velocity curve intercept
sl    = mcm_mat(par(i)%mat)%eosinput(2)	!slope coefficient
cv    = mcm_mat(par(i)%mat)%eosinput(3)	!specific heat at constant volume
t0    = mcm_mat(par(i)%mat)%eosinput(4)	!initial temperature
gama0 = mcm_mat(par(i)%mat)%eosinput(5)	!initial Grunesian gama
v0    = mcm_mat(par(i)%mat)%eosinput(6)	!initial relative volume
beta  = mcm_mat(par(i)%mat)%eosinput(7)	!beta= volume correction factor for gama0
!
! specific 'trial' value of int. energy = etry
!
etry1 = par(i)%etry / par(i)%mass
vol   = 1.0_d / par(i)%rho
ao    = par(i)%rho0 * cc**2
mu    = par(i)%rho/par(i)%rho0 - 1.0
dvol  = par(i)%mass/par(i)%rho - par(i)%mass/par(i)%rhoold
vol0  = par(i)%mass/par(i)%rho0
!
!calculate beta that makes gama thermodynamically stable
!	0.5*(gama0-(2.0/mu)<gama<(1.0+gama0)
!
!	Correction for gama according to Segletes
!gama=gama0/(1.0+beta*mu)
!	Larry's correction for gama  
gama = gama0 * par(i)%rho0 / par(i)%rho
!     
bo   = ao * (1. + 2. * (sl - 1.))
co   = ao * (2. * (sl - 1.) + 3. * (sl - 1.) ** 2)
!
f    = gama * mu
a1   = ao * (1. -       f)
b1   = bo * (2. - 1.5 * f)
c1   = co * (3. - 2.0 * f)
!
if (mu .gt. 0.0) then
   !
   ph     = ao * mu + bo * mu ** 2 + co * mu ** 3
   dpdrho = (a1 + b1 * mu + c1 * mu ** 2) / par(i)%rho0 + gama * etry1
   !
else
   !
   ph     = ao * mu
   dpdrho = ao * (1. - f) / par(i)%rho0 + gama * etry1
   !
endif
!
!Larry's expression for pressure     
!pn = ph * (1. - 0.5 * f) + rho(i) * gama * etry1
!
!Corrected value for pressure used in DYNA
pn = (ph * (1. - 0.5 * f) + par(i)%rho * gama * etry1)/(1.0+0.5*par(i)%rho*gama*dvol/vol0)
! 
! Pressure cutoff
!
!pn = max(pn,pcut(i))    
if (pn .lt. par(i)%pcut) then
   !
   pn         = par(i)%pcut
   par(i)%rho = par(i)%rho/(1.0-par(i)%tracerod*mcm_dt)
   !	h(i)=hold(i)
else
   !
   dpdi     = gama * par(i)%rho
   didrho   = pn / par(i)%rho ** 2
   dpdrho   = dpdrho + dpdi * didrho
   par(i)%e = etry1*par(i)%mass-0.5*dvol*pn
   !
endif
!     
dpdv      =  - dpdrho * par(i)%rho ** 2
!
par(i)%p  =  pn						!pressure
!
!
! *** solid temperature along shock hugoniot
!
!        if (vol .lt. v0) then
!           gama = = eosinput(5,pmat(i)) 
!           call htemp(rho(i), v0, cc, sl, gama, cv, temp, t0)
!        else
!           temp = t0
!        endif
!
! *** solid temperature assuming constant specific heat
!
!        temp = t0 + etry1 / cv
!
! *** solid temperature (Larry's routine)
!
!        if (vol .lt. v0) then
!           call xtemp (vol, etry1, v0, cc, sl, gama, cv, temp)
!        else
!           temp = t0
!        endif
!
      return
      end
! =====================================================================
!
subroutine mcm_htemp (rho, v0, c, sl, go, cv, t, t0)
!
!     This routine calculates the temperature increase along
!     the principal hugoniot. It assumes a constant specific
!     heat and a gruneisen parameter proportaional to density.
!     This routine will not work for densities greater than 
!     the limit of the hugoniot, vol/v0 < 1 - 1/sl
!
implicit none
!
real ei
real cv, c2, v0, c, sl, go, s2, s3
real t, t0, rho, vr
real one, two, three
!
vr = 1.0/(v0*rho)
c2 = c **2
s2 = sl **2
s3 = sl **3
!
!     test hugoniot limit
!
if (1. - sl*(1.-vr) .gt. 0) then
   !
   !     calculate temperature of solid along hugoniot
   !
   one = (3*sl-go)/sl*exp(go*(1-vr))
   two = ((3*sl-go)-sl*(4*sl-go)*(1-vr))/(sl*(1-sl*(1-vr))**2)
   three = (go**2-2*(go-sl)**2)/s2*exp(-go/sl*(1-sl*(1-vr)))   &
           *(ei(go/sl,1.-sl*(1.-vr))-ei(go/sl, 1.0))
   t = t0*exp(go*(1-vr)) + c2/(2*s2*cv)*(one - two + three)
   !
else
   !
   !     outside hugoniot limit, give up
   !
   t = t0
   !
end if
!
return
end
!
! =====================================================================
!
subroutine mcm_xtemp (v, e, v0, c, sl, go, cv, t)
!
implicit none
!
real ei
real t0, cv, a, c2, v0, c, sl, go, s2, s3, s4, a2
real fa, ea, ta, t, v, e
!
t0 = 300.
!
a  = go / v0
c2 = c**2
s2 = sl**2
s3 = sl**3
s4 = sl**4
a2 = a**2
!
fa = c2 * exp(a*v0) * (a*v0-3.*sl) / (2.*s3) +						&
     exp(a*v)*														&
     (c2/(2.*s2) + c2*v0*(2.*sl-a*v0)/(2.*s3*(sl*v+v0-sl*v0)))   -	&
     (c2*exp(a*(sl-1.)*v0/sl) * (2.*s2 - 4.*a*sl*v0 + a2*v0**2)) *	&
     (ei(a,v0/sl) - ei(a,v-v0+v0/sl)) / (2.*s4)
!
ea = fa * exp(-a * v)
ta = t0 * exp(go) / exp(a * v)
t  = ta + (e - ea) / cv
!
return
end
!
!==============================================================================
!
real function ei(a, y)
!
!     exponential integral function
!     Ei(x) = S_{(-infinity,x]} exp(t)/t dt (x<0)
!
implicit none
!
integer maxit
real    a, y, x, eps, euler, fpmin
parameter (eps = 6.e-8, euler = .57721566)
parameter (maxit = 100, fpmin = 1.e-30)
integer k
real fact, prev, sum, term
!
x = a * y
!
if (y .lt. 0.) then
   write (6,'(''bad arguement in ei'')')
   stop
endif
!
if (x .lt. fpmin) then
   ei = log(y)
elseif (x .le. -log(eps)) then
   sum  = 0.
   fact = 1.
   do k = 1, maxit
      fact = fact * x / k
      term = fact / k
      sum  = sum + term
      if (term .lt. eps * sum) goto 1
   enddo
   write (6,'(''series failed in ei'')')
   stop
1  ei = sum + log(y)
else
   sum  = 0.
   term = 1.
   do k = 1, maxit
      prev = term
      term = term * k / x
      if (term .lt. eps) goto 2
      if (term .lt. prev) then
         sum = sum + term
      else
         sum = sum - prev
         goto 2
      endif
   enddo
2  ei = exp(x) * (1. + sum) / x
endif
!
return
end    
