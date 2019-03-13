SUBROUTINE mcm_rep_cont
!************************************************************************
!
!		Purpose: Calculate Monaghan's Repulsive Term
!
!	  Called by: Solution
!
!	     Author: Al Colebourn, James Campbell
!
!          Date: 
!
! Last Modified: 31-10-02
!
!        Errors: 
!
!         Notes: 
!
!
!************************************************************************
USE mcm_database
!
IMPLICIT NONE
!
INTEGER :: count1, count2, i, j, l, jid
!
REAL(kind=real_acc) :: ki, kij, Epsilon, nconst, fij, havg, n, hjid
real(kind=real_acc) :: inv_w_deltap, dwdx(mcm_ndim), dwdr, ratio, inv_ratio
real(kind=real_acc) :: xj(3),massj,rhoj,hj,kj
REAL(kind=real_acc) :: othird, pj
!
! Set default values
!
othird = 1.0_d/3.0_d
!
! Calculate acceleration due to repulsive force for all particles
!
do i = mcm_svp, mcm_evp
   !
   par(i)%repulsion(1:mcm_ndim) = 0.0
   !
   if(par(i)%ncont.gt.0) then
      !
      Ki = mcm_k_cont(par(i)%mat) 
      n  = mcm_n_cont(par(i)%mat)
      !
      do j =1,par(i)%ncont + par(i)%g_ncont
	     !
		 call mcm_get_j_cont_info(i,j,xj,massj,rhoj,hj,kj)
		 !pj = -othird * (sigmaj(1,1)+sigmaj(2,2)+sigmaj(3,3))
		 !if (pj.lt.0.0) pj = 0.0
	     !
	     havg = 0.5*(par(i)%h + hj) !*0.998/1.0 !hardwired for martin sauer test - tom !/2.0
		 ! uses a slightly lower h avg to ensure no initial penetration
	     !
         ! Calculate inverse of w_delta_p
         !  WARNING: The expression here assumes:
         !             1) The B-Spline kernel is being used
         !
         ratio = hj / havg !/2.0
	     IF (ratio.LT.1) THEN
	        !
	        inv_W_deltap = 1.0/ ( (0.75*ratio - 1.5)*ratio*ratio + 1.0 )
		    !
         ELSE
	        !
	        inv_W_deltap = 1.0/ (2.0-ratio)**3
		    !
	     ENDIF
	     !
	     call mcm_gradw(dwdx, dwdr, par(i)%x, xj, havg) ! only real particles need to be checked? - tom - 01-11-02
	     !
	     ! Calculate f_ij term
	     !
	!	 IF (pj.GT.1e-4) THEN
	        call mcm_calc_fij(fij,par(i)%x,xj,havg,inv_w_deltap)
	!	 ELSE
	!	    fij = 0.0
	!	 ENDIF
	     !
	     Kij = Ki + Kj
	!	 fij = 1.0
	!	 kij = 1.0/64.0
	     !
	     do l = 1,mcm_ndim
	        !
		    par(i)%repulsion(l) = par(i)%repulsion(l) - (massj*(Kij * fij**n)*dwdx(l))
	        !
	     enddo
	     !
      enddo
      !
   endif
   !
enddo
!
RETURN
!
end subroutine mcm_rep_cont
!
!
!
subroutine mcm_calc_fij(fij,xi,xj,havg,inv_w_deltap)
!************************************************************************
!
!    Purpose: Calculate f_ij term for Monaghan's repulsive force.
!
!  Called by: repulse
!
!       Date: 31-1-02
!
!     Errors: 
!
!      Notes: Uses B-Spline kernel function.
!
!   Modified: Changed to make the force dependent on the relative 
!             velocity difference.  Tom - 10/06/02
!
!************************************************************************
USE mcm_database
!
IMPLICIT NONE
!
integer :: i, j, n
real(kind=real_acc) :: fij, inv_w_deltap, havg, xij(mcm_ndim), rad, z
REAL(kind=real_acc) :: xijdotvij, vij(mcm_ndim), xj(3), xi(3),vi(3),vj(3)
!
!havg      = 0.5*(par(i)%h+par(j)%h)
!
rad       = 0.0
xijdotvij = 0.0
do n = 1,mcm_ndim
 xij(n) = xi(n) - xj(n)
 vij(n) = vi(n) - vj(n)                ! velocity difference
 xijdotvij = xijdotvij + xij(n) * vij(n) ! calculate vij component in 
                                         !the rij direction (dot product)
 rad = rad + xij(n)**2
enddo
rad = sqrt(rad)
!
z = rad/havg
! xijdotvij = xijdotvij / rad              ! normalise rij to direction cosines
!
if(z.lt.1.0) then
 !
 fij = inv_w_deltap * ( ( (0.75 * z - 1.5) * z) * z + 1.0)
 !
else if(z.lt.2.0) then
 !
 fij = 0.25 * inv_w_deltap * (2.0 - z)**3
 ! fij = 0.0 ! hardwired - tom - martin sauer test problem
 !
else
 fij = 0.0
endif
!
!if (xijdotvij.GT.0.0) fij = 0.0
!fij = fij * (1 + 0.000*xijdotvij) ! makes force dependent on the magnitude of the 
                        ! velocity difference along the axis between i 
					    ! and j particle.
!
end subroutine mcm_calc_fij
