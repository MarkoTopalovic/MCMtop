subroutine mcm_get_j_strain_info(i,k,x,v,mass,rho,h,hold)
!************************************************************************
!
!    Purpose: return information for strain rate calculation on requested neighbour particle
!
!  Called by: interpolation routines
!
!       Date: 08-08-2002
!
!     Errors: 
!
!      Notes: This routine is intended to act as a uniform interface to 
!             the boundary plane routines, by removing the need to know if
!             the neighbour is a real or ghost particle in the main routines
!
!             It is assumed that 1 <= k <= par(i)%nnbr + par(i)%g_nnbr
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,k
real(kind=real_acc) :: mass, rho, h, hold
real(kind=real_acc), dimension(3) :: x,v
!
integer :: j,l,s1,s2,m,n
!
if(k.le.par(i)%nnbr) then
 ! The neighbour particle is a real particle, copy its data
 j = mcm_nbrlist(k,i)
 !
 x = par(j)%x
 v = par(j)%v
 mass = par(j)%mass
 h = par(j)%h
 hold = par(j)%hold
 rho = par(j)%rho
 !
else
 !
 ! The neighbour is a ghost particle, copy its data
 !
 l  = k - par(i)%nnbr
 j = mcm_g_nbrlist(l,i)
 !
 x = gpar(j)%x
 v = gpar(j)%v
 mass = gpar(j)%mass
 h = gpar(j)%h
 hold = gpar(j)%hold
 rho = gpar(j)%rho
 !
endif
!
end subroutine mcm_get_j_strain_info
