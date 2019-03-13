subroutine mcm_get_j_cont_info(i,k,x,mass,rho,h,kj)
!************************************************************************
!
!    Purpose: return all information on requested neighbour contact particle
!
!  Called by: interpolation routines
!
!       Date: 08-08-2002
!
!     Errors: 
!
!      Notes: This routine is intended to act as a uniform interface to 
!             the boundary plane routines, by removing the need to know if
!             the neighbour is a real or mirror particle in the main routines
!
!             It is assumed that 1 <= k <= par(i)%ncont + par(i)%n_symcont
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,k
real(kind=real_acc) :: mass, rho, h, kj
real(kind=real_acc), dimension(3) :: x,v
real(kind=real_acc), dimension(3,3) :: sigma, q
!
integer :: j,l,s1,s2,m,n
!
if(k.le.par(i)%ncont) then
 ! The neighbour particle is a real particle, copy the data
 j = mcm_contlist(k,i)
 !
 !
 Kj = mcm_k_cont(par(j)%mat) ! get contact parameter
 !
 x = par(j)%x
 mass = par(j)%mass
 rho = par(j)%rho
 h = par(j)%h
 !
else
 !
 ! The neighbour is a mirror particle (k > par(i)%nnbr)
 !
 l  = k - par(i)%ncont
 j = mcm_contlist(k,i)
 !
 !
 Kj = mcm_k_cont(gpar(j)%mat) ! get contact parameter
 !
 x = gpar(j)%x
 mass = gpar(j)%mass
 rho = gpar(j)%rho
 h = gpar(j)%h
 !
endif
!
end subroutine mcm_get_j_cont_info
