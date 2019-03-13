subroutine mcm_update_velocity
!************************************************************************
!
!    Purpose: update particle velocity
!
!  Called by: solution
!
!       Date: 06-07-05
!
!     Errors: 
!
!      Notes: This was originally part of the move subroutine.
!             Modified to include dynamic relaxation
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,j,jid,k
real(kind=real_acc) :: dtn,W,Volj,epsilon
!
real(kind=real_acc) :: massj, rhoj, hj, holdj
real(kind=real_acc), dimension(3) :: xj,vj
!
! calculate delta t at time n
!
dtn = 0.5_d*(mcm_dt+mcm_dtold)
!
! Update velocities
!
if(mcm_drelax) then
 ! dynamic relaxation analysis
 do i=mcm_svp,mcm_evp
  ! absolute velocity, required for time step calculation options 1,2
  par(i)%vabs = 0.0_d
  !
  do j=1,mcm_ndim
   par(i)%v(j) = mcm_drelax_scale*par(i)%v(j) + par(i)%a(j)*dtn
   par(i)%vabs = par(i)%vabs + par(i)%v(j)**2
  enddo
  par(i)%vabs = sqrt(par(i)%vabs)
 enddo
else
 do i=mcm_svp,mcm_evp
   ! absolute velocity, required for time step calculation options 1,2
   par(i)%vabs = 0.0_d
   !
   do j=1,mcm_ndim
    par(i)%v(j) = par(i)%v(j) + par(i)%a(j)*dtn
    par(i)%vabs = par(i)%vabs + par(i)%v(j)**2
   enddo
   par(i)%vabs = sqrt(par(i)%vabs)
  enddo
endif
!
! XSPH option. Correct velocity.
!  Note that as currently programmed, v and x are held at different times.
!
if(mcm_veloc_opt.eq.1) then
 !
 epsilon = 0.100_d
 do i=1,mcm_np
  par(i)%smooth_v = 0.0_d
  !
  do j=1,par(i)%nnbr + par(i)%g_nnbr
   !
   call mcm_get_j_strain_info(i,j,xj,vj,massj,rhoj,hj,holdj)
   !
   Volj=massj/rhoj
   call mcm_kernel(w,par(i)%x,xj,hj)
   !
   Volj = massj / (0.5_d*(rhoj+par(i)%rho))
   do k=1,mcm_ndim
    par(i)%smooth_v(k) = par(i)%smooth_v(k) + Volj*(vj(k)-par(i)%v(k))*W
   enddo
  enddo
  !
 enddo
 ! update v
 do i=1,mcm_np
  do j=1,mcm_ndim
   par(i)%v(j) = par(i)%v(j) + epsilon*par(i)%smooth_v(j)
  enddo
 enddo
endif
!
! End XSPH option
!
end subroutine mcm_update_velocity