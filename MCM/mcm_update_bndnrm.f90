subroutine mcm_update_bndnrm
!************************************************************************
!
!    Purpose: Identify boundary particles and calculate a normal vector
!
!  Called by: neighbour
!
!       Date: 14-08-2002
!
!     Errors: 
!
!      Notes: For boundary particles, calculate updated boundary normal
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,k,l
real(kind=real_acc) :: grad_len
real(kind=real_acc) :: grad(mcm_ndim)
!
!
!
select case (mcm_disctype)
  !
 case default
  !
  ! loop over all potential boundary particles (velocity points)
  !
  do i=mcm_svp,mcm_evp
   if(par(i)%boundary.eq.1) then
    ! get gradient of unity
    !call mcm_grad_unity(i,grad,grad_len)
    ! set normal
    do k=1,mcm_ndim
     par(i)%bndnorm(k) = grad(k)/grad_len
    enddo
   endif
   !
  enddo
  !
end select
!
end subroutine mcm_update_bndnrm