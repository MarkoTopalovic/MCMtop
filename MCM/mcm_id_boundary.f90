subroutine mcm_id_boundary
!************************************************************************
!
!    Purpose: Identify boundary particles and calculate a normal vector
!
!  Called by: initial
!
!       Date: 05-08-2002
!
!     Errors: 
!
!      Notes: The equations implemented are those used by Libersky and Randles
!             in their sphpc code. These equations are a modification of those
!             given in R&L 96. In their code L&R normalise the gradient of unity
!             using their normal B matrix. I have not yet been able to justify
!             algebraically the use of the B matrix, but from examining the 
!             method in an Excel sheet, use of the B matrix does increase the 
!             gradient at the boundary.
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
! For now identify boundary velocity points for spvp approach in input file
!
select case (mcm_disctype)
  !
 case default
  !
  ! loop over all potential boundary particles (velocity points)
  !
  do i=mcm_svp,mcm_evp
   par(i)%boundary = 0
   ! get gradient of unity
   !call mcm_grad_unity(i,grad,grad_len)
   ! check magnitude of gradient to determine if the particle is on the boundary
   if(grad_len.gt.0.1/par(i)%h) then
    !boundary particle
    par(i)%boundary = 1
    do k=1,mcm_ndim
     par(i)%bndnorm(k) = grad(k)/grad_len
    enddo
   endif
  !
  enddo
  !
end select
!
end subroutine mcm_id_boundary