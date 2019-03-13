subroutine mcm_update_ghost_stress
!************************************************************************
!
!		Purpose:Calculate ghost particle stress and artificial viscosity tensors
!
!	  Called by: momentum
!
!	     Author: James Campbell
!
!          Date: 6-10-2005
!
! Last Modified: 6-10-2005
!
!        Errors: 
!
!         Notes: THis routine has to be called after the stress update
!
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,j,k,id
real(kind=real_acc), dimension(3,3) :: factor
!
factor = -1.0_d
do k=1,3
 factor(k,k) = 1.0_d
enddo
!
do i=1,mcm_ngp
 id = gpar(i)%par
 do k=1,3
  do j=1,3
   gpar(i)%sigma(j,k) = factor(j,k)*par(id)%sigma(j,k)
   gpar(i)%q(j,k)     = factor(j,k)*par(id)%q(j,k)
  enddo
 enddo
enddo
!
end subroutine mcm_update_ghost_stress
