subroutine mcm_ghost_setup
!************************************************************************
!
!		Purpose: Set up ghost particles for current timestep
!
!	  Called by: neighbours, init_neighb
!
!	     Author: James Campbell
!
!          Date: 6-10-2005
!
! Last Modified: 6-10-2005
!
!        Errors: 
!
!         Notes: 
!
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,j,k,id,m
!
logical :: enough_mem
!
mcm_ngp = 0
enough_mem = .true.
!
do m=1,2
 !
 !------------------------------------------------------------------------
 !
 !  X min boundary
 !
 select case(mcm_boundary_code(1,1))
  case(1)
   i=2
  case(2)
   i=mcm_gridlim(1) - 1
 end select
 !
 select case(mcm_boundary_code(1,1))
  case(1,2)
   !
   ! loop over all adjacent cells and create ghost particles from real particles
   do k=1,mcm_gridlim(3)
    do j=1,mcm_gridlim(2)
     id=mcm_llgrid(i,j,k)
     if(id.gt.0) then
      ! cell contains real particles
      do
       if(id.eq.0) exit 
       call mcm_create_ghost_from_real(id,1,j,k,mcm_xmin_mult, mcm_xmin_add,enough_mem)
       ! get new id
       id = par(id)%llpointer
      enddo
     else
      ! cell is empty or contains ghost particles
      do
       if(id.eq.0) exit
       call mcm_create_ghost_from_ghost(abs(id),1,j,k,mcm_xmin_mult, mcm_xmin_add,enough_mem)
       ! get new id
       id = gpar(abs(id))%llpointer
      enddo
     endif
    enddo
   enddo
 end select
 !
 !------------------------------------------------------------------------
 !
 !  X max boundary
 !
 select case(mcm_boundary_code(2,1))
  case(1)
   i=mcm_gridlim(1) - 1
  case(2)
   i=2
 end select
 !
 select case(mcm_boundary_code(2,1))
  case(1,2)
   !
   ! loop over all adjacent cells and create ghost particles from real particles
   do k=1,mcm_gridlim(3)
    do j=1,mcm_gridlim(2)
     id=mcm_llgrid(i,j,k)
     if(id.gt.0) then
      ! cell contains real particles
      do
       if(id.eq.0) exit
	   call mcm_create_ghost_from_real(id,mcm_gridlim(1),j,k, mcm_xmax_mult, mcm_xmax_add,enough_mem)
       ! get new id
       id = par(id)%llpointer
      enddo
     else
      ! cell is empty or contains ghost particles
      do
       if(id.eq.0) exit
       call mcm_create_ghost_from_ghost(abs(id),mcm_gridlim(1),j,k,mcm_xmax_mult, mcm_xmax_add,enough_mem)
       ! get new id
       id = gpar(abs(id))%llpointer
      enddo
     endif
    enddo
   enddo
 end select
 !
 !------------------------------------------------------------------------
 !
 !  Y min boundary
 !
 select case(mcm_boundary_code(1,2))
  case(1)
   j=2
  case(2)
   j=mcm_gridlim(2) - 1
 end select
 !
 select case(mcm_boundary_code(1,2))
  case(1,2)
   !
   ! loop over all adjacent cells and create ghost particles from real particles
   do k=1,mcm_gridlim(3)
    do i=1,mcm_gridlim(1)
     id=mcm_llgrid(i,j,k)
     if(id.gt.0) then
      ! cell contains real particles
      do
       if(id.eq.0) exit
	   call mcm_create_ghost_from_real(id,i,1,k, mcm_ymin_mult, mcm_ymin_add,enough_mem)
       ! get new id
       id = par(id)%llpointer
      enddo
     else
      ! cell is empty or contains ghost particles
      do
       if(id.eq.0) exit
       call mcm_create_ghost_from_ghost(abs(id),i,1,k,mcm_ymin_mult, mcm_ymin_add,enough_mem)
       ! get new id
       id = gpar(abs(id))%llpointer
      enddo
     endif
    enddo
   enddo
 end select
 !
 !------------------------------------------------------------------------
 !
 !  Y max boundary
 !
 select case(mcm_boundary_code(2,2))
  case(1)
   j=mcm_gridlim(2) - 1
  case(2)
   j=2
 end select
 !
 select case(mcm_boundary_code(2,2))
  case(1,2)
   !
   ! loop over all adjacent cells and create ghost particles from real particles
   do k=1,mcm_gridlim(3)
    do i=1,mcm_gridlim(1)
     id=mcm_llgrid(i,j,k)
     if(id.gt.0) then
      ! cell contains real particles
      do
       if(id.eq.0) exit
	   call mcm_create_ghost_from_real(id,i,mcm_gridlim(2),k, mcm_ymax_mult, mcm_ymax_add,enough_mem)
       ! get new id
       id = par(id)%llpointer
      enddo
     else
      ! cell is empty or contains ghost particles
      do
       if(id.eq.0) exit
       call mcm_create_ghost_from_ghost(abs(id),i,mcm_gridlim(2),k,mcm_ymax_mult, mcm_ymax_add,enough_mem)
       ! get new id
       id = gpar(abs(id))%llpointer
      enddo
     endif
    enddo
   enddo
 end select
 !
 !------------------------------------------------------------------------
 !
 !  Z min boundary
 !
 select case(mcm_boundary_code(1,3))
  case(1)
   k=2
  case(2)
   k=mcm_gridlim(3) - 1
 end select
 !
 select case(mcm_boundary_code(1,3))
  case(1,2)
   !
   ! loop over all adjacent cells and create ghost particles from real particles
   do j=1,mcm_gridlim(2)
    do i=1,mcm_gridlim(1)
     id=mcm_llgrid(i,j,k)
     if(id.gt.0) then
      ! cell contains real particles
      do
       if(id.eq.0) exit
	   call mcm_create_ghost_from_real(id,i,j,1, mcm_zmin_mult, mcm_zmin_add,enough_mem)
       ! get new id
       id = par(id)%llpointer
      enddo
     else
      ! cell is empty or contains ghost particles
      do
       if(id.eq.0) exit
       call mcm_create_ghost_from_ghost(abs(id),i,j,1,mcm_zmin_mult, mcm_zmin_add,enough_mem)
       ! get new id
       id = gpar(abs(id))%llpointer
      enddo
     endif
    enddo
   enddo
 end select
 !
 !------------------------------------------------------------------------
 !
 !  Z max boundary
 !
 select case(mcm_boundary_code(2,3))
  case(1)
   k=mcm_gridlim(3) - 1
  case(2)
   k=2
 end select
 !
 select case(mcm_boundary_code(2,3))
  case(1,2)
   !
   ! loop over all adjacent cells and create ghost particles from real particles
   do j=1,mcm_gridlim(2)
    do i=1,mcm_gridlim(1)
     id=mcm_llgrid(i,j,k)
     if(id.gt.0) then
      ! cell contains real particles
      do
       if(id.eq.0) exit
	   call mcm_create_ghost_from_real(id,i,j,mcm_gridlim(3), mcm_zmax_mult, mcm_zmax_add,enough_mem)
       ! get new id
       id = par(id)%llpointer
      enddo
     else
      ! cell is empty or contains ghost particles
      do
       if(id.eq.0) exit
       call mcm_create_ghost_from_ghost(abs(id),i,j,mcm_gridlim(3),mcm_zmax_mult, mcm_zmax_add,enough_mem)
       ! get new id
       id = gpar(abs(id))%llpointer
      enddo
     endif
    enddo
   enddo
 end select
 !
 !------------------------------------------------------------------------
 if(enough_mem) exit
 !
 if(m.eq.1) then
  ! first time through loop and number of ghost particles exceeds storage
  write(*,1000)
  write(13,1000)
  call mcm_init_boundary
 else
  ! getting to this point should be impossible but call shutdown incase
  ! programming error causes the code to get here. 
  write(*,2000)
  write(13,2000)
  call mcm_shutdown(2)
 endif
 !
enddo
!
1000 format(/5x,'Extra memory required for ghost particles')
2000 format(/5x,'Fatal error in routine mcm_ghost_setup, please report.')
!
end subroutine mcm_ghost_setup
!
!========================================================================
!
subroutine mcm_create_ghost_from_real(id,i,j,k,mult,add,enough_mem)
!************************************************************************
!
!		Purpose: create ghost particle from real particle
!
!	  Called by: ghost_setup
!
!	     Author: James Campbell
!
!          Date: 7-10-2005
!
! Last Modified: 7-10-2005
!
!        Errors: 
!
!         Notes: This routine creates a ghost particle based on particle id.
!                The new ghost is added to linked list cell i,j,k.
!                mult and add are used to calculate the new coordinates.
!                Stress and viscosity are not calculated for the ghost here
!                as they are only required after the stress update
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: id, i,j,k,m
real(kind=real_acc), dimension(3) :: mult, add
logical :: enough_mem
!
if(enough_mem) then
 !
 mcm_ngp = mcm_ngp + 1
 !
 ! copy particle data
 !
 gpar(mcm_ngp)%par = id
 gpar(mcm_ngp)%mat  = par(id)%mat
 gpar(mcm_ngp)%mass = par(id)%mass
 gpar(mcm_ngp)%h    = par(id)%h
 gpar(mcm_ngp)%hold = par(id)%hold
 gpar(mcm_ngp)%rho  = par(id)%rho
 !
 do m=1,3
  gpar(mcm_ngp)%x(m) = par(id)%x(m)*mult(m) + add(m)
  gpar(mcm_ngp)%xzero(m) = par(id)%xzero(m)*mult(m) + add(m)
  gpar(mcm_ngp)%v(m) = par(id)%v(m)*mult(m)
 enddo
 !
 ! Update linked list
 !
 ! move pointer to first particle in box to the particle pointer
 !
 gpar(mcm_ngp)%llpointer = mcm_llgrid(i,j,k)
 !
 ! set the pointer to the first particle in the box
 !  the minus is used to signify ghost particles
 !
 mcm_llgrid(i,j,k) = -mcm_ngp
 !
 if(mcm_ngp.eq.mcm_max_ngp) enough_mem = .false.
endif
!
end subroutine mcm_create_ghost_from_real
!
!========================================================================
!
subroutine mcm_create_ghost_from_ghost(id,i,j,k,mult,add,enough_mem)
!************************************************************************
!
!		Purpose: create ghost particle from ghost particle
!
!	  Called by: ghost_setup
!
!	     Author: James Campbell
!
!          Date: 7-10-2005
!
! Last Modified: 7-10-2005
!
!        Errors: 
!
!         Notes: This routine creates a ghost particle based on ghost particle id.
!                The new ghost is added to linked list cell i,j,k.
!                mult and add are used to calculate the new coordinates.
!                Stress and viscosity are not calculated for the ghost here
!                as they are only required after the stress update
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: id, i,j,k,m
real(kind=real_acc), dimension(3) :: mult, add
logical :: enough_mem
!
if(enough_mem) then
 !
 mcm_ngp = mcm_ngp + 1
 !
 ! copy particle data
 !
 gpar(mcm_ngp)%par  = gpar(id)%par
 gpar(mcm_ngp)%mat  = gpar(id)%mat
 gpar(mcm_ngp)%mass = gpar(id)%mass
 gpar(mcm_ngp)%h    = gpar(id)%h
 gpar(mcm_ngp)%hold = gpar(id)%hold
 gpar(mcm_ngp)%rho  = gpar(id)%rho
 !
 do m=1,3
  gpar(mcm_ngp)%x(m)     = gpar(id)%x(m)*mult(m) + add(m)
  gpar(mcm_ngp)%xzero(m) = gpar(id)%xzero(m)*mult(m) + add(m)
  gpar(mcm_ngp)%v(m)     = gpar(id)%v(m)*mult(m)
 enddo
 !
 ! Update linked list
 !
 ! move pointer to first particle in box to the particle pointer
 !
 gpar(mcm_ngp)%llpointer = mcm_llgrid(i,j,k)
 !
 ! set the pointer to the first particle in the box
 !  the minus is used to signify ghost particles
 !
 mcm_llgrid(i,j,k) = -mcm_ngp
 !
 if(mcm_ngp.eq.mcm_max_ngp) enough_mem = .false.
endif
!
end subroutine mcm_create_ghost_from_ghost
