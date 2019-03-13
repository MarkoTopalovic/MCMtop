subroutine mcm_init_neighb
!************************************************************************
!
!    Purpose: neighbour searching.
!
!  Called by: Initialise
!
!       Date: 05-08-2002
!
!     Errors: 
!
!      Notes: call subroutines for user defined problem neighbour 
!             search method.
!         *** This routine should be identical to neighbours except for
!             first deallocation of memory and variable initialisation.
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i, error
!
! Allocate memory for neighbour list arrays as ll_neighbours assumes that the arrays are allocated
allocate (mcm_nbrlist(mcm_maxnbr,mcm_np),mcm_contlist(mcm_maxcont,mcm_np),STAT=error)  
if(error.ne.0) then
   write(*,1000)
   write(13,1000)
 write(*,1100)
 write(13,1100)
   call mcm_shutdown(2)
endif
!
! Ghost particles
!
if(mcm_boundary) then
 allocate (mcm_g_nbrlist(mcm_g_maxnbr,mcm_np),mcm_g_contlist(mcm_g_maxcont,mcm_np),STAT=error)  
 if(error.ne.0) then
  write(*,1000)
  write(13,1000)
  write(*,1200)
  write(13,1200)
  call mcm_shutdown(2)
 endif
endif
!
! Initialise variables
!
do i=1,mcm_max_np
 par(i)%nnbr = 0
 par(i)%g_nnbr = 0
enddo
!
!
select case (mcm_disctype)
 case(0)
  !
  ! Conventional SPH
  !   0 = Eulerian
  ! basic linked list search
  !
  call mcm_setuplist
  if(mcm_boundary) then
   call mcm_init_boundary ! set up boundary planes, requires linked list to be active
   call mcm_ghost_setup   ! calclate coordinates and values of ghost particles
  endif
  call mcm_ll_neighbours(0)
    !
 case(1)
  !
  ! Corrected SPH
  !   1 = Eulerian
  !   3 = Total Lagrangian
  ! basic linked list search, with i particle included
  !
  call mcm_setuplist
  if(mcm_boundary) then
   call mcm_init_boundary ! set up boundary planes, requires linked list to be active
   call mcm_ghost_setup   ! calclate coordinates and values of ghost particles
  endif
  call mcm_ll_neighbours(1)
    !
end select
!
return
!
1000 format(//5x,'Error in subroutine init_neighb')
1100 format(/5x,'Error allocating neighbour list array memory')
1200 format(/5x,'Error allocating ghost neighbour list array memory')
!
end subroutine mcm_init_neighb