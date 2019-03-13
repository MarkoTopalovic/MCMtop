subroutine mcm_neighbours
!************************************************************************
!
!    Purpose: Main neighbour searching routine
!
!  Called by: Solution
!
!       Date: 05-08-2002
!
!     Errors: 
!
!      Notes: call subroutines for user defined problem neighbour 
!             search method.
!         *** Changes to this routine should also be made to init_neighb
!             if necessary.
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i
!
integer :: j,k
character(len=22) :: thformat
!
!
select case (mcm_disctype)
  case(0)
    !
    ! Conventional SPH
    ! basic linked list search
    !
    call mcm_setuplist
  if(mcm_boundary) call mcm_ghost_setup   ! calclate coordinates and values of ghost particles
  call mcm_ll_neighbours(0)
    !
    case(1)
    !
    ! Corrected SPH
    ! basic linked list search
    !
    call mcm_setuplist
  if(mcm_boundary) call mcm_ghost_setup   ! calclate coordinates and values of ghost particles
  call mcm_ll_neighbours(1)               ! 1 means augment neighbourhood with i particle
    !
    !
end select
!
do i=1,mcm_np
 k=par(i)%nnbr
 if (k.lt.2) then
  write(unit=thformat,fmt=500) k 
 else
  write(unit=thformat,fmt=600) k
 endif
 !
 !write(130,fmt=thformat) mcm_timestep,i,k,(mcm_nbrlist(j,i),j=1,k)
 ! topalovic ovo izgleda kao neka zastarela kontrolna stampa
enddo
500 format('(3(i4,'',''),',i1,'(i4,'',''))')
600 format('(3(i4,'',''),',i2,'(i4,'',''))')
end subroutine mcm_neighbours
