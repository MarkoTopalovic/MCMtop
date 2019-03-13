subroutine mcm_ghost_ll_neighbours(enough_mem)
!************************************************************************
!
!    Purpose: find all ghost particle neighbours from linked list
!
!  Called by: neighbours
!
!       Date: 7-10-2005
!
!     Errors: 
!
!      Notes: uses linked list array set up in subroutine setuplist
!             routine is based on ll_neighbours routine
!             
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,j,k,l,m,n
integer :: id, boxcrd(3),looplim(2,3),error
real(kind=real_acc) :: dist_sq,h_sq,mindist_cp,h_avg,factor
logical :: enough_mem
integer, dimension(2,3) :: period_box
!
factor = 2.0_d
!
enough_mem = .true.
!
do i=1,3
 boxcrd(i) = 1
 looplim(1,i) = 1
 looplim(2,i) = 1
enddo
!
! loop over all particles
!
do i=1,mcm_np
 par(i)%g_nnbr = 0
 par(i)%g_ncont = 0
 mindist_cp = 1.0e20_d
 !
 ! Only do search is particle is inside the sort domain,
 !  otherwise leave nnbr = 0
 !
 if(par(i)%llpointer.gt.-1) then
  !
  ! identify which cell the particle lies in, and prevent loop over grid cells
  ! that do not exist
  !
  do j=1,mcm_ndim
   boxcrd(j)=int(mcm_gridsize(j)*(par(i)%x(j)-mcm_gridmin(j))) + 1
   looplim(1,j)=max(boxcrd(j)-1,1)
   looplim(2,j)=min(boxcrd(j)+1,mcm_gridlim(j))
  enddo
  !
  ! loop over neighbour grid cells
  !
  do n=looplim(1,3),looplim(2,3)
   do m=looplim(1,2),looplim(2,2)
    do l=looplim(1,1),looplim(2,1)
     ! get id of first particle in cell
     id=mcm_llgrid(l,m,n) 
 	 ! loop over all particles in the cell
	 do
	  ! exit if end of list has been reached, placed here so that is the 
	  ! cell is empty, the do loop is immediately exited
	  if(id.gt.-1) exit
	  id=abs(id)
	  dist_sq = 0.0_d
	  do j=1,mcm_ndim
	   dist_sq = dist_sq + (gpar(id)%x(j)-par(i)%x(j))**2
      enddo
	  h_avg = factor * 0.5_d * (gpar(id)%h+par(i)%h)
      h_sq = h_avg*h_avg
      ! check if particle is within 2h of i particle
	  if(dist_sq.lt.h_sq) then
	   ! if they are of the same material then add to the neighbour list
	   if(gpar(id)%mat.eq.par(i)%mat) then
        ! check for minimum interparticle distance, used for timestep calculation option 2
		par(i)%mindist = min(par(i)%mindist,dist_sq)
		par(i)%g_nnbr = par(i)%g_nnbr + 1
		! check if exceeded max neighbour limit
		if(par(i)%g_nnbr.le.mcm_g_maxnbr) then
		 mcm_g_nbrlist(par(i)%g_nnbr,i) = id
		else
	     enough_mem = .false.
		endif
	   else
	    !otherwise add it to the list of potential contact particles
	    ! not worrying about whether the particles are boundary particles
        ! as the contact list will be used to augment the neighbour list if
   	    ! kernel contact is being used
        !
	    ! check for minimum interparticle distance, used for timestep calculation option 2
	    mindist_cp = min(mindist_cp,dist_sq)
		par(i)%g_ncont = par(i)%g_ncont +1
		if(par(i)%g_ncont.le.mcm_g_maxcont) then
		 mcm_g_contlist(par(i)%g_ncont,i) = id
		else
	     enough_mem = .false.
		endif
	   endif
	  endif
	  ! get new id
	  id = gpar(id)%llpointer
	 enddo
    enddo
   enddo
  enddo
  !
  ! if using kernel contact then augment neighbour list with contact list
  !
  if(mcm_contacttype.eq.0) then
   !
   if(enough_mem) then
    if((par(i)%g_nnbr + par(i)%g_ncont).le.mcm_maxnbr) then
     do j=1,par(i)%g_ncont
      if ((par(i)%g_nnbr + par(i)%g_ncont).eq.0) exit
      mcm_g_nbrlist(par(i)%g_nnbr+j,i) = mcm_g_contlist(j,i)
     enddo
    else
     enough_mem = .false.
    endif
   endif
   par(i)%g_nnbr = par(i)%g_nnbr + par(i)%g_ncont
   par(i)%mindist = min(par(i)%mindist,mindist_cp)
  endif
  par(i)%mindist = sqrt(par(i)%mindist)
 endif  ! particle inside domain if  (par(i)%llpointer.gt.-1)
enddo
!
! find new maxnbr
mcm_maxnbr = 0
mcm_maxcont= 0
do i=1,mcm_np
 mcm_maxnbr = max(mcm_maxnbr,par(i)%nnbr)
 mcm_maxcont = max(mcm_maxcont,par(i)%ncont)
enddo
!
return
!
1000 format(//5x,'Error in subroutine ll_neighbours')
1010 format(/5x,'Unexpected error when deallocating linked list array memory')
!
end subroutine mcm_ghost_ll_neighbours
