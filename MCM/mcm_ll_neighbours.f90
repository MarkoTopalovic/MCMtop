subroutine mcm_ll_neighbours(ip_flag)
!************************************************************************
!
!    Purpose: find all particle neighbours from linked list
!
!  Called by: neighbours
!
!       Date: 05-08-2002
!
!     Errors: 
!
!      Notes: uses linked list array set up in subroutine setuplist
!             
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: i,j,k,l,m,n,count, ip_flag
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
! Put in safety factor on max neighbour
!
mcm_maxnbr = mcm_maxnbr + 5 + ip_flag
mcm_maxcont = mcm_maxcont + 5
!
! start of loop to allow neighbour array memory to expand if requred
!
do count = 1,2
 !
 ! deallocate neighbour list memory
 !
 deallocate(mcm_nbrlist,mcm_contlist,STAT=error)
 if(error.ne.0) then
  write(*,1000)
  write(13,1000)
  write(*,1010)
  write(13,1010)
  call mcm_shutdown(2)
 endif
 !
 ! allocate neighbour list memory with new maximum size
 !
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
  !
  deallocate(mcm_g_nbrlist,mcm_g_contlist,STAT=error)
  if(error.ne.0) then
   write(*,1000)
   write(13,1000)
   write(*,1020)
   write(13,1020)
   call mcm_shutdown(2)
  endif
  !
  ! allocate neighbour list memory with new maximum size
  !
  allocate (mcm_g_nbrlist(mcm_g_maxnbr,mcm_np),mcm_g_contlist(mcm_g_maxcont,mcm_np),STAT=error)  
  if(error.ne.0) then
   write(*,1000)
   write(13,1000)
   write(*,1120)
   write(13,1120)
   call mcm_shutdown(2)
  endif
  !
 endif
 !
 ! loop over all particles
 !
 do i=1,mcm_np
 par(i)%nnbr = 0
 par(i)%ncont = 0
  !
  par(i)%g_nnbr = 0
  par(i)%g_ncont = 0
  !
 par(i)%mindist = 1.0e20_d
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
      !
      ! Check if cell contains real or ghost particles
      !
      if(id.gt.-1) then
       !
       ! Real particles
       !
	! loop over all particles in the cell
	do
	 ! exit if end of list has been reached, placed here so that is the 
	 ! cell is empty, the do loop is immediately exited
	    if(id.lt.1) exit
	 ! ignore if current particle is the same as the i partcle
	  if(id.ne.i) then
	   dist_sq = 0.0_d
	   do j=1,mcm_ndim
	    dist_sq = dist_sq + (par(id)%x(j)-par(i)%x(j))**2
       enddo
	   h_avg = factor * 0.5_d * (par(id)%h+par(i)%h)
       h_sq = h_avg*h_avg
       ! check if particle is within 2h of i particle
	   if(dist_sq.lt.h_sq) then
		! if they are of the same material then add to the neighbour list
		if(par(id)%mat.eq.par(i)%mat) then
		 ! check for minimum interparticle distance, used for timestep calculation option 2
		 par(i)%mindist = min(par(i)%mindist,dist_sq)
		 par(i)%nnbr = par(i)%nnbr + 1
		 ! check if exceeded max neighbour limit
		 if(par(i)%nnbr.le.mcm_maxnbr) then
		  mcm_nbrlist(par(i)%nnbr,i) = id
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
		 par(i)%ncont = par(i)%ncont +1
		 if(par(i)%ncont.le.mcm_maxcont) then
		  mcm_contlist(par(i)%ncont,i) = id
		 else
	      enough_mem = .false.
		 endif
		endif
	   endif
	  endif
	 ! get new id
	 id = par(id)%llpointer
	enddo
	   !
	  else
	   !
	   ! cell contains ghost particles
	   !
	   ! loop over all particles in the cell 
	   do
	    ! exit if end of list has been reached, placed here so that is the 
	    ! cell is empty, the do loop is immediately exited
	    if(id.gt.-1) exit
	    dist_sq = 0.0_d
	    id = abs(id)
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
	   !
	   ! End loop over ghost particles
	   !
	  endif
   enddo
  enddo
 enddo
 !
 ! if using kernel contact then augment neighbour list with contact list
 !
 if(mcm_contacttype.eq.0) then
  !
  if(enough_mem) then
   if((par(i)%nnbr + par(i)%ncont).le.mcm_maxnbr) then
    do j=1,par(i)%ncont
     if ((par(i)%nnbr + par(i)%ncont).eq.0) exit
     mcm_nbrlist(par(i)%nnbr+j,i) = mcm_contlist(j,i)
    enddo
   else
    enough_mem = .false.
   endif
  endif
  par(i)%nnbr = par(i)%nnbr + par(i)%ncont
  par(i)%mindist = min(par(i)%mindist,mindist_cp)
    !
    if(mcm_boundary) then
     if(enough_mem) then
      if((par(i)%g_nnbr + par(i)%g_ncont).le.mcm_g_maxnbr) then
       do j=1,par(i)%g_ncont
        if ((par(i)%g_nnbr + par(i)%g_ncont).eq.0) exit
        mcm_g_nbrlist(par(i)%g_nnbr+j,i) = mcm_g_contlist(j,i)
       enddo
      else
       enough_mem = .false.
      endif
     endif
     par(i)%g_nnbr = par(i)%g_nnbr + par(i)%g_ncont
    endif
 endif
 par(i)%mindist = sqrt(par(i)%mindist)
  endif  ! particle inside domain if  (par(i)%llpointer.gt.-1)
 enddo
 ! find new maxnbr
 mcm_maxnbr = 0
 mcm_maxcont= 0
 do i=1,mcm_np
 mcm_maxnbr = max(mcm_maxnbr,par(i)%nnbr)
 mcm_maxcont = max(mcm_maxcont,par(i)%ncont)
 enddo
 !
 if(mcm_boundary) then
  mcm_g_maxnbr = 0
  mcm_g_maxcont= 0
  do i=1,mcm_np
   mcm_g_maxnbr = max(mcm_g_maxnbr,par(i)%g_nnbr)
   mcm_g_maxcont = max(mcm_g_maxcont,par(i)%g_ncont)
  enddo
 endif
 !
 if(ip_flag.eq.1) then
  ! if required augment neighbour list with i particle
  do i=1,mcm_np
   par(i)%nnbr = par(i)%nnbr + 1
   mcm_maxnbr = max(mcm_maxnbr,par(i)%nnbr)
   if(par(i)%nnbr.le.mcm_maxnbr) then
    mcm_nbrlist(par(i)%nnbr,i) = i
   else
    enough_mem = .false.
   endif
  enddo 
 endif
 !
 if(enough_mem) exit
 !
 if(count.eq.2) then
  write(*,1000)
  write(13,1000)
  write(*,2000)
  write(13,2000)
  call mcm_shutdown(2)
 endif
 !
enddo
! Generate transducer neighbour arrays if required
!
if(mcm_num_transducer.gt.0) call mcm_transducer_nbr
!
! deallocate linked list memory
!
deallocate(mcm_llgrid,STAT=error)
if(error.ne.0) then
 write(*,1000)
 write(13,1000)
 write(*,1200)
 write(13,1200)
 call mcm_shutdown(2)
endif
!
return
!
1000 format(//5x,'Error in subroutine ll_neighbours')
1010 format(/5x,'Unexpected error when deallocating neighbour list array memory')
1020 format(/5x,'Unexpected error when deallocating ghost neighbour list array memory')
1100 format(/5x,'Error allocating neighbour list array memory')
1120 format(/5x,'Error allocating ghost neighbour list array memory')
1200 format(/5x,'Unexpected error when deallocating linked list array memory')
!
2000 format(/5x,'Error looped through neighbour search too many times',/, &
             5x,'This error should not occur, but is present to prevent an infinite loop')
!
end subroutine mcm_ll_neighbours