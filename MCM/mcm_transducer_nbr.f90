subroutine mcm_transducer_nbr
!************************************************************************
!
!    Purpose: find all transducer neighbours from linked list
!
!  Called by: mcm_ll_neighbours
!
!       Date: 12-08-2005
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
integer :: i,j,k,l,m,n
integer :: id, boxcrd(3),looplim(2,3),error
real(kind=real_acc) :: dist_sq,h_sq,mindist_cp,h_avg,factor
!
factor = 2.0_d
!
if(allocated(mcm_tr_nbr)) then
 deallocate(mcm_tr_nbr,STAT=error)
 if(error.ne.0) then
  write(*,1000)
  write(13,1000)
  write(*,1010)
  write(13,1010)
  call mcm_shutdown(2)
 endif
endif
!
allocate (mcm_tr_nbr(mcm_maxnbr+10,mcm_num_transducer),STAT=error)  
if(error.ne.0) then
 write(*,1000)
 write(13,1000)
 write(*,1100)
 write(13,1100)
 call mcm_shutdown(2)
endif
!
do i=1,3
 boxcrd(i) = 1
 looplim(1,i) = 1
 looplim(2,i) = 1
enddo
!
! loop over all transducers
!
do i=1,mcm_num_transducer
 mcm_num_tr_nbr(i) = 0
 !
 ! identify which cell the transducer lies in, and prevent loop over grid cells
 ! that do not exist
 !
 do j=1,mcm_ndim
  boxcrd(j)=int(mcm_gridsize(j)*(mcm_trans_x(j,i)-mcm_gridmin(j))) + 1
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
	 if(id.eq.0) exit
	  dist_sq = 0.0_d
	  do j=1,mcm_ndim
	   dist_sq = dist_sq + (par(id)%x(j)-mcm_trans_x(j,i))**2
      enddo
	  h_avg = factor * par(id)%h
      h_sq = h_avg*h_avg
      ! check if particle is within 2h of i particle
	  if(dist_sq.lt.h_sq) then
	   ! if they are of the correct material then add to the neighbour list
	   if(mcm_trmat(i).eq.par(id)%mat) then
	   mcm_num_tr_nbr(i) = mcm_num_tr_nbr(i) + 1
	   ! check if exceeded max neighbour limit
	   if(mcm_num_tr_nbr(i).le.mcm_maxnbr+10) then
		mcm_tr_nbr(mcm_num_tr_nbr(i),i) = id
	   endif
      endif
	 endif
	 ! get new id
	 id = par(id)%llpointer
	enddo
   enddo
  enddo
 enddo
enddo
!
return
!
1000 format(//5x,'Error in subroutine transducer_nbr')
1010 format(/5x,'Unexpected error when deallocating array memory')
1100 format(/5x,'Error allocating transducer neighbour array memory')
!
end subroutine mcm_transducer_nbr
