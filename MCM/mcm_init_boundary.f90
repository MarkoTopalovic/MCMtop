subroutine mcm_init_boundary
!************************************************************************
!
!    Purpose: Set up values for boundary ghost particles
!
!  Called by: init_neighb, ghost_setup
!
!       Date: 10/2005
!
!     Errors: 
!
!      Notes: This routine is called if symmetry or periodic planes are
!             active. It requires that the linked list array is allocated.
!             The linked list array is allocated in routine mcm_setuplist
!             and deallocated at the end of mcm_ll_neighbours.
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer :: error, i, j, k, l, count, id, cell_count
integer :: xcount, ycount, zcount, temp, multiplier
!
! Set boundary multiplyers
!
select case(mcm_boundary_code(1,1))    ! x min boundary
 case(1) ! symmetry
   mcm_xmin_mult = (/ -1.0_d,  1.0_d,  1.0_d /)
   !
   mcm_xmin_add  = (/  0.0_d,  0.0_d,  0.0_d /)
   mcm_xmin_add(1) = 2.0_d*mcm_boundary_x(1,1)
   !
 case(2) ! periodic
   mcm_xmin_mult = (/  1.0_d,  1.0_d,  1.0_d /)
   !
   mcm_xmin_add  = (/  0.0_d,  0.0_d,  0.0_d /)
   mcm_xmin_add(1) = -(mcm_boundary_x(2,1)-mcm_boundary_x(1,1))
   !
end select
!
select case(mcm_boundary_code(2,1))    ! x max boundary
 case(1) ! symmetry
   mcm_xmax_mult = (/ -1.0_d,  1.0_d,  1.0_d /)
   !
   mcm_xmax_add  = (/  0.0_d,  0.0_d,  0.0_d /)
   mcm_xmax_add(1) = 2.0_d*mcm_boundary_x(2,1)
   !
 case(2) ! periodic
   mcm_xmax_mult = (/  1.0_d,  1.0_d,  1.0_d /)
   !
   mcm_xmax_add  = (/  0.0_d,  0.0_d,  0.0_d /)
   mcm_xmax_add(1) = mcm_boundary_x(2,1)-mcm_boundary_x(1,1)
   !
end select
!
select case(mcm_boundary_code(1,2))    ! y min boundary
 case(1) ! symmetry
   mcm_ymin_mult = (/  1.0_d, -1.0_d,  1.0_d /)
   !
   mcm_ymin_add  = (/  0.0_d,  0.0_d,  0.0_d /)
   mcm_ymin_add(2) = 2.0_d*mcm_boundary_x(1,2)
   !
 case(2) ! periodic
   mcm_ymin_mult = (/  1.0_d,  1.0_d,  1.0_d /)
   !
   mcm_ymin_add  = (/  0.0_d,  0.0_d,  0.0_d /)
   mcm_ymin_add(2) = -(mcm_boundary_x(2,2)-mcm_boundary_x(1,2))
   !
end select
!
select case(mcm_boundary_code(2,2))    ! y max boundary
 case(1) ! symmetry
   mcm_ymax_mult = (/  1.0_d, -1.0_d,  1.0_d /)
   !
   mcm_ymax_add  = (/  0.0_d,  0.0_d,  0.0_d /)
   mcm_ymax_add(2) = 2.0_d*mcm_boundary_x(2,2)
   !
 case(2) ! periodic
   mcm_ymax_mult = (/  1.0_d,  1.0_d,  1.0_d /)
   !
   mcm_ymax_add  = (/  0.0_d,  0.0_d,  0.0_d /)
   mcm_ymax_add(2) = mcm_boundary_x(2,2)-mcm_boundary_x(1,2)
   !
end select
!
select case(mcm_boundary_code(1,3))    ! z min boundary
 case(1) ! symmetry
   mcm_zmin_mult = (/  1.0_d,  1.0_d, -1.0_d /)
   !
   mcm_zmin_add  = (/  0.0_d,  0.0_d,  0.0_d /)
   mcm_zmin_add(3) = 2.0_d*mcm_boundary_x(1,3)
   !
 case(2) ! periodic
   mcm_zmin_mult = (/  1.0_d,  1.0_d,  1.0_d /)
   !
   mcm_zmin_add  = (/  0.0_d,  0.0_d,  0.0_d /)
   mcm_zmin_add(3) = -(mcm_boundary_x(2,3)-mcm_boundary_x(1,3))
   !
end select
!
select case(mcm_boundary_code(2,3))    ! z max boundary
 case(1) ! symmetry
   mcm_zmax_mult = (/  1.0_d,  1.0_d, -1.0_d /)
   !
   mcm_zmax_add  = (/  0.0_d,  0.0_d,  0.0_d /)
   mcm_zmax_add(3) = 2.0_d*mcm_boundary_x(2,3)
   !
 case(2) ! periodic
   mcm_zmax_mult = (/  1.0_d,  1.0_d,  1.0_d /)
   !
   mcm_zmax_add  = (/  0.0_d,  0.0_d,  0.0_d /)
   mcm_zmax_add(3) = mcm_boundary_x(2,3)-mcm_boundary_x(1,3)
   !
end select
!
! Check if ghost particle array is allocated, if so then deallocate
!
if(allocated(gpar)) then
 deallocate (gpar,STAT=error)  
 if(error.ne.0) then
  write(*,1000)
  write(13,1000)
  write(*,1100)
  write(13,1100)
  call mcm_shutdown(2)
 endif
endif
!
! Now estimate how many ghost particles are required for the calculation
!
count = 0
!
! Loop over all linked list cells
!
do k=1,mcm_gridlim(3)
 do j=1,mcm_gridlim(2)
  do i = 1,mcm_gridlim(1)
   id = mcm_llgrid(i,j,k)  ! get pointer to first particle in cell
   if(id.gt.0) then        ! only look at cell if it contains a particle
    !
    ! Work out if the cell needs to be counted, and how many times
	!
	xcount = 0
	ycount = 0
	zcount = 0
	!
	select case(mcm_boundary_code(1,1))    ! x min boundary
	 case(1,2)
      if(i.eq.2) xcount = xcount + 1
	end select
	!
	select case(mcm_boundary_code(2,1))    ! x max boundary
	 case(1,2)
      if(i.eq.mcm_gridlim(1)-1) xcount = xcount + 1
	end select
	!
	select case(mcm_boundary_code(1,2))    ! y min boundary
	 case(1,2)
      if(j.eq.2) ycount = ycount + 1
	end select
	!
	select case(mcm_boundary_code(2,2))    ! y max boundary
	 case(1,2)
      if(j.eq.mcm_gridlim(2)-1) ycount = ycount + 1
	end select
	!
	select case(mcm_boundary_code(1,3))    ! z min boundary
	 case(1,2)
      if(k.eq.2) zcount = zcount + 1
	end select
	!
	select case(mcm_boundary_code(2,3))    ! z max boundary
	 case(1,2)
      if(k.eq.mcm_gridlim(3)-1) zcount = zcount + 1
	end select
	!
	temp = xcount + ycount*(xcount + 1)
    multiplier = temp + zcount*(temp + 1)
	!
	if(multiplier.gt.0) count = count + multiplier * cell_count(i,j,k)
	!
   endif
  enddo
 enddo
enddo
!
! count is now the number of ghost particles needed at the start of the simulation
!  for safety add 10% of the number of particles to the total to allow for changes
!
mcm_max_ngp = count+ int(0.1*mcm_np)
!
! Allocatet ghost particle memory
!
allocate (gpar(mcm_max_ngp),stat=error)
if(error.ne.0) then
 write(*,1000)
 write(13,1000)
 write(*,1200)
 write(13,1200)
 call mcm_shutdown(2)
endif 
!
write(*,2000) mcm_max_ngp
write(13,2000) mcm_max_ngp
!
!------------------------------------------------------------------------
!
1000 format(//5x,'Error in subroutine mcm_init_boundary')
1100 format(/5x,'Error deallocating ghost particle memory')
1200 format(/5x,'Error allocating ghost particle memory')
!
2000 format(/5x,'Allocated memory for ',i6,' ghost particles.')
!
end subroutine mcm_init_boundary
!
!
integer function cell_count(i,j,k)
!
! Function to return the number of particles in box i,j,k
!
use mcm_database
!
implicit none
!
integer :: id, count, i, j, k
!
count = 0
! get id of first particle in cell
id=mcm_llgrid(i,j,k)
! loop over all particles in the cell
do
 ! exit if end of list has been reached, placed here so that is the 
 ! cell is empty, the do loop is immediately exited
 if(id.eq.0) exit
 count = count + 1
 ! get new id
 id = par(id)%llpointer
enddo
!
cell_count = count
!
end function cell_count
