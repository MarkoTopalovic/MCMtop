    subroutine mcm_inp_symmetry_v2
!************************************************************************
!
!    Purpose: Read boundary plane information from version 2 input file
!
!  Called by: getinput
!
!       Date: 01-08-2002
!
!     Errors: 
!
!      Notes: Covers symmetry, periodic and neighbour limit boundaries
!
!************************************************************************
!
use mcm_database
!
implicit none

character(len=80) :: errmsg
character(len=80) :: textline
!
logical :: found
integer :: lcount, n, i
!
integer, dimension(6) :: mcm_sym_flag
real(kind=real_acc) :: mcm_sym_xmin, mcm_sym_xmax, mcm_sym_ymin, mcm_sym_ymax, mcm_sym_zmin, mcm_sym_zmax
!
errmsg='Error reading boundary plane flag card'
call mcm_gttxsg(textline,lcount)
!
select case(mcm_ndim)
 case(1)
  read(Unit=textline,fmt=100,err=500) mcm_sym_flag(1),mcm_sym_flag(2)
 case(2)
  read(Unit=textline,fmt=100,err=500) mcm_sym_flag(1),mcm_sym_flag(2),mcm_sym_flag(3),mcm_sym_flag(4)
 case(3)
  read(Unit=textline,fmt=100,err=500) mcm_sym_flag(1),mcm_sym_flag(2),mcm_sym_flag(3),mcm_sym_flag(4), &
                                      mcm_sym_flag(5),mcm_sym_flag(6)
end select
!
errmsg='Error reading boundary plane coordinate card'
call mcm_gttxsg(textline,lcount)
!
select case(mcm_ndim)
 case(1)
   read(Unit=textline,fmt=200,err=500) mcm_sym_xmin,mcm_sym_xmax
 case(2)
   read(Unit=textline,fmt=200,err=500) mcm_sym_xmin, mcm_sym_xmax, mcm_sym_ymin, mcm_sym_ymax 	   	 
 case(3)
   read(Unit=textline,fmt=200,err=500) mcm_sym_xmin, mcm_sym_xmax, mcm_sym_ymin, mcm_sym_ymax, &
                                       mcm_sym_zmin, mcm_sym_zmax
end select
!
! Write input to logfile
!
write(13,1000)
!
if(mcm_sym_flag(1).eq.1) write(13,1100) mcm_sym_xmin
if(mcm_sym_flag(2).eq.1) write(13,1200) mcm_sym_xmax
if(mcm_sym_flag(3).eq.1) write(13,1300) mcm_sym_ymin
if(mcm_sym_flag(4).eq.1) write(13,1400) mcm_sym_ymax
if(mcm_sym_flag(5).eq.1) write(13,1500) mcm_sym_zmin
if(mcm_sym_flag(6).eq.1) write(13,1600) mcm_sym_zmax
!
if(mcm_sym_flag(1).eq.2) write(13,2100) mcm_sym_xmin
if(mcm_sym_flag(2).eq.2) write(13,2200) mcm_sym_xmax
if(mcm_sym_flag(3).eq.2) write(13,2300) mcm_sym_ymin
if(mcm_sym_flag(4).eq.2) write(13,2400) mcm_sym_ymax
if(mcm_sym_flag(5).eq.2) write(13,2500) mcm_sym_zmin
if(mcm_sym_flag(6).eq.2) write(13,2600) mcm_sym_zmax
!
if(mcm_sym_flag(1).eq.3) write(13,2100) mcm_sym_xmin
if(mcm_sym_flag(2).eq.3) write(13,2200) mcm_sym_xmax
if(mcm_sym_flag(3).eq.3) write(13,2300) mcm_sym_ymin
if(mcm_sym_flag(4).eq.3) write(13,2400) mcm_sym_ymax
if(mcm_sym_flag(5).eq.3) write(13,2500) mcm_sym_zmin
if(mcm_sym_flag(6).eq.3) write(13,2600) mcm_sym_zmax
!
! Check for bc type
!
! Type 1: symmetry
found = .false.
do i=1,6
 if(mcm_sym_flag(i).eq.1) found = .true.
enddo
!
! Type 2: periodic
do i=1,5,2
 if(mcm_sym_flag(i).eq.2) then
  found = .true.
  if(mcm_sym_flag(i+1).ne.2) then
   write(*,300) i
   write(13,300) i
   call mcm_shutdown(2)
  endif
 endif
enddo
!
do i=2,6,2
 if(mcm_sym_flag(i).eq.2) then
  found = .true.
  if(mcm_sym_flag(i-1).ne.2) then
   write(*,310) i
   write(13,310) i
   call mcm_shutdown(2)
  endif
 endif
enddo
!
if(.not.found) mcm_boundary = .false.
!
! Store boundary codes
!
mcm_boundary_code(1,1) = mcm_sym_flag(1)
mcm_boundary_code(2,1) = mcm_sym_flag(2)
mcm_boundary_code(1,2) = mcm_sym_flag(3)
mcm_boundary_code(2,2) = mcm_sym_flag(4)
mcm_boundary_code(1,3) = mcm_sym_flag(5)
mcm_boundary_code(2,3) = mcm_sym_flag(6)
!
! Copy boundary coordinates into main array
!
mcm_boundary_x(1,1) = mcm_sym_xmin
mcm_boundary_x(2,1) = mcm_sym_xmax
mcm_boundary_x(1,2) = mcm_sym_ymin
mcm_boundary_x(2,2) = mcm_sym_ymax
mcm_boundary_x(1,3) = mcm_sym_zmin
mcm_boundary_x(2,3) = mcm_sym_zmax
!
! Set boundary type
!
do i=1,3
 mcm_boundary_type(i) = 1
 if(mcm_boundary_code(1,i).gt.0) mcm_boundary_type(i) = mcm_boundary_type(i) + 1
 if(mcm_boundary_code(2,i).gt.0) mcm_boundary_type(i) = mcm_boundary_type(i) + 2
enddo
!
return
!
500 call mcm_termin (textline,errmsg,lcount,0)
!
100 format(6i10)
200 format(6e10.0)
300 format(' Error in definition of periodic boundary conditions',/ &
            '  The corresponding coordiante maximum plane to plane number ',i1, &
			'  has not been defined as a periodic plane.')
310 format(' Error in definition of periodic boundary conditions',/ &
            '  The corresponding coordiante minimum plane to plane number ',i1, &
			'  has not been defined as a periodic plane.')
!
1000 format(//,'B O U N D A R Y   P L A N E S',//)
1100 format(' Xmin symmetry plane active at x coordinate:',e11.3,/)
1200 format(' Xmax symmetry plane active at x coordinate:',e11.3,/)
1300 format(' Ymin symmetry plane active at y coordinate:',e11.3,/)
1400 format(' Ymax symmetry plane active at y coordinate:',e11.3,/)
1500 format(' Zmin symmetry plane active at z coordinate:',e11.3,/)
1600 format(' Zmax symmetry plane active at z coordinate:',e11.3,/)
2100 format(' Xmin periodic plane active at x coordinate:',e11.3,/)
2200 format(' Xmax periodic plane active at x coordinate:',e11.3,/)
2300 format(' Ymin periodic plane active at y coordinate:',e11.3,/)
2400 format(' Ymax periodic plane active at y coordinate:',e11.3,/)
2500 format(' Zmin periodic plane active at z coordinate:',e11.3,/)
2600 format(' Zmax periodic plane active at z coordinate:',e11.3,/)
3100 format(' Xmin problem limit plane active at x coordinate:',e11.3,/)
3200 format(' Xmax problem limit plane active at x coordinate:',e11.3,/)
3300 format(' Ymin problem limit plane active at y coordinate:',e11.3,/)
3400 format(' Ymax problem limit plane active at y coordinate:',e11.3,/)
3500 format(' Zmin problem limit plane active at z coordinate:',e11.3,/)
3600 format(' Zmax problem limit plane active at z coordinate:',e11.3,/)

!
!
end subroutine mcm_inp_symmetry_v2