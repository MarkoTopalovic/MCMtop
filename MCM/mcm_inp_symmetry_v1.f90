    subroutine mcm_inp_symmetry_v1
!************************************************************************
!
!    Purpose: Read symmetry plane information from version 1 input file
!
!  Called by: getinput
!
!       Date: 01-08-2002
!
!     Errors: 
!
!      Notes: 
!
!************************************************************************
!
use mcm_database
!
implicit none

character(len=80) :: errmsg
character(len=80) :: textline
!
integer :: lcount, n, i, j
integer :: xbmin,xbmax,ybmin,ybmax,zbmin,zbmax
!
integer, dimension(6) :: mcm_sym_flag
real(kind=real_acc) :: mcm_sym_xmin, mcm_sym_xmax, mcm_sym_ymin, mcm_sym_ymax, mcm_sym_zmin, mcm_sym_zmax
!
errmsg='Error reading Symetry plans cards'
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
errmsg='Error reading Reflection plans cards'
call mcm_gttxsg(textline,lcount)
!
select case(mcm_ndim)
 case(1)
  read(Unit=textline,fmt=100,err=500) xbmin,xbmax
 case(2)
  read(Unit=textline,fmt=100,err=500) xbmin,xbmax,ybmin,ybmax   	 
 case(3)
  read(Unit=textline,fmt=100,err=500) xbmin,xbmax,ybmin,ybmax,zbmin,zbmax
end select
!
write(*,*) 'Reflective_Plans',xbmin, xbmax, ybmin, ybmax, zbmin, zbmax
write(*,*) 'Domain_Limits',mcm_sym_xmin, mcm_sym_xmax, mcm_sym_ymin, mcm_sym_ymax, &
                           mcm_sym_zmin, mcm_sym_zmax
!
! convert to new data structure
!
do i=1,3
 do j=1,2
  mcm_boundary_code(j,i) = 0
 enddo
enddo
!
if(xbmin.eq.0) mcm_boundary_code(1,1) = 1
if(xbmax.eq.0) mcm_boundary_code(2,1) = 1
if(ybmin.eq.0) mcm_boundary_code(1,2) = 1
if(ybmax.eq.0) mcm_boundary_code(2,2) = 1
if(zbmin.eq.0) mcm_boundary_code(1,3) = 1
if(zbmax.eq.0) mcm_boundary_code(2,3) = 1
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
100 format(6i5)
200 format(6e10.0)
!
!
end subroutine mcm_inp_symmetry_v1