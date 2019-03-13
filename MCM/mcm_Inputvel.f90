subroutine mcm_inputvel
!************************************************************************
!
!    Purpose: Read initial velocity cards from input file
!
!  Called by: getinput
!
!       Date: 01-08-2002
!
!     Errors: Errors in format or content of input file cards
!             The individual errors are described where they are checked for
!
!      Notes: Able to handle 1D, 2D and 3D input files
!
!************************************************************************
!
use mcm_database
!
implicit none
!
character*80 errmsg
character*80 textline
integer :: nid,igenflag,lcount
integer :: i,j
!
! Initialise all velocities to zero
!
do i=1,mcm_np
 do j=1,mcm_ndim
  par(i)%v(j) = 0.0_d
 enddo
enddo
!
! If initial velocities are not definied in input file then return
!
if(mcm_init_v_opt.ne.1) return
!
! Read initial velocities from input file
!
errmsg='Error reading initial velocity cards'
call mcm_gttxsg(textline,lcount)
select case (mcm_ndim)
 case(1)
  read(unit=textline,fmt=100,err=500) nid, par(mcm_svp)%v(1)
 case(2)
  read(unit=textline,fmt=100,err=500) nid, par(mcm_svp)%v(1), par(mcm_svp)%v(2)
 case(3)
  read(unit=textline,fmt=100,err=500) nid, par(mcm_svp)%v(1), par(mcm_svp)%v(2), par(mcm_svp)%v(3)
end select
!
! error termination if first velocity is not for first velocity node
if(nid.ne.mcm_svp) then
 write(*,1000)
 write(13,1000)
 write(*,2000) mcm_svp
 write(13,2000) mcm_svp
 call mcm_shutdown(2)
endif
!
!  now loop over remaining velocity particles, allow for gaps in the data to be 
!  automaticlly filled in, provided the velocity values each side of the 
!  gap are the same. igenflag is used to sign whether to read in the
!  next line (=0) or auto generate(=1)
!
igenflag=0
do i=mcm_svp+1,mcm_evp
 if(igenflag.eq.0) then
  call mcm_gttxsg(textline,lcount)
  select case (mcm_ndim)
   case(1)
    read(unit=textline,fmt=100,err=500) nid, par(i)%v(1)
   case(2)
    read(unit=textline,fmt=100,err=500) nid, par(i)%v(1), par(i)%v(2)
   case(3)
    read(unit=textline,fmt=100,err=500) nid, par(i)%v(1), par(i)%v(2), par(i)%v(3)
  end select
  ! if a block of nodes with the same velocity is being specified by
  ! not including them individually in the input file, but specifying the first
  ! and last nodes with the velocity
  if(nid.gt.i) then
   ! if the velocity is for a node which does not exist, error
   if(nid.gt.mcm_np) then
    write(*,1000)
	write(13,1000)
    write(*,1010)
    write(13,1010)
    call mcm_termin (textline,errmsg,lcount,1)
   endif
   ! if the 1 velocity is not constant for all, error
   if(par(i)%v(1).ne.par(i-1)%v(1)) then
    write(*,1000)
	write(13,1000)
    write(*,1020)
    write(13,1020)
    call mcm_termin (textline,errmsg,lcount,1)
   endif
   ! if the 2 velocity is not constant for all, error
   if(par(i)%v(2).ne.par(i-1)%v(2)) then
    write(*,1000)
	write(13,1000)
    write(*,1030)
    write(13,1030)
    call mcm_termin (textline,errmsg,lcount,1)
   endif
   ! if the 3 velocity is not constant for all, error
   if(par(i)%v(3).ne.par(i-1)%v(3)) then
    write(*,1000)
	write(13,1000)
    write(*,1040)
    write(13,1040)
    call mcm_termin (textline,errmsg,lcount,1)
   endif
   igenflag=1
  else
   ! if the nodes have not been specified in increasing order, error
   if(nid.lt.i) then
    write(*,1000)
	write(13,1000)
    write(*,1050)
    write(13,1050)
    call mcm_termin (textline,errmsg,lcount,1)
   endif
  endif
 else
  ! if incrementing then assign correct velocity to node
  do j=1,mcm_ndim
   par(i)%v(j)=par(i-1)%v(j)
  enddo
  if(i.eq.nid) igenflag=0
 endif
enddo
!
! write out velocities to log file
!
write(13,3000)
select case (mcm_ndim)
case (1)
  write(13,3010)
  do i=mcm_svp,mcm_evp
   write(13,3040) i, par(i)%v(1)
  enddo
case (2)
  write(13,3020)
  do i=mcm_svp,mcm_evp
   write(13,3040) i, par(i)%v(1), par(i)%v(2)
  enddo
case (3)
  write(13,3030)
  do i=mcm_svp,mcm_evp
   write(13,3040) i, par(i)%v(1), par(i)%v(2), par(i)%v(3)
  enddo
end select
!
return
!
500 call mcm_termin (textline,errmsg,lcount,1)
!
100 format(i8,3e10.0)
!
1000 format(4x,'Error in subroutine inputvel')
1010 format(4x,'Nodal id exceeds number of particles')
1020 format(4x,'x velocity for auto generation is not constant')
1030 format(4x,'y velocity for auto generation is not constant')
1040 format(4x,'z velocity for auto generation is not constant')
1050 format(4x,'Nodal id is less than previous card')
!
2000 format(4x,'Error, first nodal velocity card is not for node ',i5)
!
3000 format(//,'I N I T I A L   P A R T I C A L   V E L O C I T I E S')
3010 format(//,' Node id      x velocity'/)
3020 format(//,' Node id      x velocity      y velocity'/)
3030 format(//,' Node id      x velocity      y velocity      z velocity'/)
3040 format (i8,3e16.8)
!
end