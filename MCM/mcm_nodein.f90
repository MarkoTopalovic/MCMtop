subroutine mcm_nodein
!************************************************************************
!
!    Purpose: Read nodal coordinate cards
!
!  Called by: getinput
!
!       Date: 31-07-2002
!
!     Errors: node has zero material, node card out of order.
!
!      Notes: Will handle 1D, 2D and 3D data plus co-located and non-colocated
!             cases.
!
!************************************************************************ 
!
use mcm_database
!
!	use mat_prop
!
implicit none 
integer :: i,j,lcount,n,error    
character(len=80) :: mssg
character(len=80) :: txts
!
real(kind=real_acc), dimension(:,:),allocatable :: pardata
!
! loop over all nodes
!
do j=1,mcm_np
 !
 !  read node card
 !
 call mcm_gttxsg(txts,lcount)
 mssg='Error in nodal poistion cards'
 select case (mcm_ndim)
  case(1)
   read(unit=txts,fmt=1100,err=500) n, par(j)%dispbc, par(j)%x(1), par(j)%mat
  case(2)
   read(unit=txts,fmt=1200,err=500) n, par(j)%dispbc, par(j)%x(1), par(j)%x(2), par(j)%mat
  case(3)
   read(unit=txts,fmt=1300,err=500) n, par(j)%dispbc, par(j)%x(1), par(j)%x(2), par(j)%x(3), par(j)%mat
 end select
 !
 !  check that node card is in order
 !
 if(j.ne.n) then
  write(*,3000) n,j
  write(13,3000) n,j
  call mcm_shutdown(2)
 endif
 !
 !  check that node has non-zero material id
 !
 if(par(j)%mat.eq.0) then
  write(*,3100) j
  write(13,3100) j
  call mcm_shutdown(2)
 endif
 !
enddo
!
! now write out nodal coordinates to log file
!
select case (mcm_disctype)
 case (0,1)
 ! colocated
  write(13,4000)
  select case (mcm_ndim)
   case(1) !1D
    write(13,4100)
    do j=1,mcm_np
     write(13,5100) j, par(j)%dispbc, par(j)%x(1), par(j)%mat
	enddo
   case(2) !2D
    write(13,4200)
    do j=1,mcm_np
     write(13,5200) j, par(j)%dispbc, par(j)%x(1), par(j)%x(2), par(j)%mat
	enddo
   case(3) !3D
    write(13,4300)
    do j=1,mcm_np
     write(13,5300) j, par(j)%dispbc, par(j)%x(1), par(j)%x(2), par(j)%x(3), par(j)%mat
	enddo
  end select
end select
!
! If required read in additional particle information cards
!
if(mcm_massopt.eq.1.or.mcm_init_h_opt.eq.1.or.mcm_init_rhoe.eq.1) then
 !
 allocate (pardata(4,mcm_np),STAT=error)  
  if(error.ne.0) then
   write(*,1000)
   write(13,1000)
   write(*,2000)
   write(13,2000)
   call mcm_shutdown(2)
  endif 
 do j=1,mcm_np
  !
  !  read node card
  !
  call mcm_gttxsg(txts,lcount)
  mssg='Error in additional particle information cards'
  read(unit=txts,fmt=1400,err=500) n, pardata(1,j),pardata(2,j),pardata(3,j),pardata(4,j)
  !
  !  check that node card is in order
  !
  if(j.ne.n) then
   write(*,3200) n,j
   write(13,3200) n,j
   call mcm_shutdown(2)
  endif
  !
 enddo
 !
 if(mcm_massopt.eq.1) then
  do j=1,mcm_np
   par(j)%mass = pardata(1,j)
  enddo
 endif
 !
 if(mcm_init_h_opt.eq.1) then
  do j=1,mcm_np
   par(j)%h = pardata(2,j)
  enddo
 endif
 !
 if(mcm_init_rhoe.eq.1) then
  do j=1,mcm_np
   par(j)%rho = pardata(3,j)
   par(j)%e = pardata(4,j)
  enddo
 endif
 !
 deallocate(pardata,STAT=error)
  if(error.ne.0) then
   write(*,1000)
   write(13,1000)
   write(*,2010)
   write(13,2010)
   call mcm_shutdown(2)
  endif
endif
!
return
!
500 call mcm_termin(txts,mssg,lcount,0)
!
1000 format(4x,'Error in subroutine nodein')
!
1100 format(i8,i5,e20.0,i7)
1200 format(i8,i5,2e20.0,i7)
1300 format(i8,i5,3e20.0,i7)
1400 format(i8,4e20.0)
!
2000 format(4x,'Error allocating temporary memory for additional particle data')
2010 format(4x,'Error deallocating temporary memory for additional particle data')
!
3000 format(4x,'Card for node ',i8,' is out of order. Should be ',i8)
3100 Format(4x,'Node ',i8,' has zero material type')
3200 format(4x,'Additional information card for node ',i8,' is out of order. Should be ',i8)
!
4000 format(//,'N O D E   I N P U T   C A R D S')
4010 format(/,'Stress Points')
4020 format(/,'Velocity Points')
!
4100 format(' Node id Node bc    x coordinate Node material',/)
4200 format(' Node id Node bc    x coordinate    y coordinate Node material',/)
4300 format(' Node id Node bc    x coordinate    y coordinate    z coordinate Node material',/)
!
5100 format(i8,i8,e16.8,i5)
5200 format(i8,i8,2e16.8,i5)
5300 format(i8,i8,3e16.8,i5)
!
end subroutine mcm_nodein
