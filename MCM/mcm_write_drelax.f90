subroutine mcm_write_drelax
!************************************************************************
!
!    Purpose: Write out particle data from dynamic relaxation analysis
!
!  Called by: main               
!
!       Date: 06-07-2005
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
!
logical :: fexist
integer :: i
character(len=16) :: filetest
!
if(.not.mcm_drelax) return
!
!Open files
!
filetest=mcm_fileout(1:mcm_filelen(2))//'.drlx'
!
! delete old file if it exists then open new state plot file
!
inquire(file=filetest,exist=fexist)
if(fexist) then
 open(unit=101,file=filetest,status='old',form='formatted')
 close(unit=101,status='delete')
endif
open(unit=101,file=filetest,status='new',form='formatted')
!
write(101,100)
!
! Write particle coordinates
!
select case(mcm_ndim)
 case(1)
  do i=1,mcm_np
   write(101,11) i,par(i)%dispbc,par(i)%x(1),par(i)%mat
  enddo
 case(2)
  do i=1,mcm_np
   write(101,12) i,par(i)%dispbc,par(i)%x(1),par(i)%x(2),par(i)%mat
  enddo
 case(3)
  do i=1,mcm_np
   write(101,13) i,par(i)%dispbc,par(i)%x(1),par(i)%x(2),par(i)%x(3),par(i)%mat
  enddo
end select
!
! Write particle data
!
write(101,200)
!
do i=1,mcm_np
 write(101,20) i,par(i)%mass,par(i)%h,par(i)%rho,par(i)%e
enddo
!
11 format(i8,i5,e20.13,i7)
12 format(i8,i5,2e20.13,i7)
13 format(i8,i5,3e20.13,i7)
!
20 format(i8,4e20.13)
!
100 format('*',/ &
           '* Node Coordinates and Materials',/ &
		   '*')
200 format('*',/ &
           '* Node Data',/ &
		   '*     ID  Mass                Smoothing Length    Density             Internal Energy')
!
end subroutine mcm_write_drelax
