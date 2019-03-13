subroutine mcm_sets10(cm,prop)
!************************************************************************
!
!    Purpose: sets initial material properties
!
!  Called by: matin
!
!       Date: 
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
real(kind=real_acc) :: cm(*), prop(48)
integer :: i
!
do i=1,48
 prop(i)=cm(i)
enddo
!
if (prop(4).eq.0.0) prop(4)=-1.e30
if (prop(7).eq.0.0) prop(7)=1.000000001
!
do i=1,48
 cm(i)=prop(i)
enddo
!
do i=2,16
 if(prop(i+15).lt.prop(i+16)) cm(16)=i
enddo
!
end subroutine mcm_sets10
