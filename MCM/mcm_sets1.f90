subroutine mcm_sets1(cm,prop)
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
integer:: i,j
real(kind=real_acc) :: cm(5,*),q1,q2,q3,prop(48)
!
prop(1)=cm(1,1)
prop(2)=cm(1,2)
prop(3)=cm(1,3)
prop(4)=cm(1,4)
prop(5)=cm(1,5)
prop(6)=cm(2,1)
prop(7)=cm(3,1)
do i=2,5
 do j=1,5
  cm(i,j)=0.0
 enddo
enddo
q1=cm(1,1)*cm(1,2)/((1.0+cm(1,2))*(1.0-2.0*cm(1,2)))
q2=cm(1,1)*0.5/(1.0+cm(1,2))
q3=q1+2.0*q2
cm(2,1)=q3
cm(3,2)=q3
cm(4,3)=q3
cm(5,4)=q2
cm(2,2)=q1
cm(2,3)=q1
cm(3,1)=q1
cm(3,3)=q1
cm(4,1)=q1
cm(4,2)=q1
cm(2,6)=cm(1,6)
cm(1,6)=cm(1,1)*cm(1,4)/(cm(1,1)-cm(1,4))
cm(1,7)=prop(6)
cm(2,7)=prop(7)
!
end subroutine mcm_sets1
