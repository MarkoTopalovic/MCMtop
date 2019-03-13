subroutine mcm_sets3(cm,prop)
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
!      Notes: Model 3 originally used routine sets1. This routine was created
!             when the effective plastic strain failure criterion was added
!             and only calulates the parameters used by model 3.
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
prop(6)=cm(1,6)
prop(7)=cm(2,6)
do i=2,5
 do j=1,6
  cm(i,j)=0.0
 enddo
enddo
q1=cm(1,1)*cm(1,2)/((1.0+cm(1,2))*(1.0-2.0*cm(1,2)))
q2=cm(1,1)*0.5/(1.0+cm(1,2))
q3=q1+2.0*q2
cm(2,1)=q3
cm(5,4)=q2
cm(2,4)=cm(1,1)*cm(1,4)/(cm(1,1)-cm(1,4))
cm(2,6)=prop(7)
!
end subroutine mcm_sets3
