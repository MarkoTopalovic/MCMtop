SUBROUTINE mcm_write_vector(unit,array,number,begin)
!*****************************************************************
!
!    Purpose : Allows to write the "number" values of the vector contained
!              in "array" (for example x,v or a) in the "unit" file 
!             
!
!  Called by : write_file and write_file_3d
!
!       Date : 07-08-2002
!
!     Errors : 
!
!      Notes :
!
!*****************************************************************
!
use mcm_database
!
implicit none
!
INTEGER :: i,j
REAL(kind=real_acc), DIMENSION(3,6) :: array
INTEGER, INTENT(IN) :: unit,number,begin
!
do i=1,number
 do j=1,3
  if(abs(array(j,i)).lt.1.0E-50_d) array(j,i) = 0.0_d
 enddo
enddo
!
select case(number)
 !
 case(1,2)
  write(unit,100) ( (array(j,i),j=1,3),i=1,number)
  !
 case(3,4)
 write(unit,100) ( (array(j,i),j=1,3),i=1,2)
 write(unit,100) ( (array(j,i),j=1,3),i=3,number)
  !
 case(5,6)
 write(unit,100) ( (array(j,i),j=1,3),i=1,2)
 write(unit,100) ( (array(j,i),j=1,3),i=3,4)
 write(unit,100) ( (array(j,i),j=1,3),i=5,number)
 !
end select
!
100 format(6(e12.5))
!
end subroutine mcm_write_vector