SUBROUTINE mcm_WRITE_SCALAR(unit,array,number,begin)
!*****************************************************************
!
!    Purpose : Allows to write the "number" values of "array"
!              in the "unit" file 
!             
!
!  Called by : write_file and write_file_3d
!
!       Date : 07-08-2002
!
!     Errors : 
!
!      Notes : Each time before to write, value if checked and if 
!              smaller than 1E-50 replaced by zero
!
!*****************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
INTEGER :: i
REAL(kind=real_acc), DIMENSION(6) :: array
INTEGER, INTENT(IN) :: unit,number,begin
!
do i=1,number
 if(abs(array(i)).lt.1.0E-50_d) array(i) = 0.0_d
enddo
!
write(unit,100) (array(i),i=1,number)
!
100 format(6(e12.5))
END SUBROUTINE mcm_write_scalar
