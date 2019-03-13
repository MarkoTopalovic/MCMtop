SUBROUTINE mcm_write_tensor(unit,array,dim1,dim2,number,begin)
!*****************************************************************
!
!    Purpose : Writes [dim1,dim2] component of tensor array to 
!              file
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
!*****************************************************************use mcm_database
!
use mcm_database
!
implicit none
!
INTEGER :: i,j
REAL(kind=real_acc), DIMENSION(3,3,6) :: array
INTEGER, INTENT(IN) :: unit,number,begin,dim1,dim2
!
do i=1,number
 if(abs(array(dim1,dim2,i)).lt.1.0E-50_d) array(dim1,dim2,i) = 0.0_d
enddo
!
write(unit,100) (array(dim1,dim2,i),i=1,number)
!
100 format(6(e12.5))
!
end subroutine mcm_write_tensor


