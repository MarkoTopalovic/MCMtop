module mcm_math_func
!******************************************************************************
!
! Purpose: Contains functions that are used by more than one section of the code
!
!    Date: 05-08-2002
!
!******************************************************************************
!
use mcm_database
!
implicit none
!
contains
 !
 function inverse(matrix)
 !
 ! invert an ndim x ndim matrix
 !
 implicit none
 !
 integer :: i
 real(kind=real_acc) :: determinant,invdet
 real(kind=real_acc) :: matrix(3,3),inverse(3,3),subdet(3,3)
 !
 select case (mcm_ndim)
   !
   !--------------------------------------------------------------
  case(1)
   inverse(1,1) = 1.0_d/matrix(1,1)
   !
   !--------------------------------------------------------------
  case(2)
   determinant = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)
   ! if matrix is singular then return the identity matrix
   if(abs(determinant).lt.1.0e-08_d) then
    inverse(1,1) = 1.0_d
    inverse(2,2) = 1.0_d
    inverse(1,2) = 0.0_d
    inverse(2,1) = 0.0_d
   else
    invdet = 1.0_d/determinant
    inverse(1,1) =  matrix(2,2)*invdet
    inverse(2,2) =  matrix(1,1)*invdet
    inverse(1,2) = -matrix(1,2)*invdet
    inverse(2,1) = -matrix(2,1)*invdet
   endif
   !
   !--------------------------------------------------------------
  case(3)
   subdet(1,1) = matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)
   subdet(1,2) = matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)
   subdet(1,3) = matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1)
   subdet(2,1) = matrix(1,2)*matrix(3,3) - matrix(1,3)*matrix(3,2)
   subdet(2,2) = matrix(1,1)*matrix(3,3) - matrix(1,3)*matrix(3,1)
   subdet(2,3) = matrix(1,1)*matrix(3,2) - matrix(1,2)*matrix(3,1)
   subdet(3,1) = matrix(1,2)*matrix(2,3) - matrix(1,3)*matrix(2,2)
   subdet(3,2) = matrix(1,1)*matrix(2,3) - matrix(1,3)*matrix(2,1)
   subdet(3,3) = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)
   !
    determinant = matrix(1,1)*subdet(1,1) - matrix(1,2)*subdet(1,2) + matrix(1,3)*subdet(1,3)
   !
   if(abs(determinant).lt.1.0e-08_d) then
    inverse(1:3,1:3) = 0.0_d
    do i=1,3
     inverse(i,i) = 1.0_d
    enddo
   else
    invdet = 1.0_d/determinant
    inverse(1,1) =  subdet(1,1)*invdet
    inverse(1,2) = -subdet(2,1)*invdet
    inverse(1,3) =  subdet(3,1)*invdet
    inverse(2,1) = -subdet(1,2)*invdet
    inverse(2,2) =  subdet(2,2)*invdet
    inverse(2,3) = -subdet(3,2)*invdet
    inverse(3,1) =  subdet(1,3)*invdet
    inverse(3,2) = -subdet(2,3)*invdet
    inverse(3,3) =  subdet(3,3)*invdet
   endif
 end select
 !
 end function inverse
 !
 function determinant(matrix)
 !
 ! calculate the determinant of a 3 x 3 matrix
 !
 implicit none
 !
 integer :: i
 real(kind=real_acc) :: determinant
 real(kind=real_acc) :: matrix(3,3),subdet(3)
 !
 subdet(1) = matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)
 subdet(2) = matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)
 subdet(3) = matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1)
 !
 determinant = matrix(1,1)*subdet(1) - matrix(1,2)*subdet(2) + matrix(1,3)*subdet(3)
 !
 end function determinant
 !
 function inv3(matrix)
 !
 ! invert a 3 x 3 matrix
 !
 implicit none
 !
 integer :: i
 real(kind=real_acc) :: determinant,invdet
 real(kind=real_acc) :: matrix(3,3),inv3(3,3),subdet(3,3)
 !
 subdet(1,1) = matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)
 subdet(1,2) = matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)
 subdet(1,3) = matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1)
 subdet(2,1) = matrix(1,2)*matrix(3,3) - matrix(1,3)*matrix(3,2)
 subdet(2,2) = matrix(1,1)*matrix(3,3) - matrix(1,3)*matrix(3,1)
 subdet(2,3) = matrix(1,1)*matrix(3,2) - matrix(1,2)*matrix(3,1)
 subdet(3,1) = matrix(1,2)*matrix(2,3) - matrix(1,3)*matrix(2,2)
 subdet(3,2) = matrix(1,1)*matrix(2,3) - matrix(1,3)*matrix(2,1)
 subdet(3,3) = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)
 !
 determinant = matrix(1,1)*subdet(1,1) - matrix(1,2)*subdet(1,2) + matrix(1,3)*subdet(1,3)
 !
 if(abs(determinant).lt.1.0e-08_d) then
  inv3(1:3,1:3) = 0.0_d
  do i=1,3
   inv3(i,i) = 1.0_d
  enddo
 else
  invdet = 1.0_d/determinant
  inv3(1,1) =  subdet(1,1)*invdet
  inv3(1,2) = -subdet(2,1)*invdet
  inv3(1,3) =  subdet(3,1)*invdet
  inv3(2,1) = -subdet(1,2)*invdet
  inv3(2,2) =  subdet(2,2)*invdet
  inv3(2,3) = -subdet(3,2)*invdet
  inv3(3,1) =  subdet(1,3)*invdet
  inv3(3,2) = -subdet(2,3)*invdet
  inv3(3,3) =  subdet(3,3)*invdet
 endif
 !
 end function inv3
 !
end module mcm_math_func