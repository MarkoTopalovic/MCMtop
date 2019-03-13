SUBROUTINE mcm_kernel(w,xi,xj,havg)
!************************************************************************
!
!    Purpose: COntrol routine for calculation of SPH kernel function
!
!  Called by: 
!
!       Date: 08-08-2002
!
!     Errors: 
!
!      Notes:
!
!************************************************************************
!
use mcm_database
!
IMPLICIT NONE
!
REAL(kind=real_acc)    :: W, havg
real(kind=real_acc), dimension(3) :: xi, xj
!
SELECT CASE (mcm_krtype)
 !
 CASE (1) ! B-Spline                                                                                          
  CALL mcm_kernel1(w, xi, xj, havg)
  !
 CASE DEFAULT
  WRITE(*,1000)
  WRITE(13,1010) mcm_krtype
  !
END SELECT
!
RETURN
!
1000 FORMAT(//5x,'Error in subroutine kernel')
1010 FORMAT(/5x,'Kernel type ',I5,' does not exist')
!
END SUBROUTINE mcm_kernel
!
!
SUBROUTINE mcm_kernel1(W,xi,xj,havg)
!************************************************************************
!
!    Purpose: Calculate B-Spline kernel function
!
!  Called by: 
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
IMPLICIT NONE
!
INTEGER :: i, j, n
REAL(kind=real_acc) :: W, havg, const1, const2, rad, z
real(kind=real_acc), dimension(3) :: xi, xj, xij
!
SELECT CASE (mcm_ndim)
 !
 CASE (3)
  const1 = 1.0_d  / ( pi*havg*havg*havg )
 CASE (2)
  const1 = 10.0_d / ( 7.0_d*pi*havg*havg )
 CASE (1)
  const1 = 2.0_d  / ( 3.0_d*havg )
END SELECT
!
const2 = const1 / 4.0_d
!
rad = 0.0_d
!
DO n = 1,mcm_ndim
 xij(n) = xi(n) - xj(n)
 rad = rad + xij(n)**2
ENDDO
!
rad = SQRT(rad)
z = rad/havg
!
IF (z .lt. 1.0_d) THEN
 W = const1 * ( ( (0.75_d * z - 1.5_d) * z) * z + 1.0_d)
ELSEIF (z .lt. 2.0_d) THEN
 W = const2 * (2.0_d - z)**3
ELSE
 W = 0.0_d
ENDIF

end subroutine mcm_kernel1
