SUBROUTINE mcm_gradw(dwdx,dwdr,xi,xj,havg)
!************************************************************************
!
!    Purpose: Control routine for calculation of gradient of SPH kernel function
!
!  Called by: 
!
!       Date: 05-08-2002
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
REAL(kind=real_acc) :: dwdx(mcm_ndim),dwdr,havg,xi(3),xj(3)
!
SELECT CASE (mcm_krtype)
 ! 
 CASE (1) ! B-Spline                                                                                          
  CALL mcm_grad1(dwdx,dwdr,xi,xj,havg)
  !
 CASE DEFAULT
  WRITE(*,1000)
  WRITE(13,1010) mcm_krtype
  !
END SELECT
!
RETURN
!
1000 FORMAT(//5x,'Error in subroutine gradw')
1010 FORMAT(/5x,'Kernel type ',I5,' does not exist')
!
END SUBROUTINE mcm_gradw
!
!
SUBROUTINE mcm_grad1(dwdx,dwdr,xi,xj,havg)
!************************************************************************
!
!    Purpose: Calculate gradient of B-spline kernel
!
!  Called by: 
!
!       Date: 05-08-2002
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
INTEGER :: n, i, j
REAL(kind=real_acc)    :: dw, dwdx(mcm_ndim),dwdr, xi(3),xj(3)
REAL(kind=real_acc)    :: const1, const2, havg, xij(mcm_ndim), rad, z
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
const2 = const1 * 0.25_d
!
rad = 0.0_d
!
DO n = 1,mcm_ndim
   xij(n) = xi(n) - xj(n)
   rad = rad + xij(n)**2
ENDDO
rad = SQRT(rad)
z = rad/havg
!
IF (z .lt. 1.0_d) THEN
   dw = const1 * ( (9.0_d/4.0_d * z - 3.0_d) * z )
ELSEIF (z .lt. 2.0_d) THEN
   dw = -3.0_d * const2 * (2.0_d - z)**2
ELSE
   dw = 0.0_d
ENDIF
!
DO n = 1,mcm_ndim
 IF (rad .eq. 0.0_d) THEN
  dwdx(n) = 0.0_d
  dwdr = 0.0_d
 ELSE
  dwdx(n) = dw*xij(n) / (rad*havg)
  dwdr = dw/havg
 ENDIF
ENDDO
!
END
