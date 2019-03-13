SUBROUTINE mcm_contact_in
!************************************************************************
!
!    Purpose: reads and writes contact definitions
!
!  Called by: mcm_getinput
!
!       Date: 28-10-02
!
!     Errors: 
!
!      Notes: 
!
!************************************************************************                         
!
USE mcm_database
!
IMPLICIT NONE                              
!
CHARACTER (len=80):: mssg,txts
INTEGER :: i,i1,j,j1,lcount,last,last1,l,l1,k,k1
INTEGER :: thnodeid_temp(10)
!
WRITE (unit=13,fmt=200)
!
mssg =' error reading contact definitions'
CALL mcm_gttxsg (txts,lcount)
READ (unit=txts,fmt=100,err=400) mcm_ncontmats, mcm_cont_opt
WRITE (unit=13,fmt=100) mcm_ncontmats, mcm_cont_opt
IF (mcm_ncontmats.EQ.-1) THEN
   ! All materials in contact with eachother
   DO i = 1,mcm_nummat  
      !
      CALL mcm_gttxsg (txts,lcount)	
      READ  (unit=txts,fmt=110,err=400) mcm_k_cont(i), mcm_n_cont(i)
	  WRITE (unit=13, fmt=115)
      WRITE (unit=13, fmt=120) i, mcm_k_cont(i), mcm_n_cont(i)
      !
   ENDDO
   !
ELSE
   ! Only a selection of materials in contact
   WRITE (unit=13,fmt=300)
   WRITE (unit=13,fmt=310)
   !
ENDIF
!
RETURN
!
400 CALL mcm_termin (txts,mssg,lcount,1)
!
100 FORMAT(2i5)
101 FORMAT(2i5,i8)
110 FORMAT(2E10.5)
115 FORMAT(//5x,'Material       K           n      ')
120 FORMAT(I10,2E12.5)
200 FORMAT(    //' C O N T A C T   D E F I N I T I O N S '///)
300 FORMAT(//5x,'Error Reading Contact Definitions',i8)
310 FORMAT(//5x,'This option is not implemented',i8)
!
END SUBROUTINE mcm_contact_in
