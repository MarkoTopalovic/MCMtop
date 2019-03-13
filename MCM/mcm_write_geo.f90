SUBROUTINE mcm_write_geo
!**********************************************************************
!
!    Purpose : This module writes out the .geo file
!
!  Called by : init_state_output
!
!       Date : 22-07-99
!
!     Errors : 
!
!      Notes :  The .geo file contain the coordinates of all nodes
!               and the repartition of differents parts ( that is to
!               say the different pieces )
!
!**********************************************************************
!
USE mcm_database
!
IMPLICIT NONE
!
LOGICAL mcm_fexist
INTEGER :: mcm_ind=1, i, j, k, n, num, mat
REAL(kind=real_acc) :: mcm_xnode,mcm_ynode,mcm_znode,mcm_zero=0
CHARACTER(LEN=25), EXTERNAL :: mcm_concatene
CHARACTER(LEN=25) :: mcm_file_write
INTEGER, DIMENSION(10) :: mcm_material
!
!
! if old time file exists then delete it
!
mcm_file_time = mcm_fileout(1:mcm_filelen(2))//"_time.txt"
inquire(file=mcm_file_time, exist=mcm_fexist)
!
if(mcm_fexist) then
   !
   open(unit=33, file=mcm_file_time, status='old', form='formatted')
   close(unit=33, status='delete')
   !
endif
!
mcm_file_write = mcm_fileout(1:mcm_filelen(2))//".geo"
open(unit=32, file=mcm_file_write, status='unknown', form='formatted') !,POSITION='append')
!
! The following lines must always be exactly the same
!
write(32,'(A15)') "BEGIN TIME STEP"
write(32,'(A6)')  "Title1"
write(32,'(A6)')  "Title2"
write(32,'(A13)') "node id given"
Write(32,'(A16)') "element id given"
write(32,'(A11)') "coordinates"
write(32,fmt='(I8)') mcm_np  !write number of nodes
!
!
! Next loop writes the coordinates of all nodes
!
!
! MCM
!
DO j = 1,mcm_np
   !
   mat = real(par(j)%mat)
   WRITE(32,fmt='(I8,E12.5,E12.5,E12.5)') j,par(j)%x(1),par(j)%x(2),par(j)%x(3)
   mcm_material(mat) = mcm_material(mat)+1
   !
END DO
!
! MCM
!
do j=1,mcm_nummat
	write(32,fmt='(A5,I8)') "part ",j
	write(32,fmt='(A22,I8)')"SPH particles in part ",j
	write(32,'(A5)')"point"
	write(32,fmt='(I8)') mcm_material(j)
	do i=1,mcm_np
	   IF (par(i)%mat.EQ.j) THEN
		write(32,fmt='(I8,I8)') i, i
	   ENDIF 
	end do
	mcm_ind = mcm_ind + mcm_material(j)
end do
!
write(32,'(A13)')"END TIME STEP"
!
close(unit=32)
!
RETURN
!
END SUBROUTINE mcm_write_geo
