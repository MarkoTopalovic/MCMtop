SUBROUTINE mcm_write_case_2d
!*****************************************************************
!
!    Purpose : This module writes out all the files for all variables              
!              for the 2D co-located SPH simulation
!
!  Called by : state_output
!
!       Date : 22-07-99
!
!     Errors : 
!
!      Notes : 
!
!*****************************************************************
!
USE mcm_database
!
IMPLICIT NONE
!
INTEGER r,l,rest,i,mcm_begin,j,pos
!
REAL, DIMENSION(mcm_np) :: mcm_bound
REAL, DIMENSION(mcm_np) :: mcm_energy
REAL, DIMENSION(mcm_np) :: mcm_neighbours
REAL, DIMENSION(mcm_np) :: mcm_interpol
REAL, DIMENSION(mcm_ndim,mcm_np) :: mcm_disp
! 
! write the time in the time.txt file
! 
open(unit=33,file=mcm_file_time,status='unknown',position='append',form='formatted')
write(33,fmt='(E13.6)') mcm_ptime
close(unit=33)

if (mcm_istate.lt.10) then
	write(unit=mcm_extension,fmt=1000) mcm_istate
else if (mcm_istate.lt.100) then
	write(unit=mcm_extension,fmt=1010) mcm_istate
else
	write(unit=mcm_extension,fmt=1020) mcm_istate
end if
!________________________________________
!
! Open all the files and write the title
!________________________________________
!
mcm_filename_sig = mcm_fileout(1:mcm_filelen(2))//"_sig.ens"
mcm_filename_vel = mcm_fileout(1:mcm_filelen(2))//"_vel.ens"
mcm_filename_acc = mcm_fileout(1:mcm_filelen(2))//"_acc.ens"
mcm_filename_dis = mcm_fileout(1:mcm_filelen(2))//"_dis.ens"
mcm_filename_eps = mcm_fileout(1:mcm_filelen(2))//"_eps.ens"
mcm_filename_egy = mcm_fileout(1:mcm_filelen(2))//"_egy.ens"
mcm_filename_rho = mcm_fileout(1:mcm_filelen(2))//"_rho.ens"
mcm_filename_snd = mcm_fileout(1:mcm_filelen(2))//"_snd.ens"
mcm_filename_mas = mcm_fileout(1:mcm_filelen(2))//"_mas.ens"
mcm_filename_bnd = mcm_fileout(1:mcm_filelen(2))//"_bnd.ens"
mcm_filename_nbr = mcm_fileout(1:mcm_filelen(2))//"_nbr.ens"
mcm_filename_tmp = mcm_fileout(1:mcm_filelen(2))//"_tmp.ens"
!
IF (mcm_istate.EQ.1) THEN
   !
   open(unit=32,file=mcm_filename_dis,status='unknown',form='formatted')
   open(unit=35,file=mcm_filename_vel,status='unknown',form='formatted')
   open(unit=36,file=mcm_filename_acc,status='unknown',form='formatted')
   open(unit=37,file=mcm_filename_rho,status='unknown',form='formatted')
   open(unit=39,file=mcm_filename_sig,status='unknown',form='formatted')
   open(unit=45,file=mcm_filename_eps,status='unknown',form='formatted')
   open(unit=46,file=mcm_filename_egy,status='unknown',form='formatted')
   open(unit=47,file=mcm_filename_snd,status='unknown',form='formatted')
   open(unit=48,file=mcm_filename_mas,status='unknown',form='formatted')
   open(unit=49,file=mcm_filename_bnd,status='unknown',form='formatted')
   open(unit=50,file=mcm_filename_nbr,status='unknown',form='formatted')
   open(unit=51,file=mcm_filename_tmp,status='unknown',form='formatted')
   !
ELSE
   !
   open(unit=32,file=mcm_filename_dis,status='unknown',position='append',form='formatted')
   open(unit=35,file=mcm_filename_vel,status='unknown',position='append',form='formatted')
   open(unit=36,file=mcm_filename_acc,status='unknown',position='append',form='formatted')
   open(unit=37,file=mcm_filename_rho,status='unknown',position='append',form='formatted')
   open(unit=39,file=mcm_filename_sig,status='unknown',position='append',form='formatted')
   open(unit=45,file=mcm_filename_eps,status='unknown',position='append',form='formatted')
   open(unit=46,file=mcm_filename_egy,status='unknown',position='append',form='formatted')
   open(unit=47,file=mcm_filename_snd,status='unknown',position='append',form='formatted')
   open(unit=48,file=mcm_filename_mas,status='unknown',position='append',form='formatted')
   open(unit=49,file=mcm_filename_bnd,status='unknown',position='append',form='formatted')
   open(unit=50,file=mcm_filename_nbr,status='unknown',position='append',form='formatted')
   open(unit=51,file=mcm_filename_tmp,status='unknown',position='append',form='formatted')
   !
ENDIF
!
write(32,'(A15)')"BEGIN TIME STEP"
write(35,'(A15)')"BEGIN TIME STEP"
write(36,'(A15)')"BEGIN TIME STEP"
write(37,'(A15)')"BEGIN TIME STEP"
write(39,'(A15)')"BEGIN TIME STEP"
write(45,'(A15)')"BEGIN TIME STEP"
write(46,'(A15)')"BEGIN TIME STEP"
write(47,'(A15)')"BEGIN TIME STEP"
write(48,'(A15)')"BEGIN TIME STEP"
write(49,'(A15)')"BEGIN TIME STEP"
write(50,'(A15)')"BEGIN TIME STEP"
write(51,'(A15)')"BEGIN TIME STEP"
!
write(32,'(A12)')"displacement"
write(35,'(A8)')"velocity"
write(36,'(A12)')"acceleration"
write(37,'(A7)')"density"
write(39,'(A6)')"stress"
write(45,'(A4)')"efps"
write(46,'(A6)')"energy"
write(47,'(A14)')"speed of sound"
write(48,'(A4)')"mass"
write(49,'(A18)')"boundary particles"
write(50,'(A20)')"number of neighbours"
write(51,'(A11)')"temperature"
!_________________________
!
! Calculate Displacements
!_________________________
!
DO i = 1,mcm_np
DO j = 1,mcm_ndim
   !
   mcm_disp(j,i) = par(i)%x(j) - par(i)%xzero(j)
   !
ENDDO
ENDDO
!_______________________________________
!
! Write nodal variables in each file
!_______________________________________
! CHECK IF mcm_np IS ODD OR EVEN - TOM
!
IF (mod(mcm_np,2).EQ.0) THEN
   ! EVEN
   DO i = 1,mcm_np,2
      write(UNIT=32,fmt='(6(E12.5))') mcm_disp(1,i),mcm_disp(2,i),0.0,mcm_disp(1,i+1),mcm_disp(2,i+1),0.0
      write(UNIT=35,fmt='(6(E12.5))') par(i)%v(1),par(i)%v(2),0.0,par(i+1)%v(1),par(i+1)%v(2),0.0
      write(UNIT=36,fmt='(6(E12.5))') par(i)%a(1),par(i)%a(2),0.0,par(i+1)%a(1),par(i+1)%a(2),0.0
   ENDDO
   !
ELSE
   ! ODD
   DO i = 1,mcm_np-1,2
      write(UNIT=32,fmt='(6(E12.5))') mcm_disp(1,i),mcm_disp(2,i),0.0,mcm_disp(1,i+1),mcm_disp(2,i+1),0.0
      write(UNIT=35,fmt='(6(E12.5))') par(i)%v(1),par(i)%v(2),0.0,par(i+1)%v(1),par(i+1)%v(2),0.0
      write(UNIT=36,fmt='(6(E12.5))') par(i)%a(1),par(i)%a(2),0.0,par(i+1)%a(1),par(i+1)%a(2),0.0
   ENDDO
   write(UNIT=32,fmt='(3(E12.5))') mcm_disp(1,mcm_np),mcm_disp(2,mcm_np),0.0
   write(UNIT=35,fmt='(3(E12.5))') par(mcm_np)%v(1),par(mcm_np)%v(2),0.0
   write(UNIT=36,fmt='(3(E12.5))') par(mcm_np)%a(1),par(mcm_np)%a(2),0.0
   !
ENDIF
!
pos = 0
DO i=1,mcm_np
   !
   pos = pos + 1
   IF (pos.LT.6) THEN
      ! write efps - ensight - tom - 16-9-2002
      WRITE(48,101) par(i)%mass
      WRITE(49,101) par(i)%boundary
      WRITE(50,101) par(i)%nnbr
      !
   ELSE
      !
      WRITE(48,'(E12.5)') par(i)%mass
      WRITE(49,'(E12.5)') par(i)%boundary
      WRITE(50,'(E12.5)') par(i)%nnbr
      pos = 0
   ENDIF 
   ! 
ENDDO
!
IF (pos.NE.0) THEN
   !
   WRITE(48,'(/)')
   WRITE(49,'(/)')
   WRITE(50,'(/)')
   !
ENDIF
!______________________________________________
!
! write the element variables in each file
!______________________________________________
!
do j=1,mcm_nummat
	write(37,fmt='(A5,I8)') "part ",j
	write(39,fmt='(A5,I8)') "part ",j
	write(45,fmt='(A5,I8)') "part ",j
	write(46,fmt='(A5,I8)') "part ",j
	write(47,fmt='(A5,I8)') "part ",j
	write(51,fmt='(A5,I8)') "part ",j
!	write(32,fmt='(A22,I8)')"SPH particles in part ",j+nmmat
	write(37,'(A5)')"point"
	write(39,'(A5)')"point"
	write(45,'(A5)')"point"
	write(46,'(A5)')"point"
	write(47,'(A5)')"point"
	write(51,'(A5)')"point"
!	write(32,fmt='(I8)') material_MCM(j)
    pos = 0
	do i=1,mcm_np
	   IF (par(i)%mat.EQ.j) THEN
          write(UNIT=39,fmt='(6(E12.5))') par(i)%sigma(1,1), par(i)%sigma(2,2), &
                                          0.0              , par(i)%sigma(1,2), &
										  0.0              , 0.0
          !
		  pos = pos + 1
          IF (pos.LT.6) THEN
             ! write efps - ensight - tom - 16-9-2002
	         WRITE(37,101) par(i)%rho
	         WRITE(45,101) par(i)%efps
	         WRITE(46,101) par(i)%e / par(i)%mass
	         WRITE(47,101) par(i)%c
	         WRITE(51,101) par(i)%temper
	         !
          ELSE
             !
	         WRITE(37,'(E12.5)') par(i)%rho
	         WRITE(45,'(E12.5)') par(i)%efps
	         WRITE(46,'(E12.5)') par(i)%e / par(i)%mass
	         WRITE(47,'(E12.5)') par(i)%c
	         WRITE(51,'(E12.5)') par(i)%temper
	         pos = 0
          ENDIF 
	   ENDIF 
	end do
	!
	IF (pos.NE.0) THEN
	   !
	   write(37,'(/)')
	   write(45,'(/)')
	   write(46,'(/)')
	   write(47,'(/)')
	   write(51,'(/)')
	   !
	ENDIF
!	ind_MCM=ind_MCM+material_MCM(j)
end do
!
write(32,'(A13)')"END TIME STEP"
write(35,'(A13)')"END TIME STEP"
write(36,'(A13)')"END TIME STEP"
write(37,'(A13)')"END TIME STEP"
write(39,'(A13)')"END TIME STEP"
write(45,'(A13)')"END TIME STEP"
write(46,'(A13)')"END TIME STEP"
write(47,'(A13)')"END TIME STEP"
write(48,'(A13)')"END TIME STEP"
write(49,'(A13)')"END TIME STEP"
write(50,'(A13)')"END TIME STEP"
write(51,'(A13)')"END TIME STEP"
!
close(unit=32) ! disp
!close(unit=34)
close(unit=35) ! velocity
close(unit=36) ! acceleration
close(unit=37) ! density
!close(unit=38)
close(unit=39) ! stress
!close(unit=40)
!close(unit=41)
!close(unit=42)
!close(unit=55)
!close(unit=44)
close(unit=45) ! eps
close(unit=46) ! energy
close(unit=47) ! speed of sound
close(unit=48) ! mass
close(unit=49) ! boundary
close(unit=50) ! number of neighbours
close(unit=51) ! temperature
!close(unit=52)
!close(unit=53)
!close(unit=54)
!
!  istate_MCM = istate_MCM + 1
!
 101 FORMAT ('',E12.5,$)
1000 format("00",I1)
1010 format("0",I2)
1020 format("",I3)
2000 format(22(1X,E13.6))
2010 format(28(1X,E13.6))
3000 format(3(1X,E13.6))
!
end SUBROUTINE mcm_write_case_2d