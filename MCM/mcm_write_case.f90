subroutine mcm_write_case
!*****************************************************************
!
!    Purpose : This module now writes out the case file
!
!  Called by : 
!
!       Date : 28-10-02
!
!     Errors : 
!
!      Notes : 
!
! Written by : Tom De Vuyst 
!
!*****************************************************************
USE mcm_database
!
IMPLICIT NONE 
integer :: mcm_rest,l,r,i
!
CHARACTER(LEN=25) :: mcm_file_case, mcm_file_geo, mcm_file_sig, mcm_file_dis, mcm_file_acc, &
                     mcm_file_vel, mcm_file_eps, mcm_file_rho, mcm_file_egy, mcm_file_snd,  &
					 mcm_file_mas, mcm_file_bnd, mcm_file_nbr, mcm_file_tmp
!
mcm_file_case = mcm_fileout(1:mcm_filelen(2))//".case"
mcm_file_geo  = mcm_fileout(1:mcm_filelen(2))//".geo"
mcm_file_sig  = mcm_fileout(1:mcm_filelen(2))//"_sig.ens"
mcm_file_dis  = mcm_fileout(1:mcm_filelen(2))//"_dis.ens"
mcm_file_vel  = mcm_fileout(1:mcm_filelen(2))//"_vel.ens"
mcm_file_acc  = mcm_fileout(1:mcm_filelen(2))//"_acc.ens"
mcm_file_eps  = mcm_fileout(1:mcm_filelen(2))//"_eps.ens"
mcm_file_rho  = mcm_fileout(1:mcm_filelen(2))//"_rho.ens"
mcm_file_egy  = mcm_fileout(1:mcm_filelen(2))//"_egy.ens"
mcm_file_snd  = mcm_fileout(1:mcm_filelen(2))//"_snd.ens"
mcm_file_mas  = mcm_fileout(1:mcm_filelen(2))//"_mas.ens"
mcm_file_bnd  = mcm_fileout(1:mcm_filelen(2))//"_bnd.ens"
mcm_file_nbr  = mcm_fileout(1:mcm_filelen(2))//"_nbr.ens"
mcm_file_tmp  = mcm_fileout(1:mcm_filelen(2))//"_tmp.ens"
!
OPEN(unit=32,file=mcm_file_case,status='unknown',form='formatted')
!
WRITE(32,'(A6)') "FORMAT"
WRITE(32,'(A15,/)') "type:   ensight"
WRITE(32,'(A8)') "GEOMETRY"
WRITE(32,'(A41,A26,/)') "model:                                    ",mcm_file_geo
WRITE(32,'(A9)') "VARIABLES"
WRITE(32,'(A41,A26)')   "scalar per node:         1 1 mass         ",mcm_file_mas
WRITE(32,'(A41,A26)')   "scalar per node:         1 1 boundary     ",mcm_file_bnd
WRITE(32,'(A41,A26)')   "scalar per node:         1 1 nnbr         ",mcm_file_nbr
WRITE(32,'(A41,A26)')   "vector per node:         1 1 displacement ",mcm_file_dis
WRITE(32,'(A41,A26)')   "vector per node:         1 1 velocity     ",mcm_file_vel
WRITE(32,'(A41,A26)')   "vector per node:         1 1 acceleration ",mcm_file_acc
WRITE(32,'(A41,A26)')   "scalar per element:      1 1 efps         ",mcm_file_eps
WRITE(32,'(A41,A26)')   "scalar per element:      1 1 density      ",mcm_file_rho
WRITE(32,'(A41,A26)')   "scalar per element:      1 1 energy       ",mcm_file_egy
WRITE(32,'(A41,A26)')   "scalar per element:      1 1 soundspeed   ",mcm_file_snd
WRITE(32,'(A41,A26)')   "scalar per element:      1 1 temperature  ",mcm_file_tmp
WRITE(32,'(A41,A26,/)') "tensor symm per element: 1 1 stress       ",mcm_file_sig
!
WRITE(32,'(A4)')    "TIME"
WRITE(32,'(A26)')   "time set:                1"
WRITE(32,'(A21,I5)')"number of steps:     ",  mcm_istate + 1
WRITE(32,'(A26)')   "filename start number:   1"
WRITE(32,'(A26)')   "filename increment:      1"
IF (mcm_istate.EQ.0) THEN
   WRITE(32,1030) "time values:  ", mcm_ptime
ELSE
   WRITE(32,1040) "time values:  ", (mcm_timenum(i), i=1,mcm_istate), mcm_ptime
ENDIF
!
WRITE(32,'(A4)')	"FILE"
WRITE(32,'(A26)')   "file set:                1"
WRITE(32,'(A26)')   "filename index:          1"
WRITE(32,'(A21,I5)')"number of steps:     ",  mcm_istate + 1
!
mcm_istate = mcm_istate + 1
mcm_timenum(mcm_istate) = mcm_ptime
!
close(unit=32)
!
1000 format("00",I1)
1010 format("0",I2)
1020 format("",I3)
1030 FORMAT (A14,E12.5/)
1040 FORMAT (A14,<mcm_istate>(E12.5/,14X),E12.5/)
!
END SUBROUTINE mcm_write_case 

