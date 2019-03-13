subroutine WriteVTK
!************************************************************************
!
!    Purpose: vrite VTK files for Paraview 
!
!  Called by: 
!
!       Date: 0
!
!     Errors: 
!
!      Notes: 
!
!************************************************************************
!
use mcm_database
!
implicit none
!
INTEGER r,l,rest,i,mcm_begin,j,pos
!
REAL, DIMENSION(mcm_np) :: mcm_bound
REAL, DIMENSION(mcm_np) :: mcm_energy
REAL, DIMENSION(mcm_np) :: mcm_neighbours
REAL, DIMENSION(mcm_np) :: mcm_interpol
REAL, DIMENSION(mcm_ndim,mcm_np) :: mcm_disp

REAL VMstress
INTEGER mcm_np1

!
! Open all the files and write the title
!________________________________________
!

CHARACTER(len=5)   :: str, str1
character(len=25) :: extension

!integer charParNumb,charParNumbDouble
character(len=25) :: charParNumb,charParNumbDouble
WRITE(str,'(I5)') mcm_timestep
str1=ADJUSTL(str)



if (mcm_istate.lt.10) then
	write(unit=extension,fmt=1000) mcm_istate
else if (mcm_istate.lt.100) then
	write(unit=extension,fmt=1010) mcm_istate
else
	write(unit=extension,fmt=1020) mcm_istate
end if
1000 format(I1,'.vtk')
1010 format(I2,'.vtk')
1020 format(I3,'.vtk')

if (mcm_np.lt.10) then
	write(unit=charParNumb,fmt=2000) mcm_np	
else if (mcm_np.lt.100) then
		write(unit=charParNumb,fmt=2010) mcm_np
else if (mcm_np.lt.1000) then
		write(unit=charParNumb,fmt=2020) mcm_np
else if (mcm_np.lt.10000) then
		write(unit=charParNumb,fmt=2030) mcm_np
else
		write(unit=charParNumb,fmt=2040) mcm_np	
end if

if (2*mcm_np.lt.10) then
	write(unit=charParNumbDouble,fmt=2000) 2*mcm_np
else if (2*mcm_np.lt.100) then
	write(unit=charParNumbDouble,fmt=2010) 2*mcm_np
else if (2*mcm_np.lt.1000) then
	write(unit=charParNumbDouble,fmt=2020) 2*mcm_np
else if (2*mcm_np.lt.10000) then
	write(unit=charParNumbDouble,fmt=2030) 2*mcm_np
else
	write(unit=charParNumbDouble,fmt=2040) 2*mcm_np
end if

2000 format(I1)
2010 format(I2)
2020 format(I3)
2030 format(I4)
2040 format(I5)



mcm_filename_vtk = mcm_fileout(1:mcm_filelen(2))//extension


!
   open(unit=88,file=mcm_filename_vtk,status='unknown',form='formatted')
!_________________________
!
! Calculate Displacements
!
DO i = 1,mcm_np
DO j = 1,mcm_ndim
   !
   mcm_disp(j,i) = par(i)%x(j) - par(i)%xzero(j)
   !
ENDDO
ENDDO

! ovde update-ovati broj cvorova !!!
mcm_np1 = mcm_np
!charParNumb = mcm_np
!charParNumbDouble = 2*mcm_np

WRITE(88, 81)
WRITE(88, 82)
WRITE(88, 83)
WRITE(88, 84)

!WRITE(88, 881) charParNumb
WRITE(88,*)"POINTS ",charParNumb," double"
DO i = 1,mcm_np1
 WRITE(88,882) par(i)%x(1), par(i)%x(2), par(i)%x(3)
ENDDO

!WRITE(88, 883)charParNumb, charParNumbDouble
!883 FORMAT ("CELLS ",3hNumb," ",3hNumb)
WRITE(88,*)"CELLS ",charParNumb," ",charParNumbDouble
DO i = 1,mcm_np1
 WRITE(88,884) I-1
ENDDO

WRITE(88, 885) mcm_np1
DO i = 1,mcm_np1
 WRITE(88,886)
ENDDO

WRITE(88, 887) mcm_np1
WRITE(88,*)"VECTORS Velocities double"

DO i = 1,mcm_np1
 WRITE(88,882) par(i)%v(1),par(i)%v(2),par(i)%v(3)
ENDDO

!WRITE(88,888)
!WRITE(88,889)
!DO i = 1,mcm_np1
! WRITE(88,890) par(i)%e
!ENDDO

WRITE(88,891)
WRITE(88,889)
DO i = 1,mcm_np1
 WRITE(88,890) par(i)%p
ENDDO

WRITE(88,892)
WRITE(88,889)
DO i = 1,mcm_np1
 WRITE(88,890) par(i)%sigma(1,1)
ENDDO

WRITE(88,893)
WRITE(88,889)
DO i = 1,mcm_np1
 WRITE(88,890) par(i)%sigma(2,2)
ENDDO

WRITE(88,894)
WRITE(88,889)
DO i = 1,mcm_np1
 WRITE(88,890) par(i)%sigma(3,3)
ENDDO

WRITE(88,895)
WRITE(88,889)
DO i = 1,mcm_np1
 WRITE(88,890) par(i)%sigma(1,2)
ENDDO

WRITE(88,896)
WRITE(88,889)
DO i = 1,mcm_np1
 WRITE(88,890) par(i)%sigma(1,3)
ENDDO

WRITE(88,897)
WRITE(88,889)
DO i = 1,mcm_np1
 WRITE(88,890) par(i)%sigma(2,3)
ENDDO

WRITE(88,898)
WRITE(88,889)
DO i = 1,mcm_np1
 VMStress = SQRT(0.5* &
					((par(i)%sigma(1,1)-par(i)%sigma(2,2))*(par(i)%sigma(1,1)-par(i)%sigma(2,2))+ &
					(par(i)%sigma(2,2)-par(i)%sigma(3,3))*(par(i)%sigma(2,2)-par(i)%sigma(3,3))+ &
					(par(i)%sigma(3,3)-par(i)%sigma(1,1))*(par(i)%sigma(3,3)-par(i)%sigma(1,1))+ &
					6*(par(i)%sigma(1,2)*par(i)%sigma(1,2)+par(i)%sigma(2,3)*par(i)%sigma(2,3)+par(i)%sigma(1,3)*par(i)%sigma(1,3))	))

 WRITE(88,890) VMStress
ENDDO

CLOSE(UNIT=88)

  81 FORMAT ("# vtk DataFile Version 3.0")
  82 FORMAT ("TEST")
  83 FORMAT ("ASCII")
  84 FORMAT ("DATASET UNSTRUCTURED_GRID")
   ! ovde update-ovati broj cvorova !!!
 881 FORMAT ("POINTS ",3hNumb," double")
 882 FORMAT (F15.8,F15.8,F15.8)
 ! ovde update-ovati broj cvorova !!!
883 FORMAT ("CELLS ",3hNumb," ",3hNumb)
 884 FORMAT ("1 ",I6)
 885 FORMAT ("CELL_TYPES ",I5)
 886 FORMAT ("1")
 887 FORMAT ("POINT_DATA ",I5)
 888 FORMAT ("SCALARS Energy double 1")
 889 FORMAT ("LOOKUP_TABLE default")
 890 FORMAT (E24.5)
 891 FORMAT ("SCALARS Pressure double 1")
 892 FORMAT ("SCALARS Stress_X double 1")
 893 FORMAT ("SCALARS Stress_Y double 1")
 894 FORMAT ("SCALARS Stress_Z double 1")
 895 FORMAT ("SCALARS Shear_Stress_XY double 1")
 896 FORMAT ("SCALARS Shear_Stress_XZ double 1")
 897 FORMAT ("SCALARS Shear_Stress_YZ double 1")
 898 FORMAT ("SCALARS VM_Stress double 1")

!
end subroutine WriteVTK 















subroutine WriteVTKmat(matID)
!************************************************************************
!
!    Purpose: write VTK files for Paraview 
!
!  Called by: 
!
!       Date: 0
!
!     Errors: 
!
!      Notes: 
!
!************************************************************************
!
use mcm_database
!
implicit none
!
INTEGER r,l,rest,i,mcm_begin,j,pos,matID,NoOfMatParticles,cellcounter
!
REAL, DIMENSION(mcm_np) :: mcm_bound
REAL, DIMENSION(mcm_np) :: mcm_energy
REAL, DIMENSION(mcm_np) :: mcm_neighbours
REAL, DIMENSION(mcm_np) :: mcm_interpol
REAL, DIMENSION(mcm_ndim,mcm_np) :: mcm_disp

REAL VMstress
INTEGER mcm_np1

!
! Open all the files and write the title
!________________________________________
!

CHARACTER(len=5)   :: str, str1
character(len=25) :: extension, charParNumb,charParNumbDouble
character(len=1) :: prefix
WRITE(str,'(I5)') mcm_timestep
str1=ADJUSTL(str)
write(unit=prefix,fmt=100) matID
100 format(I1)
if (mcm_istate.lt.10) then
	write(unit=extension,fmt=1000) mcm_istate
else if (mcm_istate.lt.100) then
	write(unit=extension,fmt=1010) mcm_istate
else
	write(unit=extension,fmt=1020) mcm_istate
end if
1000 format(I1,'.vtk')
1010 format(I2,'.vtk')
1020 format(I3,'.vtk')

mcm_filename_vtk = prefix//mcm_fileout(1:mcm_filelen(2))//extension


!
   open(unit=88,file=mcm_filename_vtk,status='unknown',form='formatted')
   NoOfMatParticles = 0
   DO i = 1,mcm_np
    if (par(i).mat.eq.matID) then
    NoOfMatParticles =NoOfMatParticles +1
    end if

ENDDO


if (NoOfMatParticles.lt.10) then
	    write(unit=charParNumb,fmt=4000) NoOfMatParticles
else if (NoOfMatParticles.lt.100) then
		write(unit=charParNumb,fmt=4010) NoOfMatParticles	
else if (NoOfMatParticles.lt.1000) then
		write(unit=charParNumb,fmt=4020) NoOfMatParticles	
else if (NoOfMatParticles.lt.10000) then
		write(unit=charParNumb,fmt=4030) NoOfMatParticles
else
		write(unit=charParNumb,fmt=4040) NoOfMatParticles
end if

if (2*NoOfMatParticles.lt.10) then
	    write(unit=charParNumbDouble,fmt=4000) 2*NoOfMatParticles
else if (2*NoOfMatParticles.lt.100) then
	    write(unit=charParNumbDouble,fmt=4010) 2*NoOfMatParticles
else if (2*NoOfMatParticles.lt.1000) then
	    write(unit=charParNumbDouble,fmt=4020) 2*NoOfMatParticles
else if (2*NoOfMatParticles.lt.10000) then
	    write(unit=charParNumbDouble,fmt=4030) 2*NoOfMatParticles
else		
	    write(unit=charParNumbDouble,fmt=4040) 2*NoOfMatParticles
end if

4000 format(I1)
4010 format(I2)
4020 format(I3)
4030 format(I4)
4040 format(I5)



   
!_________________________
!
! Calculate Displacements
!
DO i = 1,mcm_np
DO j = 1,mcm_ndim
   !
   if (par(i).mat.eq.matID) then
   mcm_disp(j,i) = par(i)%x(j) - par(i)%xzero(j)
   end if
   !
ENDDO
ENDDO

! ovde update-ovati broj cvorova !!!
mcm_np1 = mcm_np

WRITE(88, 81)
WRITE(88, 82)
WRITE(88, 83)
WRITE(88, 84)

!WRITE(88, 881) charParNumb
WRITE(88,*)"POINTS ",charParNumb," double"
DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 WRITE(88,882) par(i)%x(1), par(i)%x(2), par(i)%x(3)
 end if
ENDDO
cellcounter = 0
!WRITE(88, 883) charParNumb, charParNumbDouble
WRITE(88,*)"CELLS ",charParNumb," ",charParNumbDouble
DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 WRITE(88,884) cellcounter
 cellcounter = cellcounter+1
 end if
ENDDO

WRITE(88, 885) NoOfMatParticles
DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 WRITE(88,886)
 end if
ENDDO

WRITE(88, 887) NoOfMatParticles
WRITE(88,*)"VECTORS Velocities double"

DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 WRITE(88,882) par(i)%v(1),par(i)%v(2),par(i)%v(3)
 end if
ENDDO

WRITE(88,891)
WRITE(88,889)
DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 WRITE(88,890) par(i)%p
 end if
ENDDO

WRITE(88,892)
WRITE(88,889)
DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 WRITE(88,890) par(i)%sigma(1,1)
 end if
ENDDO

WRITE(88,893)
WRITE(88,889)
DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 WRITE(88,890) par(i)%sigma(2,2)
 end if
ENDDO

WRITE(88,894)
WRITE(88,889)
DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 WRITE(88,890) par(i)%sigma(3,3)
 end if
ENDDO

WRITE(88,895)
WRITE(88,889)
DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 WRITE(88,890) par(i)%sigma(1,2)
 end if
ENDDO

WRITE(88,896)
WRITE(88,889)
DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 WRITE(88,890) par(i)%sigma(1,3)
 end if
ENDDO

WRITE(88,897)
WRITE(88,889)
DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 WRITE(88,890) par(i)%sigma(2,3)
 end if
ENDDO

WRITE(88,898)
WRITE(88,889)
DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 VMStress = SQRT(0.5* &
					((par(i)%sigma(1,1)-par(i)%sigma(2,2))*(par(i)%sigma(1,1)-par(i)%sigma(2,2))+ &
					(par(i)%sigma(2,2)-par(i)%sigma(3,3))*(par(i)%sigma(2,2)-par(i)%sigma(3,3))+ &
					(par(i)%sigma(3,3)-par(i)%sigma(1,1))*(par(i)%sigma(3,3)-par(i)%sigma(1,1))+ &
					6*(par(i)%sigma(1,2)*par(i)%sigma(1,2)+par(i)%sigma(2,3)*par(i)%sigma(2,3)+par(i)%sigma(1,3)*par(i)%sigma(1,3))	))

 WRITE(88,890) VMStress
 end if
ENDDO

WRITE(88,899)
WRITE(88,889)
DO i = 1,mcm_np1
if (par(i).mat.eq.matID) then
 WRITE(88,900) par(i)%ptype
 end if
ENDDO


CLOSE(UNIT=88)

  81 FORMAT ("# vtk DataFile Version 3.0")
  82 FORMAT ("TEST")
  83 FORMAT ("ASCII")
  84 FORMAT ("DATASET UNSTRUCTURED_GRID")
  ! ovde update-ovati broj cvorova !!!
 881 FORMAT ("POINTS ",3hNumb," double")
 882 FORMAT (F15.8,F15.8,F15.8)
 ! ovde update-ovati broj cvorova !!!
883 FORMAT ("CELLS ",3hNumb," ",3hNumb)
 884 FORMAT ("1 ",I6)
 885 FORMAT ("CELL_TYPES ",I5)
 886 FORMAT ("1")
 887 FORMAT ("POINT_DATA ",I5)
 888 FORMAT ("SCALARS Energy double 1")
 889 FORMAT ("LOOKUP_TABLE default")
 890 FORMAT (E24.5)
 891 FORMAT ("SCALARS Pressure double 1")
 892 FORMAT ("SCALARS Stress_X double 1")
 893 FORMAT ("SCALARS Stress_Y double 1")
 894 FORMAT ("SCALARS Stress_Z double 1")
 895 FORMAT ("SCALARS Shear_Stress_XY double 1")
 896 FORMAT ("SCALARS Shear_Stress_XZ double 1")
 897 FORMAT ("SCALARS Shear_Stress_YZ double 1")
 898 FORMAT ("SCALARS VM_Stress double 1")
 899 FORMAT ("SCALARS PTYPE double 1")
 900 FORMAT (I5)
!
end subroutine WriteVTKmat