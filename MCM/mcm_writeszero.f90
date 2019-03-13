subroutine mcm_writeszero
!************************************************************************
!
!    Purpose: Write out *.s000 file
!
!  Called by: stateplot
!
!       Date: 07-08-2002
!
!     Errors: 
!
!      Notes: If the *.s000 file already existes then it is deleted before
!             the new one is written. It is easier to write a complete new
!             file than to work out which bits require updating.
!
!************************************************************************
!
use mcm_database
!
implicit none
!
character(len=16) :: stfilename
logical :: fexist
integer :: i,coloc
!
coloc = 0
!
stfilename=mcm_fileout(1:mcm_filelen(2))//'.s000'
!
! delete old file if it exists then open new state plot file
!
inquire(file=stfilename,exist=fexist)
if(fexist) then
 open(unit=2,file=stfilename,status='old',form='formatted')
 close(unit=2,status='delete')
endif
open(unit=2,file=stfilename,status='new',form='formatted')
!
! write out main data
!
write(2,2000) mcm_title
write(2,2010) mcm_istate,mcm_ptime,mcm_np,mcm_nummat,mcm_ndim,coloc,mcm_nvelocp,mcm_nstressp
!
! write out variable names plus max/min data
! variable names are defined in initialisation subroutine in this file
!
write(2,2100) mcm_out_cols(1), mcm_out_cols(2) ! number of variables
do i=1,mcm_ndim
 write(2,2110) mcm_maxmin(1,i),mcm_maxmin(2,i)
enddo
do i=mcm_ndim+1,mcm_out_cols(1)+mcm_ndim
 write(2,2120) mcm_variablename(i)
 write(2,2110) mcm_maxmin(1,i),mcm_maxmin(2,i)
enddo
if(coloc.eq.1) then
 do i = mcm_out_cols(1)+mcm_ndim+1 , mcm_out_cols(1)+mcm_out_cols(2)+mcm_ndim
  write(2,2120) mcm_variablename(i)
  write(2,2110) mcm_maxmin(1,i),mcm_maxmin(2,i)
 enddo
endif
!
! close file
!
close(unit=2)
!
return
!
2000 format(a78)
2010 format(i3,e14.6,i7,i3,i3,i3,2i7)
2100 format(2i4)
2110 format(2e14.6)
2120 format(a40)
!
end subroutine mcm_writeszero