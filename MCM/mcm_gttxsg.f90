subroutine mcm_gttxsg(txts,lcount)
!************************************************************************
!
!    Purpose: read in a line of text from the input file
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
implicit none                               
!
character(len=80) :: txts,mssg
integer, save :: count
integer ::lcount
logical :: ok
mssg =' error reading input from file '
10   continue
read (unit=101,fmt=20,err=100) txts
count=count+1
if (txts(1:1).eq.'*'.or.txts(1:1).eq.'$') then
  go to 10
else
  lcount=count
  return
endif
100   call mcm_termin(txts,mssg,count,1)
20 format (a80)
end subroutine mcm_gttxsg
