subroutine mcm_nodehist
!************************************************************************
!
!    Purpose: reads and writes node IDs for time history plots
!
!  Called by: matin
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
implicit none                               
!
character (len=80):: mssg,txts
integer:: i,i1,j,j1,lcount,last,last1,l,l1,k,k1
integer:: thnodeid_temp(10)
!
write (unit=13,fmt=200)
!
if (mcm_thnode.le.10) then
	last=0
	last1=mcm_thnode
else
	last=mcm_thnode/10
	last1=mod(mcm_thnode,10)
endif
!
do i=0, last
 mssg =' error reading nodes for time histry'
 call mcm_gttxsg (txts,lcount)
 if(i.lt.last) then	
  ! load node id's into temprary array	
  read (unit=txts,fmt=100,err=400) (thnodeid_temp(j), j=1,10)
  write (unit=13,fmt=100) (thnodeid_temp(k), k=1,10)
  ! copy temprary array into full array
  do l=1,10
   mcm_thnodeid(i+l)=thnodeid_temp(l)
   if(mcm_thnodeid(i+l).le.0.or.mcm_thnodeid(i+l).gt.mcm_np) then
	write( *,300) i+1
	write(13,300) i+1
   endif
  enddo
 else
  ! load node id's into temprary array 
  read (unit=txts,fmt=110,err=400) (thnodeid_temp(j1), j1=1,last1)
  write (unit=13,fmt=110) (thnodeid_temp(k1), k1=1,last1)
  ! copy temprary array into full array
  do l1=1,last1
   mcm_thnodeid(i+l1)=thnodeid_temp(l1)
   if(mcm_thnodeid(i+l1).le.0.or.mcm_thnodeid(i+l1).gt.mcm_np) then
	write( *,300) i+1
	write(13,300) i+1
   endif
  enddo
 endif
enddo
!
return
!
400 call mcm_termin (txts,mssg,lcount,1)
!
100 format(10i8)
110 format(10i8)
200 format(    //' P A R T I C L E S  S E L E C T E D  F O R  T I M E  &
			H I S T O R Y  P L O T S '///)
300 format(//5x,'Time history particle ID out of range, card:  ',i8)

end subroutine mcm_nodehist
