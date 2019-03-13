subroutine mcm_termin (txts,mssg,lcount,iprint)
!************************************************************************
!
!    Purpose: write error message for error in input file
!
!  Called by: input, matin
!
!       Date: 30-07-2002
!
!     Errors: 
!
!      Notes:
!
!************************************************************************
implicit none                               
character(len=80) :: txts,mssg
integer :: lcount, iprint
!
if (iprint.eq.0) then
 write (13,30) lcount,lcount,txts
 write ( *,30) lcount,lcount,txts
else
 write (13,20) lcount,lcount,txts,mssg
 write ( *,20) lcount,lcount,txts,mssg
endif
call mcm_shutdown(2)
return
20 format(/// &
      '     line number',i7,' contains improperly formatted data',// &
     ,'************************************************line#',i7 &
     ,/,a80,/ &
     ,'************************************************************',/, &
       a80,/)
30 format(/// &
      '     line number',i7,' contains improperly formatted data',// &
     ,'************************************************line#',i7 &
     ,/,a80,/ &
     ,'************************************************************',/)
end subroutine mcm_termin
