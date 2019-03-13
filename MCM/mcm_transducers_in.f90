subroutine mcm_transducers_in
!************************************************************************
!
!    Purpose: read pressure transducer coordinates from input file
!
!  Called by: getinput
!
!       Date: 01-08-2002
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
integer :: i,lcount
character(len=80) :: mssg
character(len=80) :: txts
!
do i=1,mcm_num_transducer
 !
 call mcm_gttxsg(txts,lcount)
 mssg='Error in pressure transducer cards'
 select case (mcm_ndim)
  case(1)
   read(unit=txts,fmt=1100,err=500) mcm_trans_x(1,i), mcm_trmat(i)
  case(2)
   read(unit=txts,fmt=1200,err=500) mcm_trans_x(1,i),mcm_trans_x(2,i), mcm_trmat(i)
  case(3)
   read(unit=txts,fmt=1300,err=500) mcm_trans_x(1,i),mcm_trans_x(2,i),mcm_trans_x(3,i), mcm_trmat(i)
 end select

enddo

return
!
500 call mcm_termin(txts,mssg,lcount,0)
!
1100 format(E10.0,I5)
1200 format(2E10.0,I5)
1300 format(3E10.0,I5)
!
end subroutine mcm_transducers_in