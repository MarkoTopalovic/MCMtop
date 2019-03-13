      subroutine mcm_matin  
!************************************************************************
!
!    Purpose: reads material properties
!
!  Called by: input
!
!       Date: 30-07-02
!
!     Errors: Improper format
!
!      Notes: Correct for input versions 1 and 2
!
!************************************************************************
use mcm_database
!
implicit none
!
integer :: i,ifst,ilst,j,lcount,m,n
! ihq and qh are hourglass cofficients inherited from Dyna format. They are ignored.
integer :: ihq
real(kind=real_acc) :: qh
character(len=80) :: txts,mssg
!
real(kind=real_acc) :: prop(48)
!
do m=1,mcm_nummat
 mcm_mat(m)%model = 0
 mcm_mat(m)%eos = 0
enddo
!
!
! loop over each material and read in the necessary data
!
do m=1,mcm_nummat
 !
 !     read material control card
 !
 mssg =' error reading material control card'
 call mcm_gttxsg (txts,lcount)
 read (unit=txts,fmt=250,err=200) n,mcm_mat(n)%model,mcm_mat(n)%rho,mcm_mat(n)%eos, &
      ihq,qh,mcm_mat(n)%visc_type,mcm_mat(n)%av_q,mcm_mat(n)%av_l
 !
 ! set default artificial viscosity information if input file is zero
 !
 if(mcm_mat(n)%visc_type.eq.0) mcm_mat(n)%visc_type= 1
 if(mcm_mat(n)%av_l.eq.0.0) mcm_mat(n)%av_l=0.06_d
 if(mcm_mat(n)%av_q.eq.0.0) mcm_mat(n)%av_q=1.5_d
 !
 ! Read in MCM material options card
 !
 mssg =' error reading MCM material options card'
  call mcm_gttxsg (txts,lcount)
  read (unit=txts,fmt=260,err=200) mcm_mat(n)%mass,mcm_mat(n)%h,mcm_mat(n)%rho_min,mcm_mat(n)%rho_max
  if (mcm_mat(n)%rho_min.eq.0.0_d) mcm_mat(n)%rho_min=1.0E-10_d
  if (mcm_mat(n)%rho_max.eq.0.0_d) mcm_mat(n)%rho_max=1.0E+10_d
!
! if the supplied material number is outside the allowed range
!
 if (n.le.0.or.n.gt.mcm_nummat) then
  write(* ,9000)
  write(13,9000)
  call mcm_shutdown(2)
 endif
9000 format(///5x,'error in input - material number is out of range')
 !
 ! check to see if the material model type is supported by the code
 !
 select case (mcm_mat(n)%model)
 !
 !
 case(1,3,4)
   !  
   ! materials that do not require an equation of state, error if one is supplied
   if(mcm_mat(n)%eos.ne.0) then
    write(*,410)  n,mcm_mat(n)%model,mcm_mat(n)%eos
    write(13,410) n,mcm_mat(n)%model,mcm_mat(n)%eos
    call mcm_shutdown(2)
   endif
 !
 !
 case(9,10)
   !
   ! for materials that require an equation of state, check if it is supplied
   select case (mcm_mat(n)%eos)
    !
   case(0)
	 ! eos is not supplied, error
     write( *,400) n,mcm_mat(n)%model
     write(13,400) n,mcm_mat(n)%model
     call mcm_shutdown(2)
   !
   case(1,4,13,28,41)
     ! eos is supplied and supported
     continue
   !
   case default
   !
     9020 format(///5x,'error in input - eos type not implemented')
     write( *,9020)
     write(13,9020)
     call mcm_shutdown(2)
   end select
  !
  !
  case default
   !
   ! material number not recognised
   write(* ,9010)
   write(13,9010)
   call mcm_shutdown(2)
 end select
9010 format(///5x,'error in input - material type not implemented')
 !
 ! read material title card
 !
 call mcm_gttxsg (txts,lcount)
 read (unit=txts,fmt=230,err=200) (mcm_mat(n)%head(1,j),j=1,12)
 !
 !     read material properties
 !
 mssg =' error reading material property cards'
 !
 select case (mcm_mat(n)%model)
 !
 case(1,3,4,9)
   !
   ifst=1
   ilst=5
   do j=1,6
    call mcm_gttxsg (txts,lcount)
    read (unit=txts,fmt=220,err=200) (mcm_mat(n)%strinput(i),i=ifst,ilst)
    ifst=ifst+5
    ilst=ilst+5
   enddo
   !
 case(10)
   !
   ifst=1
   ilst=8
   do j=1,6
    call mcm_gttxsg (txts,lcount)
    read (unit=txts,fmt=220,err=200) (mcm_mat(n)%strinput(i),i=ifst,ilst)
    ifst=ifst+8
    ilst=ilst+8
   enddo
   !
 end select
 !
 ! set material constants
 !
 select case (mcm_mat(n)%model)
 case(1)
   call mcm_sets1  (mcm_mat(n)%strinput(1),prop)
 case(3)
   call mcm_sets3  (mcm_mat(n)%strinput(1),prop)
   case(4)
   call mcm_sets4  (mcm_mat(n)%strinput(1),prop)
 case(10)
   call mcm_sets10 (mcm_mat(n)%strinput(1),prop)
 end select
 !
 ! read eos data
 !
 if (mcm_mat(n)%eos.ne.0) then
  !
  ! read eos title card
  !
  mssg =' error reading eos property cards'
  call mcm_gttxsg (txts,lcount)
  read (unit=txts,fmt=230,err=200) (mcm_mat(n)%head(2,j),j=1,12)
  !
  !  read in equation of state constants
  !
  select case (mcm_mat(n)%eos)
  !
case(1)
ifst=1
   ilst=8
 do j=1,2
 call mcm_gttxsg (txts,lcount)
     read (unit=txts,fmt=220,err=200) (mcm_mat(n)%eosinput(i),i=ifst,ilst)
    ifst=ifst+8
	ilst=ilst+1
	 enddo

  case(4)
    call mcm_gttxsg (txts,lcount)
    read (unit=txts,fmt=220,err=200) (mcm_mat(n)%eosinput(i),i=1,8)
  !
  case(13)
    call mcm_gttxsg (txts,lcount)
    read (unit=txts,fmt=222,err=200) (mcm_mat(n)%eosinput(i),i=1,5)
  !
  case(28)
    call mcm_gttxsg (txts,lcount)
    read (unit=txts,fmt=222,err=200) (mcm_mat(n)%eosinput(i),i=1,2)
  !
  case(41)
    call mcm_gttxsg (txts,lcount)
    read (unit=txts,fmt=220,err=200) (mcm_mat(n)%eosinput(i),i=1,8)
  !
  end select
 endif	
 !
 ! write out material properties
 !
 call mcm_printm(n,prop)
 !
enddo
!
return
!
! error termination due to improperly formatted data
!
200 call mcm_termin (txts,mssg,lcount,1)
!
210 format(i5,e10.0,i5,2e10.0)
215 format(6e10.0)
220 format(8e10.0)
222 format(5e10.0)
230 format(12a6)
240 format(' material',i4,' density is zero-fatal')
250 format(2i5,e10.0,2i5,e10.0,i5,2e10.0)
260 format(4e10.0)
270 format(5e16.0)
400 format(//5x,'Error in material definition', &
           /10x,'material number:  ',i8, &
           /10x,'  material type:  ',i8, &
           //5x,'An equation of state is required for this material type', &
            /5x,'and no equation of state was defined.')
410 format(//5x,'Error in material definition', &
           /10x,'material number:  ',i8, &
           /10x,'  material type:  ',i8, &
           /10x,'       eos type:  ',i8, &
           //5x,'An equation of state is not permitted for this material', &
            /5x,'type but was defined')
!
end subroutine mcm_matin
