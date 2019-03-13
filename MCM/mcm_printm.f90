subroutine mcm_printm(n,prop)
!***********************************************************************
!
!    Purpose: print out material properties
!
!  Called by: input
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
integer :: i,j,n
!
real(kind=real_acc) :: prop(48)
!
! if first material then write header information to the log file
if (n.eq.1) then
 write (13,310)
 write (13,330)
 write (13,350)
 write (13,360) 1,0.06,1.5
endif
!
write (13,370) (mcm_mat(n)%head(1,j),j=1,12)
write (13,380) n, mcm_mat(n)%model, mcm_mat(n)%eos, &
               mcm_mat(n)%visc_type, mcm_mat(n)%rho, &
               mcm_mat(n)%av_q, mcm_mat(n)%av_l
!
! write initial total material mass and h if given
if(mcm_massopt.eq.0) then
 write(13,384) mcm_mat(n)%mass
endif
!
if(mcm_init_h_opt.eq.0) then
 write(13,385) mcm_mat(n)%h
endif
!
write(13,386) mcm_mat(n)%rho_min, mcm_mat(n)%rho_max
!
! write out the material specific input data to log file
!
select case (mcm_mat(n)%model)
!
case (1)
  !
  !  model = 1    l i n e a r   i s o t r o p i c
  !
  write (13,390) (prop(i),i=1,2)
!
case (3)
  !
  !  model = 3    e l a s t o p l a s t i c (von mises yield criterion)
  !
  write (13,480) (prop(i),i=1,5)
!
case (9)
  !
  !   model = 9    n u l l   m a t e r i a l
  !
  write (13,540) (mcm_mat(n)%strinput(i),i=1,2)
!
case (10)
  !
  !   model = 10   h y d r o d y n a m i c   e l a s t i c  -  p l a s t i c
  !
  write (13,550) (prop(i),i=1,7),prop(9),(prop(i),i=17,48)
!
end select
!
!     write out equation-of-state data when applicable
!     and assign the default value for the relative volume
!
if(mcm_mat(n)%eos.gt.0) write (13,650) (mcm_mat(n)%head(2,j),j=1,12)
!
select case (mcm_mat(n)%eos)
!
case (4)
  if (mcm_mat(n)%eosinput(8).eq.0.0) mcm_mat(n)%eosinput(8)=1.0_d
  write (13,690) (mcm_mat(n)%eosinput(i),i=1,8)
!
case (13)
  if (mcm_mat(n)%eosinput(4).eq.0.0) mcm_mat(n)%eosinput(4)=1.0_d
  write (13,700) (mcm_mat(n)%eosinput(i),i=1,3)
!
case (28)
  if (mcm_mat(n)%eosinput(4).eq.0.0) mcm_mat(n)%eosinput(4)=1.0_d
  write (13,705) (mcm_mat(n)%eosinput(i),i=1,2)
!
case (41)
  if (mcm_mat(n)%eosinput(8).eq.0.0) mcm_mat(n)%eosinput(8)=1.0_d
  write (13,710) (mcm_mat(n)%eosinput(i),i=1,8)
!
end select
return
!
310 format(    //' m a t e r i a l   d e f i n i t i o n '///, &
                 ' material models                       '/, &
                 '     eq.1   isotropic elastic          '/, &
                 '     eq.3   kinematic/isotropic elastic-plastic '/, &
                 '     eq.9   fluid material              '/, &
                 '     eq.10  isotropic elastic-plastic hydrodynamic'/, &
                 '     eq.15  johnson/cook elastic-plastic '/)
330 format(&
                 ' equation-of-state models             '/, &
                 '     eq.4   gruneisen                 '/, &
                 '     eq.13  perfect gas               '/, &
				 '     eq.28  murnaghan (quasi-incompr.)'/, &
				 '     eq.41  mie-gruneisen             '/)
350 format(&
                 ' bulk viscosity models                '/, &
                 '     eq.1   standard                  '////)
360 format(&
                 ' default viscosities                  '/, &
                 '     bulk viscosity type...................',i5,/, &
                 '     linear bulk viscosity coefficient.....',1pe12.4,/, &
                 '     quadratic bulk viscosity coefficient..',1pe12.4,////)
370 format(/1x,12a6/)
380 format(' material constants set number ........ ',i5//, &
       10x,' material model ............. ',i5/, &
       10x,' equation-of-state model .... ',i5/, &
       10x,' bulk viscosity model ....... ',i5//, &
        5x,'den .............................. =', 1pe12.4/, &
        5x,'quadratic bulk viscosity ......... =', 1pe12.4/, &
        5x,'linear bulk viscosity ............ =', 1pe12.4/)
384 format(5x,'Total material mass............... =', 1pe12.4)
385 format(5x,'Initial h for material............ =', 1pe12.4)
386 format(5x,'Min density limit factor for mat.. =', 1pe12.4/, &
           5x,'Max density limit factor for mat.. =', 1pe12.4/)
390 format(&
        5x,'e ................................ =', 1pe12.4/, &
        5x,'vnu .............................. =', 1pe12.4///)
480 format(&
        5x,'e ................................ =', 1pe12.4/, &
        5x,'vnu .............................. =', 1pe12.4/, &
        5x,'yield ............................ =', 1pe12.4/, &
        5x,'e (harden) ....................... =', 1pe12.4/, &
        5x,'hardening parmeter ............... =', 1pe12.4///)
540 format(&
        5x,'pressure cutoff .................. =', 1pe12.4/, &
        5x,'viscosity coefficient ............ =', 1pe12.4 )
550 format(&
        5x,'g ................................ =', 1pe12.4/, &
        5x,'yield ............................ =', 1pe12.4/, &
        5x,'plastic hardening modulus ........ =', 1pe12.4/, &
        5x,'pressure cutoff .................. =', 1pe12.4/, &
        5x,'pressure hardening coefficient a1  =', 1pe12.4/, &
        5x,'pressure hardening coefficient a2  =', 1pe12.4/, &
        5x,'spall type ....................... =', 1pe12.4/, &
        5x,'eff. plastic strain at failure ... =', 1pe12.4/, &
        5x,'    eq.0.0: no failure'/, &
        5x,'effective plastic strain ......... =',8(1pe8.1,1x),/, &
        5x,'effective plastic strain ......... =',8(1pe8.1,1x),/, &
        5x,'effective stress ................. =',8(1pe8.1,1x),/, &
        5x,'effective stress ................. =',8(1pe8.1,1x),///)
650 format(///,1x,12a6,//)
690 format(&
        5x,'c  ............................... =', 1pe12.4/, &
        5x,'s1 ............................... =', 1pe12.4/, &
        5x,'s2 ............................... =', 1pe12.4/, &
        5x,'s3 ............................... =', 1pe12.4/, &
        5x,'g0 ............................... =', 1pe12.4/, &
        5x,'a  ............................... =', 1pe12.4/, &
        5x,'e0 ............................... =', 1pe12.4/, &
        5x,'initial relative volume .......... =', 1pe12.4///)
700 format(&
        5x,'perfect gass constant ............ =', 1pe12.4/, &
        5x,'initial energy ................... =', 1pe12.4/, &
        5x,'initial relative volume .......... =', 1pe12.4/)
!
705 format(&
        5x,'B ................................ =', 1pe12.4/, &
        5x,'Gamma ............................ =', 1pe12.4/)
!
710 format(&
        5x,'c  ............................... =', 1pe12.4/, &
        5x,'s  ............................... =', 1pe12.4/, &
        5x,'cv ............................... =', 1pe12.4/, &
        5x,'t0 ............................... =', 1pe12.4/, &
        5x,'g0 ............................... =', 1pe12.4/, &
        5x,'initial relative volume .......... =', 1pe12.4/, &
        5x,'beta ............................. =', 1pe12.4/, &
        5x,'e0 not used....................... =', 1pe12.4///)
!
end subroutine mcm_printm
