subroutine mcm_initstateplot
!************************************************************************
!
!    Purpose: Initialise state plot variables - max and min plus names
!
!  Called by: initial
!
!       Date: 7-12-98
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
integer :: i
!
mcm_istate=1
do i=1,30
 mcm_maxmin(1,i)= 1.0e29_d
 mcm_maxmin(2,i)=-1.0e29_d
enddo
!
select case (mcm_disctype)
 ! colocated
 case(0,1)
  select case (mcm_ndim)
   case(1)
    mcm_out_cols(1) = 14
	mcm_out_cols(2) = 0
    mcm_variablename(2) ='Particle Material'
    mcm_variablename(3) ='x Velocity'
    mcm_variablename(4) ='x Acceleration'
    mcm_variablename(5) ='Density'
    mcm_variablename(6) ='Mass'
    mcm_variablename(7) ='Pressure'
    mcm_variablename(8) ='11 Stress'
    mcm_variablename(9) ='Effective Plastic Strain'
    mcm_variablename(10)='Specific Internal Energy'
    mcm_variablename(11)='Speed of Sound'
    mcm_variablename(12)='Particle Critical Timestep'
    mcm_variablename(13)='Temperature'
	mcm_variablename(14)='Number of Neighbours'
	mcm_variablename(15)='Smoothing Length'
   case(2)
    mcm_variablename(3) ='Particle Material'
    mcm_variablename(4) ='x Velocity'
    mcm_variablename(5) ='y Velocity'
    mcm_variablename(6) ='x Acceleration'
    mcm_variablename(7) ='y Acceleration'
    mcm_variablename(8) ='Density'
    mcm_variablename(9) ='Mass'
    mcm_variablename(10)='Pressure'
    mcm_variablename(11)='11 Stress'
    mcm_variablename(12)='22 Stress'
    mcm_variablename(13)='12 Stress'
	select case (mcm_axopt)
	 case(2)
      mcm_out_cols(1) = 19
	  mcm_out_cols(2) = 0
      mcm_variablename(14)='Effective Plastic Strain'
      mcm_variablename(15)='Specific Internal Energy'
      mcm_variablename(16)='Speed of Sound'
      mcm_variablename(17)='Particle Critical Timestep'
      mcm_variablename(18)='Temperature'
	  mcm_variablename(19)='Boundary Particles'
	  mcm_variablename(20)='Smoothing Length'
	  mcm_variablename(21)='Number of neighbours'
     case(4)
      mcm_out_cols(1) = 20
	  mcm_out_cols(2) = 0
	  mcm_variablename(14)='Hoop Stress'
      mcm_variablename(15)='Effective Plastic Strain'
      mcm_variablename(16)='Specific Internal Energy'
      mcm_variablename(17)='Speed of Sound'
      mcm_variablename(18)='Particle Critical Timestep'
      mcm_variablename(19)='Temperature'
	  mcm_variablename(20)='Boundary Particles'
	  mcm_variablename(21)='Smoothing Length'
	  mcm_variablename(22)='Number of neighbours'
    end select
    case(3)
	 write(*,1000)
	 write(*,1010)
	 write(13,1000)
	 write(13,1010)
	 call mcm_shutdown(2)
  end select
 case(2)
  select case (mcm_ndim)
   case(1)
    mcm_out_cols(1) = 4
	mcm_out_cols(2) = 11
    mcm_variablename(2)= 'Velocity Particle Material'
    mcm_variablename(3)= 'Vp x Velocity'
    mcm_variablename(4)= 'Vp x Acceleration'
    mcm_variablename(5)= 'Vp Density'
    mcm_variablename(6)= 'Stress Particle Material'
    mcm_variablename(7)= 'Sp x Velocity'
    mcm_variablename(8)= 'Sp x Acceleration'
    mcm_variablename(9)= 'Sp Density'
    mcm_variablename(10)='Pressure'
    mcm_variablename(11)='11 Stress'
    mcm_variablename(12)='Effective Plastic Strain'
    mcm_variablename(13)='Specific Internal Energy'
    mcm_variablename(14)='Speed of Sound'
    mcm_variablename(15)='Particle Critical Timestep'
    mcm_variablename(16)='Temperature'
   case(2)
    mcm_out_cols(1) = 6 
	mcm_out_cols(2) = 15
    mcm_variablename(3)= 'Velocity Particle Material'
    mcm_variablename(4)= 'Vp x Velocity'
    mcm_variablename(5)= 'Vp y Velocity'
    mcm_variablename(6)= 'Vp x Acceleration'
    mcm_variablename(7)= 'Vp y Acceleration'
    mcm_variablename(8)= 'Vp Density'
    mcm_variablename(9)= 'Stress Particle Material'
    mcm_variablename(10)='Sp x Velocity'
    mcm_variablename(11)='Sp y Velocity'
    mcm_variablename(12)='Sp x Acceleration'
    mcm_variablename(13)='Sp y Acceleration'
    mcm_variablename(14)='Sp Density'
    mcm_variablename(15)='Pressure'
    mcm_variablename(16)='11 Stress'
    mcm_variablename(17)='22 Stress'
    mcm_variablename(18)='12 Stress'
    mcm_variablename(19)='Effective Plastic Strain'
    mcm_variablename(20)='Specific Internal Energy'
    mcm_variablename(21)='Speed of Sound'
    mcm_variablename(22)='Particle Critical Timestep'
    mcm_variablename(23)='Temperature'
    case(3)
	 write(*,1000)
	 write(*,1010)
	 write(13,1000)
	 write(13,1010)
	 call mcm_shutdown(2)
  end select
end select
!
! calculate time for second state plot
!
 mcm_nextsttime = mcm_ptime + mcm_stpltime
!
return
!
1000 format(//5X, 'Error in subroutine initstateplot')
1010 format(/5X,'MCMGUI output format is only for 1D and 2D problems')
!
end subroutine mcm_initstateplot