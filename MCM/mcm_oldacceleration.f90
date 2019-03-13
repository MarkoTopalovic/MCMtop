      subroutine mcm_oldacceleration
!************************************************************************
!
!    Purpose: Calculates particle acceleration
!
!  Called by: MOMENTUM
!
!       Date: 15-01-99
!
!     Errors: 
!
!      Notes: 
!             
!             
!
!************************************************************************
!
use mcm_database
!
implicit none
!
integer:: nn,nn1,i,j,l,k,m,n
REAL(kind=real_acc)    :: Volj, dwdx(mcm_ndim),deltasig(mcm_ndim),havg,dwdr,rhoi
real(kind=real_acc), dimension(3,3) :: grad_v
real(kind=real_acc), dimension(3,3) :: sigmai, qi
!
real(kind=real_acc) :: massj, rhoj, hj, holdj
real(kind=real_acc), dimension(3) :: xj
real(kind=real_acc), dimension(3,3) :: sigmaj, qj
real(kind=real_acc) :: inv_w_deltap, epsilon, fij, ratio,npow
!
!	loop over all paticles or all velocity particles for the 
!	noncollocated discretization,   particle interactions
!
! Contact force for writing to disk
do i = 1, mcm_nummat
 mcm_mat(i)%xcf = 0.0
 mcm_mat(i)%ycf = 0.0
 mcm_mat(i)%zcf = 0.0
enddo
!
do i=mcm_svp,mcm_evp		!svp=start velocity point, evp=end velocity point
 par(i)%a = 0.0_d
 !________________________________________________________
 !
 ! Add acceleration due to contact force term if active
 !________________________________________________________
 !
 if(mcm_contacttype.GT.0) then
  do l = 1,mcm_ndim
   par(i)%a(l) = par(i)%a(l) + par(i)%repulsion(l)
  enddo
  mcm_mat(par(i)%mat)%xcf = mcm_mat(par(i)%mat)%xcf + par(i)%mass*par(i)%repulsion(1)
  mcm_mat(par(i)%mat)%ycf = mcm_mat(par(i)%mat)%ycf + par(i)%mass*par(i)%repulsion(2)
  mcm_mat(par(i)%mat)%zcf = mcm_mat(par(i)%mat)%zcf + par(i)%mass*par(i)%repulsion(3)
 endif
 !________________
 !
 ! End Contact
 !________________
 !
 rhoi = par(i)%rho
 sigmai = par(i)%sigma
 qi = par(i)%q
 !
 do k=1,par(i)%nnbr + par(i)%g_nnbr
  !
  call mcm_get_j_moment_info(i,k,xj,massj,rhoj,hj,holdj,sigmaj,qj)
  !
  Volj = massj/rhoj
  havg = 0.5_d*(par(i)%h+hj)
  deltasig = 0.0_d
  !
  call mcm_gradw(dwdx,dwdr,par(i)%x,xj,havg)
  !
  do n=1,mcm_ndim
   do m=1,mcm_ndim
    ! Sum form of basic SPH momentum equation
    deltasig(n) = deltasig(n)+ ( (sigmaj(m,n)-qj(m,n))/(rhoj**2)  +          &
	                             (sigmai(m,n)-qi(m,n))/(rhoi**2)) * dwdx(m)
   enddo
  enddo
  !
  do l=1,mcm_ndim
   ! sum form of basic SPH momentum equation
   par(i)%a(l) = par(i)%a(l) + massj * deltasig(l)
  enddo
  !
 enddo
enddo
! 
end subroutine mcm_oldacceleration