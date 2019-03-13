subroutine mcm_sets4(cm,prop)
!************************************************************************
!
!    Purpose: sets initial material properties
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
!implicit none                               
!
implicit double precision (a-h,o-z) 
integer:: i,j
real(kind=real_acc) :: cm(5,*),q1,q2,q3,prop(48), tempE, tempNi
!
      alpha = cm(1,3)
      theta = cm(1,4)
      gama = cm(1,5)
      beta = cm(1,6)
      r = cm(2,1)

      ltype = cm(2,6)
      xint = cm(2,5)
      ! E= 9KG/(3K+G)
      tempE=(9.*cm(1,1)*cm(1,2))/(3.*cm(1,1)+cm(1,2))
      !ni = (3K-2G)/2*(3K+G)
      tempNi=(3.*cm(1,1)-2.*cm(1,2))/(2.*(3*cm(1,1)+cm(1,2)))
      tempNi=(3.*cm(1,1)-tempE)/(6.*cm(1,1))
      q1=tempE*tempNi/((1.0+tempNi)*(1.0-2.0*tempNi))
    q2=tempE*0.5/(1.0+tempNi)
    q3=q1+2.0*q2
    cm(5,2) = q3

call initel (ltype,xint,alpha,gama,beta,theta,r,capint,fcut)
!
      cm(3,3) = capint
      cm(3,4) = fcut
!
end subroutine mcm_sets4

      subroutine initel(ltype,xint,alpha,gama,beta,theta,r,capint,fcut)
     implicit double precision (a-h,o-z)                                    

!c.... initialization of new cap model (newton algorithm)

      tol = 1.0e-3
      capint = 0.0
      delcap=1.e10
      delf=1.e10
      iter1=0
      iter2=0
      maxit=200
      fcut=0.0
    2 continue
      fe = alpha - gama*exp(-beta*capint) + theta*capint
      psik = xint - capint - r*fe
      err1=abs(psik)
      err2=abs(delcap)/(abs(capint)+1.e-15)
      if((err1.gt.tol).and.(err2.gt.tol)) then
!c     if( abs(psik) .gt. tol) then
      fep = beta*gama*exp(-beta*capint) + theta
      dpsik = -1.0 - r*fep
      delcap= psik/dpsik
      capint = capint - delcap
      if( iter1.gt.maxit ) then
      write(*,100)
      write(13,100)
 !     call adios(2)  ! used for printing in dyna3d when arror occurs
      endif
      iter1=iter1+1
      go to 2
      endif
    3 continue
      fe = alpha - gama*exp(-beta*fcut) + theta*fcut
      err1=abs(fe)
      err2=abs(delf)/(abs(fcut)+1.e-15)
      if( (err1.gt.tol).and.(err2.gt.tol) ) then
!c     if( abs(fe) .gt. tol) then
      fep = gama*beta*exp(-beta*fcut) + theta
      delf = fe/fep
      fcut = fcut - delf
      if( iter2.gt.maxit) then
      write(*,110)
      write(13,110)
!      call adios(2)
      endif
      iter2=iter2+1
      go to 3
      endif
      return
  100 format(///5x,'Iteration for initial kappa unconverged', /5x,'Selected value of x0 may be too large', /5x,'or other parameters inconsistent.',//)
  110 format(///5x,'Iteration for fcut unconverged', /5x,'failure surface may be too flat', /5x,'or other parameters inconsistent.',//)
  
      end subroutine initel