subroutine mcm_f3dm25 (cm,i)
!************************************************************************
!
!    Purpose: Granular material
!
!  Called by: constitutive
!
!       Date: 15-06-15
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
integer:: i
real(kind=real_acc) :: cm(5,*)

 call f3dm25a(cm,i)
return
!
end subroutine mcm_f3dm25


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!initel
       subroutine f3dm25a(cm,i) 
        
        use mcm_database

integer:: i
real(kind=real_acc) :: cm(5,*)
real(kind=real_acc) :: davg,ym,pr,gg,gdt,gd2,blk,bulk,capan
real(kind=real_acc) :: p, p_incr, othird
real(kind=real_acc) :: sig1,sig2,sig3,sig4,sig5,sig6,dd1,dd2,dd3,dd4,dd5,dd6, dd,davgm
parameter(tol=1.0e-7)
othird = 1.0_d/3.0_d
!ym=cm(1,1)		!Young's modulus
!pr=cm(1,2)		!Poison's ratio
!! gg=ym/(1.+pr)
!gg=cm(5,4)		!shear modulus
!ym = 22.65983894
!
!pr = 0.376661016
!
!gg = 8.2299995
!bulk = 30.620001
!blk = mcm_dt*3.*bulk

!
!	calculate volumetric strain increment and pressure 
!
! xm(i)=1./volo(i)
! vlrho(i)=crho*volo(i)
!blk=mcm_dt*ym/((1.0_d-2.0_d*pr))		!Bulk modulus *3 *dt
!davg=othird*par(i)%tracerod

  
!c     common/bk02/iburn,dt1,dt2,isdo                                    
!c     common/aux2/dd1(128),dd2(128),dd3(128),dd4(128),dd5(128),         
!c    1 dd6(128),wzzdt(128),wyydt(128),wxxdt(128)                        
!c     common/aux14/                                                     
!c    1 sig1(128),sig2(128),sig3(128),sig4(128),                         
!c    2 sig5(128),sig6(128),epx1(128), capa(128),                        
!c    3 epa1(128),epa2(128),epa4(128),epa5(128),epa6(128)                
!c     common/aux18/dd(128),def(128)                                     
!c     common/aux33/ix1(128),ix2(128),ix3(128),ix4(128),ix5(128),        
!c    1             ix6(128),ix7(128),ix8(128),mxt(128),nmel             
!c     common/aux35/rhoa(128),cb(128),davg(128),p(128)                   
!c     common/aux36/lft,llt                                              
!      common/dbk02/iburn,dt1,dt2,isdo                                   
!      common/daux2/dd1(128),dd2(128),dd3(128),dd4(128),dd5(128),        
!     1 dd6(128),wzzdt(128),wyydt(128),wxxdt(128)                        
!      common/daux14/                                                    
!     1 sig1(128),sig2(128),sig3(128),sig4(128),                         
!     2 sig5(128),sig6(128),epx1(128), capa(128),                        
!     3 epa1(128),epa2(128),epa4(128),epa5(128),epa6(128)                
!      common/daux18/dd(128),def(128)                                    
!      common/daux33/ix1(128),ix2(128),ix3(128),ix4(128),ix5(128),       
!     1             ix6(128),ix7(128),ix8(128),mxt(128),nmel             
!      common/daux35/rhoa(128),cb(128),davg(128),p(128)                  
!      common/daux36/lft,llt                                             
!      common/dauxplt/ el,xl,evpi,sj1,sj2,tmises,fj1,xmtype,xit,nocon    
!      common/double/iprec,ncpw,unit
!      dimension gtr11(128),gtr22(128),gtr33(128),gtr12(128),
!     1          gtr23(128),gtr31(128)
!      dimension s1(128),s2(128),s3(128),s4(128),s5(128),s6(128),
!     1      x1tr(128),
!     1      x2dtr(128),x1(128),x2d(128),tmises(128),tcrit(128)
!      dimension f1tr(128),fe(128),fep(128),fc2tr(128),f2tr(128)
!      dimension capan(128),capal(128),xk(128),
!     1          evpn(128),devp(128),
!     2          xkn(128),elcapn(128),
!     3          devpb(128),x1crt(128),evp(128),iupd(128),imode(128),
!     4          iter(128)
!       dimension dhdk(128),dfdk(128),elcap(128),omega(128),dcd(128),f(128),fac1(128),fac2(128),sdtr11(128),sdtr22(128),sdtr33(128),sdev(6,128),fbar(128)
!      dimension  cm(1)
!      data third /0.333333333333333/
!c
!c.... mtype :  0 = tension cutoff, 1 = elastic, 2 = failure,
!c....          3 = cap mode,       4= compression corner mode
!c....          5 = tension corner mode
!c....          6 = 'gap region' (cap will catch up)
!c.... ltype : 1 = soil, 2 = rock
!c

      third = 1.0/3.0
      
      unit=1.0
      geop=0.
!!      mx=48*(mxt(lft)-1)
      bulk=cm(1,1)
      shear=cm(1,2)
      alpha = cm(1,3)
      theta = cm(1,4)
      gama = cm(1,5)
      beta = cm(1,6)
      r=cm(2,1)
      dDyna = cm(2,2) ! d nije bilo zamenjeno sa dDyna u celom fajlu !
      w=cm(2,3)
      cbar=cm(3,6)
      falfac=cm(4,1)
      alpha=alpha-falfac
      z=cm(2,4)
      nplot=cm(2,5)
      ivec=1
      if (cm(3,5).ne.0.) ivec=0
                                                    !!      if (falfac .ne. 0.) then
                                                    !
                                                    !!      trs=(par(i)%sigma(1,1)+par(i)%sigma(2,2)+par(i)%sigma(3,3))/3.
                                                    !!      sdev(1,i)=par(i)%sigma(1,1)-trs
                                                    !!      sdev(2,i)=par(i)%sigma(2,2)-trs
                                                    !!      sdev(3,i)=par(i)%sigma(3,3)-trs
                                                    !!      sdev(4,i)=par(i)%sigma(1,2)
                                                    !!      sdev(5,i)=par(i)%sigma(1,3)
                                                    !!      sdev(6,i)=par(i)%sigma(2,3)
                                                    !
                                                    !!      do 510 i=lft,llt
                                                    !!      sig1(i)=sig1(i)-epa1(i)
                                                    !!      sig2(i)=sig2(i)-epa2(i)
                                                    !!      sig3(i)=sig3(i)+epa1(i)+epa2(i)
                                                    !!      sig4(i)=sig4(i)-epa4(i)
                                                    !!      sig5(i)=sig5(i)-epa5(i)
                                                    !!      sig6(i)=sig6(i)-epa6(i)
                                                    !!  510 continue
                                                    !!      endif
      ltype=cm(2,6)
      tcut=cm(3,1)
      fcut=cm(3,4)
      t=max(fcut,tcut+3.*geop)
      fet=alpha - gama*exp(-beta*t) + theta*t
      fept = gama*beta*exp(-beta*t) + theta
      
     ! par(i)%capa = cm(3,3)
      par(i)%capa = max(par(i)%capa,cm(3,3))

      sig1 = -par(i)%sigma(1,1)
      sig2 = -par(i)%sigma(2,2)
      sig3 = -par(i)%sigma(3,3)
      sig4 = -par(i)%sigma(1,2)
      sig5 = -par(i)%sigma(1,3)
      sig6 = -par(i)%sigma(2,3)
      dd1 = -par(i)%rod(1,1)
      dd2 = -par(i)%rod(2,2)
      dd3 = -par(i)%rod(3,3)
      dd4 = -par(i)%rod(1,2)/2.0
      dd5 = -par(i)%rod(1,3)/2.0
      dd6 = -par(i)%rod(2,3)/2.0
!      iter=0
!
      dd  = (dd1+dd2+dd3)
      davgm=(1./3.)*dd
      dt1 =mcm_dt
      p=3.0*bulk*dt1*davgm
      p_incr=blk*davg				!Presure increment=volume strain *bulk modulus
  
  gdt= mcm_dt*shear
  gd2=2.0*gdt
  ! 2.0*shear*dt1
   
      gtr11 = sig1 + p +gd2*(dd1 - davgm)
      gtr22 = sig2 + p +gd2*(dd2 - davgm)
      gtr33 = sig3 + p +gd2*(dd3 - davgm)
      gtr12 = sig4        +gd2* dd4
      gtr23 = sig5        +gd2* dd5
      gtr31 = sig6        +gd2* dd6
      x1tr = gtr11 + gtr22 + gtr33
      s1 = gtr11 - (1./3.)*x1tr
      s2 = gtr22 - (1./3.)*x1tr
      s3 = gtr33 - (1./3.)*x1tr
      sdtr11=s1
      sdtr22=s2
      sdtr33=s3
      x2dtr=0.5*(s1**2 + s2**2 + s3**2 + 2.*(gtr12**2 + gtr23**2 + gtr31**2) )
      x2dtr=sqrt(x2dtr)
!
!        goto 901
      fe = alpha - gama*exp(-beta*capa) + theta*capa
      xkn = capa + r*fe
      elcapn = max(par(i)%capa,0.0*unit)
      tmises=abs(xkn - elcapn)/r

!c

      fe = alpha - gama*exp(-beta*x1tr) + theta*x1tr
      f1tr = x2dtr - min(fe,tmises)
      tcrit = t - 9.*bulk/shear*(x2dtr - fet)*fept
      fe = alpha - gama*exp(-beta*elcapn) + theta*elcapn
      fep = beta*gama*exp(-beta*elcapn) + theta
      fc2tr=((xkn-elcapn)**2-(x1tr-elcapn)**2)/r**2
      f2tr=x2dtr**2 - fc2tr
      x1crt = elcapn - (x2dtr - fe) *9.*bulk*fep/shear
      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      capal = min(0.0*unit,par(i)%capa - tol*fe)
      fe = alpha - gama*exp(-beta*capal) + theta*capal
      xk = capal + r*fe
      evpn =  w*(1. - exp(dDyna*(z-xkn)))
      devpb = w*(1. - exp(dDyna*(z-xk))) - evpn
      capan = par(i)%capa
      iupd=1


      if( (x1tr.le. t) .and. (x2dtr.le.fet) .and. (ltype.eq.1).and. (capan.gt.0.0) .and. (devpb.ge.0.0) ) then
           imode=2
      elseif( (x1tr.le. t) .and. (x2dtr.le.fet) .and.  (ltype.eq.1)  .and. (capan.gt.0.0) .and. (devpb.lt.0.0) ) then
           imode=7
      elseif( (x1tr.le. t) .and. (x2dtr.le.fet) .and. ((ltype.eq.2).or.(capan.le.0.0)) ) then
           imode=3
      elseif( (x1tr.le. t) .and. (x2dtr.gt.fet) .and.  (x1tr.le.tcrit) .and. (ltype.eq.1) .and. (capan.gt.0.0) .and. (devpb.ge.0.0) ) then
           imode=4
      elseif( (x1tr.le. t) .and. (x2dtr.gt.fet) .and. (x1tr.le.tcrit) .and. (ltype.eq.1)  .and. (capan.gt.0.0) .and. (devpb.lt.0.0) ) then
           imode=6
      elseif( (x1tr.le. t) .and. (x2dtr.gt.fet) .and.  (x1tr.le.tcrit) .and. ((ltype.eq.2) .or.  (capan.le.0.0) ) ) then
           imode=5
      elseif( (x1tr.le. t) .and. (x2dtr.gt.fet) .and.  (x1tr.gt.tcrit) .and. (capan.ge.0.0) ) then
           imode=10
      elseif( (x1tr.le. t) .and. (x2dtr.gt.fet) .and.  (x1tr.gt.tcrit) .and. (capan.lt.0.0) ) then
           imode=8
      elseif( (x1tr.gt. t) .and. (x1tr.le.elcapn)  .and. (f1tr.le.tol) ) then
           imode=1
      elseif( (x1tr.gt. t) .and. (x1tr.le.elcapn)  .and. (f1tr.gt.tol) .and. (x1tr.ge. x1crt) ) then
           imode=9
      elseif( (x1tr.gt. t) .and. (x1tr.le.elcapn)  .and. (f1tr.gt.tol) .and. (x1tr.lt. x1crt)  .and. (capan.ge.0.0) ) then
           imode=10
      elseif( (x1tr.gt. t) .and. (x1tr.le.elcapn)  .and. (f1tr.gt.tol) .and. (x1tr.lt. x1crt) .and. (capan.lt.0.0) ) then
           imode=8
      elseif( (x1tr.gt. t) .and. (x1tr.gt.elcapn)  .and. (f2tr.lt.tol)) then
           imode=1
      elseif( (x1tr.gt. t) .and. (x1tr.gt.elcapn)  .and. (f2tr.ge.tol)) then
           imode=11
      else

      endif

!c900
          
      par(i)%ptype = imode
      
      if(imode.eq.1) then
      x1 = x1tr
      x2d = x2dtr
      par(i)%capa = capan
      p = third*x1 + geop
      sig1 = s1 + p
      sig2 = s2 + p
      sig3 = s3 + p
      sig4 = gtr12
      sig5 = gtr23
      sig6 = gtr31
      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      xk = par(i)%capa + r*fe
      evp = w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif

      if(imode.eq.2) then
      x1=t
      x2d=x2dtr
      par(i)%capa=0.0
      p=third*x1 + geop
      sig1 = s1 + p
      sig2 = s2 + p
      sig3 = s3 + p
      sig4 = gtr12
      sig5 = gtr23
      sig6 = gtr31
      xk = alpha - gama
      evp=w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif

      if(imode.eq.3) then
      x1=t
      x2d=x2dtr
      par(i)%capa=capan
      p=third*x1 + geop
      sig1 = s1 + p
      sig2 = s2 + p
      sig3 = s3 + p
      sig4 = gtr12
      sig5 = gtr23
      sig6 = gtr31
      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      xk = par(i)%capa + r*fe
      evp=w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif
!c
!c.... new mode
!c

      if(imode.eq.4) then
      x1=t
      x2d=fet
      par(i)%capa=0.0
      p=third*x1 + geop
      ratio = x2d/x2dtr
      sig1 = s1*ratio + p
      sig2 = s2*ratio + p
      sig3 = s3*ratio + p
      sig4 = gtr12*ratio
      sig5 = gtr23*ratio
      sig6 = gtr31*ratio
      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      xk = par(i)%capa + r*fe
      evp=w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif


      if(imode.eq.5) then
      x1=t
      x2d=fet
      par(i)%capa = capan
      p = third*x1 + geop
      ratio=x2d/x2dtr
      sig1 = s1*ratio + p
      sig2 = s2*ratio + p
      sig3 = s3*ratio + p
      sig4 = gtr12*ratio
      sig5 = gtr23*ratio
      sig6 = gtr31*ratio
      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      xk = par(i)%capa + r*fe
      evp = w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif

 
      if (imode.eq.6) then
      x1=t
      x2d=fet
      x1old = sig1 + sig2 + sig3
      devp =  dd - (x1 - x1old)/(3.*bulk)
      xhat = min(devpb,devp)
      par(i)%capa = par(i)%capa + devp*(capal - capan)/xhat
      par(i)%capa = max(par(i)%capa,0.0*unit)
      p = third*x1 + geop
      ratio=x2d/x2dtr
      sig1 = s1*ratio + p
      sig2 = s2*ratio + p
      sig3 = s3*ratio + p
      sig4 = gtr12*ratio
      sig5 = gtr23*ratio
      sig6 = gtr31*ratio
      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      xk = par(i)%capa + r*fe
      evp = w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif


      if( imode.eq.7) then
      x1=t
      x2d=x2dtr
      x1old = sig1 + sig2 + sig3
      devp =  dd - (x1 - x1old)/(3.*bulk)
      xhat = min(devpb,devp)
      par(i)%capa = par(i)%capa + devp*(capal - capan)/xhat
      par(i)%capa = max(par(i)%capa,0.0*unit)
      p = third*x1 + geop
      sig1 = s1 + p
      sig2 = s2 + p
      sig3 = s3 + p
      sig4 = gtr12
      sig5 = gtr23
      sig6 = gtr31
      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      xk = par(i)%capa + r*fe
      evp = w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif

!c
!c.... new mode
!c

      if(imode.eq.8) then
      x1 = x1tr
      fe = alpha - gama*exp(-beta*x1) + theta*x1
      x2d=min(fe,tmises)
      par(i)%capa = capan
      ratio = x2d/x2dtr
      p = third*x1 + geop
      sig1 = s1*ratio + p
      sig2 = s2*ratio + p
      sig3 = s3*ratio + p
      sig4 = gtr12*ratio
      sig5 = gtr23*ratio
      sig6 = gtr31*ratio
      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      xk = par(i)%capa + r*fe
      evp = w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif


      if(imode.eq.9) then
      x1 = elcapn
      x2d = alpha - gama*exp(-beta*x1) + theta*x1
      par(i)%capa = capan
      ratio = x2d/x2dtr
      p = third*x1 + geop
      sig1 = s1*ratio + p
      sig2 = s2*ratio + p
      sig3 = s3*ratio + p
      sig4 = gtr12*ratio
      sig5 = gtr23*ratio
      sig6 = gtr31*ratio
      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      xk = par(i)%capa + r*fe
      evp = w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif


      if(imode.eq.10 .and. ivec.eq.0) then
      sig1=gtr11
      sig2=gtr22
      sig3=gtr33
      sig4=gtr12
      sig5=gtr23
      sig6=gtr31

      call feit (capan,x1tr,x2dtr,t,alpha,gama,beta,theta,w,dDyna,r,bulk,shear,x1,par(i)%capa,sig1,sig2, sig3,sig4,sig5,sig6,iter,nocon)

      if(x1 .gt. elcapn) x1=elcapn
      if( (ltype.eq.2) .or. (capan.eq.0.0) ) par(i)%capa=capan
      if( x1 .gt. par(i)%capa) par(i)%capa = x1
      if((ltype.eq.1).and.(capan.ge.0.0))  par(i)%capa=max(par(i)%capa,0.0*unit)
      fe = alpha - gama*exp(-beta*x1) + theta*x1
      x2d=min(fe,tmises)
      ratio = x2d/x2dtr
      p = third*x1 + geop
      sig1 = s1*ratio + p
      sig2 = s2*ratio + p
      sig3 = s3*ratio + p
      sig4 = gtr12*ratio
      sig5 = gtr23*ratio
      sig6 = gtr31*ratio
      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      xk = par(i)%capa + r*fe
      evp = w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif

!c

      if(imode.eq.10 .and. ivec.eq.1) then
      sig1=gtr11
      sig2=gtr22
      sig3=gtr33
      sig4=gtr12
      sig5=gtr23
      sig6=gtr31
      iter=4
      par(i)%capa=capan
      x1=x1tr
      x2d=x2dtr
      s1=sig1-x1tr/3.
      s2=sig2-x1tr/3.
      s3=sig3-x1tr/3.
      s4=sig4
      s5=sig5
      s6=sig6
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
      dhdk=w*dDyna*exp(-dDyna*xk)*(1.+r*(theta+gama*beta     *exp(-beta*par(i)%capa)))
      omega=theta+gama*beta*exp(-beta*x1)
      dcd=shear+9.*bulk*omega**2
      f=x2d-(alpha-gama*exp(-beta*x1)+theta*x1)
      delcap=-3.*omega*f/(dcd*dhdk)
      delam=-dhdk*delcap/(3.*omega)
      par(i)%capa=par(i)%capa+delcap
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
      fac1=shear/x2d
      fac2=3.*bulk*omega
      sig1=sig1-delam*(fac1*s1-fac2)
      sig2=sig2-delam*(fac1*s2-fac2)
      sig3=sig3-delam*(fac1*s3-fac2)
      sig4=sig4-delam*fac1*s4
      sig5=sig5-delam*fac1*s5
      sig6=sig6-delam*fac1*s6
      x1=sig1+sig2+sig3
      s1=sig1-x1/3.
      s2=sig2-x1/3.
      s3=sig3-x1/3.
      s4=sig4
      s5=sig5
      s6=sig6
      x2d=sqrt(.5*(s1**2+s2**2+s3**2)     +s4**2+s5**2+s6**2)
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
      dhdk=w*dDyna*exp(-dDyna*xk)*(1.+r*(theta+gama*beta *exp(-beta*par(i)%capa)))
      omega=theta+gama*beta*exp(-beta*x1)
      dcd=shear+9.*bulk*omega**2
      f=x2d-(alpha-gama*exp(-beta*x1)+theta*x1)
      delcap=-3.*omega*f/(dcd*dhdk)
      delam=-dhdk*delcap/(3.*omega)
      par(i)%capa=par(i)%capa+delcap
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
      fac1=shear/x2d
      fac2=3.*bulk*omega
      sig1=sig1-delam*(fac1*s1-fac2)
      sig2=sig2-delam*(fac1*s2-fac2)
      sig3=sig3-delam*(fac1*s3-fac2)
      sig4=sig4-delam*fac1*s4
      sig5=sig5-delam*fac1*s5
      sig6=sig6-delam*fac1*s6
      x1=sig1+sig2+sig3
      s1=sig1-x1/3.
      s2=sig2-x1/3.
      s3=sig3-x1/3.
      s4=sig4
      s5=sig5
      s6=sig6
      x2d=sqrt(.5*(s1**2+s2**2+s3**2)  +s4**2+s5**2+s6**2)
      endif


      if(imode.eq.10 .and. ivec.eq.1) then
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
      dhdk=w*dDyna*exp(-dDyna*xk)*(1.+r*(theta+gama *beta*exp(-beta*par(i)%capa)))
      omega=theta+gama*beta*exp(-beta*x1)
      dcd=shear+9.*bulk*omega**2
      f=x2d-(alpha-gama*exp(-beta*x1)+theta*x1)
      delcap=-3.*omega*f/(dcd*dhdk)
      delam=-dhdk*delcap/(3.*omega)
      par(i)%capa=par(i)%capa+delcap
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
      fac1=shear/x2d
      fac2=3.*bulk*omega
      sig1=sig1-delam*(fac1*s1-fac2)
      sig2=sig2-delam*(fac1*s2-fac2)
      sig3=sig3-delam*(fac1*s3-fac2)
      sig4=sig4-delam*fac1*s4
      sig5=sig5-delam*fac1*s5
      sig6=sig6-delam*fac1*s6
      x1=sig1+sig2+sig3
      s1=sig1-x1/3.
      s2=sig2-x1/3.
      s3=sig3-x1/3.
      s4=sig4
      s5=sig5
      s6=sig6
      x2d=sqrt(.5*(s1**2+s2**2+s3**2)       +s4**2+s5**2+s6**2)
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
      dhdk=w*dDyna*exp(-dDyna*xk)*(1.+r*(theta+gama  *beta*exp(-beta*par(i)%capa)))
      omega=theta+gama*beta*exp(-beta*x1)
      dcd=shear+9.*bulk*omega**2
      f=x2d-(alpha-gama*exp(-beta*x1)+theta*x1)
      delcap=-3.*omega*f/(dcd*dhdk)
      delam=-dhdk*delcap/(3.*omega)
      par(i)%capa=par(i)%capa+delcap
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
      fac1=shear/x2d
      fac2=3.*bulk*omega
      sig1=sig1-delam*(fac1*s1-fac2)
      sig2=sig2-delam*(fac1*s2-fac2)
      sig3=sig3-delam*(fac1*s3-fac2)
      sig4=sig4-delam*fac1*s4
      sig5=sig5-delam*fac1*s5
      sig6=sig6-delam*fac1*s6
      x1=sig1+sig2+sig3
      s1=sig1-x1/3.
      s2=sig2-x1/3.
      s3=sig3-x1/3.
      s4=sig4
      s5=sig5
      s6=sig6
      x2d=sqrt(.5*(s1**2+s2**2+s3**2)   +s4**2+s5**2+s6**2)
      endif


      if(imode.eq.10 .and. ivec.eq.1) then
      if(x1 .gt. elcapn) x1=elcapn
      if( (ltype.eq.2) .or. (capan.eq.0.0) ) par(i)%capa=capan
      if( x1 .gt. par(i)%capa) par(i)%capa = x1
      if((ltype.eq.1).and.(capan.ge.0.0))      par(i)%capa=max(par(i)%capa,0.0*unit)
      endif

!c

      if(imode.eq.10 .and. ivec.eq.1) then
      fe = alpha - gama*exp(-beta*x1) + theta*x1
      x2d=min(fe,tmises)
      ratio = x2d/x2dtr
      p = third*x1 + geop
      sig1 = sdtr11*ratio + p
      sig2 = sdtr22*ratio + p
      sig3 = sdtr33*ratio + p
      sig4 = gtr12*ratio
      sig5 = gtr23*ratio
      sig6 = gtr31*ratio
      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      xk = par(i)%capa + r*fe
      evp = w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif

      if(imode.eq.11 .and. ivec.eq.0) then
      sig1=gtr11
      sig2=gtr22
      sig3=gtr33
      sig4=gtr12
      sig5=gtr23
      sig6=gtr31

      call fcit (capan,x1tr,x2dtr,t,alpha,gama,beta, theta,w,dDyna,r,bulk,shear,x1,par(i)%capa,sig1,sig2, sig3,sig4,sig5,sig6,iter,nocon)

      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      xk = par(i)%capa + r*fe
      evp = w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif


      if(imode.eq.11 .and. ivec.eq.1) then
      sig1=gtr11
      sig2=gtr22
      sig3=gtr33
      sig4=gtr12
      sig5=gtr23
      sig6=gtr31
      iter=4
      par(i)%capa=capan
      x1=x1tr
      x2d=x2dtr
      s1=sig1-x1tr/3.
      s2=sig2-x1tr/3.
      s3=sig3-x1tr/3.
      s4=sig4
      s5=sig5
      s6=sig6
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
      dhdk=w*dDyna*exp(-dDyna*xk)*(1.+r*(theta+gama    *beta*exp(-beta*par(i)%capa)))
      if(par(i)%capa.ge.0.0) then
      dfdk=-2.*(abs(x1-par(i)%capa)/r**2    + (alpha-gama*exp(-beta*par(i)%capa)   +theta*par(i)%capa)*(theta+gama*beta*exp(-beta*par(i)%capa)))
      else
      dfdk=-2.*xk*(1.+r*(theta+gama*beta  *exp(-beta*par(i)%capa)))/r**2
      endif
      elcap=max(par(i)%capa,0.0*unit)
      omega=min(-2.*abs(x1-elcap)/r**2,-1.e-15*unit)
      dcd=4.*shear*x2d**2 + 9.*bulk*(2.*(x1  -elcap)/r**2)**2
      f=x2d**2 + (x1-elcap)**2/r**2   - (xk-elcap)**2/r**2
      delcap=-f/(dfdk+(dcd*dhdk)/(3.*omega))
      delam=-dhdk*delcap/(3.*omega)
      par(i)%capa=par(i)%capa+delcap
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
!c     elcap=max(par(i)%capa,0.0*unit)
      fac1=2.*shear
      fac2=6.*bulk*abs(x1-elcap)/r**2
      sig1=sig1-delam*(fac1*s1+fac2)
      sig2=sig2-delam*(fac1*s2+fac2)
      sig3=sig3-delam*(fac1*s3+fac2)
      sig4=sig4-delam*fac1*s4
      sig5=sig5-delam*fac1*s5
      sig6=sig6-delam*fac1*s6
      x1=sig1+sig2+sig3
      s1=sig1-x1/3.
      s2=sig2-x1/3.
      s3=sig3-x1/3.
      s4=sig4
      s5=sig5
      s6=sig6
      x2d=sqrt(.5*(s1**2+s2**2+s3**2) +s4**2+s5**2+s6**2)
      endif


      if(imode.eq.11 .and. ivec.eq.1) then
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
      dhdk=w*dDyna*exp(-dDyna*xk)*(1.+r*(theta+gama  *beta*exp(-beta*par(i)%capa)))
      if(par(i)%capa.ge.0.0) then
      dfdk=-2.*(abs(x1-par(i)%capa)/r**2   + (alpha-gama*exp(-beta*par(i)%capa)   +theta*par(i)%capa)*(theta+gama*beta*exp(-beta*par(i)%capa)))
      else
      dfdk=-2.*xk*(1.+r*(theta+gama*beta   *exp(-beta*par(i)%capa)))/r**2
      endif
      elcap=max(par(i)%capa,0.0*unit)
      omega=min(-2.*abs(x1-elcap)/r**2,-1.e-15*unit)
      dcd=4.*shear*x2d**2 + 9.*bulk*(2.*(x1  -elcap)/r**2)**2
      f=x2d**2 + (x1-elcap)**2/r**2     - (xk-elcap)**2/r**2
      delcap=-f/(dfdk+(dcd*dhdk)/(3.*omega))
      delam=-dhdk*delcap/(3.*omega)
      par(i)%capa=par(i)%capa+delcap
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
!c     elcap=max(par(i)%capa,0.0*unit)
      fac1=2.*shear
      fac2=6.*bulk*abs(x1-elcap)/r**2
      sig1=sig1-delam*(fac1*s1+fac2)
      sig2=sig2-delam*(fac1*s2+fac2)
      sig3=sig3-delam*(fac1*s3+fac2)
      sig4=sig4-delam*fac1*s4
      sig5=sig5-delam*fac1*s5
      sig6=sig6-delam*fac1*s6
      x1=sig1+sig2+sig3
      s1=sig1-x1/3.
      s2=sig2-x1/3.
      s3=sig3-x1/3.
      s4=sig4
      s5=sig5
      s6=sig6
      x2d=sqrt(.5*(s1**2+s2**2+s3**2) +s4**2+s5**2+s6**2)
      endif


      if(imode.eq.11 .and. ivec.eq.1) then
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
      dhdk=w*dDyna*exp(-dDyna*xk)*(1.+r*(theta+gama    *beta*exp(-beta*par(i)%capa)))
      if(par(i)%capa.ge.0.0) then
      dfdk=-2.*(abs(x1-par(i)%capa)/r**2   + (alpha-gama*exp(-beta*par(i)%capa)    +theta*par(i)%capa)*(theta+gama*beta*exp(-beta*par(i)%capa)))
      else
      dfdk=-2.*xk*(1.+r*(theta+gama*beta    *exp(-beta*par(i)%capa)))/r**2
      endif
      elcap=max(par(i)%capa,0.0*unit)
      omega=min(-2.*abs(x1-elcap)/r**2,-1.e-15*unit)
      dcd=4.*shear*x2d**2 + 9.*bulk*(2.*(x1   -elcap)/r**2)**2
      f=x2d**2 + (x1-elcap)**2/r**2     - (xk-elcap)**2/r**2
      delcap=-f/(dfdk+(dcd*dhdk)/(3.*omega))
      delam=-dhdk*delcap/(3.*omega)
      par(i)%capa=par(i)%capa+delcap
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
!c     elcap=max(par(i)%capa,0.0*unit)
      fac1=2.*shear
      fac2=6.*bulk*abs(x1-elcap)/r**2
      sig1=sig1-delam*(fac1*s1+fac2)
      sig2=sig2-delam*(fac1*s2+fac2)
      sig3=sig3-delam*(fac1*s3+fac2)
      sig4=sig4-delam*fac1*s4
      sig5=sig5-delam*fac1*s5
      sig6=sig6-delam*fac1*s6
      x1=sig1+sig2+sig3
      s1=sig1-x1/3.
      s2=sig2-x1/3.
      s3=sig3-x1/3.
      s4=sig4
      s5=sig5
      s6=sig6
      x2d=sqrt(.5*(s1**2+s2**2+s3**2) +s4**2+s5**2+s6**2)
      endif

!c

      if(imode.eq.11 .and. ivec.eq.1) then
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
      dhdk=w*dDyna*exp(-dDyna*xk)*(1.+r*(theta+gama   *beta*exp(-beta*par(i)%capa)))
      if(par(i)%capa.ge.0.0) then
      dfdk=-2.*(abs(x1-par(i)%capa)/r**2   + (alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)*(theta+gama*beta*exp(-beta*par(i)%capa)))
      else
      dfdk=-2.*xk*(1.+r*(theta+gama*beta    *exp(-beta*par(i)%capa)))/r**2
      endif
      elcap=max(par(i)%capa,0.0*unit)
      omega=min(-2.*abs(x1-elcap)/r**2,-1.e-15*unit)
      dcd=4.*shear*x2d**2 + 9.*bulk*(2.*(x1 -elcap)/r**2)**2
      f=x2d**2 + (x1-elcap)**2/r**2    - (xk-elcap)**2/r**2
      delcap=-f/(dfdk+(dcd*dhdk)/(3.*omega))
      delam=-dhdk*delcap/(3.*omega)
      par(i)%capa=par(i)%capa+delcap
      xk=par(i)%capa+r*(alpha-gama*exp(-beta*par(i)%capa)+theta*par(i)%capa)
!c     elcap=max(par(i)%capa,0.0*unit)
      fac1=2.*shear
      fac2=6.*bulk*abs(x1-elcap)/r**2
      sig1=sig1-delam*(fac1*s1+fac2)
      sig2=sig2-delam*(fac1*s2+fac2)
      sig3=sig3-delam*(fac1*s3+fac2)
      sig4=sig4-delam*fac1*s4
      sig5=sig5-delam*fac1*s5
      sig6=sig6-delam*fac1*s6
      x1=sig1+sig2+sig3
      s1=sig1-x1/3.
      s2=sig2-x1/3.
      s3=sig3-x1/3.
      s4=sig4
      s5=sig5
      s6=sig6
      x2d=sqrt(.5*(s1**2+s2**2+s3**2)+s4**2+s5**2+s6**2)
      if(x1.lt.elcap) then
      x1=elcap+abs(x1-elcap)
      sig1=s1+x1/3.
      sig2=s2+x1/3.
      sig3=s3+x1/3.
      endif
      fe = alpha - gama*exp(-beta*par(i)%capa) + theta*par(i)%capa
      xk = par(i)%capa + r*fe
      evp = w*(1. - exp(dDyna*(z-xk)))
      iupd=0
      endif

  900 continue

      par(i)%sigma(1,1) = -sig1
      par(i)%sigma(2,2) = -sig2
      par(i)%sigma(3,3) = -sig3
      par(i)%sigma(1,2) = -sig4
      par(i)%sigma(1,3) = -sig5
      par(i)%sigma(2,3) = -sig6
      par(i)%rod(1,1) = -dd1
      par(i)%rod(2,2) = -dd2
      par(i)%rod(3,3) = -dd3
      par(i)%rod(1,2) = -dd4*2.0
      par(i)%rod(1,3) = -dd5*2.0
      par(i)%rod(2,3) = -dd6*2.0

! 901 continue
!        par(i)%sigma(1,1) = -gtr11
!        par(i)%sigma(2,2) = -gtr22
!        par(i)%sigma(3,3) = -gtr33
!        par(i)%sigma(1,2) = -gtr12
!        par(i)%sigma(1,3) = -gtr31
!        par(i)%sigma(2,3) = -gtr23
!        par(i)%rod(1,1)=-dd1
!        par(i)%rod(2,2)=-dd2
!        par(i)%rod(3,3)=-dd3
!        par(i)%rod(1,2)=-2.0*dd4
!        par(i)%rod(1,3)=-2.0*dd5
!        par(i)%rod(2,3)=-2.0*dd6

        par(i)%p = -othird * (par(i)%sigma(1,1) + par(i)%sigma(2,2) + par(i)%sigma(3,3))


!      idone=idone+iupd


!      if (nplot.eq.0) epx1(i)=0.0
!      if (nplot.eq.1) epx1(i)=par(i)%capa
!      if (nplot.eq.2) epx1(i)=xk(i)
!      if (nplot.eq.3) epx1(i)=evp(i)
!      if (nplot.eq.4) epx1(i)=x1(i)
!      if (nplot.eq.5) epx1(i)=x2d(i)
!      if (nplot.eq.6) epx1(i)=tmises(i)
!      if (nplot.eq.7) epx1(i)=0.0
!      if (nplot.eq.8) epx1(i)=imode(i)
!      if (nplot.eq.9) epx1(i)=iter(i)

!      if(idone.ne.0) stop 'f3dm252'

                                                                                    !      if (falfac .ne. 0.) then

                                                                                    !      epa3=-epa1-epa2
                                                                                    !      fac7= 0.5*(sig1*epa1+sig2*epa2+sig3*epa3+2.*(sig4*epa4+sig5*epa5+sig6*epa6))
                                                                                    !      fac8 = alpha - gama*exp(-beta*x1) + theta*x1
                                                                                    !      fbar=max(0.0*unit,1.-fac7/(1.*falfac*fac8))
                                                                                    !
                                                                                    !      fac7=1./(2.*shear)
                                                                                    !
                                                                                    !      dtrs=(dd1+dd2+dd3)/3.
                                                                                    !      ddev1=dd1-dtrs
                                                                                    !      ddev2=dd2-dtrs
                                                                                    !      strs=(sig1+sig2+sig3)/3.
                                                                                    !      sx1=sig1-strs
                                                                                    !      sx2=sig2-strs
                                                                                    !      epa1=epa1+cbar*fbar*(ddev1*dt1-fac7*(sx1-(sdev(1,i)-epa1)))
                                                                                    !      epa2=epa2+cbar*fbar*(ddev2*dt1-fac7*(sx2-(sdev(2,i)-epa2)))
                                                                                    !      epa4=epa4+cbar*fbar*(dd4/2.*dt1-fac7*(sig4-(sdev(4,i)-epa4)))
                                                                                    !      epa5=epa5+cbar*fbar*(dd5/2.*dt1-fac7*(sig5-(sdev(5,i)-epa5)))
                                                                                    !      epa6=epa6+cbar*fbar*(dd6/2.*dt1-fac7*(sig6-(sdev(6,i)-epa6)))
                                                                                    !
                                                                                    !!c
                                                                                    !!c.... radial return
                                                                                    !!c
                                                                                    !
                                                                                    !      epa3=-epa1-epa2
                                                                                    !      anorm=sqrt(0.5*(epa1**2+epa2**2+epa3**2+2.*(epa4**2+epa5**2+epa6**2)))
                                                                                    !      if (anorm .gt. falfac) then
                                                                                    !      fac8=falfac/anorm
                                                                                    !      epa1=epa1*fac8
                                                                                    !      epa2=epa2*fac8
                                                                                    !      epa4=epa4*fac8
                                                                                    !      epa5=epa5*fac8
                                                                                    !      epa6=epa6*fac8
                                                                                    !      endif
                                                                                    !
                                                                                    !
                                                                                    !      sig1=sig1+epa1
                                                                                    !      sig2=sig2+epa2
                                                                                    !      sig3=sig3-epa1-epa2
                                                                                    !      sig4=sig4+epa4
                                                                                    !      sig5=sig5+epa5
                                                                                    !      sig6=sig6+epa6

                                                                                    !      endif
      return
      end 
      
      
      subroutine feit(capan,xj1tr,xj2dtr,t,alpha,gama,beta,theta,w,dDyna,r,bulk,shear,xj1,capa,sig1,sig2,sig3,sig4,sig5,sig6,iter,nocon)
      use mcm_database  
      real(kind=real_acc) :: bulk,capan,capa
      real(kind=real_acc) :: sig1,sig2,sig3,sig4,sig5,sig6
      parameter (tol=1.0e-8,maxit=100)
!c     implicit double precision (a-h,o-z)  
!c
!
      iter=1
      nocon=1
      capa=capan
      xj1=xj1tr
      xj2d=xj2dtr
      s1=sig1-xj1tr/3.
      s2=sig2-xj1tr/3.
      s3=sig3-xj1tr/3.
      s4=sig4
      s5=sig5
      s6=sig6
   10 continue
      xk=capa+r*(alpha-gama*exp(-beta*capa)+theta*capa)
      dhdk=w*dDyna*exp(-dDyna*xk)*(1.+r*(theta+gama*beta*exp(-beta*capa)))
      omega=theta+gama*beta*exp(-beta*xj1)
      dcd=shear+9.*bulk*omega**2
      f=xj2d-(alpha-gama*exp(-beta*xj1)+theta*xj1)
      delcap=-3.*omega*f/(dcd*dhdk)
      delam=-dhdk*delcap/(3.*omega)
      capa=capa+delcap
      xk=capa+r*(alpha-gama*exp(-beta*capa)+theta*capa)
      fac1=shear/xj2d
      fac2=3.*bulk*omega
      sig1=sig1-delam*(fac1*s1-fac2)
      sig2=sig2-delam*(fac1*s2-fac2)
      sig3=sig3-delam*(fac1*s3-fac2)
      sig4=sig4-delam*fac1*s4
      sig5=sig5-delam*fac1*s5
      sig6=sig6-delam*fac1*s6
      xj1=sig1+sig2+sig3
      s1=sig1-xj1/3.
      s2=sig2-xj1/3.
      s3=sig3-xj1/3.
      s4=sig4
      s5=sig5
      s6=sig6
      xj2d= sqrt(.5*(s1**2+s2**2+s3**2)+s4**2+s5**2+s6**2)
      f=xj2d-(alpha-gama*exp(-beta*xj1)+theta*xj1)
      if(abs(f).le.tol) return
      if( iter.lt.maxit ) then
      iter=iter+1
      go to 10
      else
      write(13,1000) iter,f
      write( *,1000) iter,f
 1000 format(/5x,'cap model failure mode iterations unconverged',    /10x,'iterations: ',i5,5x,'fe:  ',e14.4)
      stop 'feit1'
      endif
      end
      
      subroutine fcit(capan,xj1tr,xj2dtr,t,alpha,gama,beta,theta,w,dDyna,r,bulk,shear,xj1,capa,sig1,sig2,sig3,sig4,sig5,sig6,iter,nocon)    
      use mcm_database  
      real(kind=real_acc) :: bulk,capan,capa
      real(kind=real_acc) :: sig1,sig2,sig3,sig4,sig5,sig6   
      parameter (tol=1.0e-8,maxit=100)
!c     implicit double precision (a-h,o-z)                               
!c
      unit=1.0
      iter=1
      nocon=1
      capa=capan
      xj1=xj1tr
      xj2d=xj2dtr
      s1=sig1-xj1tr/3.
      s2=sig2-xj1tr/3.
      s3=sig3-xj1tr/3.
      s4=sig4
      s5=sig5
      s6=sig6
   10 continue
      xk=capa+r*(alpha-gama*exp(-beta*capa)+theta*capa)
      dhdk=w*dDyna*exp(-dDyna*xk)*(1.+r*(theta+gama*beta*exp(-beta*capa)))
      if(capa.ge.0.0) then
      dfdk=-2.*(abs(xj1-capa)/r**2 + (alpha-gama*exp(-beta*capa)   +theta*capa)*(theta+gama*beta*exp(-beta*capa)))
      else
      dfdk=-2.*xk*(1.+r*(theta+gama*beta*exp(-beta*capa)))/r**2
      endif
      elcap=max(capa,0.0*unit)
      omega=min(-2.*abs(xj1-elcap)/r**2,-1.e-15*unit)
      dcd=4.*shear*xj2d**2 + 9.*bulk*(2.*(xj1-elcap)/r**2)**2
      f=xj2d**2 + (xj1-elcap)**2/r**2 - (xk-elcap)**2/r**2
      delcap=-f/(dfdk+(dcd*dhdk)/(3.*omega))
      delam=-dhdk*delcap/(3.*omega)
      capa=capa+delcap
      xk=capa+r*(alpha-gama*exp(-beta*capa)+theta*capa)
!c     elcap=max(capa,0.0*unit)
      fac1=2.*shear
      fac2=6.*bulk*abs(xj1-elcap)/r**2
      sig1=sig1-delam*(fac1*s1+fac2)
      sig2=sig2-delam*(fac1*s2+fac2)
      sig3=sig3-delam*(fac1*s3+fac2)
      sig4=sig4-delam*fac1*s4
      sig5=sig5-delam*fac1*s5
      sig6=sig6-delam*fac1*s6
      xj1=sig1+sig2+sig3
      s1=sig1-xj1/3.
      s2=sig2-xj1/3.
      s3=sig3-xj1/3.
      s4=sig4
      s5=sig5
      s6=sig6
      xj2d= sqrt(.5*(s1**2+s2**2+s3**2)+s4**2+s5**2+s6**2)
      f=xj2d**2 + (xj1-elcap)**2/r**2 - (xk-elcap)**2/r**2
      if(abs(f).le.tol) then
      if(xj1.lt.elcap) then
      xj1=elcap+abs(xj1-elcap)
      sig1=s1+xj1/3.
      sig2=s2+xj1/3.
      sig3=s3+xj1/3.
      endif
      return
      endif
      if( iter.lt.maxit ) then
      iter=iter+1
      go to 10
      else
      write(13,1000) iter,f
      write( *,1000) iter,f
 1000 format(/5x,'cap model cap mode iterations unconverged',       /10x,'iterations: ',i5,5x,'fe:  ',e14.4)
      stop 'fcit1'
      endif
      end