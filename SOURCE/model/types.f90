!***********************************************************************
  
      subroutine type1(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 1 : FAN or PUMP MODEL
!
!       Required units: pressures in kPa; specific heats in kJ/kg-C;
!                       densities in kg/m3; flow rate in kg/s.
!       Head curve funp(z), efficiency curve funn(z), and flow kW are
!       dimensionless.
!
!       Updated :  June 29, 1992
!       Upadted to Fortran 90  March 23, 2007 Cheol Park, NIST
!
!***********************************************************************

      use modsim_head
      implicit none
      integer, parameter                :: ni=4, no=3, np=12, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real           :: n,kw
      integer        :: i,mode
      real           :: z,funp,funn,w,p2,t1,d,cp,&
                        rho,p1,eff,t2,power
      real           :: rhoa=1.2,rhow=1000.0,cpa=1.0,cpw=4.180

      funp(z)=par_v(1)+z*(par_v(2)+z*(par_v(3)+&
              z*(par_v(4)+z*par_v(5))))
      funn(z)=par_v(6)+z*(par_v(7)+z*(par_v(8)+&
              z*(par_v(9)+z*par_v(10))))
  
      w=    xin(1)
      p2=   xin(2)
      n=    xin(3)
      t1=   xin(4)
  
      d=    par_v(11)
      mode= nint(par_v(12))
  
      if(mode == 1) then
         cp=cpa
         rho=rhoa
      else
         cp=cpw
         rho=rhow
      endif
  
      if( n > 0.0  .and.  w > 0.0 ) then
         kw=w/(rho*d**3*n)
         p1=p2-0.001*funp(kw)*rho*(d*n)**2
         eff=funn(kw)
         t2=t1+(p2-p1)*((1.0/eff)-1.0)/(rho*cp)
         power=(p2-p1)*w/(rho*eff)
      else
         p1=p2
         t2=t1
         power=0.0
      endif
  
      yout(1)= p1
      yout(2)= t2
      yout(3)= power
  
      do i=1,3
         iostat(i)=1
      end do

      return
      end subroutine type1
  
!***********************************************************************
  
      subroutine type2(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 2 : CONUIT MODEL (PIPE or DUCT)
!
!   flow:   kg/s            (Fluid mass flow rate)
!   tin:    C               (Fluid inlet temperature)
!   tamb:   C               (Ambient temperature)
!   p2:     kPa             (Pressure at conduit outlet)
!   hia:    kW/C            (Heat transfer coeff*area, inside surface)
!   hoa:    kW/C            (Heat transfer coeff*area, outside surface)
!   cm:     kJ/C            (Thermal capacitance of conduit material)
!   vol:    m3              (Volume of conduit)
!   k:      0.001/(kg-m)    (Flow Resistance)
!   h:      m               (height of outlet relative to inlet)
!                  (sign convention: positive H --> uphill flow)
!
!  The last parameter, mode, has the following meaning:
!   mode = -2 => water, detailed dynamics (one diff.eq.)
!   mode = -1 =>   air, detailed dynamics (one diff.eq.)
!   mode =  1 =>   air,   simple dynamics (no  diff.eq.)
!   mode =  2 => water,   simple dynamics (no  diff.eq.)
!
!***********************************************************************
  
      use modsim_head
      implicit none

      integer, parameter                :: ni=5, no=2, np=7, ns=11,&
                                           ndeq=1
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      logical        :: frzto
      integer        :: mode
      real           :: k
      real           :: flow,p2,tin,tamb,tout,hia,hoa,&
                        cm,vol,h,cp,rho,ua,told,p1,tss,dt,&
                        taux,t,dtindt,alpha,tauc,alag,&
                        dtodt,delay,gamma
      real           :: grav=9.80665
      real           :: rhoa=1.2,rhow=1000.0,cpa=1.0,cpw=4.180

      flow=  xin(1)
      p2=    xin(2)
      tin=   xin(3)
      tamb=  xin(4)
      tout=  xin(5)
  
      frzto=(iostat(1)<-1).or.(iostat(1)==2)

      hia=   par_v(1)
      hoa=   par_v(2)
      cm=    par_v(3)
      vol=   par_v(4)
      k=     par_v(5)
      h=     par_v(6)
      mode=  nint(par_v(7))
  
      if(mode == 1) then
         cp=cpa
         rho=rhoa
      else
         cp=cpw
         rho=rhow
      endif
  
      if((hia>0.0).and.(hoa>0.0))  then
         ua=(hia*hoa)/(hia+hoa)
      else
         ua=0.0
      endif
  
      if(itime ==1) then
         if(init==0) then
            saved_v(9)=0.0
            saved_v(11)=tin
         endif
         if(saved_v(9)>time) then
            saved_v(9)=0.0
         endif
      endif
  
      if(time>saved_v(9)) then
        saved_v(10)=saved_v(11)
      endif
      told=saved_v(10)
  
      p1=p2+(rho*grav*h*0.001)+k*flow*abs(flow)
  
!     Mmultiplication by 0.001 converts from Pa to kPa.
  
      if(flow <= 0.00001) then
         if(abs(mode)==mode)then
            yout(1)=tin
            saved_v(11)=told
         else
            yout(1)=0.0
            saved_v(11)=tin
         endif
         yout(2)=p1
         iostat(1)=1
         iostat(2)=1
         return
      endif
  
      gamma=ua/(flow*cp)
!      a=ua/(flow*cp)
      dt=tin-tamb
      if(gamma<=20.0 .and. abs(dt)>=1.0e-10) then
         tss=tamb+dt*exp(-gamma)
      else
         tss=tamb
      endif
  
      taux=rho*vol/flow
      if(abs(mode)==mode) then
  
!     Simplified dynamics (no differential equations)

250      t=tss
         if((cm>0.0) .and. (hia>0.0)) then
            gamma=flow*cp*tstep*(1.+hoa/hia)/cm
            dt=told-tss
            if(gamma<=20.0 .and. abs(dt)>=1.0e-10) then
               t=tss+dt*exp(-gamma)
            endif
         endif
         tout=delay(t,taux,saved_v(1),saved_v(7),saved_v(8),saved_v(9))
         yout(1)=tout
         saved_v(11)=t
  
      else

!     Detailed dynamics (one differential equation)

         dtindt=(tin-told)/tstep
         alpha=0.5*gamma*hia/hoa
         tauc=cm/(hia+hoa)
         alag=tss-tamb+tauc*exp(-gamma-alpha)*dtindt
         dtodt=(exp(-alpha)/tauc)*(tamb-tout+delay(alag,taux,saved_v(1),&
           saved_v(7),saved_v(8),saved_v(9)))
         yout(1)=dtodt
         saved_v(11)=tin
  
      endif
  
      yout(2)=p1
      iostat(2)=1
      iostat(1)=0

      if((abs(tss-tout)<=rtolx*abs(tss)+atolx).and.frzto) then
         iostat(1)=1
      end if
  
      return
      end subroutine type2
!***********************************************************************
  
      subroutine type3(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 3 : INLET CONDUIT MODEL ( PIPE or DUCT)
!
!   flow:   kg/s            (Fluid mass flow rate)
!   tin:    C               (Fluid inlet temperature)
!   tamb:   C               (Ambient temperature)
!   p2:     kPa             (Pressure at conduit outlet)
!   hia:    kW/C            (Heat transfer coeff*area, inside surface)
!   hoa:    kW/C            (Heat transfer coeff*area, outside surface)
!   cm:     kJ/C            (Thermal capacitance of conduit material)
!   vol:    m3              (Volume of conduit)
!   k:      0.001/(kg-m)    (Flow Resistance)
!   h:      m               (height of outlet relative to inlet)
!                  (sign convention: positive h --> uphill flow)
!
!  The last parameter, mode, has the following meaning:
!   mode = -2 => water, detailed dynamics (one diff.eq.)
!   mode = -1 =>   air, detailed dynamics (one diff.eq.)
!   mode =  1 =>   air,   simple dynamics (no  diff.eq.)
!   mode =  2 => water,   simple dynamics (no  diff.eq.)
!
!***********************************************************************
  
      use modsim_head
      implicit none
      integer, parameter                :: ni=5, no=2, np=7, ns=11,&
                                           ndeq=1
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real           :: k
      logical        :: frzto
      integer        :: mode
      real           :: p1,p2,tin,tamb,tout,hia,hoa,cm,vol,h,cp,rho,ua,&
                        told,s,dp,flow,tss,dt,tau,dtindt,alpha,tauc,alag,&
                        dtodt,t,delay,gamma
      real           :: grav=9.80665
      real           :: rhoa=1.2,rhow=1000.0,cpa=1.0,cpw=4.180

      p1=   xin(1)
      p2=   xin(2)
      tin=  xin(3)
      tamb= xin(4)
      tout= xin(5)
  
      frzto= (iostat(2)<-1).or.(iostat(2)==2)
  
      hia=  par_v(1)
      hoa=  par_v(2)
      cm=   par_v(3)
      vol=  par_v(4)
      k=    par_v(5)
      h=    par_v(6)
      mode= nint(par_v(7))
  
      if(abs(mode) == 1) then
         cp=cpa
         rho=rhoa
      else
         cp=cpw
         rho=rhow
      endif
  
      if((hia>0.0).and.(hoa>0.0))  then
         ua=(hia*hoa)/(hia+hoa)
      else
         ua=0.0
      endif
  
      if(itime <=1) then
         if(init==0) then
            saved_v(9)=0.0
            saved_v(11)=tin
         endif
         if(saved_v(9)>time) then
            saved_v(9)=0.0
         endif
      endif

 50   if(time>saved_v(9)) then
        saved_v(10)=saved_v(11)
      endif

      told=saved_v(10)
      dp=p1-p2-rhow*grav*h*0.001

!     Multiplication by 0.001 converts Pa to kPa.

      if(dp/=0.0) then
        s=sign(1.0,dp)
      else
        s=1.0
      endif

      flow=s*sqrt(abs(dp)/k)
      if(flow<=0.00001) then
         if(abs(mode)==mode)then
            yout(1)=tin
            saved_v(11)=told
         else
            yout(1)=0.0
            saved_v(11)=tin
         endif

         yout(2)=flow
         iostat(1)=1
         iostat(2)=1
         return
      endif
  
100   tss=tamb
      gamma=ua/(flow*cp)

      if(gamma<=20.0) then
        dt=tin-tamb
        if(abs(dt)>=1.0e-10) then
          tss=tamb+dt*exp(-gamma)
        endif
      endif

      tau=rho*vol/flow

      if(abs(mode)/=mode) then

!       Detailed dynamics (one differential equation)
                                                                              
        dtindt=(tin-told)/tstep                                               
        alpha=0.5*gamma*hia/hoa                                               
        tauc=cm/(hia+hoa)                                                     
        alag=tss-tamb+tauc*exp(-gamma-alpha)*dtindt                           
        dtodt=(exp(-alpha)/tauc)*(tamb-tout+delay(alag,tau,saved_v(1),&
          saved_v(7),saved_v(8),saved_v(9)))
        yout(1)=dtodt                                                         
        saved_v(11)=tin
      else

!       Simplified dynamics (no differential equations)
                                                                           
250     t=tss                                                              
        if((cm>0.0).and.(hia>0.0)) then
          gamma=flow*cp*tstep*(1.0+hoa/hia)/cm
          if(gamma<=20.0) then
            dt=told-tss
            if(abs(dt)>=1.0e-10) then
              t=tss+dt*exp(-gamma)
            endif
          endif
        endif

        tout=delay(t,tau,saved_v(1),saved_v(7),saved_v(8),saved_v(9))      
        yout(1)=tout                                                       
        saved_v(11)=t                                                      
      endif

      yout(2)=flow
      iostat(2)=1
      iostat(1)=0
      if((abs(tss-tout)<=rtolx*abs(tss)+atolx).and.frzto)&
        iostat(1)=1
  
      return
      end subroutine type3
!***********************************************************************
  
      subroutine type4(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 4 : FLOW MERGE MODEL
!
!***********************************************************************
  
      implicit none
      integer, parameter                :: ni=5, no=4, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      integer        :: i
      real           :: k
      real           :: w1,w2,p3,t1,t2,w3,p1,p2,t3

      w1=  xin(1)
      w2=  xin(2)
      p3=  xin(3)
      t1=  xin(4)
      t2=  xin(5)
  
      k=par_v(1)
  
      w3=w1+w2
      p1=p3+0.5*k*(w1*abs(w1)+w3*abs(w3))
      p2=p3+0.5*k*(w2*abs(w2)+w3*abs(w3))

      if(abs(w3)>1.0e-12) then
         t3=(w1*t1+w2*t2)/w3
      else
         t3=0.5*(t1+t2)
      endif
  
      yout(1)= w3
      yout(2)= p1
      yout(3)= p2
      yout(4)= t3
  
      do i=1,4
         iostat(i)=1
      end do
  
      return
      end subroutine type4
!***********************************************************************
  
      subroutine type5(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 5 : DAMPER or VALVE MODEL
!
!***********************************************************************
  
      implicit none
      integer, parameter                :: ni=3, no=1, np=4, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: k
      real           :: w,p2,c,xl,xm,rl,re,r,p1

      w=  xin(1)
      p2= xin(2)
      c=  xin(3)
  
      k=  par_v(1)
      xl= par_v(2)
      xm= par_v(3)
  
      if(par_v(4)/=0.0) then
        c=1.0-c
      endif
  
      if(c<0.0) then
        c=0.0
      elseif(c>1.0) then
        c=1.0
      endif

      rl=k/((1.0-xl)*c+xl)**2
      re=k*(xl**(2.0*c-2.0))
      r=xm*rl+(1.0-xm)*re
      p1=p2+r*w*abs(w)

      yout(1)=p1
      iostat(1)=1
  
      return
      end subroutine type5
  
 !***********************************************************************
  
      subroutine type6(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 6 : FLOW SPLIT MODEL
!
!***********************************************************************
  
      use modsim_head
      implicit none
      integer, parameter                :: ni=3, no=3, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      integer        :: i
      real           :: k
      real           :: w1,p2,p3,dp,dw,discr,w2,sp,w3,p1

      w1=  xin(1)
      p2=  xin(2)
      p3=  xin(3)
  
      k=   par_v(1)
  
      dp=p3-p2
      dw=2.0*dp/(k*abs(w1))
      discr=4.0*abs(dp)/k-w1**2

      if(dp/=0.0) then
        sp=sign(1.0,dp)
      else
        sp=1.0
      endif

      if(discr>w1**2) then
        w2=0.5*(w1+sp*sqrt(discr))
      else
        w2=0.5*(w1+dw)
      endif

      w3=w1-w2
      p1=p2+0.5*k*(w1*abs(w1)+w2*abs(w2))
  
      yout(1)=  w2
      yout(2)=  w3
      yout(3)=  p1
  
      do i=1,3
         iostat(i)=1
      end do
  
      return
      end subroutine type6
!***********************************************************************
!
      subroutine type7(xin,yout,par_v,saved_v,iostat)
!
! ----------------------------------------------------------------------
!
!     TYPE 7 : FIRST ORDER TEMPERATURE SENSOR MODEL
!
!***********************************************************************

      use modsim_head
      implicit none
      integer, parameter                :: ni=2, no=1, np=3, ns=3,&
                                           ndeq=1
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: ti,cto,tcon,tgain,cti,dtodt,dt,ctout,a1

      ti   =xin(1)
      cto  =xin(2)
  
      tcon =par_v(1)
      tzero=par_v(2)
      tgain=par_v(3)
!
      cti=(ti-tzero)/tgain
  
      if(tcon>=1.0) then
         dtodt  = (cti-cto)/tcon
         yout(1) = dtodt
         if(((iostat(1)<-1).or.(iostat(1)==2)) .and.&
           (abs(cti-cto)<=rtolx*abs(cto)+atolx)) then
            iostat(1)=1
         else
            iostat(1)=0
         endif
         return
      endif
  
      if(itime == 1) then
         if (init == 0 ) then
            saved_v(1) = 0.0
            saved_v(2) = cto
         endif
         if ( saved_v(1) > time ) then
            saved_v(1) = 0.0
         endif
      endif
  
      if(time > saved_v(1))  then
         saved_v(3)=saved_v(2)
      endif
      saved_v(1) = time
  
      a1 = tcon/tstep
      dt = cti - saved_v(3)
      if(a1 < 0.05 .or. abs(dt) < 1.0e-10) then
         ctout=cti-dt*exp(-1.0/a1)
      else
         ctout = cti
      endif
  
      yout(1)   = ctout
      saved_v(2) = ctout
  
      if(((iostat(1)<-1).or.(iostat(1)==2)) .and.&
        (abs(cti-ctout)<=rtolx*abs(ctout)+atolx)) then
         iostat(1)=1
      else
         iostat(1)=0
      endif
  
      return
      end subroutine type7

!***********************************************************************
  
      subroutine type8(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 8 : PROPORTIONAL PLUS INTEGRAL CONTROLLER  MODEL
!
!***********************************************************************
  
      use modsim_head
      implicit none
      integer, parameter                :: ni=4, no=2, np=3, ns=3,&
                                           ndeq=2
      real,dimension(ni),intent(in out) :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: ctm,cts,cint,c,gainp,css,gaini,tcon,e,cpss,didt,&
                        dcdt,out2,dc

      ctm=  xin(1)
      cts=  xin(2)
      cint= xin(3)
      c=    xin(4)
  
      gainp=par_v(1)
      gaini=par_v(2)
      tcon= par_v(3)
  
!     Initialize at beginning of simulation

      if(itime<=1) then
         if(init==0 .or. saved_v(1)>time) then
           saved_v(1)=0.0
         endif

         if(init==0) then
           saved_v(2)=c
         endif

         if(gaini<=0.0) then
           xin(3)=0.
         endif
      endif

      if(time>saved_v(1)) then
        saved_v(3)=saved_v(2)
      endif

      saved_v(1)=time

      e=cts-ctm
      cpss=gainp*e
      css=cpss+cint

!     Limit c and css to values from 0.0 to 1.0

      if(css<0.0) then
        css=0.0
      elseif(css>1.0) then
        css=1.0
      endif

      if(c<0.0) then
        c=0.0
      elseif(c>1.0) then
        c=1.0
      endif

      didt=0.0
      if(gaini>0.0)then
         didt=gaini*e
         if(((cint>=1.0).and.(didt>0.0)).or.&
           ((cint<=0.0).and.(didt<0.0))) then
           didt=0.0
         endif

!        **Don't integrate error if within tolerances

         if(abs(e)<rtolx*abs(cts)+atolx) then
           didt=0.0
         endif
      endif

      if(tcon>=1.0) then
         dcdt=(css-c)/tcon
         out2=dcdt
      else
         dc=css-saved_v(3)
         out2=css
         if((tcon/tstep>=0.05).and.(abs(dc)>=1.0e-10)) then
            out2=css-dc*exp(-tstep/tcon)
            saved_v(2)=out2
         endif
      endif

      yout(1)=didt
      yout(2)=out2
      iostat(1)=0
      iostat(2)=0

!     **Disable old conditions for freezing of controller
!     **        if(abs(cpss)<=rtolx) iostat(1)=1
!     **        if(abs(css-c)<=rtolx*abs(c)+atolx) iostat(2)=1
!     **Try alternate condition for freezing integral portion

      if(abs(e)<=rtolx*abs(cts)+atolx) iostat(1)=1
!     **mess with freezing conditions at your own risk.
  
      return
      end subroutine type8

!***********************************************************************
  
      subroutine type9(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 9 : LINEAR VALVE MODEL
!
!***********************************************************************
  
      use modsim_head
      implicit none

      integer, parameter                :: ni=4, no=3, np=4, ns=5,&
                                           ndeq=1
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: k
      logical        :: con
      real           :: c,cv,tcon,cl,hys,dcvdt,dc,s,cvact,p1,p2,w,&
                        hyster,a1

      p2=  xin(1)
      w=   xin(2)
      c=   xin(3)
      cv=  xin(4)
  
      k=   par_v(1)
      tcon=par_v(2)
      cl=  par_v(3)
      hys= par_v(4)
  
      con=(iostat(3)<-1).or.(iostat(3)==2)
      if(c<0.0) then
        c=0.0
      elseif(c>1.0) then
        c=1.0
      endif

      if(cv<0.0) then
        cv=0.0
      elseif(cv>1.0) then
        cv=1.0
      endif

      if(tcon>=1.0) then
        dcvdt=(c-cv)/tcon
      else
        if(itime<=1) then
          if(init==0) then                              
            saved_v(4)=cv                               
          endif                                         
          if(init==0 .or. saved_v(3)>time) then         
            saved_v(3)=0.                               
          endif                                         
        endif                                           
                                                        
3000    if(time>saved_v(3)) then                        
          saved_v(5)=saved_v(4)                         
        endif

        dcvdt=c                                         
        cv=c                                            
        a1=tcon/tstep

        if(a1>=0.05) then
          dc=c-saved_v(5)
          if(abs(dc)>=1.0e-10) then
            cv=c-dc*exp(-1./a1)
            dcvdt=cv
          endif
        endif
      endif

      if(w/=0.0) then
        s=sign(1.0,w)
      else
        s=1.0
      endif

      cvact=hyster(cv,saved_v(1),saved_v(2),saved_v(3),hys)
      p1=p2+s*k*(w/((1.-cl)*cvact+cl))**2
  
      saved_v(4)=cv
  
      yout(1)= dcvdt
      yout(2)= p1
      yout(3)= cvact
  

      if((abs(c-cv)<=rtolx*abs(cv)+atolx).and.con) then
        iostat(1)=1
      else
        iostat(1)=0
      endif

      iostat(2)=1
      iostat(3)=1
  
      return
      end subroutine type9
  
!***********************************************************************
  
      subroutine type10(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 10 : HOT WATER COIL MODEL (CROSSFLOW)
!
!***********************************************************************
  
      use modsim_head
      implicit none

      integer, parameter                :: ni=8, no=4, np=5, ns=0,&
                                           ndeq=2
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      integer        :: i
      real           :: k1,k2,ntu
      logical        :: flo1,flo2,frzde
      real           :: eta,r,e1,e2,efec,t2s,t4s,dt2dt,dt4dt,t3,t4,ua,&
                        vol,tcon,s,p1,p3,c1,c2,rate,cmin,cmax,p2,w1,t1,&
                        t2,p4,w2,a1
      real           :: rhoa=1.2,rhow=1000.0,cpa=1.0,cpw=4.180

      p2= xin(1)
      w1= xin(2)
      t1= xin(3)
      t2= xin(4)
      p4= xin(5)
      w2= xin(6)
      t3= xin(7)
      t4= xin(8)
  
      ua=  par_v(1)
      k1=  par_v(2)
      k2=  par_v(3)
      vol= par_v(4)
      tcon=par_v(5)
  
      flo1=(iostat(2)<-1).or.(iostat(2)==2)
      flo2=(iostat(3)<-1).or.(iostat(3)==2)
      frzde=flo1.and.flo2

      if(w1/=0.0) then
        s=sign(1.0,w1)
      else
        s=1.0
      endif

      p1=p2+s*k1*w1**2

      if(w2/=0.0) then
        s=sign(1.0,w2)
      else
        s=1.0
      endif

      p3=p4+s*k2*w2**2
      c1=cpa*w1
      c2=cpw*w2

      if(tcon>0.0)then
        rate=rate+1.0/tcon
      else
        rate=0.0
      endif

      if(vol>0.0) then
        rate=rate+w2/(rhow*vol)
      endif

      cmin=amin1(c1,c2)
      cmax=amax1(c1,c2)

      if((c1<=0.0).or.(c2<=0.0)) then
        t2s=t1
        t4s=t3
      else
        ntu=ua/cmin
        eta=ntu**0.22                
        r=cmin/cmax                  
        e1=0.0                       
        e2=0.0                       
        a1=r*ntu/eta                 
                                     
        if(a1<=20.0) then
          e1=exp(-a1)
        endif                        
        a1=eta*(1.0-e1)/r

        if(a1<=20.0) then
          e2=exp(-a1)
        endif                        
        efec=1.0-e2

        t2s=t1+efec*cmin*(t3-t1)/c1
        t4s=t3-c1*(t2s-t1)/c2
      endif

      if(rate<=0.0) then
        yout(1)=t2s
        yout(2)=t4s              
        yout(3)=p1               
        yout(4)=p3               
                                 
        do i=1,4                 
           iostat(1)=1           
        end do                   
      else
        dt2dt=(t2s-t2)*rate
        dt4dt=(t4s-t4)*rate                                               
                                                                          
        yout(1)=dt2dt                                                     
        yout(2)=dt4dt                                                     
        yout(3)=p1                                                        
        yout(4)=p3                                                        

        if((abs(t2s-t2)<=rtolx*abs(t2)+atolx).and.frzde) then
          iostat(1)=1
        else
          iostat(1)=0
        endif

        if((abs(t4s-t4)<=rtolx*abs(t4)+atolx).and.frzde) then
          iostat(2)=1
        else
          iostat(2)=0
        endif

        iostat(3)=1                                                       
        iostat(4)=1                                                       
      endif

      return
      end subroutine type10
  
!***********************************************************************
  
      subroutine type11(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 11 : HOT WATER COIL MODEL (CROSSFLOW)
!
!***********************************************************************
  
      use modsim_head
      implicit none

      integer, parameter                :: ni=9, no=5, np=8, ns=11,&
                                           ndeq=3
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: ka,kw,ntu,lag
      logical        :: flo1,tem1,frzde
      real           :: e2,efec,taoss,twoss,dtwidt,taux,b1,b2,b3,b4,phi,&
                        e,dphi,z,ab1,aphi,a2,a3,two,td,haa,haw,vol,cm,s,&
                        pai,pwi,rate,dtaodt,dtwodt,ca,cw,cmin,cmax,ua,&
                        eta,r,e1,tauw,pao,wa,tai,tao,pwo,ww,twi,dtddt,&
                        a1,delay
      real           :: rhoa=1.2,rhow=1000.0,cpa=1.0,cpw=4.180

      pao=  xin(1)
      wa=   xin(2)
      tai=  xin(3)
      tao=  xin(4)
      pwo=  xin(5)
      ww=   xin(6)
      twi=  xin(7)
      two=  xin(8)
      td=   xin(9)
  
      haa=par_v(1)*wa**par_v(2)
      haw=par_v(3)*ww**par_v(4)
      ka= par_v(5)
      kw= par_v(6)
      vol=par_v(7)
      cm= par_v(8)
  
      if(itime<=1) then
  
        if(init==0 .or. saved_v(9)>time) then
          saved_v(9)=0.0                          
        endif                                     
                                                  
        if(init==0) then                          
          saved_v(11)=twi                         
        endif

      endif

      if(time>saved_v(9)) then
        saved_v(10)=saved_v(11)
      endif
  
      flo1=(iostat(2)<-1).or.(iostat(2)==2)
      tem1=(iostat(3)<-1).or.(iostat(3)==2)
      frzde=flo1.and.tem1

      if(wa/=0.0) then
        s=sign(1.0,wa)
      else
        s=1.0
      endif

      pai=pao+s*ka*wa**2

      if(ww/=0.0) then
        s=sign(1.0,ww)
      else
        s=1.0
      endif

      pwi=pwo+s*kw*ww**2
      rate=1.0/tstep
  
      if(wa>0.0 .and. ww>0.0) then
        ca=cpa*wa
        cw=cpw*ww                                                           
        cmin=amin1(ca,cw)                                                   
        cmax=amax1(ca,cw)                                                   
        ua=haa*haw/(haa+haw)                                                
                                                                            
!       Calculate effectiveness: unmixed both sides, approx. equation.      
                                                                            
        ntu=ua/cmin                                                         
        eta=ntu**0.22                                                       
        r=cmin/cmax                                                         
        e1=0.0                                                              
        e2=0.0                                                              
        a1=r*ntu/eta                                                        

        if(a1<=20.0) then
          e1=exp(-a1)
        endif
        a1=eta*(1.0-e1)/r

        if(a1<=20.0) then
          e2=exp(-a1)
        endif
        efec=1.0-e2
                                                                            
!       The next 3 lines are effectiveness for fully mixed fluids           
                                                                            
!**       if(ntu<=20.0) e1=exp(-ntu)                                        
!**       if(r*ntu<=20.0) e2=exp(-r*ntu)                                    
!**       efec=ntu/(ntu/(1.-e1)+r*ntu/(1.-e2)-1.0)                          
        taoss=tai+efec*cmin*(twi-tai)/ca                                    
        twoss=twi-ca*(taoss-tai)/cw                                         
                                                                            
        dtwidt=(twi-saved_v(10))*rate                                       
        taux=rhow*vol/ww                                                    
        b1=haw/cw                                                           
        b2=taux*haa/cm                                                      
        b3=taux*haw/cm                                                      
        b4=haa/ca                                                           
        a1=b3+2.0*b2/(2.0+b4)                                               
        phi=b1-b1*b3/a1                                                     
        e=exp(-phi)                                                         

        dphi=1.0+b1*b3/(a1*a1)                                              
        z=(a1+b1)/(a1*phi)-dphi*e/(1.0-e)                                   
        e1=e/(1.0-e)                                                        
        ab1=a1+b1                                                           
        aphi=a1*phi                                                         
        a2=ab1/aphi-dphi*e1                                                 
        a3=a2*a2-(ab1*ab1-aphi)/(aphi*aphi)&                                
         +(dphi*(ab1*(2.0+phi)+phi*(2.0-phi))-2.0*phi)*e1/(2.0*aphi)        
        dtddt=(taoss-tao)/taux                                              
        dtaodt=(td-a2*tao)/(a3*taux)                                        
                                                                            
        tauw=taux*exp(b1*b3/(2.0*a1))/a1                                    
        lag=dtwidt*exp(-b1)*tauw+(twoss-tai)                                
        dtwodt=((tai-two)+delay(lag,taux,saved_v(1),saved_v(7),&            
        saved_v(8),saved_v(9)))/tauw
      else
        dtaodt=0.0
        dtddt=0.0        
        dtwodt=0.0       
      endif

      saved_v(11)=twi

      yout(1)=dtaodt
      yout(2)=dtwodt
      yout(3)=dtddt
      yout(4)=pai
      yout(5)=pwi

      iostat(4)=1
      iostat(5)=1
      iostat(1)=0
      iostat(2)=0
      iostat(3)=0

!     Either all three d.e.'s freeze, or none of them freeze.

      if((abs(dtaodt) <= atolx) .and.&
       (abs(dtddt)  <= atolx) .and.&
       (abs(dtwodt) <= atolx) .and. frzde) then
         iostat(1)=1
         iostat(2)=1
         iostat(3)=1
      endif

      return
      end subroutine type11
  
!***********************************************************************
!
      subroutine type12(xin,yout,par_v,saved_v,iostat)
!
! ----------------------------------------------------------------------
!
!     TYPE 12 : COOLING and DEHUMIDIFYING COIL
!
!
!       "Fortran IV program to simulate cooling and dehumidifying
!       finned-tube multi-row heat exchangers", A.H. Elmahdy and
!       G.P. Mitalas, Computer program no. 43, Division of Building
!       Research, National Research Council of Canada, Ottawa,
!       March 1977.
!
!       Adapted for use with modsim, January 1984, by D.R. Clark,
!       Center for Building Technology, NBS.
!
!       Dynamics added to steady state model, though wet/dry condition
!       of coil is determined from steady state conditions.
!       Numerous pointless changes made in nomenclature.
!
!-----------------------------------------------------------------------
!
!    This program models circular finned or continuous finned heat
!    exchangers with four or more rows in counterflow crossflow
!    configuration.
!
!    inputs:
!
!     1) fmr:        water mass flow rate            (kg/s)
!     2) twi:        entering water temperature      (C)
!     3) fma:        dry air mass flow rate          (kg/s)
!     4) taidb:      entering air dry bulb temp.     (C)
!     5) wai:        entering air humidity ratio     (-)
!     6) two         exit water temp. [from yout(1)] (C)
!     7) tao         exit air temp. [from yout(2)]   (C)
!     8) wao         exit humid. ratio [from yout(3)](-)
!
!    outputs (all steady state values):
!
!     1) two:        leaving water temperature       (C)
!     2) tao:        leaving air dry bulb temp.      (C)
!     3) wao:        leaving air humidity ratio      (-)
!     4) qt:         total cooling load              (kW)
!     5) qs:         sensible cooling load           (kW)
!     6) fwet:       wet surface area as a fraction
!                    of total surface area           (-)
!
!    parameters:
!
!     1) icoil:      coil type:
!                     0 = flat continuous fins
!                     1 = circular fins              (-)
!     2) ap:         primary surface area            (sq. m)
!     3) as:         secondary surface area          (sq. m)
!     4) ai:         internal surface area           (sq. m)
!     5) flfa:       min. flow area/face area        (-)
!     6) fcon:       fin material conductivity       (kW/m.K)
!     7) fa:         coil face area                  (sq. m)
!     8) fpi:        number of fins per centimeter   (1/cm)
!     9) npr:        number of tubes per row         (-)
!    10) nor:        number of rows                  (-)
!    11) tdo:        outside tube diameter           (m)
!    12) tdi:        inside tube diameter            (m)
!    13) fth:        fin thickness                   (m)
!    14) cm:         thermal capacitance of metal    (kJ/K)
!    15) sl:         tube longitudinal spacing       (m)
!    16) fd:         fin diameter (if icoil = 1)     (m)
!    17) coid:       coil depth                      (m)
!    18) tcon:       tube thermal conductivity       (kW/m.K)
!
!***********************************************************************

      use modsim_head
      implicit none
      integer, parameter                :: ni=8, no=6, np=18, ns=12,&
                                           ndeq=3
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real,dimension(10)                :: cof
      integer                           :: icount,i,icm,icoil,npr,nor
      logical                           :: wflo,aflo,wtem,atem,ahum,frzde

      real  :: derr,awold,ad,uad,cpbeta,tai,tsi,daw,hsi,c,y,hsp,tsp,&
               dtaodt,dtwodt,dwaodt,rate,hcd,effd,fid,usd1,cde1,trmd,&
               hcrd,betad,cd,taodb,taowb,wao,fwet,qs,tri,delt,err,fai,&
               effw,trmw,hcrw,rw,rf,rx,r,cwe1,uaw,e,cw,trmn,hse,hs1,&
               qt,hai,aw,frat,fmacp,renr,hcr1,tso,tse,hal,two,ua,fiw,&
               hcw,cfact,ftr,awstr,shr,awd,srei,fsao,roh,hae,fh,&
               fpi,fs,c1,c2,tafres,dia2,v,fv,dp,vm,hsl,tdo,tdi,fth,cm,&
               sl,fd,coid,tcon,fmr,twi,fma,taidb,wai,twodyn,taodyn,ao,&
               hd,ap,ai,flfa,fcon,fa,waodyn,a1,twb,as1

      real     :: pi    =3.141592653
      real     :: foulf =0.0003
      real     :: siacon=10.76391,silcon=3.28084,sikcon=577.789,&
                  siqcon=3412.14, simcon=7936.68,sifcon=2.54,&
                  sit1cn=1.8,     sit2cn=32.0,   sijcon=0.526543
      real     :: visc=0.0437,    cp=0.243,      h0=-7.55,&
                  slope=0.557
      real     :: rhoa1=0.075,    rhow1=62.4
      real     :: adp0=11.471,    adp1=6875.333, adp2=-280479.327,&
                  adp3=5025099.04

      icoil = par_v(1)
      ap =    par_v(2)*siacon
      as1 =    par_v(3)*siacon
      ai =    par_v(4)*siacon
      flfa =  par_v(5)
      fcon =  par_v(6)*sikcon
      fa =    par_v(7)*siacon
      npr =   par_v(9)
      nor =   par_v(10)
      tdo =   par_v(11)*silcon
      tdi =   par_v(12)*silcon
      fth =   par_v(13)*silcon
      cm =    par_v(14)*sijcon
      sl =    par_v(15)*silcon
      fd =    par_v(16)*silcon
      coid =  par_v(17)*silcon
      tcon =  par_v(18)*sikcon
  
      if(icoil==0) fd=sqrt(4.*fd*coid/(pi*nor*npr))
  
      fmr =   xin(1)*simcon
      if(fmr<1.)fmr=1.
      twi =   xin(2)*sit1cn+sit2cn
      if(twi<0.)twi=0.
      fma =   xin(3)*simcon
      if(fma<1.)fma=1.
      taidb = xin(4)*sit1cn+sit2cn
      wai =   xin(5)
      if(wai<0.)wai=0.
      twodyn= xin(6)*sit1cn+sit2cn
      taodyn= xin(7)*sit1cn+sit2cn
      waodyn= xin(8)
  
      wflo=(iostat(1)<-1).or.(iostat(1)==2)
      wtem=(iostat(2)<-1).or.(iostat(2)==2)
      aflo=(iostat(3)<-1).or.(iostat(3)==2)
      atem=(iostat(4)<-1).or.(iostat(4)==2)
      ahum=(iostat(5)<-1).or.(iostat(5)==2)
      frzde=wflo.and.aflo.and.wtem.and.atem.and.ahum
  
      ao=ap+as1
      hd=4.*fa*flfa*coid/ao
      srei=ao/ai
      fsao=as1/ao
      roh=tdo/fd
      hae=0.240*taidb+wai*(1061.+0.444*taidb)
      if(itime<=1) then

!  Initialize c1,c2, and dry fin efficiency coefficients
!    on first call of simulation.

         call sufed(roh,cof)
         do i=1,5
           saved_v(i)=cof(i)
         enddo
         fh=0.5*(fd-tdo)
         fpi = par_v(8)*sifcon
         fs=1./(12.*fpi)
         saved_v(6)=0.159*(fth/hd)**(-0.065)*(fth/fh)**0.141
         saved_v(7)=-0.323*(fs/fh)**.049*(fd/sl)**0.549*&
                   (fth/fs)**(-.028)
         saved_v(10)=slope
         saved_v(11)=h0
         saved_v(12)=0.
      endif
      do i=1,5
        cof(i)=saved_v(i)
      enddo
      c1=saved_v(6)
      c2=saved_v(7)
      if(time>saved_v(12)) then
         saved_v(8)=0.5*(saved_v(8)+saved_v(10))
         saved_v(9)=0.5*(saved_v(9)+saved_v(11))
         saved_v(12)=time
      endif

!  Begin cooling coil proper

      tafres=srei*(.5*(tdo-tdi)/tcon+foulf)
      dia2=(tdi*12.)**.2
      v=fmr*4./(npr*pi*tdi*tdi*rhow1*3600.)        ! <<<<<
      fv=fma/(rhoa1*fa)                            ! <<<<<

!  Calculate entering air dew point

      dp=adp0+wai*(adp1+wai*(adp2+wai*adp3))

!  Calculate mass velocity

      vm=(1.+wai)*fma/(flfa*fa)
      frat=fmr/fma
      fmacp=fma*cp
      renr=hd*vm/visc

!  Calc. h2o film coefficient

      hcr1=(150./dia2)*v**0.8

!  If entering air dp does not exceed twi by more than 5 F,
!  skip to dry coil section.

      icm=-1
      if(dp-twi <= 5.) then
         tso=dp+1.
         tse=dp+1.
!   (next 5 lines serve only to prevent compiler warning messages)
         hal=hae
         two=twi
         ua=0.
         fiw=0.
         hcw=0.
!   (there.)
         goto 33
      endif

!  Calculate hcw (wet surface coeff.)

      cfact=1.425-.00051*renr+.263*renr*renr*10.**(-6.)
      hcw=0.3*cfact*c1*renr**c2*vm

!  Calculate the wet fin efficiency

      ftr=twi+0.5
      awstr=-0.00793+0.00031*ftr+0.0000075*(ftr-53.)*(ftr-53.)
      shr=0.24*(taidb-ftr)/(.24*taidb+(1061.+.444*taidb)*wai-.24*ftr-&
       (1061.+0.444*ftr)*awstr)
      shr=abs(shr)
      awd=wai-awstr
      awd=abs(awd)
      fai=0.5*(fd-tdo)*sqrt(2.*hcw/(fcon*fth))
      if(taidb>80.) then
         effw=exp(-0.41718)*shr**(0.09471)*awd**(0.0108)*fai**(-0.50303)
      else
         effw=exp(-0.3574)*shr**(0.16081)*awd**(0.01995)*fai**(-0.52951)
      endif

!  Calculate outlet water temp and air enthalpy for wet coil

      fiw=fsao*effw+(1.-fsao)
      trmw=twi
      slope=saved_v(8)
      h0=saved_v(9)
 41   hcrw=hcr1*(1.+0.011*trmw)
      rw=srei/hcrw
      rf=(cp/slope)*(1.-fiw)/(fiw*hcw)
      rx=(tafres+rw)/(rf+cp/(slope*hcw))
      r=rx/slope
      cwe1=1./fma-slope/fmr
      uaw=ao/(slope*(srei/hcrw+tafres)+cp/(fiw*hcw))
      ua=uaw*slope
      e=uaw*cwe1

      if(e>20.) then
         two=(hae-h0-(slope-frat)*twi)/frat
      elseif(e<-20.) then
         two=(hae-h0)/slope
      else
         cw=exp(e)
         two=((1.-cw)*(hae-h0)+cw*(slope-frat)*twi)/(slope-cw*frat)
      endif

      trmn=(twi+two)/2.
      hal=hae-(fmr/fma)*(two-twi)
      tso=(two+r*(hae-h0))/(1.+slope*r)
      tse=(twi+r*(hal-h0))/(1.+slope*r)

      if(abs(trmn-trmw) > 0.01) then
         trmw=trmn
         hse=h0+slope*tse
         hsl=h0+slope*tso
         tse=twb(hse)
         tso=twb(hsl)
         slope=(hsl-hse)/(tso-tse)
         h0=hse-slope*tse
         saved_v(10)=slope
         saved_v(11)=h0
         goto 41
      endif

33    if(tso <= dp) then

!  Completely wet

         icm=0
         qt=fma*(hae-hal)
         hai=hae
         aw=ao

!  Completely wet coil resumes at statement c100

      else
         if(dp > tse) icm=1

!  Calculate hcd (dry surface coeff.)

         hcd=0.3*c1*renr**c2*vm

!  Calculate the dry fin efficiency

         fai=0.5*(fd-tdo)*sqrt(2.*hcd/(fth*fcon))
         effd=0.
         do i=1,5
           effd=effd+cof(i)*fai**(i-1)
         enddo
         fid=fsao*effd+(1.-fsao)
         usd1=1./(fid*hcd)
         cde1=(1./fmacp)-1./fmr
         if(icm == -1) then

!  Completely dry coil

            trmd=twi
 32         hcrd=hcr1*(1.+0.011*trmd)
            ua=ao/(usd1+srei/hcrd+tafres)
            e=cde1*ua
            if(e<-20.) then
               betad=frat/cp
            elseif(e>20.) then
               betad=1.
            else
               cd=exp(e)
               betad=(frat*(1.-cd))/(cp-frat*cd)
            endif
            taodb=taidb+betad*(twi-taidb)
            qt=fmacp*(taidb-taodb)
            two=twi+qt/fmr

!     Calculate new trm

            trmn=(twi+two)/2.
            if(abs(trmn-trmd) > 0.01) then
               trmd=trmn
               goto 32
            endif
            hal=hae-qt/fma
            taowb=twb(hal)
            taowb=taowb+0.0080*(taodb-taowb)
            wao=(hal-.24*taodb)/(1061.+0.444*taodb)
            fwet=0.
            qs=qt

!  Completely dry coil steady state calculation ends here

            goto 1000
         else

!  Partly wet coil

            tri=dp
            delt=tso-tse
!     Estimate wet area
            aw=ao*(dp-tse)/(tso-tse)
!     Estimate dry area
            slope=saved_v(8)
            h0=saved_v(9)
            err=0.
            derr=0.
            awold=aw
            icount=0
 90         ad=ao-aw
            trmd=(tri+two)/2.
            hcrd=hcr1*(1.+0.011*trmd)
            uad=ad/(usd1+srei/hcrd)
            trmw=(twi+tri)/2.
            hcrw=hcr1*(1.+0.011*trmw)
            rw=srei/hcrw
            rf=(cp/slope)*(1.-fiw)/(fiw*hcw)
            rx=(tafres+rw)/(rf+cp/(slope*hcw))
            r=rx/slope
            cwe1=1./fma-slope/fmr
            uaw=aw/(slope*(srei/hcrw+tafres)+cp/(fiw*hcw))
            ua=uad+uaw*slope
            e=cde1*uad

            if(e<-20.) then
               betad=frat/cp
            elseif(e>20.) then
               betad=1.
            else
               cd=exp(e)
               betad=(frat*(1.-cd))/(cp-frat*cd)
            endif

            cpbeta=cp*betad
            e=cwe1*uaw

            if(e>20.) then
               tri=((slope-frat)*twi-(hae-h0-cpbeta*taidb))&
                   /(cpbeta-frat)
            elseif(e<-20.) then
               tri=(hae-h0-cpbeta*taidb)/(slope-cpbeta)
            else
               cw=exp(e)
               tri=(cw*(slope-frat)*twi+(1.-cw)*(hae-h0-cpbeta*taidb))/&
                   (slope-cw*frat-(1.-cw)*cpbeta)
            endif

            tai=taidb+betad*(tri-taidb)
            hai=hae-cp*(taidb-tai)
            tsi=(tri+r*(hai-h0))/(1.+slope*r)
            hal=hai-frat*(tri-twi)
            tse=(twi+r*(hal-h0))/(1.+slope*r)
            qt=fma*(hae-hal)
            two=twi+qt/fmr
            if(icount>=1) derr=dp-tsi-err
            err=dp-tsi

            if(abs(err) > 0.01) then
               if(icount<1 .or. abs(derr)<1.e-6) then
                  awold=aw
                  aw=aw+0.7*ao*err/(tso-tse)
               else
                  daw=aw-awold
                  awold=aw
                  aw=aw-err*daw/derr
               endif
               icount=icount+1
               if(aw>ao)aw=ao
               if(aw<=0.)aw=0.001*ao
               hse=h0+slope*tse
               hsi=h0+slope*tsi
               tse=twb(hse)
               tsi=twb(hsi)
               slope=(hsi-hse)/(tsi-tse)
               h0=hse-slope*tse
               saved_v(10)=slope
               saved_v(11)=h0
               goto 90
            endif
         endif
      endif
!  Partly or completely wet coil

      taowb=twb(hal)

!    Calculate dry bulb temperature

      c=(hcw*aw)/(cp*fma)
      c=c*fiw
      y=0.
      if(c<20.)y=exp(-c)
      hsp=hai-(hai-hal)/(1.-y)
      tsp=twb(hsp)
      taodb=tsp+(taidb-tsp)*y
      wao=(hal-.24*taodb)/(1061.+0.444*taodb)
      qs=0.243*fma*(taidb-taodb)
      fwet=aw/ao
 1000 continue

!  End of steady state model; tack on simple dynamics
!  steady state output if par_v(14) < 1.0

      if(cm<=sijcon) then
         dtaodt=taodb-sit2cn
         dtwodt=two-sit2cn
         dwaodt=wao
         iostat(1)=1
         iostat(2)=1
         iostat(3)=1
      else
         rate=ua/(3600.*cm)
         dtaodt=rate*(taodb-taodyn)
         dwaodt=rate*(wao-waodyn)
         dtwodt=rate*(two-twodyn)
         iostat(1)=0
         iostat(2)=0
         iostat(3)=0
         if((abs(two-twodyn)<=rtolx*abs(two)+atolx) .and.&
           (abs(taodb-taodyn)<=rtolx*abs(taodb)+atolx) .and.&
           (abs(wao-waodyn)<=rtolx*abs(wao)+atolx).and.frzde) then

            iostat(1)=1
            iostat(2)=1
            iostat(3)=1
         endif
      endif

!  Convert units for output (all coil conditions)

      yout(1)=dtwodt/sit1cn
      yout(2)=dtaodt/sit1cn
      yout(3)=dwaodt
      yout(4)=qt/siqcon
      yout(5)=qs/siqcon
      yout(6)=fwet

      iostat(4)=1
      iostat(5)=1
      iostat(6)=1

      return
      end subroutine type12
! ----------------------------------------------------------------------
!
      real function twb(h)
!
!   The following function defines wet bulb temperature
!   as a function of log(enthalpy):
!
      if(h<11.83)then
         twb=32.
      else
         alogh=log(h)
         twb=30.9185+alogh*(-39.682+alogh*(20.5841-1.758*alogh))
      endif
      return
      end function twb

!***********************************************************************

      subroutine type13(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 13 : THREE-WAY VALVE MODEL
!
!     c:   control signal from controller
!     cv:  valve position before hysteresis is considered
!     y1:  actual valve position
!
!     Port 1 is closed when C is zero.
!
!     Head loss (k) parameters are related to valve constant (av) by:
!       k=1/(rhow*av)**2
!
!     isel ([Par(6)] determines characteristics of Port 1 / Port 2 :
!         isel=0 -> Linear/Linear
!         isel=1 -> Linear/Exponential
!         isel=2 -> Exponential/Exponential
!
!***********************************************************************
  
      use modsim_head
      implicit none
      integer, parameter                :: ni=7, no=6, np=6, ns=5,&
                                           ndeq=1
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: k1,k2
      logical        :: con
      integer        :: isel
      real           :: hys,dcvdt,dc,y1,y2,w3,s1,s2,p1,p2,t3,w1,w2,p3,&
                        t1,t2,c,cv,cl,tcon,hyster,a1

      w1=  xin(1)
      w2=  xin(2)
      p3=  xin(3)
      t1=  xin(4)
      t2=  xin(5)
      c=   xin(6)
      cv=  xin(7)
  
      k1=  par_v(1)
      k2=  par_v(2)
      cl=  par_v(3)
      tcon=par_v(4)
      hys= par_v(5)
      isel=int(par_v(6)+0.1)
  
      con=(iostat(6)<-1).or.(iostat(6)==2)

      if(c<0.0) then
        c=0.0
      elseif(c>1.0) then
        c=1.0
      endif

      if(cv<0.0) then
        cv=0.0
      elseif(cv>1.0) then
        cv=1.0
      endif

      if(tcon>=1.0) then
        dcvdt=(c-cv)/tcon
      else
        if(itime<=1) then
          if(init==0) then                           
            saved_v(4)=cv                            
          endif                                      
          if(init==0 .or. saved_v(3)>time) then      
            saved_v(3)=0.0                           
          endif                                      
        endif                                        
                                                     
        if(time>saved_v(3)) then                     
          saved_v(5)=saved_v(4)                      
        endif                                        
                                                     
        dcvdt=c                                      
        cv=c                                         
        a1=tcon/tstep

        if(a1>=0.05) then
          dc=c-saved_v(5)
          if(abs(dc)>=1.0e-10) then
            cv=c-dc*exp(-1.0/a1)
            dcvdt=cv
          endif
        endif
      endif

2000  y1=hyster(cv,saved_v(1),saved_v(2),saved_v(3),hys)
      y2=1.0-y1
      w3=w1+w2

      if(w1/=0.0) then
        s1=sign(1.0,w1)
      else
        s1=1.0
      endif

      if(w2/=0.0) then
        s2=sign(1.0,w2)
      else
        s2=1.0
      endif

      select case(isel)
        case(0)
          p1=p3+k1*s1*(w1/((1.0-cl)*y1+cl))**2
          p2=p3+k2*s2*(w2/((1.0-cl)*y2+cl))**2
        case(1)
          p1=p3+k1*s1*(w1/(cl**(1.0-y1)))**2
          p2=p3+k2*s2*(w2/((1.0-cl)*y2+cl))**2
        case default
          p1=p3+k1*s1*(w1/(cl**(1.0-y1)))**2
          p2=p3+k2*s2*(w2/(cl**(1.0-y2)))**2
      end select

      if(w3>0.0) then
        t3=((w1*t1)+(w2*t2))/w3
      else
        t3=0.5*(t1+t2)
      endif

      saved_v(4)=cv

      yout(1)= dcvdt
      yout(2)= w3
      yout(3)= p1
      yout(4)= p2
      yout(5)= t3
      yout(6)= y1

      if((abs(c-cv)<=rtolx*abs(cv)+atolx).and.con) then
        iostat(1)=1
      else
        iostat(1)=0
      endif

      iostat(2)=1
      iostat(3)=1
      iostat(4)=1
      iostat(5)=1
      iostat(6)=1
  
      return
      end subroutine type13
!***********************************************************************
  
      subroutine type14(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 14 : EVAPORATIVE HUMIDIFIER MODEL
!
!       This model is developed  based on the model by J. Chi.
!
!***********************************************************************
  
      use modsim_head
      implicit none

      integer, parameter                :: ni=3, no=3, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real           :: tai,wai,flwai,ha,hai,twtb,waos,wao,tao,flwao,psats
      real           :: cps=1.8723, hfg=2501.0, cpa=1.0
      real           :: rmol=0.62198, tref=-20.0, patm=101.325
  
      tai=  xin(1)
      wai=  xin(2)
      flwai=xin(3)
  
      ha=par_v(1)
  
      hai=(cps*(tai-tref)+hfg)*wai/(1.0+wai)+cpa*(tai-tref)
  
!     Convert inlet enthalpy to btu/lbm for the purpose of
!     calculating the wet bulb temperature.
  
      hai=hai/2.326
      twtb=30.9185+log(hai)*(-39.682+log(hai)*(20.5841-&
       1.758*log(hai)))
      twtb=(twtb-32.0)*5.0/9.0
      waos=rmol*psats(twtb)/(patm-psats(twtb))
      wao=wai+(waos-wai)*(1.0-exp(-ha/(cpa*flwai)))
      tao=tai-(wao-wai)*hfg/cpa
      flwao=flwai*(1.0+wao)/(1.0+wai)
  
      yout(1)=wao
      yout(2)=tao
      yout(3)=flwao
  
      iostat(1)=1
      iostat(2)=1
      iostat(3)=1
  
      return
      end subroutine type14
  
!***********************************************************************
  
      subroutine type15(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!       TYPE 15 : ROOM MODEL
!
!***********************************************************************

      use modsim_head
      implicit none
      integer, parameter                :: ni=8, no=6, np=6, ns=9,&
                                           ndeq=4
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      integer                           :: i
      real           :: to,tbar,tau2,tau1,dt,e,t2,e1a,tr,dtm1dt,dtm2dt,&
                        dt2bdt,q1,q2,vol,wcap,rcap,h1,h2,f,f1,acap,dt1dt,&
                        ts,wa,ti,t1,tm1,tm2,t2b,delay,a1
      real           :: rhoa=1.2,rhow=1000.0,cpa=1.0,cpw=4.180

      wa=  xin(1)
      ti=  xin(2)
      t1=  xin(3)
      tm1= xin(4)
      tm2= xin(5)
      t2b= xin(6)
      q1=  xin(7)
      q2=  xin(8)
  
      vol= par_v(1)
      wcap=par_v(2)
      rcap=par_v(3)
      h1=  par_v(4)
      h2=  par_v(5)
      f=   par_v(6)
  
      if(saved_v(9)>time) saved_v(9)=0.
      f1=1.-f
  
      do i=1,6
        iostat(i)=0
      enddo
      acap=vol*rhoa*cpa
      dt1dt=(q2/f+h1*(tm1-t1)+h2*(tm2-t1))/acap
      ts=(h1*tm1+h2*tm2)/(h1+h2)
      if(wa>1.e-6) then
1000    tau1=vol*rhoa/wa                                                   
        tau2=f1*tau1                                                       
        tau1=f*tau1                                                        
        dt=ts-t1
        if(abs(dt)<1.e-10) dt=0.

        a1=f1*(h1+h2)/(wa*cpa)                                             

        if(a1<20.) then
          e=exp(-a1)
        else
          e=0.
        endif

        t2=ts-dt*e                                                         
        to=delay(t2,tau2,saved_v(1),saved_v(7),saved_v(8),saved_v(9))      
        dt1dt=dt1dt+(ti-t1)/tau1                                           

        if(e<1.0) then
          e1a=(1.0-e)/a1
        else
          e1a=1.
        endif

        tbar=ts-dt*e1a                                                     
      else
        to=ts
        tbar=ts
        tau2=f1*(h1+h2)/acap
      endif

2000  tr=f*t1+f1*t2b
      dtm1dt=(h1*(tr-tm1)+q1)/wcap
      dtm2dt=h2*(tr-tm2)/rcap
      dt2bdt=(tbar-t2b)/tau2
  
      yout(1)= dt1dt
      yout(2)= dtm1dt
      yout(3)= dtm2dt
      yout(4)= dt2bdt
      yout(5)= tr
      yout(6)= to
  
      return
      end subroutine type15

!***********************************************************************
  
      subroutine type16(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 16 : STICKY PROPORTIONAL CONTROLLER
!
!       NOTE: this component it to be used ONLY to pass control signals
!       between   superblocks. The variable defined by its output
!       should not be solved simultaneously!
!
!***********************************************************************

      use modsim_head
      implicit none

      integer, parameter                :: ni=2, no=1, np=3, ns=3
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: gainp,band,e,enew,de,outt,ctm,cts,offset

      ctm=  xin(1)
      cts=  xin(2)
  
      gainp= par_v(1)
      band=  par_v(2)
      offset=par_v(3)
  
      e=cts-ctm

      if(itime<=1) then
        if(init==0) then
          saved_v(2)=e
        endif

        if(init==0 .or. saved_v(1)>time) then
          saved_v(1)=0.0
        endif
      endif

      if(time>saved_v(1)) then
        saved_v(3)=saved_v(2)
      endif

      saved_v(1)=time
      enew=saved_v(3)
      de=e-saved_v(3)

      if(abs(de)>=band) then
        if(abs(de)<=2.0*band) then
          enew=e
        else
          enew=saved_v(3)+sign(band,de)
        endif
      endif

      saved_v(2)=enew
      outt=offset+gainp*enew

      if(outt>1.0) then
        outt=1.0
      elseif(outt<0.0) then
        outt=0.0
      endif
  
      yout(1)=outt
      iostat(1)=0
  
      return
      end subroutine type16
  
!***********************************************************************
  
      subroutine type17(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 17 : MIXING DAMPERS and MERGE MODEL
!
!***********************************************************************
  
      implicit none
      integer, parameter                :: ni=6, no=4, np=3, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      integer        :: i
      real           :: k
      real           :: w1,w2,w3,t3,t2,c,cl,wf,dp,s,r1,r2,p1,p2,p3,t1

      p1= xin(1)
      p2= xin(2)
      p3= xin(3)
      t1= xin(4)
      t2= xin(5)
      c=  xin(6)
  
      k=  par_v(1)
      cl= par_v(2)
      wf= par_v(3)
  
      if(c<0.0) then
        c=0.0
      elseif(c>1.0) then
        c=1.0
      endif

      dp=p1-p3

      if(dp/=0.0) then
        s=sign(1.0,dp)
      else
        s=1.0
      endif

      r1=wf*k/((1.-cl)*c+cl)**2 + (1.0-wf)*k*cl**(2.0*c-2.0)
      r2=wf*k/((1.-cl)*(1.-c)+cl)**2 + (1.0-wf)*k*cl**(-2.0*c)
      w1=s*sqrt(abs(dp)/r1)
      dp=p2-p3

      if(dp/=0.0) then
        s=sign(1.0,dp)
      else
        s=1.0
      endif

      w2=s*sqrt(abs(dp)/r2)
      w3=w1+w2
      t3=0.5*(t1+t2)

      if(w3>0.0) then
        t3=(w1*t1+w2*t2)/w3
      endif
  
      yout(1)= w1
      yout(2)= w2
      yout(3)= w3
      yout(4)= t3
  
      do i=1,4
         iostat(i)=1
      end do
  
      return
      end subroutine type17
  
!***********************************************************************
  
      subroutine type18(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 18 : PLENUM MODEL
!
!***********************************************************************
  
      implicit none
      integer, parameter                :: ni=10, no=1, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      integer        :: n,i
      real           :: w

      w=0.0
      n=ifix(par_v(1))
      do i=1,n
         w=w+xin(i)
      end do
  
      yout(1)=w
      iostat(1)=1
  
      return
      end subroutine type18
  
!***********************************************************************
  
      subroutine type19(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 19 : FLOW BALANCE CONTROL
!
!       Used to equalize flows calculated independently in different
!       superblocks.
!
!***********************************************************************
  
      use modsim_head
      implicit none

      integer, parameter                :: ni=3, no=1, np=1, ns=3
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: k
      real           :: w,p2,wp,dk,s,p1

      w=  xin(1)
      p2= xin(2)
      wp= xin(3)
  
      if(itime<=1) then  !  goto 1000
        if(init==0 .or. saved_v(1)>time) then
          saved_v(1)=0.0
        endif
        if(init==0) then
          saved_v(2)=par_v(1)
        endif
      endif

1000  if(time>saved_v(1)) then
        saved_v(3)=saved_v(2)
      endif

      saved_v(1)=time

      if(w>0.0) then
        dk=2.0*saved_v(3)*(1.0-wp/w)
      else
        dk=0.0
      endif

      k=saved_v(3)+dk

      if(w/=0.0) then
        s=sign(1.0,w)
      else
        s=1.0
      endif

      p1=p2+s*k*w**2
      saved_v(2)=k
  
      iostat(1)=1
      yout(1)=p1
  
      return
      end subroutine type19
  
!***********************************************************************
  
      subroutine type20(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 20 : HIGH/LOW LIMIT CONTROLLER
!
!***********************************************************************
  
      implicit none
      integer, parameter                :: ni=2, no=1, np=2, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: ctm,c,gain,ctl,e

      ctm= xin(1)
      c=   xin(2)
  
      gain= par_v(1)
      ctl=abs(par_v(2))
  
      e=ctl-ctm
      if(par_v(2)>=0.0) then  
        if(ctm>ctl) then
          c=c+gain*e
        endif
      else
        if(ctm<ctl) then
          c=c+gain*e
        endif
      endif

      if(c>1.0) then
        c=1.0
      elseif(c<0.0) then
        c=0.0
      endif
  
      yout(1)=c
      iostat(1)=1
  
      return
      end subroutine type20
  
!***********************************************************************
  
      subroutine type21(xin,yout,par_v,saved_v,iostat)
  
! ---------------------------------------------------------------------
!
!     TYPE 21 :  CLAMPED SPLIT
!
!***********************************************************************
  
      implicit none
      integer, parameter                :: ni=2, no=4, np=2, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      integer        :: i
      real           :: k
      real           :: p2,p3,p0,s1,s2,s3,dp2,dp3,w1,w2,w3,p1
  
      p2=  xin(1)
      p3=  xin(2)
  
      p0=  par_v(1)
      k=   par_v(2)
  
      dp2=p0-p2
      dp3=p0-p3

      if(dp2/=0.0) then
        s2=sign(1.0,dp2)
      else
        s2=1.0
      endif

      if(dp3/=0.0) then
        s3=sign(1.0,dp3)
      else
        s3=1.0
      endif

      w2=s2*sqrt(2.0*abs(dp2)/k)
      w3=s3*sqrt(2.0*abs(dp3)/k)
      w1=w2+w3

      if(w1/=0.0) then
        s1=sign(1.0,w1)
      else
        s1=1.0
      endif

      p1=p0+0.5*k*s1*w1**2
  
      yout(1)= p1
      yout(2)= w1
      yout(3)= w2
      yout(4)= w3
  
      do i=1,4
         iostat(i)=1
      end do
  
      return
      end subroutine type21
  
!***********************************************************************
  
      subroutine type22(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 22 : STEAM SPRAY HUMIDIFIER
!                 based on mathematical model by J.Chi
!
!***********************************************************************
  
      use modsim_head
      implicit none

      integer, parameter                :: ni=5, no=3, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      integer                           :: i
      real           :: ts,tai,wai,flwai,flws,eta,cpai,tao,wao,&
                        waos,psats,flwao
      real           :: cps=1.8723, rmol=0.62198, patm=101.325, cpa=1.0

      ts=   xin(1)
      tai=  xin(2)
      wai=  xin(3)
      flwai=xin(4)
      flws= xin(5)
  
      eta=par_v(1)
  
      cpai=(cpa+wai*cps)/(1.0+wai)
      tao=(flws*cps*ts+flwai*cpai*tai)/(flws*cps+flwai*cpai)
      wao=wai+flws*(1.0+wai)/flwai
      waos=rmol*psats(tao)/(patm-psats(tao))

      if(wao>eta*waos) then
        wao=eta*waos
      endif

      flwao=flwai*(1.0+wao)/(1.0+wai)
      yout(1)=tao
      yout(2)=flwao
      yout(3)=wao
  
      do i=1,3
         iostat(i)=1
      end do
  
      return
      end subroutine type22
  
!***********************************************************************
  
      subroutine type23(xin,yout,par_v,saved_v,iostat)
  
! ---------------------------------------------------------------------
!
!     TYPE 23 : ISENTROPIC STEAM NOZZLE
!                  using real gas properties
!
!***********************************************************************
  
      implicit none
      integer, parameter                :: ni=3, no=2, np=2, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: pzero,tin,pin,pout,tsin,tsats,s,ssats,hin,&
                        hsats,ss,hs,tout,tpss,tsat,sg,hout,vout,vs,&
                        sf,ssatw,hf,hsatw,vf,vsatw,vg,vsats,hg,x,w,a1

      a1=   par_v(1)
      pzero=par_v(2)
  
      tin= xin(1)
      pin= xin(2)+pzero
      pout=xin(3)+pzero
  
      tsin=tsats(pin)

      if(tin<=tsin) then

!     Entering steam is saturated

         s=ssats(tsin)
         hin=hsats(tsin)
      else

!     Entering steam is superheated

         s=ss(pin,tin)
         hin=hs(pin,tin)
      endif

      tout=tpss(pout,s)
      tsat=tsats(pout)
      sg=ssats(tsat)

      if(s>sg)then
         hout=hs(pout,tout)
         vout=vs(pout,tout)
      else
         sf=ssatw(tsat)
         hf=hsatw(tsat)
         vf=vsatw(tsat)
         vg=vsats(pout,tsat)
         hg=hsats(tsat)
         x=(s-sf)/(sg-sf)
         vout=vf+x*(vg-vf)
         hout=hf+x*(hg-hf)
      endif

      w=a1/vout*sqrt(2.0*1000.0*(hin-hout))
  
      yout(1)= tout
      yout(2)= w
  
      iostat(1)=1
      iostat(2)=1
  
      return
      end subroutine type23
  
!***********************************************************************
  
      subroutine type24(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 24 : IDEAL GAS FLOW NOZZLE
!
!***********************************************************************
  
      implicit none
      integer, parameter                :: ni=3, no=2, np=4, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: pu,pd,t,ac,po,gama,r,pcrit,ratio,c3,w,texit

      pu=  xin(1)
      pd=  xin(2)
      t=   xin(3)
  
      ac=  par_v(1)
      po=  par_v(2)
      gama=par_v(3)
      r=   par_v(4)
  
      t=t+273.15
      pu=pu+po
      pd=pd+po
      pcrit=(2.0/(1.0+gama))**(gama/(gama-1.0))
      ratio=pd/pu
      c3=(ratio**(2.0/gama)-ratio**((1.0+gama)/gama))*2.0/(gama-1.0)
      if(ratio<pcrit)c3=(2.0/(1.0+gama))**((1.0+gama)/(gama-1.0))
      w=1000.0*ac*pu*sqrt(gama*c3/(r*t))
      texit=t*ratio**((gama-1.0)/gama)-273.15
  
      yout(1)= w
      yout(2)= texit
  
      iostat(1)=1
      iostat(2)=1
  
      return
      end subroutine type24
  
!***********************************************************************
  
      subroutine type25(xin,yout,par_v,saved_v,iostat)
  
! ---------------------------------------------------------------------
!
!     TYPE 25 : STEAM TO AIR HEATING COIL
!
!  Parameters:
!
!  1) tsc: degrees of subcooling of condensate (C)
!  2) aa:  minimum air flow area (m2)
!  3) ais: inside heat transfer area (m2)
!  4) aop: primary (tube) outside area (m2)
!  5) aos: secondary (fin) heat transfer area (m2)
!  6) di:  inside tube diameter (m)
!  7) do:  outside tube diameter (m)
!  8) akt: tube thermal conductivity [kJ/(kg K)]
!  9) df:  fin diameter (m)
! 10) fcon: fin thermal conductivity [kJ/(kg K)]
! 11) fth: fin thickness (m)
! 12) fpm: number of fins per meter
! 13) cm:  thermal capacitance of fins and tubes (kJ/K)
!
!***********************************************************************
  
      use modsim_head
      implicit none

      integer, parameter                :: ni=5, no=4, np=13, ns=10,&
                                           ndeq=1
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real                              :: ka,kw,mua,muw,nusslt,ntu
      real                              :: grav=9.80665
      real           :: tai,tao,wa,tsc,aa,ais,aop,aos,di,do,akt,df,&
                        fcon,fth,fpm,cm,tair,tssat,roh,rt,x,cpair,re,&
                        pr,ffo,htcos,fai,fineff,ro,dh,hsl,rr,ts,htcis,&
                        ri,rsum,rrr,q,ws,taoss,two,z,tau,dtaodt,pso,&
                        tsi,tsats,aka,visca,hfg,hsatw,wmu,wk,hs,a1
      real           :: rhoa=1.2,rhow=1000.0,cpa=1.0,cpw=4.180

      pso=  xin(1)
      tsi=  xin(2)
      tai=  xin(3)
      tao=  xin(4)
      wa=   xin(5)
  
      tsc=  par_v(1)
      aa=   par_v(2)
      ais=  par_v(3)
      aop=  par_v(4)
      aos=  par_v(5)
      di=   par_v(6)
      do=   par_v(7)
      akt=  par_v(8)
      df=   par_v(9)
      fcon= par_v(10)
      fth=  par_v(11)
      fpm=  par_v(12)
      cm=   par_v(13)
  
      tair=0.5*(tai+tao)
      tssat=tsats(pso)
      roh=do/df
      rt=(do-di)/(2.0*akt*ais)
  
!     Determine outside heat transfer coefficient for air, using
!     correlation from H.Hausen, "Heat transfer in counterflow,
!     parallel flow and cross flow," McGraw-Hill, 1983.
  
      ka=aka(tair)
      mua=visca(tair)
      call cpcva(tair,cpair,x,x,x)
      re=wa*do/(mua*aa)
      pr=mua*cpair/ka

!     ffo is the ratio of finned tube surface area to the surface area of
!     a similar but unfinned tube.

      ffo=(0.5*(df*df-do*do)+fth*(df-do)+do*(1.0/fpm-fth))*fpm/do
      nusslt=0.3*re**0.625*ffo**(-0.375)*pr**(1.0/3.0)
      htcos=nusslt*ka/do

!     Determine dry fin efficiency using coefficients on first call

      if(itime<=1) then
        call sufed(roh,saved_v)
      endif

      fai=0.5*(df-do)*sqrt(2.0*htcos/(fth*fcon))
      fineff=saved_v(1)+fai*(saved_v(2)+fai*(saved_v(3)+fai*(saved_v(4)&
       +fai*saved_v(5))))
      ro=1./(htcos*(aos*fineff+aop))

!     If inlet steam temp. is below saturation, ignore input pressure
!     and assume inlet steam is sat. vapor at inlet temp.

      if(tsi<=tssat) then
        tssat=tsi
      endif

      dh=hfg(tssat)
      hsl=hsatw(tssat)
      muw=wmu(tssat)
      kw=wk(tssat)
!**     rhow=wrho(tssat)

!     Iterate for inside surface temperature and inside film coefficient

      rr=0.1
  50  ts=tssat-rr*(tssat-tair)
      htcis=0.612*(kw**3*rhow*rhow*dh*grav/(muw*di*(tssat-ts)))**0.25
      ri=1.0/(htcis*ais)
      rsum=ri+rt+ro
      rrr=ri/rsum
      if(abs(rrr-rr)/rrr>rtolx)then
         rr=rrr
         goto 50
      endif
      rr=(ro+rt)/ri
      ntu=1.0/(rsum*wa*cpa)

!     Find mass flow rate of steam and condensate

      q=htcis*ais*(tssat-ts)
      ws=q/(hs(pso,tsi)-(hsl-cpw*tsc))
      taoss=tai+q/(wa*cpair)
      two=tssat-tsc

!     Air temperature dynamics, based on Myers et al., Trans. ASME, J. Heat
!     Transfer, May 1970, pp. 269-275.

      z=0.5*ntu*(1.+rr)/rr
      a1=2.0*ntu*(1.+rr)*exp(-z)*sinh(z)/((sinh(z)+cosh(z))*exp(-z)&
       -exp(-ntu))
      tau=cm/(a1*wa*cpair)
      dtaodt=(taoss-tao)/tau
      q=wa*cpair*(tao-tai)

      iostat(1)=0
      if((abs(taoss-tao)<rtolx*taoss+atolx) .and.&
        (iostat(1) < -1 .or. iostat(1) == 2) .and.&
        (iostat(2) < -1 .or. iostat(2) == 2) .and.&
        (iostat(3) < -1 .or. iostat(3) == 2)) then

        iostat(1)=1
        iostat(2)=1  
        iostat(3)=1  
        iostat(4)=1  
      end if

      yout(1)=dtaodt
      yout(2)=two
      yout(3)=ws
      yout(4)=q
  
      return
      end subroutine type25
  
!***********************************************************************
  
      subroutine type26(xin,yout,par_v,saved_v,iostat)
  
! ---------------------------------------------------------------------
!
!     TYPE 26 : CONTROL SIGNAL INVERTER
!                 Cout = 1.0 - Cin ; 0.0 < Cout < 1.0
!       Revised : Sept. 24, 1986 C.P.
!
!***********************************************************************
  
      implicit none
      integer, parameter                :: ni=1, no=1, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      yout(1)=par_v(1)*min(max(1.0-xin(1),0.0),1.0)
      iostat(1)=1
  
      return
      end subroutine type26
  
!***********************************************************************
  
      subroutine type27(xin,yout,par_v,saved_v,iostat)
!
! ----------------------------------------------------------------------
!
!     TYPE 27 : MOIST AIR FLOW MERGE MODEL
! 
!***********************************************************************
  
      implicit none
      integer, parameter                :: ni=7, no=5, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      integer        :: i
      real           :: k
      real           :: w1,w2,p3,t1,t2,ahum1,ahum2,w3,p1,p2,t3,ahum3

      w1=   xin(1)
      w2=   xin(2)
      p3=   xin(3)
      t1=   xin(4)
      t2=   xin(5)
      ahum1=xin(6)
      ahum2=xin(7)
  
      k=    par_v(1)
  
      w3=w1+w2
      p1=p3+0.5*k*(w1*abs(w1)+w3*abs(w3))
      p2=p3+0.5*k*(w2*abs(w2)+w3*abs(w3))
      if(abs(w3)<1.0e-12) then
         t3=0.5*(t1+t2)
         ahum3=0.5*(ahum1+ahum2)
      else
         t3=(w1*t1+w2*t2)/w3
         ahum3=(w1*ahum1+w2*ahum2)/w3
      endif
  
      yout(1)= w3
      yout(2)= p1
      yout(3)= p2
      yout(4)= t3
      yout(5)= ahum3
  
      do i=1,5
         iostat(i)=1
      end do
  
      return
      end subroutine type27
  
!***********************************************************************
  
      subroutine type28(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 28 : CONSTANT FLOW RESISTANCE MODEL
!
!***********************************************************************
  
      implicit none
      integer, parameter                :: ni=2, no=1, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: k
      real           :: w,p2,p1

      w=   xin(1)
      p2=  xin(2)
  
      k=   par_v(1)
  
      p1=p2+k*w*abs(w)
  
      yout(1)=p1
      iostat(1)=1
  
      return
      end subroutine type28
  
!************************************************************************
  
      subroutine type29(xin,yout,par_v,saved_v,iostat)
  
! ----------------------------------------------------------------------
!
!     TYPE 29 : INLET CONSATANT FLOW RESISTANCE MODEL
!
!***********************************************************************
  
      implicit none
      integer, parameter                :: ni=2, no=1, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real           :: k
      real           :: p1,p2,dpk,w

      p1=  xin(1)
      p2=  xin(2)
  
      k=   par_v(1)
  
      dpk=(p1-p2)/k

      if(dpk<abs(dpk)) then
        w=-w
      else
        w=sqrt(abs(dpk))
      endif
  
      yout(1)=   w
      iostat(1)=1
  
      return
      end subroutine type29
