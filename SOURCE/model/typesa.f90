!***********************************************************************
!
      subroutine type30(xin,yout,par_v,saved_v,iostat)
!
! ----------------------------------------------------------------------
!
!     TYPE 30 :   COOLING COIL DUTY FROM INLET CONDITIONS
!
!     J. S. Nizet, University of Liege, Belgium
!
!     Modified at National Bureau of Standards
!
!     Revised : July 7, 1987 Cheol Park
!     Upadted to Fortran 90  February 9, 2007 Cheol Park, NIST
!
!
!   input===>  1.inlet air mass flow rate (kg/sec)
!              2.inlet air dry temperature (C)
!              3.inlet air humidity ratio (kg/kg)
!              4.inlet water temperature (C)
!              5.water mass flow rate (kg/sec)
!              6.switch for cooling coil valve (-)
!
!   output===> 1.outlet air dry temperature (C)
!              2.outlet air relative humidity (%)
!              3.outlet air humidity ratio (kg/kg)
!              4.outlet air entalpy (kJ/kg)
!              5.outlet water temperature on coil side (C)
!              6.total heat transfer = coil duty (kW)
!              7.sensible heat transfer (kW)
!              8.removal water from air flow (kg/h)
!              9.sensible heat ratio
!             10.coil effectiveness
!             11.coil bypass factor
!             12.equivalent au value (kW/K)
!
!   param.===> 1.atmospheric pressure (Pa)
!              2.free flow area (m**2)
!              3.air side area (m**2)
!              4.water side area (m**2)
!              5.hydraulic diameter (m)
!              6.number of water circuits
!              7.fin thickness (m)
!              8.fin hight (m)
!              9.fin area (m**2)
!             10.tube inside diameter (m)
!             11.tube outside diameter (m)
!             12.colburn constant c1
!             13.colburn constant c2
!             14.amlow  : mimimum air flow rate (kg/s)
!
!***********************************************************************

      use modsim_head
      implicit none

      integer, parameter                :: ni=6, no=12, np=14, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real      :: rifin,rofin,cflow,cpi,hi,viscmu,gair,&
                   reynol,ha,x1,x2,x3,efin,efingl,ra,rm,&
                   uwater,hw,rw,cair,cwat,xnum1,bf,crr,&
                   rdry,antu,xnum2,eo,effect,qt,shr,ho,&
                   hs,ts,tao,r,xnum3,qtnew,two,qsensi,go,&
                   at,pws,rho,au,rewat,tai,gi,twi,wm,swccv,&
                   patm,afree,aext,aint,hydiam,ncirc,thifin,&
                   heifin,afin,diaint,diaext,cair1,cair2,&
                   amlow,pi,aratio,bratio,fhair,calar,visca,&
                   ftsat,fwha,am1

      real      :: c1=-5674.5359,    c2=6.3925247,&
                   c3=-0.9677843e-2, c4=0.62215701e-6,&
                   c5=0.20747825e-8, c6=0.9484024e-12,&
                   c7=4.1635019,     c8=-5800.2206,&
                   c9=1.3914993,     c10=-0.04860239,&
                   c11=0.41764768e-4,c12=-0.14452093e-7,&
                   c13=6.5459673

      real      :: rhowat=1000.0,cpg=1.805,cpa=1.0,cpw=4.18,&
                   prandl=0.71,confin=0.204

      real      :: eps=1.0e-3
      integer   :: i,niter=50,ifile=1,k
      ! MAG change k from real to int

!     Inputs, parameters and geometrical characteristics

      am1=   xin(1)
      tai=   xin(2)
      gi=    xin(3)
      twi=   xin(4)
      wm=    xin(5)
      swccv= xin(6)

      patm=  par_v(1)
      afree= par_v(2)
      aext=  par_v(3)
      aint=  par_v(4)
      hydiam=par_v(5)
      ncirc= par_v(6)
      thifin=par_v(7)
      heifin=par_v(8)
      afin=  par_v(9)
      diaint=par_v(10)
      diaext=par_v(11)
      cair1= par_v(12)
      cair2= par_v(13)
      amlow= par_v(14)

      pi=4.0*atan(1.0)
      aratio=afin/aext
      bratio=aext/aint
      rifin=diaext/2.0
      rofin=rifin+heifin
      cflow=ncirc*(pi*diaint**2/4)

!     Calculation of air inlet entalpy

      cpi=cpa+cpg*gi
      hi=fhair(tai,gi)

!     If swccv is in the range of [0, 0.5], then bypass the cooling
!     coil.

      if(swccv>=0.0 .and. swccv<=0.5) then
        goto 1000
      endif

!     Is the cooling required ?

!     If the air flow rate is lower than the minimum air flow rate,
!     no cooling coiled is needed.

      if((am1<amlow) .or. (wm<eps)) goto 1000

!     Air side thermal resistance :  colburn=cair1*reynolds**cair2
!     ra in k*m**2/kw

!      viscmu=9.90e-8*(tref+tai)-9.87e-6

      viscmu=visca(tai)

      gair=am1/afree
      reynol=gair*hydiam/viscmu
      ha=gair*cpi*cair1*reynol**cair2*prandl**(-2.0/3.0)

!     Fin effectiveness

      x1=sqrt(2.0*ha/confin/thifin)
      x2=1.0+0.35*alog(rofin/rifin)
      x3=x1*x2*heifin
      efin=tanh(x3)/x3
      efingl=1.0-aratio*(1.0-efin)
      ra=1.0/(efingl*ha)
      rm=0.0

!     Water side thermal resistence : nu=0.023*re**0.8**pr**0.33
!     rw in k*m**2/kw

      uwater=wm/(rhowat*cflow)
      hw=1.429*(1.0+0.0146*twi)*diaint**(-0.2)*uwater**0.8
      rw=bratio/hw

!     bf:     bypass factor
!     crr:    thermal capacity rate ratio (dry)
!     antu:   number of transfer units (dry)
!     effect: coil effectiveness (dry)

      cair=cpi*am1
      cwat=cpw*wm
      xnum1=aext/(cair*ra)

      if(xnum1>100.0) then
        bf=0.0
      else
        bf=exp(-xnum1)
      endif

      crr=cair/cwat

!     Initialization for coil duty (dry coil assumption)

      rdry=ra+rm+rw
      antu=aext/(cair*rdry)
      xnum2=antu*(1.0-crr)

      if(xnum2>100.0) then
        eo=0.0
      else
        eo=exp(-xnum2)
      endif

      effect=(1.0-eo)/(1.0-crr*eo)
      qt=effect*cair*(tai-twi)

      shr=1.0

!     Iteration for coil duty

      do k=1,niter

!       Outlet air enthalpy
                                                                                  
          ho=hi-effect*cpi*(tai-twi)/shr                                          
!         ho=hi-qt/am1
                                                                                  
!       Entalpy of saturated air and corresponding                                
!       effective surface temperature                                             
                                                                                  
        hs=(ho-bf*hi)/(1-bf)                                                      
        ts=ftsat(hs,patm)                                                         
                                                                                  
!       Outlet air dry bulb temperature                                           
                                                                                  
        tao=bf*(tai-ts)+ts                                                        
                                                                                  
!       Air side sensible heat ratio                                              
                                                                                  
        shr=(tai-tao)*cpi/(hi-ho)                                                 
        if(shr>0.999) then
          shr=1.0
        endif
                                                                                  
!       Cooling and dehumidifying coil model:                                     
!  ---  ----------------------------------                                        
!       r      : equivalent overall thermal resistance                            
!       antu   : equivalent number of transfer units                              
!       crr    : equivalent thermal capacity rate ratio                           
!       effect : equivalent air side effectiveness                                
                                                                                  
        r=shr*ra+rm+rw                                                            
        antu=aext*shr/(cair*r)                                                    
        crr=cair/cwat/shr                                                         
        xnum3=antu*(1-crr)

        if(xnum3>100.0) then
          eo=0.0                                                                  
        else                                                                      
          eo=exp(-xnum3)                                                          
        endif

        effect=(1.0-eo)/(1.0-crr*eo)
        qtnew=cair*effect*(tai-twi)/shr                                           
                                                                                  
        if(abs(qtnew-qt)<eps) go to 50
                                                                                  
        qt=qtnew                                                                  
      end do

!     No convergence apres niter iterations

      write(ifile,9001)time,niter
9001  format(' At time',f10.2,'  no convergence after',i5,' iterations')
!*      write(ifile,9001)time,qtnew,shr
!*9001  format(' At time',f6.2,'no convergence after niter iterations',&
!*       //,' total heat transfer =',f6.2,'shr =',f6.2)
      go to 1000

!     Results printing

50    two=twi+qtnew/cwat
      ho=hi-cpi*effect*(tai-twi)/shr
      qsensi=shr*qtnew
      tao=tai-qsensi/cair

      if(tao>200.0) then
        tao=200.0
      elseif(tao<0.0) then
        tao=0.0
      endif

      hs=(ho-bf*hi)/(1.0-bf)
      ts=ftsat(hs,patm)
      go=fwha(tao,ho)
      at=tao+273.15

      if(tao<=0.0) then
        pws=exp(c1/at+c2+at*(c3+at*(c4+at*(c5+at*c6)))+c7*alog(at))
      else
        pws=exp(c8/at+c9+at*(c10+at*(c11+at*c12))+c13*alog(at))
      endif

      rho=100.0*(patm*go)/(pws*(0.62198+go))

      if(rho.gt.100.0) then
        write(ifile,9000)time,rho
9000    format(' At time ',f10.2,' rho = ',f6.2,' higher than 100%')
!*       '  imperfection of the simulation method, ',
!*       'out(3)=rho is not put at 100%')
      endif

      au=aext/r
      rewat=3.6*am1*(gi-go)

      yout(1)= tao
      yout(2)= rho
      yout(3)= go
      yout(4)= ho
      yout(5)= two
      yout(6)= qtnew
      yout(7)= qsensi
      yout(8)= rewat
      yout(9)= shr
      yout(10)=effect
      yout(11)=bf
      yout(12)=au

      do i=1,12
        iostat(i)=1
      end do

      return

!     The cooling coil is not required

1000  yout(1)= tai
      yout(2)= 0.0
      yout(3)= gi
      yout(4)= hi
      yout(5)= twi
      yout(6)= 0.0
      yout(7)= 0.0
      yout(8)= 0.0
      yout(9)= 0.0
      yout(10)=0.0
      yout(11)=0.0
      yout(12)=0.0

      do i=1,12
        iostat(i)=1
      end do

      return
      end subroutine type30

!***********************************************************************
!
        subroutine type33(xin,yout,par_v,saved_v,iostat)
!
! ----------------------------------------------------------------------
!
!       TYPE 33 : MOIST AIR MIXING DAMPERS AND MERGER MODEL
!
!***********************************************************************

      implicit none
      integer, parameter                :: ni=8, no=5, np=3, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      integer                           :: i
      real                              :: k
      real                              :: p3,t1,t2,h1,h2,c,cl,wf,dp,&
                                           s,r1,r2,w1,p2,w3,t3,h3,p1,w2
      p1=  xin(1)
      w2=  xin(2)
      p3=  xin(3)
      t1=  xin(4)
      t2=  xin(5)
      h1=  xin(6)
      h2=  xin(7)
      c=   xin(8)

      k=   par_v(1)
      cl=  par_v(2)
      wf=  par_v(3)

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

      r1=wf*k/((1.-cl)*c+cl)**2 + (1.-wf)*k*cl**(2.0*c-2.0)
      r2=wf*k/((1.-cl)*(1.-c)+cl)**2 + (1.-wf)*k*cl**(-2.0*c)
      w1=s*sqrt(abs(dp)/r1)

      if(w2/=0.0) then
        s=sign(1.0,w2)
      else
        s=1.0
      endif

      p2=p3+s*r2*w2*w2
      w3=w1+w2

      if(w3>0.0) then
        t3=(w1*t1+w2*t2)/w3
        h3=(w1*h1+w2*h2)/w3
      else
        t3=0.5*(t1+t2)
        h3=0.5*(h1+h2)
      endif

      yout(1)=w1
      yout(2)=p2
      yout(3)=w3
      yout(4)=t3
      yout(5)=h3

      do i=1,5
        iostat(i)=1
      end do

      return
      end subroutine type33

!***********************************************************************
!
      subroutine type34 (xin,yout,par_v,saved_v,iostat)
!
! ----------------------------------------------------------------------
!
!      TYPE 34 : MULTIPLIER
!                state variables are multiplied by given factors
!
!***********************************************************************

      implicit none
      integer, parameter                :: ni=1, no=1, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      yout(1)=xin(1)*par_v(1)
      iostat(1)=1

      return
      end subroutine type34

!***********************************************************************
!
      subroutine type35 (xin,yout,par_v,saved_v,iostat)
!
! ----------------------------------------------------------------------
!
!     TYPE 35 :  MEAN VALUES OF TEMPERATURES AND HUMIDITY RATIO
!
!***********************************************************************

      implicit none
      integer, parameter                :: ni=9, no=2, np=3, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real                              :: w1,w2,w3,t1,t2,t3,h1,h2,&
                                           h3,f1,f2,f3,tmean,hmean     
      real                              :: eps=1.0e-4

!     Inputs

      w1=xin(1)
      w2=xin(2)
      w3=xin(3)
      t1=xin(4)
      t2=xin(5)
      t3=xin(6)
      h1=xin(7)
      h2=xin(8)
      h3=xin(9)

      f1=par_v(1)
      f2=par_v(2)
      f3=par_v(3)

!     Compute mean values using mass flow rate as a weighting factor

      if(w1>=eps .and. w2>=eps .and. w3>=eps) then
        tmean=(f1*w1*t1+f2*w2*t2+f3*w3*t3)/(f1*w1+f2*w2+f3*w3)
        hmean=(f1*w1*h1+f2*w2*h2+f3*w3*h3)/(f1*w1+f2*w2+f3*w3)
      else
        tmean=(f1*t1+f2*t2+f3*t3)/(f1+f2+f3)
        hmean=(f1*h1+f2*h2+f3*h3)/(f1+f2+f3)
      endif

!     Outputs

      yout(1)=tmean
      yout(2)=hmean

!     Output status

      iostat(1)=1
      iostat(2)=1

      return
      end subroutine type35

!***********************************************************************
!
        subroutine type36(xin,yout,par_v,saved_v,iostat)
!
! ----------------------------------------------------------------------
!
!       TYPE 36 : RESET SCHEDULE
!
!       Summer schedule for resetting the set-point of chilled water
!       valve controller  of air handling units.  the contol scheme
!       to cutoff the chilled water valve during nighttime is also
!       implemented.
!
! ---------------------------------------------------------------------
!
!       ti      zone air dry-bulb temperature (c)
!       fansp   fan rotational speed (1/s)
!       ts      supply air dry-bulb temperature (c)
!       cntf    conversion factor for a controller using ts (-)
!       cntset  set-point value of a controller (-)
!       tstart  start time of implementing summer schedule of reset (s)
!       fansp0  change-over fan rotational speed (1/s)
!
!       Note:   It is recommended to implement this schedule at tstart,
!               which is the time when resonable convergence is obtained
!               during a simulation.
!
!       December 22, 1986 Cheol Park
!
!       March 11, 1987    Revised to include the nighttime cutoff scheme
!
!***********************************************************************

      use modsim_head
      implicit none
      
      integer, parameter                :: ni=2, no=1, np=7, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real                              :: ti,fansp,ti1,ts1,ti2,ts2,&
                                           cntf,tstart,fansp0,ts,cntset

!     Inputs

      ti     =xin(1)
      fansp  =xin(2)

      ti1    =par_v(1)
      ts1    =par_v(2)
      ti2    =par_v(3)
      ts2    =par_v(4)
      cntf   =par_v(5)
      tstart =par_v(6)
      fansp0 =par_v(7)

!     Reset schedule of summer mode operation of air handling units

      if(time>tstart) then
        if(ti<ti1) then
          ts=ts1
        elseif(ti>=ti1 .and. ti<=ti2) then
          ts=ts1+(ti-ti1)*(ts2-ts1)/(ti2-ti1)
        else
          ts=ts2
        endif
        cntset=ts*cntf
      else
        cntset=ts1*cntf
      endif

!     Cut off the chilled water valve during nighttime if the
!     rotational speed of fan is set to be lower than a specified
!     value, fansp0.

      if(fansp<fansp0) then
        cntset=0.5
      endif

!     Outputs

      yout(1)=cntset
      iostat(1)=1

      return
      end subroutine type36

!***********************************************************************
!
        subroutine type39(xin,yout,par_v,saved_v,iostat)
!
! ----------------------------------------------------------------------
!
!       TYPE 39   TIME-OF-DAY CONTROL  WITH ZONE DEMAND RESET
!
!       The cooling coil valves, fans, and dampers are controlled by
!       the schedule given in boundary.dat file. the control status
!       is indicated by 0 for off and 1 for on, respectively.
!
!       May 19, 1987  Cheol Park
!       Revised:  August 7, 1987 C.P.
!
! ---------------------------------------------------------------------
!
!       ti      zone air dry-bulb temperature (C)
!       ti1     lower limit of zone air dry-bulb temperature (C)
!       ti2     higher limit of zone air dry-bulb temperature (C)
!       ts      supply air dry-bulb temperature (C)
!       ts1     higher set point of  supply air dry-bulb temperature (C)
!       ts2     lower set point of supply air dry-bulb temperature (C)
!       cntf    conversion factor for a controller using ts (-)
!       cntset  set-point value of a controller (-)
!       tstart  start time of implementing summer schedule  (s)
!       fansup  supply fan rotational speed (1/s)
!       fansu1  supply fan rotational speed at full power (1/s)
!       fanrt   return fan rotational speed (1/s)
!       fanrt1  return fan rotational speed at full power (1/s)
!       swfan   switch for fan, 0 for off and 1 for on
!       swccv   switch for cooling coil valve, 0 for off and 1 for on
!       swdamp  switch for dampers, 0 for regular control and 1 for
!               purging cycle
!       toa     outdoor dry-bulb temperature (C)
!       tdamp   change-over temperature for damper control (C)
!       damper  outdoor air damper opening
!       tai     inlet air temperature of the cooling coil (C)
!
!
!       Note:   It is recommended to implement this schedule at tstart,
!               which is the time when resonable convergence is obtained
!               during a simulation.
!
!***********************************************************************

      use modsim_head
      implicit none

      integer, parameter                :: ni=6, no=4, np=9, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      integer      :: i
      real         :: ti,swfan,swccv,swdamp,tai,ti1,ts1,ti2,ts2,&
                      cntf,tstart,fansu1,fanrt1,tdamp,fansup,fanrt,&
                      ts,cntset,damper

!     Inputs

      ti     =xin(1)
      swfan  =xin(2)
      swccv  =xin(3)
      swdamp =xin(4)
      toa    =xin(5)
      tai    =xin(6)

!     Parameters

      ti1    =par_v(1)
      ts1    =par_v(2)
      ti2    =par_v(3)
      ts2    =par_v(4)
      cntf   =par_v(5)
      tstart =par_v(6)
      fansu1 =par_v(7)
      fanrt1 =par_v(8)
      tdamp  =par_v(9)

!     Supply and exhaust fan control  which allows a minimum
!     flow rate of 5 % of full flow rate for simulation purpose.

      if(swfan>0.5 .and. swfan<=1.0) then
        fansup=fansu1*swfan
        fanrt=fanrt1*swfan
      elseif(swfan>=0.0 .and. swfan<=0.5) then
        fansup=(0.05+0.9*swfan)*fansu1
        fanrt=(0.05+0.9*swfan)*fanrt1
      endif

!     Cooling coil valve control

      if(swccv>0.5 .and. swccv<=1.0) then

!     Reset schedule of summer mode operation of air handling units

        if(time>tstart) then
          if(ti<ti1) then
            ts=ts1
          elseif(ti>=ti1 .and. ti<=ti2) then
            ts=ts1+(ti-ti1)*(ts2-ts1)/(ti2-ti1)
          else
            ts=ts2
          endif

!     When the inlet air temperature of the cooling coil is lower
!     than the higher set point of supply air temperature, close the
!     cooling coil valve by raising the set point.

          if(tai<ts1) then
            cntset=3.0*ts1*cntf*swccv
          else
            cntset=ts*cntf*swccv
          endif
        else
          cntset=ts1*cntf*swccv
        endif

!       Set high value of set point to close a cooling coil valve
!       by multiplying a factor of 3.

      elseif(swccv>=0.0 .and. swccv<=0.5) then
        cntset=3.0*ts1*cntf*(1.0-swccv)
      endif

!     Damper control

!     For purging cycle

      if(swdamp>0.5 .and. swdamp<=1.0) then
        if(toa<=tdamp .and. toa>ts2) then
          damper=1.0*swdamp
        else
          damper=0.0
        endif

!     For regular economizer cycle based on NIST Admin. Building control
!     scheme.

      elseif(swdamp>=0.0 .and. swdamp<=0.5) then
        if(toa<=ts2 .or. toa>=tdamp) then
          damper=swdamp
        elseif(toa<tdamp .and. toa>=ts1) then
          damper=1.0
        elseif(toa<ts1 .and. toa>ts2) then
          damper=(toa-ts2)/(ts1-ts2)
        endif

!       If the cooling coil is off, the damper opening is set to
!       minimum regardless the damper control scheme.

        if(swccv>=0.0 .and. swccv<=0.5) then
          damper=swdamp
        endif
      endif

!     Outputs

      yout(1)=fansup
      yout(2)=fanrt
      yout(3)=cntset
      yout(4)=damper

      do i=1,4
        iostat(i)=1
      end do

      return
      end subroutine type39

