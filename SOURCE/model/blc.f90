! ***********************************************************************

        subroutine type62(xin,yout,par_v,saved_v,iostat)

! -----------------------------------------------------------------------
!
!     TYPE 62: HOT WATER BOILER WITH A DOMESTIC HOT WATER HEATING COIL
!
!     January 21, 1988     Cheol Park
!     Updated:    June 22, 1988
!     Updated:    Convert to Fortran 90 code, April 27, 2007 C.P.
!
!     INPUTS:
!       tstk   Stack gas temperature (C)
!       tra    Room air temperature (C)
!       toa    Outdoor air temperature (C)
!       tblw   Boiler water temperature (C)
!       psw    Supply boiler water pressure (kPa)
!       wsw    Supply boiler water flowrate (kg/s)
!       trw    Return boiler water temperature (C)
!       qcw    Heat flow rate to the tankless coil (kW)
!       swch   Switch for controlling on/off
!              1 for burner on and 0 for burner off
!
!     OUTPUTS:
!       tstk   Stack gas temperature (C)     --- Diff. eq.
!       tblw   Boiler water temperature (C)  --- Diff. eq.
!       tsw    Supply water temperature (C)
!       qhwt   Heat flow to water (kW)
!       qstk   Heat loss from the stack (kW)
!       qinput Input heat flow (kW)
!       qltnt  Latent heat loss (kW)
!       qjakt  Heat loss from the boiler jacket (kW)
!       qss    Heat flow to boiler water (kW)
!
!     PARAMETERS:
!
!       vgf    Volume of gas in the boiler fire-box (m3)
!       asf    Boiler fire-box effective radiation heat transfer area (m2)
!       arf    Boiler fire-box refractory surface area (m2)
!       lgpf   Gas-path length in the fire-box (m)
!       emisf  Fire-box surface emissivity in fraction
!       vwb    Volume of water in the boiler (m3)
!       ajakt  Boiler jacket surface area (m2)
!       ujakt  Boiler jacket U-factor (kW/m2-C)
!       carb   Atomic ratio of carbon in fuel
!       hydr   Atomic ratio of hydrogen in fuel
!       oxyg   Atomic ratio of oxygen in fuel
!       xntr   Atomic ratio of nitrogen in fuel
!       sulf   Atomic ratio of sulfur in fuel
!       hfuel  Fuel higher heating value (kJ/kg)
!       cpfuel Fuel specific heat value (kJ/kg-C)
!       wfuel  Fuel supply rate (kg/s)
!       tfuel  Fuel temperature (C)
!       xair   Excess air for combustion (-)
!       effyss Steady-state boiler efficiency in fraction (-)
!       tblws  Boiler water temperature at steady state (C)
!       tsgs   Stack gas temperature at full load (C)
!       icfgon Integration control factor for on-period stack gas temp.
!       icfgof Integration control factor for off-period stack gas temp.
!       icfw   Integration control factor for water temperature
!       dtfg   Short start period for burner on/off cycle (s)
!       icfgns Integration control for short start on-period
!       icfgfs Integration control for short start off-period
!       dtbw   Minimum time interval for updating QSS (s)
!
! ***********************************************************************
!
!     REFERENCE:  Park, C. and G.E. Kelly, A Study on the Performance of
!                 Residential Boilers for Space and Domestic Hot Water
!                 Heating, NISTIR 89-4104, NIST, June 1989.
!
! ***********************************************************************
!
      use modsim_head
      use blc_head
      implicit none
      integer,parameter                 :: ni=9, no=9, np=28, ns=5,&
                                           ndeq=2
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      logical     :: comint=.true.
      real        :: icfgon,icfgof,icfw,icfgns,icfgfs
      real        :: tstk,tra,tblw,psw,wsw,trw,qcw,swch,carb,hydr,oxyg,&
                     xntr,sulf,xair,effyss,tblws,dtfg,dtbw,tsw,tsws,&
                     wgon,texgf,qfc,wgoff,qhxc,dummy,qf,ionoff,capg,&
                     tafc,qfcr,cpf,capw,wcp,qinput,qltnt,qjakt,qhwt,&
                     qstk,qss,qss1,iburn,time0,delt,time1,deltb,dtstk,&
                     tinf,dtblw
      integer     :: i                
      real        :: wtfac=0.8,tblref=97.5,cpa=1.0

!       namelist  /nam1/ qinp
!       namelist  /nam2/ time,swch,qss,tstk,tblw,qf,qhxc,&
!                         dtstk,dtblw,texgf,wgoff,wttl
!       namelist  /nam3/ time,ionoff,time0,delt,icfgon,icfgof
!       namelist  /nam4/ time,time1,deltb,qss,qss1
!       namelist  /nam5/ time,tblw,tstk,qcw
!       namelist  /nam6/ time,qinput,qltnt,qstk,qjakt,qhwt,qss

!     Inputs:

      tstk    =xin(1)
      tra     =xin(2)
      toa     =xin(3)
      tblw    =xin(4)
      psw     =xin(5)
      wsw     =xin(6)
      trw     =xin(7)
      qcw     =xin(8)
      swch    =xin(9)

!     Parameters:

      vgf     =par_v(1)
      asf     =par_v(2)
      arf     =par_v(3)
      lgpf    =par_v(4)
      emisf   =par_v(5)
      vwb     =par_v(6)
      ajakt   =par_v(7)
      ujakt   =par_v(8)
      carb    =par_v(9)
      hydr    =par_v(10)
      oxyg    =par_v(11)
      xntr    =par_v(12)
      sulf    =par_v(13)
      hfuel   =par_v(14)
      cpfuel  =par_v(15)
      wfuel   =par_v(16)
      tfuel   =par_v(17)
      xair    =par_v(18)
      effyss  =par_v(19)
      tblws   =par_v(20)
      tsgs    =par_v(21)
      icfgon  =par_v(22)
      icfgof  =par_v(23)
      icfw    =par_v(24)
      dtfg    =par_v(25)
      icfgns  =par_v(26)
      icfgfs  =par_v(27)
      dtbw    =par_v(28)

!     Assume that the boiler supply water temperature is the same
!     as the boiler water temperature, i.e. tsw=tblw.

      tsw=tblw

!     Initial conditions at steady state.   Input data at
!     beginning of simulation  must be entered using a steady-state
!     values.

      if(itime==1 .and. comint) then
        pt=1.
        tstk=tsgs
        tblw=tblws
        tsws=tblws
        call  blinit(carb,hydr,oxyg,xntr,sulf,xair,&
          tra,tblws,tsws,effyss,wgon)
        comint=.false.
        saved_v(1)=time
        saved_v(2)=wgon
        saved_v(3)=0.0
        saved_v(4)=1
!       print nam1
      else
        wgon=saved_v(2)
      endif

!     Burner off

      if(time > saved_v(1)) then

        if(swch<=0.5) then
          call blfoff(tstk,tra,tblw,wgon,wgoff,qfc,texgf)
          call blhx(tsw,texgf,tstk,wgoff,qhxc)
          call cpcva(tra,cpa,dummy,dummy,dummy)
          qf=qfc
          ionoff=0
          capg=wgoff*cpa

!     Burner on
                                                               
        else
          call blfon(tra,tblw,wgon,qfcr,tafc,texgf)
          call blhx(tsw,texgf,tstk,wgon,qhxc)
          qf=qfcr
          ionoff=1
          capg=wgon*cpf(tra,tstk)
        endif
      endif

!     Heat fluxes

      capw=wsw*wcp(tblw)

      qinput=wfuel*hfuel*ionoff
      qltnt=2442.*wfuel*rwtf*ionoff
      qjakt=ajakt*ujakt*(tblw-tra)
      qhwt=capw*(tsw-trw)+qcw
      qstk=qinput-qltnt-qf-qhxc
      qss=qf+qhxc-qjakt-qcw
      qss1=saved_v(3)

!     Setting up a short time interval of rapid change of gas
!     temperature  and  using different integration constants

      if(time > saved_v(1)) then
        iburn=saved_v(4)+0.001
        if(ionoff /= iburn) then
          time0 = time
        else
          delt = time-time0
        endif
        saved_v(4) = ionoff
      endif

      if(delt < dtfg) then
        icfgon = icfgns
        icfgof = icfgfs
      endif

!     Introducing time-delay to the boiler water temperature,
!     qss is updated every time-interval of deltb, if deltb is
!     greater than or equal to the specified minimum value of dtbw.
!     At the same time, the reference time,time1,
!     is also changed.  It is assumed that qss holds a steady value
!     in the time-interval.
!
!     Note that deltb depends upon the specified maximum time step,
!     tmax.  Due to this fact, the use of different tmax's could result
!     in different simulation results.

      if(itime == 1) then
        time1=time
        saved_v(5)=time1
      endif

      if(time > saved_v(1)) then
        time1=saved_v(5)
        deltb = time-time1
        if(deltb >= dtbw) then
          saved_v(3)=qss
          time1=time
          saved_v(5)=time1
        endif
      endif

!     Derivatives  of stack gas and boiler water temperatures

      if(ionoff == 1) then
        dtstk=(tsgs-tstk)/(icfgon*taug)
      elseif(ionoff == 0) then
        tinf=wtfac*tblref+(1.-wtfac)*tra
        dtstk=(tinf-tstk)/(icfgof*taug)
      endif

      dtblw=(capw*(trw-tsw)+qss1)/(icfw*tauw)

!     Outputs

      yout(1)  =dtstk
      yout(2)  =dtblw
      yout(3)  =tsw
      yout(4)  =qhwt
      yout(5)  =qstk
      yout(6)  =qinput
      yout(7)  =qltnt
      yout(8)  =qjakt
      yout(9)  =qss

      if(time>saved_v(1)) then
!         print nam2
!         print nam3
!         print nam4
!         print nam5
!         print nam6
      endif
      saved_v(1)=time

      do i=1,no
        if(i<4) then
          iostat(i)=0
        else
          iostat(i)=1
        endif
      end do

      return
      end subroutine type62
! **********************************************************************

      subroutine blinit(carb,hydr,oxyg,xntr,sulf,xair,&
        tra,tblw,tsw,effyss,wgon)

! --------------------------------------------------------------------
!
!       cntu     Modified boiler HX heat transfer number at full load
!       taug     Modified boiler stack gas time constant
!       tauw     Modified boiler water time constant
!
! **********************************************************************

      use blc_head
      implicit none
      real    :: wgon,effyss,tsw,tblw,tra,xair,sulf,xntr,oxyg,hydr,&
                 carb,texgf ,tafc,qfcr,qjakt ,tinhx ,qfns,t1g,t2g,cpg,&
                 cpf,capg

       namelist  /namini/ wgon,qfcr,tafc,texgf,qjakt,qfns,taug,tauw

!     Combustion gas properties

      call prdpp(carb,hydr,oxyg,xntr,sulf,xair)

!     Estimate the gas path lenth in the boiler fire-box by matching
!     the calculated and the measured gas temperatures at the exit of
!     the fire-box.

      call blfon(tra,tblw,wgon,qfcr,tafc,texgf)

!     Steady-state condition

      qjakt=ajakt*ujakt*(tblw-tra)
      tinhx=texgf
      qfns=qfcr-qjakt
      call blfld(wgon,qfns,tinhx,effyss,tsw)

!     Time constants of gas and water

      if(taug<1.e-10) then
        t1g=(tfuel+ratf*tra)/rptf
        t2g=tsgs
        cpg=cpf(t1g,t2g)
        capg=cpg*wgon
        taug=1.0/capg
      endif
      if(tauw<=1.e-10) then
        tauw=4200.*vwb
      endif
!      print namini

      return
      end subroutine blinit

! **********************************************************************

      subroutine blfld(wgon,qfns,tinhx,effyss,tsw)

! ---------------------------------------------------------------------
!
!     Boiler heat exchanger performance at full load.  Calculation
!     of the number of heat transfer unit.
!
! **********************************************************************
!
      use blc_head
      implicit none
      real              :: k,mu,ntu
      real              :: tsw,effyss,tinhx,qfns,wgon,cp,cpf,qhxss,&
                           entu,tave

!       namelist  /namfld/ qhxss,ntu,cntu

!     The specific heat of gas passing through the heat
!     exchanger by using the measured stack gas temperature at
!     a steady state.

      cp=cpf(tinhx,tsgs)

!     The number of transfer units at full load

      qhxss=effyss*qinp-qfns
      entu=1.-qhxss/(cp*wgon*(tinhx-tsw))
      if(entu>1.e-20) then
        ntu=log(1./entu)
      else
        ntu=50.
        cp=cpf(tinhx,tsw)
        tave=0.5*(tinhx+tsw)
        print *,' --- too large efficiency value ---'
      endif
      call prdpr(tave,mu,k)
      cntu=ntu*wgon**0.2*mu**0.4*(cp/k)**0.6
!       print namfld

      return
      end subroutine blfld

! *********************************************************************

      subroutine blhx(tsw,tinhx,texhx,wg,qhxc)

! ---------------------------------------------------------------------
!
!      Boiler heat-exchanger part load performance
!
! *********************************************************************

      use blc_head
      implicit none
      real          :: k,mu,ntu
      real          :: qhxc,wg,texhx,tinhx,tsw,tave,cp,cpf,dummy,cpa

!     The number of transfer unit for boiler heat exchanger
!     at part load condition.

      tave=0.5*(tinhx+texhx)
      if(abs(tinhx-texhx) > 0.001) then
        cp=cpf(tinhx,texhx)
      else
        call cpcva(tave,cpa,dummy,dummy,dummy)
      endif

      call prdpr(tave,mu,k)

!     Convective heat transfer rate

      if(wg> 1.0e-6) then
        ntu=cntu/(wg**0.2*mu**0.4*(cp/k)**0.6)
        qhxc=cp*wg*(tinhx-tsw)*(1.-exp(-ntu))
      else
        qhxc=cp*wg*(tinhx-tsw)
      endif

      return
      end subroutine blhx

! **********************************************************************

       subroutine blfon(tra,tblw,wgon,qfcr,tafc,texgf)

! ----------------------------------------------------------------------
!
!     Simulation of boiler fire-box performance during on-period
!
!       pt     Gas pressure (atm)
!       pco2   Number of moles of CO2
!       ph2o   Number of moles of H2O
!       pso2   Number of moles of SO2
!       po2    Number of moles of O2
!       pn2    Number of moles of N2
!       ratf   Mass ratio of air to fuel
!       rwtf   Mass ratio of water to fuel
!       rptf   Mass ratio of combustion product to fuel
!       ts     Fire-box surface temperature (C)
!       qfcr   Fire-box heat transfer rate (kW)
!       tafc   Adiabatic flame temperature (C)
!       texgf  Fire-box exit gas temperature (C)
!
! *********************************************************************

      use blc_head
      implicit none
      real           :: k,mu
      real           :: sigma=5.670e-11,ckelvn=273.15,pi=3.14159
      real           :: texgf,tafc,qfcr,wgon,tblw,tra,hhv,ts,tks,tka,&
                        wh2o,hf,atotal,cs,al,taff,taf,tg2,tg2c,tgave,&
                        cpoe,cpf,cpae,cpoa,tt,pr,ac,dh,re,hcov,hcovf,&
                        ge,gef,ags,gs,tratio,capwoe,capwoa,ccoat,ccoe,&
                        tgcal,dtg2,dtg1,tg1,tg
      integer        :: itr
!     Fuel input

      hhv=hfuel+cpfuel*(tfuel-tra)
      qinp=wfuel*hhv

!     Set the fire-box surface tmperature to be the same as the boiler
!     water temperature.

      ts=tblw

!     Convert temperature unit into absolute unit

      tks=ts+ckelvn
      tka=tra+ckelvn

!     The gas and water vapor flow rates during on-cycle

      wgon=wfuel*rptf
      wh2o=wfuel*rwtf

      hf=wfuel*hhv-2442.*wh2o
      atotal=asf+arf
      cs=asf/atotal
      al=3.5*vgf/atotal

!     Adiabatic flame temperature

      tafc=taff(hhv,tra)
      taf=tafc+ckelvn

!     The heat transfer rate due to radiation and
!     convection using the well-stirred furnace theory

      tg2=taf-250
      itr=0
10    itr=itr+1
      tg2c=tg2-ckelvn
      tgave=0.5*(taf+tg2)
      cpoe=cpf(tra,tg2c)
      cpae=cpf(tg2c,tafc)
      cpoa=cpf(tra,tafc)
      tt=tgave-ckelvn
      call prdpr(tt,mu,k)
      pr=mu*cpae/k
      ac=vgf/lgpf
      dh=sqrt(4.0*ac/pi)
      re=wgon*dh/(mu*ac)
      hcov=hcovf(re,pr,k,dh)
      ge=gef(pco2,ph2o,al,tt)
      ags=gs(atotal,cs,ge,emisf)
      tratio=tks/tgave
      ags=ags*(1.-tratio**3)/(1.-tratio**4)
      ags=ags+2.*asf*hcov/(sigma*(tgave+tks)**3)
      qfcr=ags*sigma*(tgave**4-tks**4)

!     A new estimate of gas temperature by using Newton's method

      capwoe=cpoe*wgon
      capwoa=cpoa*wgon
      ccoat=capwoa*(taf-tka)*(hf-qfcr)
      ccoe=capwoe*hf
      tgcal=ccoat/ccoe+tka
      dtg2=tgcal-tg2
      if(abs(dtg2)>10.) then
        if(itr>1) then
          tg=tg1-dtg1*(tg2-tg1)/(dtg2-dtg1)
          if(abs(dtg1)>abs(dtg2)) then
            dtg1=dtg2
            tg1=tg2
          endif
          tg2=tg
        else
          dtg1=dtg2
          tg1=tg2
          tg2=tg2+100.
          if(dtg1<0.) then
            tg2=tg2-200.
          endif
        endif
        goto 10
      endif

!     Gas temperature at boiler fire-box during on-period

      texgf=tg2-ckelvn

      return
      end subroutine blfon

! **********************************************************************

      subroutine blfoff(tstk,tra,tblw,wgon,wgoff,qfc,texgf)

! --------------------------------------------------------------------
!
!     Boiler fire-box performance during off-period
!
! **********************************************************************

      use blc_head
      implicit none
      real            :: k,mu,ntu
      real            :: ckelvn=273.15,pi=3.14159
      real            :: texgf,qfc,wgoff,wgon,tblw,tra,tstk,tave,&
                         dummy,cpa,pr,ac,dh,re,hcov,hcovf,capa

!     Mass flow rate of air during off-period

      if(tstk-tra > 0.01) then
        wgoff=wgon*((tstk-tra)/(tsgs-tra))**0.56&
          *((tsgs+ckelvn)/(tstk+ckelvn))**1.19
      else
        wgoff=0.0
      endif

!     The properties of gas in the boiler fire-box:
!     the dynamic viscosity, thermal conductivity, specific heat
!     capacity, Prantl number, and Reynolds number based on the
!     hydraulic diameter.

      tave=0.5*(tra+texgf)
      call prdpr(tave,mu,k)
      call cpcva(texgf,cpa,dummy,dummy,dummy)
      pr=mu*cpa/k
      ac=vgf/lgpf
      dh=sqrt(4.0*ac/pi)
      re=wgoff*dh/(mu*ac)

!     The convective heat transfer coefficient.

      hcov=hcovf(re,pr,k,dh)

!     The off-cycle heat transfer rate using the
!     effectiveness method of a heat exchanger.
!     the gas temperature at the exit of the boiler fire-box.

      capa=cpa*wgoff
      if(wgoff > 0.001) then
        ntu=asf*hcov/capa
        qfc=capa*(tra-tblw)*(1.-exp(-ntu))
        texgf=tra-qfc*(1./capa)
      else
        qfc=0.0
        texgf=(tra+tblw)/2.
      endif

      return
      end subroutine blfoff


! **********************************************************************

        subroutine type63(xin,yout,par_v,saved_v,iostat)

! ----------------------------------------------------------------------
!
!     TYPE 63: HOT WATER COIL WITH CONSTANT WALL TEMPERATURE
!
!     February 16, 1988 Cheol Park
!     Revised : April 15, 1988
!     Updated: Convert to Fortran 90 code, April 27, 2007 C.P.
!
!     INPUTS:
!       tblw   Boiler water temperature (C)
!       wcw    Water flow rate through the coil (kg/s)
!       ticw   Inlet coil water temperature (C)
!
!     OUTPUTS:
!       tocw   Outlet coil water temperature (C)
!       qcw    Heat flow rate from boiler water to coil water (kW)
!
!     PARAMETERS:
!       dcoil  Diameter of coil
!       lcoil  Length of coil
!
! **********************************************************************
!
      use modsim_head
      implicit none
      integer,parameter                 :: ni=3, no=2, np=2, ns=2
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real         :: k,lcoil,lmtd,mu
      real         :: pi=3.14159
      real         :: tblw,wcw,ticw,dcoil,ts,tave,tocw,wmu,wk,cp,&
                      wcp,asx,ac,re,pr,hcov,hcovf,ax,qcw

!       namelist /nam1/ hcov,tocw,lmtd,qcw
!       namelist /nam2/ time,tblw,wcw,ticw,tocw,qcw

!     Inputs:

      tblw   =xin(1)
      wcw    =xin(2)
      ticw   =xin(3)

!     Parameters:

      dcoil  =par_v(1)
      lcoil  =par_v(2)

!     Temporally set the wall temperature of the water coil equal to
!     the boiler water temperature.

      ts=tblw

!     Water properities:  Dynamic viscosity, thermal conductivity, and
!     specific heat of water from the subroutine watpr with average
!     temperature except at the beginning.

      if(itime==1) then
        tave=0.5*(tblw+ticw)
        saved_v(1)=time
        saved_v(2)=ticw
      else
        tocw=saved_v(2)
        tave=0.5*(tocw+ticw)
      endif
      mu=wmu(tave)
      k=wk(tave)
      cp=wcp(tave)

!     the coil surface area, and wetted cross-sectinal area.

      asx=pi*dcoil*lcoil
      ac=pi*dcoil*dcoil/4.0

!     reynolds number and prandltl number of the coil water

      re=wcw*dcoil/(ac*mu)
      pr=mu*cp/k
      hcov=hcovf(re,pr,k,dcoil)

!     the outlet temperature and the heat flow rate when
!     the coil surface temperature is constant.

      if(wcw> 0.0001) then
        ax=hcov*asx/(wcw*cp)
        tocw=ts-(ts-ticw)*exp(-ax)
        lmtd=(tocw-ticw)/log((ts-ticw)/(ts-tocw))
        qcw=hcov*asx*lmtd
      else
        tocw=ts
        qcw=0.0
      endif

!     outputs:

      yout(1)  =tocw
      yout(2)  =qcw
      saved_v(2)=tocw

      if(time > saved_v(1)) then
!         print nam1
!         print nam2
      endif
      saved_v(1)=time

      iostat(1)=0
      iostat(2)=0

      return
      end subroutine type63

! **********************************************************************

        subroutine type64(xin,yout,par_v,saved_v,iostat)

! -----------------------------------------------------------------------
!
!     TYPE 64: BOILER BURNER AND CIRCULATING PUMP CONTROL
!
!     February 17, 1988 Cheol Park
!     Revised : June 6, 1988
!
!     INPUTS:
!       tburn  Burner control temperature (C)
!       tblw   Boiler water temperature (C)
!       tload  Load temperature (C)
!       ctherm Thermostat control indicator for the burner
!               1 for thermostat control, 0 for manual control (-)
!       wcw    Domestic hot water flow rate (kg/s)
!
!     OUTPUTS:
!       swch   Burner control on/off signal (-)
!       wsw    Boiler circulating water flow rate (kg/s)
!
!     PARAMETERS:
!       tburon Boiler burner-on temperature (C)
!       tburof Boiler burner-off temperature (C)
!       tpmpon Boiler water circulating pump-on temperature (C)
!       tpmpof Boiler water circulating pump-off temperature (C)
!       wpump  Circulating water flow rate (kg/s)
!       rateof Minimum water flow rate ratio
!                                  during off-period of pump (-)
!
! **********************************************************************

      use modsim_head
      implicit none
      integer,parameter                 :: ni=5, no=2, np=6, ns=1
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat
      real         :: tburn,tblw,tload,ctherm,wcw,tburon,tburof,&
                      tpmpon,tpmpof,wpump,rateof,swch,wsw

       namelist  /nam1/ time,tburof,tburon,tblw,swch,wcw

!     Inputs:

      tburn  =xin(1)
      tblw   =xin(2)
      tload  =xin(3)
      ctherm =xin(4)
      wcw    =xin(5)

!     Parameters:

      tburon =par_v(1)
      tburof =par_v(2)
      tpmpon =par_v(3)
      tpmpof =par_v(4)
      wpump  =par_v(5)
      rateof =par_v(6)

      if(itime == 1) then
        saved_v(1)=time
      endif

!     Boiler burner control

!     The burner is controlled by the thermostat sensing the boiler
!     water temperature if the indicator, ctherm, is greater than or
!     equal to 0.5.  Otherwise tburn controls the burner.
!     The values of ctherm and tburn must be present in the boundary
!     data file.

      if(ctherm >= 0.5) then
        tburn=tblw
      endif

!     Turn on and off the burner and the pump only when the time step
!     varies to prevent false action due to numerical instability
!     during the large change of a variable.

      if(time > saved_v(1)) then
        if(tburn >= tburof) then
          swch=0.0
        elseif(tburn <= tburon) then
          swch=1.0
        endif

!     Boiler water circulating pump control

        if(tload <= tpmpon) then
          wsw=wpump
        elseif(tload >= tpmpof) then
          wsw=rateof*wpump
        endif
      endif

!     print nam1
!     Outputs:

      yout(1)  =swch
      yout(2)  =wsw

      iostat(1)=0
      iostat(2)=0

!       if(time > saved_v(1)) then
!         print nam1
!       endif
      saved_v(1)=time

      return
      end subroutine type64

