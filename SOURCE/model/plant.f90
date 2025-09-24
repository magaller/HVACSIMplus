! ********************************************************************

	subroutine type122(xin,yout,par_v,saved_v,iostat)

! **********************************************************************
! *
! * STATIC BOILER                                                      *
! * ref : "boiler specification"         IEA Annex 10   an10 870408-07 *
! *
! **********************************************************************
! * M.Dachelet May 1987                                                *
! **********************************************************************
! * mc      burner fuel mass flow rate                                 *
! * qu      useful power                                               *
! * qc      consumed power                                             *
! * co2     co2 content in gas                                         *
! * eta     efficiency of the boiler                                   *
! * etaon   nominal efficiency                                         *
! * tetaon  fraction of cycle when burner is "on"                      *
! * tetasb  fraction of cycle when burner is "on" during stand-by      *
! * eff     efficiency of equivalent heat exchangers                   *
! * au      au of equivalent heat exchanger                            *
! * mw      water mass flow rate                                       *
! * mw0     nominal water mass flow rate                               *
! * mc      fuel mass flow rate                                        *
! * mc0     nominal fuel mass flow rate                                *
! * twsu    water supply temperature                                   *
! * twex    water exhaust temperature                                  *
! * tamb    ambient temperature                                        *
! * tc      fuel temperature                                           *
! * pci     low calorific value of the fuel                            *
! * ntu     number of transfert units                                  *
! *                                                                    *
! * mfor more information on variables and nomenclature :              *
! * "boiler specifications"     (IEA - Annex 10)                       *
! **********************************************************************
! * note : units                                                       *
! *                                                                    *
! *         For correct run , please use coherent unit system          *
! *                                                                    *
! * examples:                                                          *
! *                                                                    *
! * time           s         s         h         h                     *
! * power          W         kW        kJ/h      kW                    *
! * mass           kg        kg        kg        kg                    *
! * volume         m3        m3        m3        m3                    *
! * energy         J         kJ        kJ        kWh                   *
! *                                                                    *
! **********************************************************************
! * inputs :  1 occupancy ( 0=off  1=on )                              *
! *           2 water supply temperature                               *
! *           3 water mass flow rate                                   *
! *           4 ambient temperature                                    *
! *           5 mean aquastat temperature                              *
! *                                                                    *
! * outputs : 1 water exhaust temperature                              *
! *           2 useful power                                           *
! *           3 consumed power                                         *
! *           4 mean fuel consumption                                  *
! *           5 tetaon                                                 *
! *           6 tetasb                                                 *
! *           7 efficiency                                             *
! *           8 efficiency 100% on                                     *
! *           9 efficiency of the equivalent heat exchanger            *
! *          10 equivalent au                                          *
! *                                                                    *
! * param. :  1 burner fuel mass flow rate                             *
! *           2 concentration of CO2 in gas                            *
! *           3,4,5   au0, k1, k2                                      *
! *           6,7,8   yw, dyw, kw                                      *
! *           9,10,11 mc0, mw0, co20                                   *
! *           12,13   c1,c2                                            *
! *           14,15   c3,c4                                            *
! *           16      cpw                                              *
! *           17,18   cpc,pci                                          *
! *                                                                    *
! **********************************************************************
!                                                                      *
!   proposed values :   c1         3.294      kJ/K kg                  *
!                       c2         2.105      kJ/K kg                  *
!                       c3         0.982      kJ/K kg                  *
!                       c4         2.125      kJ/K kg                  *
!                       cpw        4.187      kJ/K kg                  *
!                       cpc        1.880      kJ/K kg                  *
!                       pci      42875.0      kJ/kg                    *
!                                                                      *
! **********************************************************************
!
!      revised:  Dec. 10, 2002  Cheol Park, NIST
!      updated : January 25, 2007
!            Converted Fortran 77 code into Fortran 90
!
! **********************************************************************

      use modsim_head
      implicit none
      integer, parameter                :: ni=5, no=10, np=18, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real        :: c,o2,u0,yw,dyw,co20,c1,c2,c3,c4,cpc,&
                     pci,twsu,tamb,taqua,tc,cw,twex,qu,qc,tetaon,&
                     tetasb,eta,etaon,eff,au,qb,co2,cpg,cpg0,&
                     cg,cg0,cw0,au0,omega,um0,tav,qumax,cc,cd,cpw,cpa

      real        :: ntu,mw,mc,mw0,mc0,k1,k2,&
                     kb,kd,ke,kw
      integer     :: gamma,i
      logical     :: occupancy, retolog

! ********************
! * Values of inputs *
! ********************

      occupancy = retolog(xin(1))
      twsu  = xin(2)
      mw    = xin(3)
      tamb  = xin(4)
      taqua = xin(5)

! ************************
! * Values of parameters *
! ************************

      mc   = par_v(1)
      co2  = par_v(2)
      au0  = par_v(3)
      k1   = par_v(4)
      k2   = par_v(5)
      yw   = par_v(6)
      dyw  = par_v(7)
      kw   = par_v(8)
      mc0  = par_v(9)
      mw0  = par_v(10)
      co20 = par_v(11)
      c1   = par_v(12)
      c2   = par_v(13)
      c3   = par_v(14)
      c4   = par_v(15)
      cpw  = par_v(16)
      cpc  = par_v(17)
      pci  = par_v(18)
!      if(occupancy) then
         gamma = 1.0
!      else
!         gamma = 0.0
!      endif

      tc=tamb

! **********************
! * No water flow rate *
! **********************

      if (mw <= 0.) then
         yout(1) = twsu
         do i=2,no
            yout(i) = 0.
         enddo
         go to 20
      endif

! ***************************
! * Checks if boiler is off *
! ***************************

      if ((gamma == 0) .or. (twsu >= taqua)) then
         cw   = mw*cpw
         twex = (twsu-(twsu*yw*(1.-kw)/cw)+(yw*tamb/cw))/(1.+(yw*kw/cw))
         qu   = cw*(twex-twsu)
         qc   = 0.
         tetaon = 0.
         tetasb = 0.
         eta    = 0.
         etaon  = 0.
         eff    = 0.
         au     = 0.
         go to 10
      endif

! ************************
! * General calculations *
! ************************

      qb    = mc*pci
      cpg   = c1+c2/co2
      cpg0  = c1+c2/co20
      cpa   = c3+c4/co2
      cg    = mc*cpg
      cg0   = mc0*cpg0
      cw    = mw*cpw
      cw0   = mw0*cpw
      au    = au0/(1.+k1*(cg0-cg)/cg+k2*(cw0-cw)/cw)
      ntu   = au/cg
      omega = cg/cw
      um0   = 1.-omega
      eff   = (1.-exp(-ntu*um0))/(1.-omega*exp(-ntu*um0))
      kb    = (cpg*twsu-(cpa*tamb+cpc*tc))/pci
      twex  = taqua
      qu    = cw*(twex-twsu)

      tav    = kw*twex+(1.-kw)*twsu
      kd     = (dyw*(tav-tamb))/qb
      ke     = (yw*(tav-tamb))/qb
      etaon  = eff*(1.-kb)-kd-ke
      qumax  = qb*etaon
      tetaon = ((qu/qb)+ke)/(eff*(1.-kb)-kd)
      tetasb = ke/(eff*(1.-kb)-kd)

! ************************************
! * Calculation in case of full load *
! ************************************

      if(qu > qumax) then

        cc = cw+((yw+dyw)*kw)
        cd = ((cw*twsu)+(qb*eff*(1.-kb))&
             -((yw+dyw)*(((1.-kw)*twsu)-tamb)))
        twex   = cd/cc
        qu     = cw*(twex-twsu)
        tetaon = 1.
        etaon  = qu/qb

      endif

      qc  = qb*tetaon
      eta = qu/qc

! *********************
! * Values of outputs *
! *********************

   10 yout(1)  = twex
      yout(2)  = qu
      yout(3)  = qc
      yout(4)  = qc/pci
      yout(5)  = tetaon
      yout(6)  = tetasb
      yout(7)  = eta
      yout(8)  = etaon
      yout(9)  = eff
      yout(10) = au
      if(time > 700000.0) then                      ! debug 12/10/02
!         print boiler_data
      endif

   20   do i=1,no
         iostat(i)=1
      enddo

      return
      end subroutine type122

! ********************************************************************

	subroutine type124(xin,yout,par_v,saved_v,iostat)

! ********************************************************************
! *                   *
! *      CHILLER      *
! *                   *
! *********************
! ***********************************************************
! * parameters                                              *
! * ----------                                              *
! *  1) type of compressor          1 = mono-screw          *
! *                                 2 = piston              *
! *                                 3 = centrifugal         *
! *  2) minimum partial cooling load (0-1)                  *
! *  3) design cooling capacity [kW]                        *
! *  4) design compressor power [kW]                        *
! *  5) design chilled water mass flow rate [kg/s]          *
! *  6) design evaporator outlet water temperature [C]      *
! *  7) design final temperature difference between         *
! *     exhaust water and refrigerant in the evaporator [C] *
! *  8) design temperature difference between               *
! *     exhaust water and supply water in the condenser [C] *
! *  9) design condenser inlet water temperature [C]        *
! * 10) design condensation temperature [C]                 *
! * 11) maximum condensation temperature [C]                *
! * 12) proportional band around set point                  *
! *     of evaporator outlet water temperature [C]          *
! *                                                         *
! * inputs                                                  *
! * ------                                                  *
! * 1) occupancy (0/1)                                      *
! * 2) evaporator inlet water temperature [C]               *
! * 3) evaporator water mass flow rate [kg/s]               *
! * 4) condenser inlet water temperature [C]                *
! * 5) condenser mass flow rate [kg/s]                      *
! * 6) evaporator outlet water temperature set point [C]    *
! *                                                         *
! * outputs                                                 *
! * -------                                                 *
! *  1) evaporator oulet water temperature [C]              *
! *  2) condenser outlet water temperature [C]              *
! *  3) evaporation temperature [C]                         *
! *  4) condensation temperature [C]                        *
! *  5) heat transfer in the evaporator [kW]                *
! *  6) compressor power [kw]                               *
! *  7) heat transfer in the condenser [kW]                 *
! *  8) cop                                                 *
! *  9) compressor working status (0,1,2,3)                 *
! * 10) operation time fraction                             *
! **********************************************************************
!
!      revised:  Dec. 10, 2002  Cheol Park, NIST
!      updated : January 25, 2007
!            Converted Fortran 77 code into Fortran 90
!
! **********************************************************************

      use modsim_head
      use cool_plant
      implicit none
      integer, parameter                :: ni=6, no=10, np=12, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real           :: twiev,wmfev,twicd,wmfcd,tset,gamma,twievd,&
                        twocdd,tevd,copcar,etacar,copd,qdcond,wmfcdd,&
                        etaevd,auev,etacdd,aucd,qmax,qmin,qmed,qset,&    
                        twoev,qreg,status,qevap,optime,etaev,tev,etacd,& 
                        tcd,ct,tcds,pl,f,twocd,cop,power,qcond,twoev2,&  
                        twocd2,dtreg,crtemp,plfact,stop_fault,&
                        start_fault,capmin,capch,dpower,wmfevd,&
                        twoevd,deltev,deltc0,twicdd,tcdd,tcdmax
      real           :: kk
      integer        :: i,iter,icompr,itermax=20
      logical        :: retolog, occupancy
      real,dimension(3)      :: eta0 = (/0.63,0.57,0.60/),&                        
                                cap0 = (/0.275,0.25,0.25/)                          
      real                   :: tref = 273.16, tol =0.005
      real                   :: cpw=4.187

      if(time > start_chiller .and. time < stop_chiller) then   ! 12/3/02
         cooling = .true.                                                          
      else                                                                          
         cooling = .false.                                                          
      endif                                                                         

!*************************
! * Values of parameters *
! ************************
      icompr = nint(par_v(1))
      capmin = par_v(2)              
      capch  = par_v(3)              
      dpower = par_v(4)              
      wmfevd = par_v(5)              
      twoevd = par_v(6)              
      deltev = par_v(7)              
      deltc0 = par_v(8)              
      twicdd = par_v(9)              
      tcdd   = par_v(10)             
      tcdmax = par_v(11)             
      dtreg  = par_v(12)             
! ****************
! * Input values *
! ****************
      occupancy = retolog(xin(1))
      twiev = xin(2)
      wmfev = xin(3)                                 
      twicd = xin(4)                                 
      wmfcd = xin(5)                                 
      tset  = xin(6)                                 

      if(cooling) then                ! 12/3/02
         gamma = 1.0
      else
         gamma = 0.0
      endif

      if (time < 3600.) then          ! 12/6/02
        goto  20
      end if

! ********************************
! * Checks if chiller is working *
! ********************************
      if(gamma < 0.5) goto 20
      if(wmfcd<=0.0 .or. wmfev<=0.0 .or. twiev<=tset) goto 20
! ************************************
! * Calculation of design conditions *
! ************************************
      twievd = capch/wmfevd/cpw + twoevd
      twocdd = twicdd + deltc0                                                
      tevd   = twoevd - deltev                                                
      if (capmin<=0.0 .or. capmin>=1.) capmin = cap0(icompr)                  
      copcar = (tevd+tref)/(tcdd-tevd)               ! 12/17/02               
!      copcar = atemp(tevd)/(tcdd-tevd)
      if (dpower<=0) then
         etacar = eta0(icompr)
         copd   = copcar*etacar
         dpower = capch/copd
      else
         copd   = capch/dpower
         etacar = copd/copcar
      endif
      qdcond = capch + dpower
      wmfcdd = qdcond/cpw/deltc0
      etaevd = (twievd - twoevd)/(twievd - tevd)
      auev   = -wmfevd*cpw*alog(1.-etaevd)
      etacdd = (twocdd - twicdd)/(tcdd - twicdd)
      aucd   = -wmfcdd*cpw*alog(1.-etacdd)
! ***********************
! * Temperature control *
! ***********************
!
!         qreg »
!              |
!              |          ______qmax
!              |         /
!              |        /
!              |       /
!              |      /
!         qmin |_____/
!              |
!              -------------------> twoev
!                      tset
!
      qmax   = capch
      qmin   = qmax*capmin
      qmed   = (qmin + qmax)/2.

      qset   = wmfev*cpw*(twiev - tset)
      qreg   = qmed + ((qmax-qmin)/dtreg)*(twoev-tset)      !  12/17/02  Jeff Schein

!      qreg   = wmfev*cpw*(twiev-tset+qmed*dtreg/(qmax-qmin))/&
!               (1.+wmfev*cpw*dtreg/(qmax-qmin))
      if (qreg<qmin) qreg = wmfev*cpw*(twiev-tset+dtreg/2)
      if (qreg<=0.) then
      		    status  = 0.
      		    qevap   = 0.
      		    optime  = 0.
      else
         if (qreg<qmin) then
      		    status  = 1.
      		    qevap   = qmin
      		    optime  = qreg/qmin
         else
      	    if (qreg<qmax) then
      		    status  = 2.
      		    qevap   = qreg
      		    optime  = 1.
      	    else
      		    status  = 3.
      		    qevap   = qmax
      		    optime  = 1.
      	    endif
         endif
      endif
      twoev = twiev - qevap/wmfev/cpw
      if (status==0.0) go to 20
! *********************
! * The chiller is on *
! *********************
! Evaporator
! ----------
      etaev = 1.-exp(-auev/wmfev/cpw)
      tev   = twiev - (twiev - twoev)/etaev
!
! Condenser
! ---------
      etacd  = 1.-exp(-aucd/wmfcd/cpw)
! *************
! * part load *
! *************
! Loop on partial load correction
! -------------------------------
! (takes into account the temperature influence on the partial load factor)
      iter = 0
      tcd = 0.
      ct = 1.
200   continue
      if(iter>0) then
        ct = crtemp(icompr,tcd,tev)
      endif
      tcds = tcd
      pl   = qevap/qmax*ct
      f    = plfact(icompr,pl)
! ***************************************************************
! * tcd is calculated by solving this system of equations:      *
! *   -------------------------------------------------------   *
! *   | QCOND = qevap + POWER                               |   *
! *   | QCOND = (TWOCD - twicd) * cpw * wmfcd               |   *
! *   | qevap = POWER * etacar * f * (tev+tref) / (tev-TCD) |   * 12/17/02
! *   | etacd = (TWOCD - twicd) / (TCD -twicd)              |   *
! *   -------------------------------------------------------   *
! * The capital letters represent the not known variables       *
! ***************************************************************
      kk     = wmfcd*cpw*etacd/qevap
      tcd    = (1.-tev/f/etacar/(tev+tref)+kk*twicd)&
                /(kk-1./f/etacar/(tev+tref))
!      tcd    = (1.-tev/f/etacar/atemp(tev)+kk*twicd)&
!                 /(kk-1./f/etacar/atemp(tev))
      iter   = iter + 1                                                      
      if (abs((tcd-tcds)/tcd)>tol.and.iter<itermax) then   ! itermax = 20
         goto 200
!      else
!         write(*,*) 'type 124: iter exceeds itermax'
!         write(*,*) 'condensor outlet temp: tcd =', tcd
!         stop
      endif

      twocd  = etacd*tcd + twicd*(1.-etacd)                                  
      cop    = f*etacar*(tev+tref)/(tcd-tev)                                 
      power  = qevap/cop                                                     
      qcond  = qevap + power                                                 
      twoev2 = twiev - qevap*optime/wmfev/cpw                                
      twocd2 = twicd + qcond*optime/wmfcd/cpw                                
      go to 30                                                               
! **********************
! * The chiller is off *
! **********************
20    twoev  = twiev
!if(cooling) then
!  print *,' type124: time, twicd, twoev, twocd, cooling =',&
!        time, twicd, twoev, twocd, cooling
!endif
      twocd  = twicd
      twoev2 = twoev          
      twocd2 = twocd          
      tev    = 0.0            
      tcd    = 0.0            
      qevap  = 0.0            
      power  = 0.0            
      qcond  = 0.0            
      cop    = 0.0            
      status = 0.0            
      optime = 0.0            
! *****************
! * Output values *
! *****************
   30 continue

      yout(1) = twoev2
      yout(2) = twocd2                                               
      yout(3) = tev                                                  
      yout(4) = tcd                                                  
      yout(5) = qevap*optime                                         
      yout(6) = power*optime                                         
      yout(7) = qcond*optime                                         
      yout(8)=  cop                                                  
      yout(9)=  status                                               
      yout(10)= optime                                               

      do i=1,no
        iostat(i)=1                                                  
      enddo                                                          

      if(twoev2 > 50.0) then
         print *, ' type124: evaporator inlet water temp. is too high'
         print *, ' twoev, twoev2 =', twoev, twoev2
         stop
      endif

      return
      end subroutine type124

! ************************************
! * Full load temperature correction *
! ************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real function crtemp(icompr,tcd,tev)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      real                  :: tcd,tev,ct,tetae,tetac
      real                  :: teref=2.0, tcref=36.0
      real                  :: a=-0.119363098e-01,&
                               b=0.403618892e-01,&
                               c=-0.270317305e-03
      integer               :: icompr

      select case (icompr)
        case (1)
          ct = 1./(0.0536*tev-0.0130*tcd-0.0001*tev*tcd+1.7887)
        case (2)
          ct = 1./(0.0798*tev-0.0168*tcd-0.0008*tev*tcd+2.0108)
        case (3)
!         values are not available for centrifugal compressor
          ct = 1.
        case (4)
!         correlation cgwa-105r de trane
          tetac = tcd - tcref
          tetae = tev- teref                                     
          ct = 1./(1. + a*tetac + b*tetae + c*tetae*tetac)       
        case default
      end select
      crtemp = ct

      return
      end function crtemp

! ********************
! * part load factor *
! ********************
!+++++++++++++++++++++++++++++++++++++
	real function plfact(icompr,pl)
!+++++++++++++++++++++++++++++++++++++
      integer               :: icompr
      real                  :: pl,f

      select case (icompr)
         case (1)
            if (pl<=0.275) then
               f = 0.46004
            else
               f = pl/(0.7469-1.1368*pl+2.4548*pl**2-1.0649*pl**3)
            endif
         case (2)
            if (pl<=0.25) then
               f = 0.6667
            else
               f = 1.2*pl/(pl+0.2)
            endif
         case (3)
            if (pl<=0.1) then
               f = pl/(0.5*pl+0.2)
            else
               if (pl<=0.5) then
            	f = 8*pl/(5*pl+1.5)
               else
            	f=1.
               endif
            endif
         case (4)
            if (pl<=0.25) then
               f = 0.6667
            else
               f = 1.2*pl/(pl+0.2)
            endif
         case default
      end select

      plfact = f
      return
      end function plfact


! ********************************************************************

	subroutine type143(xin,yout,par_v,saved_v,iostat)

! ********************************************************************
!
! *************************
! *                       *
! *     COOLING TOWER     *
! *                       *
! *************************
!
! Modified by P. Haves, University of Oxford, 6/90, to include three
! port valve, variable speed fan, and sump.  Assumes wind induced
! flow gives 10% duty with fans off at design conditions.
!
! *******************************************************************
! * 1) Webb, R. L., "A unified theoretical treatment for thermal    *
! *    analysis of cooling towers, evaporative condensers and       *
! *    fluid coolers", ASHRAE Transactions, vol. 90, part 2, 1984.  *
! * 2) Webb, R. L., Villacres, A., "Algorithms for performance      *
! *    simulation of cooling towers, evaporative condensers and     *
! *    fluid coolers", ASHRAE Transactions, vol. 90, part 2, 1984.  *
! * 3) Modified and improved by J.L. Nizet, Laboratoire de          *
! *    Physique du Batiment et de Thermodynamique,                  *
! *    University of Liege, March 1985, and V.Giaretto, A. Mazza    *
! *    Dipartimento di Energetica del Politecnico di Torino, Feb.'87*
! * 4) Simplified method contributed by H.C.Pettsman, C.Kwakernaak  *
! *    and H.J.Nicolaas, TNO-Institute of Applied Phisic,           *
! *    June 1986.                                                   *
! *******************************************************************
!
! ********************************************************************
! * parameters                                                       *
! * 1) design supply air wet bulb temperature (C)                    *
! * 2) design supply water temperature (C)                           *
! * 3) design exhaust water temperature (C)                          *
! * 4) design air volume flow rate (m3/s)                            *
! * 5) design water volume flow rate (m3/s)                          *
! * 6) design fan power (kW)                                         *
! * 7) valve authority (parabolic/linear valve) (0-1)                *
! * 8) valve leakage (0-1)                                           *
! * 9) sump volume (m3)                                              *
! *                                                                  *
! * inputs                                                           *
! * 1) on/off switch for pump (0/1)                                  *
! * 2) inlet water temperature (C)                                   *
! * 3) supply water mass flow rate from chiller (C)                  *
! * 4) ambient wet bulb temperature (C)                              *
! * 5) atmospheric pressure (Pa)                                     *
! * 6) fan speed control (fraction of design air flow) (0-1)         *
! * 7) valve stem position (0-1)                                     *
! * 8) temperature of water exiting sump (C)                         *
! *                                                                  *
! * outputs                                                          *
! * 1) temperature of water exiting sump (C)                         *
! * 2) mixed water temperature returned to sump (C)                  *
! * 3) tower outlet water temperature (C)                            *
! * 4) water mass flow rate through tower (kg/s)                     *
! * 5) air mass flow rate (kg/s)                                     *
! * 6) heat transfer to the air (kW)                                 *
! * 7) fan electrical power consumption (kW)                         *
! * 8) water side effectiveness                                      *
! * 9) air side effectiveness                                        *
! **********************************************************************
!
!      revised:  Dec. 10, 2002  Cheol Park, NIST
!      updated : January 25, 2007
!            Converted Fortran 77 code into Fortran 90
!
! **********************************************************************

      use modsim_head
      use cool_plant
      implicit none
      integer, parameter                :: ni=8, no=9, np=9, ns=2,&
                                           ndeq=1
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real              :: twexd,avfld,wvfld,qfan0,avalve,cl,vsump,&
                           msump,etactd,ws,h,wsb,hb,rho,amfld,wmfld,&      
                           rlmind,rlod,rlambd,ck,gamma,twsu,wmfls,&        
                           twbamb,fansp,y1,tsump,y2,aa,hr1,hr2,ct,&        
                           cb,wmfl,avfl,amfl,wvfl,rlo,rlmin,rlamb,&        
                           effw,twex,twmix,qload,effa,qfan,dtsumpdt,&
                           stop_fault,start_fault,w,ta,patm,&
                           ro,twbambd,twsud,fwphi,fhsat,ftwb
      logical           :: onof,wtem,wflo,awbt,totp,fans,valv,&
                           frzde
      real              :: patmd=101325.0,rhow=1000.0,cpw=4.187
      integer           :: i

      ro(patm,ta,w) = patm/(273.16+ta)/(287.06+461.52*w)

      if(time > start_chiller .and. time < stop_chiller) then       !  1/27/03
         cooling = .true.
      else
         cooling = .false.
      endif

! ************************
! * Values of parameters *
! ************************

      twbambd= par_v(1)
      twsud  = par_v(2)
      twexd  = par_v(3)
      avfld  = par_v(4)
      wvfld  = par_v(5)
      qfan0  = par_v(6)
      avalve = par_v(7)
      cl     = par_v(8)
      vsump  = par_v(9)

      msump  = vsump*rhow

! *****************************************************
! * Calculation of the rlo-value at design conditions *
! *****************************************************

      if (itime == 1) then
                          ! first timestep, store ck for use at later
                          ! timesteps
         etactd = (twsud-twexd)/(twsud-twbambd)

         ws  = fwphi(twsud,100.,patmd)
         h   = fhsat(twsud,patmd)
         wsb = fwphi(twbambd,100.,patmd)
         hb  = fhsat(twbambd,patmd)

         rho   = ro(patmd,twbambd,wsb)
         amfld = avfld*rho
         wmfld = wvfld*998.

         rlmind = (cpw*(twsud-twbambd))/(h-hb-cpw*twbambd*(ws-wsb))
         rlod   = amfld/wmfld
         rlambd = rlod/rlmind
         
         if(rlambd > 6.0) then     ! limiting rlambd    5/18/07
            rlambd = 6.0
         elseif(rlambd < 0.0) then            
            rlambd = 0.0
         endif   

         ck = etactd/(1.-exp(-rlambd))
         saved_v(1)=rlod
         saved_v(2)=ck
      else
                          ! use value of ck calculated at first timestep
          rlod=saved_v(1)
          ck=saved_v(2)
      endif

! ****************
! * Input values *
! ****************

      gamma  = xin(1)
      twsu   = xin(2)
      wmfls  = xin(3)
      twbamb = xin(4)
      patm   = xin(5)*1000.
      fansp  = xin(6)
      y1     = xin(7)
      tsump  = xin(8)

      twbamb = ftwb(toa, woa, patm) ! Calculate ambient wet-bulb temp. using
                                    ! data from common block of weather 12/18/02

      onof =(iostat(1)<-1).or.(iostat(1)==2)
      wtem =(iostat(2)<-1).or.(iostat(2)==2)
      wflo =(iostat(3)<-1).or.(iostat(3)==2)
      awbt =(iostat(4)<-1).or.(iostat(4)==2)
      totp =(iostat(5)<-1).or.(iostat(5)==2)
      fans =(iostat(6)<-1).or.(iostat(6)==2)
      valv =(iostat(7)<-1).or.(iostat(7)==2)
      frzde=onof.and.wtem.and.wflo.and.awbt.and.totp.and.fans.and.valv

! ***********************************
! * Check pump on and water flowing
! ***********************************

      if ((wmfls > 1.e-10).and.(gamma > 0)) then
                          ! water flow from chiller, calculate water flow
                          ! through tower
                          ! y1: actual valve position including crank angle and
                          ! hysteresis
      		y2=1.0-y1
          if (avalve>1.e-10) then
              aa=(1.0-avalve)/avalve
          else
      	      stop 'type143: valve authority = 0'
          endif
                            ! hydraulic resistances through valve and
                            ! conductances of each path
                            ! port 1 (through tower) is parabolic,
                            ! port 2 is linear
      		hr1=1.0/((1.0-cl)*y1*y1+cl)**2
      		hr2=1.0/((1.0-cl)*y2+cl)**2
      		ct=(aa+hr1)**(-0.5)
      		cb=(aa+hr2)**(-0.5)
                            ! flow through tower
      		wmfl=wmfls*ct/(ct+cb)
        
          if (wmfl > 1.e-10) then

! ************************************************
! * Calculation of the cooling tower performance *
! ************************************************

             ws  = fwphi(twsu,100.,patm)
             h   = fhsat(twsu,patm)
             wsb = fwphi(twbamb,100.,patm)
             hb  = fhsat(twbamb,patm)
             rho  = ro(patm,twbamb,wsb)
                          ! set effective air speed to be no less than 6.5%
                          ! of design value to account for wind driven heat
                          ! loss when fans are off
      	     avfl = avfld*max(fansp,0.065)
             amfl = avfl*rho
             wvfl = wmfl/rhow
             rlo  = rlod*(avfl/avfld)/(wvfl/wvfld)
             rlmind = h-hb-cpw*twbamb*(ws-wsb)
             if(abs(rlmind) < 1.0e-4) then
                return
             endif    
             rlmin = (cpw*(twsu-twbamb))/(h-hb-cpw*twbamb*(ws-wsb))
             rlamb = rlo/rlmin

!if(time > 349000.0) then
!   print *,'time, rlo, rlmin, rlamb:', time, rlo, rlmin, rlamb
!endif   
 
             if(rlamb > 6.0) then     ! limiting rlamb   5/21/07
                rlamb = 6.0
             elseif(rlamb < 0.0) then            
                rlamb = 0.0
             endif   
             effw  = ck*(1.-(exp(-rlamb)))
             twex = twsu-effw*(twsu-twbamb)
             twmix = (twex*wmfl + twsu*(wmfls-wmfl))/wmfls
             qload = wmfl*cpw*(twsu-twex)

             if (amfl > 0.01) then                 ! preventing divide zero 1/24/03
                effa = (qload/amfl)/(cpw*twsu-hb)
             end if

          else
                          ! flow is all through by-pass, set exit to supply
              twex = twbamb
              twmix = twsu
              qload = 0.0
              effw = 1.0
              effa = 0.0
          endif
      else
                          ! no flow from chiller, set outlet temperatures to
                          ! wet bulb
          twex = twbamb
          twmix = twbamb
          qload = 0.0
          effw = 1.0
          effa = 0.0
      endif

                          ! fan power
      if (fansp>0.0) then
          qfan = qfan0*(avfl/avfld)**3.
      else
          qfan = 0.0
      endif

                          ! Determine sump temperature, calculate derivative
                          ! if sump volume > 0
      if(vsump > 0.) then
        dtsumpdt = (twmix-tsump)*wmfls/msump
        yout(1)=dtsumpdt

        if ((abs(twmix-tsump) <= rtolx*abs(tsump)+atolx).and. frzde) then
      	    iostat(1)=1
        else
            iostat(1)=0
        endif
      else
        yout(1)=twmix
        iostat(1)=1
      endif

! *****************
! * Output values *
! *****************

      if(.not. cooling) then                      ! no cooling condition  1/27/03
         twmix = twbamb
         twex = twbamb
         wmfl = 0.0
         amfl = 0.0
         qload = 0.0
         qfan  = 0.0
         effw = 1.0
         effa = 0.0
      endif

      yout(2) = twmix
      yout(3) = twex
      yout(4) = wmfl
      yout(5) = amfl
      yout(6) = qload
      yout(7) = qfan
      yout(8) = effw
      yout(9) = effa

      do i=2,no
         iostat(i)=1
      enddo

      return
      end subroutine type143

! ********************************************************************

	subroutine type144(xin,yout,par_v,saved_v,iostat)

! ********************************************************************
!
! *************************
! *     CHILLER SUMP      *
! *************************
!
! ********************************************************************
! * inputs                                                           *
! * 1) inlet water temperature (C)                                   *
! * 2) supply water mass flow rate from chiller (C)                  *
! * 3) temperature of water exiting sump (C)                         *
! *                                                                  *
! * outputs                                                          *
! * 1) temperature of water exiting sump (C)                         *
! *
! * parameters                                                       *
! * 1) sump volume (m3)                                              *
! *                                                                  *
! **********************************************************************
!
!      created:  Dec. 17, 2002  Cheol Park, NIST
!      updated : January 25, 2007
!            Converted Fortran 77 code into Fortran 90
!
! **********************************************************************

      use modsim_head
      use cool_plant
      implicit none
      integer, parameter                :: ni=3, no=1, np=1, ns=1,&
                                           ndeq=1
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real            :: twsu,wmfls,tsump,vsump,msump,dtsumpdt
      real            :: rhow=1000.0

! ****************
! * input values *
! ****************

      twsu   = xin(1)
      wmfls  = xin(2)
      tsump  = xin(3)

      vsump  = par_v(1)
      msump  = vsump*rhow

                          ! Determine sump temperature, calculate derivative
                          ! if sump volume > 0
      if(msump>0.0) then
         dtsumpdt = (twsu-tsump)*wmfls/msump
      else
         print *,' Chiller sump volume is zero. '
         stop
      endif

! *****************
! * Output values *
! *****************
      yout(1)=dtsumpdt
      iostat(1)=0

      return
      end subroutine type144

! *******************************************************************
!
!      COOLING TOWER
!
! Modified by P. Haves, University of Oxford, 6/90, to include three
! port valve, variable speed fan, and sump.  Assumes wind induced
! flow gives 10% duty with fans off at design conditions.
!
! *******************************************************************
! * 1) Webb, R. L., "A unified theoretical treatment for thermal    *
! *    analysis of cooling towers, evaporative condensers and       *
! *    fluid coolers", ASHRAE Transactions, vol. 90, part 2, 1984.  *
! * 2) Webb, R. L., Villacres, A., "Algorithms for performance      *
! *    simulation of cooling towers, evaporative condensers and     *
! *    fluid coolers", ASHRAE Transactions, vol. 90, part 2, 1984.  *
! * 3) Modified and improved by J.L. Nizet, Laboratoire de          *
! *    Physique du Batiment et de Thermodynamique,                  *
! *    University of Liege, March 1985, and V.Giaretto, A. Mazza    *
! *    Dipartimento di Energetica del Politecnico di Torino, Feb.'87*
! * 4) Simplified method contributed by H.C.Pettsman, C.Kwakernaak  *
! *    and H.J.Nicolaas, TNO-Institute of Applied Phisic,           *
! *    June 1986.                                                   *
! *******************************************************************
!
! ********************************************************************
! *                                                                  *
! * inputs                                                           *
! * ======
! * 1) gamma:   on/off switch for pump (0/1)                         *
! * 2) twsu:    inlet water temperature (C)                          *
! * 3) wmfls:   supply water mass flow rate from chiller (C)         *
! * 4) twbamb:  ambient wet bulb temperature (C)                     *
! * 5) patm:    atmospheric pressure (Pa)                            *
! * 6) fansp:   fan speed control (fraction of design air flow)(0-1) *
! * 7) y1:      valve stem position (0-1)                            *
! *                                                                  *
! * outputs                                                          *
! * =======                                                          *
! * 1) tsump:   temperature of water exiting sump (C)                *
! * 2) twmix:   mixed water temperature returned to sump (C)         *
! * 3) twex:    tower outlet water temperature (C)                   *
! * 4) wmfl:    water mass flow rate through tower (kg/s)            *
! * 5) amfl:    air mass flow rate (kg/s)                            *
! * 6) qload:   heat transfer to the air (kW)                        *
! * 7) qfan:    fan electrical power consumption (kW)                *
! * 8) effw:    water side effectiveness                             *
! * 9) effa:    air side effectiveness                               *
! *                                                                  *
! * parameters                                                       *
! * ==========                                                       *
! * 1) twbambd: design supply air wet bulb temperature (C)           *
! * 2) twsud:   design supply water temperature (C)                  *
! * 3) twexd:   design exhaust water temperature (C)                 *
! * 4) avfld:   design air volume flow rate (m3/s)                   *
! * 5) wvfld:   design water volume flow rate (m3/s)                 *
! * 6) qfan0:   design fan power (kW)                                *
! * 7) avalve:  valve authority (parabolic/linear valve) (0-1)       *
! * 8) cl:      valve leakage (0-1)                                  *
! * 9) vsump:   sump volume (m3)                                     *
! ********************************************************************
!
!      revised:  Dec. 10, 2002  Cheol Park, NIST
!      updated : April 17, 2007
!            Converted Fortran 77 code into Fortran 90
!            Subroutine diffeq is used to solve the diff. eq. 
!
! ********************************************************************

    Subroutine type145(xin,yout,par_v,saved_v,iostat)

    use modsim_head
    use cool_plant
    implicit none
    integer, parameter                :: ni=7,no=9,np=9,&
                                         ndiffeq=1,ns=3+ndiffeq*2
    real,dimension(ni),intent(in)     :: xin
    real,dimension(no),intent(out)    :: yout
    real,dimension(np),intent(in)     :: par_v
    real,dimension(ns),intent(in out) :: saved_v
    integer,dimension(no),intent(out) :: iostat

    real              :: twexd,avfld,wvfld,qfan0,avalve,cl,vsump,&
                         msump,etactd,ws,h,wsb,hb,rho,amfld,wmfld,&
                         rlmind,rlod,rlambd,ck,gamma,twsu,wmfls,&
                         twbamb,fansp,y1,tsump,y2,aa,hr1,hr2,ct,&
                         cb,wmfl,avfl,amfl,wvfl,rlo,rlmin,rlamb,&
                         effw,twex,twmix,qload,effa,qfan,dtsumpdt,&
                         stop_fault,start_fault,w,ta,patm,&
                         ro,twbambd,twsud,fwphi,fhsat,ftwb
    real              :: tsumpp,aa1,bb1,tsumpmax,tsumpmin,tsumpbar
    logical           :: onof,wtem,wflo,awbt,totp,fans,valv,frzde
    real              :: patmd=101325.0,rhow=1000.0,cpw=4.187
    integer           :: i

    ro(patm,ta,w) = patm/(273.16+ta)/(287.06+461.52*w)

    if(time > start_chiller .and. time < stop_chiller) then  !  1/27/03
       cooling = .true.
    else
       cooling = .false.
    endif

!   Input values

    gamma  = xin(1)
    twsu   = xin(2)
    wmfls  = xin(3)
    twbamb = xin(4)
    patm   = xin(5)*1000.
    fansp  = xin(6)
    y1     = xin(7)

!   Parameters

    twbambd  = par_v(1)
    twsud    = par_v(2)
    twexd    = par_v(3)
    avfld    = par_v(4)
    wvfld    = par_v(5)
    qfan0    = par_v(6)
    avalve   = par_v(7)
    cl       = par_v(8)
    vsump    = par_v(9)

    msump  = vsump*rhow

!   Calculation of the rlo-value at design conditions

    if (itime == 1) then
                        ! first timestep, store ck for use at later
                        ! timesteps
       etactd = (twsud-twexd)/(twsud-twbambd)

       ws  = fwphi(twsud,100.,patmd)
       h   = fhsat(twsud,patmd)
       wsb = fwphi(twbambd,100.,patmd)
       hb  = fhsat(twbambd,patmd)

       rho   = ro(patmd,twbambd,wsb)
       amfld = avfld*rho
       wmfld = wvfld*998.

       rlmind = (cpw*(twsud-twbambd))/(h-hb-cpw*twbambd*(ws-wsb))
       rlod   = amfld/wmfld
       rlambd = rlod/rlmind

       ck = etactd/(1.-exp(-rlambd))

       saved_v(1)= -99999.0
       if (init==0) then
                        ! tsump is in equilibrium with environment
                        ! at start of simulation
          saved_v(2) = twsu
       endif
       saved_v(4)=rlod
       saved_v(5)=ck
    else
                        ! use value of ck calculated at first timestep
       rlod=saved_v(4)
       ck=saved_v(5)
    endif

    if (time>saved_v(1)) then
                        ! First call of timestep - update value of sump
                        ! temperature from previous time-step
       saved_v(3) = saved_v(2)
    endif
                        !   Update previous values
    tsumpp = saved_v(3)

    onof =(iostat(1)<-1).or.(iostat(1)==2)
    wtem =(iostat(2)<-1).or.(iostat(2)==2)
    wflo =(iostat(3)<-1).or.(iostat(3)==2)
    awbt =(iostat(4)<-1).or.(iostat(4)==2)
    totp =(iostat(5)<-1).or.(iostat(5)==2)
    fans =(iostat(6)<-1).or.(iostat(6)==2)
    valv =(iostat(7)<-1).or.(iostat(7)==2)
    frzde=onof.and.wtem.and.wflo.and.awbt.and.totp.and.fans.and.valv

!   Check pump on and water flowing

    if ((wmfls > 1.e-10).and.(gamma > 0)) then
                        ! water flow from chiller, calculate water flow
                        ! through tower
                        ! y1: actual valve position including crank angle and
                        ! hysteresis
       y2=1.0-y1
       if (avalve>1.e-10) then
           aa=(1.0-avalve)/avalve
       else
             stop 'type145: valve authority = 0'
       endif
                        ! hydraulic resistances through valve and
                        ! conductances of each path
                        ! port 1 (through tower) is parabolic,
                        ! port 2 is linear
       hr1=1.0/((1.0-cl)*y1*y1+cl)**2
       hr2=1.0/((1.0-cl)*y2+cl)**2               
       ct=(aa+hr1)**(-0.5)                       
       cb=(aa+hr2)**(-0.5)                       
                        ! flow through tower
       wmfl=wmfls*ct/(ct+cb)                     

       if (wmfl > 1.e-10) then

!   Calculation of the cooling tower performance

          ws  = fwphi(twsu,100.,patm)
          h   = fhsat(twsu,patm)
          wsb = fwphi(twbamb,100.,patm)
          hb  = fhsat(twbamb,patm)
          rho  = ro(patm,twbamb,wsb)
                        ! set effective air speed to be no less than 6.5%
                        ! of design value to account for wind driven heat
                        ! loss when fans are off
    	  avfl = avfld*max(fansp,0.065)
          amfl = avfl*rho
          wvfl = wmfl/rhow
          rlo  = rlod*(avfl/avfld)/(wvfl/wvfld)
          rlmin = (cpw*(twsu-twbamb))/(h-hb-cpw*twbamb*(ws-wsb))
          rlamb = rlo/rlmin
          effw  = ck*(1.-(exp(-rlamb)))
          twex = twsu-effw*(twsu-twbamb)
          twmix = (twex*wmfl + twsu*(wmfls-wmfl))/wmfls
          qload = wmfl*cpw*(twsu-twex)
                        ! preventing divide zero 1/24/03
          if (amfl > 0.01) then
             effa = (qload/amfl)/(cpw*twsu-hb)
          end if
       else
                        ! flow is all through by-pass, set exit to supply
          twex = twbamb
          twmix = twsu
          qload = 0.0
          effw = 1.0
          effa = 0.0
       endif
    else
                        ! no flow from chiller, set outlet temperatures to
                        ! wet bulb
       twex = twbamb
       twmix = twbamb
       qload = 0.0
       effw = 1.0
       effa = 0.0
    endif
                        ! fan power
    if (fansp>0.0) then
       qfan = qfan0*(avfl/avfld)**3.
    else
       qfan = 0.0
    endif

!   Determine sump temperature, calculate derivative, if sump volume > 0

    if(vsump > 0.) then
                        ! Integrate analytically using diffeq
       aa1 = -wmfls/msump
       bb1 = wmfls*twmix/msump                                                     
       call diffeq(time,aa1,bb1,tsumpp,tsump,tsumpbar)                             
                        ! Save current value to be used in updating at
                        ! next time step
       saved_v(2) = tsump                                                          
                        ! Save time of current call
       saved_v(1) = time                                                           
                                                                                   
!    Output
                                                                                   
       yout(1) = tsump

       if ((abs(twmix-tsump) <= rtolx*abs(tsump)+atolx).and. frzde) then           
       	  iostat(1)=1
       else                                                                        
          iostat(1)=0                                                              
       endif                                                                       
    else
       yout(1)=twmix
       iostat(1)=1
    endif

!   Output values
                        ! no cooling condition  1/27/03
    if(.not. cooling) then
       twmix = twbamb
       twex = twbamb
       wmfl = 0.0
       amfl = 0.0
       qload = 0.0
       qfan  = 0.0
       effw = 1.0
       effa = 0.0
    endif

    yout(2) = twmix
    yout(3) = twex
    yout(4) = wmfl
    yout(5) = amfl
    yout(6) = qload
    yout(7) = qfan
    yout(8) = effw
    yout(9) = effa

    do i=2,no
       iostat(i)=1
    enddo

    return
    end subroutine type145

! *********************************************************************
!
!  TYPE146 - Chiller sump temperature
!            Subroutine diffeq is used to solve the diff. eq.
!
!   April 12, 2007  Cheol Park, NIST
!
! *********************************************************************
!  inputs
!  ======
!   twsu    : inlet water temperature to chiller sump           (C)
!   wmfls      : water mass flow rate from cooling coil            (kg/s)
!
!  outputs
!  =======
!   tsump   : water temperaute exiting sump                     (C)
!
!  parameters
!  ==========
!   vsump   : sump volume                                       (m3)
!   tsumpmax: upper limit of water temperaute exiting sump      (C)
!   tsumpmin: lower limit of water temperaute exiting sump      (C)
!
! *********************************************************************

    subroutine type146(xin,yout,par_v,saved_v,iostat)

    use modsim_head
    use cool_plant
    implicit none
    integer,parameter                 :: ni=2,no=1,np=1,&
                                         ndiffeq=1,ns=1+ndiffeq*2
    real,dimension(ni),intent(in)     :: xin
    real,dimension(no),intent(out)    :: yout
    real,dimension(np),intent(in)     :: par_v
    real,dimension(ns),intent(in out) :: saved_v
    integer,dimension(no),intent(out) :: iostat

    real         :: twsu,wmfls,vsump,tsumpmax,tsumpmin,aa,bb,&
                    tsump,tsumpp,tsumpbar
    real         :: rhow=1000.0

    if(time > start_chiller .and. time < stop_chiller) then  !  1/27/03
       cooling = .true.
    else
       cooling = .false.
    endif

!   Read in input

    twsu  = xin(1)
    wmfls = xin(2)

!   Read in parameters

    vsump    = par_v(1)

!   Initialize at beginning of simulation

    if (itime<=1) then
       if (init==0 .or. saved_v(1)>time) then
          saved_v(1) = -999999.0
       endif
       if (init==0) then
!   tsump is in equilibrium with environment at start of simulation
            saved_v(2) = twsu
       endif
    endif

    if (time>saved_v(1)) then

!   First call of timestep - update value of sump temperature from
!   previous time-step

       saved_v(3) = saved_v(2)
    endif

!   Update previous values

    tsumpp = saved_v(3)
    if(vsump<=0.0) then
       tsump = twsu
    else

!   Integrate analytically

       aa = -wmfls/(rhow*vsump)
       bb = wmfls*twsu/(rhow*vsump)
       call diffeq(time,aa,bb,tsumpp,tsump,tsumpbar)
    endif

!   Save current value to be used in updating at next step time

    saved_v(2) = tsump

!   Save time of current call

    saved_v(1) = time

!   Output

    if(cooling) then
       yout(1)=tsump
    else
       yout(1)=twsu
    endif

!   Determine whether to allow freezing
!   Freezing of the output is allowed only if the input is a constant
!   or a boundary variable and the change in the output is small

    if (((iostat(1)<-1).or.(iostat(1)==2))&
       .and. ((abs(twsu - tsump))<=(rtolx * abs(tsump)+atolx))) then
       iostat(1) = 1
    else
       iostat(1) = 0
    endif

    return
    end subroutine type146

!******************************************************************

	  subroutine type179(xin,yout,par_v,saved_v,iostat)

! ********************************************************************
!
!  COOLING TOWER CONTROLLER for IEA Annex 17 exercise c1 (trend format)
!
!  assumes two speed fan
!---------------------------------------------------------------
!
!	input(1) = water temperature sensor (C)
!	input(2) = water temperature setpoint (C)
!	input(3) = status: (0=plant off, 1=on)
!
!	output(1)= valve demanded position (0-1)
!	output(2)= fan speed (0-1)
!
!	saved_v(1) = time of previous call
!	saved_v(2) = time of previous controller execution
!	saved_v(3) = time of previous sample
!	saved_v(4) = integral term from previous call
!	saved_v(5) = integral term from previous sample
!	saved_v(6) = output from previous call
!	saved_v(7) = output from previous sample
!	saved_v(8) = fan low hysteresis output from previous call
!	saved_v(9) = fan low hysteresis output from previous sample
!	saved_v(10) = fan high hysteresis output from previous call
!	saved_v(11) = fan high hysteresis output from previous sample
!	saved_v(12) = valve demanded position from previous call
!	saved_v(13) = fan speed from previous call
!
!       par_v(1) = proportional gain (%/C)
!	par_v(2) = integral time (sec)
!	par_v(3) = reschedule time (sec)
!       par_v(4) = number of times entered in sequence table
!       par_v(5) = controller number (parameter file='contn.par', n > 0)
!
! **********************************************************************
!
!      revised:  Dec. 10, 2002  Cheol Park, NIST
!      updated : January 25, 2007
!            Converted Fortran 77 code into Fortran 90
!
! **********************************************************************

      use modsim_head
      use cool_plant
      implicit none
      integer, parameter                :: ni=3, no=2, np=5, ns=13
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real       :: stop_fault,start_fault,propg,tint,rsec,t,&
                    tsp,tsamp,ti,pl1,dif,pid,fansp,valve,rsclto
      real       :: l1
      integer    :: status,fanlow,fanlowp,fanhigh,fanhighp,&
                    icont,is,ihyst,nseq

      if(time > start_chiller .and. time < stop_chiller) then
         cooling = .true.
      else
         cooling = .false.
      endif
                          !  parameters
      propg = par_v(1)
      tint  = par_v(2)
      rsec  = par_v(3)
      nseq  = nint(par_v(4))
      icont = nint(par_v(5))
                          ! inputs - Check for out of range
      t     = xin(1)
      tsp   = xin(2)
      status= nint(xin(3))
                          ! Initialize at beginning of simulation
      if (itime<=1) then
      	if (init==0 .or. saved_v(1)>time) then
      	   saved_v(1)=0.
      	   saved_v(2)=0.
        endif

      	if (init==0) then
      	   do is=4,10,2
      	      saved_v(is)=0.0
           enddo
      	   saved_v(8)=0.0
      	   saved_v(9)=0.0
        endif
      endif
                          ! Run controller if a sample instant
      tsamp=rsec/nseq

      if (time>=(saved_v(3)+tsamp)) then
                        ! first call of timestep - Update previous sample
                        ! instant values
          if (time>saved_v(1)) then
             do is=2,10,2
      	       saved_v(is+1)=saved_v(is)
             enddo
          endif
                        ! cooling demand (l1)
          ti=saved_v(5)
          pl1=saved_v(7)
          dif = 0.0                                                ! 11/19/02
          l1=pid(t,tsp,tsp,0.0,0,status,-propg,tint,0.0,&
                  ti,dif,pl1,0.,rsec,nseq)                        ! 11/19/02
          saved_v(4)=ti
          saved_v(6)=l1
                        ! fan off/low (fan on low at 50 rising, off at 10
                        ! falling)
          fanlowp=nint(saved_v(9))
          fanlow=ihyst(l1,30.0,20.0,fanlowp)
          saved_v(8)=float(fanlow)
                        ! fan low/high (fan on high at 90 rising, on low at
                        ! 50 falling)
          fanhighp=nint(saved_v(11))
          fanhigh=ihyst(l1,70.0,20.0,fanhighp)
          saved_v(10)=float(fanhigh)
                        ! fan speed (0=off, 0.5=low, 1.0=high)
          if (fanlow==0 .and. fanhigh==0) then
!      =======
!      Fan off
!      =======
	      fansp=0.0
                        ! valve (0-50 -> 0-100)
	      valve=rsclto(0.0,50.0,l1)
	      saved_v(12)=valve
	  elseif (fanlow==1 .and. fanhigh==0) then
!      =======
!      Fan low
!      =======
	      fansp=0.5
                        ! valve (0-90 -> 0-100)
	      valve=rsclto(0.0,90.0,l1)
	      saved_v(12)=valve
	  elseif (fanlow==1 .and. fanhigh==1) then
!      =======
!      Fan high
!      ========
      	      fansp=1.0
                        ! valve (0-90 -> 0-100)
              valve=rsclto(0.0,90.0,l1)
              saved_v(12)=valve
          else
              stop 'type179: illegal combination of fan speeds'
          endif
                        ! Save fan speed for non-sample-times
          saved_v(13)=fansp
                        ! Save time of controller execution
          saved_v(2)=time
      else
                        ! Not a sample instant, set output to
                        ! value from prev sample instant
          valve=saved_v(12)
          fansp=saved_v(13)
      endif

      saved_v(1)=time

      if(.not. cooling) then
          valve = 0.0
          fansp = 0.0
      endif

      yout(1)=valve/100.0
      yout(2)=fansp

      iostat(1)=1
      iostat(2)=1

      return
      end subroutine type179

! ******************************************************************

      subroutine type200(xin,yout,par_v,saved_v,iostat)            ! 8/30/02

! ******************************************************************
!     Mixing of water flows
! ******************************************************************
!
!      revised:  Dec. 10, 2002  Cheol Park, NIST
!      updated : January 25, 2007
!            Converted Fortran 77 code into Fortran 90
!
! **********************************************************************

      use modsim_head
      implicit none
      integer, parameter                :: maxmix = 12
      integer, parameter                :: ni=maxmix*2, no=2, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real, dimension(maxmix)           :: t, m
      real                              :: ma,to
      integer                           :: nmix,ind

      namelist /mixing/ nmix, to, ma               ! debug

! inputs
      nmix=nint(par_v(1))
      do ind=1,nmix
	t(ind) = xin(2*ind-1)
	m(ind) = xin(2*ind)
      end do
!  mass flow rate
      ma = 0.
      to = 0.
      do ind=1,nmix
         ma = ma + m(ind)
      end do
! humidity ratio and temperature
! avoid division by zero or overflow, but allow calculation of jacobian
      if (ma > 1.e-10) then
         do ind=1,nmix
            to = to + m(ind)*t(ind)
         end do
         to = to/ma
      else
! all flows zero, take arithmetic mean as reasonable starting point
         do ind=1,nmix
            to = to + t(ind)
         end do
            to = to/nmix
      endif
! outputs
!     if(nmix == 3)   print mixing                                ! debug
      yout(1) = to
      yout(2) = ma

      iostat(1)=1
      iostat(2)=1

      return
      end subroutine type200

! *****************************************************************************

	subroutine type201(xin,yout,par_v,saved_v,iostat)

! ********************************************************************
!
!       Adding electric powers
!
! **********************************************************************
!
!  inputs
!
!   1.  boiler total power			        (kW)
!   2.  chiller electric power			        (kW)
!   3.  cooling tower fan power			        (kW)
!   4.  ahu1 supply fan power		                (kW)
!   5.  ahu1 return fan power				(kW)
!   6.  zone electric power (zone 1)     		(kW)
!   7.  zone electric power (zone 2)     		(kW)
!   8.  zone electric power (zone 3)     		(kW)
!   9.  ahu2 supply fan power		                (kW)
!  10.  ahu2 return fan power			        (kW)
!  11.  zone electric power (zone 4)     		(kW)
!  12.  zone electric power (zone 5)     		(kW)
!  13.  zone electric power (zone 6)     		(kW)
!  14.  ahu3 supply fan power		                (kW)
!  15.  ahu3 return fan power				(kW)
!  16.  zone electric power (zone 7)     		(kW)
!  17.  zone electric power (zone 8)     		(kW)
!  18.  zone electric power (zone 9)     		(kW)
!
!  outputs
!
!   1.  boiler energy				        (kW)
!   2.  chiller & cooling tower energy                  (kW)
!   3.  ahu fan energy		                        (kW)
!   4.  zone electric consumption                       (kW)
!
!  pamrameter
!
!    1. (dummy)
!
! **********************************************************************
!
!      created:  August 30, 2002  Cheol Park, NIST
!      updated : January 25, 2007
!            Converted Fortran 77 code into Fortran 90
!
! **********************************************************************
!
      use modsim_head
!      implicit none
      integer, parameter                :: ni=18, no=4, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      namelist /electric_power/boiler, chiller, yout
                                                                             
      boiler           = xin(1)            ! input powers in kw              
      chiller          = xin(2)                                              
      towfan           = xin(3)                                              
      ahu1_sfan        = xin(4)                                              
      ahu1_rfan        = xin(5)                                              
      zone1_electric   = xin(6)                                              
      zone2_electric   = xin(7)                                              
      zone3_electric   = xin(8)                                              
      ahu2_sfan        = xin(9)                                              
      ahu2_rfan        = xin(10)                                             
      zone4_electric   = xin(11)                                             
      zone5_electric   = xin(12)                                             
      zone6_electric   = xin(13)                                             
      ahu3_sfan        = xin(14)                                             
      ahu3_rfan        = xin(15)                                             
      zone7_electric   = xin(16)                                             
      zone8_electric   = xin(17)                                             
      zone9_electric   = xin(18)                                             
                                                                             
      yout(1) = boiler                                                       
      yout(2) = chiller + towfan                                             
      yout(3) = ahu1_sfan + ahu1_rfan + ahu2_sfan&                           
             + ahu2_rfan + ahu3_sfan + ahu3_rfan                             
      yout(4) = zone1_electric + zone2_electric + zone3_electric&            
             + zone4_electric + zone5_electric + zone6_electric&             
             + zone7_electric + zone8_electric + zone9_electric              
                                                                             
      do i=1, no
          iostat(i)=1                                                        
      end do                                                                 
                                                                             
!       print electric_power
                                                                             
      return                                                                 
      end subroutine type201

! *****************************************************************************

	subroutine type202(xin,yout,par_v,saved_v,iostat)

! ********************************************************************
!
!       Integrator of thermal power for energy
!
! *****************************************************************************
!
!      created:   August 30, 2002  Cheol Park, NIST
!      updated : January 25, 2007
!            Converted Fortran 77 code into Fortran 90
!
! **********************************************************************
!
!  inputs
!
!   1.  ahu1 heating coil power        (kW)
!   2.  ahu1 cooling coil power        (kW)
!   3.  reheat coil power (zone 1)     (kW)
!   4.  reheat coil power (zone 2)     (kW)
!   5.  reheat coil power (zone 3)     (kW)
!   6.  ahu2 supply fan power	       (kW)
!   7.  ahu2 return fan power	       (kW)
!   8.  reheat coil power (zone 4)     (kW)
!   9.  reheat coil power (zone 5)     (kW)
!  10.  reheat coil power (zone 6)     (kW)
!  11.  ahu3 supply fan power	       (kW)
!  12.  ahu3 return fan power	       (kW)
!  13.  reheat coil power (zone 7)     (kW)
!  14.  reheat coil power (zone 8)     (kW)
!  15.  reheat coil power (zone 9)     (kW)
!
!  outputs
!
!   1. heating energy transfered from coils to air (kJ)
!   2. cooling energy transfered from coils to air (kJ)
!   3. total energy transfered from coils to air   (kJ)
!
!  pamrameter
!
!   1. (dummy)
!
! *****************************************************************************

      use modsim_head
!!      implicit none
      integer, parameter                :: ni=15, no=3, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      ahu1_hc      = xin(1)
      ahu1_cc      = xin(2)                                  
      zone1_rhc    = xin(3)                                  
      zone2_rhc    = xin(4)                                  
      zone3_rhc    = xin(5)                                  
      ahu2_hc      = xin(6)                                  
      ahu2_cc      = xin(7)                                  
      zone4_rhc    = xin(8)                                  
      zone5_rhc    = xin(9)                                  
      zone6_rhc    = xin(10)                                 
      ahu3_hc      = xin(11)                                 
      ahu3_cc      = xin(12)                                 
      zone7_rhc    = xin(13)                                 
      zone8_rhc    = xin(14)                                 
      zone9_rhc    = xin(15)                                 
                                                             
      dahu1_hc     = ahu1_hc/1000.0                          
      dahu1_cc     = ahu1_cc/1000.0                          
      dzone1_rhc   = zone1_rhc/1000.0                        
      dzone2_rhc   = zone2_rhc/1000.0                        
      dzone3_rhc   = zone3_rhc/1000.0                        
      dahu2_hc     = ahu2_hc/1000.0                          
      dahu2_cc     = ahu2_cc/1000.0                          
      dzone4_rhc   = zone4_rhc/1000.0                        
      dzone5_rhc   = zone5_rhc/1000.0                        
      dzone6_rhc   = zone6_rhc/1000.0                        
      dahu3_hc     = ahu3_hc/1000.0                          
      dahu3_cc     = ahu3_cc/1000.0                          
      dzone7_rhc   = zone7_rhc/1000.0                        
      dzone8_rhc   = zone8_rhc/1000.0                        
      dzone9_rhc   = zone9_rhc/1000.0                        
                                                             
                                                             
      ahu_hc_all  = dahu1_hc + dahu2_hc + dahu3_hc           
      ahu_cc_all  = dahu1_cc + dahu2_cc + dahu3_cc           
      rhc_all = dzone1_rhc + dzone2_rhc + dzone3_rhc&        
              + dzone4_rhc + dzone5_rhc + dzone6_rhc&        
              + dzone7_rhc + dzone8_rhc + dzone9_rhc         
                                                             
      yout(1) =  ahu_hc_all + rhc_all                        
      yout(2) =  ahu_cc_all                                  
      yout(3) =  ahu_hc_all + rhc_all - ahu_cc_all           
                                                             
      do i=1, no                                             
          iostat(i)=1                                        
      end do                                                 
                                                             
      return                                                 
      end subroutine type202

! *****************************************************************************

	subroutine type203(xin,yout,par_v,saved_v,iostat)

! ********************************************************************
!
!       Relative humidity and wet-bulb temprature
!
! **********************************************************************
!
!       created:   November 13, 2002 Cheol Park, NIST
!       Updated : January 25, 2007  Cheol Park, NIST
!            Converted Fortran 77 code into Fortran 90
!
! **********************************************************************
!
!  inputs
!
!   1. dry-bulb temperature           (C)
!   2. humidity ratio                 (-)
!   3. atmospheric pressure           (kPa)
!
!  outputs
!
!   1. relative humidity              (%)
!   2. wet-bulb temperature           (C)
!   3. dew-point temperature          (C)
!
!  pamrameter
!
!   1. (dummy)
!
! *****************************************************************************

      use modsim_head
      implicit none
      integer, parameter                :: ni=3, no=3, np=1, ns=0
      real,dimension(ni),intent(in)     :: xin
      real,dimension(no),intent(out)    :: yout
      real,dimension(np),intent(in)     :: par_v
      real,dimension(ns),intent(in out) :: saved_v
      integer,dimension(no),intent(out) :: iostat

      real                              :: tdb,w,patm,rh,fphi,twb,&
                                           ftwb,dp,ftdew
      integer                           :: i

      tdb   = xin(1)
      w     = xin(2)                                    
      patm  = xin(3)* 1000.0                            

      rh  = fphi(tdb, w, patm)                          
      twb = ftwb(tdb, w, patm)                          
      dp  = ftdew(w, patm)                              
                                                        
      yout(1) = rh
      yout(2) = twb                                     
      yout(3) = dp                                      
                                                        
      do i=1, no                                        
          iostat(i)=1                                   
      end do                                            
                                                        
      return                                            
      end subroutine type203

