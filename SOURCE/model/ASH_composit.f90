! **********************************************************************
! Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:  Dynamic or steady state coil and two port valve
! *              - calculates water-side inlet pressure from flow rate
! *
! * PURPOSE:     Calculate the outlet air and water conditions for a
! *              heating or cooling coil
! **********************************************************************
! * INPUTS
! * ======
! *  1. tai     : inlet air dry bulb temperature                     (C)
! *  2. gi      : inlet air humidity ratio                       (kg/kg)
! *  3. pao     : outlet air pressure                              (kPa)
! *  4. wa      : dry air mass flow rate                          (kg/s)
! *  5. twi     : inlet water temperature                            (C)
! *  6. ww      : water mass flow rate                            (kg/s)
! *  7. pwo     : outlet water pressure                            (kPa)
! *  8. vstem   : valve stem position                                (-)
! *  9. tsdyn   : effective coil surface temperature                 (C)
! *
! * OUTPUTS
! * =======
! *  1. ts      : effective coil surface temperature                 (C)
! *  2. tao     : outlet air dry bulb temperature                    (C)
! *  3. go      : outlet air humidity ratio                      (kg/kg)
! *  4. pai     : inlet air pressure                               (kPa)
! *  5. two     : outlet water temperature                           (C)
! *  6. pwi     : inlet water pressure                             (kPa)
! *  7. qa      : total heat transfer to the air                    (kW)
! *  8. shr     : sensible heat ratio                                (-)
! *  9. effect  : coil effectiveness                                 (-)
! * 10. bf      : coil by-pass factor                                (-)
! * 11. ho      : outlet air specific enthalpy                   (kJ/kg)
! * 12. rho     : outlet air relative humidity                       (-)
! * 13. twbo    : outlet air wet-bulb temperature                    (C)
! *
! * PARAMETERS
! * ==========
! *  1. dynamic : 0 = steady state, 1 = dynamic                      (-)
! *  2. ifault  : 0 = no faults, 1 = parallel flow (cooling coils)   (-)
! *  3. psychro : FALSE = no psychrometric output calcs, TRUE = calcs(-)
! *  4. nrow    : number of rows of tubes                            (-)
! *  5. ntpr    : number of tubes per row                            (-)
! *  6. ncir    : number of parallel water circuits                  (-)
! *  7. lcoil   : length of finned section in direction of flow      (m)
! *  8. hcoil   : heigth of finned section                           (m)
! *  9. wcoil   : width of finned section                            (m)
! * 10. dotube  : tube outside diameter                              (m)
! * 11. thitube : tube wall thickness                                (m)
! * 12. watube  : tube material (Al=1,Cu=2,Fe=3,CaCO3=4)             (-)
! * 13. spafin  : fin spacing (pitch)                                (m)
! * 14. thifin  : fin thickness                                      (m)
! * 15. wafin   : fin material (Al=1,Cu=2,Fe=3)                      (-)
! * 16. fra     : flow resistance on air side               (0.001 kg.m)
! * 17. frwcoil : coil water flow resistance                (0.001 kg.m)
! * 18. kv      : valve capacity (Kv)                    (m3/hr @ 1 bar)
! * 19. eqpchar : valve curvature parameter (0=linear)               (-)
! * 20. sv      : valve rangability                                  (-)
! * 21. cl      : valve leakage (fractional flow)                    (-)
! *
! * SAVED (dynamic mode only)
! * =====
! *  1. time    : time of previous call
! *  2. dt      : solution of differential equation from previous call
! *  3. dtp     : solution of DE from previous step time
! *  4. aa      : A coefficent in dT/dt=A*T+B from previous call
! *  5. aap     : A coefficent in dT/dt=A*T+B from previous step time
! *  6. bb      : B coefficent in dT/dt=A*T+B from previous call
! *  7. bbp     : B coefficent in dT/dt=A*T+B from previous step time
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  Assumes coil is all dry or all wet
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 14, 1995
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  LICOILDY,LICOILSS
!   FUNCTIONS  CALLED:   RLINPORT,REQPPORT
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           Cooling Coil Models to be used in Transient
!                        and/or Wet Regimes -- Theorical Analysis and
!                        Experimental Validation.
!                        X.DING, J-P EPPE,J.LEBRUN,M.WASACZ.
!                        I.E.A. ANNEX 17 document AN17-901019-01.
!                        University of Liege, Belgium.
!
!                        (adapted from the cooling coil model written
!                        by V. Corrado - Polytechnic of Turin, Italy)
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! *
! * Material Properties
! * -------------------
! * VOLCPM  : volumic specific heat of Al, Cu and mild steel     (kJ/m3)
! * CONMA   : thermal conductivity of Al, Cu and mild steel   (kW/(m.K))
! *
! * Geometrical characteristics
! * ---------------------------
! * fracfin : fraction of face area blocked by fins                  (-)
! * acgross : gross face area                                       (m2)
! * ntube   : total number of tubes                                  (-)
! * ttubelen: total tube length                                      (m)
! * axo     : gross outside tube surface area                       (m2)
! * ap      : net outside tube surface area, allowing for fins      (m2)
! * as      : net air side fin area                                 (m2)
! * ditube  : tube inside diameter                                   (m)
! * ai      : inside surface area of tubes                          (m2)
! * afree   : exchanger minimum free-flow area on air side          (m2)
! * aext    : exchanger total heat transfer area on air side        (m2)
! * aratio  : ratio of total fin area to total heat transfer area    (-)
! * bratio  : ratio of air side total area : water side internal area(-)
! * rifin   : fin inside radius                                      (m)
! * rofin   : effective fin radius                                   (m)
! * heifin  : effective fin height                                   (m)
! * aflow   : flow area on water side                               (m2)
! * hydiam  : coil hydraulic diameter                                (m)
! * confin  : thermal conductivity of fin material            (kW/(m.K))
! * contube : thermal conductivity of tube material           (kW/(m.K))
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type521(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndtdt=0,ndiffeq=1,&
                                             ni=9,no=13,np=21,&
                                             ns=1+ndiffeq*6
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: lcoil,kv
        integer      :: matube,mafin
        logical      :: dynamic,psychro
! **** Thermophysical constants
        real,dimension(4)   :: volcpm=(/2455.,3475.,3915.,2250./),&
                               conma=(/0.221,0.393,0.0453,0.0013/)
! **** Miscellaneous constants (pi, atmospheric pressure, minimum
! **** water-side resistance
        real         :: pi=3.14159,patmpa=101325.,frsmall=1.0e-8
! **** Break-points in valve characteristics (see documentation for
! **** type329 two port valve)
        real         :: xlin=0.01,xeqp=0.33
        real         :: eqpchar,sv,cl,fracfin,acgross,ttubelen,axo,asx,&
                        ap,ditube,ai,afree,aext,aratio,bratio,rifin,&
                        rofin,heifin,aflow,hydiam,confin,contube,vcpfin,&
                        vcptube,frwvalve,frwtot,pai,pwi,tsp,aap,bbp,&
                        twbo,rho,ho,bf,shr,qa,two,go,tao,bb,aa,aabar,&
                        bbbar,tsbar,ts,effect,tai,gi,gsat,fwphi,pao,&
                        wa,twi,ww,pwo,vstem,tsdyn,hcoil,wcoil,&
                        dotube,thitube,spafin,thifin,fra,frwcoil,&
                        rlinport,reqpport,dpqudlin
        integer      :: i,icountfl,ntube,is,ifault,nrow,ntpr,ncir

! **** Read in inputs
        tai     = xin(1)
        gi      = xin(2)
! **** Limit humidity ratio to saturated value
        gsat    = fwphi(tai,100.0,patmpa)
        gi      = min(gi,gsat)
        pao     = xin(3)
        wa      = xin(4)
        twi     = xin(5)
        ww      = xin(6)
        pwo     = xin(7)
        vstem   = xin(8)
        tsdyn   = xin(9)
! **** Read in parameters
! **** Dynamic or steady state
        if (nint(par_v(1))==1) then
            dynamic = .true.
        else
            dynamic = .false.
        endif
! **** Fault type
        ifault = nint(par_v(2))
! **** Peform additional pyschrometric calculations?
        if (nint(par_v(3))/=0) then
            psychro = .true.
        else
            psychro = .false.
        endif
! **** Coil geometry and materials
        nrow    = nint(par_v(4))
        ntpr    = nint(par_v(5))
        ncir    = nint(par_v(6))
        lcoil   = par_v(7)
        hcoil   = par_v(8)
        wcoil   = par_v(9)
        dotube  = par_v(10)
        thitube = par_v(11)
        matube  = nint(par_v(12))
        spafin  = par_v(13)
        thifin  = par_v(14)
        mafin   = nint(par_v(15))
! **** Air and water flow resistances
        fra     = par_v(16)
        frwcoil = par_v(17)
! **** Valve characteristics
        kv      = par_v(18)
        eqpchar = par_v(19)
        sv      = par_v(20)
        cl      = par_v(21)
! **** Heat exchanger geometry
        if (nrow<=2) then
! **** One or two rows -> cross flow
            icountfl=0
        elseif (nrow>2) then
            if (ifault/=1) then
! **** Three or more rows and no connection fault -> counterflow
                icountfl=1
            else
! **** Three or more rows and connection fault -> parallel flow
                icountfl=-1
            endif
        endif
! **** Heat transfer areas
        fracfin = thifin/spafin
        acgross = hcoil*wcoil
        ntube   = nrow*ntpr
        ttubelen= ntube*wcoil
        axo     = pi*dotube*dotube/4.0
        ap      = pi*dotube*ttubelen*(1.0-fracfin)
        asx      = 2.0*(hcoil*lcoil-ntube*axo)*wcoil/spafin
        ditube  = dotube-2.0*thitube
        ai      = pi*ditube*ttubelen
        afree   = acgross*(1.0-ntpr*dotube/hcoil)*(1.0-fracfin)
        aext    = asx+ap
        aratio  = asx/aext
        bratio  = aext/ai
        rifin   = dotube/2.0
        rofin   = sqrt(hcoil*lcoil/(pi*ntube))
        heifin  = rofin-rifin
        aflow   = ncir*pi*ditube*ditube/4.0
        hydiam  = 4.0*afree*lcoil/aext
        confin  = conma(mafin)
        contube = conma(matube)
        vcpfin  = volcpm(mafin)
        vcptube = volcpm(matube)
! **** Valve type and resistance
! **** Limit stem position to range 0-1
        vstem   = max( 0.0, min(1.0,vstem) )
! **** Select valve type and calculate resistance
        if (eqpchar<1.0) then
! **** Linear characteristic - two segments (linear and close-off)
            frwvalve=rlinport(vstem,kv,sv,xlin,cl)
        else
! **** Equal percentage characteristic - three segments (exponential,
! **** linear and close-off)
            frwvalve=reqpport(vstem,kv,xeqp,eqpchar,xlin,sv,cl)
        endif
! **** Total water-side resistance - coil + valve
        frwtot  = frwcoil+frwvalve
        if (frwtot<frsmall) then
! **** Total resistance almost zero
            write(*,*)&
            'type521: water-side flow resistance must not be < ',frsmall
            stop
        endif
! **** Calculate air and water inlet pressures from flow resistances
        pai     = pao+dpqudlin(fra,wa)
        pwi     = pwo+dpqudlin(frwtot,ww)
! **** Calculate coil duty and outlet conditions
        if (dynamic) then
! **** Initialize at beginning of simulation
            if (itime<=1) then
                if (init==0 .or. saved_v(1)>time) then
                    saved_v(1) = -999999.
               endif
                if (init==0) then
! **** Capacitance of coil initially at inlet air temperature (coil "off")
                    saved_v(2) = tai
                    saved_v(4) = 0.0
                    saved_v(6) = 0.0
                endif
            endif
! **** Every time-step
            if (time>saved_v(1)) then
! **** First call of time-step - update value of temperature rise, aa and bb
! **** from previous time-step
                do is = 2,ns-1,2
                    saved_v(is+1) = saved_v(is)
                enddo
            endif
! **** Update previous values
            tsp = saved_v(3)
            aap = saved_v(5)
            bbp = saved_v(7)
            call licoilab(tai,gi,wa,twi,ww,tsdyn,&
                          psychro,icountfl,lcoil,wcoil,dotube,thitube,&
                          thifin,ttubelen,asx,ditube,afree,aext,&
                          aratio,bratio,rifin,rofin,heifin,aflow,&
                          hydiam,confin,contube,vcpfin,vcptube,&
                          aa,bb,tao,go,two,qa,shr,bf,ho,rho,twbo)
! **** Integrate differential equation for surface temperature
            aabar = (aa + aap)/2.0
            bbbar = (bb + bbp)/2.0
            call diffeq(time,aabar,bbbar,tsp,ts,tsbar)
! **** Save time of current call
            saved_v(1) = time
! **** Save provisional value to be used in updating at next step time
            saved_v(2) = ts
            saved_v(4) = aa
            saved_v(6) = bb
        else
! **** Steady state model - output is value of surface temperature
            call licoilss(tai,gi,wa,twi,ww,&
                          psychro,icountfl,lcoil,wcoil,thitube,&
                          thifin,ditube,afree,aext,&
                          aratio,bratio,rifin,rofin,heifin,aflow,&
                          hydiam,confin,contube,&
                          ts,tao,go,two,qa,shr,effect,bf,ho,rho,twbo)
        endif
! **** Assign output values
        yout(1)  = ts
        yout(2)  = tao
        yout(3)  = go
        yout(4)  = pai
        yout(5)  = two
        yout(6)  = pwi
        yout(7)  = qa
        yout(8)  = shr
        if (dynamic) then
! **** Effectiveness not defined for dynamic case
            yout(9)  = -99.9
        else
            yout(9)  = effect
        endif
        yout(10) = bf
        if (psychro) then
! **** Auxilliary outputs
            yout(11) = ho
            yout(12) = rho/100.
            yout(13) = twbo
        else
! **** Auxilliary outputs not calculated - set values to flag this
            yout(11) = -99.9
            yout(12) = -99.9
            yout(13) = -99.9
        endif
! **** Allow freezing of algebraic variables
        if (.not.dynamic) then
            iostat(1)=1
        else
            iostat(1)=0
        endif
        do i=2,no
            iostat(i)=1
        enddo
        return
        end subroutine type521

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:  Dynamic or steady state coil and two port valve
! *              - calculates water-side flow rate from inlet pressure
! *
! * PURPOSE:     Calculate the outlet air and water conditions for a
! *              heating or cooling coil
! **********************************************************************
! * INPUTS
! * ======
! *  1. tai     : inlet air dry bulb temperature                     (C)
! *  2. gi      : inlet air humidity ratio                       (kg/kg)
! *  3. pao     : outlet air pressure                              (kPa)
! *  4. wa      : dry air mass flow rate                          (kg/s)
! *  5. twi     : inlet water temperature                            (C)
! *  6. pwi     : inlet water pressure                             (kPa)
! *  7. pwo     : outlet water pressure                            (kPa)
! *  8. vstem   : valve stem position                                (-)
! *  9. tsdyn   : effective coil surface temperature                 (C)
! *
! * OUTPUTS
! * =======
! *  1. ts      : effective coil surface temperature                 (C)
! *  2. tao     : outlet air dry bulb temperature                    (C)
! *  3. go      : outlet air humidity ratio                      (kg/kg)
! *  4. pai     : inlet air pressure                               (kPa)
! *  5. two     : outlet water temperature                           (C)
! *  6. ww      : water mass flow rate                            (kg/s)
! *  7. qa      : total heat transfer to the air                    (kW)
! *  8. shr     : sensible heat ratio                                (-)
! *  9. effect  : coil effectiveness                                 (-)
! * 10. bf      : coil by-pass factor                                (-)
! * 11. ho      : outlet air specific enthalpy                   (kJ/kg)
! * 12. rho     : outlet air relative humidity                       (-)
! * 13. twbo    : outlet air wet-bulb temperature                    (C)
! *
! * PARAMETERS
! * ==========
! *  1. dynamic : 0 = steady state, 1 = dynamic                      (-)
! *  2. ifault  : 0 = no faults, 1 = parallel flow (cooling coils)   (-)
! *  3. psychro : FALSE = no psychrometric output calcs, TRUE = calcs(-)
! *  4. nrow    : number of rows of tubes                            (-)
! *  5. ntpr    : number of tubes per row                            (-)
! *  6. ncir    : number of parallel water circuits                  (-)
! *  7. lcoil   : length of finned section in direction of flow      (m)
! *  8. hcoil   : heigth of finned section                           (m)
! *  9. wcoil   : width of finned section                            (m)
! * 10. dotube  : tube outside diameter                              (m)
! * 11. thitube : tube wall thickness                                (m)
! * 12. watube  : tube material (Al=1,Cu=2,Fe=3,CaCO3=4)             (-)
! * 13. spafin  : fin spacing (pitch)                                (m)
! * 14. thifin  : fin thickness                                      (m)
! * 15. wafin   : fin material (Al=1,Cu=2,Fe=3)                      (-)
! * 16. fra     : flow resistance on air side               (0.001 kg.m)
! * 17. frwcoil : coil water flow resistance                (0.001 kg.m)
! * 18. kv      : valve capacity (Kv)                    (m3/hr @ 1 bar)
! * 19. eqpchar : valve curvature parameter (0=linear)               (-)
! * 20. sv      : valve rangability                                  (-)
! * 21. cl      : valve leakage (fractional flow)                    (-)
! *
! * SAVED (dynamic mode only)
! * =====
! *  1. time    : time of previous call
! *  2. dt      : solution of differential equation from previous call
! *  3. dtp     : solution of DE from previous step time
! *  4. aa      : A coefficent in dT/dt=A*T+B from previous call
! *  5. aap     : A coefficent in dT/dt=A*T+B from previous step time
! *  6. bb      : B coefficent in dT/dt=A*T+B from previous call
! *  7. bbp     : B coefficent in dT/dt=A*T+B from previous step time
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  Assumes coil is all dry or all wet
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 14, 1995
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  LICOILDY,LICOILSS
!   FUNCTIONS  CALLED:   RLINPORT,REQPPORT
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           Cooling Coil Models to be used in Transient
!                        and/or Wet Regimes -- Theorical Analysis and
!                        Experimental Validation.
!                        X.DING, J-P EPPE,J.LEBRUN,M.WASACZ.
!                        I.E.A. ANNEX 17 document AN17-901019-01.
!                        University of Liege, Belgium.
!
!                        (adapted from the cooling coil model written
!                        by V. Corrado - Polytechnic of Turin, Italy)
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! *
! * Material Properties
! * -------------------
! * volcpm  : volumic specific heat of Al, Cu and mild steel     (kJ/m3)
! * conma   : thermal conductivity of Al, Cu and mild steel   (kW/(m.K))
! *
! * Geometrical characteristics
! * ---------------------------
! * fracfin : fraction of face area blocked by fins                  (-)
! * acgross : gross face area                                       (m2)
! * ntube   : total number of tubes                                  (-)
! * ttubelen: total tube length                                      (m)
! * axo     : gross outside tube surface area                       (m2)
! * ap      : net outside tube surface area, allowing for fins      (m2)
! * as      : net air side fin area                                 (m2)
! * ditube  : tube inside diameter                                   (m)
! * ai      : inside surface area of tubes                          (m2)
! * afree   : exchanger minimum free-flow area on air side          (m2)
! * aext    : exchanger total heat transfer area on air side        (m2)
! * aratio  : ratio of total fin area to total heat transfer area    (-)
! * bratio  : ratio of air side total area : water side internal area(-)
! * rifin   : fin inside radius                                      (m)
! * rofin   : effective fin radius                                   (m)
! * heifin  : effective fin height                                   (m)
! * aflow   : flow area on water side                               (m2)
! * hydiam  : coil hydraulic diameter                                (m)
! * confin  : thermal conductivity of fin material            (kW/(m.K))
! * contube : thermal conductivity of tube material           (kW/(m.K))
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type522(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndtdt=0,ndiffeq=1,&
                                             ni=9,no=13,np=21,&
                                             ns=1+ndiffeq*6
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: lcoil,kv
        integer      :: matube,mafin
        logical      :: dynamic,psychro
! **** Thermophysical constants
        real,dimension(4)   :: volcpm=(/2455.,3475.,3915.,2250./),&
                               conma=(/0.221,0.393,0.0453,0.0013/)
! **** Miscellaneous constants (pi, atmospheric pressure, minimum
! **** water-side resistance
        real         :: pi=3.14159,patmpa=101325.,frsmall=1.0e-8
! **** Break-points in valve characteristics (see documentation for
! **** type329 two port valve)
        real         :: xlin=0.01,xeqp=0.33
        real         :: eqpchar,sv,cl,fracfin,acgross,ttubelen,axo,asx,&
                        ap,ditube,ai,afree,aext,aratio,bratio,rifin,&
                        rofin,heifin,aflow,hydiam,confin,contube,vcpfin,&
                        vcptube,frwvalve,frwtot,pai,pwi,tsp,aap,bbp,&
                        twbo,rho,ho,bf,shr,qa,two,go,tao,bb,aa,aabar,&
                        bbbar,tsbar,ts,effect,tai,gi,gsat,fwphi,pao,&
                        wa,twi,ww,pwo,vstem,tsdyn,hcoil,wcoil,&
                        dotube,thitube,spafin,thifin,fra,frwcoil,&
                        rlinport,reqpport,dpqudlin,dpw,wqudlin
        integer      :: i,icountfl,ntube,is,ifault,nrow,ntpr,ncir

 ! **** Read in inputs
        tai     = xin(1)
        gi      = xin(2)
! **** Limit humidity ratio to saturated value
        gsat    = fwphi(tai,100.0,patmpa)
        gi      = min(gi,gsat)
        pao     = xin(3)
        wa      = xin(4)
        twi     = xin(5)
        pwi     = xin(6)
        pwo     = xin(7)
        vstem   = xin(8)
        tsdyn   = xin(9)
! **** Read in parameters
! **** Dynamic or steady state
        if (nint(par_v(1))==1) then
            dynamic = .true.
        else
            dynamic = .false.
        endif
! **** Fault type
        ifault = nint(par_v(2))
! **** Peform additional pyschrometric calculations?
        if (nint(par_v(3))/=0) then
            psychro = .true.
        else
            psychro = .false.
        endif
! **** Coil geometry and materials
        nrow    = nint(par_v(4))
        ntpr    = nint(par_v(5))
        ncir    = nint(par_v(6))
        lcoil   = par_v(7)
        hcoil   = par_v(8)
        wcoil   = par_v(9)
        dotube  = par_v(10)
        thitube = par_v(11)
        matube  = nint(par_v(12))
        spafin  = par_v(13)
        thifin  = par_v(14)
        mafin   = nint(par_v(15))
! **** Air and water flow resistances
        fra     = par_v(16)
        frwcoil = par_v(17)
! **** Valve characteristics
        kv      = par_v(18)
        eqpchar = par_v(19)
        sv      = par_v(20)
        cl      = par_v(21)
! **** Heat exchanger geometry
        if (nrow<=2) then
! **** One or two rows -> cross flow
            icountfl=0
        elseif (nrow>2) then
            if (ifault/=1) then
! **** Three or more rows and no connection fault -> counterflow
                icountfl=1
            else
! **** Three or more rows and connection fault -> parallel flow
                icountfl=-1
            endif
        endif
! **** Heat transfer areas
        fracfin = thifin/spafin
        acgross = hcoil*wcoil
        ntube   = nrow*ntpr
        ttubelen= ntube*wcoil
        axo     = pi*dotube*dotube/4.0
        ap      = pi*dotube*ttubelen*(1.0-fracfin)
        asx     = 2.0*(hcoil*lcoil-ntube*axo)*wcoil/spafin
        ditube  = dotube-2.0*thitube
        ai      = pi*ditube*ttubelen
        afree   = acgross*(1.0-ntpr*dotube/hcoil)*(1.0-fracfin)
        aext    = asx+ap
        aratio  = asx/aext
        bratio  = aext/ai
        rifin   = dotube/2.0
        rofin   = sqrt(hcoil*lcoil/(pi*ntube))
        heifin  = rofin-rifin
        aflow   = ncir*pi*ditube*ditube/4.0
        hydiam  = 4.0*afree*lcoil/aext
        confin  = conma(mafin)
        contube = conma(matube)
        vcpfin  = volcpm(mafin)
        vcptube = volcpm(matube)
! **** Valve type and resistance
! **** Limit stem position to range 0-1
        vstem   = max( 0.0, min(1.0,vstem) )
! **** Select valve type and calculate resistance
        if (eqpchar<1.0) then
! **** Linear characteristic - two segments (linear and close-off)
            frwvalve = rlinport(vstem,kv,sv,xlin,cl)
        else
! **** Equal percentage characteristic - three segments (exponential,
! **** linear and close-off)
            frwvalve = reqpport(vstem,kv,xeqp,eqpchar,xlin,sv,cl)
        endif
! **** Total water-side resistance - coil + valve
        frwtot  = frwcoil+frwvalve
        if (frwtot<frsmall) then
! **** Total resistance almost zero
            write(*,*)&
            'type522: water-side flow resistance must not be < ',frsmall
            stop
        endif
! **** Calculate air and water inlet pressures from flow resistances
        pai = pao+dpqudlin(fra,wa)
        dpw = pwi-pwo
        ww  = wqudlin(frwtot,dpw)
! **** Calculate coil duty and outlet conditions
        if (dynamic) then
! **** Initialize at beginning of simulation
            if (itime<=1) then
                if (init==0 .or. saved_v(1)>time) then
                    saved_v(1) = -999999.
               endif
                if (init==0) then
! **** Capacitance of coil initially at inlet air temperature (coil "off")
                    saved_v(2) = tai
                    saved_v(4) = 0.0
                    saved_v(6) = 0.0
                endif
            endif
! **** Every time-step
            if (time>saved_v(1)) then
! **** First call of time-step - update value of temperature rise, aa and bb
! **** from previous time-step
                do is = 2,ns-1,2
                    saved_v(is+1) = saved_v(is)
                enddo
            endif
! **** Update previous values
            tsp = saved_v(3)
            aap = saved_v(5)
            bbp = saved_v(7)
            call licoilab(tai,gi,wa,twi,ww,tsdyn,&
                          psychro,icountfl,lcoil,wcoil,dotube,thitube,&
                          thifin,ttubelen,asx,ditube,afree,aext,&
                          aratio,bratio,rifin,rofin,heifin,aflow,&
                          hydiam,confin,contube,vcpfin,vcptube,&
                          aa,bb,tao,go,two,qa,shr,bf,ho,rho,twbo)
! **** Integrate differential equation for surface temperature
            aabar = (aa + aap)/2.0
            bbbar = (bb + bbp)/2.0
            call diffeq(time,aabar,bbbar,tsp,ts,tsbar)
! **** Save time of current call
            saved_v(1) = time
! **** Save provisional value to be used in updating at next step time
            saved_v(2) = ts
            saved_v(4) = aa
            saved_v(6) = bb
        else
! **** Steady state model - output is value of surface temperature
            call licoilss(tai,gi,wa,twi,ww,&
                          psychro,icountfl,lcoil,wcoil,thitube,&
                          thifin,ditube,afree,aext,&
                          aratio,bratio,rifin,rofin,heifin,aflow,&
                          hydiam,confin,contube,&
                          ts,tao,go,two,qa,shr,effect,bf,ho,rho,twbo)
        endif
! **** Assign output values
        yout(1)  = ts
        yout(2)  = tao
        yout(3)  = go
        yout(4)  = pai
        yout(5)  = two
        yout(6)  = pwi
        yout(7)  = qa
        yout(8)  = shr
        if (dynamic) then
! **** Effectiveness not defined for dynamic case
            yout(9)  = -99.9
        else
            yout(9)  = effect
        endif
        yout(10) = bf
        if (psychro) then
! **** Auxilliary outputs
            yout(11) = ho
            yout(12) = rho/100.
            yout(13) = twbo
        else
! **** Auxilliary outputs not calculated - set values to flag this
            yout(11) = -99.9
            yout(12) = -99.9
            yout(13) = -99.9
        endif
! **** Allow freezing of algebraic variables
        if (.not.dynamic) then
            iostat(1)=1
        else
            iostat(1)=0
        endif
        do i=2,no
            iostat(i)=1
        enddo
        return
        end subroutine type522
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:  Dynamic or steady state coil and three port valve
! *
! * PURPOSE:     Calculate the outlet air and water conditions for a
! *              heating or cooling coil with a three port valve.
! *              Calculates inlet pressure and coil water flow from
! *              primary circuit flow rate and outlet pressure
! **********************************************************************
! * INPUTS
! * ======
! *  1. tai     : inlet air dry bulb temperature                     (C)
! *  2. gi      : inlet air humidity ratio                       (kg/kg)
! *  3. pao     : outlet air pressure                              (kPa)
! *  4. wa      : dry air mass flow rate                          (kg/s)
! *  5. twi     : inlet water temperature                            (C)
! *  6. wwprim  : primary circuit water mass flow rate            (kg/s)
! *  7. pwo     : outlet water pressure                            (kPa)
! *  8. vstem   : valve stem position                                (-)
! *  9. tsdyn   : effective coil surface temperature                 (C)
! *
! * OUTPUTS
! * =======
! *  1. ts      : effective coil surface temperature                 (C)
! *  2. tao     : outlet air dry bulb temperature                    (C)
! *  3. go      : outlet air humidity ratio                      (kg/kg)
! *  4. pai     : inlet air pressure                               (kPa)
! *  5. two     : outlet water temperature                           (C)
! *  6. pwi     : inlet water pressure                             (kPa)
! *  7. ww      : water mass flow rate though coil                (kg/s)
! *  8. tret    : mixed return water temperature                     (C)
! *  9. qa      : total heat transfer to the air                    (kW)
! * 10. shr     : sensible heat ratio                                (-)
! * 11. effect  : coil effectiveness                                 (-)
! * 12. bf      : coil by-pass factor                                (-)
! * 13. ho      : outlet air specific enthalpy                   (kJ/kg)
! * 14. rho     : outlet air relative humidity                       (-)
! * 15. twbo    : outlet air wet-bulb temperature                    (C)
! *
! * PARAMETERS
! * ==========
! *  1. dynamic : 0 = steady state, 1 = dynamic                      (-)
! *  2. ifault  : 0 = no faults, 1 = parallel flow (cooling coils)   (-)
! *  3. psychro : FALSE = no psychrometric output calcs, TRUE = calcs(-)
! *  4. nrow    : number of rows of tubes                            (-)
! *  5. ntpr    : number of tubes per row                            (-)
! *  6. ncir    : number of parallel water circuits                  (-)
! *  7. lcoil   : length of finned section in direction of flow      (m)
! *  8. hcoil   : height of finned section                           (m)
! *  9. wcoil   : width of finned section                            (m)
! * 10. dotube  : tube outside diameter                              (m)
! * 11. thitube : tube wall thickness                                (m)
! * 12. watube  : tube material (Al=1,Cu=2,Fe=3,CaCO3=4)             (-)
! * 13. spafin  : fin spacing (pitch)                                (m)
! * 14. thifin  : fin thickness                                      (m)
! * 15. wafin   : fin material (Al=1,Cu=2,Fe=3)                      (-)
! * 16. fra     : flow resistance on air side               (0.001 kg.m)
! * 17. frwcoil : coil water flow resistance                (0.001 kg.m)
! * 18. frwbypas: by-pass water flow resistance             (0.001 kg.m)
! * 19. ivaltype: valve type: 0=lin/lin, 1=eq%(flow)/lin(byp), 2=lin/eq%
! * 20. kv      : valve capacity (Kv)                    (m3/hr @ 1 bar)
! * 21. eqpchar : valve curvature parameter (equal percentage port)  (-)
! * 22. sv      : valve rangability                                  (-)
! * 23. cl      : valve leakage (fractional flow)                    (-)
! *
! * SAVED (dynamic mode only)
! * =====
! *  1. time    : time of previous call
! *  2. dt      : solution of differential equation from previous call
! *  3. dtp     : solution of DE from previous step time
! *  4. aa      : A coefficent in dT/dt=A*T+B from previous call
! *  5. aap     : A coefficent in dT/dt=A*T+B from previous step time
! *  6. bb      : B coefficent in dT/dt=A*T+B from previous call
! *  7. bbp     : B coefficent in dT/dt=A*T+B from previous step time
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  Assumes coil is all dry or all wet
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 14, 1995
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  LICOILDY,LICOILSS
!   FUNCTIONS  CALLED:   RLINPORT,REQPPORT
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           Cooling Coil Models to be used in Transient
!                        and/or Wet Regimes -- Theorical Analysis and
!                        Experimental Validation.
!                        X.DING, J-P EPPE,J.LEBRUN,M.WASACZ.
!                        I.E.A. ANNEX 17 document AN17-901019-01.
!                        University of Liege, Belgium.
!
!                        (adapted from the cooling coil model written
!                        by V. Corrado - Polytechnic of Turin, Italy)
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! *
! * Material Properties
! * -------------------
! * volcpm  : volumic specific heat of Al, Cu and mild steel     (kJ/m3)
! * conma   : thermal conductivity of Al, Cu and mild steel   (kW/(m.K))
! *
! * Geometrical characteristics
! * ---------------------------
! * fracfin : fraction of face area blocked by fins                  (-)
! * acgross : gross face area                                       (m2)
! * ntube   : total number of tubes                                  (-)
! * ttubelen: total tube length                                      (m)
! * axo     : gross outside tube surface area                       (m2)
! * ap      : net outside tube surface area, allowing for fins      (m2)
! * as      : net air side fin area                                 (m2)
! * ditube  : tube inside diameter                                   (m)
! * ai      : inside surface area of tubes                          (m2)
! * afree   : exchanger minimum free-flow area on air side          (m2)
! * aext    : exchanger total heat transfer area on air side        (m2)
! * aratio  : ratio of total fin area to total heat transfer area    (-)
! * bratio  : ratio of air side total area : water side internal area(-)
! * rifin   : fin inside radius                                      (m)
! * rofin   : effective fin radius                                   (m)
! * heifin  : effective fin height                                   (m)
! * aflow   : flow area on water side                               (m2)
! * hydiam  : coil hydraulic diameter                                (m)
! * confin  : thermal conductivity of fin material            (kW/(m.K))
! * contube : thermal conductivity of tube material           (kW/(m.K))
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type523(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndtdt=0,ndiffeq=1,&
                                             ni=9,no=15,np=23,&
                                             ns=1+ndiffeq*6
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: lcoil,kv
        integer      :: matube,mafin
        logical      :: dynamic,psychro
! **** Thermophysical constants
        real,dimension(4)   :: volcpm=(/2455.,3475.,3915.,2250./),&
                               conma=(/0.221,0.393,0.0453,0.0013/)
! **** Miscellaneous constants (pi, atmospheric pressure, minimum
! **** water-side resistance
        real         :: pi=3.14159,patmpa=101325.,frsmall=1.0e-8
! **** Break-points in valve characteristics (see documentation for
! **** type329 two port valve)
        real         :: xlin=0.01,xeqp=0.33
        real         :: eqpchar,sv,cl,fracfin,acgross,ttubelen,axo,asx,&
                        ap,ditube,ai,afree,aext,aratio,bratio,rifin,&
                        rofin,heifin,aflow,hydiam,confin,contube,vcpfin,&
                        vcptube,frwvalve,frwtot,pai,pwi,tsp,aap,bbp,&
                        twbo,rho,ho,bf,shr,qa,two,go,tao,bb,aa,aabar,&
                        bbbar,tsbar,ts,effect,tai,gi,gsat,fwphi,pao,&
                        wa,twi,ww,pwo,vstem,tsdyn,hcoil,wcoil,&
                        dotube,thitube,spafin,thifin,fra,frwcoil,&
                        rlinport,reqpport,dpqudlin,wwprim,frwbypas,&
                        vstembyp,frvalflo,frvalbyp,frtotflo,frtotbyp,&
                        frcombin,dp,dpqud,wqud,tret
        integer      :: i,icountfl,ntube,is,ifault,nrow,ntpr,ncir,ivaltype

! **** Read in inputs
        tai     = xin(1)
        gi      = xin(2)
! **** Limit humidity ratio to saturated value
        gsat    = fwphi(tai,100.0,patmpa)
        gi      = min(gi,gsat)
        pao     = xin(3)
        wa      = xin(4)
        twi     = xin(5)
        wwprim  = xin(6)
        pwo     = xin(7)
        vstem   = xin(8)
        tsdyn   = xin(9)
! **** Read in parameters
! **** Dynamic or steady state
        if (nint(par_v(1))==1) then
            dynamic = .true.
        else
            dynamic = .false.
        endif
! **** Fault type
        ifault = nint(par_v(2))
! **** Peform additional pyschrometric calculations?
        if (nint(par_v(3))/=0) then
            psychro = .true.
        else
            psychro = .false.
        endif
! **** Coil geometry and materials
        nrow    = nint(par_v(4))
        ntpr    = nint(par_v(5))
        ncir    = nint(par_v(6))
        lcoil   = par_v(7)
        hcoil   = par_v(8)
        wcoil   = par_v(9)
        dotube  = par_v(10)
        thitube = par_v(11)
        matube  = nint(par_v(12))
        spafin  = par_v(13)
        thifin  = par_v(14)
        mafin   = nint(par_v(15))
! **** Air and water flow resistances
        fra     = par_v(16)
        frwcoil = par_v(17)
        frwbypas= par_v(18)
! **** Valve characteristics
        ivaltype= par_v(19)
        kv      = par_v(20)
        eqpchar = par_v(21)
        sv      = par_v(22)
        cl      = par_v(23)
! **** Heat exchanger geometry
        if (nrow<=2) then
! **** One or two rows -> cross flow
            icountfl=0
        elseif (nrow>2) then
            if (ifault/=1) then
! **** Three or more rows and no connection fault -> counterflow
                icountfl=1
            else
! **** Three or more rows and connection fault -> parallel flow
                icountfl=-1
            endif
        endif
! **** Heat transfer areas
        fracfin = thifin/spafin
        acgross = hcoil*wcoil
        ntube   = nrow*ntpr
        ttubelen= ntube*wcoil
        axo     = pi*dotube*dotube/4.0
        ap      = pi*dotube*ttubelen*(1.0-fracfin)
        asx     = 2.0*(hcoil*lcoil-ntube*axo)*wcoil/spafin
        ditube  = dotube-2.0*thitube
        ai      = pi*ditube*ttubelen
        afree   = acgross*(1.0-ntpr*dotube/hcoil)*(1.0-fracfin)
        aext    = asx+ap
        aratio  = asx/aext
        bratio  = aext/ai
        rifin   = dotube/2.0
        rofin   = sqrt(hcoil*lcoil/(pi*ntube))
        heifin  = rofin-rifin
        aflow   = ncir*pi*ditube*ditube/4.0
        hydiam  = 4.0*afree*lcoil/aext
        confin  = conma(mafin)
        contube = conma(matube)
        vcpfin  = volcpm(mafin)
        vcptube = volcpm(matube)
! **** Valve types and resistances
! **** Limit stem position to range 0-1
        vstem   = max( 0.0, min(1.0,vstem) )
        vstembyp= 1.-vstem
! **** Flow port - select valve type and calculate resistance
        if (ivaltype==0 .or. ivaltype==2) then
! **** Linear characteristic - two segments (linear and close-off)
            frvalflo=rlinport(vstem,kv,sv,xlin,cl)
        elseif (ivaltype==1) then
! **** Equal percentage characteristic - three segments (exponential,
! **** linear and close-off)
            frvalflo=reqpport(vstem,kv,xeqp,eqpchar,xlin,sv,cl)
        else
            stop 'type523: illegal valve type'
        endif
! **** By-pass port - select valve type and calculate resistance
        if (ivaltype==0 .or. ivaltype==1) then
! **** Linear characteristic - two segments (linear and close-off)
            frvalbyp=rlinport(vstembyp,kv,sv,xlin,cl)
        else
! **** Equal percentage characteristic - three segments (exponential,
! **** linear and close-off)
            frvalbyp=reqpport(vstembyp,kv,xeqp,eqpchar,xlin,sv,cl)
        endif
! **** Total water-side resistances - coil + valve, by-pass + valve
        frtotflo  = frwcoil+frvalflo
        frtotbyp  = frwbypas+frvalbyp
        if (frtotflo<frsmall .or. frtotbyp<frsmall) then
! **** Total resistance almost zero
            write(*,*)&
            'type523: water-side flow resistance must not be < ',frsmall
            stop
        endif
! **** Parallel combination of total flow and by-pass resistances
        frcombin  = frtotflo*frtotbyp/(frtotflo+frtotbyp)
! **** Inlet pressure from combined resistance and primary circuit flow
! **** treat as turbulent at low flows for consistency with the way valve
! **** leakage is modeled
        dp      = dpqud(frcombin,wwprim)
        pwi     = pwo+dp
! **** Flow rate through coil
        ww      = wqud(frtotflo,dp)
! **** Calculate air inlet pressure from flow resistance
        pai     = pao+dpqudlin(fra,wa)
! **** Calculate coil duty and outlet conditions
        if (dynamic) then
! **** Initialize at beginning of simulation
            if (itime<=1) then
                if (init==0 .or. saved_v(1)>time) then
                    saved_v(1) = -999999.
               endif
                if (init==0) then
! **** Capacitance of coil initially at inlet air temperature (coil "off")
                    saved_v(2) = tai
                    saved_v(4) = 0.0
                    saved_v(6) = 0.0
                endif
            endif
! **** Every time-step
            if (time>saved_v(1)) then
! **** First call of time-step - update value of temperature rise, aa and bb
! **** from previous time-step
                do is = 2,ns-1,2
                    saved_v(is+1) = saved_v(is)
                enddo
            endif
! **** Update previous values
            tsp = saved_v(3)
            aap = saved_v(5)
            bbp = saved_v(7)
            call licoilab(tai,gi,wa,twi,ww,tsdyn,&
                          psychro,icountfl,lcoil,wcoil,dotube,thitube,&
                          thifin,ttubelen,asx,ditube,afree,aext,&
                          aratio,bratio,rifin,rofin,heifin,aflow,&
                          hydiam,confin,contube,vcpfin,vcptube,&
                          aa,bb,tao,go,two,qa,shr,bf,ho,rho,twbo)
! **** Integrate differential equation for surface temperature
            aabar = (aa + aap)/2.0
            bbbar = (bb + bbp)/2.0
            call diffeq(time,aabar,bbbar,tsp,ts,tsbar)
! **** Save time of current call
            saved_v(1) = time
! **** Save provisional value to be used in updating at next step time
            saved_v(2) = ts
            saved_v(4) = aa
            saved_v(6) = bb
        else
! **** Steady state model - output is value of surface temperature
            call licoilss(tai,gi,wa,twi,ww,&
                          psychro,icountfl,lcoil,wcoil,thitube,&
                          thifin,ditube,afree,aext,&
                          aratio,bratio,rifin,rofin,heifin,aflow,&
                          hydiam,confin,contube,&
                          ts,tao,go,two,qa,shr,effect,bf,ho,rho,twbo)
        endif
! **** Return water temperature
        if (wwprim>1.0e-10) then
            tret    = (ww*two + (wwprim-ww)*twi)/wwprim
        else
            tret    = two
        endif
! **** Assign output values
        yout(1)  = ts
        yout(2)  = tao
        yout(3)  = go
        yout(4)  = pai
        yout(5)  = two
        yout(6)  = pwi
        yout(7)  = ww
        yout(8)  = tret
        yout(9)  = qa
        yout(10) = shr
        if (dynamic) then
! **** Effectiveness not defined for dynamic case
            yout(11)  = -99.9
        else
            yout(11)  = effect
        endif
        yout(12) = bf
        if (psychro) then
! **** Auxilliary outputs
            yout(13) = ho
            yout(14) = rho/100.
            yout(15) = twbo
        else
! **** Auxilliary outputs not calculated - set values to flag this
            yout(13) = -99.9
            yout(14) = -99.9
            yout(15) = -99.9
        endif
! **** Allow freezing of algebraic variables
        if (.not.dynamic) then
            iostat(1)=1
        else
            iostat(1)=0
        endif
        do i=2,no
            iostat(i)=1
        enddo
        return
        end subroutine type523

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:  Dynamic or steady state coil and three port valve
! *
! * PURPOSE:     Calculate the outlet air and water conditions for a
! *              heating or cooling coil with a three port valve.
! *              Calculates the flow through the coil and in the
! *              primary circuit from the inlet and outlet pressures.
! **********************************************************************
! * INPUTS
! * ======
! *  1. tai     : inlet air dry bulb temperature                     (C)
! *  2. gi      : inlet air humidity ratio                       (kg/kg)
! *  3. pao     : outlet air pressure                              (kPa)
! *  4. wa      : dry air mass flow rate                          (kg/s)
! *  5. twi     : inlet water temperature                            (C)
! *  6. pwi     : inlet water pressure                             (kPa)
! *  7. pwo     : outlet water pressure                            (kPa)
! *  8. vstem   : valve stem position                                (-)
! *  9. tsdyn   : effective coil surface temperature                 (C)
! *
! * OUTPUTS
! * =======
! *  1. ts      : effective coil surface temperature                 (C)
! *  2. tao     : outlet air dry bulb temperature                    (C)
! *  3. go      : outlet air humidity ratio                      (kg/kg)
! *  4. pai     : inlet air pressure                               (kPa)
! *  5. two     : outlet water temperature                           (C)
! *  6. wwprim  : primary circuit water mass flow rate            (kg/s)
! *  7. ww      : water mass flow rate though coil                (kg/s)
! *  8. tret    : mixed return water temperature                     (C)
! *  9. qa      : total heat transfer to the air                    (kW)
! * 10. shr     : sensible heat ratio                                (-)
! * 11. effect  : coil effectiveness                                 (-)
! * 12. bf      : coil by-pass factor                                (-)
! * 13. ho      : outlet air specific enthalpy                   (kJ/kg)
! * 14. rho     : outlet air relative humidity                       (-)
! * 15. twbo    : outlet air wet-bulb temperature                    (C)
! *
! * PARAMETERS
! * ==========
! *  1. dynamic : 0 = steady state, 1 = dynamic                      (-)
! *  2. ifault  : 0 = no faults, 1 = parallel flow (cooling coils)   (-)
! *  3. psychro : FALSE = no psychrometric output calcs, TRUE = calcs(-)
! *  4. nrow    : number of rows of tubes                            (-)
! *  5. ntpr    : number of tubes per row                            (-)
! *  6. ncir    : number of parallel water circuits                  (-)
! *  7. lcoil   : length of finned section in direction of flow      (m)
! *  8. hcoil   : height of finned section                           (m)
! *  9. wcoil   : width of finned section                            (m)
! * 10. dotube  : tube outside diameter                              (m)
! * 11. thitube : tube wall thickness                                (m)
! * 12. watube  : tube material (Al=1,Cu=2,Fe=3,CaCO3=4)             (-)
! * 13. spafin  : fin spacing (pitch)                                (m)
! * 14. thifin  : fin thickness                                      (m)
! * 15. wafin   : fin material (Al=1,Cu=2,Fe=3)                      (-)
! * 16. fra     : flow resistance on air side               (0.001 kg.m)
! * 17. frwcoil : coil water flow resistance                (0.001 kg.m)
! * 18. frwbypas: by-pass water flow resistance             (0.001 kg.m)
! * 19. ivaltype: valve type: 0=lin/lin, 1=eq%(flow)/lin(byp), 2=lin/eq%
! * 20. kv      : valve capacity (Kv)                    (m3/hr @ 1 bar)
! * 21. eqpchar : valve curvature parameter (equal percentage port)  (-)
! * 22. sv      : valve rangability                                  (-)
! * 23. cl      : valve leakage (fractional flow)                    (-)
! *
! * SAVED (dynamic mode only)
! * =====
! *  1. time    : time of previous call
! *  2. dt      : solution of differential equation from previous call
! *  3. dtp     : solution of DE from previous step time
! *  4. aa      : A coefficent in dT/dt=A*T+B from previous call
! *  5. aap     : A coefficent in dT/dt=A*T+B from previous step time
! *  6. bb      : B coefficent in dT/dt=A*T+B from previous call
! *  7. bbp     : B coefficent in dT/dt=A*T+B from previous step time
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  Assumes coil is all dry or all wet
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 14, 1995
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  LICOILDY,LICOILSS
!   FUNCTIONS  CALLED:   RLINPORT,REQPPORT
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           Cooling Coil Models to be used in Transient
!                        and/or Wet Regimes -- Theorical Analysis and
!                        Experimental Validation.
!                        X.DING, J-P EPPE,J.LEBRUN,M.WASACZ.
!                        I.E.A. ANNEX 17 document AN17-901019-01.
!                        University of Liege, Belgium.
!
!                        (adapted from the cooling coil model written
!                        by V. Corrado - Polytechnic of Turin, Italy)
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! *
! * Material Properties
! * -------------------
! * volcpm  : volumic specific heat of Al, Cu and mild steel     (kJ/m3)
! * conma   : thermal conductivity of Al, Cu and mild steel   (kW/(m.K))
! *
! * Geometrical characteristics
! * ---------------------------
! * fracfin : fraction of face area blocked by fins                  (-)
! * acgross : gross face area                                       (m2)
! * ntube   : total number of tubes                                  (-)
! * ttubelen: total tube length                                      (m)
! * axo     : gross outside tube surface area                       (m2)
! * ap      : net outside tube surface area, allowing for fins      (m2)
! * as      : net air side fin area                                 (m2)
! * ditube  : tube inside diameter                                   (m)
! * ai      : inside surface area of tubes                          (m2)
! * afree   : exchanger minimum free-flow area on air side          (m2)
! * aext    : exchanger total heat transfer area on air side        (m2)
! * aratio  : ratio of total fin area to total heat transfer area    (-)
! * bratio  : ratio of air side total area : water side internal area(-)
! * rifin   : fin inside radius                                      (m)
! * rofin   : effective fin radius                                   (m)
! * heifin  : effective fin height                                   (m)
! * aflow   : flow area on water side                               (m2)
! * hydiam  : coil hydraulic diameter                                (m)
! * confin  : thermal conductivity of fin material            (kW/(m.K))
! * contube : thermal conductivity of tube material           (kW/(m.K))
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type524(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndtdt=0,ndiffeq=1,&
                                             ni=9,no=15,np=23,&
                                             ns=1+ndiffeq*6
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: lcoil,kv
        integer      :: matube,mafin
        logical      :: dynamic,psychro
! **** Thermophysical constants
        real,dimension(4)   :: volcpm=(/2455.,3475.,3915.,2250./),&
                               conma=(/0.221,0.393,0.0453,0.0013/)
! **** Miscellaneous constants (pi, atmospheric pressure, minimum
! **** water-side resistance
        real         :: pi=3.14159,patmpa=101325.,frsmall=1.0e-8
! **** Break-points in valve characteristics (see documentation for
! **** type329 two port valve)
        real         :: xlin=0.01,xeqp=0.33
        real         :: eqpchar,sv,cl,fracfin,acgross,ttubelen,axo,asx,&
                        ap,ditube,ai,afree,aext,aratio,bratio,rifin,&
                        rofin,heifin,aflow,hydiam,confin,contube,vcpfin,&
                        vcptube,frwvalve,frwtot,pai,pwi,tsp,aap,bbp,&
                        twbo,rho,ho,bf,shr,qa,two,go,tao,bb,aa,aabar,&
                        bbbar,tsbar,ts,effect,tai,gi,gsat,fwphi,pao,&
                        wa,twi,ww,pwo,vstem,tsdyn,hcoil,wcoil,&
                        dotube,thitube,spafin,thifin,fra,frwcoil,&
                        rlinport,reqpport,dpqudlin,frwbypas,&
                        vstembyp,frvalflo,frvalbyp,frtotflo,frtotbyp ,&
                        dp,wqud,wwbypas,wwprim,frprim,frtotcoi,frcoibyp ,&
                        frtot,dpprim,dpqud,dpcoibyp,tret
        integer      :: i,icountfl,ntube,is,ifault,nrow,ntpr,ncir,ivaltype

! **** Read in inputs
        tai     = xin(1)
        gi      = xin(2)
! **** Limit humidity ratio to saturated value
        gsat    = fwphi(tai,100.0,patmpa)
        gi      = min(gi,gsat)
        pao     = xin(3)
        wa      = xin(4)
        twi     = xin(5)
        pwi     = xin(6)
        pwo     = xin(7)
        vstem   = xin(8)
        tsdyn   = xin(9)
! **** Read in parameters
! **** Dynamic or steady state
        if (nint(par_v(1))==1) then
            dynamic = .true.
        else
            dynamic = .false.
        endif
! **** Fault type
        ifault = nint(par_v(2))
! **** Peform additional pyschrometric calculations?
        if (nint(par_v(3))/=0) then
            psychro = .true.
        else
            psychro = .false.
        endif
! **** Coil geometry and materials
        nrow    = nint(par_v(4))
        ntpr    = nint(par_v(5))
        ncir    = nint(par_v(6))
        lcoil   = par_v(7)
        hcoil   = par_v(8)
        wcoil   = par_v(9)
        dotube  = par_v(10)
        thitube = par_v(11)
        matube  = nint(par_v(12))
        spafin  = par_v(13)
        thifin  = par_v(14)
        mafin   = nint(par_v(15))
! **** Air and water flow resistances
        fra     = par_v(16)
        frwcoil = par_v(17)
        frwbypas= par_v(18)
! **** Valve characteristics
        ivaltype= par_v(19)
        kv      = par_v(20)
        eqpchar = par_v(21)
        sv      = par_v(22)
        cl      = par_v(23)
! **** Heat exchanger geometry
        if (nrow<=2) then
! **** One or two rows -> cross flow
            icountfl=0
        elseif (nrow>2) then
            if (ifault/=1) then
! **** Three or more rows and no connection fault -> counterflow
                icountfl=1
            else
! **** Three or more rows and connection fault -> parallel flow
                icountfl=-1
            endif
        endif
! **** Heat transfer areas
        fracfin = thifin/spafin
        acgross = hcoil*wcoil
        ntube   = nrow*ntpr
        ttubelen= ntube*wcoil
        axo     = pi*dotube*dotube/4.0
        ap      = pi*dotube*ttubelen*(1.0-fracfin)
        asx     = 2.0*(hcoil*lcoil-ntube*axo)*wcoil/spafin
        ditube  = dotube-2.0*thitube
        ai      = pi*ditube*ttubelen
        afree   = acgross*(1.0-ntpr*dotube/hcoil)*(1.0-fracfin)
        aext    = asx+ap
        aratio  = asx/aext
        bratio  = aext/ai
        rifin   = dotube/2.0
        rofin   = sqrt(hcoil*lcoil/(pi*ntube))
        heifin  = rofin-rifin
        aflow   = ncir*pi*ditube*ditube/4.0
        hydiam  = 4.0*afree*lcoil/aext
        confin  = conma(mafin)
        contube = conma(matube)
        vcpfin  = volcpm(mafin)
        vcptube = volcpm(matube)
! **** Valve types and resistances
! **** Limit stem position to range 0-1
        vstem   = max( 0.0, min(1.0,vstem) )
        vstembyp= 1.-vstem
! **** Water flow rates through coil and bypass
        if (ivaltype>=0 .and. ivaltype<=2) then
! **** Common port carries common flow
! **** flow port - select valve type and calculate resistance
            if (ivaltype==0 .or. ivaltype==2) then
! **** Linear characteristic - two segments (linear and close-off)
                frvalflo=rlinport(vstem,kv,sv,xlin,cl)
            elseif (ivaltype==1) then
! **** Equal percentage characteristic - three segments (exponential,
! **** linear and close-off)
                frvalflo=reqpport(vstem,kv,xeqp,eqpchar,xlin,sv,cl)
            endif
! **** By-pass port - select valve type and calculate resistance
            if (ivaltype==0 .or. ivaltype==1) then
! **** Linear characteristic - two segments (linear and close-off)
                frvalbyp=rlinport(vstembyp,kv,sv,xlin,cl)
            else
! **** Equal percentage characteristic - three segments (exponential,
! **** linear and close-off)
                frvalbyp=reqpport(vstembyp,kv,xeqp,eqpchar,xlin,sv,cl)
            endif
! **** Total water-side resistances - coil + valve, by-pass + valve
            frtotflo  = frwcoil+frvalflo
            frtotbyp  = frwbypas+frvalbyp
            if (frtotflo<frsmall .or. frtotbyp<frsmall) then
! **** Total resistance almost zero
                write(*,*)&
                'type524: water-side flow resistance must not be < ',&
                frsmall
                stop
            endif
! **** Pressure drop
            dp      = pwi-pwo
! **** Flow rate through coil
            ww      = wqud(frtotflo,dp)
! **** Flow rate through by-pass
            wwbypas = wqud(frtotbyp,dp)
! **** Flow rate in primary circuit
            wwprim  = ww+wwbypas
        elseif (ivaltype==3 .or. ivaltype==4) then
! **** Faulty installation
            if (ivaltype==3) then
! **** Flow and common ports swopped
! **** Resistance of flow port - equal percentage characteristic
                frprim   = reqpport(vstem,kv,xeqp,eqpchar,xlin,sv,cl)
! **** Resistance of bypass port - linear characteristic
                frvalbyp = rlinport(vstembyp,kv,sv,xlin,cl)
! **** Resistance of coil branch (common port has zero resistance)
                frtotcoi = frwcoil
! **** Resistance of bypass branch (balancing valve + bypass port)
                frtotbyp = frwbypas + frvalbyp
! **** Parallel resistance of coil and bypass branches
                frcoibyp = (frtotcoi**(-0.5)+frtotbyp**(-0.5))**2
! **** Total resistance of subsystem
                frtot    = frcoibyp + frprim
! **** Pressure drop across subsystem
                dp       = pwi-pwo
! **** Flow rate in primary circuit
                wwprim   = wqud(frtot,dp)
                dpprim   = dpqud(frprim,wwprim)
                dpcoibyp = dp - dpprim
! **** Flow rate through coil
                ww       = wqud(frtotcoi,dpcoibyp)
! **** Flow rate through by-pass
                wwbypas  = wqud(frtotbyp,dpcoibyp)
            elseif (ivaltype==4) then
! **** Bypass and common ports swopped
! **** Resistance of flow port - equal percentage characteristic
                frvalflo = reqpport(vstem,kv,xeqp,eqpchar,xlin,sv,cl)
! **** Resistance of bypass port - linear characteristic
                frprim   = rlinport(vstembyp,kv,sv,xlin,cl)
! **** Resistance of coil branch (flow port + coil)
                frtotcoi = frvalflo + frwcoil
! **** Resistance of bypass branch (common port has zero resistance)
                frtotbyp = frwbypas
! **** Parallel resistance of coil and bypass branches
                frcoibyp = (frtotcoi**(-0.5)+frtotbyp**(-0.5))**2
! **** Total resistance of subsystem
                frtot    = frcoibyp + frprim
! **** Pressure drop across subsystem
                dp       = pwi-pwo
! **** Flow rate in primary circuit
                wwprim   = wqud(frtot,dp)
                dpprim   = dpqud(frprim,wwprim)
                dpcoibyp = dp - dpprim
! **** Flow rate through coil
                ww       = wqud(frtotcoi,dpcoibyp)
! **** Flow rate through by-pass
                wwbypas  = wqud(frtotbyp,dpcoibyp)
            endif
        else
            stop 'valve type out of range'
        endif
! **** Calculate air inlet pressure from flow resistance
        pai     = pao+dpqudlin(fra,wa)
! **** Calculate coil duty and outlet conditions
        if (dynamic) then
! **** Initialize at beginning of simulation
            if (itime<=1) then
                if (init==0 .or. saved_v(1)>time) then
                    saved_v(1) = -999999.
               endif
                if (init==0) then
! **** Capacitance of coil initially at inlet air temperature (coil "off")
                    saved_v(2) = tai
                    saved_v(4) = 0.0
                    saved_v(6) = 0.0
                endif
            endif
! **** Every time-step
            if (time>saved_v(1)) then
! **** First call of time-step - update value of temperature rise, aa and bb
! **** from previous time-step
                do is = 2,ns-1,2
                    saved_v(is+1) = saved_v(is)
                enddo
            endif
! **** Update previous values
            tsp = saved_v(3)
            aap = saved_v(5)
            bbp = saved_v(7)
            call licoilab(tai,gi,wa,twi,ww,tsdyn,&
                          psychro,icountfl,lcoil,wcoil,dotube,thitube,&
                          thifin,ttubelen,asx,ditube,afree,aext,&
                          aratio,bratio,rifin,rofin,heifin,aflow,&
                          hydiam,confin,contube,vcpfin,vcptube,&
                          aa,bb,tao,go,two,qa,shr,bf,ho,rho,twbo)
! **** Integrate differential equation for surface temperature
            aabar = (aa + aap)/2.0
            bbbar = (bb + bbp)/2.0
            call diffeq(time,aabar,bbbar,tsp,ts,tsbar)
! **** Save time of current call
            saved_v(1) = time
! **** Save provisional value to be used in updating at next step time
            saved_v(2) = ts
            saved_v(4) = aa
            saved_v(6) = bb
        else
! **** Steady state model - output is value of surface temperature
            call licoilss(tai,gi,wa,twi,ww,&
                          psychro,icountfl,lcoil,wcoil,thitube,&
                          thifin,ditube,afree,aext,&
                          aratio,bratio,rifin,rofin,heifin,aflow,&
                          hydiam,confin,contube,&
                          ts,tao,go,two,qa,shr,effect,bf,ho,rho,twbo)
        endif
! **** Return water temperature
        if (wwprim>1.0e-10) then
            tret    = (ww*two + (wwprim-ww)*twi)/wwprim
        else
            tret    = two
        endif
! **** Assign output values
        yout(1)  = ts
        yout(2)  = tao
        yout(3)  = go
        yout(4)  = pai
        yout(5)  = two
        yout(6)  = wwprim
        yout(7)  = ww
        yout(8)  = tret
        yout(9)  = qa
        yout(10) = shr
        if (dynamic) then
! **** Effectiveness not defined for dynamic case
            yout(11)  = -99.9
        else
            yout(11)  = effect
        endif
        yout(12) = bf
        if (psychro) then
! **** Auxilliary outputs
            yout(13) = ho
            yout(14) = rho/100.
            yout(15) = twbo
        else
! **** Auxilliary outputs not calculated - set values to flag this
            yout(13) = -99.9
            yout(14) = -99.9
            yout(15) = -99.9
        endif
! **** Allow freezing of algebraic variables
        if (.not.dynamic) then
            iostat(1)=1
        else
            iostat(1)=0
        endif
        do i=2,no
            iostat(i)=1
        enddo
        return
        end subroutine type524

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     Motorized pressure-independent VAV box
! *
! * PURPOSE:        Calculate the inlet pressure of a motor-driven,
! *                 pressure-independent VAV box from the outlet
! *                 pressure, the mass flow rate and the demanded
! *                 normalized flow rate.
! *
! *********************************************************************
! * INPUTS
! * ======
! *  1. wair    : air mass flow rate                             (kg/s)
! *  2. pout    : outlet pressure                                 (kPa)
! *  3. fflowset: demanded flowrate (normalized to nominal)       (0-1)
! *
! * OUTPUTS
! * =======
! *  1. pin     : inlet pressure                                  (kPa)
! *  2. damppos : damper position (0=closed, 1=open)                (-)
! *  3. fsped   : fractional motor velocity (0-1)                   (-)
! *  4. tssrev  : number of stops/starts/reversals                  (-)
! *  5. ttrav   : total distance travelled by valve/damper          (-)
! *
! * PARAMETERS
! * ==========
! *  1. volumnom: nominal volumetric flow rate                   (m3/s)
! *  2. farea   : face area of damper(s)                           (m2)
! *  3. pdopen  : pressure drop at nominal flow with damper open  (kPa)
! *  4. ttran   : travel time (0-90 deg)                            (s)
! *  5. fspedmin: minimum fractional motor speed                    (-)
! *  6. hys     : hysteresis                                        (-)
! *  7. kp      : controller gain (frac speed/frac error)           (-)
! *
! *  SAVED
! *  =====
! *  1. time    : time of previous call for rate limit
! *  2. actpos  : actuator position at previous call
! *  3. actposp : actuator position at previous time-step
! *  4. fsped   : actuator speed at previous call
! *  5. fspedp  : actuator speed at prev time-step
! *  6. isped   : +1=forward, -1=backward, 0=stopped at previous call
! *  7. ispedp  : +1=forward, -1=backward, 0=stopped at prev time-step
! *  8. tssrev  : total no. of stop/starts/reversals at previous call
! *  9. tssrevp : total no. of stop/starts/reversals at prev time-step
! * 10. ttrav   : total distance travelled at previous call
! * 11. ttravp  : total distance travelled at previous timestep
! * 12. damppos : position of final control element at previous call
! * 13. dampposp: position of final control element at previous step time
! *
! *********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                June 16, 1994
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   HYSTER2
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           Haves, P., Component-Based Modelling of VAV
!                        Systems, Proc. System Simulation in Buildings
!                        '94, Liege, Belgium, December 1994
!
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type525(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndtdt=0,ndiffeq=0,&
                                             ni=3,no=5,np=7,&
                                             ns=13
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: kp
        logical      :: con
        real         :: rhoa=1.2
! **** Damper coefficients (ln(k)=a+b*theta) from legg
        real,dimension(0:1)   :: alegg=(/-1.51,-1.51/),&
                                 blegg=(/0.105,0.0842/)
! **** Loss coefficient of fully open damper
        real         :: kopen=0.2
! **** Assumed fractional flow leakage of closed damper
        real         :: fleak=1.0e-2
! **** Assume opposed or single blade damper
        integer      :: ipx=0
! **** Assume damper is fully open at start
        real         :: startpos=1.0
        real         :: volumnom,farea,pdopen,ttran,fspedmin,hys,actposp,&
                        fspedp,tssrevp,ttravp,dampposp,fflow,fspeddem,&
                        fsped,fspedbar,actpos,tssrev,ttrav,damppos,ropen,&
                        alegg1,blegg1,r,radd,rtot,pin,wair,pout,fflowset,&
                        spedlim,clip,rdamper,dpqudlin
        integer      :: i,is,ispedp,isped

! **** Read in inputs
        wair     = xin(1)
        pout     = xin(2)
        fflowset = xin(3)
! **** Read in par_vameters
        volumnom = par_v(1)
        if (volumnom<=0.0) then
            stop 'type 525: nominal flow rate must be greater than zero'
        endif
        farea    = par_v(2)
        pdopen   = par_v(3)
        ttran    = par_v(4)
        fspedmin = par_v(5)
        hys      = par_v(6)
        kp       = par_v(7)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
                saved_v(2)  = startpos
                do is = 4,(ns-3),2
                    saved_v(is) = 0.0
                enddo
                saved_v(12) = startpos
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update previous sample instant values
            do is=2,(ns-1),2
                saved_v(is+1) = saved_v(is)
            enddo
        endif
! **** Determine curent position based on situation at previous step time
! **** Update previous values
        actposp  = saved_v(3)
        fspedp   = saved_v(5)
        ispedp   = saved_v(7)
        tssrevp  = saved_v(9)
        ttravp   = saved_v(11)
        dampposp = saved_v(13)
! **** Limit control signal
        if (fflowset<0.0) fflowset = 0.0
        if (fflowset>1.0) fflowset = 1.0
! **** Measured fractional flow at end of time-step
        fflow    = (wair/rhoa)/volumnom
! **** Fractional speed of motor at end of time-step - proportional to
! **** control error and limited to range fspedmin - 1.0 or zero
        fspeddem = kp*(fflowset-fflow)
        fsped    = spedlim(fspeddem,fspedmin)
! **** Average speed over time-step
        fspedbar = (fsped+fspedp)/2.0
! **** New actuator position
        actpos   = actposp+fspedbar*tstep/ttran
        actpos   = clip(actpos,0.0,1.0)
! **** Operational statistics
! **** Stops/starts/reversals
! **** Current motor direction
        if (fsped>0.0) then
            isped = +1
        elseif (fsped<0.0) then
            isped = -1
        else
            isped = 0
        endif
! **** Previous stops/starts/reversals
        if (isped/=ispedp) then
            tssrev = tssrevp+1.0
        else
            tssrev = tssrevp
        endif
! **** Distance travelled
        ttrav = ttravp+abs(actpos-actposp)
! ****  hysteresis (no correction for reduction in range)
        call hystrsis(actpos,dampposp,hys,damppos)
! **** Damper resistance
! **** Resistance when fully open
        ropen = 1.0e-3*kopen/(2.0*rhoa*farea*farea)
! **** Resistance at position damppos
        alegg1 = alegg(0)                 ! added 12/23/99
        blegg1 = blegg(0)                 ! added 12/23/99
        r = rdamper(damppos,ropen,fleak,farea,alegg1,blegg1,ipx)
! **** Resistance of box, diffuser etc, excluding damper
        radd = pdopen/(rhoa*volumnom*rhoa*volumnom) - ropen
        rtot = r+radd
! **** Pressure drop and inlet pressure
        pin = pout+dpqudlin(rtot,wair)
! **** Save time of current call
        saved_v(1)  = time
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = actpos
        saved_v(4)  = fsped
        saved_v(6)  = float(isped)
        saved_v(8)  = tssrev
        saved_v(10) = ttrav
        saved_v(12) = damppos
! **** Output
        yout(1) = pin
        yout(2) = damppos
        yout(3) = fsped
        yout(4) = tssrev
        yout(5) = ttrav
! **** Determine whether to allow freezing
! **** Freezing allowed if motor stopped and demanded flow no changing
        con = (iostat(3)<-1).or.(iostat(3)==2)
        if (fspedbar==0.0.and.con) then
            do i=1,no
                iostat(i) = 1
            enddo
        else
            do i=1,no
                iostat(i) = 1
            enddo
        endif
! **** Return
        return
        end subroutine type525

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     Pressure-independent VAV box - calculates flow
! *
! * PURPOSE:        Calculate the flow rate through a motor-driven,
! *                 pressure-independent VAV box from the inlet
! *                 pressure, the outlet pressure and the demanded
! *                 normalized flow rate.
! *
! *********************************************************************
! * INPUTS
! * ======
! *  1. pin     : inlet pressure                                  (kPa)
! *  2. pout    : outlet pressure                                 (kPa)
! *  3. fflowset: demanded flowrate (normalized to nominal)       (0-1)
! *  4. wair    : air mass flow rate                             (kg/s)
! *
! * OUTPUTS
! * =======
! *  1. wair    : air mass flow rate                             (kg/s)
! *  2. damppos : damper position (0=closed, 1=open)                (-)
! *  3. fsped   : fractional motor velocity (0-1)                   (-)
! *  4. tssrev  : number of stops/starts/reversals                  (-)
! *  5. ttrav   : total distance travelled by valve/damper          (-)
! *
! * PARAMETERS
! * ==========
! *  1. volumnom: nominal volumetric flow rate                   (m3/s)
! *  2. farea   : face area of damper(s)                           (m2)
! *  3. pdopen  : pressure drop at nominal flow with damper open  (kPa)
! *  4. ttran   : travel time (0-90 deg)                            (s)
! *  5. fspedmin: minimum fractional motor speed                    (-)
! *  6. hys     : hysteresis                                        (-)
! *  7. kp      : controller gain (frac speed/frac error)           (-)
! *
! *  SAVED
! *  =====
! *  1. time    : time of previous call for rate limit
! *  2. actpos  : actuator position at previous call
! *  3. actposp : actuator position at previous time-step
! *  4. fsped   : actuator speed at previous call
! *  5. fspedp  : actuator speed at prev time-step
! *  6. isped   : +1=forward, -1=backward, 0=stopped at previous call
! *  7. ispedp  : +1=forward, -1=backward, 0=stopped at prev time-step
! *  8. tssrev  : total no. of stop/starts/reversals at previous call
! *  9. tssrevp : total no. of stop/starts/reversals at prev time-step
! * 10. ttrav   : total distance travelled at previous call
! * 11. ttravp  : total distance travelled at previous timestep
! * 12. damppos : position of final control element at previous call
! * 13. dampposp: position of final control element at previous step time
! *
! *********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                June 16, 1994
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   HYSTER2
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           Haves, P., Component-Based Modelling of VAV
!                        Systems, Proc. System Simulation in Buildings
!                        '94, Liege, Belgium, December 1994
!
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type526(xin,yout,par_v,saved_v,iostat)

         use modsim_head
        implicit none
        integer,parameter                 :: ndtdt=0,ndiffeq=0,&
                                             ni=4,no=5,np=7,&
                                             ns=13
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: kp
        logical      :: con
        real         :: rhoa=1.2
! **** Damper coefficients (ln(k)=a+b*theta) from legg
        real,dimension(0:1)   :: alegg=(/-1.51,-1.51/),&
                                 blegg=(/0.105,0.0842/)
! **** Loss coefficient of fully open damper
        real         :: kopen=0.2
! **** Assumed fractional flow leakage of closed damper
        real         :: fleak=1.0e-2
! **** Assume opposed or single blade damper
        integer      :: ipx=0
! **** Assume damper is fully open at start
        real         :: startpos=1.0
        real         :: pin,pout,fflowset,wair,volumnom,farea,pdopen,&
                        ttran,fspedmin,hys,actposp,fspedp,&
                        tssrevp,ttravp,dampposp,fflow,fspeddem,fsped,&
                        spedlim,fspedbar,actpos,clip,tssrev,ttrav,&
                        damppos,ropen,alegg1,blegg1,r,rdamper,radd,&
                        rtot,dp,wqudlin
        integer      :: i,is,ispedp,isped

! **** Read in inputs
        pin      = xin(1)
        pout     = xin(2)
        fflowset = xin(3)
        wair     = xin(4)
! **** Read in parameters
        volumnom = par_v(1)
        if (volumnom<=0.0) then
            stop 'type 526: nominal flow rate must be greater than zero'
        endif
        farea    = par_v(2)
        pdopen   = par_v(3)
        ttran    = par_v(4)
        fspedmin = par_v(5)
        hys      = par_v(6)
        kp       = par_v(7)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
                saved_v(2)  = startpos
                do is = 4,(ns-3),2
                    saved_v(is) = 0.0
                enddo
                saved_v(12) = startpos
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update previous sample instant values
            do is=2,(ns-1),2
                saved_v(is+1) = saved_v(is)
            enddo
        endif
! **** Determine curent position based on situation at previous step time
! **** Update previous values
        actposp  = saved_v(3)
        fspedp   = saved_v(5)
        ispedp   = saved_v(7)
        tssrevp  = saved_v(9)
        ttravp   = saved_v(11)
        dampposp = saved_v(13)
! **** Limit control signal
        if (fflowset<0.0) fflowset = 0.0
        if (fflowset>1.0) fflowset = 1.0
! **** Measured fractional flow at end of time-step
        fflow    = (wair/rhoa)/volumnom
! **** Fractional speed of motor at end of time-step - proportional to
! **** control error and limited to range fspedmin - 1.0 or zero
        fspeddem = kp*(fflowset-fflow)
        fsped    = spedlim(fspeddem,fspedmin)
! **** Average speed over time-step
        fspedbar = (fsped+fspedp)/2.0
!       fspedbar = fspedp
! **** New actuator position
        actpos   = actposp+fspedbar*tstep/ttran
        actpos   = clip(actpos,0.0,1.0)
! **** Operational statistics
! **** Stops/starts/reversals
! **** Current motor direction
        if (fsped>0.0) then
            isped = +1
        elseif (fsped<0.0) then
            isped = -1
        else
            isped = 0
        endif
! **** Previous stops/starts/reversals
        if (isped/=ispedp) then
            tssrev = tssrevp+1.0
        else
            tssrev = tssrevp
        endif
! **** Distance travelled
        ttrav = ttravp+abs(actpos-actposp)
! ****  hysteresis (no correction for reduction in range)
        call hystrsis(actpos,dampposp,hys,damppos)
! **** Damper resistance
! **** Resistance when fully open
        ropen = 1.0e-3*kopen/(2.0*rhoa*farea*farea)
! **** Resistance at position damppos
        alegg1 = alegg(0)                 ! added 12/23/99
        blegg1 = blegg(0)                 ! added 12/23/99
        r = rdamper(damppos,ropen,fleak,farea,alegg1,blegg1,ipx)
! **** Resistance of box, diffuser etc, excluding damper
        radd = pdopen/(rhoa*volumnom*rhoa*volumnom) - ropen
        rtot = r+radd
! **** Pressure drop and mass flow rate
        dp   = pin-pout
        wair = wqudlin(rtot,dp)
! **** Save time of current call
        saved_v(1)  = time
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = actpos
        saved_v(4)  = fsped
        saved_v(6)  = float(isped)
        saved_v(8)  = tssrev
        saved_v(10) = ttrav
        saved_v(12) = damppos
! **** Output
        yout(1) = wair
        yout(2) = damppos
        yout(3) = fsped
        yout(4) = tssrev
        yout(5) = ttrav
! **** Determine whether to allow freezing
! **** Freezing allowed if motor stopped and demanded flow no changing
        con = (iostat(3)<-1).or.(iostat(3)==2)
        if (fspedbar==0.0.and.con) then
            do i=1,no
                iostat(i) = 1
            enddo
        else
            do i=1,no
                iostat(i) = 1
            enddo
        endif
! **** Return
        return
        end subroutine type526
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     Pressure-independent VAV box (implicit flow)
! *
! * PURPOSE:        Calculate the flow rate through a motor-driven,
! *                 pressure-independent VAV box from the design
! *                 pressure drop and the demanded normalized flow rate.
! *
! *********************************************************************
! * INPUTS
! * ======
! *  1. fflowset: demanded flowrate (normalized to nominal)       (0-1)
! *  2. wair    : air mass flow rate                             (kg/s)
! *  3. sfstatus: supply fan status (1 = on, 0 = off)               (-)
! *
! * OUTPUTS
! * =======
! *  1. wair    : mass flow rate                                 (kg/s)
! *  2. damppos : damper position (0=closed, 1=open)                (-)
! *  3. fsped   : fractional motor velocity (0-1)                   (-)
! *  4. tssrev  : number of stops/starts/reversals                  (-)
! *  5. ttrav   : total distance travelled by valve/damper          (-)
! *
! * PARAMETERS
! * ==========
! *  1. volumnom: nominal volumetric flow rate                   (m3/s)
! *  2. farea   : face area of damper(s)                           (m2)
! *  3. pdopen  : pressure drop at nominal flow with damper open  (kPa)
! *  4. ttran   : travel time (0-90 deg)                            (s)
! *  5. fspedmin: minimum fractional motor speed                    (-)
! *  6. hys     : hysteresis                                        (-)
! *  7. kp      : controller gain (frac speed/frac error)           (-)
! *
! *  SAVED
! *  =====
! *  1. time    : time of previous call for rate limit
! *  2. actpos  : actuator position at previous call
! *  3. actposp : actuator position at previous time-step
! *  4. fsped   : actuator speed at previous call
! *  5. fspedp  : actuator speed at prev time-step
! *  6. isped   : +1=forward, -1=backward, 0=stopped at previous call
! *  7. ispedp  : +1=forward, -1=backward, 0=stopped at prev time-step
! *  8. tssrev  : total no. of stop/starts/reversals at previous call
! *  9. tssrevp : total no. of stop/starts/reversals at prev time-step
! * 10. ttrav   : total distance travelled at previous call
! * 11. ttravp  : total distance travelled at previous timestep
! * 12. damppos : position of final control element at previous call
! * 13. dampposp: position of final control element at previous step time
! *
!**********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                June 16, 1994
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   HYSTER2
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           Haves, P., Component-Based Modelling of VAV
!                        Systems, Proc. System Simulation in Buildings
!                        '94, Liege, Belgium, December 1994
!
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type527(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndtdt=0,ndiffeq=0,&
                                             ni=3,no=5,np=7,&
                                             ns=13
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: kp
        logical      :: con,sfstatus,retolog
        real         :: rhoa=1.2
! **** Damper coefficients (ln(k)=a+b*theta) from legg
        real,dimension(0:1)   :: alegg=(/-1.51,-1.51/),&
                                 blegg=(/0.105,0.0842/)
! **** Loss coefficient of fully open damper
        real         :: kopen=0.2
! **** Assumed fractional flow leakage of closed damper
        real         :: fleak=1.0e-2
! **** Assume opposed or single blade damper
        integer      :: ipx=0
! **** Assume damper is fully open at start
        real         :: startpos=1.0
        real         :: fflowset,wair,volumnom,farea,pdopen,ttran,&
                        fspedmin,hys,actposp,fspedp,tssrevp,ttravp,&
                        dampposp,fflow,fspeddem,fsped,spedlim,fspedbar,&
                        actpos,clip,tssrev,ttrav,damppos,ropen,alegg1,&
                        blegg1,r,rdamper,radd,rtot,wqudlin
        integer      :: i,is,ispedp,isped

! **** Read in inputs
        fflowset = xin(1)
        wair     = xin(2)
        sfstatus = retolog(xin(3))
! **** Read in parameters
        volumnom = par_v(1)
        if (volumnom<=0.0) then
            stop 'type 527: nominal flow rate must be greater than zero'
        endif
        farea    = par_v(2)
        pdopen   = par_v(3)
        ttran    = par_v(4)
        fspedmin = par_v(5)
        hys      = par_v(6)
        kp       = par_v(7)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1)  = 0.
            endif
            if (init==0) then
                saved_v(2)  = startpos
                do is = 4,(ns-3),2
                    saved_v(is) = 0.0
                enddo
                saved_v(12) = startpos
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update previous sample instant values
            do is=2,(ns-1),2
                saved_v(is+1) = saved_v(is)
            enddo
        endif
! **** Determine curent position based on situation at previous step time
! **** Update previous values
        actposp  = saved_v(3)
        fspedp   = saved_v(5)
        ispedp   = saved_v(7)
        tssrevp  = saved_v(9)
        ttravp   = saved_v(11)
        dampposp = saved_v(13)
! **** Limit control signal
        if (fflowset<0.0) fflowset = 0.0
        if (fflowset>1.0) fflowset = 1.0
! **** Measured fractional flow at end of time-step
        fflow    = (wair/rhoa)/volumnom
! **** Fractional speed of motor at end of time-step - proportional to
! **** control error and limited to range fspedmin - 1.0 or zero
        fspeddem = kp*(fflowset-fflow)
        fsped    = spedlim(fspeddem,fspedmin)
! **** Average speed over time-step
        fspedbar = (fsped+fspedp)/2.0
! **** New actuator position
        actpos   = actposp+fspedbar*tstep/ttran
        actpos   = clip(actpos,0.0,1.0)
! **** Operational statistics
! **** Stops/starts/reversals
! **** Current motor direction
        if (fsped>0.0) then
            isped = +1
        elseif (fsped<0.0) then
            isped = -1
        else
            isped = 0
        endif
! **** Previous stops/starts/reversals
        if (isped/=ispedp) then
            tssrev = tssrevp+1.0
        else
            tssrev = tssrevp
        endif
! **** Distance travelled
        ttrav = ttravp+abs(actpos-actposp)
! ****  Hysteresis (no correction for reduction in range)
        call hystrsis(actpos,dampposp,hys,damppos)
! **** Calculate damper resistance and flow rate if supply fan on
        if (sfstatus) then
! **** Resistance when fully open
            ropen = 1.0e-3*kopen/(2.0*rhoa*farea*farea)
! **** Resistance at position damppos
            alegg1 = alegg(0)            ! added 12/23/99
            blegg1 = blegg(0)            ! added 12/23/99
            r = rdamper(damppos,ropen,fleak,farea,alegg1,blegg1,ipx)
! **** Resistance of box, diffuser etc, excluding damper
            radd = pdopen/(rhoa*volumnom*rhoa*volumnom) - ropen
            rtot = r+radd
! **** Calculate flow rate from total resistance and external pressure drop
           wair = wqudlin(rtot,pdopen)
        else
! **** Supply fan off - set flow rate to zero
            wair = 0.0
        endif
! **** Save time of current call
        saved_v(1)  = time
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = actpos
        saved_v(4)  = fsped
        saved_v(6)  = float(isped)
        saved_v(8)  = tssrev
        saved_v(10) = ttrav
        saved_v(12) = damppos
! **** Output
        yout(1) = wair
        yout(2) = damppos
        yout(3) = fsped
        yout(4) = tssrev
        yout(5) = ttrav
! **** Determine whether to allow freezing
! **** Freezing allowed if motor stopped and demanded flow not changing
        con = (iostat(1)<-1).or.(iostat(1)==2)
        if (fspedbar==0.0.and.con) then
            do i=1,no
                iostat(i) = 1
            enddo
        else
            do i=1,no
                iostat(i) = 1
            enddo
        endif
! **** Return
        return
        end subroutine type527

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     Flow split and pressure-independent VAV box
! *
! * PURPOSE:        Calculate the outlet flow rates and the inlet pressure
! *                 for a flow split with a motorized pressure-independent
! *                 VAV box in one branch, given the outlet pressures, the
! *                 inlet flow rate and the demanded normalized flow rate
! *                 through the VAV box.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. wi      : inlet mass flow rate                            (kg/s)
! *  2. pomain  : pressure at main outlet                          (kPa)
! *  3. povav   : pressure at VAV box outlet                       (kPa)
! *  4. fflowset: demanded VAV flowrate (normalized to nominal)    (0-1)
! *  5. wovav   : mass flow rate through VAV box                  (kg/s)
! *
! * OUTPUTS
! * =======
! *  1. pi      : inlet pressure                                   (kPa)
! *  2. womain  : mass flow rate at main outlet                   (kg/s)
! *  3. wovav   : mass flow rate through VAV box                  (kg/s)
! *  4. damppos : damper position (0=closed, 1=open)                 (-)
! *  5. fsped   : fractional motor velocity (0-1)                    (-)
! *  6. tssrev  : number of stops/starts/reversals                   (-)
! *  7. ttrav   : total distance travelled by valve/damper           (-)
! *
! * PARAMETERS
! * ==========
! *  1. rin     : inlet resistance                           (0.001/k.m)
! *  2. rmain   : resistance of main outlet                  (0.001/k.m)
! *  3. rbranch : resistance of branch to VAV box            (0.001/k.m)
! *  4. volumnom: nominal volumetric flow rate for VAV box        (m3/s)
! *  5. farea   : face area of damper(s)                            (m2)
! *  6. pdopen  : box pressure drop at nom. flow with damper open  (kPa)
! *  7. ttran   : actuator travel time (0-90 deg)                    (s)
! *  8. fspedmin: minimum fractional motor speed                     (-)
! *  9. hys     : hysteresis                                         (-)
! * 10. kp      : controller gain (frac speed/frac error)            (-)
! *
! *  SAVED
! *  =====
! *  1. time    : time of previous call for rate limit
! *  2. actpos  : actuator position at previous call
! *  3. actposp : actuator position at previous time-step
! *  4. fsped   : actuator speed at previous call
! *  5. fspedp  : actuator speed at prev time-step
! *  6. isped   : +1=forward, -1=backward, 0=stopped at previous call
! *  7. ispedp  : +1=forward, -1=backward, 0=stopped at prev time-step
! *  8. tssrev  : total no. of stop/starts/reversals at previous call
! *  9. tssrevp : total no. of stop/starts/reversals at prev time-step
! * 10. ttrav   : total distance travelled at previous call
! * 11. ttravp  : total distance travelled at previous timestep
! * 12. damppos : position of final control element at previous call
! * 13. dampposp: position of final control element at previous step time
! *
! *********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                June 16, 1994
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  FLOWSPLT
!   FUNCTIONS  CALLED:   HYSTRSIS
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           Haves, P., Component-Based Modelling of VAV
!                        Systems, Proc. System Simulation in Buildings
!                        '94, Liege, Belgium, December 1994
!
! **********************************************************************
! *
! * INTERNAL VARIABLES
! * ==================
! * wcritu  : flow threshold, lower limit of pure quadratic solution
! * wcritl  : flow threshold, upper limit of pure linear solution
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type528(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndtdt=0,ndiffeq=0,&
                                             ni=5,no=7,np=10,&
                                             ns=13
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: kp
        logical      :: con
        real         :: rhoa=1.2
! **** Damper coefficients (ln(k)=a+b*theta) from legg
        real,dimension(0:1)   :: alegg=(/-1.51,-1.51/),&
                                 blegg=(/0.105,0.0842/)
! **** Loss coefficient of fully open damper
        real         :: kopen=0.2
! **** Assumed fractional flow leakage of closed damper
        real         :: fleak=1.0e-2
! **** Assume opposed or single blade damper
        integer      :: ipx=0
! **** Assume damper is fully open at start
        real         :: startpos=1.0
        real         :: wcritl=3.e-2, wcritu=6.e-2
        real         :: wi,pomain,povav,fflowset,wovav,rin,rmain,&
                        rbranch,volumnom,farea,pdopen,ttran,fspedmin,&
                        hys,actposp,fspedp,tssrevp,ttravp,dampposp,&
                        fflow,fspeddem,fsped,spedlim,fspedbar,actpos,&
                        clip,tssrev,ttrav,damppos,ropen,alegg1,blegg1,&
                        r,rdamper,radd,rvav,womain,pi
        integer      :: i,is,ispedp,isped,ifail

! **** Read in inputs
        wi       = xin(1)
        pomain   = xin(2)
        povav    = xin(3)
        fflowset = xin(4)
        wovav    = xin(5)
! **** Check for reverse inlet flow at last iteration
        if (wi<-0.01 .and. lastit==1)&
            write(*,*) 'warning:', ' type 528: ',&
                       'reverse inlet flow'
! **** Read in parameters
        rin      = par_v(1)
        rmain    = par_v(2)
        rbranch  = par_v(3)
        volumnom = par_v(4)
        if (volumnom<=0.0) then
            stop 'type 528: nominal flow rate must be greater than zero'
        endif
        farea    = par_v(5)
        pdopen   = par_v(6)
        ttran    = par_v(7)
        fspedmin = par_v(8)
        hys      = par_v(9)
        kp       = par_v(10)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
                saved_v(2)  = startpos
                do is = 4,(ns-3),2
                    saved_v(is) = 0.0
                enddo
                saved_v(12) = startpos
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update previous sample instant values
            do is=2,(ns-1),2
                saved_v(is+1) = saved_v(is)
            enddo
        endif
! **** Determine curent position based on situation at previous step time
! **** Update previous values
        actposp  = saved_v(3)
        fspedp   = saved_v(5)
        ispedp   = saved_v(7)
        tssrevp  = saved_v(9)
        ttravp   = saved_v(11)
        dampposp = saved_v(13)
! **** Limit control signal
        if (fflowset<0.0) fflowset = 0.0
        if (fflowset>1.0) fflowset = 1.0
! **** Measured fractional flow at end of time-step
        fflow    = (wovav/rhoa)/volumnom
! **** Fractional speed of motor at end of time-step - proportional to
! **** Control error and limited to range fspedmin - 1.0 or zero
        fspeddem = kp*(fflowset-fflow)
        fsped    = spedlim(fspeddem,fspedmin)
! **** Average speed over time-step
        fspedbar = (fsped+fspedp)/2.0
! **** New actuator position
        actpos   = actposp+fspedbar*tstep/ttran
        actpos   = clip(actpos,0.0,1.0)
! **** Operational statistics
! **** Stops/starts/reversals
! **** Current motor direction
        if (fsped>0.0) then
            isped = +1
        elseif (fsped<0.0) then
            isped = -1
        else
            isped = 0
        endif
! **** Increment stops/starts/reversals
        if (isped/=ispedp) then
            tssrev = tssrevp+1.0
        else
            tssrev = tssrevp
        endif
! **** Increment distance travelled
        ttrav = ttravp+abs(actpos-actposp)
! ****  Hysteresis (no correction for reduction in range)
        call hystrsis(actpos,dampposp,hys,damppos)
! **** Damper resistance
! **** Resistance when fully open
        ropen = 1.0e-3*kopen/(2.0*rhoa*farea*farea)
! **** Resistance at position damppos
        alegg1 = alegg(0)               ! added 12/23/1999
        blegg1 = blegg(0)               ! added 12/23/99
        r = rdamper(damppos,ropen,fleak,farea,alegg1,blegg1,ipx)
! **** Resistance of box, diffuser etc, excluding damper
        radd = pdopen/(rhoa*volumnom*rhoa*volumnom) - ropen
        rvav = r+radd
! **** Call flow split routine
        call flowsplt(wi,pomain,povav,rin,rmain,rvav,wcritl,wcritu,&
                      rtolx,pi,womain,wovav,ifail)
! **** Check for unsuccesful completion
        if (ifail==1) then
! **** One or more resistances negative
            write(*,*) ' type 528: ',&
                       'resistances must not be negative'
            stop
        elseif (ifail==2) then
! **** Zero resistance for both outlet branches
            write(*,*) ' type 528: ',&
                       'both outlet resistances cannot be zero'
            stop
        endif
! **** Save time of current call
        saved_v(1)  = time
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = actpos
        saved_v(4)  = fsped
        saved_v(6)  = float(isped)
        saved_v(8)  = tssrev
        saved_v(10) = ttrav
        saved_v(12) = damppos
! **** Output
        yout(1) = pi
        yout(2) = womain
        yout(3) = wovav
        yout(4) = damppos
        yout(5) = fsped
        yout(6) = tssrev
        yout(7) = ttrav
! **** Determine whether to allow freezing
! **** Freezing allowed if motor stopped and demanded flow no changing
        con = (iostat(4)<-1).or.(iostat(4)==2)
        if (fspedbar==0.0.and.con) then
            do i=1,no
                iostat(i) = 1
            enddo
        else
            do i=1,no
                iostat(i) = 1
            enddo
        endif
! **** Return
        return
        end subroutine type528


