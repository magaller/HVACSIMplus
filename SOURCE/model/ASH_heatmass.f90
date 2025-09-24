! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:  Dynamic or steady state heating and cooling coil
! *
! * PURPOSE:     Calculate the outlet air and water conditions for a
! *              heating or cooling coil
! **********************************************************************
! * INPUTS
! * ======
! *  1. tai     : inlet air dry bulb temperature                     (C)
! *  2. gi      : inlet air humidity ratio                       (kg/kg)
! *  3. pao     : outlet air gauge pressure                        (kPa)
! *  4. wa      : dry air mass flow rate                          (kg/s)
! *  5. twi     : inlet water temperature                            (C)
! *  6. ww      : water mass flow rate                            (kg/s)
! *  7. pwo     : outlet water gauge pressure                      (kPa)
! *  8. tsdyn   : effective coil surface temperature                 (C)
! *
! * OUTPUTS
! * =======
! *  1. ts      : effective coil surface temperature                 (C)
! *  2. tao     : outlet air dry bulb temperature                    (C)
! *  3. go      : outlet air humidity ratio                      (kg/kg)
! *  4. pai     : inlet air gauge static pressure                  (kPa)
! *  5. two     : outlet water temperature                           (C)
! *  6. pwi     : inlet water gauge pressure                       (kPa)
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
! *  3. psychro : 0 = no psychrometric output calcs, 1 = calcs       (-)
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
! * 17. frw     : flow resistance on water side             (0.001 kg.m)
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
!   SUBROUTINES CALLED:  LICOILAB,LICOILSS
!   FUNCTIONS  CALLED:
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
! * cpw     : water specific heat                            (kJ/(kg.K))
! * cpa     : dry air specific heat                          (kJ/(kg.K))
! * cpg     : water vapor specific heat                      (kJ/(kg/K))
! * volcpm  : volumic specific heat of Al, Cu and mild steel     (kJ/m3)
! * conma   : thermal conductivity of Al, Cu and mild steel   (kW/(m.K))
! * rhoa    : density of air                                     (kg/m3)
! * rhow    : water density                                      (kg/m3)
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

        subroutine type362(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=8,no=13,np=17,&
                                             ndiffeq=1,ns=1+ndiffeq*6
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: lcoil
        integer      :: matube,mafin
        logical      :: dynamic,psychro
! **** Thermophysical constants
        real,dimension(4)     :: volcpm=(/2455.,3475.,3915.,2250./),&
                                 conma =(/0.221,0.393,0.0453,0.0013/)
! **** Miscellaneous constants
        real         :: pi=3.14159, patmpa=101325.
        real         :: spafin,thifin,fracfin,acgross,ttubelen,axo,&
                        ap,ditube,ai,afree,aext,aratio,ratio,rifin,&
                        rofin,heifin,aflow,hydiam,confin,contube,&
                        vcpfin,vcptube,tsp,aap,bbp,twbo,rho,ho,bf,&
                        shr,qa,two,go,tao,bb,aa,aabar,bbbar,tsbar,ts,&
                        effect,tai,gi,gsat,fwphi,pao,wa,twi,ww,pwo,&
                        tsdyn,fra,pai,dpqudlin,frw,pwi,hcoil,wcoil,&
                        dotube,thitube,bratio,asx
        integer      :: i,is,icountfl,ntube,ifault,nrow,ntpr,ncir

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
        tsdyn   = xin(8)
! **** Calculate air and water inlet pressures from flow resistances
        pai     = pao+dpqudlin(fra,wa)
        pwi     = pwo+dpqudlin(frw,ww)
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
        frw     = par_v(17)
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
! **** Return
        return
        end subroutine type362


! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Ideal heating of fluid stream
! *
! * PURPOSE:    Treats either instantaneous sensible heaters or heat
! *             addition by other components, e.g. fans
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. w       : mass flow rate                                  (kg/s)
! *  2. ti      : inlet temperature                                  (C)
! *  3. qa      : rate of heat addition                             (kW)
! *
! * OUTPUT
! * ======
! *  1. to      : outlet temperature                                 (C)
! *
! * PARAMETERS
! * ==========
! *  1. mode    : fluid: 1 = air, any other value + water            (-)
! *  2. tcon    : time constant                                      (s)
! *
! *  SAVED
! *  =====
! *  1. time    : time of previous call
! *  2. dt      : solution of differential equation from previous call
! *  3. dtp     : solution of differential equation from previous step
! *               time
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  None
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                October 20, 1993
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   WQUDLIN
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type366(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=2,&
                                             ndiffeq=1,ns=1+ndiffeq*2
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: cpa=1.0, cpw=4.18
        real         :: w,ti,qa,tcon,cp,dtp,dt0,aa,bb,dtba,dt,to,dtbar
        integer      :: mode

! **** Read in inputs
        w    = xin(1)
        ti   = xin(2)
        qa   = xin(3)
! **** Read in parameters
        mode = nint(par_v(1))
        tcon = par_v(2)
! **** Set specific heat for appropriate fluid
        if (mode==1) then
            cp=cpa
        else
            cp=cpw
        endif
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Source of heat (e.g.\ fan motor) off at start of simulation
                saved_v(2) = 0.0
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update value of temperature rise from
! **** previous time-step
                saved_v(3) = saved_v(2)
        endif
! **** Update previous values
        dtp = saved_v(3)
! **** Determine temperature rise
! **** Equilibrium temperature rise and coefficients of differential
! **** equation
        if (w>0.0) then
            dt0 = qa/(w*cp)
            aa  = -1./tcon
            bb  = dt0/tcon
        else
! **** Zero flow - set derivative of temperature rise to zero until flow
! **** restarts
            dt0 = 0.0              ! added  12/22/99
            aa  = 0.0
            bb  = 0.0
        endif
        call diffeq(time,aa,bb,dtp,dt,dtbar)
! **** Calculate outlet temperature
        to = ti+dt
! **** Save time of current call
        saved_v(1) = time
! **** Save provisional value to be used in updating at next step time
        saved_v(2) = dt
! **** Assign output
        yout(1) = to
! **** Determine whether to allow freezing
! **** Freezing of the output is allowed only if the input is a constant
! **** or a boundary variable and the change in the output is small
        if (((iostat(1)<-1).or.(iostat(1)==2)) .and.&
                ((abs(dt0 - dt))<=(rtolx * abs(dt)+atolx))) then
            iostat(1) = 1
        else
            iostat(1) = 0
        endif
! **** Return
        return
        end subroutine type366

! *********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Mixing of two moist air streams
! *
! * PURPOSE:    Calculates the mixed temperature, humidity and mass
! *             flow rate resulting from mixing two air streams
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. t(1)    : Stream 1 dry bulb temperature                      (C)
! *  2. g(1)    : Stream 1 humidity ratio                        (kg/kg)
! *  3. w(1)    : Stream 1 mass flow rate                         (kg/s)
! *  4. t(2)    : Stream 2 dry bulb temperature                      (C)
! *  5. g(2)    : Stream 2 humidity ratio                        (kg/kg)
! *  6. w(2)    : Stream 2 mass flow rate                         (kg/s)
! *
! * OUTPUTS
! * =======
! *  1. tmix    : Mixed air dry bulb temperature                     (C)
! *  2. gmix    : Mixed air humidity ratio                       (kg/kg)
! *  3. wmix    : Mixed air mass flow rate                        (kg/s)
! *
! * PARAMETERS
! * ==========
! * (none)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 8, 1995
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  MOISTMIX
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type367(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=6,no=3,np=1,ns=0
        integer,parameter                 :: nstream=2
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nstream)           :: t,gx,w
        real         :: wmix,gmix,tmix
        integer      :: i


! **** read in inputs
        t(1) = xin(1)
        gx(1) = xin(2)
        w(1) = xin(3)
        t(2) = xin(4)
        gx(2) = xin(5)
        w(2) = xin(6)
! **** calculate mixed air conditions
! **** moistmix also tests for reverse flows and sets the mixed condition
! **** to the forward flow condition, or to an average of the flow conditions
! **** if both are reverse
        call moistmix(t,gx,w,nstream,tmix,gmix,wmix)
! **** assign outputs
        yout(1) = tmix
        yout(2) = gmix
        yout(3) = wmix
! **** allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo

        return
        end subroutine type367
! *********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Mixing of six moist air streams
! *
! * PURPOSE:    Calculates the mixed temperature, humidity and mass
! *             flow rate resulting from mixing six air streams
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. t(1)    : Stream 1 dry bulb temperature                      (C)
! *  2. g(1)    : Stream 1 humidity ratio                        (kg/kg)
! *  3. w(1)    : Stream 1 mass flow rate                         (kg/s)
! *  4. t(2)    : Stream 2 dry bulb temperature                      (C)
! *  5. g(2)    : Stream 2 humidity ratio                        (kg/kg)
! *  6. w(2)    : Stream 2 mass flow rate                         (kg/s)
! *  7. t(3)    : Stream 3 dry bulb temperature                      (C)
! *  8. g(3)    : Stream 3 humidity ratio                        (kg/kg)
! *  9. w(3)    : Stream 3 mass flow rate                         (kg/s)
! * 10. t(4)    : Stream 4 dry bulb temperature                      (C)
! * 11. g(4)    : Stream 4 humidity ratio                        (kg/kg)
! * 12. w(4)    : Stream 4 mass flow rate                         (kg/s)
! * 13. t(5)    : Stream 5 dry bulb temperature                      (C)
! * 14. g(5)    : Stream 5 humidity ratio                        (kg/kg)
! * 15. w(5)    : Stream 5 mass flow rate                         (kg/s)
! * 16. t(6)    : Stream 6 dry bulb temperature                      (C)
! * 17. g(6)    : Stream 6 humidity ratio                        (kg/kg)
! * 18. w(6)    : Stream 6 mass flow rate                         (kg/s)
! *
! * OUTPUTS
! * =======
! *  1. tmix    : Mixed air dry bulb temperature                     (C)
! *  2. gmix    : Mixed air humidity ratio                       (kg/kg)
! *  3. wmix    : Mixed air mass flow rate                        (kg/s)
! *
! * PARAMETERS
! * ==========
! * (none)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 8, 1995
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  MOISTMIX
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type368(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=18,no=3,np=1,ns=0
        integer,parameter                 :: nstream=6
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nstream)           :: t,gx,w
        real         :: wmix,gmix,tmix
        integer      :: i

! **** Read in inputs
        t(1) = xin(1)
        gx(1) = xin(2)
        w(1) = xin(3)
        t(2) = xin(4)
        gx(2) = xin(5)
        w(2) = xin(6)
        t(3) = xin(7)
        gx(3) = xin(8)
        w(3) = xin(9)
        t(4) = xin(10)
        gx(4) = xin(11)
        w(4) = xin(12)
        t(5) = xin(13)
        gx(5) = xin(14)
        w(5) = xin(15)
        t(6) = xin(16)
        gx(6) = xin(17)
        w(6) = xin(18)
! **** Calculate mixed air conditions
! **** Moistmix also tests for reverse flows and sets the mixed condition
! **** to the forward flow condition, or to an average of the flow conditions
! **** if both are reverse
        call moistmix(t,gx,w,nstream,tmix,gmix,wmix)
! **** assign outputs
        yout(1) = tmix
        yout(2) = gmix
        yout(3) = wmix
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo

        return
        end subroutine type368

! *********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Supply and return air flow rates
! *
! * PURPOSE:    Determines supply air flow rate by summing up to six
! *             zone flow rates, determines return air flow rate from
! *             supply flow rate by subtracting specified flow rate
! *             difference if OA damper is open. Sets both flows to
! *             zero if supply fan status is 'OFF'. Determines return
! *             flow rate from each zone, assuming equal fractional
! *             leakage from each zone.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. ws(1)   : Zone 1 supply mass flow rate                    (kg/s)
! *  2. ws(2)   : Zone 2 supply mass flow rate                    (kg/s)
! *  3. ws(3)   : Zone 3 supply mass flow rate                    (kg/s)
! *  4. ws(4)   : Zone 4 supply mass flow rate                    (kg/s)
! *  5. ws(5)   : Zone 5 supply mass flow rate                    (kg/s)
! *  6. ws(6)   : Zone 6 supply mass flow rate                    (kg/s)
! *  7. sfstatus: Supply fan status (1 = on, 0 = off)                (-)
! *  8. oastatus: Outside air damper status (0=closed, 1=open)       (-)
! *
! * OUTPUTS
! * =======
! *  1. wsup    : Supply air mass flow rate                       (kg/s)
! *  2. wret    : Return air mass flow rate                       (kg/s)
! *  3. wr(1)   : Zone 1 return mass flow rate                    (kg/s)
! *  4. wr(2)   : Zone 2 return mass flow rate                    (kg/s)
! *  5. wr(3)   : Zone 3 return mass flow rate                    (kg/s)
! *  6. wr(4)   : Zone 4 return mass flow rate                    (kg/s)
! *  7. wr(5)   : Zone 5 return mass flow rate                    (kg/s)
! *  8. wr(6)   : Zone 6 return mass flow rate                    (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. deltw   : Supply-return flow rate diff if OA damper open  (kg/s)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 8, 1995
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  MOISTMIX
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type369(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nstream=6
        integer,parameter                 :: ni=nstream+2,no=nstream+2,&
                                             np=1,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nstream)           :: ws,wr

        logical      :: sfstatus,oastatus,retolog
        real         :: deltw,wsup,wret,rfrac
        integer      :: i

        if (itime<=1) then                          ! 12/21/99
            if (init==0) then
                do i = 1,nstream
                    wr(i) = 0.0
                end do
            endif
        endif
! **** Read in inputs
        do i=1,nstream
            ws(i) = xin(i)
        enddo
        sfstatus  = retolog(xin(7))
        oastatus  = retolog(xin(8))
! **** Read in parameter
        deltw     = par_v(1)
! **** Calculate supply and return flow rates
        if (sfstatus) then
! **** Supply fan on - sum zone flow rates to get supply flow rate
            wsup = 0.0
            do i=1,nstream
                wsup = wsup+ws(i)
            enddo
            if (oastatus) then
! **** (Minimum) outside air damper open - use specified flow rate
! **** difference to calculate return flow rate
                wret = max((wsup-deltw),0.0)
            else
! **** (Minimum) outside air damper closed - return flow rate equals
! **** supply flow rate
                wret = wsup
            endif
            if (wsup>0.0) then
                rfrac = wret/wsup
                do i=1,nstream
                    wr(i) = rfrac*ws(i)
                enddo
            endif
        else
! ****  Supply fan off - set  flow rates to zero
            wsup = 0.0
            wret = 0.0
            do i=1,nstream
                wr(i) = 0.0
            enddo
        endif
! **** Assign outputs
        yout(1) = wsup
        yout(2) = wret
        do i=1,nstream
            yout(2+i) = wr(i)
        enddo
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo

        return
        end subroutine type369
        
