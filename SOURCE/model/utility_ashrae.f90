! ======================================================================

!    Utility_ashrae.f90
!
!    Updated:  January 19, 2007 Cheol Park, NIST
!              Fortran 77 code was converted into Fortran 90

! ======================================================================

! **********************************************************************
!
        subroutine diffeq(dummy,aa,bb,ti,tf,tbar)
!
! **********************************************************************
!   SUBROUTINE:  Analytical integration
!
!   PURPOSE:     Integrate dT/dt = A*T + B over the current time-step
!                Adapted for HVACSIM+ from the TRNSYS routine DIFFEQ.
!                See the TRNSYS manual for further details.
!  
!                The main difference is that there is no need for special
!                code to deal with the first step time, since initial 
!                values in HVACSIM+ apply to the start of the time-step
!                and the routine calculates the value of T at the end of
!                the time-step.
!  
! **********************************************************************
!   INPUTS                                                            
!   ======  
!   dummy   : (not need here - retained for compatibility with TRNSYS)
!   aa      : A coefficent, averaged over the time-step
!   bb      : B coefficent, averaged over the time-step
!   ti      : value of T at the beginning of the time-step  
!  
!   OUTPUTS
!   =======
!   tf      : value of T at the end of the time-step
!   tbar    : value of T averaged over the time-step
!  
! **********************************************************************
!
        use precision1
        use modsim_head
        implicit none
        real(kind=pp)       :: tbar, tf, ti, aa, bb, dummy

! ** Solution depends on whether the A coefficient is zero
        if(abs(aa) == 0.0) then
! ** A coefficient is zero - linear solution

            tf = bb*tstep + ti
            tbar = (tf + ti)/2.0
        else
! ** A coefficient is non-zero - exponential solution
            tf = (ti + bb/aa)*exp(aa*tstep) - bb/aa
            tbar = (ti + bb/aa)/aa/tstep*(exp(aa*tstep)- 1.0) - bb/aa
        endif
        return
        end subroutine diffeq

! **********************************************************************
!  
        subroutine licoilab(tai,gi,wa,twi,ww,tsdyn,                     &
                           psychro,icountfl,lcoil,wcoil,dotube,thitube, &
                           thifin,ttubelen,as,ditube,afree,aext,        &
                           aratio,bratio,rifin,rofin,heifin,aflow,      &
                           hydiam,confin,contube,vcpfin,vcptube,        &
                           aa,bb,tao,go,two,qa,shr,bf,ho,rho,twbo)
!  
! ********************************************************************** 
!
!   SUBROUTINE:  Dynamic heating and cooling coil
!  
!   PURPOSE:     Calculate the dynamic outlet air and water conditions
!                for a heating or cooling coil 
! **********************************************************************
!   INPUTS                                                            
!   ======                                                            
!   tai     : inlet air dry bulb temperature                         (C)
!   gi      : inlet air humidity ratio                           (kg/kg)
!   wa      : dry air mass flow rate                              (kg/s)
!   twi     : inlet water temperature                                (C)
!   ww      : water mass flow rate                                (kg/s)
!   tsdyn   : effective coil surface temperature                     (C)
!                                                                     
!   OUTPUTS                                                           
!   =======                                                           
!   aa      : A coefficent in DTS/DT=A*TS+B where DTS/DT is the 
!             derivative of effective coil surface temperature      (/s)
!   bb      : B coefficent in DTS/DT=A*TS+B where DTS/DT is the 
!             derivative of effective coil surface temperature     (K/s)
!   tao     : outlet air dry bulb temperature                        (C)
!   go      : outlet air humidity ratio                          (kg/kg)
!   two     : outlet water temperature                              (C)
!   qa      : total heat transfer to the air                        (kW)
!   shr     : sensible heat ratio                                    (-)
!   bf      : coil by-pass factor                                    (-)
!   ho      : outlet air specific enthalpy                       (kJ/kg)
!   rho     : outlet air relative humidity                           (-)
!   twbo    : outlet air wet-bulb temperature                        (C)
!                                                                        
!   PARAMETERS                                                          
!   ==========                                                          
!   icountfl: 0 = cross-flow, 1 = counterflow, -1 = parallel flow    (-)
!   psychro : FALSE = no psychrometric output calcs, TRUE = calcs    (-)
!   lcoil   : length of finned section in direction of flow          (m)   
!   wcoil   : width of finned section in direction of flow           (m)   
!   dotube  : tube outside diameter                                  (m)
!   thitube : tube thickness                                         (m)
!   thifin  : fin thickness                                          (m)
!   ttubelen: total tube length                                      (m)
!   as      : net air side tube area                                (m2)
!   ditube  : tube inside diameter                                   (m)
!   afree   : exchanger minimum free-flow area on air side          (m2)
!   aext    : exchanger total heat transfer area on air side        (m2) 
!   aratio  : ratio of total fin area to total heat transfer area    (-)
!   bratio  : ratio of air side total area : water side intrnl area  (-)
!   rifin   : fin inside radius                                      (m)
!   rofin   : effective fin radius                                   (m)
!   heifin  : effective fin height                                   (m)
!   aflow   : flow area on water side                               (m2)
!   hydiam  : coil hydraulic diameter                                (m) 
!   confin  : thermal conductivity of fin material            (kW/(m.K)) 
!   contube : thermal conductivity of tube material           (kW/(m.K)) 
!   vcpfin  : volumetric heat capacity of fin material       (kJ/(m3.K)) 
!   vcptube : volumetric heat capacity of tube material      (kJ/(m3.K))  
!  
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
!   SUBROUTINES CALLED:  DRYOUTCO, TSHBYPAS, WETOUTCO
!   FUNCTIONS  CALLED:   AIRCOEFF, AIRFINRES, BYPASFAC, COILCAP,
!                        QAIRSIDE, QEFFNTU, QWATSIDE, WATERRES
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
!   INTERNAL VARIABLES                                           
!   ==================                                         
!  
!   Material Properties
!   -------------------
!   rhow    : water density                                      (kg/m3)  
!   cpw     : water specific heat                            (kJ/(kg.K))  
!   cpa     : dry air specific heat                          (kJ/(kg.K))  
!   cpg     : water vapor specific heat                      (kJ/(kg/K)) 
!   tref    : reference temperature                                  (K)  
!   kw      : water thermal conductivity                      (kW/(m.K))  
!   viscw   : water dynamic viscosity                         (kg/(m.s))
!   rhoa    : density of air                                     (kg/m3)
!   rhow    : density of water                                   (kg/m3)
!   cpai    : inlet air mass specific heat                   (kJ/(kg.C)) 
!  
!   Heat Transfer Variables
!   -----------------------
!   va      : air velocity (referred to free area)                 (m/s)
!   bf      : coil by-pass factor                                    (-)
!   cair    : air capacity rate                               (kJ/(s.K))  
!   cwat    : water capacity rate                             (kJ/(s.K))
!   effdry  : dry coil effectiveness                                 (-)
!   hadry   : dry heat transfer coefficient on air side      (kW/(K.m2))
!   radry   : dry air side thermal resistance                  (K.m2/kW) 
!   rw      : water side thermal resistance                    (K.m2/kW)
!   rm      : metal thermal resistance                         (K.m2/kW) 
!   rtdry   : global thermal resistance                        (K.m2/kW)
!   rtwet   : equivalent wet global thermal resistance         (K.m2/kW)
!   antuwet : number of equivalent transfer units (wet conditions)   (-)   
!   cs      : average saturation specific heat               (kJ/(kg.K))
!   
!   Miscellaneous
!   -------------
!   vsmall  : threshold velocity for significant heat transfer     (m/s)
!   patmpa  : outlet air absolute pressure in Pascals for 
!             psychrometric routines                                (Pa)
!   tdpai   : inlet air dew point temperature                        (C)
! **********************************************************************        
        use precision1
        implicit none
        real(kind=pp)    :: tai,gi,wa,twi,ww,tsdyn,                  &
                            lcoil,wcoil,dotube,thitube,     &
                            thifin,ttubelen,as,ditube,afree,aext,    &
                            aratio,bratio,rifin,rofin,heifin,aflow,  &
                            hydiam,confin,contube,vcpfin,vcptube,    &
                            aa,bb,tao,go,two,qa,shr,bf,ho,rho,twbo
        real(kind=pp)    :: t,csat,cpai,cair,cwat,va,vw,gair,hadry,  &
                            radry,rw,rm,rtdry,aircoeff,airfinres,    &
                            waterres,ftdew,qadry,qairside,qw,qwatside, &
                            !dryoutco,rst1,rst1wet,bfwet,hi,fhair, &
                            !wetoutco,rtwet,tdpai,w,twbi,ftwb,tsm,cs, &
                            rst1,rst1wet,bfwet,hi,fhair, &
                            rtwet,tdpai,w,twbi,ftwb,tsm,cs, &
                            hawet,rawet,cairsat,qawet,bfdry,hccoil, &
                            dpai,antuw,cst2,rst1dry,rst2,coilcap, &
                            fwphi,fphi,gsat,effect,a=100.0
        integer          :: icountfl

        logical          :: psychro
! ** Thermophysical constants
        real(kind=pp)    :: rhoa = 1.2, rhow = 1000., cpa = 1., cpw = 4.18, &
                            cpg = 1.805, patmpa = 101325.
! ** Miscellaneous constants
        real(kind=pp)    :: vsmall = 1.e-3
! ** Saturation specific heat as a function of temperature
        csat(t)=1.786+2.27e-2*t+2.96565e-3*t*t
! ** Moist air specific heat and air and water capacity rates and velocities
        cpai    = cpa+cpg*gi
        cair    = cpai*wa
        cwat    = cpw*ww
        va      = wa/(rhoa*afree)
        vw      = ww/(rhow*aflow)
! ** Calculation of coil performance depends on whether there is
! ** significant flow on either side. The four possible cases are
! ** considered separately.
        if (va>vsmall.and.vw>vsmall) then
! ** General case: Significant flow on both air side and water side
! ** performance depends on whether the coil is wet or dry. 
! ** In this model, the coil is taken to be completely
! ** dry or completely wet, depending on which has the higher duty.
! ** Calculate thermal resistances referred to unit area on the air side
            gair   = wa/afree
            hadry  = aircoeff(tai,gair,hydiam,lcoil)
            radry  = airfinres(hadry,confin,thifin,rofin,rifin,&
                               heifin,aratio)
            rw     = waterres(ww,aflow,twi,ditube,bratio,wcoil)
            rm     = bratio*thitube/contube
            rtdry  = radry+rm+rw
! ** Inlet dew point temperature
            tdpai   = ftdew(gi,patmpa)
! ** Dynamic model - Calculate separate resistances on air and water sides,
! ** each referred to the inlet temperature (see printed documentation).
! ** Associate a heat capacity with the surface equal to that of the
! ** metal plus the water.
! ** Dry coil air side heat transfer. Resistance includes the fins
            qadry  = qairside(tai,tsdyn,cair,aext,radry,rst1dry,&
                              bfdry)
! ** Compare surface temperature to inlet dew point temperature
            if (tsdyn>=tdpai) then
! ** Either coil is completely dry or dry (cooling) duty must exceed
! ** wet (cooling) duty. Calculate water side heat transfer for dry coil
                qw     = qwatside(aext,cair,cwat,rtdry,&
                                  icountfl,rst1dry,twi,tsdyn,rst2)
! ** Air side duty, resistance and by-pass factor
                qa     = qadry      
                rst1   = rst1dry 
                bf     = bfdry
! ** Outlet conditions
                call dryoutco(tai,gi,qa,cair,wa,psychro,patmpa,tao,&
                              go,shr,ho,rho,twbo)
            else
! ** Coil may be wet or dry. Calculate air side heat transfer for wet coil
! ** use wet bulb as driving temperature and enhance air side heat transfer 
! ** coefficient by the ratio of the saturation specific heat to the 
! ** sensible specific heat. This is a good approximation to using enthalpy
! ** as the driving potential.
                twbi   = ftwb(tai,gi,patmpa)
! ** Saturation specific heat at average of driving temperatures
                tsm    = (twbi+twi)/2.0
                cs     = csat(tsm)
! ** Enhanced air side heat transfer coefficient and reduced resistance
                hawet  = hadry*cs/cpai
                rawet  = airfinres(hawet,confin,thifin,rofin,rifin,&
                            heifin,aratio)
! ** Air side heat transfer
                cairsat=cs*wa
                qawet  = qairside(twbi,tsdyn,cairsat,aext,rawet,&
                                  rst1wet,bfwet)
! ** Compare dry and wet cooling duties
                if (-qadry>=-qawet) then
! ** Dry coil - Calculate outlet air condition and water side heat transfer
! ** Air side heat transfer, resistance and bypass factor
                    qa     = qadry
                    rst1   = rst1dry
                    bf     = bfdry
! ** Air outlet condition
                    call dryoutco(tai,gi,qa,cair,wa,psychro,patmpa,&
                                  tao,go,shr,ho,rho,twbo)
! ** Water side heat transfer
                    qw     = qwatside(aext,cair,cwat,rtdry,&
                                      icountfl,rst1dry,twi,tsdyn,rst2)
                else
! ** Wet coil - calculate outlet air condition and water side heat transfer
! ** Air side duty, resistance and by-pass factor
                    qa     = qawet
                    rst1   = rst1wet
                    bf     = bfwet
! ** Outlet enthalpy
                    hi     = fhair(tai,gi)
                    ho     = hi+qa/wa
! ** Outlet air dry bulb and humidity ratio
                    call wetoutco(tai,qa,ho,cair,bf,tsdyn,patmpa, &
                                  psychro,shr,tao,go,rho,twbo)
! ** Calculate water side resistance from total resistance and
! ** air side resistance under wet conditons
                    rtwet  = rawet+rw+rm
                    qw     = qwatside(aext,cairsat,cwat,rtwet, &
                                      icountfl,rst1wet,twi,tsdyn,rst2)
                endif
! ** End of wet/dry model
            endif
! ** Outlet water temperature
            two    = twi-qw/cwat
! ** Heat capacity of fins, tubes and water
            hccoil = coilcap(vcpfin,as,thifin,vcptube,&
                            ttubelen,dotube,thitube,ditube)
! ** Coefficients of dtsdyn/dt=aa*tsdyn+bb (nb dtsdyn/dt=(qw-qa)/hccoil)
            aa     = -(rst1+rst2)/(rst1*rst2*hccoil)
            bb     = (tai/rst1+twi/rst2)/hccoil
! ** End of significant flows on both sides
        elseif ((va>vsmall.and.vw<=vsmall)) then
! ** Water flow rate is negligible and air flow rate is significant
! ** calculate rate of change of surface temperature from heat  
! ** balance on coil heat capacity.
! ** Air side resistances and capacity rate
            gair   = wa/afree
            hadry  = aircoeff(tai,gair,hydiam,lcoil)
            radry  = airfinres(hadry,confin,thifin,rofin,rifin,&
                            heifin,aratio)
! ** Air side duty assuming dry coil
            qadry  = qairside(tai,tsdyn,cair,aext,radry,rst1dry,&
                              bfdry)
! ** Inlet dew point temperature
            tdpai   = ftdew(gi,patmpa)
            if (tsdyn>tdpai) then
! ** Coil is dry or dry coil heat flow exceeds wet coil heat flow
! ** Air side heat transfer, resistance and bypass factor
                qa     = qadry      
                rst1   = rst1dry
                bf     = bfdry
! ** Outlet air conditions
                call dryoutco(tai,gi,qa,cair,wa,psychro,patmpa,&
                              tao,go,shr,ho,rho,twbo)
! ** End dry model
            else
! ** Coil may be wet or dry
! ** Calculate air side heat transfer for wet coil
! ** Air side resistance
                twbi   = ftwb(tai,gi,patmpa)
                tsm    = (twbi+tsdyn)/2.0
                cs     = csat(tsm)
                hawet  = hadry*cs/cpai
                rawet  = airfinres(hawet,confin,thifin,rofin,rifin,&
                           heifin,aratio)
! ** Air side heat transfer
                cairsat=cs*wa
                qawet  = qairside(twbi,tsdyn,cairsat,aext,rawet,&
                                  rst1wet,bfwet)
                if (-qadry>=-qawet) then
! ** Coil effectively dry
! ** Air side heat transfer, resistance and bypass factor
                    qa     = qadry      
                    rst1   = rst1dry
                    bf     = bfdry
! ** Outlet air conditions
                    call dryoutco(tai,gi,qa,cair,wa,psychro,patmpa,&
                                  tao,go,shr,ho,rho,twbo)
                else
! ** Wet coil - calculate outlet air condition and water side heat transfer
! ** Air side duty, resistance and by-pass factor
                    qa     = qawet
                    rst1   = rst1wet
                    bf     = bfwet
! ** Outlet enthalpy
                    hi     = fhair(tai,gi)
                    ho     = hi+qa/wa
! ** Outlet air dry bulb and humidity ratio
                    call wetoutco(tai,qa,ho,cair,bf,tsdyn,patmpa,&
                                  psychro,shr,tao,go,rho,twbo)
! ** End wet model
                endif
! ** End wet/dry model
            endif
! ** No water side heat flow - outlet water at temperature of surface
            two    = tsdyn
! ** Heat capacity of fins, tubes and water
            hccoil = coilcap(vcpfin,as,thifin,vcptube,&
                            ttubelen,dotube,thitube,ditube)
! ** Coefficients of dtsdyn/dt=aa*tsdyn+bb (nb dtsdyn/dt=-qa/hccoil)
            aa     = -1./(rst1*hccoil)
            bb     = tai/(rst1*hccoil)
! ** End negligible water flow, significant air flow
        elseif (va<=vsmall.and.vw>vsmall) then
! ** air flow rate negligible and water flow rate significant.
! ** Calculate rate of change of surface temperature from heat
! ** balance on coil heat capacity (this yields a particularly poor 
! ** approximation to the dynamics when the air flow is negligible, but
! ** continues to maintain an energy balance until the air flow becomes
! ** significant again).
! ** Water side duty - Calculate from ntu-effectiveness relationship
! ** assuming infinite capacity rate for the coil surface
            rw     = waterres(ww,aflow,twi,ditube,bratio,wcoil)
            antuw  = aext/(rw*ww*cpw)      
            cst2   = (1.-exp(-antuw))*cwat
            qw     = cst2*(twi-tsdyn)
! ** Heat capacity of fins, tubes and water
            hccoil = coilcap(vcpfin,as,thifin,vcptube,&
                            ttubelen,dotube,thitube,ditube)
! ** Coefficients of dtsdyn/dt=aa*tsdyn+bb (nb dtsdyn/dt=qw/hccoil)
            aa     = -cst2/hccoil
            bb     = twi*cst2/hccoil
! ** Outlet air is in equilibrium with surface
            tao    = tsdyn
! ** Heat balance on water side
            two    = twi-qw/cwat
! ** Outlet humidity ratio is equal to inlet value or saturated value,
! ** whichever is lower

            go     = min(gi,fwphi(tao,a,patmpa))
!            go     = min(gi,fwphi(tao,100.,patmpa))
            qa     = 0.
            effect = 1.
            bf     = 0.
! ** Sensible heat ratio is indeterminate
            shr    = 0.
            if (psychro) then
                ho     = fhair(tao,go)
                rho    = fphi(tao,go,patmpa)
                twbo   = ftwb(tao,go,patmpa)
            endif
! ** End of case of negligible air flow
        else
! ** Both flow rates negligible
! ** set outlet conditions equal to effective surface conditions
! ** as reasonable starting point for solution when one or other flow
! ** rate becomes significant
            tao    = tsdyn
            gsat   = fwphi(tao,a,patmpa)
            go     = min(gi,gsat)
            two    = tsdyn
! ** No heat transfer on either side - no change in surface temperature
            aa     = 0.
            bb     = 0.
            qa     = 0.
            effect = 1.
            bf     = 0.
            shr    = 0.
            if (psychro) then
                ho     = fhair(tao,go)
                rho    = fphi(tao,go,patmpa)
                twbo   = ftwb(tao,go,patmpa)
            endif
        endif

        return
        end subroutine licoilab

! **********************************************************************
!  
        subroutine licoildy(tai,gi,wa,twi,ww,tsdyn,                     &
                           psychro,icountfl,lcoil,wcoil,dotube,thitube, &
                           thifin,ttubelen,as,ditube,afree,aext,        &
                           aratio,bratio,rifin,rofin,heifin,aflow,      &
                           hydiam,confin,contube,vcpfin,vcptube,        &
                           dtsdt,tao,go,two,qa,shr,bf,ho,rho,twbo)
!  
! **********************************************************************
!   SUBROUTINE:  Dynamic heating and cooling coil
!  
!   PURPOSE:     Calculate the dynamic outlet air and water conditions
!                for a heating or cooling coil 
! **********************************************************************
!   INPUTS                                                            
!   ======                                                            
!   tai     : inlet air dry bulb temperature                         (C)
!   gi      : inlet air humidity ratio                           (kg/kg)
!   wa      : dry air mass flow rate                              (kg/s)
!   twi     : inlet water temperature                                (C)
!   ww      : water mass flow rate                                (kg/s)
!   tsdyn   : effective coil surface temperature                     (C)
!                                                                     
!   OUTPUTS                                                           
!   =======                                                           
!   dtsdt   : derivative of effective coil surface temperature       (C)
!   tao     : outlet air dry bulb temperature                        (C)
!   go      : outlet air humidity ratio                          (kg/kg)
!   two     : outlet water temperature                              (C)
!   qa      : total heat transfer to the air                        (kW)
!   shr     : sensible heat ratio                                    (-)
!   bf      : coil by-pass factor                                    (-)
!   ho      : outlet air specific enthalpy                       (kJ/kg)
!   rho     : outlet air relative humidity                           (-)
!   twbo    : outlet air wet-bulb temperature                        (C)
!                                                                        
!   PARAMETERS                                                          
!   ==========                                                          
!   icountfl: 0 = cross-flow, 1 = counterflow, -1 = parallel flow    (-)
!   psychro : FALSE = no psychrometric output calcs, TRUE = calcs    (-)
!   lcoil   : length of finned section in direction of flow          (m)   
!   wcoil   : width of finned section in direction of flow           (m)   
!   dotube  : tube outside diameter                                  (m)
!   thitube : tube thickness                                         (m)
!   thifin  : fin thickness                                          (m)
!   ttubelen: total tube length                                      (m)
!   as      : net air side tube area                                (m2)
!   ditube  : tube inside diameter                                   (m)
!   afree   : exchanger minimum free-flow area on air side          (m2)
!   aext    : exchanger total heat transfer area on air side        (m2) 
!   aratio  : ratio of total fin area to total heat transfer area    (-)
!   bratio  : ratio of air side total area : water side intrnl area  (-)
!   rifin   : fin inside radius                                      (m)
!   rofin   : effective fin radius                                   (m)
!   heifin  : effective fin height                                   (m)
!   aflow   : flow area on water side                               (m2)
!   hydiam  : coil hydraulic diameter                                (m) 
!   confin  : thermal conductivity of fin material            (kW/(m.K)) 
!   contube : thermal conductivity of tube material           (kW/(m.K)) 
!   vcpfin  : volumetric heat capacity of fin material       (kJ/(m3.K)) 
!   vcptube : volumetric heat capacity of tube material      (kJ/(m3.K))  
!  
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
!   SUBROUTINES CALLED:  DRYOUTCO, TSHBYPAS, WETOUTCO
!   FUNCTIONS  CALLED:   AIRCOEFF, AIRFINRES, BYPASFAC, COILCAP,
!                        QAIRSIDE, QEFFNTU, QWATSIDE, WATERRES
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
!   INTERNAL VARIABLES
!   ==================
!
!   Material Properties
!   -------------------
!   rhow    : water density                                      (kg/m3)  
!   cpw     : water specific heat                            (kJ/(kg.K))  
!   cpa     : dry air specific heat                          (kJ/(kg.K))  
!   cpg     : water vapor specific heat                      (kJ/(kg/K)) 
!   tref    : reference temperature                                  (K)  
!   kw      : water thermal conductivity                      (kW/(m.K))  
!   viscw   : water dynamic viscosity                         (kg/(m.s))
!   rhoa    : density of air                                     (kg/m3)
!   rhow    : density of water                                   (kg/m3)
!   cpai    : inlet air mass specific heat                   (kJ/(kg.C)) 
!  
!   Heat Transfer Variables
!   -----------------------
!   va      : air velocity (referred to free area)                 (m/s)
!   bf      : coil by-pass factor                                    (-)
!   cair    : air capacity rate                               (kJ/(s.K))  
!   cwat    : water capacity rate                             (kJ/(s.K))
!   effdry  : dry coil effectiveness                                 (-)
!   hadry   : dry heat transfer coefficient on air side      (kW/(K.m2))
!   radry   : dry air side thermal resistance                  (K.m2/kW) 
!   rw      : water side thermal resistance                    (K.m2/kW)
!   rm      : metal thermal resistance                         (K.m2/kW) 
!   rtdry   : global thermal resistance                        (K.m2/kW)
!   rtwet   : equivalent wet global thermal resistance         (K.m2/kW)
!   antuwet : number of equivalent transfer units (wet conditions)   (-)   
!   cs      : average saturation specific heat               (kJ/(kg.K))
!   
!   Miscellaneous
!   -------------
!   vsmall  : threshold velocity for significant heat transfer     (m/s)
!   patmpa  : outlet air absolute pressure in Pascals for 
!             psychrometric routines                                (Pa)
!   tdpai   : inlet air dew point temperature                        (C)
! **********************************************************************
        use precision1
        implicit none
        real(kind=pp)   :: tai,gi,wa,twi,ww,tsdyn,                      &
                           lcoil,wcoil,dotube,thitube,         &
                           thifin,ttubelen,as,ditube,afree,aext,        &
                           aratio,bratio,rifin,rofin,heifin,aflow,      &
                           hydiam,confin,contube,vcpfin,vcptube,        &
                           dtsdt,tao,go,two,qa,shr,bf,ho,rho,twbo
        real(kind=pp)   :: t,csat,cpai,cair,cwat,va,vw,gair,hadry,  &
                           radry,rw,rm,rtdry,aircoeff,airfinres,    &
                           waterres,ftdew,qadry,qairside,qw,qwatside, &
                           rst1,rst1wet,bfwet,hi,fhair, &
                           rtwet,tdpai,w,twbi,ftwb,tsm,cs, &
                           hawet,rawet,cairsat,qawet,bfdry,hccoil, &
                           dpai,antuw,cst2,rst1dry,rst2,coilcap, &
                           fwphi,fphi,gsat,effect,a=100.0
        integer         :: icountfl
!  
        logical         :: psychro
! ** Thermophysical constants
        real(kind=pp)   :: rhoa = 1.2, rhow = 1000., cpa = 1., cpw = 4.18, &
                           cpg = 1.805, patmpa = 101325.
! ** Miscellaneous constants
        real(kind=pp)   :: vsmall = 1.e-3
! ** Saturation specific heat as a function of temperature
        csat(t)=1.786+2.27e-2*t+2.96565e-3*t*t
! ** Moist air specific heat and air and water capacity rates and velocities
        cpai    = cpa+cpg*gi
        cair    = cpai*wa
        cwat    = cpw*ww
        va      = wa/(rhoa*afree)
        vw      = ww/(rhow*aflow)
! ** Calculation of coil performance depends on whether there is
! ** significant flow on either side. The four possible cases are
! ** considered separately.
        if (va>vsmall.and.vw>vsmall) then
! ** General case: Significant flow on both air side and water side
! ** performance depends on whether the coil is wet or dry. 
! ** In this model, the coil is taken to be completely
! ** dry or completely wet, depending on which has the higher duty.
! ** Calculate thermal resistances referred to unit area on the air side
            gair   = wa/afree
            hadry  = aircoeff(tai,gair,hydiam,lcoil)
            radry  = airfinres(hadry,confin,thifin,rofin,rifin,&
                               heifin,aratio)
            rw     = waterres(ww,aflow,twi,ditube,bratio,wcoil)
            rm     = bratio*thitube/contube
            rtdry  = radry+rm+rw
! ** Inlet dew point temperature
            tdpai   = ftdew(gi,patmpa)
! ** Dynamic model - Calculate separate resistances on air and water sides,
! ** each referred to the inlet temperature (see printed documentation).
! ** Associate a heat capacity with the surface equal to that of the
! ** metal plus the water.
! ** Dry coil air side heat transfer. resistance includes the fins
            qadry  = qairside(tai,tsdyn,cair,aext,radry,rst1dry,&
                              bfdry)
! ** Compare surface temperature to inlet dew point temperature
            if (tsdyn>=tdpai) then
! ** Either coil is completely dry or dry (cooling) duty must exceed
! ** wet (cooling) duty. Calculate water side heat transfer for dry coil
                qw     = qwatside(aext,cair,cwat,rtdry,&
                                  icountfl,rst1dry,twi,tsdyn,rst2)
! ** Air side duty, resistance and by-pass factor
                qa     = qadry      
                rst1   = rst1dry
                bf     = bfdry
! ** Outlet conditions
                call dryoutco(tai,gi,qa,cair,wa,psychro,patmpa,tao,&
                              go,shr,ho,rho,twbo)
            else
! ** Coil may be wet or dry. Calculate air side heat transfer for wet coil
! ** use wet bulb as driving temperature and enhance air side heat transfer 
! ** coefficient by the ratio of the saturation specific heat to the 
! ** sensible specific heat. This is a good approximation to using enthalpy
! ** as the driving potential.
                twbi   = ftwb(tai,gi,patmpa)
! ** Saturation specific heat at average of driving temperatures
                tsm    = (twbi+twi)/2.0
                cs     = csat(tsm)
! ** Enhanced air side heat transfer coefficient and reduced resistance
                hawet  = hadry*cs/cpai
                rawet  = airfinres(hawet,confin,thifin,rofin,rifin,&
                           heifin,aratio)
! ** Air side heat transfer
                cairsat=cs*wa
                qawet  = qairside(twbi,tsdyn,cairsat,aext,rawet,&
                                  rst1wet,bfwet)
! ** Compare dry and wet cooling duties
                if (-qadry>=-qawet) then
! ** Dry coil - calculate outlet air condition and water side heat transfer
! ** Air side heat transfer, resistance and bypass factor
                    qa     = qadry
                    rst1   = rst1dry
                    bf     = bfdry
! ** Air outlet condition
                    call dryoutco(tai,gi,qa,cair,wa,psychro,patmpa,&
                                  tao,go,shr,ho,rho,twbo)
! ** Water side heat transfer
                    qw     = qwatside(aext,cair,cwat,rtdry,&
                                      icountfl,rst1dry,twi,tsdyn,rst2)
                else
! ** Wet coil - calculate outlet air condition and water side heat transfer
! ** Air side duty, resistance and by-pass factor
                    qa     = qawet
                    rst1   = rst1wet
                    bf     = bfwet
! ** Outlet enthalpy
                    hi     = fhair(tai,gi)
                    ho     = hi+qa/wa
! ** Outlet air dry bulb and humidity ratio
                    call wetoutco(tai,qa,ho,cair,bf,tsdyn,patmpa, &
                                  psychro,shr,tao,go,rho,twbo)
! ** Calculate water side resistance from total resistance and
! ** air side resistance under wet conditons
                    rtwet  = rawet+rw+rm
                    qw     = qwatside(aext,cairsat,cwat,rtwet, &
                                      icountfl,rst1wet,twi,tsdyn,rst2)
                endif
! ** End of wet/dry model
            endif
! ** Outlet water temperature
            two    = twi-qw/cwat
! ** Heat capacity of fins, tubes and water
            hccoil = coilcap(vcpfin,as,thifin,vcptube,&
                             ttubelen,dotube,thitube,ditube)
! ** Calculate rate of change of surface temperature from heat balance
            dtsdt  = (qw-qa)/hccoil
! ** End of significant flows on both sides
        elseif ((va>vsmall.and.vw<=vsmall)) then
! ** Water flow rate is negligible and air flow rate is significant
! ** calculate rate of change of surface temperature from heat  
! ** balance on coil heat capacity.
! ** Air side resistances and capacity rate
            gair   = wa/afree
            hadry  = aircoeff(tai,gair,hydiam,lcoil)
            radry  = airfinres(hadry,confin,thifin,rofin,rifin,&
                               heifin,aratio)
! ** Air side duty assuming dry coil
            qadry  = qairside(tai,tsdyn,cair,aext,radry,rst1dry,&
                              bfdry)
! ** Inlet dew point temperature
            tdpai   = ftdew(gi,patmpa)
            if (tsdyn>tdpai) then
! ** Coil is dry or dry coil heat flow exceeds wet coil heat flow
! ** Air side heat transfer, resistance and bypass factor
                qa     = qadry      
                rst1   = rst1dry
                bf     = bfdry
! ** Outlet air conditions
                call dryoutco(tai,gi,qa,cair,wa,psychro,patmpa,&
                              tao,go,shr,ho,rho,twbo)
! ** End dry model
            else
! ** Coil may be wet or dry
! ** Calculate air side heat transfer for wet coil
! ** Air side resistance
                twbi   = ftwb(tai,gi,patmpa)
                tsm    = (twbi+tsdyn)/2.0
                cs     = csat(tsm)
                hawet  = hadry*cs/cpai
                rawet  = airfinres(hawet,confin,thifin,rofin,rifin,&
                         heifin,aratio)
! ** Air side heat transfer
                cairsat=cs*wa
                qawet  = qairside(twbi,tsdyn,cairsat,aext,rawet,&
                                  rst1wet,bfwet)
                if (-qadry>=-qawet) then
! ** Coil effectively dry
! ** Air side heat transfer, resistance and bypass factor
                    qa     = qadry      
                    rst1   = rst1dry
                    bf     = bfdry
! ** Outlet air conditions
                    call dryoutco(tai,gi,qa,cair,wa,psychro,patmpa,&
                                  tao,go,shr,ho,rho,twbo)
                else
! ** Wet coil - calculate outlet air condition and water side heat transfer
! ** Air side duty, resistance and by-pass factor
                    qa     = qawet
                    rst1   = rst1wet
                    bf     = bfwet
! ** Outlet enthalpy
                    hi     = fhair(tai,gi)
                    ho     = hi+qa/wa
! ** Outlet air dry bulb and humidity ratio
                    call wetoutco(tai,qa,ho,cair,bf,tsdyn,patmpa,&
                                  psychro,shr,tao,go,rho,twbo)
! ** End wet model
                endif
! ** End wet/dry model
            endif
! ** No water side heat flow - outlet water at temperature of surface
            two    = tsdyn
! ** Heat capacity of fins, tubes and water
            hccoil = coilcap(vcpfin,as,thifin,vcptube,&
                             ttubelen,dotube,thitube,ditube)
! ** Heat balance on coil heat capacity
            dtsdt  = -qadry/hccoil
! ** End negligible water flow, significant air flow
        elseif (va.le.vsmall.and.vw>vsmall) then
! ** Air flow rate negligible and water flow rate significant.
! ** Calculate rate of change of surface temperature from heat
! ** balance on coil heat capacity (this yields a particularly poor 
! ** approximation to the dynamics when the air flow is negligible, but
! ** continues to maintain an energy balance until the air flow becomes
! ** significant again).
! ** Water side duty - Calculate from ntu-effectiveness relationship
! ** assuming infinite capacity rate for the coil surface
            rw     = waterres(ww,aflow,twi,ditube,bratio,wcoil)
            antuw  = aext/(rw*ww*cpw)      
            cst2   = (1.-exp(-antuw))*cwat
            qw     = cst2*(twi-tsdyn)
! ** Heat capacity of fins, tubes and water
            hccoil = coilcap(vcpfin,as,thifin,vcptube,&
                             ttubelen,dotube,thitube,ditube)
! ** Heat balance on coil heat capacity
            dtsdt  = qw/hccoil
! ** Outlet air is in equilibrium with surface
            tao    = tsdyn
! ** Heat balance on water side
            two    = twi-qw/cwat
! ** Outlet humidity ratio is equal to inlet value or saturated value,
! ** whichever is lower

            go     = min(gi,fwphi(tao,a,patmpa))
!            go     = min(gi,fwphi(tao,100.,patmpa))
            qa     = 0.
            effect = 1.
            bf     = 0.
! ** Sensible heat ratio is indeterminate
            shr    = 0.
            if (psychro) then
                ho     = fhair(tao,go)
                rho    = fphi(tao,go,patmpa)
                twbo   = ftwb(tao,go,patmpa)
            endif
! ** End of case of negligible air flow
        else
! ** Both flow rates negligible
! ** set outlet conditions equal to effective surface conditions
! ** as reasonable starting point for solution when one or other flow
! ** rate becomes significant
            tao    = tsdyn
            gsat   = fwphi(tao,a,patmpa)
!            gsat   = fwphi(tao,100.0,patmpa)
            go     = min(gi,gsat)
            two    = tsdyn
! ** No heat transfer on either side - no change in surface temperature
            dtsdt  = 0.
            qa     = 0.
            effect = 1.
            bf     = 0.
            shr    = 0.
            if (psychro) then
                ho     = fhair(tao,go)
                rho    = fphi(tao,go,patmpa)
                twbo   = ftwb(tao,go,patmpa)
            endif
        endif
!
        return
        end subroutine licoildy

! **********************************************************************
!  
        subroutine licoilss(tai,gi,wa,twi,ww,                      &
                           psychro,icountfl,lcoil,wcoil,thitube,   &
                           thifin,ditube,afree,aext,               &
                           aratio,bratio,rifin,rofin,heifin,aflow, &
                           hydiam,confin,contube,                  &
                           ts,tao,go,two,qa,shr,effect,bf,ho,rho,twbo)
!  
! **********************************************************************
!   SUBROUTINE:  Steady state Liege heating and cooling coil
!  
!   PURPOSE:     Calculate the steady state outlet air and water
!                conditions for a heating or cooling coil 
! **********************************************************************
!   INPUTS                                                            
!   ======                                                            
!   tai      : inlet air dry bulb temperature                        (C)
!   gi       : inlet air humidity ratio                          (kg/kg)
!   wa       : dry air mass flow rate                             (kg/s)
!   twi      : inlet water temperature                               (C)
!   ww       : water mass flow rate                               (kg/s)
!                                                                     
!   OUTPUTS                                                           
!   =======                                                           
!   ts       : effective coil surface temperature                    (C)
!   tao      : outlet air dry bulb temperature                       (C)
!   go       : outlet air humidity ratio                         (kg/kg)
!   two      : outlet water temperature                              (C)
!   qa       : total heat transfer to the air                       (kW)
!   shr      : sensible heat ratio                                   (-)
!   effect   : heat exchanger effectiveness                          (-)
!   bf       : coil by-pass factor                                   (-)
!   ho       : outlet air specific enthalpy                      (kJ/kg)
!   rho      : outlet air relative humidity                          (-)
!   twbo     : outlet air wet-bulb temperature                       (C)
!                                                                        
!   PARAMETERS                                                          
!   ==========                                                          
!   icountfl : 0 = cross flow, 1 = counterflow, -1 = parallel flow   (-)
!   psychro  : FALSE = no psychrometric output calcs, TRUE = calcs   (-)
!   lcoil    : length of finned section in direction of flow         (m)   
!   wcoil    : width of finned section in direction of flow          (m)   
!   thitube  : tube thickness                                        (m)
!   thifin   : fin thickness                                         (m)
!   ditube   : tube inside diameter                                  (m)
!   afree    : exchanger minimum free-flow area on air side         (m2)
!   aext     : exchanger total heat transfer area on air side       (m2) 
!   aratio   : ratio of total fin area to total heat transfer area   (-)
!   bratio   : ratio of air side total area : water side intrnl area (-)
!   rifin    : fin inside radius                                     (m)
!   rofin    : effective fin radius                                  (m)
!   heifin   : effective fin height                                  (m)
!   aflow    : flow area on water side                              (m2)
!   hydiam   : coil hydraulic diameter                               (m) 
!   confin   : thermal conductivity of fin material           (kW/(m.K)) 
!   contube  : thermal conductivity of tube material          (kW/(m.K)) 
!  
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
!   SUBROUTINES CALLED:  DRYOUTCO, TSHBYPAS, WETOUTCO
!   FUNCTIONS  CALLED:   AIRCOEFF, AIRFINRES, BYPASFAC, COILCAP,
!                        QAIRSIDE, QEFFNTU, QWATSIDE, WATERRES
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
!   INTERNAL VARIABLES                                           
!   ==================                                         
!  
!   Material Properties
!   -------------------
!   rhow    : water density                                      (kg/m3)  
!   cpw     : water specific heat                            (kJ/(kg.K))  
!   cpa     : dry air specific heat                          (kJ/(kg.K))  
!   cpg     : water vapor specific heat                      (kJ/(kg/K)) 
!   tref    : reference temperature                                  (K)  
!   kw      : water thermal conductivity                      (kW/(m.K))  
!   viscw   : water dynamic viscosity                         (kg/(m.s))
!   rhoa    : density of air                                     (kg/m3)
!   rhow    : density of water                                   (kg/m3)
!   cpai    : inlet air mass specific heat                   (kJ/(kg.C)) 
!  
!   Heat transfer Variables
!   -----------------------
!   va      : air velocity (referred to free area)                 (m/s)
!   bfdry   : coil by-pass factor                                    (-)
!   cair    : air capacity rate                               (kJ/(s.K))  
!   cwat    : water capacity rate                             (kJ/(s.K))
!   antuw   : number of transfer units on water side                 (-)
!   effdry  : dry coil effectiveness                                 (-)
!   hadry   : dry heat transfer coefficient on air side      (kW/(K.m2))
!   radry   : dry air side thermal resistance                  (K.m2/kW) 
!   rw      : water side thermal resistance                    (K.m2/kW)
!   rm      : metal thermal resistance                         (K.m2/kW) 
!   rtdry   : global thermal resistance                        (K.m2/kW)
!   rtwet   : equivalent wet global thermal resistance         (K.m2/kW)
!   antuwet : number of equivalent transfer units (et conditions)    (-)
!   cs      : average saturation specific heat               (kJ/(kg.K))
!   
!   Miscellaneous
!   -------------
!   vsmall  : threshold velocity for significant heat transfer     (m/s)
!   patmpa  : outlet air absolute pressure in Pascals for 
!             psychrometric routines                                (Pa)
!   tdpai   : inlet air dew point temperature                        (C)
! **********************************************************************

        use precision1
        implicit none

        real(kind=pp)   :: tai,gi,wa,twi,ww,                       &
                           lcoil,wcoil,thitube,           &
                           thifin,ditube,afree,aext,               &
                           aratio,bratio,rifin,rofin,heifin,aflow, &
                           hydiam,confin,contube,                  &
                           ts,tao,go,two,qa,shr,effect,bf,ho,rho,twbo
        real(kind=pp)   :: t,csat,cpai,cair,cwat,va,vw,gair,hadry,  &
                           radry,rw,rm,rtdry,aircoeff,airfinres,    &
                           waterres,ftdew,qadry,qairside,qw,qwatside, &
                           rst1,rst1wet,bfwet,hi,fhair, &
                           rtwet,tdpai,w,twbi,ftwb,tsm,cs, &
                           hawet,rawet,cairsat,qawet,bfdry,hccoil, &
                           dpai,antuw,cst2,rst1dry,rst2,coilcap, &
                           fwphi,fphi,gsat,qdry,qeffntu,effdry, &
                           bypasfac,effwet,qwet
        integer         :: icountfl
        logical         :: psychro
! ** Thermophysical constants
        real(kind=pp)   :: rhoa = 1.2, rhow = 1000., cpa = 1., cpw = 4.18, &
                           cpg = 1.805, patmpa = 101325.
! ** Miscellaneous constants
        real(kind=pp)   :: vsmall = 1.e-3
! ** Saturation specific heat as a function of temperature
        csat(t)=1.786+2.27e-2*t+2.96565e-3*t*t
! ** Moist air specific heat and air and water capacity rates and velocities
        cpai    = cpa+cpg*gi
        cair    = cpai*wa
        cwat    = cpw*ww
        va      = wa/(rhoa*afree)
        vw      = ww/(rhow*aflow)
! ** Calculation of coil performance depends on whether there is
! ** significant flow on either side. The four possible cases are
! ** considered separately.
        if (va>vsmall.and.vw>vsmall) then
! ** General case: Significant flow on both air side and water side
! ** the performance depends on whether the coil is wet or dry. 
! ** In this model, the coil is taken to be completely dry or
! ** completely wet, depending on which has the higher duty.
! ** Calculate thermal resistances referred to unit area on the air side
            gair   = wa/afree
            hadry  = aircoeff(tai,gair,hydiam,lcoil)
            radry  = airfinres(hadry,confin,thifin,rofin,rifin,&
                               heifin,aratio)
            rw     = waterres(ww,aflow,twi,ditube,bratio,wcoil)
            rm     = bratio*thitube/contube
            rtdry  = radry+rm+rw
! ** Inlet dew point temperature
            tdpai   = ftdew(gi,patmpa)
! ** Steady state - Calculate dry coil duty from ntu-effectiveness
! ** relationship corresponding to coil configuration
            qdry = qeffntu(cair,cwat,tai,twi,aext,rtdry,icountfl,&
                           effdry)
! ** Test for completely dry coil by comparing inlet dew point and inlet
! ** water temperature
            if (twi>=tdpai) then
! ** Completely dry coil - Calculate outlet conditions and return
! ** Air side duty, by-pass factor, sensible heat ratio and effectiveness
                qa     = qdry
                bf     = bypasfac(cair,aext,radry)
                shr    = 1.0
                effect = effdry
! ** Outlet air conditions
                call dryoutco(tai,gi,qa,cair,wa,psychro,patmpa,tao,&
                              go,shr,ho,rho,twbo)
! ** Surface temperature
                ts     = (tao-bf*tai)/(1.-bf)
            else
! ** Coil is (partly) wet - Calculate completely wet duty and compare with
! ** completely dry duty, taking the larger. use wet bulb as the driving
! ** temperature and enhance the air side heat transfer coefficient by the
! ** ratio of the saturation specific heat to the sensible specific heat.
! ** This is a good approximation to the use of enthalpy as the driving
! ** potential.
                twbi   = ftwb(tai,gi,patmpa)
! ** Saturation specific heat at average of driving temperatures used to
! ** calculate enhanced air side heat transfer coefficient and reduced 
! ** resistance
                tsm    = (twbi+twi)/2.0
                cs     = csat(tsm)
                hawet  = hadry*cs/cpai
                rawet  = airfinres(hawet,confin,thifin,rofin,rifin,&
                         heifin,aratio)
                rtwet  = rawet+rw+rm
! ** Capacity rate of saturated air
                cairsat=cs*wa
! ** Calculate wet duty from appropriate ntu-effectiveness relationship
                qwet   = qeffntu(cairsat,cwat,twbi,twi,aext,&
                                 rtwet,icountfl,effwet)
! ** Compare completely dry and completely wet duties.
! ** duties are assumed to be negative (cooling)
                if (-qdry>=-qwet) then
! ** Dry coil duty is greater
! ** air side duty, by-pass factor, sensible heat ratio and effectiveness 
                    qa     = qdry
                    bf     = bypasfac(cair,aext,radry)
                    shr    = 1.0
                    effect = effdry
! ** Calculate outlet conditions from dry duty
                    call dryoutco(tai,gi,qa,cair,wa,psychro,patmpa,&
                                  tao,go,shr,ho,rho,twbo)
! ** Surface temperature
                    ts     = (tao-bf*tai)/(1.-bf)
                else
! ** Wet coil duty is greater - Calculate outlet conditions from wet duty
! ** air side duty and effectiveness
                    qa     = qwet
                    effect = effwet
! ** Outlet enthalpy, bypass factor and estimate of surface temperature
                    call tshobf(tai,gi,qa,wa,cairsat,aext,rawet,&
                                patmpa,ts,ho,bf)
! ** Outlet air conditions
                    call wetoutco(tai,qa,ho,cair,bf,ts,patmpa,&
                                  psychro,shr,tao,go,rho,twbo)
                endif
! ** End of (partly) wet model
            endif
! ** Outlet water temperature
            two    = twi-qa/cwat
! ** End of significant flows on both sides
        elseif ((va>vsmall.and.vw.le.vsmall)) then
! ** Water flow rate is negligible and air flow rate is significant
! ** Set outlet and surface temperatures to inlet air
! ** temperature. No duty, so no change in humidity ratio
            qa     = 0.
            effect = 0.
! ** By-pass factor and sensible heat ratio are indeterminate
            bf     = 0.
            shr    = 0.
            ts     = tai
! ** Outlet conditions
            two    = tai
            tao    = tai
            go     = gi
            effect = 1.
            if (psychro) then
! ** Auxilliary outputs
                ho     = fhair(tao,go)
                rho    = fphi(tao,go,patmpa)
                twbo   = ftwb(tao,go,patmpa)
            endif
! ** End negligible water flow, significant air flow
        elseif (va.le.vsmall.and.vw>vsmall) then
! ** Air flow rate negligible and water flow rate significant.
! ** Set outlet and surface temperatures to inlet water
! ** temperature
            ts     = twi
            tao    = twi
! ** No duty, so outlet water temperature equals inlet water temperature
            two    = twi
! ** Outlet humidity ratio is equal to inlet value or saturated value,
! ** whichever is lower
            go     = min(gi,fwphi(tao,100.,patmpa))
            qa     = 0.
            effect = 1.
            bf     = 0.
! ** Sensible heat ratio is indeterminate
            shr    = 0.
            if (psychro) then
                ho     = fhair(tao,go)
                rho    = fphi(tao,go,patmpa)
                twbo   = ftwb(tao,go,patmpa)
            endif
! ** End of case of negligible air flow
        else
! ** Both flow rates negligible
! ** Set outlet conditions equal to inlet conditions
! ** as reasonable starting point for solution when one or other flow
! ** rate becomes significant.
            tao    = tai
            go     = gi
            two    = twi
            ts     = (tai+twi)/2.0
            qa     = 0.
            effect = 0.
            bf     = 0.
            shr    = 0.
            if (psychro) then
                ho     = fhair(tao,go)
                rho    = fphi(tai,gi,patmpa)
                twbo   = ftwb(tao,go,patmpa)
            endif
! ** End both flows negligible
        endif
!
        return
        end subroutine licoilss

! **********************************************************************
!
!   LIBRARY:        Routines for use in coil models
!  
!   DEVELOPER:      Philip Haves and Li Mei
!                   Loughborough University of Technology
!
!   REVISION DATE:       
!
!   REFERENCE:      ASHRAE 825-RP Final Report
!
!=======================================================================
!
        real function effectiv(iconfig,crr,antu)
!
! **********************************************************************
!  
!   PURPOSE:  Calculate the effectiveness of different heat exchanger
!             configurations
!  
! **********************************************************************
!   INPUTS
!   ======
!   iconfig : configuration: +1 = counterflow                        (-)
!                             0 = cross-flow (both sides unmixed)    (-)
!                            -1 = parallel flow                      (-)
!                            +2 = cross-flow (maximum fluid unmixed) (-)
!                            -2 = cross-flow (minimum fluid unmixed) (-)
!   crr     : ratio of minimum to maximum capacity rate              (-)
!   antu    : total number of transfer units                         (-)
!  
!   OUTPUT
!   ======
!   effectiv: heat exchanger effectiveness                
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)             :: crr,antu,omcrr,e0,eta,a, &
                                     e1,b,e2
        integer                   :: iconfig
        real(kind=pp)             :: small = 1.e-15
! ** Trap special cases
! ** Zero ntu; no duty
        if (antu<small) effectiv=0.
! ** Zero capacity rate ratio; effectiveness is independent of
! ** configuration
        if (crr<small) effectiv=1.-exp(-antu)
! ** Calculate effectiveness for different configurations
        if (iconfig==1) then
! ** Counterflow
            omcrr=1.-crr
            if (abs(omcrr)>1.0e-5) then
                e0=exp(-antu*(omcrr))
                effectiv=(1.-e0)/(1.-crr*e0)
            else
! ** Capacity rates approx equal, avoid indeterminacy
                effectiv=antu/(1. + antu)
            endif
        elseif (iconfig==-1) then
! ** Parallel flow
            e0=exp(-antu*(1.+crr))
            effectiv=(1.-e0)/(1.+crr)
        elseif (iconfig==0) then
! ** Cross-flow (both sides unmixed)
            eta=antu**0.22
            a=crr*antu/eta
            if(a<20.) then
                e1=exp(-a)
            else
                e1=0.
            endif
            b=eta*(1.-e1)/crr
            if(b<20.) then
                e2=exp(-b)
            else
                e2=0.
            endif
            effectiv=1.-e2
        elseif (iconfig==2) then
! ** Cross-flow - maximum fluid unmixed
            effectiv=1.-exp(-(1.-exp(-antu*crr))/crr)
        elseif (iconfig==-2) then
! ** Cross-flow - minimum fluid unmixed
            effectiv=(1.-exp(-crr*(1.-exp(-antu))))/crr
        else
            stop 'effectiv: iconfig out of range'
        endif
!
        return
        end function effectiv
!
!=======================================================================
!
        real function aircoeff(tai,gair,hydiam,lcoil)
!
! **********************************************************************
!  
!   PURPOSE:    Calculate the convective heat transfer coefficient on the
!               air side
!  
!   WRITTEN BY: Laboratoire de Thermodynamique, University of Liege
!  
! **********************************************************************
!   INPUTS
!   ======
!   tai     : inlet air dry bulb temperature                         (C)
!   gair    : mass flux                                      (kg/(s.m2))
!   hydiam  : hydraulic diameter                                     (m)
!   lcoil   : length of coil in direction of air flow                (m)
!  
!   OUTPUT                                                            
!   ======
!   aircoeff: heat transfer coefficient on air side          (kW/(K*m2))
!                                                                       
!   INTERNAL VARIABLES
!   ==================
!   prair   : Prantl number for air                                  (-)
!   kair    : conductivity of air                              (W/(K.m))
!   viscmu  : air dynamic viscosity                           (kg/(m*s)) 
!   reyair  : Reynolds number                                        (-)
!   nusair  : Nusselt number                                         (-)
!  
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)             :: tai,gair,hydiam, &
                                     viscmu,reyair,gz
        real(kind=pp)             :: lcoil,nusair
        real(kind=pp)             :: prair = 0.71, kair = .0257e-3
! ** Calculate heat transfer coefficent from an empirical correlation
! ** for the Nusselt number
! ** viscosity and Reynolds number
        viscmu = 1.458e-6*((273.16+tai)**1.5/(383.4+tai))
        reyair = gair*hydiam/viscmu
        gz     = reyair*prair*hydiam/lcoil
! ** Nusselt number
        nusair = 3.66+(0.0668*gz/(1+0.04*gz**(2/3)))
! ** Heat transfer coefficient
        aircoeff = nusair*kair/hydiam
!
        return
        end function aircoeff
!
!=======================================================================
!
        real function airfinres(ha,confin,thifin,rofin,rifin,heifin,aratio)
!
! **********************************************************************
!  
!   PURPOSE:  Calculate the heat exchange resistance on the air side of
!             a fin-tube heat exchanger
!  
! **********************************************************************
!   INPUTS
!   ======
!   ha      : heat transfer coefficient on air side          (kW/(K*m2))
!   confin  : thermal conductivity of fin material            (kW/(m.K)) 
!   thifin  : fin thickness                                          (m)
!   rofin   : effective fin radius                                   (m)
!   rifin   : fin inside radius                                      (m)
!   heifin  : effective fin height                                   (m)
!   aratio  : ratio of total fin area to total heat transfer area    (-)
!  
!   output                                                            
!   ======
!   airfinres : air side thermal resistance                    (K*m2/kW)
!                                                                       
!   INTERNAL VARIABLES
!   ==================
!   efin    : fin efficiency                                         (-)
!   efingl  : fin global effectiveness                               (-)
!  
! **********************************************************************

        use precision1
        implicit none
        real(kind=pp)             :: ha,confin,thifin,rofin,rifin, &
                                     heifin,aratio
        real(kind=pp)             :: x1, x2, x3, efin, efingl

! ** Fin efficiency
        x1     = sqrt(2.*ha/confin/thifin)
        x2     = 1.+0.35*log(rofin/rifin)
        x3     = x1*x2*heifin
        efin   = tanh(x3)/x3
! ** Fin effectiveness
        efingl = 1.-aratio*(1.-efin)
! ** Air side resistance
        airfinres = 1./(efingl*ha)
!
        return
        end function airfinres
!
!=======================================================================
!
        real function viscwat(t)
!
! **********************************************************************
!  
!   PURPOSE:  Calculate the dynamic viscosity of water as a function of
!             temperature by interpolating values given in the Handbook 
!             of Physics and Chemistry                   
!  
! **********************************************************************
!   INPUT
!   =====
!   t       : temperature                                            (C)
!  
!   OUTPUT                                                            
!   ======
!   viscwat : dynamic viscosity of water                          (Pa.s)
!   
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)             :: t

! ** Interpolate Handbook values
        if (t>=0.0 .and. t<20.0) then
            viscwat = (1.787e-3*(20.0-t) + 1.002e-3*t)/20.0
        elseif (t>=20.0 .and. t<40.0) then
            viscwat = (1.002e-3*(40.0-t) + 0.6529e-3*(t-20.0))/20.0
        elseif (t>=40.0 .and. t<60.0) then
            viscwat = (0.6529e-3*(60.0-t) + 0.4665e-3*(t-40.0))/20.0
        elseif (t>=60.0 .and. t<85.0) then
            viscwat = (0.4665e-3*(80.0-t) + 0.3540e-3*(t-60.0))/20.0
        else 
            write(*,*) 'viscwat: temperature out of range ',t
        endif
!
        return
        end function viscwat
!
!=======================================================================
!
      real function waterres(ww,aflow,twi,ditube,bratio,wcoil)
!
! **********************************************************************
!  
!   PURPOSE:     Calculate water-side heat transfer coefficient and 
!                thermal resistance
!  
!   WRITTEN BY:  Laboratoire de Thermodynamique, University of Liege
!  
!   MODIFIED BY: P. Haves, 7.3.92                                      
!                                                                     
!             1) Reynolds number break points changed from 2300 -> 2245.7   
!                and 6000 -> 6109.3                                       
!             2) Take greater of laminar and turbulent Nusselt numbers      
!  
! **********************************************************************
!   INPUTS
!   ======
!   ww      : water mass flow rate                                (kg/s)
!   aflow   : flow area on water side                               (m2)
!   twi     : inlet water temperature                                (C)
!   ditube  : tube inside diameter                                   (m)
!   bratio  : ratio of air side total area : water side intrnl area  (-)
!   wcoil   : width of finned section in direction of flow           (m)   
!  
!   OUTPUT
!   ======
!   waterres: water-side thermal resistance                    (K*m2/kW)
!  
!   INTERNAL VARIABLES
!   ==================
!   prwat   : Prantl number for water                                (-)
!   kwat    : conductivity of water                              W/(K.m)
!   visc    : dynamic viscosity as a functiion of temperature (kg/(m*s)) 
!   reywat  : Reynolds number                                        (-)
!   nuwatlam: laminar Nusselt number                                 (-)
!   nuwattur: turbulent Nusselt number                               (-)
!   hw      : heat transfer coefficient on water side        (kW/(K*m2))
!  
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)           :: ww,aflow,twi,ditube,bratio,wcoil
        real(kind=pp)           :: nuwat,nuwatlam,nuwattur,kwat,gwat,visc, &
                                   reywat,prwat,pw,hw,viscwat
        real(kind=pp)           :: cpw = 4.18

! ** Calculate heat transfer coefficient from the larger of the laminar
! ** and turbulent Nusselt numbers
! ** Reynolds and Prantl numbers
        gwat    = ww/aflow
        visc    = viscwat(twi)
        reywat  = gwat*ditube/visc
        kwat    = 5.53e-4+twi*2.41e-6-twi*twi*1.32e-8+ &
                  twi*twi*twi*1.71e-11
        prwat   = cpw*visc/kwat
! ** Laminar Nusselt number
        nuwatlam = 1.86*(reywat*prwat*ditube/wcoil)**0.33
! ** Empirical correlations for the turbulent Nusselt number
        if (reywat<2245.7) then
            pw       = 1.8*prwat**0.3-0.8
            nuwattur = 0.013978*(reywat**0.8-36.4)*pw
        elseif (reywat>=2245.7.and.reywat<6109.3) then
            pw       = 1.8*prwat**0.3-0.8
            nuwattur = 0.027384*(reywat**0.8-253.48)*pw     
        else
            pw       = 1.8*prwat**0.3-0.8
            nuwattur = 0.02264*(reywat**0.8-82.69)*pw
        endif
! ** Heat transfer coefficient and thermal resistance
        nuwat    = max(nuwatlam,nuwattur) 
        hw       = nuwat*kwat/ditube
        waterres = bratio/hw
!
        return
        end function waterres
!
!======================================================================
!
      real function bypasfac(cair,aext,ra)
!
! **********************************************************************
!   
!   PURPOSE:    Calculate the by-pass factor on the air side.
!    
! **********************************************************************
!   INPUTS
!   ======
!   cair    : air capacity rate                                   (kW/K)   
!   aext    : area for air side heat transfer                       (m2)
!   ra      : air side thermal resistance                      (K.m2/kW)
!  
!   OUTPUT
!   ======
!   bypasfac: by-pass factor                                         (-)
!  
! **********************************************************************
!
      use precision1
      implicit none
      real(kind=pp)           :: cair,aext,ra
      bypasfac = exp(-aext/(cair*ra))
!
      return
      end function bypasfac
!
!=======================================================================
!
      subroutine capratio(cair,cwat,crr,cmin)
!
! **********************************************************************
!       
!   PURPOSE:    Calculate the capacity rate ratio and the minimum
!               capacity rate from the capacity rates on the air and 
!               water side
!    
! **********************************************************************
!   INPUTS
!   ======
!   cair    : air capacity rate                                   (kW/K)     
!   cwat    : water capacity rate                                 (kW/K)   
!  
!   OUTPUTS
!   =======
!   crr     : capacity rate ratio                                    (-)
!   cmin    : minimum capacity rate                               (kW/K)
!    
! **********************************************************************
!
      use precision1
      implicit none
      real(kind=pp)           :: cair,cwat,crr,cmin

! ** Minimum capacity rate and capacity rate ratio
        if (cair<=cwat) then
! ** Air is minimum fluid
           crr = cair/cwat
           cmin = cair
        else      
! ** Water is minimum fluid
           crr = cwat/cair
           cmin = cwat
        endif
!
        return
        end subroutine capratio
!
!=======================================================================
!
      subroutine limithr(tao,go,ho,tai,cair,qa,patmpa,shr)
!
! **********************************************************************
!  
!   PURPOSE:    Check that the humidity ratio at the outlet of a coil
!               does not exceed saturation value and recalculate outlet 
!               conditions if necessary
!  
! **********************************************************************
!   INPUTS
!   ======
!   tao     : outlet air dry bulb temperature                        (C)
!   go      : outlet air humidity ratio                          (kg/kg)
!   ho      : outlet air specific enthalpy                       (kJ/kg)
!   tai     : inlet air dry bulb temperature                         (C)
!   cair    : air capacity rate                               (kJ/(s.K))  
!   qa      : total heat transfer to the air                        (kW)
!   patmpa  : outlet air absolute pressure in Pascals for 
!             psychrometric routines                                (Pa)
!  
!   OUTPUTS
!   =======
!   tao     : outlet air dry bulb temperature                        (C)
!   go      : outlet air humidity ratio                          (kg/kg)
!   shr     : sensible heat ratio                                    (-)
!  
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)           :: tao,go,ho,tai,cair,qa,patmpa,shr, &
                                   gsat,fwphi,ftair,qasens,a=100.0
! ** Compare actual and saturation humidity ratios
        gsat   = fwphi(tao,a,patmpa)
        if (go>gsat) then
! ** Limit humidity ratio to saturation value while conserving enthalpy
            go     = gsat
            tao    = ftair(ho,go) 
            qasens = (tao-tai)/(cair)
            shr    = qasens/qa
        endif   
!
        return
        end subroutine limithr
!
!=======================================================================
!
        subroutine wetoutco(tai,qa,ho,cair,bf,ts,patmpa,psychro, &
                            shr,tao,go,rho,twbo)
!
! **********************************************************************
!  
!   PURPOSE:    Calculate the condition of the outlet air from a wet 
!               coil, limiting the outlet humidity ratio to the 
!               saturation value while conserving enthalpy
!  
! **********************************************************************
!   INPUTS
!   ======
!   tai     : inlet air dry bulb temperature                         (C)
!   qa      : total heat transfer to the air                        (kW)
!   ho      : outlet air specific enthalpy                       (kJ/kg)
!   cair    : air capacity rate                               (kJ/(s.K))  
!   bf      : coil by-pass factor                                    (-)
!   ts      : effective coil surface temperature                     (C)
!   patmpa  : outlet air absolute pressure in Pascals for 
!             psychrometric routines                                (Pa)
!   psychro : FALSE = no psychrometric output calcs, TRUE = calcs    (-)
!   
!   OUTPUTS
!   =======
!   shr     : sensible heat ratio                                    (-)
!   tao     : outlet air dry bulb temperature                        (C)
!   go      : outlet air humidity ratio                          (kg/kg)
!   rho     : outlet air relative humidity                           (-)
!   twbo    : outlet air wet-bulb temperature                        (C)
!   
! ***********************************************************************

        use precision1
        implicit none
        real(kind=pp)           :: tai,qa,ho,cair,bf,ts,patmpa, &
                                   shr,tao,go,rho,twbo,teo,qasens,fwha, &
                                   fphi,ftwb
        logical                 :: psychro

! ** Wet coil outlet conditions
! ** Estimate of outlet dry bulb temperature
        teo    = bf*(tai-ts)+ts
! ** Sensible heat ratio - limited to 1
        shr    = cair*(teo-tai)/qa
        if (shr>1.0) shr = 1.0
        qasens = shr*qa
! ** Outlet air dry bulb and humidity ratio
        tao    = tai+qasens/cair 
        go     = fwha(tao,ho)
! ** Check outlet air not supersaturated
        call limithr(tao,go,ho,tai,cair,qa,patmpa,shr)
        if (psychro) then
! ** Auxilliary outputs
            rho     = fphi(tao,go,patmpa)
            twbo    = ftwb(tao,go,patmpa)
        endif
!
        return
        end subroutine wetoutco
!
!=======================================================================
!
        subroutine dryoutco(tai,gi,qa,cair,wa,psychro,patmpa,tao,go,shr, &
                            ho,rho,twbo)
!
! **********************************************************************
!  
!   PURPOSE:  Calculate the outlet air condition for a dry coil
!             
! **********************************************************************
!  
!   INPUTS
!   ======
!   tai     : inlet air dry bulb temperature                         (C)
!   gi      : inlet air humidity ratio                           (kg/kg)
!   qa      : total heat transfer to the air                        (kW)
!   cair    : air capacity rate                               (kJ/(s.K))  
!   wa      : dry air mass flow rate                              (kg/s)
!   psychro : FALSE = no psychrometric output calcs, TRUE = calcs    (-)
!   patmpa  : outlet air absolute pressure in Pascals for 
!             psychrometric routines                                (Pa)
!  
!   OUTPUTS
!   =======
!   tao     : outlet air dry bulb temperature                        (C)
!   go      : outlet air humidity ratio                          (kg/kg)
!   shr     : sensible heat ratio                                    (-)
!   ho      : outlet air specific enthalpy                       (kJ/kg)
!   rho     : outlet air relative humidity                           (-)
!   twbo    : outlet air wet-bulb temperature                        (C)
!  
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)           :: tai,gi,qa,cair,wa,patmpa,tao,go,shr, &
                                   ho,rho,twbo,hi,fhair,fphi,ftwb
        logical                 :: psychro

! ** Outlet air temperature
        tao=tai+qa/cair
! ** No latent duty
        go=gi
        shr=1.0
        if (qa<0.0) then
! ** Cooling - Check outlet air not supersaturated
            hi=fhair(tai,gi)
            ho=hi+qa/wa
            call limithr(tao,go,ho,tai,cair,qa,patmpa,shr)
        else
! ** Heating - Calculate outlet enthalpy if required
            if (psychro) ho=fhair(tao,go)
        endif
        if (psychro) then
! ** Auxilliary outputs
            rho=fphi(tao,go,patmpa)
            twbo=ftwb(tao,go,patmpa)
        endif            
        return
        end subroutine dryoutco
!
!=======================================================================
!
        subroutine tshobf(tai,gi,qa,wa,cairsat,aext,rawet,patmpa, &
                          ts,ho,bf)
! **********************************************************************
!  
!   PURPOSE:  Calculate the effective surface temperature, outlet 
!             enthalpy and bypass factor for a wet coil from the inlet 
!             conditions and the duty
!  
! **********************************************************************
!   INPUTS
!   ======
!   tai     : inlet air dry bulb temperature                         (C)
!   gi      : inlet air humidity ratio                           (kg/kg)
!   qa      : total heat transfer to the air                        (kW)
!   wa      : dry air mass flow rate                              (kg/s)
!   cairsat : saturated air capacity rate                     (kJ/(s.K))  
!   aext    : area for air side heat transfer                       (m2)
!   rawet   : air side thermal resistance under wet conditions (K.m2/kW)
!   patmpa  : outlet air absolute pressure in Pascals for 
!             psychrometric routines                                (Pa)
!  
!   OUTPUTS
!   =======
!   ts      : effective coil surface temperature                     (C)
!   ho      : outlet air specific enthalpy                       (kJ/kg)
!   bf      : coil by-pass factor                                    (-)
!  
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)           :: tai,gi,qa,wa,cairsat,aext,rawet,patmpa, &
                                   ts,ho,bf,hi,fhair,bypasfac,hs,ftsat
! ** Outlet enthalpy
        hi=fhair(tai,gi)
        ho=hi+qa/wa
! ** By-pass factor
        bf=bypasfac(cairsat,aext,rawet)
! ** Enthalpy of saturated air at surface temperature
        hs=(ho-bf*hi)/(1.-bf)
! ** Surface temperature
        ts=ftsat(hs,patmpa)
        return
        end subroutine tshobf
!
!=======================================================================
!
        real function qeffntu(cair,cwat,tai,twi,aext,rt,icountfl,eff)
!
! **********************************************************************
!  
!   PURPOSE:      Calculate duty for air-water heat exchangers.
!                 Sensible or total duty is calculated depending
!                 on whether dry bulb or wet bulb is used as the
!                 driving potential (for total duty, the air side
!                 capacity rate must be calculated from the 
!                 saturation specific heat and the air side 
!                 resistance reduced accordingly.
!  
! **********************************************************************
!   INPUTS
!   ======
!   cair    : air side capacity rate                              (kW/K)
!   cwat    : water side capacity rate                            (kW/K)
!   tai     : inlet air temperature (dry bulb or wet bulb)           (C)
!   twi     : inlet water temperature                                (C)
!   aext    : air side heat transfer area                           (m2)
!   rt      : total resistance, referred to unit area on air 
!           :                                             side (K.m2/kW)
!   icountfl: coil geometry (1=counterflow, 0=cross-flow - air side
!             unmixed and water side mixed, -1= parallel flow        (-)
!  
!   OUTPUT 
!   ======
!   qeffntu : coil duty                                             (kW) 
!   eff     : heat exchanger effectiveness                           (-) 
!  
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)           :: cair,cwat,tai,twi,aext,rt,eff, &
                                   crr,cmin,antu,effectiv
        integer                 :: iconfig,icountfl

! ** Set heat transfer configuration from geometry and capacity rates
        if (icountfl==1) then
! ** Counterflow
            iconfig=1
        elseif (icountfl==0) then
! ** Cross-flow - air side unmixed, water side mixed
            if (cair>cwat) then
! ** Maximum capacity rate side is unmixed
                iconfig=2
            else
! ** Minimum capacity rate side is unmixed
                iconfig=-2
            endif
        elseif (icountfl==-1) then
! ** Parallel flow
            iconfig=-1
        else
            write(*,*) 'function qeffntu: icountfl = ',icountfl, &
                       ' is out of range'
            stop
        endif
! ** Capacity rate ratio, minimum capacity and number of transfer units
        call capratio(cair,cwat,crr,cmin)
        antu = aext/(cmin*rt)
! ** Effectiveness and duty (+ve duty <-> heating of the air)
        eff  = effectiv(iconfig,crr,antu)
        qeffntu = eff*cmin*(twi-tai)
!
        return
        end function qeffntu
!
!======================================================================
!
        real function qairside(tai,ts,cair,aext,ra,rst1,bf)

! **********************************************************************
!  
!   PURPOSE:  Calculate the air-side resistance, by-pass factor and
!             heat transfer
!  
! **********************************************************************
!   INPUTS
!   ======
!   tai     : inlet air dry bulb temperature                         (C)
!   ts      : effective coil surface temperature                     (C)
!   cair    : air capacity rate                               (kJ/(s.K))  
!   aext    : area for air side heat transfer                       (m2)
!   ra      : air side thermal resistance                      (K.m2/kW)

!   OUTPUTS
!   =======
!   qairside: air side heat transfer                                (kW)
!   rst1    : air side thermal resistance referred to inlet    (K.m2/kW)
!   bf      : coil by-pass factor                                    (-)
!  
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)           :: tai,ts,cair,aext,ra,rst1,bf, &
                                   bypasfac
        real(kind=pp)           :: rsmall = 1.0e-3
! ** By-pass factor and resistance referred to the inlet temperature
        bf=bypasfac(cair,aext,ra)
        rst1=1./(cair*(1.-bf))
! ** Limit minimum value of rst1 to avoid division by zero
        rst1=max(rsmall,rst1)
! ** Air-side duty
        qairside=(ts-tai)/rst1
!
        return
        end function qairside
!
!======================================================================
!
        real function qwatside(aext,cair,cwat,rt,icountfl,rst1,twi,ts,rst2)

! **********************************************************************
!  
!   PURPOSE:  Calculate the water-side resistance from the total 
!             resistance and the air-side resistance and calculate
!             the water-side heat transfer
!  
! **********************************************************************
!   INPUTS
!   ======
!   aext    : area for air side heat transfer                       (m2)
!   cair    : air capacity rate                               (kJ/(s.K))  
!   cwat    : water capacity rate                             (kJ/(s.K))  
!   rt      : total thermal resistance                         (K.m2/kW)
!   icountfl: 0 = cross-flow, 1 = counterflow, -1 = parallel flow    (-)
!   rst1    : air side thermal resistance referred to inlet    (K.m2/kW)
!   twi     : inlet water temperature                                (C)
!   tsdyn   : effective coil surface temperature                     (C)
!  
!   OUTPUTS
!   =======
!   qwatside: water side heat transfer                              (kW)
!   rst2    : water side thermal resistance referred to inlet  (K.m2/kW)
!  
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)           :: aext,cair,cwat,rt,rst1,twi,ts,rst2
        real(kind=pp)           :: qeffntu,tai,eff,crr,cmin,antu,effectiv
        integer                 :: iconfig,icountfl
        real(kind=pp)           :: rsmall = 1.0e-3

! ** Set heat transfer configuration from geometry and capacity rates
        if (icountfl==1) then
! ** Counterflow
            iconfig=1
        elseif (icountfl==0) then
! ** Cross-flow - air side unmixed, water side mixed
            if (cair>cwat) then
! ** Maximum capacity rate side is unmixed
                iconfig=2
            else
! ** Minimum capacity rate side is unmixed
                iconfig=-2
            endif
        elseif (icountfl==-1) then
! ** Parallel flow
            iconfig=-1
        else
            write(*,*) 'function qeffntu: icountfl = ',icountfl,&
                       ' is out of range'
            stop
        endif
! ** Capacity rate ratio, minimum capacity and number of transfer units
        call capratio(cair,cwat,crr,cmin)
        antu= aext/(cmin*rt)
        eff = effectiv(iconfig,crr,antu)
        rst2= 1./(eff*cmin)-rst1
! ** Limit minimum value of rst2 to avoid division by zero
        rst2=max(rsmall,rst2)
! ** Water side heat transfer
        qwatside=(twi-ts)/rst2
        return
        end function qwatside
 
!======================================================================

        real function coilcap(vcpfin,as,thifin,vcptube,ttubelen,dotube, &
                        thitube,ditube)
!
! **********************************************************************
!  
!   PURPOSE:  Calculate the heat capacity of the fins, tubes and water
!             in a fin-tube heat exchanger
!  
! **********************************************************************
!   INPUTS
!   ======
!   vcpfin  : volumetric heat capacity of fin material       (kJ/(m3.K)) 
!   as      : net air side tube area                                (m2)
!   thifin  : fin thickness                                          (m)
!   vcptube : volumetric heat capacity of tube material      (kJ/(m3.K))  
!   ttubelen: total tube length                                      (m)
!   dotube  : tube outside diameter                                  (m)
!   thitube : tube thickness                                         (m)
!   ditube  : tube inside diameter                                   (m)
!  
!   OUTPUTS
!   =======
!   coilcap : coil heat capacity                                  (kJ/K)
!  
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)           :: vcpfin,as,thifin,vcptube,ttubelen, &
                                   dotube,thitube,ditube,hcfin,hctub,hcwat
        real(kind=pp)           :: pi = 3.14159, cpw=4.18, rhowat=1000.

! ** Fins
        hcfin  = vcpfin*as*thifin/2.0
! ** Tubes
        hctub  = vcptube*ttubelen*pi*dotube*thitube
! ** Water in tubes
        hcwat  = rhowat*cpw*ttubelen*pi*ditube*ditube/4.
! ** Total heat capacity
        coilcap = hcfin+hctub+hcwat
!
        return
        end function coilcap
!
! **********************************************************************
!
!   LIBRARY:        Routines for use in control strategies
!
!   DEVELOPER:      Philip Haves and Li Mei
!                   Loughborough University of Technology
!
!   REVISION DATE:  January 29, 1996
!
!   REFERENCE:      ASHRAE 825-RP Final Report
!
! **********************************************************************
!
! ======================================================================
!
!         Math Functions
!
! ======================================================================
!
        real function absvalue(a)

! Absolute value

        use precision1
        implicit none
        real(kind=pp)           :: a
 
        absvalue = abs(a)

        return
        end function absvalue
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        real function sqroot(a)

! Square root - value is negative if argument is negative

        use precision1
        implicit none
        real(kind=pp)           :: a

        if (a>=0) then
            sqroot = sqrt(a)
        else 
            sqroot = -sqrt(-a)
        endif

        return
        end function sqroot
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        real function sum(a,b)

! Sum first and second arguments
!
        use precision1
        implicit none
        real(kind=pp)           :: a, b

        sum = a+b
!
        return
        end function sum
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        real function diff(a,b)
!
! Subtract second argument from first argument
!
        use precision1
        implicit none
        real(kind=pp)           :: a, b
        diff = a-b

        return
        end function diff
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        real function product(a,b)
!
! Product of first and second arguments
!
        use precision1
        implicit none
        real(kind=pp)           :: a, b
        product = a*b

        return
        end function product
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        real function quotient(a,b)
!
! Quotient of first and second arguments
!
        use precision1
        implicit none
        real(kind=pp)           :: a, b
        if (b.ne.0) then
            quotient = a/b
        else 
            stop 'quotient: divide by zero'
        endif

        return
        end function quotient
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        real function average(a,n)
!
! Average of first n elements of a
!
        use precision1
        implicit none
        real(kind=pp)           :: suma
        integer                 :: i, n
        real(kind=pp), dimension(n)  :: a

        if (n<1) then
            stop 'function average: n<1'
        else
            suma = 0.0
            do 10 i=1,n
                suma = suma+a(i)
 10         continue
            average = suma/float(n)
        endif

        return
        end function average
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        real function sumarray(a,n)
!
! Sum of first n elements of a
!
        use precision1
        implicit none
        integer                 :: i, n
        real(kind=pp), dimension(n)  :: a

        if (n<1) then
            stop 'function sumarray: n<1'
        else
            sumarray = 0.0
            do i=1,n
                sumarray = sumarray+a(i)
            enddo
        endif

        return
        end function sumarray
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
        real function smallest(a,n)
!
! Minimum of first n elements of a
!
        use precision1
        implicit none
        real(kind=pp)                 :: large = 1.e+30
        integer                       :: i, n
        real(kind=pp), dimension(n)   :: a

        if (n<1) then
            stop 'function smallest: n<1'
        else
            smallest = large
            do i=1,n
                smallest = min(smallest,a(i))
            enddo
        endif

        return
        end function smallest
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
        real function biggest(a,n)
!
! Maximum of first n elements of a
!
        use precision1
        implicit none
        real(kind=pp)                 :: neglarge = -1.e+30
        integer                       :: i, n
        real(kind=pp), dimension(n)   :: a

        if (n<1) then
            stop 'function biggest: n<1'
        else
            biggest = neglarge
            do i=1,n
                biggest = max(biggest,a(i))
            enddo
        endif

        return
        end function biggest
!
! ======================================================================
!  
!       logic functions
!
! ======================================================================
!
        logical function logicnot(p)
!
! Logical compliment
!
        logical      :: p

        logicnot = .not. p

        return
        end function logicnot
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        logical function logicand(p,q)
!
! Logical AND of first and second arguments
!
        logical          :: p,q

        logicand = p .and. q

        return
        end function logicand
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        logical function logicor(p,q)
!
! Logical OR of first and second arguments
!
        logical          :: p,q

        logicor = p .or. q

        return
        end function logicor
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      
!
        logical function logicxor(p,q)
!
! Exclusive OR of first and second arguments (XOR)
!
        logical          :: p,q

        if (p .and. q) then
            logicxor = .false.
        else
            logicxor = p .or. q
        endif

        return
        end function logicxor
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        logical function retolog(a)
!
! Convert a real to a boolean
!
        if (a > 0.5) then
            retolog = .true.
        else
            retolog = .false.
        endif

        return
        end function retolog
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        real function logtore(p)
!
! Convert a boolean to a real
!
        logical           :: p

        if (p) then
           logtore = 1.0
        else
           logtore = 0.0
        endif

        return
        end function logtore
!
! ======================================================================
!  
!         psychrometric functions
!
! ======================================================================
!
         real function enthalpy(tdb,phi)
!
! Calculate the specific enthalpy of moist air (kJ/kg) at standard pressure
!
!   inputs
!   ======
!    1. tdb     : dry bulb temperature                               (C)
!    2. phi     : relative humidity (%)                              (-)

         use precision1
         implicit none
         real(kind=pp)                 :: tdb, phi, patm = 101325.
         real(kind=pp)                 :: w, fwphi, fhair

         w=fwphi(tdb,phi,patm)
         enthalpy=fhair(tdb,w)

         return
         end function enthalpy
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
         real function wetbulb(tdb,phi)
!
! Calculate the wet bulb temperature of moist air (C) at standard pressure
!
!   inputs
!   ======
!    1. tdb     : dry bulb temperature                               (C)
!    2. phi     : relative humidity (%)                              (-)

         use precision1
         implicit none
         real(kind=pp)                 :: tdb, phi, patm = 101325.
         real(kind=pp)                 :: w, ftwb, fwphi

         w=fwphi(tdb,phi,patm)
         wetbulb=ftwb(tdb,w,patm)

         return
         end function wetbulb
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
         real function dewpoint(tdb,phi)
!
! Calculate the dew point temperature of moist air (C) at standard pressure
!
!   inputs
!   ======
!    1. tdb     : dry bulb temperature                               (C)
!    2. phi     : relative humidity (%)                              (-)

         use precision1
         implicit none
         real(kind=pp)                 :: tdb, phi, patm = 101325.
         real(kind=pp)                 :: w, fwphi, ftdew

         w=fwphi(tdb,phi,patm)
         dewpoint=ftdew(w,patm)
         return
         end function dewpoint
!
! ======================================================================
!  
!       controller and sequencing functions
!
! ======================================================================
!
        real function clip(x,a,b)
!
! Limit the first argument to the range defined by the second and third
! arguments
!
        use precision1
        implicit none
        real(kind=pp)                 :: x, a, b

        if (a<b) then
            clip = min(b,max(a,x))
        else
            clip = min(a,max(b,x))
        endif

        return
        end function clip
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
         real function switch(p,a,b)
!
! Select the second argument if the first argument is true, else select
! the third argument
!
         use precision1
         implicit none
         real(kind=pp)                 :: a, b
         logical                       :: p

         if (p) then
            switch = a
         else
            switch = b
         endif

         return
         end function switch
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        real function span(x,a,b,c,d)
!
!  Rescale x from a-b to c-d
!
!   inputs
!   ======
!    1. x       : variable to be rescaled
!    2. a       : lower limit of input range
!    3. b       : upper limit of input range
!    4. c       : lower limit of output range
!    5. d       : upper limit of output range

        use precision1
        implicit none
        real(kind=pp)                 :: x, a, b, c, d, clip
        real(kind=pp)                 :: small = 1.e-15

        if (abs(a-b)>small) then
            span = (x-a)*(d-c)/(b-a) + c
        else
            stop &
       'function span: upper and lower limit of input range must differ'
        endif
        span=clip(span,c,d)

        return
        end function span
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
        subroutine deadband(y,w,db,wdb)
!
! Determine the effective setpoint in the presence of a deadband
!
!   inputs
!   ======
!    1. y       : measured variable
!    2. w       : setpoint
!    3. db      : half width of deadband (deadband extends from w-db to w+db)
!    4. wdb     : effective setpoint


        use precision1
        implicit none
        real(kind=pp)                 :: y, w, db, wdb

        if (y>=(w + db)) then
             wdb = w + db
        elseif (y<=(w - db)) then
             wdb = w - db
        else
             wdb = y
        endif

        return
        end subroutine deadband
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        logical function compare(a,b)
!
! Compare two real numbers - true if second argument is less than
! first argument
!
        use precision1
        implicit none
        real(kind=pp)                 :: a, b

        if (a>=b) then
           compare = .true.
        else 
           compare = .false.
        endif

        return
        end function compare
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        logical function comphys(a,b,db,comphysp)
!
! Compare two real numbers using a deadband - true if the first argument
! exceeds the second argument by more than the deadband or if the 
! magnitude of the difference is less than the deadband and the previous
! output was true
!
!   inputs
!   ======
!    1. a       : first number to be compared
!    2. b       : second number to be compared
!    3. db      : half width of deadband (deadband extends from w-db to w+db)
!    4. comphysp: previous output (stored in calling routine)

        use precision1
        implicit none
        real(kind=pp)                 :: a, b, db
        logical                       :: comphysp

        if (a>=(b+db)) then
           comphys = .true.
        elseif (a<=(b-db)) then
           comphys = .false.
        else
            comphys=comphysp
        endif

        return
        end function comphys
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        real function pidcont(y,w,pb,ti,td,intp,difp,pidp,ep,tsamp,umax,umin)
!
! Discrete-time Proportional plus Integral plus Derivative controller -
! integral action is implemented in feedback form to avoid "wind-up", and
! the effective derivative time is limited to ten sampling intervals.
!
! Ref: Clarke, D W, "PID algorithms and their computer implementation",
!     Department of Engineering Science, University of Oxford,
!     Report No. OUEL 1482/83, 1983.
!
!   inputs
!   ======
!    1. y       : process variable
!    2. w       : setpoint
!    3. pb      : proportional band, 0 -> no proportional action
!    4. ti      : integral time (sec), 0 -> no integral action
!    5. td      : derivative time (sec) 
!    6. intp    : integral term from previous sample instant (entry), 
!                 current value of integral term (return)
!    7. difp    : derivative term from previous sample instant (entry), 
!                 current value of derivative term (return)
!    8. pidp    : controller output from previous sample instant
!    9. ep      : error from prev sample instant (derivative action only)
!   10. tsamp   : sample time (sec) 
!   11. umax    : maximum output 
!   12. umin    : minimum output 
!   
!   output
!   ======
!      pidcont  = controller output


        use precision1
        implicit none
        real(kind=pp)                 :: y,w,pb,ti,td,intp,difp,pidp,ep,tsamp,umax,umin
        real(kind=pp)                 :: iterm,kp,pterm,dterm,clip,beta,d1,d2,e
        logical                       :: zeroi,zerod
        real(kind=pp)                 :: small = 1.e-15
!  Control error
        e = w - y
!  Determine proportional gain
        if (abs(pb)<small) then
            stop 'function pid: proportional band must not be zero'
        else
            kp = 1.0/pb
        endif
!  Test for no integral action
        zeroi = (ti<=0.0).or.((tsamp/ti)<1.0e-3)
!  Evaluate proportional and integral terms
        if (zeroi) then
!  No integral action
            pterm = kp*e
            iterm = 0.0
        else
!  Integral action
            if (ti<(tsamp/2.)) then
!  Integral gain too high for sampling interval
                stop 'function pidcont: ti < tsamp/2'
            else
!  Calculate proportional and integral terms using feedback form of
!  discretization
                pterm = kp*e*(1.0 + 0.5*tsamp/ti) 
                beta  = (2.0*ti - tsamp)/(2.0*ti + tsamp)
                iterm = beta*intp + (1.0 - beta)*pidp
            endif
        endif
!  Derivative term - check for no derivative action
        zerod = (td/tsamp)<1.0e-3
        if (zerod) then
!  No derivative action
            dterm = 0.0
        else
!  Derviative action - limit the effective derivative gain
            d1 = (0.2*td - tsamp)/(0.2*td + tsamp)
            if (d1<0.0) write(*,*) 'func pidcont: td < 5*tsamp'
            d2 = 2.0*td/(0.2*td + tsamp)
            dterm = d1*difp + d2*kp*(e - ep)
        endif
!  Controller output - limit to allowed range
        pidcont = pterm + iterm + dterm
        pidcont = clip(pidcont,umax,umin)
!  Update previous stored values
        intp = iterm
        difp = dterm
        ep   = e
        pidp = pidcont

        return
        end function pidcont

! **********************************************************************
!
!   LIBRARY:        Routines for use flow/pressure models
!  
!   DEVELOPER:      Philip Haves
!                   Loughborough University of Technology
!
!   REVISION DATE:  January 29, 1996 
!
!   REFERENCE:      ASHRAE 825-RP Final Report
!
!=======================================================================
!
        subroutine moistmix(t,g,wa,n,tmix,gmix,wmix)
!
! **********************************************************************
!  
!   PURPOSE:  Calculate the temperature and humidity ratio produced by
!             mixing up to six moist air streams
!  
! **********************************************************************
!   INPUTS
!   ======
!   t(1)    : dry bulb temperature of stream 1                       (C)
!   g(1)    : humidity ratio of stream 1                         (kg/kg)
!   wa(1)   : mass flow rate of stream 1                          (kg/s)
!   t(2)    : dry bulb temperature of stream 2                       (C)
!   g(2)    : humidity ratio of stream 2                         (kg/kg)
!   wa(2)   : mass flow rate of stream 2                          (kg/s)
!   ...
!  
!   OUTPUT
!   ======
!   tmix    : dry bulb temperature of mixed stream                   (C)
!   gmix    : humidity ratio of mixed stream                     (kg/kg)
!   wmix    : mass flow rate of mixed streaM(s)                   (kg/s)
!
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)                 :: tmix,gmix,wmix,sdwa,sdwv,sdh,sg,st
        real(kind=pp)                 :: cpa = 1., cpg = 1.805
        integer                       :: n, i
        real(kind=pp), dimension(n)   :: t, g, wa
        real(kind=pp), dimension(n)   :: dwa
        logical                       :: zeroi, zerod
        
! ** Perform moisture and enthalpy balances, checking for reverse flow
! ** sum positive inward dry air mass flow rates, moisture mass flow
! ** rates and enthalpy flow rates 
        sdwa = 0.0
        sdwv = 0.0
        sdh  = 0.0

        do i=1,n
            dwa(i) = max(0.0,wa(i))
            sdwa   = sdwa+dwa(i)
            sdwv   = sdwv+dwa(i)*g(i)
            sdh    = sdh+dwa(i)*t(i)*(cpa+cpg*g(i))
        enddo
        if (sdwa>0.0) then
! ** Some inward flow - Calculate mixed temperature and humidity ratio
            gmix   = sdwv/sdwa
            tmix   = sdh/(sdwa*(cpa+cpg*gmix))
        else
! ** No inlet flow, take mixed temperature and humidity ratio to be
! ** an average of inlet conditions so as to provide a
! ** reasonable guess for subsequent iterations, if required
            sg = 0.0
            st = 0.0
            do i=1,n
                sg = sg+g(i)
                st = st+t(i)
            enddo
            gmix = sg/float(n)
            tmix = st/float(n)
        endif
        wmix = sdwa
 
        return
        end subroutine moistmix
!=======================================================================
!
        real function dpqudlin(r,w)
!
! **********************************************************************
!  
!   PURPOSE:  Calculate pressure drop from flow rate and quadratic 
!             resistance - use linear flow relationship at low flow to 
!             avoid zero Jacobian. See ASHRAE 825-RP Final Report
!  
!             Uses (mass) flow rate as a very approximate proxy for 
!             Reynolds number (i.e. assumes a hydraulic diameter of 
!             0.1 - 1 m for air and 0.01 - 0.1 m for water)
!  
! **********************************************************************
!   INPUTS
!   ======
!   r       : (quadratic) flow resistance                   (0.001/kg.m)
!   w       : mass flow rate                                      (kg/s)
!  
!   OUTPUT
!   ======
!   dpqudlin: pressure drop                                        (kPa)
!  
! **********************************************************************
!   INTERNAL VARIABLES              
!   ==================
!   wcrit   : critical flow, lower limit of purely quadratic flow (kg/s)
!  
! **********************************************************************

        use precision1
        implicit none
        real(kind=pp)                 :: r, w
        real(kind=pp)                 :: dpqud, dplin
        real(kind=pp)                 :: wcrit = 6.e-2

! ** The relationship used depends on the magnitude of the flow
        if (abs(w)>=wcrit) then
! ** Above the critical bound - fully quadratic flow
            dpqudlin = dpqud(r,w)
        else
! ** Below the critical bound - use linear relationship
            dpqudlin = dplin(r,wcrit,w)
        endif

        return
        end function dpqudlin
!
!=======================================================================
!
        real function dplin(r,wcrit,w)
!
! **********************************************************************
!  
!   PURPOSE:  Calculate pressure drop from flow rate and quadratic 
!             resistance at low flow - use linear flow relationship to 
!             avoid zero Jacobian. See ASHRAE 825-RP Final Report
!  
!             Uses (mass) flow rate as a very approximate proxy for 
!             Reynolds number (i.e. assumes a hydraulic diameter of 
!             0.1 - 1 m for air and 0.01 - 0.1 m for water)
!  
! **********************************************************************
!   INPUTS
!   ======
!   r       : (quadratic) flow resistance                   (0.001/kg.m)
!   wcrit   : critical flow, lower limit of purely quadratic flow (kg/s)
!   w       : mass flow rate                                      (kg/s)
!  
!   OUTPUT
!   ======
!   dplin   : pressure drop                                        (kPa)
!
! **********************************************************************
!
! ** Linear dependence of pressure drop on flow rate
        use precision1
        implicit none
        real(kind=pp)                 :: r, w, wcrit

        dplin = r*wcrit*w

        return
        end function dplin
!
!=======================================================================
!
        real function dpqud(r,w)
!
! **********************************************************************
!  
!   PURPOSE:  Calculate pressure drop from flow rate and quadratic 
!             resistance
!  
! **********************************************************************
!   INPUTS
!   ======
!   r       : (quadratic) flow resistance                   (0.001/kg.m)
!   w       : mass flow rate                                      (kg/s)
!  
!   OUTPUT
!   ======
!   dpqud: pressure drop                                        (kPa)
!  
! **********************************************************************
!
! ** Quadratic dependence of pressure drop on flow
! ** rate, treats negative flows and resistances
        use precision1
        implicit none
        real(kind=pp)                 :: r, w

        dpqud = r*w*abs(w)

        return
        end function dpqud
!
!=======================================================================
!
        real function wqudlin(r,dp)
!
! **********************************************************************
!  
!   PURPOSE:  Calculate flow rate from pressure drop and turbulent 
!             resistance - use laminar flow relationship at low flow to 
!             avoid zero Jacobian. See ASHRAE 825-RP Final Report
!  
!             Uses (mass) flow rate as a very approximate proxy for 
!             Reynolds number (i.e. assumes a hydraulic diameter of 
!             0.1 - 1 m for air and 0.01 - 0.1 m for water)
!  
! **********************************************************************
!   INPUTS
!   ======
!   r       : (quadratic) flow resistance                   (0.001/kg.m)
!   dp      : pressure drop                                        (kPa)
!  
!   OUTPUT
!   ======
!   wqudlin : mass flow rate                                      (kg/s)
!  
! **********************************************************************
!   INTERNAL VARIABLES              
!   ==================
!   wcrit   : critical flow, lower limit of purely quadratic flow (kg/s)
!
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)                 :: r, dp, dpcrit
        real(kind=pp)                 :: wqud, wlin, wcrit = 6.e-2

! ** The relationship used depends on the magnitude of the pressure drop
! ** critical pressure drop corresponding to critical flow rate
        dpcrit = abs(r*wcrit*wcrit)
        if (abs(dp)>dpcrit) then
! ** Above the critical bound - fully quadratic flow
            wqudlin=wqud(r,dp)
        else
! ** Below the critical bound - use linear relationship
            wqudlin=wlin(r,wcrit,dp)
        endif

        return
        end function wqudlin
!
!=======================================================================
!
        real function wlin(r,wcrit,dp)
!
! **********************************************************************
!  
!   PURPOSE:  Calculate flow rate from pressure drop and turbulent 
!             resistance at low flow - use laminar flow relationship to 
!             avoid zero Jacobian. See ASHRAE 825-RP Final Report
!  
!             Uses (mass) flow rate as a very approximate proxy for 
!             Reynolds number (i.e. assumes a hydraulic diameter of 
!             0.1 - 1 m for air and 0.01 - 0.1 m for water)
!  
! **********************************************************************
!   INPUTS
!   ======
!   r       : (quadratic) flow resistance                   (0.001/kg.m)
!   wcrit   : critical flow, lower limit of purely quadratic flow (kg/s)
!   dp      : pressure drop                                        (kPa)
!  
!   OUTPUT
!   ======
!   wlin    : mass flow rate                                      (kg/s)
!
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)                 :: r, wcrit, dp

! * Trap zero resistance
        if (r==0.0) then
! ** Zero resistance - potentially infinite flow rate
            stop 'function wlin: zero resistance'
        endif
! * Trap zero critical flow
        if (wcrit==0.0) then
! ** Zero critical flow - potentially infinite flow rate
            stop 'function wlin: zero critical flow'
        endif
! ** Linear dependence of pressure drop on flow rate
        wlin = dp/(r*wcrit)

        return
        end function wlin
!
!=======================================================================
!
        real function wqud(r,dp)
!
! **********************************************************************
!  
!   PURPOSE:  Calculate flow rate from pressure drop and turbulent 
!             resistance 
!  
! **********************************************************************
!   INPUTS
!   ======
!   r       : (quadratic) flow resistance                   (0.001/kg.m)
!   dp      : pressure drop                                        (kPa)
!  
!   OUTPUT
!   ======
!   wqud    : mass flow rate                                      (kg/s)
!
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)                 :: r, dp

! * Trap zero resistance
        if (r==0.0) then
! ** Zero resistance - potentially infinite flow rate
            stop 'function wqud: zero resistance'
        endif
! ** Quadratic dependence of pressure drop on flow
! ** rate, treats negative flows and resistances
        wqud = sign(1.0,(dp/r))*sqrt(abs(dp/r))

        return
        end function wqud
!
!======================================================================
!
        real function rlinport(x,kv,sv,xlin,cl)
!
! *********************************************************************
!  
!   PURPOSE:  Calculate hydraulic resistance of a linear port in a
!             control valve. Valve characteristic has two segments:    
!  
!             1. Linear from (0,cl) to (xlin,1/sv)                           
!             2. Linear from (xlin,1/sv) to (1,1)
!                                                                     
! *********************************************************************
!   INPUT
!   =====
!   x       : valve stem position (assumed limited to range 0-1)    (-)
!   kv      : valve capacity                   (flow in m3/hr at 1 bar)  
!   sv      : rangability - ratio of highest to lowest 
!             controllable flow                                     (-)
!   xlin    : stem position at which close-off starts               (-)
!   cl      : fractional leakage flow when fully closed             (-)
!   
!   OUTPUT 
!   ======
!   rlinport: hydraulic resistance of valve              (0.001/(kg.m))
!
! *********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)                 :: x, kv, sv, xlin, &
                                         cl, recipsv, fracflow, slope
        real(kind=pp)                 :: small = 1.e-8

! ** Determine the segment of the characteristic and then determine
! ** the resistance
! ** Check that the characteristic parameters are consistent
        recipsv = 1.0/sv
        if (cl > recipsv) stop    &
           'rlinport: leakage > 1/rangability'
        if (cl < small) stop &
           'rlinport: leakage is too small'
        if (kv < small) stop &
           'rlinport: valve capacity is too small'
        if (xlin < small) stop &
           'rlinport: close-off/linear break-point is too small'
        if (xlin > 1.0) stop &
           'rlinport: close-off/linear break-point is too large'
! ** Two segments - determine fractional flow
        if (x <= xlin) then
! ** Shut-off
            fracflow=((recipsv-cl)/xlin)*x + cl
! ** Linear
        else
            slope = (1.0-recipsv) / (1.0-xlin)
            fracflow = slope*(x-xlin) + recipsv
        endif
! ** Resistance from kv and fractional flow
        rlinport = 1296./(kv*kv*fracflow*fracflow)
        return
        end function rlinport
!
!======================================================================
!
        real function reqpport(x,kv,xeqp,eqpchar,xlin,sv,cl)
!
! *********************************************************************
!  
!   PURPOSE:  Calculate hydraulic resistance of an equal percentage 
!             port in a control valve. Valve characteristic has three
!             segments:    
!  
!             1. Linear from (0,CL) to (XLIN,1/SV)                           
!             2. Linear from (XLIN,1/SV) to (XEQP,YEQP)
!             3. Exponential from (XEQP,YEQP) to (1,1)
!  
!             where: YEQP = EQPCHAR**(XEQP-1)
!                                                                     
! *********************************************************************
!   INPUT
!   =====
!   y       : valve stem position (assumed limited to range 0-1)    (-)
!   kv      : valve capacity                   (flow in m3/hr at 1 bar)  
!   xeqp    : stem position of break-point between linear and equal 
!             percentage characteristics                            (-)
!   eqpchar : curvature parameter for equal percentage 
!             characteristic                                        (-)
!   xlin    : stem position at which close-off starts               (-)
!   sv      : rangability - ratio of highest to lowest 
!             controllable flow                                     (-)
!   cl      : fractional leakage flow when fully closed             (-)
!   
!   OUTPUT 
!   ======
!   reqpport: hydraulic resistance of valve              (0.001/(kg.m))
!
! *********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)                 :: x,kv,xeqp,eqpchar,xlin,sv,cl
        real(kind=pp)                 :: small = 1.e-8, yeqp, recipsv,slope,fracflow

! ** Determine the segment of the characteristic and then determine
! ** the resistance
! ** Check that the characteristic parameters are consistent
        recipsv = 1.0/sv
        if (cl > recipsv) stop &
           'reqpport: leakage > 1/rangability'
        if (cl < small) stop &
           'reqpport: leakage is too small'
        if (kv < small) stop &
           'reqpport: valve capacity is too small'
        if (xlin < small) stop &
           'reqpport: close-off/linear break-point is too small'
        if (xlin > 1.0) stop &
           'reqpport: close-off/linear break-point is too large'
        if (xlin > xeqp) stop &
           'reqpport: break-points in wrong order'
           yeqp = eqpchar**(xeqp-1.0)
        if (yeqp < recipsv) stop &
           'reqpport: linear segment has negative slope'
! ** Three segments:
        if (x <= xlin) then
! ** Shut-off
            fracflow=((recipsv-cl)/xlin)*x + cl
        elseif (x <= xeqp) then
! ** Linear
            slope = (yeqp - recipsv) / (xeqp-xlin)
            fracflow = slope*(x-xlin) + recipsv
        else
! ** Exponential characteristic
            fracflow = eqpchar**(x-1.0)
        endif
! ** Resistance from kv and fractional flow
        reqpport = 1296./(kv*kv*fracflow*fracflow)
        return
        end function reqpport
        
!======================================================================
!
        real function rdamper(c,ropen,fleak,area,a,b,ip)
!
! *********************************************************************
!  
!   PURPOSE:  Calculate resistance of damper using Legg's correlation
!             for the log of the loss factor, ln(K)=A+B*THETA, in the 
!             range 15< THETA < 55 deg (65 for parallel). Use quadratic     
!             interpolation for 0 < THETA < 15 and 55(65) < THETA < 90.
!  
! *********************************************************************
!   INPUT
!   =====
!   c       : damper position (0=closed, 1=open)                    (-)
!   ropen   : open resistance of damper                  (0.001/(kg.m))
!   fleak   : leakage, expressed as a fraction of full flow         (-)  
!   area    : face area of damper                                   (-)  
!   a       : A coefficient in ln(K)=A+B*THETA                      (-)  
!   b       : B coefficient in ln(K)=A+B*THETA                      (-)  
!   ip      : damper type (0=opposed/single, 1=parallel)            (-)
!   
!   OUTPUT 
!   ======
!   rdamper : resistance of damper                        (0.001/(kg.m)
!     
! *********************************************************************
!   INTERNAL VARIABLES
!   ==================
!   theta   : damper position (90=closed, 0=open)                 (deg)
!   thetal  : lower limit of validity of correlation              (deg)
!   thetau  : upper limit of validity of correlation              (deg)
!   rhoa    : density of air                                    (kg/m3)
!   rtok    : conversion from resistance to loss factor     (1000.kg.m)  
!   thetae  : damper position at end point of correlation         (deg)  
!   lke     : nat log of loss factor at end point of correlation    (-)  
!   dlkedthe: slope of Legg correlation at end point                (-)  
!   thetap  : position of damper at extreme position (0 or 90 deg)(deg)  
!   lkp     : nat log of loss factor at THETAP                      (-)  
!   aa      : coefficient of x**2 in interpolation function         (-)
!   bb      : coefficient of x**1 in interpolation function         (-)
!   cc      : coefficient of x**0 in interpolation function         (-)
!   lk      : nat log of loss factor at THETA                       (-)  
!
! *********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)                 :: c,ropen,fleak,area,a,b
        real(kind=pp)                 :: lk,lke,lkp,rtok,theta,thetae,dlkedthe, &
                                         aa,bb,cc,yquad,thetap
        real(kind=pp), dimension(0:1) :: thetau = (/55.0,65.0/)
        real(kind=pp)                 :: rhoa = 1.2, small = 1.e-8, thetal = 15.0
        integer                       :: ip

! ** Check for effectively zero area
        if (area<small) then
            write(*,*) 'function rdamper: area < ', small
            stop
        endif
! ** Conversion factor from resistance to loss coefficient
        rtok=1000.*2.*rhoa*area*area
!        rtok=2.*rhoa*area*area
! ** Convert control signal to degrees
! ** (0 deg <-> open <-> c=1, 90 deg <-> closed <-> c=0)
        theta=90.*(1.-c)
! ** Limit damper position
        if (theta<0.0) then
            theta=0.0
        elseif (theta>90.) then
            theta=90.
        endif
! ** Calculate damper loss factor and resistance
        if (theta==0.) then
! ** Damper fully open
            rdamper=ropen
        elseif (theta<thetal) then
! ** Nearly fully open - quadratic interpolation between
! ** correlation and open resistance
! ** Position, nat log of loss factor and slope at lower limit of
! ** correlation
            thetae=thetal
            lke=a+b*thetal
            dlkedthe=b
! ** Position and nat log of loss factor at 0 deg
            thetap=0.
            lkp=log(ropen*rtok)
! ** Determine coefficients of interpolation function for 0<theta<thetal
            call interpar(thetae,lke,dlkedthe,thetap,lkp,aa,bb,cc)
! ** Evaluate interpolation function at theta
            lk=yquad(theta,aa,bb,cc)
! ** Resistance of damper from nat log of loss coefficient
            rdamper=exp(lk)/rtok
        elseif (theta.le.thetau(ip)) then
! ** Mid range - legg correlation valid
! ** nat log of loss coefficient
            lk=(a+b*theta)
! ** Resistance of damper
            rdamper=exp(lk)/rtok
        elseif (theta<90.) then
! ** Nearly closed - quadratic interpolation between correlation and
! ** leakage resistance
! ** Position, nat log of loss factor and slope at upper limit of
! ** correlation
            thetae=thetau(ip)
            lke=a+b*thetau(ip)
            dlkedthe=b
! ** Position and nat log of loss factor at 90 deg
            thetap=90.
            lkp=log(rtok*ropen/(fleak*fleak))
! ** dDtermine coefficients of interpolation function for thetau(ip)<theta<90
            call interpar(thetae,lke,dlkedthe,thetap,lkp,aa,bb,cc)
! ** Evaluate interpolation function at theta
            lk=yquad(theta,aa,bb,cc)
! ** Resistance of damper from nat log of loss coefficient
            rdamper=exp(lk)/rtok
        else
! ** Damper fully closed - calculate leakage resistance
            rdamper=ropen/(fleak*fleak)
        endif
        return
        end function rdamper
!======================================================================

        subroutine interpar(xe,ye,dyedxe,xp,yp,a,b,c)
!
! *********************************************************************
!  
!   PURPOSE:  Determine coefficients of line <-> point interpolation 
!             function
!  
!             The interpolation function connects the end point of the 
!             line, E, to the point P. The gradient is continuous at E.
!             The interpolation function is Y = A*X**2 + B*X + C, so
!             that the gradient varies linearly with X. 
!  
! *********************************************************************
!   INPUT
!   =====
!   xe      : x co-ordinate of end of line                          (-)
!   ye      : y co-ordinate of end of line                          (-)
!   dyedxe  : gradient of line evaluated at end point               (-)  
!   xp      : x co-ordinate of point                                (-)  
!   yp      : y co-ordinate of point                                (-)  
!   
!   OUTPUT 
!   ======
!   a       : a coefficient of interpolation function               (-)  
!   b       : b coefficient of interpolation function               (-)  
!   c       : c coefficient of interpolation function               (-)  
!
! *********************************************************************
        use precision1
        implicit none
        real(kind=pp)                 :: xe,ye,dyedxe,xp,yp,a,b,c
        real(kind=pp)                 :: dx
!
! ** Check the point and the end of the line are separated in the
! ** x direction
        if (xp==xe) then
            stop 'subroutine interpar: xp=xe'
        endif
! ** Calculate the coefficients of the interpolation function
        dx = xp-xe
        a = ( (yp-ye) - dyedxe*dx ) / (dx*dx)
        b = dyedxe - 2.0*a*xe
        c = ye - a*xe*xe - b*xe
        return
        end subroutine interpar
!======================================================================
!
        real function yquad(x,a,b,c)
!
! *********************************************************************
!
!   PURPOSE:  Evaluate quadratic function  A*X**2 + B*X + C
!
! *********************************************************************
!   INPUT
!   =====
!   x       : independent variable                                  (-)
!   a       : coefficient of X**2                                   (-)  
!   b       : coefficient of X**1                                   (-)  
!   c       : coefficient of X**0                                   (-)  
!   
!   OUTPUT 
!   ======
!   yquad   : dependent variable                                    (-)
!
! *********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)                 :: x,a,b,c

        yquad   = a*x*x + b*x + c
        return
        end function yquad
!
!======================================================================
!
        subroutine flowmerg(wi1,wi2,po,ri1,ri2,ro,pi1,pi2,wo)
!
! *********************************************************************
!  
!   PURPOSE:     Calculate the outlet flow rate and the inlet pressures 
!                for a flow merge, given the outlet pressure and the 
!                inlet flow rates
!             
!   LIMITATIONS: 
!  
!   REFERENCE:   Haves, P., Component-Based Modelling of VAV Systems, 
!                Proc. System Simulation in Buildings '94, Liege, 
!                Belgium, December 1994
!  
! **********************************************************************
!   INPUTS                                       
!   ======                                       
!   wi1     : Inlet flow rate 1                                   (kg/s)
!   wi2     : Inlet flow rate 2                                   (kg/s)
!   po      : Pressure at outlet                                   (kPa)
!                                                
!   parameters
!   ==========
!   ri1     : Inlet resistance 1                             (0.001/k.m)
!   ri2     : Inlet resistance 2                             (0.001/k.m)
!   ro      : Resistance of outlet                           (0.001/k.m)
!  
!   outputs                                      
!   =======                                     
!   wo      : Mass flow rate at outlet                            (kg/s)
!   pi1     : Inlet pressure 1                                     (kPa)
!   pi2     : Inlet pressure 2                                     (kPa)
!  
! **********************************************************************
        use precision1
        implicit none
        real(kind=pp)                 :: wi1,wi2,po,ri1,ri2,ro,pi1,pi2,wo, &
                                         dpqudlin, dpo, dri1, dri2,dpi1,dpi2

! *** Calculate outlet flow rate and inlet pressures
! *** Outlet flow rate is sum of inlet flow rates
        wo=wi1+wi2
! *** Pressure drops in each branch - quadratic dependence on flow at
! *** higher flow rates, linear at low flow rates
        dpi1=dpqudlin(ri1,wi1)
        dpi2=dpqudlin(ri2,wi2)
        dpo=dpqudlin(ro,wo)
! *** Inlet pressures - outlet pressure plus pressure drops in branches
        pi1=po+dpi1+dpo
        pi2=po+dpi2+dpo

        return
        end subroutine flowmerg
!
!======================================================================
!
        subroutine flowsplt(wi,po1,po2,ri,ro1,ro2,dummy,wcrit,tol, &
                            pi,wo1,wo2,ifail)
!
! *********************************************************************
!  
!   PURPOSE:     Calculate the outlet flow rates and the inlet pressure 
!                for a flow split, given the outlet pressures and the 
!                inlet flow rate
!             
!   LIMITATIONS: Resistances must not be negative. At least one outlet
!                resistance must be non-zero.
!  
!   REFERENCE:   Haves, P., Component-Based Modelling of VAV Systems, 
!                Proc. System Simulation in Buildings '94, Liege, 
!                Belgium, December 1994
!  
! **********************************************************************
!   INPUTS
!   ======
!   wi      : Inlet flow rate                                     (kg/s)
!   po1     : Pressure at outlet 1                                 (kPa)
!   po2     : Pressure at outlet 2                                 (kPa)
!   ri      : Inlet resistance                               (0.001/k.m)
!   ro1     : Resistance of outlet 1                         (0.001/k.m)
!   ro2     : Resistance of outlet 2                         (0.001/k.m)
!   wcrit   : Critical flow, lower limit of purely quadratic solution
!             also used to linearise resistances                  (kg/s)
!   tol     : Relative tolerance for flow rates (TOL**2 for pressures)
!                                                                    (-)
!  
!   OUTPUTS
!   =======
!   pi      : Inlet pressure                                       (kPa)
!   wo1     : Mass flow rate at outlet 1                          (kg/s)
!   wo2     : Mass flow rate at outlet 2                          (kg/s)
!   ifail   : 0=successful completion, 1=negative resistance(s),
!             2=both outlet resistances zero                         (-)
!  
! **********************************************************************
!   INTERNAL VARIABLES
!   ==================
!   dp21    : difference between the outlet pressures
!   wo21    : proxy for the difference between the outlet flows
!   wmax    : estimate of the largest of the three flows
!   wcrit   : flow threshold, lower limit of pure quadratic solution
!   sri     : "signed" inlet resistance
!   sr01    : "signed" resistance for outlet 1
!   sr02    : "signed" resistance for outlet 2
!   rtot    : sum of the three branch resistances
!   smallres: threshold for resistances to be considered negligible
!   abswo1  : absolute value of outlet flow 1
!   abswo2  : absolute value of outlet flow 2
!   b2m4acd4: (b**2-4ac)/4 for the quadratic equation for the general case
!   sqr     : the square root of (b**2-4ac)/4
!   wo1p    : the value of WO1 corresponding to the positive root
!   wo1n    : the value of WO1 corresponding to the negative root
!   wo2p    : the value of WO2 corresponding to the positive root
!   wo2n    : the value of WO2 corresponding to the negative root
!   posiposs: true if the positive root is a possible solution
!   negaposs: true if the negative root is a possible solution
!   piq     : the inlet pressure corresponding to the quadratic solution
!   ril     : linearised inlet resistance
!   ro1l    : linearised resistance for outlet 1
!   ro2l    : linearised resistance for outlet 2
!
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)                :: wi,po1,po2,ri,ro1,ro2,dummy,wcrit,tol, &
                                        pi,wo1,wo2,abswo2,wo2q,wo1q, &
                                        b2m4acd4,wo1p,wo1n,wo2p,wo2n,sqr, &
                                        ril,ro1l,ro2l,dp21,wo21,wmax, &
                                        sri,sro1,sro2,rtot,smallres,abswo1
        integer                      :: ifail

        logical                      :: posiposs, negaposs

! ** Trap negative branch resistance
        if (ri<0.0 .or. ro1<0.0 .or. ro2<0.0) then
            ifail=1
            return
        endif
! ** Trap zero resistance for both outlet branches
        if (ro1==0.0 .and. ro2==0.0) then
            ifail=2
            return
        endif
! ** Determine whether to use quadratic or linear pressure vs. flow
! ** relationships for the outlet branches based on the magnitudes of the
! ** flows in the outlet branches. Use linear relationships if, and only
! ** if, both flows are small. Both flows are deemed to be small if their
! ** algebraic sum (i.e. the inlet flow) is small and a quantity related
! ** to their difference is also small.
! ** Calculate wo21, the flow that would occur between the outlets if the
! ** inlet flow were zero, which is taken as a proxy for the difference
! ** between the actual outlet flows
        dp21=po2-po1
        wo21=sqrt(abs(dp21)/(ro1+ro2))
! ** Use the greater of the inlet flow and wo21 as an estimate of the
! ** largest flow, to be used to test whether to use the quadratic or the
! ** linear solution
        wmax=max(abs(wi),wo21)
! ** Two thresholds are defined. if the estimate of the largest flow is
! ** above the upper threshold, use the quadratic solution. if it is
! ** below the lower threshold, use the linear solution. if it is between
! ** the thresholds, use a linear combination of the quadratic and linear
! ** solutions.
        if (wmax>wcrit) then
! ** Above the threshold - Calculate the quadratic solution
! ** Calculate "signed" resistances (the products of the flow directions
! ** and the resistances)
            sri=sign(ri,wi)
            sro1=sign(ro1,(ro2*wi*abs(wi)+dp21))
            sro2=sign(ro2,(ro1*wi*abs(wi)-dp21))
! ** Trap cases where one resistance is very small
            rtot=ri+ro1+ro2
            smallres=rtot*tol*tol
            if (ro1<smallres) then
! ** Ro1 is very small, so dp21=po2-po1 is pressure drop across ro2
                abswo2=sqrt(abs(dp21/ro2))
! ** Sign of sro2 indicates direction of flow through ro2
                wo2q=sign(abswo2,sro2)
                wo1q=wi-wo2q
            elseif (ro2<smallres) then
! ** Ro2 is very small, so dp21=po2-po1 is pressure drop across ro1
                abswo1=sqrt(abs(dp21/ro1))
! ** Sign of sro1 indicates direction of flow through ro1
                wo1q=sign(abswo1,sro1)
                wo2q=wi-wo1q
! ** Trap cases where one outlet flow is zero
            elseif (sro1==0.0) then
! ** Outlet flow in branch 1 is zero (ro1=0 already trapped)
                wo1q=0.0
                wo2q=wi
            elseif (sro2==0.0) then
! ** Outlet flow in branch 2 is zero (ro2=0 already trapped)
                wo2q=0.0
                wo1q=wi
            elseif (abs(sro1-sro2)<smallres) then
! ** Special case - outlet resistances equal and flows have same sign
                if (wi==0.0) then
! ** Inlet flow zero and outlet flows have same sign: outlet flows zero
                    wo1q=0.0
                    wo2q=0.0
                else
! ** General quadratic equation reduces to a linear equation (see printed
! ** documentation). n.b. sro2=0 already trapped
                    wo1q=wi/2.+dp21/(2.*sro2*wi)
                    wo2q=wi-wo1q
                endif
            else
! ** General case
! ** Determine roots of quadratic (see printed documentation)
! ** Calculate (b**2 - 4ac)/4
                b2m4acd4=(sro1*sro2*wi*wi-(sro2-sro1)*dp21)
                if (b2m4acd4.lt.0.0) then
!                    print *,'subroutine flowsplt: complex roots'
                    stop
                endif
! ** Calculate outlet flows corresponding to positive and negative roots
                sqr=sqrt(b2m4acd4)
                wo1p=(sro2*wi+sqr)/(sro2-sro1)
                wo1n=(sro2*wi-sqr)/(sro2-sro1)
                wo2p=wi-wo1p
                wo2n=wi-wo1n
! ** Check for solution(s) that are consistent with known flow directions
                if (wo1p*sro1>=0.0 .and. wo2p*sro2>=0.0) then
                    posiposs=.true.
                else
                    posiposs=.false.
                endif
                if (wo1n*sro1>=0.0 .and. wo2n*sro2>=0.0) then
                    negaposs=.true.
                else
                    negaposs=.false.
                endif
!  Select consistent root - stop if no roots or two roots consistent
                if (posiposs .and. .not.negaposs) then
                    wo1=wo1p
                    wo2=wo2p
                elseif (negaposs .and. .not.posiposs) then
                    wo1=wo1n
                    wo2=wo2n
                else
                    stop &
                    'subroutine flowsplt: ambiguous or inconsistent roots'
                endif
            endif
            pi=po1+sri*wi*wi+sro1*wo1*wo1
        elseif (wmax<=wcrit) then
! ** Small flows, calculate linear solution
            ril=ri*wcrit
            ro1l=ro1*wcrit
            ro2l=ro2*wcrit
            wo1=(ro2l*wi+dp21)/(ro1l+ro2l)
            wo2=wi-wo1
            pi=po1+ril*wi+ro1l*wo1
        endif
! ** Successful completion
        ifail=0
        return
        end subroutine flowsplt
!
! *********************************************************************
!
        real function spedlim(fspeddem,fspedmin)
!
        use precision1
        implicit none
        real(kind=pp)                 :: fspeddem,fspedmin
        real(kind=pp)                 :: clip,a=1.0,b=-2.0

! ** Apply upper and lower limits to actuator speed
! ** Limit to maximum speed
        fspeddem = clip(fspeddem,a,b)
!        fspeddem = clip(fspeddem,1.0,-1.0)
! ** Limit to minimum speed, or if demanded speed is less than half the
! ** minimum speed, set speed to zero
        if (abs(fspeddem)>=fspedmin) then
            spedlim = fspeddem
        elseif (abs(fspeddem)>=(fspedmin/2.0)) then
            spedlim = sign(fspedmin,fspeddem)
        else
            spedlim = 0.0
        endif
        return
        end function spedlim
!
! ********************************************************************** 
!
        subroutine hystrsis(posin,posoutp,hys,posout)
!
! ********************************************************************** 
!  
!   PURPOSE:  Models hysteresis of actuators etc.
!  
! **********************************************************************
!   INPUTS
!   ======
!   posin   : position of input (e.g. actuator) (0<=POSIN<=1)
!   posoutp : position of output at previous step time (0<=POSOUTP<=1-HYS)
!   hys     : hysteresis (fraction of full range)     
!  
!   OUTPUT
!   ======
!   posout  : position of output (0<=POSOUT<=1-HYS, i.e. not rescaled)
!  
! **********************************************************************
!
        use precision1
        implicit none
        real(kind=pp)                 :: posin,posoutp,hys,posout,clip
        real(kind=pp)                 :: a, b
! ** Limit input position to 0-1
        a=0.0
        b=1.0
        posin=clip(posin,a,b)           ! <<<< 5/5/05

!        posin=clip(posin,0.0,1.0)
! ** Determine new output position
! ** Rescale previous output position
        if ( posin>(posoutp+hys) ) then
! ** New input is above slack region
            posout=posin-hys
        elseif ( posin<posoutp )  then
! ** New input is below slack region
            posout=posin
        else
! ** New input is within slack region
            posout=posoutp
        endif
        return
        end subroutine hystrsis

! **********************************************************************
!   Psychrometric functions as a function of the barometric pressure - *
!   based of the ASHRAE 1981 Fundamentals Handbook -  Chapter 5             *
!                                                                      *
!   units : kJ, kg, Pa, degree C                                       *
!                                                                      *
!       function name        returned value                            *
!   ___________________________________________________________________*
!    1. fpws(tdb)            pws  = saturated vapour pressure(Pa)      *
!    2. ftdew(w,patm)        tdew = dew point temp. or satur. temp.(C) *
!    3. fpww(w,patm)         pw   = vapor pressure(Pa)                 *
!    4. fwpw(pw,patm)        w    = humidity ratio(kg moist./kg dry air*
!    5. fwphi(tdb,phi,patm)  w    = humidity ratio(kg moist./kg dry air*
!    6. fwtwb(tdb,twb,patm)  w    = humidity ratio(kg moist./kg dry air*
!    7. fwha(tdb,ha)         w    = humidity ratio(kg moist./kg dry air*
!    8. ftdb(w,ha)           tdb  = air dry bulb temperature(C)        *
!    9. fphi(tdb,w,patm)     phi  = relative humidity(%)               *
!   10. fhair(tdb,w)         hair = enthalpy of moist air(kJ/kg dry air*
!   11. fhsat(tsat,patm)     hsat = saturated enthalpy(kJ/kg dry air)  *
!   12. ftsat(hs,patm)       tsat = saturated temperature(C)           *
!   13. ftwb(tdb,w,patm)     twb  = wet bulb temperature(C)            *
!   14. ftair(ha,w)          tair = temperature of moist air(C)        *
!                                                                      *
!   ___________________________________________________________________*
!   The specific heats are not varying with the temperature            *
!   cpa = 1.006 kJ/kg*K, specific heat of dry air                      *
!   cpg = 1.805 kJ/kg*K, specific heat of water vapor                  *
!   hfg = 2501. kJ/kg,   latent heat of vaporisation                   *
!   cpw = 4.187 kJ/kg*K, specific heat of water                        *
! **********************************************************************
!
! *********************
!   list of variables *
! *********************
!
! ha   = h, air enthalpy (kJ/kg dry air)
! hsat = hs, air saturation enthapy (kJ/kg dry air)
! patm = atmospheric pressure (Pa)
! phi  = relative humidity (%)
! ps   = pws, water vapor saturation pressure (Pa)
! pw   = partial pressure of water vapor (Pa)
! tdb  = air dry bulb temperature (C)
! tdew = dew point temperature (C)
! tsat = ts, saturation temperature (C)
! twb  = air wet bulb temperature (C)
! w    = humidity ratio (kg moisture/kg dry air)
! ws   = humidity ratio at saturation conditions (kg moisture/kg dry air
!
!
! ****************************************************************
!   1. Saturation pressure over liquid water for the temperature *
!   range of 0 deg C to 200 deg C                                *
!   reference: ASHRAE Fundamentals Handbook - Chapter 6          *
! ****************************************************************
!
      real function fpws(tdb)

      use precision1
      implicit none
      real(kind=pp)                 :: tdb
      real(kind=pp), dimension(10)  :: a = (/52.26,71.1,95.5,126.68,166.08, &
                                           215.38,276.36,351.16,441.94,551.36/), &
                                       b = (/611.2,517.,273.,-194.7,-982.7, &
                                           -2215.2,-4044.6,-6662.6,-10293.8,-15217.7/)

      real(kind=pp)                 :: c1 = -5800.2206, c2 = 1.3914993, c3 = -0.048640239, &
                                       c4 = 0.41764768e-4, c5 = -0.14452093e-7, &
                                       c6 = 6.5459673
      real(kind=pp)                 :: td                                 
      integer                       :: it

      it=int(tdb/5.0) + 1
      if (it>0 .and. it<10) then
!  Temperature in range 0 - 50 C
!  Use linear approximations produced by GEK from ASHRAE, Chap 6, Table 1
          fpws=a(it)*tdb + b(it)
      elseif (it>=10 .and. it<40) then
!  50 < t < 200, outside range of piece-wise linear fit, use Chap 6, Eqn 4
          td=tdb+273.15
          fpws=exp(c1/td+c2+c3*td+c4*td**2.+c5*td**3.+c6*log(td))
      elseif (it<=0) then
!  t < 0 , out of range, use t=0 value and write warning
	  fpws=611.2
!	  write(*,*) "function fpws: temperature out of range, t=",tdb
      else
!  t > 200 , out of range, use t=200 value and write warning
	  fpws=1555074.
!	  write(*,*) "function fpws: temperature out of range, t=",tdb
      endif
      return
      end function fpws
!
! ********************************************************************
!   2. Dew point temperature or saturation temperature as a function *
!   of the humidity ratio
!   Range of 0 deg C to 70 deg C                                     *
!   reference: ASHRAE 1981 Fundamentals Handbook - Chapter 5              *
! ********************************************************************
!
      real function ftdew(w,patm)
      use precision1
      implicit none
      real(kind=pp)                 :: w,patm
      real(kind=pp)                 :: fpww, ps, alpha
      real(kind=pp)                 :: c0 = -35.957, c1 = -1.8726, c2 = 1.1689

      ps=fpww(w,patm)

      if (ps.lt.1.e-3) then
        ftdew=0.
        return
      endif

      alpha=log(ps)
      ftdew=c0+c1*alpha+c2*alpha**2

      return
      end function ftdew
!
! **********************************************************************
!   3. Partial pressure of water vapor as a function of the humidity   *
!   ratio and of the barometric pressure                               *
! **********************************************************************
!
      real function fpww(w,patm)
      use precision1
      implicit none
      real(kind=pp)                 :: w,patm

      fpww=patm*w/(0.62198+w)
      return
      end function fpww
!
! **********************************************************************
!   4. Humidity ratio as a function of the vapor pressure and the      *
!   barometric pressure                                                *
! **********************************************************************
!
      real function fwpw(pw,patm)
      use precision1
      implicit none
      real(kind=pp)                 :: pw,patm

      fwpw=0.62198*pw/(patm-pw)
      return
      end function fwpw
! **********************************************************************
!   5. Humidity ratio as a function of the air dry bulb temperature and*
!   the relative humidity                                              *
! **********************************************************************
!
      real function fwphi(tdb,phi,patm)
      use precision1
      implicit none
      real(kind=pp)                 :: tdb,phi,patm
      real(kind=pp)                 :: ps,pw,fpws,fwpw

      ps=fpws(tdb)
      pw=0.01*phi*ps
      fwphi=fwpw(pw,patm)
      return
      end function fwphi
!
! *********************************************************************
!   6. Humidity ratio as a function of the air dry bulb temperature   *
!   and of the air wet bulb temperature                               *
! *********************************************************************
!
      real function fwtwb(tdb,twb,patm)
      use precision1
      implicit none
      real(kind=pp)                 :: tdb,twb,patm
      real(kind=pp)                 :: ws,fwpw,pstwb,fpws
      real(kind=pp)                 :: cpa=1.006,cpg=1.805,hfg=2501.,cpw=4.187

      pstwb=fpws(twb)
      ws=fwpw(pstwb,patm)
      fwtwb=(ws*(hfg+(cpg-cpw)*twb)-cpa*(tdb-twb))/(hfg+cpg*tdb-cpw*twb)

      return
      end function fwtwb
!
! *********************************************************************
!   7. Humidity ratio as a function of the air dry bulb temperature   *
!   and of the air enthalpy                                            *
! *********************************************************************
!
      real function fwha(tdb,ha)
      use precision1
      implicit none
      real(kind=pp)                 :: tdb,ha
      real(kind=pp)                 :: cpa=1.006,cpg=1.805,hfg=2501.

      fwha=(ha-cpa*tdb)/(cpg*tdb+hfg)

      return
      end function fwha
!
! *********************************************************************
!   8. Air dry bulb temperature as a function of the humidity ratio   *
!   and of the air enthalpy                                            *
! *********************************************************************
!
      real function ftdb(w,ha)
      use precision1
      implicit none
      real(kind=pp)                 :: w,ha
      real(kind=pp)                 :: cpa=1.006,cpg=1.805,hfg=2501.

      ftdb=(ha-hfg*w)/(cpa+cpg*w)

      return
      end function ftdb
!
! **********************************************************************
!   9. Relative humidity as a function of the dry air temperature and  *
!   the humidity ratio                                              *
! **********************************************************************
!
      real function fphi(tdb,w,patm)
      use precision1
      implicit none
      real(kind=pp)                 :: tdb,w,patm
      real(kind=pp)                 :: pw,ps,fpww,fpws

      pw=fpww(w,patm)
      ps=fpws(tdb)

      if (ps==0.) then
        fphi=0.
      endif

      fphi=100*pw/ps
      return
      end function fphi
!
! *****************************
!   10. Enthalpy of moist air *
! *****************************
!
      real function fhair(tdb,w)
      use precision1
      implicit none
      real(kind=pp)                 :: tdb,w
      real(kind=pp)                 :: cpa=1.006,cpg=1.805,hfg=2501.

      fhair=cpa*tdb+w*(cpg*tdb+hfg)

      return
      end function fhair
!
! *********************************************************************
!   11. Air saturation enthalpy as function of saturation temperature *
! *********************************************************************
!
      real function fhsat(tsat,patm)
      use precision1
      implicit none
      real(kind=pp)                 :: tsat,patm
      real(kind=pp)                 :: ps, fpws, ws, fwpw, fhair

      ps=fpws(tsat)
      ws=fwpw(ps,patm)
      fhsat=fhair(tsat,ws)
      return
      end function fhsat
!
! **********************************************************************
!   12. Saturation temperature as a function of the air saturation     *
!   enthalpy                                                           *
!                                                                      *
!   Use piecewise linear approximation if 0 < ts < 47 and atmospheric  *
!   pressure within +/- 0.5% of standard atmospheric pressure          *
! **********************************************************************
!
      real function ftsat(hs,patm)
      use precision1
      implicit none
      real(kind=pp)                 :: hs,patm
      real(kind=pp)                 :: dp,ts1,ts2,hs1,hs2,delth1,delth2,fhsat,ts
      integer                       :: i

      real(kind=pp)                 :: c0=-6.0055,c1=0.68510,c2=-0.0056978,c3=3.5344e-5, &
                                       c4=-1.2891e-7,c5=2.0165e-10

      dp=abs(patm-101325.0)
      if (hs>9.473 .and. hs.lt.236.759 .and. dp.le.500.0) then
!
!  Temperature and pressure range OK for piecewise linear approximation
!  produced from ASHRAE Chap 6, Table 1 by GEK
!
          if (hs >= 9.4730 .and. hs .lt. 18.639) then
!-------- ( 0 < ts < 5 )
              ftsat = 0.545494*hs - 5.16747
          elseif (hs >= 18.639 .and. hs .lt. 31.724) then
!-------- ( 5 < ts < 11 )
              ftsat = 0.45854*hs - 3.546733
          elseif (hs >= 31.724 .and. hs .lt. 48.0) then
!-------- ( 11 < ts < 17 )
              ftsat = 0.36844*hs - 0.694765
          elseif (hs >= 48.0 .and. hs .lt. 68.44) then
!-------- ( 17 < ts < 23 )
              ftsat = 0.293542*hs + 2.90998
          elseif (hs >= 68.44 .and. hs .lt. 94.878) then
!-------- ( 23 < ts < 29 )
              ftsat = 0.226946*hs + 7.467811
          elseif (hs >= 94.878 .and. hs .lt. 129.455) then
!-------- ( 29 < ts < 35 )
              ftsat = 0.173526*hs + 12.536224
          elseif (hs >= 129.455 .and. hs .lt. 175.265) then
!-------- ( 35 < ts < 41 )
              ftsat = 0.130976*hs + 18.04453
          elseif (hs >= 175.265 .and. hs .lt. 236.759) then
!-------- ( 41 < ts < 47 )
              ftsat = 0.0975705*hs + 23.89931
	  endif
      elseif (hs < 9.473) then
!
!  Outside range of correlations
!
!	  write(*,*) 'function ftsat: enthalpy out of range, hs =',hs
          ftsat = 0.
      else
!
!  Temperature > 47 or pressure not within 0.5% of atmospheric, use
!  secant method to find inverse of fhsat
!
! ******************
!   Initial values *
! *****************************************************************
!   First value : fit on the ASHRAE values at barometric pressure*
! *****************************************************************
!
          ts1=c0+hs*(c1+hs*(c2+hs*(c3+hs*(c4+hs*c5))))
          hs1=fhsat(ts1,patm)
          delth1=hs1-hs
          if (abs(delth1)<1.e-3) then
              ftsat=ts1
          else
!   Iterate (why is guess temperature reduced by 5 C ?)
              ts2=ts1-5.
!
! ******************************
!   starting of the iterations *
! ******************************
!
              do i=1,50
                  hs2=fhsat(ts2,patm)
                  delth2=hs2-hs

                  if (abs(delth2)<1.e-3) exit

                  ts=ts1-delth1*(ts2-ts1)/(delth2-delth1)
                  ts1=ts2
                  hs1=hs2
                  delth1=delth2
                  ts2=ts
! error message if convergence fails
!	          if (i==50) write(*,*) &
!                      'function ftsat: convergence failure'
              enddo
              ftsat=ts2
          endif
      endif
      return
      end function ftsat
!
! **********************************************************************
!   13. Air wet bulb temperatures as a function of the air dry bulb    *
!   temperature,the humidity ratio and the barometric pressure         *
!                                                                      *
!   Determination of twb by the secant method                          *
! **********************************************************************
!
      real function ftwb(tdb,w,patm)

      use precision1
      implicit none
      real(kind=pp)          :: tdb,w,patm
      real(kind=pp)          :: h,fhair,twb1,ftsat,pws1,fpws,ws1,fwpw,h1,delth, &
                                twb2,pws2,ws2,h2,twb,delth1,delth2
      integer                :: i
      real(kind=pp)          :: cpw = 4.187
!
! *********************************************************************
!   Initial values - twb2 = saturation temperature corresponding to h *
!                    twb1 = twb2 - 1                                  *
! ******************************************************************** *

      h=fhair(tdb,w)

      twb1=ftsat(h,patm)-5.
      pws1=fpws(twb1)
      ws1=fwpw(pws1,patm)
      h1=fhair(twb1,ws1)-cpw*twb1*(ws1-w)
      delth1=h1-h

      twb2=twb1+5.
!
! ******************************
!   Starting of the iterations *
! ******************************
!
      do i=1,50
          pws2=fpws(twb2)
          ws2=fwpw(pws2,patm)
          h2=fhair(twb2,ws2)-cpw*twb2*(ws2-w)
          delth2=h2-h

          if(abs(delth2)<1.e-3) exit

          twb=twb1-delth1*(twb2-twb1)/(delth2-delth1)
          twb1=twb2
          h1=h2
          delth1=delth2
          twb2=twb
! error message if convergence fails
!	          if (i==50) write(*,*) &
!                     'function ftwb: convergence failure'
      enddo

      ftwb=twb2
      return
      end function ftwb
!
! ********************************
!   14. temperature of moist air *
! ********************************
!
      real function ftair(ha,w)
      use precision1
      implicit none
      real(kind=pp)                 :: ha,w
      real(kind=pp)                 :: cpa=1.006,cpg=1.805,hfg=2501.

      ftair=(ha-w*hfg)/(cpa+w*cpg)

      return
      end function ftair
! **********************************************************************
! * Copyright ASHRAE. Control Simulation Testbed
! **********************************************************************
!
!   LIBRARY:        Routines for use in communications 
!  
!   DEVELOPER:      Philip Haves and Li Mei
!                   Loughborough University of Technology
!
!   REVISION DATE:  January 29, 1996  
!
!   REFERENCE:      ASHRAE 825-RP Final Report
!
! ======================================================================
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        real function elapsed(iss,im,is0,im0)

        use precision1
        implicit none
        integer                :: iss,is0,im,im0
!
! *** Calculate elapsed time (real) from current and initial second and
! *** millisecond values (integer)
!
        elapsed=real(iss-is0) + 0.001*real(im-im0)
        return
        end function elapsed
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        subroutine tcorrect(is,im,t)

        use precision1
        implicit none
        real(kind=pp)          :: t
        integer                :: is,im,its,itsm,itm
!
! *** Correct time expressed in integer second (is) and millisecond (im) by
! *** subtracting a real time interval (t)
!
        its=int(t)
        itsm=nint(t*1000.0)
        itm=mod(itsm,1000)
        if (im>=itm) then
! *** No need to "carry"
            is=is-its
            im=im-itm
        else
! *** "Carry"
            is=is-its-1
            im=im-itm+1000
        endif
        return
        end subroutine tcorrect
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        subroutine rfile(nfile,fnbase,nval,itype,store)
!
        use precision1
        implicit none

! *** Read values of simulation variables from a file
!
! *** nval numbers are read from nval lines of file fnbase//nfile//.par
! *** into array store
!
        integer                         :: nfile,nval,itype,nlines,io,i
        real(kind=pp),dimension(nval)   :: store
        character(len=*)                :: fnbase
        character(len=10)               :: parfile
        character(len=1)                :: onedigit
        character(len=2)                :: twodigit
!
! *** No read if filenumber zero, else build filename
        if (nfile>0) then
            if (nfile<10) then
                write(onedigit,fmt='(i1)') nfile
                parfile=fnbase//onedigit//'.par '
            elseif (nfile<100) then
                write(twodigit,fmt='(i2)') nfile
                parfile=fnbase//twodigit//'.par'
            else
                write(*,*) 'rfile: type',itype,&
                           ' file number = ',nfile,' is out of range'
                stop
            endif

! *** Close file as a precaution, then attempt to open it
            close(20)
            open(unit=20,file=parfile,status='old',iostat=io)
            if(io>0) then
               write(*,*) 'rfile: type',itype,' file ',parfile,&
                          ' cannot be opened'
               stop
            endif   

! *** Read number of data and check for consistency
            read(20,*) nlines
            if (nlines /= nval) then
                write(*,*) 'rfile: ',parfile,&
                           ': wrong number of parameters in file'
                stop
            endif
! *** Read data
            do i=1,nval
                read(20,*) store(i)               
            enddo
            close(20)
        endif
        return
        end subroutine rfile

