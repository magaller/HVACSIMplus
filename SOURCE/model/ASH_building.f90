! ***********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! ***********************************************************************
! * SUBROUTINE:     Two time constant room model (no plenum)
! *
! * PURPOSE:        Calculate the temperature and humidity ratio in the
! *                 room by performing heat balances on the aggregated
! *                 structure node and the air node a moisture balance
! *                 on the air node. Include the effects of air flows
! *                 from adjacent zones and any (inward) air leakage
! *                 and reverse return flow inaddition to the HVAC supply
! *                 air flow.
! *
! ***********************************************************************
! * INPUTS
! * ======
! *  1. tsup    : Supply air dry bulb temperature                    (C)
! *  2. gsup    : Supply air humidity ratio                      (kg/kg)
! *  3. wsup    : Supply dry air mass flow rate                   (kg/s)
! *  4. tinz1   : Interzone 1 air dry bulb temperature               (C)
! *  5. ginz1   : Interzone 1 air humidity ratio                 (kg/kg)
! *  6. winz1   : Interzone 1 dry air mass flow rate              (kg/s)
! *  7. tinz2   : Interzone 2 air dry bulb temperature               (C)
! *  8. ginz2   : Interzone 2 air humidity ratio                 (kg/kg)
! *  9. winz2   : Interzone 2 dry air mass flow rate              (kg/s)
! * 10. wret    : Extract dry air mass flow rate                  (kg/s)
! * 11. tamb    : Ambient dry bulb temperature                       (C)
! * 12. gamb    : Ambient humidity ratio                         (kg/kg)
! * 13. wleak   : Leakage air mass flow rate (positive out)      (kg/kg)
! * 14. tsolairr: Equivalent "sol-air" outdoor temperature for room  (C)
! * 15. fracocc : Fractional occupancy                               (-)
! * 16. fraclig : Fractional lighting heat gain                      (-)
! * 17. fracpow : Fractional equipment heat gain                     (-)
! * 18. qsduct  : Heat gain from supply duct                        (kW)
! * 19. qeduct  : Heat gain from return duct                        (kW)
! * 20. troom   : Room temperature                                   (C)
! * 21. tsroom  : Room structure temperature                         (C)
! * 22. groom   : Room humidity ratio                            (kg/kg)
! *
! * OUTPUTS
! * =======
! *  1. troomn  : Room temperature                                   (C)
! *  2. tsroomn : Room structure temperature                         (C)
! *  3. groomn  : Room humidity ratio                            (kg/kg)
! *  4. tret    : Return temperature                                 (C)
! *  5. qsensr  : Sensible heat gains of room                       (kW)
! *  6. qsensp  : Sensible heat gains of plenum                     (kW)
! *  7. wvapr   : Water vapour gains of room                      (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. xcap    : Room air capacity multiplier                       (-)
! *  2. rwsr    : Direct resistance room air node <-> ambient     (K/kW)
! *  3. risr    : Resistance room air node <-> room mass node     (K/kW)
! *  4. rosr    : Resistance ambient <-> room mass node           (K/kW)
! *  5. csr     : Capacitance of room mass node                   (kJ/K)
! *  6. crair   : Capacitance of room air node (unmodified)       (kJ/K)
! *  7. vroom   : Volume of room                                    (m3)
! *  8. noccup  : Number of occupants                                (-)
! *  9. qlite   : Lighting heat gain                                (kW)
! * 10. flpln   : Fraction of lighting heat gain to return air       (-)
! * 11. qequp   : Equipment heat gain                               (kW)
! * 12. nfile   : Zone number (parameter file='ZONEn.PAR', n > 0)    (-)
! *
! * SAVED
! * =====
! *  1. time    : Time of previous call of TYPE
! *  2. troom   : Room temperature from previous call
! *  3. troomp  : Room temperature from previous step time
! *  4. tsroom  : Room structure temperature from previous call
! *  5. tsroomp : Room structure temperature from previous step time
! *  6. groom   : Room humidity ratio from previous call
! *  7. groomp  : Room humidity ratio from previous step time
! *  8. qartif  : Heat input for commissioning (read from file)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 29, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  RFILE, DIFFEQ
!   FUNCTIONS  CALLED:
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! * mroom   : mass of air in room
! * mplen   : mass of air in plenum
! * cpa     : specific heat of dry air
! * cpg     : specific heat of water vapor
! * hfg     : latent heat of water vapor
! * rhoair  : density of air
! * qsperson: sensible heat gain from a person
! * qlperson: latent heat gain from a person
! * wsmall  : threshold for significant air flow rate
! * diwsup  : value of supply air flow rate if flowing in, zero if out
! * diwinz1 : value of air flow rate from Adjacent Zone 1 if flowing in
! * diwinz2 : value of air flow rate from Adjacent Zone 2 if flowing in
! * diwleak : value of leakage air flow rate if flowing in, zero if out
! * csup    : value of supply air capacity rate if flowing in, zero if out
! * cinz1   : value of air capacity rate from Adjacent Zone 1 if flowing in
! * cinz2   : value of air capacity rate from Adjacent Zone 2 if flowing in
! * cleak   : value of leakage air capacity rate if flowing in, zero if out
! * cret    : extract air capacity rate
! * qoccs   : sensible gains from occupants
! * wvapr   : latent gains from occupants
! * qrlig   : heat gain to room from lights
! * qelig   : heat gain to extract air stream from lights
! * qplig   : heat gain to plenum from lights
! * qpow    : heat gain to room from equipment
! * aa      : A coefficent in dT/dt = A*T + B
! * bb      : B coefficent in dT/dt = A*T + B
! * xxxxxb  : average over time-step value of integrated variable (not used)
! *
! **********************************************************************
!
!   Updated to Fortran 90 April 27, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type401(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndiffeq=3,nsv=(1+ndiffeq*2),&
                                             ni=22,no=7,np=12,nfp=1,&
                                             ns=nsv+nfp
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

        real         :: mroom,noccup
        real         :: cpa=1.006,cpg=1.805,hfg=2501.,rhoair=1.2
        real         :: qsperson=0.095,qlperson=0.045,wsmall=1.e-4
        integer      :: itype=401

        real         :: tsup,troom,tsroom,groom,xcap,rwsr,risr,rosr,&
                        csr,crair,vroom,qlite,flpln,qequp,cr,&
                        qartif,troomp,tsroomp,groomp,wfan,diwsup,&
                        diwinz1,diwinz2,diwleak,diwret,diwfan,csup,&
                        cinz1,cinz2,cleak,cret,cfan,cretl,qoccs,wvapr,&
                        qrlig,qelig,qpow,qsensr,aa,bb,troomb,troomn,&
                        tsroomb,tsroomn,groomb,groomn,tret,gsup,wsup,&
                        tinz1,ginz1,winz1,tinz2,ginz2,winz2,wret,tamb,&
                        gamb,wleak,tsolair,fracocc,fraclig,fracpow,&
                        qsduct,qeduct,tsolairr
        integer      :: i,is,nfile

! **** Read in inputs
        tsup     = xin(1)
        gsup     = xin(2)
        wsup     = xin(3)
        tinz1    = xin(4)
        ginz1    = xin(5)
        winz1    = xin(6)
        tinz2    = xin(7)
        ginz2    = xin(8)
        winz2    = xin(9)
        wret     = xin(10)
        tamb     = xin(11)
        gamb     = xin(12)
        wleak    = xin(13)
        tsolairr = xin(14)
        fracocc  = xin(15)
        fraclig  = xin(16)
        fracpow  = xin(17)
        qsduct   = xin(18)
        qeduct   = xin(19)
        troom    = xin(20)
        tsroom   = xin(21)
        groom    = xin(22)
! **** Read in parameters
        xcap     = par_v(1)
        rwsr     = par_v(2)
        risr     = par_v(3)
        rosr     = par_v(4)
        csr      = par_v(5)
        crair    = par_v(6)
        vroom    = par_v(7)
        noccup   = par_v(8)
        qlite    = par_v(9)
        flpln    = par_v(10)
        qequp    = par_v(11)
        nfile    = nint(par_v(12))
        cr       = crair*xcap
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Initialize room conditions to 20 deg c, 50% rh
                do is = 2,4,2
                    saved_v(is) = 20.0
                enddo
                do is = 6,6,2
                    saved_v(is) = 0.0074
                enddo
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update previous sample instant values
            do is = 2,nsv-1,2
                saved_v(is+1) = saved_v(is)
            enddo
!    Open zone heat gain/loss file - error -> zero extra gain
            call rfile(nfile,'zone',nfp,itype,fpar)
            do i=1,nfp
                saved_v(nsv+i)=fpar(i)
            enddo
        endif
        qartif = saved_v(nsv+1)
! **** Update previous values
        troomp  = saved_v(3)
        tsroomp = saved_v(5)
        groomp  = saved_v(7)
! **** Set up and solve heat and moisture balances
! **** Calculate local extract from mass balance
        wfan = wsup + winz1 - winz2 - wleak - wret
! **** Calculate capacity rates entering room and plenum.
! **** nb wsup and winz1 are positive entering the room, winz2, wleak, wret and
! **** wfan are positive leaving the room. Capacity rates are non-zero if the
! **** flow is into the zone.
        diwsup  = max(0.0,wsup)
        diwinz1 = max(0.0,winz1)
        diwinz2 = max(0.0,-winz2)
        diwleak = max(0.0,-wleak)
        diwret  = max(0.0,-wret)
        diwfan  = max(0.0,-wfan)
        csup    = diwsup*(cpa+gsup*cpg)
        cinz1   = diwinz1*(cpa+ginz1*cpg)
        cinz2   = diwinz2*(cpa+ginz2*cpg)
        cleak   = diwleak*(cpa+gamb*cpg)
        cret    = diwret*(cpa+gamb*cpg)
        cfan    = diwfan*(cpa+gamb*cpg)
        cretl   = wret*(cpa+groom*cpg)
! **** Calculate internal gains
! **** Sensible gains from occupants
        qoccs = fracocc*noccup*qsperson
! **** Latent gains from occupants
        wvapr = (fracocc*noccup*qlperson)/hfg
        if (wret>wsmall) then
! **** Flow through luminaire - divide heat between room and return air
            qrlig = fraclig*(1.-flpln)*qlite
            qelig = fraclig*flpln*qlite
        else
! **** No flow through luminaire - 100% to room
            qrlig = fraclig*qlite
            qelig = 0.0
        endif
! **** Equipment
        qpow   = fracpow*qequp
! **** Sensible gains to room and plenum
        qsensr = qoccs+qrlig+qpow+qartif+qsduct+qeduct
! **** Heat balance on room air - reverse flow in return is assumed to be from
! **** ambient
        aa = -(1./cr)*(1./risr+1./rwsr+csup+cinz1+cinz2+cleak+cret+cfan)
        bb =  (1./cr)*(tsroom/risr+tsolairr/rwsr+&
                       tsup*csup+tinz1*cinz1+tinz2*cinz2+&
                       tamb*(cleak+cret+cfan)+qsensr)
        call diffeq(time,aa,bb,troomp,troomn,troomb)
! **** Heat balance on room structure
        aa = -(1./csr)*(1./risr+1./rosr)
        bb =  (1./csr)*(troom/risr+tsolairr/rosr)
        call diffeq(time,aa,bb,tsroomp,tsroomn,tsroomb)
! **** Moisture balance on room air - reverse flow in return is assumed to be ! ** <<<<<<
        mroom = vroom*rhoair
        aa = -(1./mroom)*(diwsup+diwinz1+diwinz2+diwleak+diwret+diwfan)
        bb =  (1./mroom)*(gsup*diwsup+ginz1*diwinz1+ginz2*diwinz2+&
                          gamb*(diwleak+diwret+diwfan)+wvapr)
        call diffeq(time,aa,bb,groomp,groomn,groomb)
! **** Extract temperature
        if (wret>wsmall) then
! **** Flow through luminaire - calculate heat pick-up
            tret = troom+qelig/cretl
        else
! **** No flow - temperature indeterminate
            tret = troom
        endif
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = troomn
        saved_v(4)  = tsroomn
        saved_v(6)  = groomn
! **** Save time of current call
        saved_v(1)  = time
! **** Output
        yout(1) = troomn
        yout(2) = tsroomn
        yout(3) = groomn
        yout(4) = tret
        yout(5) = qsensr
        yout(6) = wvapr
! **** Disallow freezing of dynamic variables
        do i=1,3
            iostat(i) = 0
        enddo
! **** Allow freezing of algebraic variables
        do i=4,6
            iostat(i) = 1
        enddo
! **** Return
        return
        end subroutine type401
! ***********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! ***********************************************************************
! * SUBROUTINE:     Room with plenum return
! *
! * PURPOSE:        Calculate the temperature and humidity ratio in the
! *                 room and in the plenum by performing heat
! *                 balances on the aggregated structure node and the
! *                 air node for each space and moisture balances on the
! *                 air node in each space. Include the effects of air
! *                 flows from adjacent zones and any (inward) air leakage
! *                 and reverse return flow in addition to the HVAC supply
! *                 air flow. Include supply duct heat gains/losses
! *                 from/to the plenum.
! *
! ***********************************************************************
! * INPUTS
! * ======
! *  1. tsup    : Supply air dry bulb temperature                    (C)
! *  2. gsup    : Supply air humidity ratio                      (kg/kg)
! *  3. wsup    : Supply dry air mass flow rate                   (kg/s)
! *  4. tinz1   : Interzone 1 air dry bulb temperature               (C)
! *  5. ginz1   : Interzone 1 air humidity ratio                 (kg/kg)
! *  6. winz1   : Interzone 1 dry air mass flow rate              (kg/s)
! *  7. tinz2   : Interzone 2 air dry bulb temperature               (C)
! *  8. ginz2   : Interzone 2 air humidity ratio                 (kg/kg)
! *  9. winz2   : Interzone 2 dry air mass flow rate              (kg/s)
! * 10. wret    : Extract dry air mass flow rate                  (kg/s)
! * 11. tamb    : Ambient dry bulb temperature                       (C)
! * 12. gamb    : Ambient humidity ratio                         (kg/kg)
! * 13. wleak   : Leakage air mass flow rate (positive out)      (kg/kg)
! * 14. tsolairr: Equivalent "sol-air" outdoor temperature for room  (C)
! * 15. tsolairp: Equivalent "sol-air" outdoor temperature for plenum(C)
! * 16. fracocc : Fractional occupancy                               (-)
! * 17. fraclig : Fractional lighting heat gain                      (-)
! * 18. fracpow : Fractional equipment heat gain                     (-)
! * 19. qsduct  : Heat gain from supply duct                        (kW)
! * 20. troom   : Room temperature                                   (C)
! * 21. tsroom  : Room structure temperature                         (C)
! * 22. tplen   : Plenum temperature                                 (C)
! * 23. tsplen  : Plenum structure temperature                       (C)
! * 24. groom   : Room humidity ratio                            (kg/kg)
! * 25. gplen   : Plenum humidity ratio                          (kg/kg)
! *
! * OUTPUTS
! * =======
! *  1. troom   : Room temperature                                   (C)
! *  2. tsroom  : Room structure temperature                         (C)
! *  3. tplen   : Plenum temperature                                 (C)
! *  4. tsplen  : Plenum structure temperature                       (C)
! *  5. groom   : Room humidity ratio                            (kg/kg)
! *  6. gplen   : Plenum humidity ratio                          (kg/kg)
! *  7. tret    : Return temperature                                 (C)
! *  8. qsensr  : Sensible heat gains of room                       (kW)
! *  9. qsensp  : Sensible heat gains of plenum                     (kW)
! * 10. wvapr   : Water vapour gains of room                      (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. xcap    : Room air capacity multiplier                       (-)
! *  2. rwsr    : Direct resistance room air node <-> ambient     (K/kW)
! *  3. risr    : Resistance room air node <-> room mass node     (K/kW)
! *  4. rosr    : Resistance ambient <-> room mass node           (K/kW)
! *  5. rwsp    : Direct resistance plenum air node <-> ambient   (K/kW)
! *  6. risp    : Resistance plenum air node <-> plenum mass node (K/kW)
! *  7. rosp    : Resistance ambient <-> plenum mass node         (K/kW)
! *  8. rr      : Resistance room air node <-> plenum air node    (K/kW)
! *  9. csr     : Capacitance of room mass node                   (kJ/K)
! * 10. crair   : Capacitance of room air node (unmodified)       (kJ/K)
! * 11. csp     : Capacitance of plenum mass node                 (kJ/K)
! * 12. cp      : Capacitance of plenum air node                  (kJ/K)
! * 13. vroom   : Volume of room                                    (m3)
! * 14. vplen   : Volume of plenum                                  (m3)
! * 15. noccup  : Number of occupants                                (-)
! * 16. qlite   : Lighting heat gain                                (kW)
! * 17. flpln   : Fraction of lighting heat gain to extract air      (-)
! * 18. qequp   : Equipment heat gain                               (kW)
! * 19. nfile   : Zone number (parameter file='zoneN.par', n > 0)    (-)
! *
! * SAVED
! * =====
! *  1. time    : Time of previous call of TYPE
! *  2. troom   : Room temperature from previous call
! *  3. troomp  : Room temperature from previous step time
! *  4. tsroom  : Room structure temperature from previous call
! *  5. tsroomp : Room structure temperature from previous step time
! *  6. tplen   : Plenum temperature from previous call
! *  7. tplenp  : Plenum temperature from previous step time
! *  8. tsplen  : Plenum structure temperature from previous call
! *  9. tsplenp : Plenum structure temperature from previous step time
! * 10. groom   : Room humidity ratio from previous call
! * 11. groomp  : Room humidity ratio from previous step time
! * 12. gplen   : Plenum humidity ratio from previous call
! * 13. gplenp  : Plenum humidity ratio from previous step time
! * 14. qartif  : Heat input for commissioning (read from file)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 29, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  RFILE, DIFFEQ
!   FUNCTIONS  CALLED:
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! * mroom   : mass of air in room
! * mplen   : mass of air in plenum
! * cpa     : specific heat of dry air
! * cpg     : specific heat of water vapor
! * hfg     : latent heat of water vapor
! * rhoair  : density of air
! * qsperson: sensible heat gain from a person
! * qlperson: latent heat gain from a person
! * wsmall  : threshold for significant air flow rate
! * diwsup  : value of supply air flow rate if flowing in, zero if out
! * diwinz1 : value of air flow rate from Adjacent Zone 1 if flowing in
! * diwinz2 : value of air flow rate from Adjacent Zone 2 if flowing in
! * diwleak : value of leakage air flow rate if flowing in, zero if out
! * csup    : value of supply air capacity rate if flowing in, zero if out
! * cinz1   : value of air capacity rate from Adjacent Zone 1 if flowing in
! * cinz2   : value of air capacity rate from Adjacent Zone 2 if flowing in
! * cleak   : value of leakage air capacity rate if flowing in, zero if out
! * cret    : extract air capacity rate
! * qoccs   : sensible gains from occupants
! * wvapr   : latent gains from occupants
! * qrlig   : heat gain to room from lights
! * qelig   : heat gain to extract air stream from lights
! * qplig   : heat gain to plenum from lights
! * qpow    : heat gain to room from equipment
! * aa      : A coefficent in dT/dt = A*T + B
! * bb      : B coefficent in dT/dt = A*T + B
! * xxxxxb  : average over time-step value of integrated variable (not used)
! *
! **********************************************************************
!
!   Updated to Fortran 90 April 27, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type402(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndiffeq=6,nsv=(1+ndiffeq*2),&
                                             ni=25,no=10,np=19,nfp=1,&
                                             ns=nsv+nfp
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

        real         :: mroom,noccup
        real         :: cpa=1.006,cpg=1.805,hfg=2501.,rhoair=1.2
        real         :: qsperson=0.095,qlperson=0.045,wsmall=1.e-4
        integer      :: itype=402

        real         :: qlite,flpln,qequp,cr,qartif,troomp,&
                        tsroomp,tplenp,tsplenp,groomp,gplenp,wfan,&
                        diwsup,diwinz1,diwinz2,diwleak,diwret,diwfan,&
                        csup,cinz1,cinz2,cleak,cret,cfan,cretp,qoccs,&
                        wvapr,qrlig,qelig,qplig,qpow,qsensr,qsensp,&
                        aa,bb,troomb,troomn,tsroomb,tsroomn,tplenb,&
                        tplenn,tsplenb,tsplenn,groomb,groomn,mplen,&
                        gplenb,gplenn,tret,tsup,gsup,wsup,tinz1,ginz1,&
                        winz1,tinz2,ginz2,winz2,wret,tamb,gamb,wleak,&
                        tsolairr,tsolairp,fracocc,fraclig,fracpow,&
                        qsduct,troom,tsroom,tplen,tsplen,groom,gplen,&
                        xcap,rwsr,risr,rosr,rwsp,risp,rosp,rrp,csr,&
                        cp,vroom,vplen,crair,csp,qeduct
        integer      :: i,is,nfile


! **** Read in inputs
        tsup     = xin(1)
        gsup     = xin(2)
        wsup     = xin(3)
        tinz1    = xin(4)
        ginz1    = xin(5)
        winz1    = xin(6)
        tinz2    = xin(7)
        ginz2    = xin(8)
        winz2    = xin(9)
        wret     = xin(10)
        tamb     = xin(11)
        gamb     = xin(12)
        wleak    = xin(13)
        tsolairr = xin(14)
        tsolairp = xin(15)
        fracocc  = xin(16)
        fraclig  = xin(17)
        fracpow  = xin(18)
        qsduct   = xin(19)
        troom    = xin(20)
        tsroom   = xin(21)
        tplen    = xin(22)
        tsplen   = xin(23)
        groom    = xin(24)
        gplen    = xin(25)
! **** Read in parameters
        xcap     = par_v(1)
        rwsr     = par_v(2)
        risr     = par_v(3)
        rosr     = par_v(4)
        rwsp     = par_v(5)
        risp     = par_v(6)
        rosp     = par_v(7)
        rrp      = par_v(8)
        csr      = par_v(9)
        crair    = par_v(10)
        csp      = par_v(11)
        cp       = par_v(12)
        vroom    = par_v(13)
        vplen    = par_v(14)
        noccup   = par_v(15)
        qlite    = par_v(16)
        flpln    = par_v(17)
        qequp    = par_v(18)
        nfile    = nint(par_v(19))
        cr       = crair*xcap
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Initialize room conditions to 20 deg c, 50% rh
                do is = 2,8,2
                    saved_v(is) = 20.0
                enddo
                do is = 10,12,2
                    saved_v(is) = 0.0074
                enddo
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update previous sample instant values
            do is = 2,nsv-1,2
                saved_v(is+1) = saved_v(is)
            enddo
!    Open zone heat gain/loss file - error -> zero extra gain
            call rfile(nfile,'zone',nfp,itype,fpar)
            do i=1,nfp
                saved_v(nsv+i)=fpar(i)
            enddo
        endif
        qartif = saved_v(nsv+1)
! **** Update previous values
        troomp  = saved_v(3)
        tsroomp = saved_v(5)
        tplenp  = saved_v(7)
        tsplenp = saved_v(9)
        groomp  = saved_v(11)
        gplenp  = saved_v(13)
! **** Set up and solve heat and moisture balances
! **** Calculate local extract from mass balance
        wfan = wsup + winz1 - winz2 - wleak - wret
! **** Calculate capacity rates entering room and plenum.
! **** nb wsup and winz1 are positive entering the room, winz2, wleak, wret and
! **** wfan are positive leaving the room. capacity rates are non-zero if
! **** the flow is into the zone. All leakage is assumed to be from the zone
! **** rather than the plenum.
        diwsup  = max(0.0,wsup)
        diwinz1 = max(0.0,winz1)
        diwinz2 = max(0.0,-winz2)
        diwleak = max(0.0,-wleak)
        diwret  = max(0.0,-wret)
        diwfan  = max(0.0,-wfan)
        csup    = diwsup*(cpa+gsup*cpg)
        cinz1   = diwinz1*(cpa+ginz1*cpg)
        cinz2   = diwinz2*(cpa+ginz2*cpg)
        cleak   = diwleak*(cpa+gamb*cpg)
        cret    = diwret*(cpa+gamb*cpg)
        cfan    = diwfan*(cpa+gamb*cpg)
        cretp   = wret*(cpa+groom*cpg)
! **** Calculate internal gains
! **** Sensible gains from occupants
        qoccs = fracocc*noccup*qsperson
! **** Latent gains from occupants
        wvapr = (fracocc*noccup*qlperson)/hfg
        if (wret>wsmall) then
! **** Flow through luminaire - divide heat between room and extract air
            qrlig = fraclig*(1.-flpln)*qlite
            qelig = fraclig*flpln*qlite
            qplig = 0.0
        else
! **** No flow through luminaire - 50% to room, 50% to plenum
            qrlig = fraclig*0.5*qlite
            qplig = fraclig*0.5*qlite
            qelig = 0.0
        endif
! **** Equipment
        qpow   = fracpow*qequp
! **** Sensible gains to room and plenum
        qsensr = qoccs+qrlig+qpow+qartif
        qsensp = qsduct+qplig
! **** Heat balance on room air
        aa = -(1./cr)*(1./risr+1./rwsr+1./rrp+csup+cinz1+cinz2+cleak+&
                       cret+cfan)
        bb =  (1./cr)*(tsroom/risr+tsolairr/rwsr+tplen/rrp+&
                       tsup*csup+tinz1*cinz1+tinz2*cinz2+&
                       tamb*(cleak+cret+cfan)+qsensr)
        call diffeq(time,aa,bb,troomp,troomn,troomb)
! **** Heat balance on room structure
        aa = -(1./csr)*(1./risr+1./rosr)
        bb =  (1./csr)*(troom/risr+tsolairr/rosr)
        call diffeq(time,aa,bb,tsroomp,tsroomn,tsroomb)
! **** Heat balance on plenum air
        aa = -(1./cp)*(1./risp+1./rwsp+1./rrp+cretp)
        bb =  (1./cp)*(tsplen/risp+tsolairp/rwsp+troom/rrp+troom*cretp+&
                       qsensp)
        call diffeq(time,aa,bb,tplenp,tplenn,tplenb)
! **** Heat balance on plenum structure
        aa = -(1./csp)*(1./risp+1./rosp)
        bb =  (1./csp)*(tplen/risp+tsolairp/rosp)
        call diffeq(time,aa,bb,tsplenp,tsplenn,tsplenb)
! **** Moisture balance on room air
        mroom = vroom*rhoair
        aa = -(1./mroom)*(diwsup+diwinz1+diwinz2+diwleak+diwret+diwfan)
        bb =  (1./mroom)*(gsup*diwsup+ginz1*diwinz1+ginz2*diwinz2+&
                          gamb*(diwleak+diwret+diwfan)+wvapr)
        call diffeq(time,aa,bb,groomp,groomn,groomb)
! **** Moisture balance on plenum air
        mplen = vplen*rhoair
        aa = -(1./mplen)*wret
        bb =  (1./mplen)*groom*wret
        call diffeq(time,aa,bb,gplenp,gplenn,gplenb)
! **** Extract temperature
        if (wret>wsmall) then
! **** Flow through luminaire - calculate heat pick-up
            tret = troom+qelig/cretp
        else
! **** No flow - temperature indeterminate
            tret = troom
        endif
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = troomn
        saved_v(4)  = tsroomn
        saved_v(6)  = tplenn
        saved_v(8)  = tsplenn
        saved_v(10) = groomn
        saved_v(12) = gplenn
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1)  = troomn
        yout(2)  = tsroomn
        yout(3)  = tplenn
        yout(4)  = tsplenn
        yout(5)  = groomn
        yout(6)  = gplenn
        yout(7)  = tret
        yout(8)  = qsensr
        yout(9)  = qsensp
        yout(10) = wvapr
! **** Disallow freezing of dynamic variables
        do i=1,6
            iostat(i) = 0
        enddo
! **** Allow freezing of algebraic variables
        do i=7,10
            iostat(i) = 1
        enddo
! **** Return
        return
        end subroutine type402

! ***********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! ***********************************************************************
! * SUBROUTINE:     Room with plenum and ducted return (analytical
! *                 integration)
! *
! * PURPOSE:        Calculate the temperature and humidity ratio in the
! *                 room and in the plenum by performing heat
! *                 balances on the aggregated structure node and the
! *                 air node for each space and moisture balances on the
! *                 air node in each space. Include the effects of air
! *                 flows from adjacent zones and any (inward) air leakage
! *                 and reverse return flow in addition to the HVAC supply
! *                 air flow. Include duct heat gains/losses from/to the
! *                 plenum. Use DIFFEQ to peform analytical integration.
! *
! ***********************************************************************
! * INPUTS
! * ======
! *  1. tsup    : Supply air dry bulb temperature                    (C)
! *  2. gsup    : Supply air humidity ratio                      (kg/kg)
! *  3. wsup    : Supply dry air mass flow rate                   (kg/s)
! *  4. tinz1   : Interzone 1 air dry bulb temperature               (C)
! *  5. ginz1   : Interzone 1 air humidity ratio                 (kg/kg)
! *  6. winz1   : Interzone 1 dry air mass flow rate              (kg/s)
! *  7. tinz2   : Interzone 2 air dry bulb temperature               (C)
! *  8. ginz2   : Interzone 2 air humidity ratio                 (kg/kg)
! *  9. winz2   : Interzone 2 dry air mass flow rate              (kg/s)
! * 10. wret    : Extract dry air mass flow rate                  (kg/s)
! * 11. tamb    : Ambient dry bulb temperature                       (C)
! * 12. gamb    : Ambient humidity ratio                         (kg/kg)
! * 13. wleak   : Leakage air mass flow rate (positive out)      (kg/kg)
! * 14. tsolairr: Equivalent "sol-air" outdoor temperature for room  (C)
! * 15. tsolairp: Equivalent "sol-air" outdoor temperature for plenum(C)
! * 16. fracocc : Fractional occupancy                               (-)
! * 17. fraclig : Fractional lighting heat gain                      (-)
! * 18. fracpow : Fractional equipment heat gain                     (-)
! * 19. qsduct  : Heat gain from supply duct                        (kW)
! * 20. qeduct  : Heat gain from return duct                        (kW)
! * 21. troom   : Room temperature                                   (C)
! * 22. tsroom  : Room structure temperature                         (C)
! * 23. tplen   : Plenum temperature                                 (C)
! * 24. tsplen  : Plenum structure temperature                       (C)
! * 25. groom   : Room humidity ratio                            (kg/kg)
! *
! * OUTPUTS
! * =======
! *  1. troomn  : Room temperature                                   (C)
! *  2. tsroomn : Room structure temperature                         (C)
! *  3. tplenn  : Plenum temperature                                 (C)
! *  4. tsplenn : Plenum structure temperature                       (C)
! *  5. groomn  : Room humidity ratio                            (kg/kg)
! *  6. tret    : Return temperature                                 (C)
! *  7. qsensr  : Sensible heat gains of room                       (kW)
! *  8. qsensp  : Sensible heat gains of plenum                     (kW)
! *  9. wvapr   : Water vapour gains of room                      (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. xcap    : Room air capacity multiplier                       (-)
! *  2. rwsr    : Direct resistance room air node <-> ambient     (K/kW)
! *  3. risr    : Resistance room air node <-> room mass node     (K/kW)
! *  4. rosr    : Resistance ambient <-> room mass node           (K/kW)
! *  5. rwsp    : Direct resistance plenum air node <-> ambient   (K/kW)
! *  6. risp    : Resistance plenum air node <-> plenum mass node (K/kW)
! *  7. rosp    : Resistance ambient <-> plenum mass node         (K/kW)
! *  8. rr      : Resistance room air node <-> plenum air node    (K/kW)
! *  9. csr     : Capacitance of room mass node                   (kJ/K)
! * 10. crair   : Capacitance of room air node (unmodified)       (kJ/K)
! * 11. csp     : Capacitance of plenum mass node                 (kJ/K)
! * 12. cp      : Capacitance of plenum air node                  (kJ/K)
! * 13. vroom   : Volume of room                                    (m3)
! * 14. vplen   : Volume of plenum                                  (m3)
! * 15. noccup  : Number of occupants                                (-)
! * 16. qlite   : Lighting heat gain                                (kW)
! * 17. flpln   : Fraction of lighting heat gain to return air       (-)
! * 18. qequp   : Equipment heat gain                               (kW)
! * 19. nfile   : Zone number (parameter file='zoneN.par', N > 0)    (-)
! *
! * SAVED
! * =====
! *  1. time    : Time of previous call of TYPE
! *  2. troomn  : Room temperature from previous call
! *  3. troomp  : Room temperature from previous step time
! *  4. tsroom  : Room structure temperature from previous call
! *  5. tsroomp : Room structure temperature from previous step time
! *  6. tplenn  : Plenum temperature from previous call
! *  7. tplenp  : Plenum temperature from previous step time
! *  8. tsplen  : Plenum structure temperature from previous call
! *  9. tsplenp : Plenum structure temperature from previous step time
! * 10. groom   : Room humidity ratio from previous call
! * 11. groomp  : Room humidity ratio from previous step time
! * 12. qartif  : Heat input for commissioning (read from file)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 29, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  RFILE, DIFFEQ
!   FUNCTIONS  CALLED:
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! * mroom   : mass of air in room
! * mplen   : mass of air in plenum
! * cpa     : specific heat of dry air
! * cpg     : specific heat of water vapor
! * hfg     : latent heat of water vapor
! * rhoair  : density of air
! * qsperson: sensible heat gain from a person
! * qlperson: latent heat gain from a person
! * wsmall  : threshold for significant air flow rate
! * diwsup  : value of supply air flow rate if flowing in, zero if out
! * diwinz1 : value of air flow rate from Adjacent Zone 1 if flowing in
! * diwinz2 : value of air flow rate from Adjacent Zone 2 if flowing in
! * diwleak : value of leakage air flow rate if flowing in, zero if out
! * csup    : value of supply air capacity rate if flowing in, zero if out
! * cinz1   : value of air capacity rate from Adjacent Zone 1 if flowing in
! * cinz2   : value of air capacity rate from Adjacent Zone 2 if flowing in
! * cleak   : value of leakage air capacity rate if flowing in, zero if out
! * cret    : return air capacity rate
! * qoccs   : sensible gains from occupants
! * wvapr   : latent gains from occupants
! * qrlig   : heat gain to room from lights
! * qelig   : heat gain to return air stream from lights
! * qplig   : heat gain to plenum from lights
! * qpow    : heat gain to room from equipment
! * aa      : A coefficent in dT/dt = A*T + B
! * bb      : B coefficent in dT/dt = A*T + B
! * xxxxxb  : average over time-step value of integrated variable (not used)
! *
! **********************************************************************
!
!   Updated to Fortran 90 April 27, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type403(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndiffeq=5,nsv=(1+ndiffeq*2),&
                                             ni=25,no=9,np=19,nfp=1,&
                                             ns=nsv+nfp
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

        real         :: mroom,noccup
        real         :: cpa=1.006,cpg=1.805,hfg=2501.,rhoair=1.2
        real         :: qsperson=0.095,qlperson=0.045,wsmall=1.e-4
        integer      :: itype=403

        real         :: vroom,vplen,qlite,flpln,qequp,cr,qartif,&
                        troomp,tsroomp,tplenp,tsplenp,groomp,wfan,&
                        diwsup,diwinz1,diwinz2,diwleak,diwret,diwfan,&
                        csup,cinz1,cinz2,cleak,cret,cfan,cretp,qoccs,&
                        wvapr,qrlig,qelig,qplig,qpow,qsensr,qsensp,&
                        aa,bb,troomb,troomn,tsroomb,tsroomn,tplenb,&
                        tplenn,tsplenb,tsplenn,groomb,groomn,tret,tsup,&
                        gsup,wsup,tinz1,ginz1,winz1,tinz2,ginz2,winz2,&
                        wret,tamb,gamb,wleak,tsolairr,tsolairp,fracocc ,&
                        fraclig,fracpow ,qsduct,qeduct,troom,tsroom,&
                        tplen,tsplen,groom,xcap,rwsr,risr,rosr,rwsp,risp,&
                        rosp,rrp,csr,crair,csp,cp
        integer      :: i,is,nfile

! **** Read in inputs
        tsup     = xin(1)
        gsup     = xin(2)
        wsup     = xin(3)
        tinz1    = xin(4)
        ginz1    = xin(5)
        winz1    = xin(6)
        tinz2    = xin(7)
        ginz2    = xin(8)
        winz2    = xin(9)
        wret     = xin(10)
        tamb     = xin(11)
        gamb     = xin(12)
        wleak    = xin(13)
        tsolairr = xin(14)
        tsolairp = xin(15)
        fracocc  = xin(16)
        fraclig  = xin(17)
        fracpow  = xin(18)
        qsduct   = xin(19)
        qeduct   = xin(20)
        troom    = xin(21)
        tsroom   = xin(22)
        tplen    = xin(23)
        tsplen   = xin(24)
        groom    = xin(25)
! **** Read in parameters
        xcap     = par_v(1)
        rwsr     = par_v(2)
        risr     = par_v(3)
        rosr     = par_v(4)
        rwsp     = par_v(5)
        risp     = par_v(6)
        rosp     = par_v(7)
        rrp      = par_v(8)
        csr      = par_v(9)
        crair    = par_v(10)
        csp      = par_v(11)
        cp       = par_v(12)
        vroom    = par_v(13)
        vplen    = par_v(14)
        noccup   = par_v(15)
        qlite    = par_v(16)
        flpln    = par_v(17)
        qequp    = par_v(18)
        nfile    = nint(par_v(19))
        cr       = crair*xcap
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Initialize room conditions to 20 deg c, 50% rh
                do is = 2,8,2
                    saved_v(is) = 20.0
                enddo
                do is = 10,10,2
                    saved_v(is) = 0.0074
                enddo
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update previous sample instant values
            do is = 2,nsv-1,2
                saved_v(is+1) = saved_v(is)
            enddo
!    Open zone heat gain/loss file - error -> zero extra gain
            call rfile(nfile,'zone',nfp,itype,fpar)
            do i=1,nfp
                saved_v(nsv+i)=fpar(i)
            enddo
        endif
        qartif = saved_v(nsv+1)
! **** Update previous values
        troomp  = saved_v(3)
        tsroomp = saved_v(5)
        tplenp  = saved_v(7)
        tsplenp = saved_v(9)
        groomp  = saved_v(11)
! **** Set up and solve heat and moisture balances
! **** Calculate local extract from mass balance
        wfan = wsup + winz1 - winz2 - wleak - wret
! **** Calculate capacity rates entering room and plenum.
! **** nb wsup and winz1 are positive entering the room, winz2, wleak, wret and
! **** wfan are positive leaving the room. Capacity rates are non-zero if the
! **** flow is into the zone.  All leakage is assumed to be from the zone rather
! **** than the plenum.
        diwsup  = max(0.0,wsup)
        diwinz1 = max(0.0,winz1)
        diwinz2 = max(0.0,-winz2)
        diwleak = max(0.0,-wleak)
        diwret  = max(0.0,-wret)
        diwfan  = max(0.0,-wfan)
        csup    = diwsup*(cpa+gsup*cpg)
        cinz1   = diwinz1*(cpa+ginz1*cpg)
        cinz2   = diwinz2*(cpa+ginz2*cpg)
        cleak   = diwleak*(cpa+gamb*cpg)
        cret    = diwret*(cpa+gamb*cpg)
        cfan    = diwfan*(cpa+gamb*cpg)
        cretp   = wret*(cpa+groom*cpg)
! **** Calculate internal gains
! **** Sensible gains from occupants
        qoccs = fracocc*noccup*qsperson
! **** Latent gains from occupants
        wvapr = (fracocc*noccup*qlperson)/hfg
        if (cret>wsmall) then
! **** Flow through luminaire - divide heat between room and return air
            qrlig = fraclig*(1.-flpln)*qlite
            qelig = fraclig*flpln*qlite
            qplig = 0.0
        else
! **** No flow through luminaire - 50% to room, 50% to plenum
            qrlig = fraclig*0.5*qlite
            qplig = fraclig*0.5*qlite
            qelig = 0.0
        endif
! **** Equipment
        qpow   = fracpow*qequp
! **** Sensible gains to room and plenum
        qsensr = qoccs+qrlig+qpow+qartif
        qsensp = qsduct+qeduct+qplig
! **** Heat balance on room air
        aa = -(1./cr)*(1./risr+1./rwsr+1./rrp+csup+cinz1+cinz2+cleak+&
                       cret+cfan)
        bb =  (1./cr)*(tsroom/risr+tsolairr/rwsr+tplen/rrp+&
                       tsup*csup+tinz1*cinz1+tinz2*cinz2+&
                       tamb*(cleak+cret+cfan)+qsensr)
        call diffeq(time,aa,bb,troomp,troomn,troomb)
! **** Heat balance on room structure
        aa = -(1./csr)*(1./risr+1./rosr)
        bb =  (1./csr)*(troom/risr+tsolairr/rosr)
        call diffeq(time,aa,bb,tsroomp,tsroomn,tsroomb)
! **** Heat balance on plenum air
        aa = -(1./cp)*(1./risp+1./rwsp+1./rrp)
        bb =  (1./cp)*(tsplen/risp+tsolairp/rwsp+troom/rrp+&
                       qsensp)
        call diffeq(time,aa,bb,tplenp,tplenn,tplenb)
! **** Heat balance on plenum structure
        aa = -(1./csp)*(1./risp+1./rosp)
        bb =  (1./csp)*(tplen/risp+tsolairp/rosp)
        call diffeq(time,aa,bb,tsplenp,tsplenn,tsplenb)
! **** Moisture balance on room air
        mroom = vroom*rhoair
        aa = -(1./mroom)*(diwsup+diwinz1+diwinz2+diwleak+diwret+diwfan)
        bb =  (1./mroom)*(gsup*diwsup+ginz1*diwinz1+ginz2*diwinz2+&
                          gamb*(diwleak+diwret+diwfan)+wvapr)
        call diffeq(time,aa,bb,groomp,groomn,groomb)
! **** Extract temperature
        if (wret>wsmall) then
! **** Flow through luminaire - calculate heat pick-up
            tret = troomn+qelig/cretp
        else
! **** No flow - temperature indeterminate
            tret = troomn
        endif
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = troomn
        saved_v(4)  = tsroomn
        saved_v(6)  = tplenn
        saved_v(8)  = tsplenn
        saved_v(10) = groomn
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = troomn
        yout(2) = tsroomn
        yout(3) = tplenn
        yout(4) = tsplenn
        yout(5) = groomn
        yout(6) = tret
        yout(7) = qsensr
        yout(8) = qsensp
        yout(9) = wvapr
! **** Disallow freezing of dynamic variables
        do i=1,5
            iostat(i) = 0
        enddo
! **** Allow freezing of algebraic variables
        do i=6,9
            iostat(i) = 1
        enddo
! **** Return
        return
        end subroutine type403

! ***********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! ***********************************************************************
! * SUBROUTINE:     Room with plenum and ducted return (numerical
! *                 integration)
! *
! * PURPOSE:        Calculate the temperature and humidity ratio in the
! *                 room and in the plenum by performing heat
! *                 balances on the aggregated structure node and the
! *                 air node for each space and moisture balances on the
! *                 air node in each space. Include the effects of air
! *                 flows from adjacent zones and any (inward) air leakage
! *                 and reverse return flow in addition to the HVAC supply
! *                 air flow. Include duct heat gains/losses from/to the
! *                 plenum. Calculate derivatives for numerical integration.
! *
! ***********************************************************************
! * INPUTS
! * ======
! *  1. tsup    : Supply air dry bulb temperature                    (C)
! *  2. gsup    : Supply air humidity ratio                      (kg/kg)
! *  3. wsup    : Supply dry air mass flow rate                   (kg/s)
! *  4. tinz1   : Interzone 1 air dry bulb temperature               (C)
! *  5. ginz1   : Interzone 1 air humidity ratio                 (kg/kg)
! *  6. winz1   : Interzone 1 dry air mass flow rate              (kg/s)
! *  7. tinz2   : Interzone 2 air dry bulb temperature               (C)
! *  8. ginz2   : Interzone 2 air humidity ratio                 (kg/kg)
! *  9. winz2   : Interzone 2 dry air mass flow rate              (kg/s)
! * 10. wret    : Extract dry air mass flow rate                  (kg/s)
! * 11. tamb    : Ambient dry bulb temperature                       (C)
! * 12. gamb    : Ambient humidity ratio                         (kg/kg)
! * 13. wleak   : Leakage air mass flow rate (positive out)      (kg/kg)
! * 14. tsolairr: Equivalent "sol-air" outdoor temperature for room  (C)
! * 15. tsolairp: Equivalent "sol-air" outdoor temperature for plenum(C)
! * 16. fracocc : Fractional occupancy                               (-)
! * 17. fraclig : Fractional lighting heat gain                      (-)
! * 18. fracpow : Fractional equipment heat gain                     (-)
! * 19. qsduct  : Heat gain from supply duct                        (kW)
! * 20. qeduct  : Heat gain from return duct                        (kW)
! *
! * OUTPUTS
! * =======
! *  1. troomn  : Room temperature                                   (C)
! *  2. tsroomn : Room structure temperature                         (C)
! *  3. tplenn  : Plenum temperature                                 (C)
! *  4. tsplenn : Plenum structure temperature                       (C)
! *  5. groomn  : Room humidity ratio                            (kg/kg)
! *  6. tret    : Return temperature                                 (C)
! *  7. qsensr  : Sensible heat gains of room                       (kW)
! *  8. qsensp  : Sensible heat gains of plenum                     (kW)
! *  9. wvapr   : Water vapour gains of room                      (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. xcap    : Room air capacity multiplier                       (-)
! *  2. rwsr    : Direct resistance room air node <-> ambient     (K/kW)
! *  3. risr    : Resistance room air node <-> room mass node     (K/kW)
! *  4. rosr    : Resistance ambient <-> room mass node           (K/kW)
! *  5. rwsp    : Direct resistance plenum air node <-> ambient   (K/kW)
! *  6. risp    : Resistance plenum air node <-> plenum mass node (K/kW)
! *  7. rosp    : Resistance ambient <-> plenum mass node         (K/kW)
! *  8. rr      : Resistance room air node <-> plenum air node    (K/kW)
! *  9. csr     : Capacitance of room mass node                   (kJ/K)
! * 10. crair   : Capacitance of room air node (unmodified)       (kJ/K)
! * 11. csp     : Capacitance of plenum mass node                 (kJ/K)
! * 12. cp      : Capacitance of plenum air node                  (kJ/K)
! * 13. vroom   : Volume of room                                    (m3)
! * 14. vplen   : Volume of plenum                                  (m3)
! * 15. noccup  : Number of occupants                                (-)
! * 16. qlite   : Lighting heat gain                                (kW)
! * 17. flpln   : Fraction of lighting heat gain to return air       (-)
! * 18. qequp   : Equipment heat gain                               (kW)
! * 19. nfile   : Zone number (parameter file='zoneN.par', N > 0)    (-)
! *
! * SAVED
! * =====
! *  1. time    : Time of previous call of TYPE
! *  2. troomn  : Room temperature from previous call
! *  3. troomp  : Room temperature from previous step time
! *  4. tsroomn : Room structure temperature from previous call
! *  5. tsroomp : Room structure temperature from previous step time
! *  6. tplenn  : Plenum temperature from previous call
! *  7. tplenp  : Plenum temperature from previous step time
! *  8. tsplenn : Plenum structure temperature from previous call
! *  9. tsplenp : Plenum structure temperature from previous step time
! * 10. groomn  : Room humidity ratio from previous call
! * 11. groomp  : Room humidity ratio from previous step time
! * 12. qartif  : Heat input for commissioning (read from file)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 29, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  RFILE
!   FUNCTIONS  CALLED:
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! * mroom   : mass of air in room
! * mplen   : mass of air in plenum
! * cpa     : specific heat of dry air
! * cpg     : specific heat of water vapor
! * hfg     : latent heat of water vapor
! * rhoair  : density of air
! * qsperson: sensible heat gain from a person
! * qlperson: latent heat gain from a person
! * wsmall  : threshold for significant air flow rate
! * diwsup  : value of supply air flow rate if flowing in, zero if out
! * diwinz1 : value of air flow rate from Adjacent Zone 1 if flowing in
! * diwinz2 : value of air flow rate from Adjacent Zone 2 if flowing in
! * diwleak : value of leakage air flow rate if flowing in, zero if out
! * csup    : value of supply air capacity rate if flowing in, zero if out
! * cinz1   : value of air capacity rate from Adjacent Zone 1 if flowing in
! * cinz2   : value of air capacity rate from Adjacent Zone 2 if flowing in
! * cleak   : value of leakage air capacity rate if flowing in, zero if out
! * cret    : return air capacity rate
! * qoccs   : sensible gains from occupants
! * wvapr   : latent gains from occupants
! * qrlig   : heat gain to room from lights
! * qelig   : heat gain to return air stream from lights
! * qplig   : heat gain to plenum from lights
! * qpow    : heat gain to room from equipment
! * aa      : A coefficent in dT/dt = A*T + B
! * bb      : B coefficent in dT/dt = A*T + B
! * xxxxxb  : average over time-step value of integrated variable (not used)
! *
! **********************************************************************
!
!   Updated to Fortran 90 April 27, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type404(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndiffeq=5,nsv=(1+ndiffeq*2),&
                                             ni=20,no=9,np=19,nfp=1,&
                                             ns=nsv+nfp
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

        real         :: mroom,noccup
        real         :: cpa=1.006,cpg=1.805,hfg=2501.,rhoair=1.2
        real         :: qsperson=0.095,qlperson=0.045,wsmall=1.e-4
        integer      :: itype=403

        real         :: rrp,csr,crair,csp,cp,vroom,vplen,qlite,flpln,&
                        qequp,cr,qartif,troomp,tsroomp,tplenp,&
                        tsplenp,groomp,wfan,diwsup,diwinz1,diwinz2,&
                        diwleak,diwret,diwfan,csup,cinz1,cinz2,cleak,&
                        cret,cfan,qoccs,wvapr,qrlig,qelig,qplig,qpow,&
                        qsensr,qsensp,dtroomdt,troomn,dtsroomdt,tsroomn,&
                        dtplendt,tplenn,dtsplendt,tsplenn,dgroomdt,&
                        groomn,tret,tsup,gsup,wsup,tinz1,ginz1,winz1,&
                        tinz2,ginz2,winz2,wret,tamb,gamb,wleak,tsolairr,&
                        tsolairp,fracocc,fraclig,fracpow,qsduct,&
                        rwsr,risr,rosr,rwsp,risp,rosp,qeduct,xcap

        integer      :: i,is,nfile

! **** Read in inputs
        tsup     = xin(1)
        gsup     = xin(2)
        wsup     = xin(3)
        tinz1    = xin(4)
        ginz1    = xin(5)
        winz1    = xin(6)
        tinz2    = xin(7)
        ginz2    = xin(8)
        winz2    = xin(9)
        wret     = xin(10)
        tamb     = xin(11)
        gamb     = xin(12)
        wleak    = xin(13)
        tsolairr = xin(14)
        tsolairp = xin(15)
        fracocc  = xin(16)
        fraclig  = xin(17)
        fracpow  = xin(18)
        qsduct   = xin(19)
        qeduct   = xin(20)
! **** Read in parameters
        xcap     = par_v(1)
        rwsr     = par_v(2)
        risr     = par_v(3)
        rosr     = par_v(4)
        rwsp     = par_v(5)
        risp     = par_v(6)
        rosp     = par_v(7)
        rrp      = par_v(8)
        csr      = par_v(9)
        crair    = par_v(10)
        csp      = par_v(11)
        cp       = par_v(12)
        vroom    = par_v(13)
        vplen    = par_v(14)
        noccup   = par_v(15)
        qlite    = par_v(16)
        flpln    = par_v(17)
        qequp    = par_v(18)
        nfile    = nint(par_v(19))
        cr       = crair*xcap
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Initialize room conditions to 20 deg c, 50% rh
                do is = 2,8,2
                    saved_v(is) = 20.0
                enddo
                do is = 10,10,2
                    saved_v(is) = 0.0074
                enddo
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update previous sample instant values
            do is = 2,nsv-1,2
                saved_v(is+1) = saved_v(is)
            enddo
!     open zone heat gain/loss file - error -> zero extra gain
            call rfile(nfile,'zone',nfp,itype,fpar)
            do i=1,nfp
                saved_v(nsv+i)=fpar(i)
            enddo
        endif
        qartif = saved_v(nsv+1)
! **** Update previous values
        troomp  = saved_v(3)
        tsroomp = saved_v(5)
        tplenp  = saved_v(7)
        tsplenp = saved_v(9)
        groomp  = saved_v(11)
! **** Set up and solve heat and moisture balances
! **** Calculate local extract from mass balance
        wfan = wsup + winz1 - winz2 - wleak - wret
! **** Calculate capacity rates entering room and plenum.
! **** nb wsup and winz1 are positive entering the room, winz2, wleak and wret
! **** are positive leaving the room. Capacity rates are non-zero only if the
! **** flow is into the zone. the hvac and local extract are assumed always to
! **** flow out. All leakage is assumed to be from the zone rather than the
! **** plenum.
        diwsup  = max(0.0,wsup)
        diwinz1 = max(0.0,winz1)
        diwinz2 = max(0.0,-winz2)
        diwleak = max(0.0,-wleak)
        diwret  = max(0.0,-wret)
        diwfan  = max(0.0,-wfan)
        csup    = diwsup*(cpa+gsup*cpg)
        cinz1   = diwinz1*(cpa+ginz1*cpg)
        cinz2   = diwinz2*(cpa+ginz2*cpg)
        cleak   = diwleak*(cpa+gamb*cpg)
        cret    = diwret*(cpa+groomp*cpg)
        cfan    = diwfan*(cpa+gamb*cpg)
! **** Calculate internal gains
! **** Sensible gains from occupants
        qoccs = fracocc*noccup*qsperson
! **** Latent gains from occupants
        wvapr = (fracocc*noccup*qlperson)/hfg
        if (wret>wsmall) then
! **** Flow through luminaire - divide heat between room and return air
            qrlig = fraclig*(1.-flpln)*qlite
            qelig = fraclig*flpln*qlite
            qplig = 0.0
        else
! **** No flow through luminaire - 50% to room, 50% to plenum
            qrlig = fraclig*0.5*qlite
            qplig = fraclig*0.5*qlite
            qelig = 0.0
        endif
! **** Equipment
        qpow   = fracpow*qequp
! **** Sensible gains to room and plenum
        qsensr = qoccs+qrlig+qpow+qartif
        qsensp = qsduct+qeduct+qplig
! **** Heat balance on room air
        dtroomdt = (1./cr)*((tsroomp-troomp)/risr+(tplenp-troomp)/rrp+&
                   (tsolairr-troomp)/rwsr+csup*(tsup-troomp)+&
                   cinz1*(tinz1-troomp)+cinz2*(tinz2-troomp)+&
                   (cleak+cret+cfan)*(tamb-troomp)+qsensr)
        troomn = troomp + dtroomdt*tstep
! **** Heat balance on room structure
        dtsroomdt = (1./csr)*((troomp-tsroomp)/risr+&
                              (tsolairr-tsroomp)/rosr)
        tsroomn = tsroomp + dtsroomdt*tstep
! **** Heat balance on plenum air
        dtplendt = (1./cp)*((troomp-tplenp)/rrp+(tsplenp-tplenp)/risp+&
                   (tsolairp-tplenp)/rwsp+qsensp)
        tplenn = tplenp + dtplendt*tstep
! **** Heat balance on plenum structure
        dtsplendt = (1./csp)*((tplenp-tsplenp)/risp+&
                             (tsolairp-tsplenp)/rosp)
        tsplenn = tsplenp + dtsplendt*tstep
! **** Moisture balance on room air
        mroom = vroom*rhoair
        dgroomdt = (1./mroom)*(diwsup*(gsup-groomp)+&
                    diwinz1*(ginz1-groomp)+diwinz2*(ginz2-groomp)+&
                    (diwleak+diwret+diwfan)*(gamb-groomp)+wvapr)
        groomn = groomp + dgroomdt*tstep
! **** Extract temperature
        if (cret>wsmall) then
! **** Flow through luminaire - calculate heat pick-up
            tret = troomn+qelig/cret
        else
! **** No flow - temperature indeterminate
            tret = troomn
        endif
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = troomn
        saved_v(4)  = tsroomn
        saved_v(6)  = tplenn
        saved_v(8)  = tsplenn
        saved_v(10) = groomn
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = troomn
        yout(2) = tsroomn
        yout(3) = tplenn
        yout(4) = tsplenn
        yout(5) = groomn
        yout(6) = tret
        yout(7) = qsensr
        yout(8) = qsensp
        yout(9) = wvapr
! **** Disallow freezing of dynamic variables
        do i=1,5
            iostat(i) = 0
        enddo
! **** Allow freezing of algebraic variables
        do i=6,9
            iostat(i) = 1
        enddo
! **** Return
        return
        end subroutine type404
! ***********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! ***********************************************************************
! * SUBROUTINE:     Room (no plenum, no interzone flows)
! *
! * PURPOSE:        Calculate the temperature and humidity ratio in the
! *                 room by performing heat balances on the aggregated
! *                 structure node and the air node and a moisture
! *                 balance on the air node. Include any (inward) air
! *                 leakage and reverse flow local extract in addition
! *                 to the HVAC supply air flow.
! *
! ***********************************************************************
! * INPUTS
! * ======
! *  1. tsup    : Supply air dry bulb temperature                    (C)
! *  2. gsup    : Supply air humidity ratio                      (kg/kg)
! *  3. wsup    : Supply dry air mass flow rate                   (kg/s)
! *  4. wret    : Extract dry air mass flow rate                  (kg/s)
! *  5. tamb    : Ambient dry bulb temperature                       (C)
! *  6. gamb    : Ambient humidity ratio                         (kg/kg)
! *  7. tsolairr: Equivalent "sol-air" outdoor temperature for room  (C)
! *  8. fracocc : Fractional occupancy                               (-)
! *  9. fraclig : Fractional lighting heat gain                      (-)
! * 10. fracpow : Fractional equipment heat gain                     (-)
! * 11. qsduct  : Heat gain from supply duct                        (kW)
! * 12. qeduct  : Heat gain from return duct                        (kW)
! * 13. troom   : Room temperature                                   (C)
! * 14. tsroom  : Room structure temperature                         (C)
! * 15. groom   : Room humidity ratio                            (kg/kg)
! *
! * OUTPUTS
! * =======
! *  1. troomn  : Room temperature                                   (C)
! *  2. tsroomn : Room structure temperature                         (C)
! *  3. groomn  : Room humidity ratio                            (kg/kg)
! *  4. tret    : Return temperature                                 (C)
! *  5. qsensr  : Sensible heat gains of room                       (kW)
! *  6. qsensp  : Sensible heat gains of plenum                     (kW)
! *  7. wvapr   : Water vapour gains of room                      (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. xcap    : Room air capacity multiplier                       (-)
! *  2. rwsr    : Direct resistance room air node <-> ambient     (K/kW)
! *  3. risr    : Resistance room air node <-> room mass node     (K/kW)
! *  4. rosr    : Resistance ambient <-> room mass node           (K/kW)
! *  5. csr     : Capacitance of room mass node                   (kJ/K)
! *  6. crair   : Capacitance of room air node (unmodified)       (kJ/K)
! *  7. vroom   : Volume of room                                    (m3)
! *  8. noccup  : Number of occupants                                (-)
! *  9. qlite   : Lighting heat gain                                (kW)
! * 10. flpln   : Fraction of lighting heat gain to return air       (-)
! * 11. qequp   : Equipment heat gain                               (kW)
! * 12. nfile   : Zone number (parameter file='ZONEn.PAR', n > 0)    (-)
! *
! * SAVED
! * =====
! *  1. time    : Time of previous call of TYPE
! *  2. troom   : Room temperature from previous call
! *  3. troomp  : Room temperature from previous step time
! *  4. tsroom  : Room structure temperature from previous call
! *  5. tsroomp : Room structure temperature from previous step time
! *  6. groom   : Room humidity ratio from previous call
! *  7. groomp  : Room humidity ratio from previous step time
! *  8. qartif  : Heat input for commissioning (read from file)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 29, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  RFILE, DIFFEQ
!   FUNCTIONS  CALLED:
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! * mroom   : mass of air in room
! * mplen   : mass of air in plenum
! * cpa     : specific heat of dry air
! * cpg     : specific heat of water vapor
! * hfg     : latent heat of water vapor
! * rhoair  : density of air
! * qsperson: sensible heat gain from a person
! * qlperson: latent heat gain from a person
! * wsmall  : threshold for significant air flow rate
! * diwsup  : value of supply air flow rate if flowing in, zero if out
! * diwleak : value of leakage air flow rate if flowing in, zero if out
! * csup    : value of supply air capacity rate if flowing in, zero if out
! * cleak   : value of leakage air capacity rate if flowing in, zero if out
! * cretl   : return air capacity rate
! * qoccs   : sensible gains from occupants
! * wvapr   : latent gains from occupants
! * qrlig   : heat gain to room from lights
! * qelig   : heat gain to return air stream from lights
! * qplig   : heat gain to plenum from lights
! * qpow    : heat gain to room from equipment
! * aa      : A coefficent in dT/dt = A*T + B
! * bb      : B coefficent in dT/dt = A*T + B
! * xxxxxb  : average over time-step value of integrated variable (not used)
! *
! **********************************************************************
!
!   Updated to Fortran 90 April 27, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type411(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndiffeq=3,nsv=(1+ndiffeq*2),&
                                             ni=15,no=7,np=12,nfp=1,&
                                             ns=nsv+nfp
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

        real         :: mroom,noccup
        real         :: cpa=1.006,cpg=1.805,hfg=2501.,rhoair=1.2
        real         :: qsperson=0.095,qlperson=0.045,wsmall=1.e-4
        integer      :: itype=411
        real         :: wsup,wret,tamb,gamb,tsolairr,fracocc,fraclig, &
                        fracpow,qsduct,qeduct,troom,tsroom,groom,xcap,&
                        rwsr,risr,rosr,csr,crair,vroom,qlite,flpln,qequp,&
                        cr,qartif,troomp,tsroomp,groomp,wleak,&
                        diwsup,diwleak,csup,cleak,cretl,qoccs,wvapr,&
                        qrlig,qelig,qpow,qsensr,aa,bb,troomb,troomn,&
                        tsroomb,tsroomn,groomb,groomn,tret,tsup,gsup
        integer      :: i,is,nfile

! **** Read in inputs
        tsup     = xin(1)
        gsup     = xin(2)
        wsup     = xin(3)
        wret     = xin(4)
        tamb     = xin(5)
        gamb     = xin(6)
        tsolairr = xin(7)
        fracocc  = xin(8)
        fraclig  = xin(9)
        fracpow  = xin(10)
        qsduct   = xin(11)
        qeduct   = xin(12)
        troom    = xin(13)
        tsroom   = xin(14)
        groom    = xin(15)
! **** Read in parameters
        xcap     = par_v(1)
        rwsr     = par_v(2)
        risr     = par_v(3)
        rosr     = par_v(4)
        csr      = par_v(5)
        crair    = par_v(6)
        vroom    = par_v(7)
        noccup   = par_v(8)
        qlite    = par_v(9)
        flpln    = par_v(10)
        qequp    = par_v(11)
        nfile    = nint(par_v(12))
        cr       = crair*xcap
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Initialize room conditions to 20 deg c, 50% rh
                do is = 2,4,2
                    saved_v(is) = 20.0
                enddo
                do is = 6,6,2
                    saved_v(is) = 0.0074
                enddo
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update previous sample instant values
            do is = 2,nsv-1,2
                saved_v(is+1) = saved_v(is)
            enddo
!    open zone heat gain/loss file - error -> zero extra gain
            call rfile(nfile,'zone',nfp,itype,fpar)
            do i=1,nfp
                saved_v(nsv+i)=fpar(i)
            enddo
        endif
        qartif = saved_v(nsv+1)
! **** Update previous values
        troomp  = saved_v(3)
        tsroomp = saved_v(5)
        groomp  = saved_v(7)
! **** Set up and solve heat and moisture balances
! **** Calculate leakage from mass balance
        wleak   = wsup - wret
! **** Calculate capacity rates entering room and plenum.
! **** nb wsup is positive entering the room, wleak and wret
! **** are positive leaving the room. Capacity rates are non-zero if the
! **** flow is into the zone. The hvac return is assumed always to flow out.
        diwsup  = max(0.0,wsup)
        diwleak = max(0.0,-wleak)
        csup    = diwsup*(cpa+gsup*cpg)
        cleak   = diwleak*(cpa+gamb*cpg)
        cretl   = wret*(cpa+groom*cpg)
! **** Calculate internal gains
! **** Sensible gains from occupants
        qoccs = fracocc*noccup*qsperson
! **** Latent gains from occupants
        wvapr = (fracocc*noccup*qlperson)/hfg
        if (wret>wsmall) then
! **** Flow through luminaire - divide heat between room and return air
            qrlig = fraclig*(1.-flpln)*qlite
            qelig = fraclig*flpln*qlite
        else
! **** No flow through luminaire - 100% to room
            qrlig = fraclig*qlite
            qelig = 0.0
        endif
! **** Equipment
        qpow   = fracpow*qequp
! **** Sensible gains to room and plenum
        qsensr = qoccs+qrlig+qpow+qartif+qsduct+qeduct
! **** Heat balance on room air
        aa = -(1./cr)*(1./risr+1./rwsr+csup+cleak)
        bb =  (1./cr)*(tsroom/risr+tsolairr/rwsr+&
                       tsup*csup+tamb*cleak+qsensr)
        call diffeq(time,aa,bb,troomp,troomn,troomb)
! **** Heat balance on room structure
        aa = -(1./csr)*(1./risr+1./rosr)
        bb =  (1./csr)*(troom/risr+tsolairr/rosr)
        call diffeq(time,aa,bb,tsroomp,tsroomn,tsroomb)
! **** Moisture balance on room air
        mroom = vroom*rhoair
        aa = -(1./mroom)*(diwsup+diwleak)
        bb =  (1./mroom)*(gsup*diwsup+gamb*diwleak+wvapr)
        call diffeq(time,aa,bb,groomp,groomn,groomb)
! **** Extract temperature
        if (wret>wsmall) then
! **** Flow through luminaire - calculate heat pick-up
            tret = troom+qelig/cretl
        else
! **** No flow - temperature indeterminate
            tret = troom
        endif
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = troomn
        saved_v(4)  = tsroomn
        saved_v(6)  = groomn
! **** Save time of current call
        saved_v(1)  = time
! **** Output
        yout(1) = troomn
        yout(2) = tsroomn
        yout(3) = groomn
        yout(4) = tret
        yout(5) = qsensr
        yout(6) = wvapr
! **** Disallow freezing of dynamic variables
        do i=1,3
            iostat(i) = 0
        enddo
! **** Allow freezing of algebraic variables
        do i=4,6
            iostat(i) = 1
        enddo
! **** Return
        return
        end subroutine type411
! ***********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! ***********************************************************************
! * SUBROUTINE:     Room with plenum return (no interzone flows)
! *
! * PURPOSE:        Calculate the temperature and humidity ratio in the
! *                 room and in the plenum by performing heat
! *                 balances on the aggregated structure node and the
! *                 air node for each space and moisture balances on the
! *                 air node in each space. Include any (inward) air
! *                 leakage and reverse flow local extract in addition
! *                 to the HVAC supply air flow. Include supply duct heat
! *                 gains/losses from/to the plenum.
! *
! ***********************************************************************
! * INPUTS
! * ======
! *  1. tsup    : Supply air dry bulb temperature                    (C)
! *  2. gsup    : Supply air humidity ratio                      (kg/kg)
! *  3. wsup    : Supply dry air mass flow rate                   (kg/s)
! *  4. wret    : Extract dry air mass flow rate                  (kg/s)
! *  5. tamb    : Ambient dry bulb temperature                       (C)
! *  6. gamb    : Ambient humidity ratio                         (kg/kg)
! *  7. tsolairr: Equivalent "sol-air" outdoor temperature for room  (C)
! *  8. tsolairp: Equivalent "sol-air" outdoor temperature for plenum(C)
! *  9. fracocc : Fractional occupancy                               (-)
! * 10. fraclig : Fractional lighting heat gain                      (-)
! * 11. fracpow : Fractional equipment heat gain                     (-)
! * 12. qsduct  : Heat gain from supply duct                        (kW)
! * 13. troom   : Room temperature                                   (C)
! * 14. tsroom  : Room structure temperature                         (C)
! * 15. tplen   : Plenum temperature                                 (C)
! * 16. tsplen  : Plenum structure temperature                       (C)
! * 17. groom   : Room humidity ratio                            (kg/kg)
! * 18. gplen   : Plenum humidity ratio                          (kg/kg)
! *
! * OUTPUTS
! * =======
! *  1. troom   : Room temperature                                   (C)
! *  2. tsroom  : Room structure temperature                         (C)
! *  3. tplen   : Plenum temperature                                 (C)
! *  4. tsplen  : Plenum structure temperature                       (C)
! *  5. groom   : Room humidity ratio                            (kg/kg)
! *  6. gplen   : Plenum humidity ratio                          (kg/kg)
! *  7. tret    : Return temperature                                 (C)
! *  8. qsensr  : Sensible heat gains of room                       (kW)
! *  9. qsensp  : Sensible heat gains of plenum                     (kW)
! * 10. wvapr   : Water vapour gains of room                      (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. xcap    : Room air capacity multiplier                       (-)
! *  2. rwsr    : Direct resistance room air node <-> ambient     (K/kW)
! *  3. risr    : Resistance room air node <-> room mass node     (K/kW)
! *  4. rosr    : Resistance ambient <-> room mass node           (K/kW)
! *  5. rwsp    : Direct resistance plenum air node <-> ambient   (K/kW)
! *  6. risp    : Resistance plenum air node <-> plenum mass node (K/kW)
! *  7. rosp    : Resistance ambient <-> plenum mass node         (K/kW)
! *  8. rr      : Resistance room air node <-> plenum air node    (K/kW)
! *  9. csr     : Capacitance of room mass node                   (kJ/K)
! * 10. crair   : Capacitance of room air node (unmodified)       (kJ/K)
! * 11. csp     : Capacitance of plenum mass node                 (kJ/K)
! * 12. cp      : Capacitance of plenum air node                  (kJ/K)
! * 13. vroom   : Volume of room                                    (m3)
! * 14. vplen   : Volume of plenum                                  (m3)
! * 15. noccup  : Number of occupants                                (-)
! * 16. qlite   : Lighting heat gain                                (kW)
! * 17. flpln   : Fraction of lighting heat gain to extract air      (-)
! * 18. qequp   : Equipment heat gain                               (kW)
! * 19. nfile   : Zone number (parameter file='zoneN.par', n > 0)    (-)
! *
! * SAVED
! * =====
! *  1. time    : Time of previous call of TYPE
! *  2. troom   : Room temperature from previous call
! *  3. troomp  : Room temperature from previous step time
! *  4. tsroom  : Room structure temperature from previous call
! *  5. tsroomp : Room structure temperature from previous step time
! *  6. tplen   : Plenum temperature from previous call
! *  7. tplenp  : Plenum temperature from previous step time
! *  8. tsplen  : Plenum structure temperature from previous call
! *  9. tsplenp : Plenum structure temperature from previous step time
! * 10. groom   : Room humidity ratio from previous call
! * 11. groomp  : Room humidity ratio from previous step time
! * 12. gplen   : Plenum humidity ratio from previous call
! * 13. gplenp  : Plenum humidity ratio from previous step time
! * 14. qartif  : Heat input for commissioning (read from file)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 29, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  RFILE, DIFFEQ
!   FUNCTIONS  CALLED:
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! * mroom   : mass of air in room
! * mplen   : mass of air in plenum
! * cpa     : specific heat of dry air
! * cpg     : specific heat of water vapor
! * hfg     : latent heat of water vapor
! * rhoair  : density of air
! * qsperson: sensible heat gain from a person
! * qlperson: latent heat gain from a person
! * wsmall  : threshold for significant air flow rate
! * wlocal  : leakage plus local extract flow rate
! * diwsup  : value of supply air flow rate if flowing in, zero if out
! * diwlocal: value of leakage+local air flow rate if flowing in, zero if out
! * csup    : value of supply air capacity rate if flowing in, zero if out
! * clocal  : value of leakage+local air cap rate if flowing in, zero if out
! * cret    : extract air capacity rate
! * qoccs   : sensible gains from occupants
! * wvapr   : latent gains from occupants
! * qrlig   : heat gain to room from lights
! * qelig   : heat gain to extract air stream from lights
! * qplig   : heat gain to plenum from lights
! * qpow    : heat gain to room from equipment
! * aa      : A coefficent in dT/dt = A*T + B
! * bb      : B coefficent in dT/dt = A*T + B
! * xxxxxb  : average over time-step value of integrated variable (not used)
! *
! **********************************************************************
!
!   Updated to Fortran 90 April 27, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type412(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndiffeq=6,nsv=(1+ndiffeq*2),&
                                             ni=18,no=10,np=19,nfp=1,&
                                             ns=nsv+nfp
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

        real         :: mroom,noccup
        real         :: cpa=1.006,cpg=1.805,hfg=2501.,rhoair=1.2
        real         :: qsperson=0.095,qlperson=0.045,wsmall=1.e-4
        integer      :: itype=412

        real         :: rosp,rrp,csr,crair,csp,cp,vroom,vplen,qlite,&
                        flpln,qequp,cr,qartif,troomp,tsroomp,&
                        tplenp,tsplenp,groomp,gplenp,wlocal,diwsup,&
                        diwlocal,csup,clocal,cretp,qoccs,wvapr,qrlig,&
                        qelig,qplig,qpow,qsensr,qsensp,aa,bb,troomb,&
                        troomn,tsroomb,tsroomn,tplenb,tplenn,tsplenb,&
                        tsplenn,groomb,groomn,mplen,gplenb,gplenn,tret,&
                        tsup,gsup,wsup,wret,tamb,gamb,tsolair,tsolairr,&
                        fracocc,fraclig,fracpow,qsduct ,troom,tsroom ,&
                        tplen,tsplen ,groom,gplen,xcap,rwsr,risr,rosr,&
                        rwsp,risp,tsolairp
        integer      :: i,is,nfile

! **** Read in inputs
        tsup     = xin(1)
        gsup     = xin(2)
        wsup     = xin(3)
        wret     = xin(4)
        tamb     = xin(5)
        gamb     = xin(6)
        tsolairr = xin(7)
        tsolairp = xin(8)
        fracocc  = xin(9)
        fraclig  = xin(10)
        fracpow  = xin(11)
        qsduct   = xin(12)
        troom    = xin(13)
        tsroom   = xin(14)
        tplen    = xin(15)
        tsplen   = xin(16)
        groom    = xin(17)
        gplen    = xin(18)
! **** Read in parameters
        xcap     = par_v(1)
        rwsr     = par_v(2)
        risr     = par_v(3)
        rosr     = par_v(4)
        rwsp     = par_v(5)
        risp     = par_v(6)
        rosp     = par_v(7)
        rrp      = par_v(8)
        csr      = par_v(9)
        crair    = par_v(10)
        csp      = par_v(11)
        cp       = par_v(12)
        vroom    = par_v(13)
        vplen    = par_v(14)
        noccup   = par_v(15)
        qlite    = par_v(16)
        flpln    = par_v(17)
        qequp    = par_v(18)
        nfile    = nint(par_v(19))
        cr       = crair*xcap
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Initialize room conditions to 20 deg c, 50% rh
                do is = 2,8,2
                    saved_v(is) = 20.0
                enddo
                do is = 10,12,2
                    saved_v(is) = 0.0074
                enddo
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update previous sample instant values
            do is = 2,nsv-1,2
                saved_v(is+1) = saved_v(is)
            enddo
!     open zone heat gain/loss file - error -> zero extra gain
            call rfile(nfile,'zone',nfp,itype,fpar)
            do i=1,nfp
                saved_v(nsv+i)=fpar(i)
            enddo
        endif
        qartif = saved_v(nsv+1)
! **** Update previous values
        troomp  = saved_v(3)
        tsroomp = saved_v(5)
        tplenp  = saved_v(7)
        tsplenp = saved_v(9)
        groomp  = saved_v(11)
        gplenp  = saved_v(13)
! **** Set up and solve heat and moisture balances
! **** Calculate leakage plus local extract flow rate from mass balance
        wlocal = wsup - wret
! **** Calculate capacity rates entering room and plenum.
! **** nb wsup are positive entering the room, wlocal and wret
! **** are positive leaving the room. Capacity rates are non-zero if
! **** the flow is into the zone. All leakage is assumed to be from the zone
! **** rather than the plenum.
        diwsup   = max(0.0,wsup)
        diwlocal = max(0.0,-wlocal)
        csup     = diwsup*(cpa+gsup*cpg)
        clocal   = diwlocal*(cpa+gamb*cpg)
        cretp    = wret*(cpa+groom*cpg)
! **** Calculate internal gains
! **** Sensible gains from occupants
        qoccs = fracocc*noccup*qsperson
! **** Latent gains from occupants
        wvapr = (fracocc*noccup*qlperson)/hfg
        if (wret>wsmall) then
! **** Flow through luminaire - divide heat between room and extract air
            qrlig = fraclig*(1.-flpln)*qlite
            qelig = fraclig*flpln*qlite
            qplig = 0.0
        else
! **** No flow through luminaire - 50% to room, 50% to plenum
            qrlig = fraclig*0.5*qlite
            qplig = fraclig*0.5*qlite
            qelig = 0.0
        endif
! **** Equipment
        qpow   = fracpow*qequp
! **** Sensible gains to room and plenum
        qsensr = qoccs+qrlig+qpow+qartif
        qsensp = qsduct+qplig
! **** Heat balance on room air
        aa = -(1./cr)*(1./risr+1./rwsr+1./rrp+csup+clocal)
        bb =  (1./cr)*(tsroom/risr+tsolairr/rwsr+tplen/rrp+&
                       tsup*csup+tamb*clocal+qsensr)
        call diffeq(time,aa,bb,troomp,troomn,troomb)
! **** Heat balance on room structure
        aa = -(1./csr)*(1./risr+1./rosr)
        bb =  (1./csr)*(troom/risr+tsolairr/rosr)
        call diffeq(time,aa,bb,tsroomp,tsroomn,tsroomb)
! **** Heat balance on plenum air
        aa = -(1./cp)*(1./risp+1./rwsp+1./rrp+cretp)
        bb =  (1./cp)*(tsplen/risp+tsolairp/rwsp+troom/rrp+troom*cretp+&
                       qsensp)
        call diffeq(time,aa,bb,tplenp,tplenn,tplenb)
! **** Heat balance on plenum structure
        aa = -(1./csp)*(1./risp+1./rosp)
        bb =  (1./csp)*(tplen/risp+tsolairp/rosp)
        call diffeq(time,aa,bb,tsplenp,tsplenn,tsplenb)
! **** Moisture balance on room air
        mroom = vroom*rhoair
        aa = -(1./mroom)*(diwsup+diwlocal)
        bb =  (1./mroom)*(gsup*diwsup+gamb*diwlocal+wvapr)
        call diffeq(time,aa,bb,groomp,groomn,groomb)
! **** Moisture balance on plenum air
        mplen = vplen*rhoair
        aa = -(1./mplen)*wret
        bb =  (1./mplen)*groom*wret
        call diffeq(time,aa,bb,gplenp,gplenn,gplenb)
! **** Extract temperature
        if (wret>wsmall) then
! **** Flow through luminaire - calculate heat pick-up
            tret = troom+qelig/cretp
        else
! **** No flow - temperature indeterminate
            tret = troom
        endif
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = troomn
        saved_v(4)  = tsroomn
        saved_v(6)  = tplenn
        saved_v(8)  = tsplenn
        saved_v(10) = groomn
        saved_v(12) = gplenn
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1)  = troomn
        yout(2)  = tsroomn
        yout(3)  = tplenn
        yout(4)  = tsplenn
        yout(5)  = groomn
        yout(6)  = gplenn
        yout(7)  = tret
        yout(8)  = qsensr
        yout(9)  = qsensp
        yout(10) = wvapr
! **** Disallow freezing of dynamic variables
        do i=1,6
            iostat(i) = 0
        enddo
! **** Allow freezing of algebraic variables
        do i=7,10
            iostat(i) = 1
        enddo
! **** Return
        return
        end subroutine type412
! ***********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! ***********************************************************************
! * SUBROUTINE:     Room with plenum and ducted return (analytical
! *                 integration)
! *
! * PURPOSE:        Calculate the temperature and humidity ratio in the
! *                 room and in the plenum by performing heat
! *                 balances on the aggregated structure node and the
! *                 air node for each space and moisture balances on the
! *                 air node in each space. Include any (inward) air
! *                 leakage and reverse flow local return in addition
! *                 to the HVAC supply air flow. Include duct heat
! *                 gains/losses from/to the plenum. Use DIFFEQ to
! *                 perform analytical integration.
! *
! ***********************************************************************
! * INPUTS
! * ======
! *  1. tsup    : Supply air dry bulb temperature                    (C)
! *  2. gsup    : Supply air humidity ratio                      (kg/kg)
! *  3. wsup    : Supply dry air mass flow rate                   (kg/s)
! *  4. wret    : Extract dry air mass flow rate                  (kg/s)
! *  5. tamb    : Ambient dry bulb temperature                       (C)
! *  6. gamb    : Ambient humidity ratio                         (kg/kg)
! *  7. tsolairr: Equivalent "sol-air" outdoor temperature for room  (C)
! *  8. tsolairp: Equivalent "sol-air" outdoor temperature for plenum(C)
! *  9. fracocc : Fractional occupancy                               (-)
! * 10. fraclig : Fractional lighting heat gain                      (-)
! * 11. fracpow : Fractional equipment heat gain                     (-)
! * 12. qsduct  : Heat gain from supply duct                        (kW)
! * 13. qeduct  : Heat gain from return duct                        (kW)
! * 14. troom   : Room temperature                                   (C)
! * 15. tsroom  : Room structure temperature                         (C)
! * 16. tplen   : Plenum temperature                                 (C)
! * 17. tsplen  : Plenum structure temperature                       (C)
! * 18. groom   : Room humidity ratio                            (kg/kg)
! *
! * OUTPUTS
! * =======
! *  1. troomn  : Room temperature                                   (C)
! *  2. tsroomn : Room structure temperature                         (C)
! *  3. tplenn  : Plenum temperature                                 (C)
! *  4. tsplenn : Plenum structure temperature                       (C)
! *  5. groomn  : Room humidity ratio                            (kg/kg)
! *  6. tret    : Return temperature                                 (C)
! *  7. qsensr  : Sensible heat gains of room                       (kW)
! *  8. qsensp  : Sensible heat gains of plenum                     (kW)
! *  9. wvapr   : Water vapour gains of room                      (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. xcap    : Room air capacity multiplier                       (-)
! *  2. rwsr    : Direct resistance room air node <-> ambient     (K/kW)
! *  3. risr    : Resistance room air node <-> room mass node     (K/kW)
! *  4. rosr    : Resistance ambient <-> room mass node           (K/kW)
! *  5. rwsp    : Direct resistance plenum air node <-> ambient   (K/kW)
! *  6. risp    : Resistance plenum air node <-> plenum mass node (K/kW)
! *  7. rosp    : Resistance ambient <-> plenum mass node         (K/kW)
! *  8. rr      : Resistance room air node <-> plenum air node    (K/kW)
! *  9. csr     : Capacitance of room mass node                   (kJ/K)
! * 10. crair   : Capacitance of room air node (unmodified)       (kJ/K)
! * 11. csp     : Capacitance of plenum mass node                 (kJ/K)
! * 12. cp      : Capacitance of plenum air node                  (kJ/K)
! * 13. vroom   : Volume of room                                    (m3)
! * 14. vplen   : Volume of plenum                                  (m3)
! * 15. noccup  : Number of occupants                                (-)
! * 16. qlite   : Lighting heat gain                                (kW)
! * 17. flpln   : Fraction of lighting heat gain to return air       (-)
! * 18. qequp   : Equipment heat gain                               (kW)
! * 19. nfile   : Zone number (parameter file='zoneN.par', N > 0)    (-)
! *
! * SAVED
! * =====
! *  1. time    : Time of previous call of TYPE
! *  2. troomn  : Room temperature from previous call
! *  3. troomp  : Room temperature from previous step time
! *  4. tsroom  : Room structure temperature from previous call
! *  5. tsroomp : Room structure temperature from previous step time
! *  6. tplenn  : Plenum temperature from previous call
! *  7. tplenp  : Plenum temperature from previous step time
! *  8. tsplen  : Plenum structure temperature from previous call
! *  9. tsplenp : Plenum structure temperature from previous step time
! * 10. groom   : Room humidity ratio from previous call
! * 11. groomp  : Room humidity ratio from previous step time
! * 12. qartif  : Heat input for commissioning (read from file)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 29, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  RFILE, DIFFEQ
!   FUNCTIONS  CALLED:
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! * mroom   : mass of air in room
! * mplen   : mass of air in plenum
! * cpa     : specific heat of dry air
! * cpg     : specific heat of water vapor
! * hfg     : latent heat of water vapor
! * rhoair  : density of air
! * qsperson: sensible heat gain from a person
! * qlperson: latent heat gain from a person
! * wsmall  : threshold for significant air flow rate
! * wlocal  : leakage plus local extract flow rate
! * diwsup  : value of supply air flow rate if flowing in, zero if out
! * diwlocal: value of leakage+local air flow rate if flowing in, zero if out
! * csup    : value of supply air capacity rate if flowing in, zero if out
! * clocal  : value of leakage+local air cap rate if flowing in, zero if out
! * cretp   : return air capacity rate
! * qoccs   : sensible gains from occupants
! * wvapr   : latent gains from occupants
! * qrlig   : heat gain to room from lights
! * qelig   : heat gain to return air stream from lights
! * qplig   : heat gain to plenum from lights
! * qpow    : heat gain to room from equipment
! * aa      : A coefficent in dT/dt = A*T + B
! * bb      : B coefficent in dT/dt = A*T + B
! * xxxxxb  : average over time-step value of integrated variable (not used)
! *
! **********************************************************************
!
!   Updated to Fortran 90 April 27, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type413(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndiffeq=5,nsv=(1+ndiffeq*2),&
                                             ni=18,no=9,np=19,nfp=1,&
                                             ns=nsv+nfp
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

        real         :: mroom,noccup
        real         :: cpa=1.006,cpg=1.805,hfg=2501.,rhoair=1.2
        real         :: qsperson=0.095,qlperson=0.045,wsmall=1.e-4
        integer      :: itype=413

        real         :: risr,rosr,rwsp,risp,rosp,rrp,csr,crair,csp,cp,&
                        vroom,vplen,qlite,flpln,qequp,cr,qartif,troomp,&
                        tsroomp,tplenp,tsplenp,groomp,wlocal,diwsup,&
                        diwlocal,csup,clocal,cretp,qoccs,wvapr,qrlig,&
                        qelig,qplig,qpow,qsensr,qsensp,aa,bb,troomb,&
                        troomn,tsroomb,tsroomn,tplenb,tplenn,tsplenb,&
                        tsplenn,groomb,groomn,tret,tsup,gsup,wsup,wret,&
                        tamb,gamb,tsolairr,tsolairp,fracocc,fraclig,&
                        fracpow,qsduct,qeduct,troom,tsroom,tplen,tsplen,&
                        groom,xcap,rwsr
        integer      :: i,is,nfile

! **** Read in inputs
        tsup     = xin(1)
        gsup     = xin(2)
        wsup     = xin(3)
        wret     = xin(4)
        tamb     = xin(5)
        gamb     = xin(6)
        tsolairr = xin(7)
        tsolairp = xin(8)
        fracocc  = xin(9)
        fraclig  = xin(10)
        fracpow  = xin(11)
        qsduct   = xin(12)
        qeduct   = xin(13)
        troom    = xin(14)
        tsroom   = xin(15)
        tplen    = xin(16)
        tsplen   = xin(17)
        groom    = xin(18)
! **** Read in parameters
        xcap     = par_v(1)
        rwsr     = par_v(2)
        risr     = par_v(3)
        rosr     = par_v(4)
        rwsp     = par_v(5)
        risp     = par_v(6)
        rosp     = par_v(7)
        rrp      = par_v(8)
        csr      = par_v(9)
        crair    = par_v(10)
        csp      = par_v(11)
        cp       = par_v(12)
        vroom    = par_v(13)
        vplen    = par_v(14)
        noccup   = par_v(15)
        qlite    = par_v(16)
        flpln    = par_v(17)
        qequp    = par_v(18)
        nfile    = nint(par_v(19))
        cr       = crair*xcap
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Initialize room conditions to 20 deg c, 50% rh
                do is = 2,8,2
                    saved_v(is) = 20.0
                enddo
                do is = 10,10,2
                    saved_v(is) = 0.0074
                enddo
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update previous sample instant values
            do is = 2,nsv-1,2
                saved_v(is+1) = saved_v(is)
            enddo
!     open zone heat gain/loss file - error -> zero extra gain
            call rfile(nfile,'zone',nfp,itype,fpar)
            do i=1,nfp
                saved_v(nsv+i)=fpar(i)
            enddo
        endif
        qartif = saved_v(nsv+1)
! **** Update previous values
        troomp  = saved_v(3)
        tsroomp = saved_v(5)
        tplenp  = saved_v(7)
        tsplenp = saved_v(9)
        groomp  = saved_v(11)
! **** Set up and solve heat and moisture balances
! **** Calculate leakage plus local extract flow rate from mass balance
        wlocal = wsup - wret
! **** Calculate capacity rates entering room and plenum.
! **** nb wsup are positive entering the room, wlocal and wret
! **** are positive leaving the room. Capacity rates are non-zero if
! **** the flow is into the zone. All leakage is assumed to be from the zone
! **** rather than the plenum.
        diwsup   = max(0.0,wsup)
        diwlocal = max(0.0,-wlocal)
        csup     = diwsup*(cpa+gsup*cpg)
        clocal   = diwlocal*(cpa+gamb*cpg)
        cretp    = wret*(cpa+groom*cpg)
! **** Calculate internal gains
! **** Sensible gains from occupants
        qoccs = fracocc*noccup*qsperson
! **** Latent gains from occupants
        wvapr = (fracocc*noccup*qlperson)/hfg
        if (wret>wsmall) then
! **** Flow through luminaire - divide heat between room and return air
            qrlig = fraclig*(1.-flpln)*qlite
            qelig = fraclig*flpln*qlite
            qplig = 0.0
        else
! **** No flow through luminaire - 50% to room, 50% to plenum
            qrlig = fraclig*0.5*qlite
            qplig = fraclig*0.5*qlite
            qelig = 0.0
        endif
! **** Equipment
        qpow   = fracpow*qequp
! **** Sensible gains to room and plenum
        qsensr = qoccs+qrlig+qpow+qartif
        qsensp = qsduct+qeduct+qplig
! **** Heat balance on room air
        aa = -(1./cr)*(1./risr+1./rwsr+1./rrp+csup+clocal)
        bb =  (1./cr)*(tsroom/risr+tsolairr/rwsr+tplen/rrp+&
                       tsup*csup+tamb*clocal+qsensr)
        call diffeq(time,aa,bb,troomp,troomn,troomb)
! **** Heat balance on room structure
        aa = -(1./csr)*(1./risr+1./rosr)
        bb =  (1./csr)*(troom/risr+tsolairr/rosr)
        call diffeq(time,aa,bb,tsroomp,tsroomn,tsroomb)
! **** Heat balance on plenum air
        aa = -(1./cp)*(1./risp+1./rwsp+1./rrp)
        bb =  (1./cp)*(tsplen/risp+tsolairp/rwsp+troom/rrp+&
                       qsensp)
        call diffeq(time,aa,bb,tplenp,tplenn,tplenb)
! **** Heat balance on plenum structure
        aa = -(1./csp)*(1./risp+1./rosp)
        bb =  (1./csp)*(tplen/risp+tsolairp/rosp)
        call diffeq(time,aa,bb,tsplenp,tsplenn,tsplenb)
! **** Moisture balance on room air
        mroom = vroom*rhoair
        aa = -(1./mroom)*(diwsup+diwlocal)
        bb =  (1./mroom)*(gsup*diwsup+gamb*diwlocal+wvapr)
        call diffeq(time,aa,bb,groomp,groomn,groomb)
! **** Extract temperature
        if (wret>wsmall) then
! **** Flow through luminaire - calculate heat pick-up
            tret = troomn+qelig/cretp
        else
! **** No flow - temperature indeterminate
            tret = troomn
        endif
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = troomn
        saved_v(4)  = tsroomn
        saved_v(6)  = tplenn
        saved_v(8)  = tsplenn
        saved_v(10) = groomn
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = troomn
        yout(2) = tsroomn
        yout(3) = tplenn
        yout(4) = tsplenn
        yout(5) = groomn
        yout(6) = tret
        yout(7) = qsensr
        yout(8) = qsensp
        yout(9) = wvapr
! **** Disallow freezing of dynamic variables
        do i=1,5
            iostat(i) = 0
        enddo
! **** Allow freezing of algebraic variables
        do i=6,9
            iostat(i) = 1
        enddo
! **** Return
        return
        end subroutine type413
! ***********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! ***********************************************************************
! * SUBROUTINE:     Room with plenum and ducted return (numerical
! *                 integration)
! *
! * PURPOSE:        Calculate the temperature and humidity ratio in the
! *                 room and in the plenum by performing heat
! *                 balances on the aggregated structure node and the
! *                 air node for each space and moisture balances on the
! *                 air node in each space. Include any (inward) air
! *                 leakage and reverse flow local return in addition
! *                 to the HVAC supply air flow. Include duct heat
! *                 gains/losses from/to the plenum. Calculate
! *                 derivatives and integrate by forward difference.
! *
! ***********************************************************************
! * INPUTS
! * ======
! *  1. tsup    : Supply air dry bulb temperature                    (C)
! *  2. gsup    : Supply air humidity ratio                      (kg/kg)
! *  3. wsup    : Supply dry air mass flow rate                   (kg/s)
! *  4. wret    : Extract dry air mass flow rate                  (kg/s)
! *  5. tamb    : Ambient dry bulb temperature                       (C)
! *  6. gamb    : Ambient humidity ratio                         (kg/kg)
! *  7. tsolairr: Equivalent "sol-air" outdoor temperature for room  (C)
! *  8. tsolairp: Equivalent "sol-air" outdoor temperature for plenum(C)
! *  9. fracocc : Fractional occupancy                               (-)
! * 10. fraclig : Fractional lighting heat gain                      (-)
! * 11. fracpow : Fractional equipment heat gain                     (-)
! * 12. qsduct  : Heat gain from supply duct                        (kW)
! * 13. qeduct  : Heat gain from return duct                        (kW)
! *
! * OUTPUTS
! * =======
! *  1. troomn  : Room temperature                                   (C)
! *  2. tsroomn : Room structure temperature                         (C)
! *  3. tplenn  : Plenum temperature                                 (C)
! *  4. tsplenn : Plenum structure temperature                       (C)
! *  5. groomn  : Room humidity ratio                            (kg/kg)
! *  6. tret    : Return temperature                                 (C)
! *  7. qsensr  : Sensible heat gains of room                       (kW)
! *  8. qsensp  : Sensible heat gains of plenum                     (kW)
! *  9. wvapr   : Water vapour gains of room                      (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. xcap    : Room air capacity multiplier                       (-)
! *  2. rwsr    : Direct resistance room air node <-> ambient     (K/kW)
! *  3. risr    : Resistance room air node <-> room mass node     (K/kW)
! *  4. rosr    : Resistance ambient <-> room mass node           (K/kW)
! *  5. rwsp    : Direct resistance plenum air node <-> ambient   (K/kW)
! *  6. risp    : Resistance plenum air node <-> plenum mass node (K/kW)
! *  7. rosp    : Resistance ambient <-> plenum mass node         (K/kW)
! *  8. rr      : Resistance room air node <-> plenum air node    (K/kW)
! *  9. csr     : Capacitance of room mass node                   (kJ/K)
! * 10. crair   : Capacitance of room air node (unmodified)       (kJ/K)
! * 11. csp     : Capacitance of plenum mass node                 (kJ/K)
! * 12. cp      : Capacitance of plenum air node                  (kJ/K)
! * 13. vroom   : Volume of room                                    (m3)
! * 14. vplen   : Volume of plenum                                  (m3)
! * 15. noccup  : Number of occupants                                (-)
! * 16. qlite   : Lighting heat gain                                (kW)
! * 17. flpln   : Fraction of lighting heat gain to return air       (-)
! * 18. qequp   : Equipment heat gain                               (kW)
! * 19. nfile   : Zone number (parameter file='zoneN.par', N > 0)    (-)
! *
! * SAVED
! * =====
! *  1. time    : Time of previous call of TYPE
! *  2. troomn  : Room temperature from previous call
! *  3. troomp  : Room temperature from previous step time
! *  4. tsroomn : Room structure temperature from previous call
! *  5. tsroomp : Room structure temperature from previous step time
! *  6. tplenn  : Plenum temperature from previous call
! *  7. tplenp  : Plenum temperature from previous step time
! *  8. tsplenn : Plenum structure temperature from previous call
! *  9. tsplenp : Plenum structure temperature from previous step time
! * 10. groomn  : Room humidity ratio from previous call
! * 11. groomp  : Room humidity ratio from previous step time
! * 12. qartif  : Heat input for commissioning (read from file)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 29, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  RFILE
!   FUNCTIONS  CALLED:
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! * mroom   : mass of air in room
! * mplen   : mass of air in plenum
! * cpa     : specific heat of dry air
! * cpg     : specific heat of water vapor
! * hfg     : latent heat of water vapor
! * rhoair  : density of air
! * qsperson: sensible heat gain from a person
! * qlperson: latent heat gain from a person
! * wsmall  : threshold for significant air flow rate
! * wlocal  : leakage plus local extract flow rate
! * diwsup  : value of supply air flow rate if flowing in, zero if out
! * diwlocal: value of leakage+local air flow rate if flowing in, zero if out
! * csup    : value of supply air capacity rate if flowing in, zero if out
! * clocal  : value of leakage+local air cap rate if flowing in, zero if out
! * cret    : return air capacity rate
! * qoccs   : sensible gains from occupants
! * wvapr   : latent gains from occupants
! * qrlig   : heat gain to room from lights
! * qelig   : heat gain to return air stream from lights
! * qplig   : heat gain to plenum from lights
! * qpow    : heat gain to room from equipment
! * aa      : A coefficent in dT/dt = A*T + B
! * bb      : B coefficent in dT/dt = A*T + B
! * xxxxxb  : average over time-step value of integrated variable (not used)
! *
! **********************************************************************
!
!   Updated to Fortran 90 April 27, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type414(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndiffeq=5,nsv=(1+ndiffeq*2),&
                                             ni=13,no=9,np=19,nfp=1,&
                                             ns=nsv+nfp
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

        real         :: mroom,noccup
        real         :: cpa=1.006,cpg=1.805,hfg=2501.,rhoair=1.2
        real         :: qsperson=0.095,qlperson=0.045,wsmall=1.e-4
        integer      :: itype=414

        real         :: rwsr,risr,rosr,rwsp,risp,rosp,rrp,csr,crair,&
                        csp,cp,vroom,vplen,qlite,flpln,qequp,&
                        cr,qartif,troomp,tsroomp,tplenp,tsplenp,&
                        groomp,wlocal,diwsup,diwlocal,csup,clocal,&
                        cret,qoccs,wvapr,qrlig,qelig,qplig,qpow,qsensr,&
                        qsensp,dtroomdt,troomn,dtsroomd,tsroomn,dtplendt,&
                        tplenn,dtsplend,tsplenn,dgroomdt,groomn,tret,&
                        tsup,gsup,wsup,wret,tamb,gamb,tsolairr,tsolairp,&
                        fracocc,fraclig,fracpow,qsduct,qeduct,xcap,&
                        dtsroomdt,dtsplendt
        integer      :: i,is,nfile

! **** Read in inputs
        tsup     = xin(1)
        gsup     = xin(2)
        wsup     = xin(3)
        wret     = xin(4)
        tamb     = xin(5)
        gamb     = xin(6)
        tsolairr = xin(7)
        tsolairp = xin(8)
        fracocc  = xin(9)
        fraclig  = xin(10)
        fracpow  = xin(11)
        qsduct   = xin(12)
        qeduct   = xin(13)
! **** Read in parameters
        xcap     = par_v(1)
        rwsr     = par_v(2)
        risr     = par_v(3)
        rosr     = par_v(4)
        rwsp     = par_v(5)
        risp     = par_v(6)
        rosp     = par_v(7)
        rrp      = par_v(8)
        csr      = par_v(9)
        crair    = par_v(10)
        csp      = par_v(11)
        cp       = par_v(12)
        vroom    = par_v(13)
        vplen    = par_v(14)
        noccup   = par_v(15)
        qlite    = par_v(16)
        flpln    = par_v(17)
        qequp    = par_v(18)
        nfile    = nint(par_v(19))
        cr       = crair*xcap
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = 999999.
            endif
            if (init==0) then
! **** Initialize room conditions to 20 deg c, 50% rh
                do is = 2,8,2
                    saved_v(is) = 20.0
                enddo
                do is = 10,10,2
                    saved_v(is) = 0.0074
                enddo
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update previous sample instant values
            do is = 2,nsv-1,2
                saved_v(is+1) = saved_v(is)
            enddo
!    Open zone heat gain/loss file - error -> zero extra gain
            call rfile(nfile,'zone',nfp,itype,fpar)
            do i=1,nfp
                saved_v(nsv+i)=fpar(i)
            enddo
        endif
        qartif = saved_v(nsv+1)
! **** Update previous values
        troomp  = saved_v(3)
        tsroomp = saved_v(5)
        tplenp  = saved_v(7)
        tsplenp = saved_v(9)
        groomp  = saved_v(11)
! **** Set up and solve heat and moisture balances
! **** Calculate leakage plus local extract flow rate from mass balance
        wlocal = wsup - wret
! **** Calculate capacity rates entering room and plenum.
! **** nb wsup are positive entering the room, wlocal and wret
! **** are positive leaving the room. Capacity rates are non-zero if
! **** the flow is into the zone. All leakage is assumed to be from the zone
! **** rather than the plenum.
        diwsup   = max(0.0,wsup)
        diwlocal = max(0.0,-wlocal)
        csup     = diwsup*(cpa+gsup*cpg)
        clocal   = diwlocal*(cpa+gamb*cpg)
        cret     = wret*(cpa+groomp*cpg)
! **** Calculate internal gains
! **** Sensible gains from occupants
        qoccs = fracocc*noccup*qsperson
! **** Latent gains from occupants
        wvapr = (fracocc*noccup*qlperson)/hfg
        if (wret>wsmall) then
! **** Flow through luminaire - divide heat between room and return air
            qrlig = fraclig*(1.-flpln)*qlite
            qelig = fraclig*flpln*qlite
            qplig = 0.0
        else
! **** No flow through luminaire - 50% to room, 50% to plenum
            qrlig = fraclig*0.5*qlite
            qplig = fraclig*0.5*qlite
            qelig = 0.0
        endif
! **** Equipment
        qpow   = fracpow*qequp
! **** Sensible gains to room and plenum
        qsensr = qoccs+qrlig+qpow+qartif
        qsensp = qsduct+qeduct+qplig
! **** Heat balance on room air
        dtroomdt = (1./cr)*((tsroomp-troomp)/risr+(tplenp-troomp)/rrp+&
                   (tsolairr-troomp)/rwsr+csup*(tsup-troomp)+&
                   clocal*(tamb-troomp)+qsensr)
        troomn = troomp + dtroomdt*tstep
! **** Heat balance on room structure
        dtsroomdt = (1./csr)*((troomp-tsroomp)/risr+&
                              (tsolairr-tsroomp)/rosr)
        tsroomn = tsroomp + dtsroomdt*tstep
! **** Heat balance on plenum air
        dtplendt = (1./cp)*((troomp-tplenp)/rrp+(tsplenp-tplenp)/risp+&
                   (tsolairp-tplenp)/rwsp+qsensp)
        tplenn = tplenp + dtplendt*tstep
! **** Heat balance on plenum structure
        dtsplendt = (1./csp)*((tplenp-tsplenp)/risp+&
                             (tsolairp-tsplenp)/rosp)
        tsplenn = tsplenp + dtsplendt*tstep
! **** Moisture balance on room air
        mroom = vroom*rhoair
        dgroomdt = (1./mroom)*(diwsup*(gsup-groomp)+&
                    diwlocal*(gamb-groomp)+wvapr)
        groomn = groomp + dgroomdt*tstep
! **** Extract temperature
        if (wret>wsmall) then
! **** Flow through luminaire - calculate heat pick-up
            tret = troomn+qelig/cret
        else
! **** No flow - temperature indeterminate
            tret = troomn
        endif
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = troomn
        saved_v(4)  = tsroomn
        saved_v(6)  = tplenn
        saved_v(8)  = tsplenn
        saved_v(10) = groomn
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = troomn
        yout(2) = tsroomn
        yout(3) = tplenn
        yout(4) = tsplenn
        yout(5) = groomn
        yout(6) = tret
        yout(7) = qsensr
        yout(8) = qsensp
        yout(9) = wvapr
! **** Disallow freezing of dynamic variables
        do i=1,5
            iostat(i) = 0
        enddo
! **** Allow freezing of algebraic variables
        do i=6,9
            iostat(i) = 1
        enddo
! **** Return
        return
        end subroutine type414
       
