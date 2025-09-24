! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     Motor-driven actuator
! *
! * PURPOSE:        Calculate the position of a constant speed,
! *                 motor-driven, actuator, accounting for hysteresis,
! *                 crank geometry and range mismatch. Reverse action
! *                 and "stuck" faults are also treated.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. posdem  : demanded position                                  (-)
! *
! * OUTPUTS
! * =======
! *  1. poscoel : position of final control element                  (-)
! *  2. posmot  : motor position                                     (-)
! *  3. tssrev  : number of stop/starts/reversals                    (-)
! *  4. ttrav   : total distance travelled by final control element  (-)
! *
! * PARAMETERS
! * ==========
! *  1. directn : 1=forward, -1=reverse, 0=stuck
! *  2. startpos: starting position (0-1)
! *  3. ttran   : travel time (lim-lim) (s)
! *  4. restart : minimum change in demanded position for movement   (-)
! *  5. hys     : hysteresis (-)
! *  6. crang   : crank travel angle (0 for linear)
! *  7. poscoelu: upper limit of control element range on actuator scale (-)
! *  8. poscoell: lower limit of control element range on actuator scale (-)
! *
! *  SAVED
! *  =====
! *  1. time    : time of previous call
! *  2. posmot  : motor position at previous call
! *  3. posmotp : motor position at previous timestep
! *  4. posdem  : demanded position at previous call
! *  5. posdemp : demanded position at previous timestep
! *  6. modir   : +1=forward, -1=backward, 0=stopped at previous call
! *  7. modirp  : +1=forward, -1=backward, 0=stopped at previous step time
! *  8. tssrev  : total no. of stop/starts/reversals at previous call
! *  9. tssrevp : total no. of stop/starts/reversals at prev step time
! * 10. ttrav   : total distance travelled at previous call
! * 11. ttravp  : total distance travelled at previous step time
! * 12. posact  : position of actuator at previous call
! * 13. posactp : position of actuator at previous step time
! * 14. poscoel : position of final control element at previous call
! * 15. poscoelp: position of final control element at previous step time
! *
! **********************************************************************
! 
!   MAJOR RESTRICTIONS:  Assumes actuator stops when position error is
!                        exactly zero. Assumes motor starts and stops
!                        instanteously. Assumes demanded position has
!                        the same value throughout the time-step, i.e.
!                        the value of the input changes at the beginning
!                        of the time-step and remains valid through the
!                        time-step.
! 
! 
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
! 
!   DATE:                June 16, 1994
!   LAST MODIFIED:       January 18, 1996
! 
!   INCLUDE FILES:
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   HYSTRSIS
!
!   REVISION HISTORY:    None
! 
!   REFERENCE:           Haves, P., Component-Based Modelling of VAV
!                        Systems, Proc. System Simulation in Buildings
!                        '94, Liege, Belgium, December 1994
! 
! **********************************************************************
! * INTERNAL VAIABLES
! * =================
! * quick    : true if actuator moves essentially instantaneously
! * dposmax  : maximum distance actuator can travel in one time-step
! * dposdem  : difference between demanded and current positions at
! *            start of new time-step
! * dposdemp : difference between demanded and current positions at
! *            start of previous time-step
! * modirp   : motor direction at start of previous time-step
! * modire   : motor direction at end of previous time-step
! * modir    : motor direction at start of new time-step
! * poscrank : linear position of crank
! * posact   : position of actuator after linkage slack
! * con      : true if demanded position constant (frozen or fixed boundary
! *            condition) (HVACSIM+ only)
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type321(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=1,no=4,np=8,ns=15
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        logical      :: con
        integer      :: directn
        logical      :: quick
        real         :: dtr=0.01745,small=1.e-6
        real         :: posdem,startpos,ttran,restart,hys,crang,&
                        poscoelu,poscoell,posmotp,posdemp,tssrevp,&
                        ttravp,posactp,poscoelp,posmot,tssrev,dposmax, &
                        dposdemp,dposdem,range,poscrank,posact,poscoel,&
                        span,ttrav
        integer      :: i,is,modirp,modire,modirdem,modir

        quick = .true.      ! added   12/6/1999

! **** Read in inputs
        posdem   = xin(1)
! **** Read in parameters
        directn  = nint(par_v(1))
        startpos = par_v(2)
        if (startpos<0.0 .or. startpos>1.0) then
            stop 'type 321: starting position out of range'
        endif
        ttran    = par_v(3)
        restart  = par_v(4)
        hys      = par_v(5)
        crang    = par_v(6)
        poscoelu = par_v(7)
        poscoell = par_v(8)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.0
            endif
            if (init==0) then
                saved_v(2) = startpos
                saved_v(4) = startpos
                do is = 6,ns-1,2
                    saved_v(is) = 0.0
                enddo
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - Update previous sample instant values
            do is=2,ns-1,2
                saved_v(is+1) = saved_v(is)
            enddo
        endif
! **** Previous values
        posmotp  = saved_v(3)
        posdemp  = saved_v(5)
        modirp   = nint(saved_v(7))
        tssrevp  = saved_v(9)
        ttravp   = saved_v(11)
        posactp  = saved_v(13)
        poscoelp = saved_v(15)
! **** Limit control signal and invert if reverse acting
        posdem=max(0.,min(posdem,1.))
        if (directn==-1) posdem=1.-posdem
! **** Determine curent position based on situation at previous step time
        if (modirp==0) then
! **** Motor off throughout previous time-step
            posmot = posmotp
            modire = 0
            tssrev = tssrevp
        else
! **** Motor on at beginning of previous time-step - Determine if demanded
! **** Position attained
            if (ttran<=(tstep*small)) then
                quick   = .true.
            else
                quick   = .false.
                dposmax = tstep/ttran
            endif
! **** Demanded change in position during previous time-step
            dposdemp = posdemp-posmotp
            if (quick .or. abs(dposdemp)<=dposmax) then
! **** Demanded position attained - Actuator stopped moving at some point
! **** during previous time-step
                posmot = posdemp
                modire = 0.0
                tssrev = tssrevp+1.
            else
! **** Demanded position not attained - Actuator was on continuously
                posmot = posmotp+sign(dposmax,dposdemp)
                modire = nint(sign(1.,dposdemp))
                tssrev = tssrevp
            endif
        endif
! **** Determine response to current control signal
        dposdem  = posdem-posmot
! **** State of motor required to minimise position error
        if (dposdem==0.0) then
            modirdem = 0
        else
            modirdem = nint(sign(1.,dposdem))
        endif
! **** First consider case of motor not in desired state
        if (modire/=modirdem) then
! **** Motor was off or moving in `wrong' direction - Determine if motor starts
            if (abs(dposdem)>restart) then
! **** Demanded and actual positions sufficiently different to start motor
                modir = modirdem
                if (modir/=modire) then
! **** Actuator motion different from that at end of previous timestep
                    tssrev = tssrev+1.
                else
! *****Actuator motion unchanged from previous timestep
                    tssrev = tssrev
                endif
            else
! **** Within deadband - Motor stays off
                modir  = 0
                tssrev = tssrev
            endif
        else
! **** Motor in desired state at end of previous timestep - No change in
! **** state
            modir  = modire
            tssrev = tssrev
        endif
        if (directn==0) then
! **** Stuck/disconnected - Set to initial position
            posmot = startpos
        endif
! **** Non-linearity due to crank, and hysteresis
        if (crang>0.0) then
            range    = 2.*sin(dtr*crang/2.)
            poscrank = 0.5+sin(dtr*crang*(posmot-0.5))/range
        else
            poscrank = posmot
        endif
! **** Hysteresis due to slack in linkage
        call hystrsis(poscrank,posactp,hys,posact)
! **** Range mismatch - linear transformation, then limit to range 0-1
        poscoel = span(posact,poscoell,poscoelu,0.0,1.0)
! **** Distance travelled by control element
        ttrav   = ttravp+abs(poscoel-poscoelp)
! **** Save time of current call
        saved_v(1)=time
! **** Save provisional values to be used in updating at next step time
        saved_v(2)  = posmot
        saved_v(4)  = posdem
        saved_v(6)  = float(modir)
        saved_v(8)  = tssrev
        saved_v(10) = ttrav
        saved_v(12) = posact
        saved_v(14) = poscoel
! **** Output
        yout(1) = poscoel
        yout(2) = posmot
        yout(3) = tssrev
        yout(4) = ttrav
! **** Determine whether to allow freezing
! **** Freezing allowed if position error small and demanded position
! **** constant or if actuator responds instantly
        con = (iostat(1)<-1).or.(iostat(1)==2)
        if ((abs(dposdem)<=max(restart,rtolx).and.con).or.quick) then
            do i=1,no
                iostat(i) = 1
            enddo
        else
            do i=1,no
                iostat(i) = 0
            enddo
        endif
! **** Return
        return
        end subroutine type321

! *********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Damper
! *
! * PURPOSE:    Calculates the pressure drop across a damper from the
! *             mass flow rate. Uses a linear relationship between the
! *             position angle of the blade(s) (THETA) and the logarithm
! *             of the loss coefficient in the range 15 deg < THETA <
! *             55 deg (opposed/single) or 65 deg (parallel). Uses
! *             quadratic interpolation functions in the ranges
! *             0 < THETA < 15 deg and 55/65 deg < THETA < 90 deg.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. w       : Mass flow rate                                  (kg/s)
! *  2. pout    : Outlet pressure                                  (kPa)
! *  3. c       : Damper position (0=closed, 1=open)                 (-)
! *
! * OUTPUTS
! * =======
! *  1. pin     : Inlet pressure                                   (kPa)
! *
! * PARAMETERS
! * ==========
! *  1. ipar    : Damper geometry: opposed/single (0) or parallel (1)(-)
! *  2. ropen   : Open resistance                            (0.001/k.m)
! *  3. farea   : Face area                                         (m2)
! *  4. fleak   : Leakage (fraction of full flow)                    (-)
! *  5. a       : a coefficient in ln(K)=a+b*THETA (-1 for Legg)     (-)
! *  6. b       : b coefficient in ln(K)=a+b*THETA (-1 for Legg)     (-)
! *
! **********************************************************************
! 
!   MAJOR RESTRICTIONS:  Uses fixed limits on range of validity of
!                        logoarithmic relationship
! 
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
! 
!   DATE:                November 8, 1995
! 
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   RDAMPER
! 
!   REVISION HISTORY:    None
! 
!   REFERENCE:           825-RP Final Report
! 
! **********************************************************************
! *
! * INTERNAL VARIABLES
! * ==================
! * dp      : pressure drop across fresh air branch
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type322(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=6,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real,dimension(0:1)   :: alegg=(/-1.51,-1.51/),&
                                 blegg=(/0.105,0.0842/)
        real         :: w,pout,c,ropen,farea,fleak,ax,b,r,rdamper,dp,&
                        dpqudlin,pin
        integer      :: i,ipar

! **** Read in inputs
        w    = xin(1)
        pout = xin(2)
        c    = xin(3)
! **** Read and check parameters
        ipar    = nint(par_v(1))
        ropen   = par_v(2)
        if (ropen<=0.0) stop 'type322: ropen not greater than 0'
        farea   = par_v(3)
        fleak   = par_v(4)
        ax      = par_v(5)
! **** Default is value recommended by legg
        if (ax<0.0) ax = alegg(ipar)
        b       = par_v(6)
! **** Default is value recommended by legg
        if (b<0.0) b = blegg(ipar)
! **** Calculate resistance of damper
        r = rdamper(c,ropen,fleak,farea,ax,b,ipar)
! **** Calculate pressure drop and inlet pressure
        dp = dpqudlin(r,w)
        pin = pout+dp
! **** Assign outputs
        yout(1) = pin
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo
        return
        end subroutine type322
! *********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Damper - calculates flow rate
! *
! * PURPOSE:    Calculates the mass flow rate through a damper from the
! *             presuure drop. Uses a linear relationship between the
! *             position angle of the blade(s) (THETA) and the logarithm
! *             of the loss coefficient in the range 15 deg < THETA <
! *             55 deg (opposed/single) or 65 deg (parallel). Uses
! *             quadratic interpolation functions in the ranges
! *             0 < THETA < 15 deg and 55/65 deg < THETA < 90 deg.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. pin     : Inlet pressure                                   (kPa)
! *  2. pout    : Outlet pressure                                  (kPa)
! *  3. c       : Damper position (0=closed, 1=open)                 (-)
! *
! * OUTPUTS
! * =======
! *  1. w       : Mass flow rate                                  (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. ipar    : Damper geometry: opposed/single (0) or parallel (1)(-)
! *  2. ropen   : Open resistance                            (0.001/k.m)
! *  3. farea   : Face area                                         (m2)
! *  4. fleak   : Leakage (fraction of full flow)                    (-)
! *  5. a       : a coefficient in ln(K)=a+b*THETA (-1 for Legg)     (-)
! *  6. b       : b coefficient in ln(K)=a+b*THETA (-1 for Legg)     (-)
! *
! **********************************************************************
! 
!   MAJOR RESTRICTIONS:  Uses fixed limits on range of validity of
!                        logoarithmic relationship
! 
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
! 
!   DATE:                November 8, 1995
! 
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   RDAMPER
! 
!   REVISION HISTORY:    None
! 
!   REFERENCE:           825-RP Final Report
! 
! **********************************************************************
! *
! * INTERNAL VARIABLES
! * ==================
! * dp      : pressure drop across fresh air branch
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type323(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=6,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real,dimension(0:1)   :: alegg=(/-1.51,-1.51/),&
                                 blegg=(/0.105,0.0842/)
        real         :: pin,pout,c,ropen,farea,fleak,ax,b,r,rdamper,&
                        dp,w,wqudlin
        integer      :: i,ipar

! **** Read in inputs
        pin  = xin(1)
        pout = xin(2)
        c    = xin(3)
! **** Read and check parameters
        ipar    = nint(par_v(1))
        ropen   = par_v(2)
        if (ropen<=0.0) stop 'type323: ropen not greater than 0'
        farea   = par_v(3)
        fleak   = par_v(4)
        ax      = par_v(5)
! **** Default is value recommended by legg
        if (ax<0.0) ax = alegg(ipar)
        b       = par_v(6)
! **** Default is value recommended by legg
        if (b<0.0) b = blegg(ipar)
! **** Calculate resistance of damper
        r = rdamper(c,ropen,fleak,farea,ax,b,ipar)
! **** Calculate pressure drop and inlet pressure
        dp = pin-pout
        w = wqudlin(r,dp)
! **** Assign output
        yout(1) = w
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo

        return
        end subroutine type323
! *********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Mixing box (implicit flow)
! *
! * PURPOSE:    Calculates the mixed air temperature and humidity and
! *             the mixed air and extract air pressures. Uses Legg's
! *             correlations for the resistance of parallel and opposed
! *             blade dampers
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. tout    : Outside air dry bulb temperature                   (C)
! *  2. gout    : Outside air humidity ratio                     (kg/kg)
! *  3. text    : Extract air dry bulb temperature                   (C)
! *  4. gext    : Extract air humidity ratio                     (kg/kg)
! *  5. pout    : Outside air intake pressure                      (kPa)
! *  6. pexh    : Exhaust air outlet pressure                      (kPa)
! *  7. wmix    : Mixed air mass flow rate                        (kg/s)
! *  8. wext    : Extract air mass flow rate                      (kg/s)
! *  9. cout    : Outside air damper position (0=closed, 1=open)     (-)
! * 10. crec    : Recirc air damper position (0=open if PAR(16)=0)   (-)
! * 11. cexh    : Exhaust air damper position (0=closed, 1=open)     (-)
! *
! * OUTPUTS
! * =======
! *  1. tmix    : Mixed air dry bulb temperature                     (C)
! *  2. gmix    : Mixed air humidity ratio                       (kg/kg)
! *  3. pmix    : Mixed air pressure                               (kPa)
! *  4. pext    : Extract air pressure                             (kPa)
! *  5. wout    : Outside air mass flow rate                      (kg/s)
! *  6. wrec    : Recirc air mass flow rate                       (kg/s)
! *  7. wexh    : Exhaust air mass flow rate                      (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. iparout : Outside air damper: opposed (0) or parallel (1)    (-)
! *  2. iparrec : Recirc air damper: opposed (0) or parallel (1)     (-)
! *  3. iparexh : Exhaust air damper: opposed (0) or parallel (1)    (-)
! *  4. ropenout: Open resist. of outside air damper         (0.001/k.m)
! *  5. ropenrec: Open resist. of recirc air damper          (0.001/k.m)
! *  6. ropenexh: Open resist. of exhaust air damper         (0.001/k.m)
! *  7. fareaout: Face area of outside air damper                   (m2)
! *  8. farearec: Face area of recirc air damper                    (m2)
! *  9. fareaexh: Face area of exhaust air damper                   (m2)
! * 10. fleakout: Leakage for outside air damper (fraction of full flow)
! * 11. fleakrec: Leakage for recirc air damper  (fraction of full flow)
! * 12. fleakexh: Leakage for exhaust air damper (fraction of full flow)
! * 13. rfixout : Fixed resistance in outside air branch     (0.001/k.m)
! * 14. rfixrec : Fixed resistance in recirc air branch      (0.001/k.m)
! * 15. rfixexh : Fixed resistance in exhaust air branch     (0.001/k.m)
! * 16. noninver: 0=invert recirc air damper, 1=not inverted         (-)
! *
! **********************************************************************
! 
!   MAJOR RESTRICTIONS:  Assumes positive mixed and extract flow rates
! 
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
! 
!   DATE:                November 8, 1995
! 
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  FLOWSPLT, MOISTMIX
!   FUNCTIONS  CALLED:   None
! 
!   REVISION HISTORY:    None
! 
!   REFERENCE:
! 
! **********************************************************************
! *
! * INTERNAL VARIABLES
! * ==================
! * rout    : total resistance of outside air branch
! * rrec    : total resistance of recirculation branch
! * rexh    : total resistance of exhaust branch
! * dpexhout: difference between exhaust and outside pressures
! * dpout   : pressure drop across outside branch
! * dprec   : pressure drop across recirculation branch
! * dpexh   : pressure drop across exhaust branch
! * srout   : "signed" resistance of outside branch
! * srrec   : "signed" resistance of recirculation branch
! * srexh   : "signed" resistance of exhaust branch
! * diff    : (b^2 - 4*a*c)/(2*a)
! * wxxxp   : flow in xxx branch corresponding to positive root
! * wxxxn   : flow in xxx branch corresponding to negative root
! * rootpos : true if positive root consistent with known directions
! * rootneg : true if negative root consistent with known directions
! * wsmall  : threshold for significant flow rate
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type324(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=11,no=7,np=16,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real,dimension(20)    :: t,gx,w
        real,dimension(0:1)   :: alegg=(/-1.51,-1.51/),&
                                 blegg=(/0.105,0.0842/)
        logical      :: rootpos,rootneg
        real         :: wsmall=1.0e-6
        real         :: ropenout,ropenrec,ropenexh,fareaout,farearec,&
                        fareaexh,fleakout,fleakrec,fleakexh,rfixout,&
                        rfixrec,rfixexh,rout,rrec,rexh,dpexhout,&
                        dpout,srout,dprec,srrec,dpexh,srexh,wout,wrec,&
                        wexh,aa,bb,cc,bbb,ccc,bbbsq,diff,z,wrecp,wrecn,&
                        wexhp,woutp,wexhn,woutn,pmix,pext,tmix,gmix,&
                        wexhrev,wdummy,grec,trec,tout,gout,text,gext,&
                        pout,pexh,wmix,wext,cout,crec,cexh,rdamper
        integer      :: i,iparrec,iparexh,iparout,noninver

! **** Read and check parameters
        iparout = nint(par_v(1))
        iparrec = nint(par_v(2))
        iparexh = nint(par_v(3))
        ropenout= par_v(4)
        ropenrec= par_v(5)
        ropenexh= par_v(6)
        fareaout= par_v(7)
        farearec= par_v(8)
        fareaexh= par_v(9)
        fleakout= par_v(10)
        fleakrec= par_v(11)
        fleakexh= par_v(12)
        rfixout= par_v(13)
        rfixrec= par_v(14)
        rfixexh= par_v(15)
        noninver = nint(par_v(16))

! **** Read in inputs
        tout = xin(1)                 !  Moved on Feb. 6, 2007
        gout = xin(2)
        text = xin(3)
        gext = xin(4)
        pout = xin(5)
        pexh = xin(6)
        wmix = xin(7)
        wext = xin(8)
        cout = xin(9)
        if (noninver/=1) then
            crec = 1.0-xin(10)
        else
            crec = xin(10)
        endif
        cexh = xin(11)

! **** Calculate resistance of each branch
! **** Resistance is sum of damper resistance and fixed resistance
        rout = rdamper(cout,ropenout,fleakout,fareaout,alegg(iparout),&
                       blegg(iparout),iparout) + rfixout
        if (rout<=0.0) stop &
            'type324: resistance in outside branch not greater than 0'
        rrec = rdamper(crec,ropenrec,fleakrec,farearec,alegg(iparrec),&
                       blegg(iparrec),iparrec) + rfixrec
        if (rrec<=0.0) stop &
            'type324: resistance in recirc branch not greater than 0'
        rexh = rdamper(cexh,ropenexh,fleakexh,fareaexh,alegg(iparexh),&
                       blegg(iparexh),iparexh) + rfixexh
        if (rexh<=0.0) stop &
            'type324: resistance in exhaust branch not greater than 0'
! **** Calculate "signed" resistances (see printed documentation)
! **** Test for direction of the outside, return and exhaust flows by
! **** determining the sign of the pressure difference across the
! **** corresponding resistance assuming the flow in question is zero.
        dpexhout = pexh - pout
! **** Outside:
        dpout = rrec*wmix*abs(wmix) - rexh*(wext-wmix)*abs(wext-wmix) -&
                dpexhout
        srout = sign(rout,dpout)
! **** Return:
        dprec = rexh*wext*abs(wext) + rout*wmix*abs(wmix) + dpexhout
        srrec = sign(rrec,dprec)
! **** Exhaust:
        dpexh = rrec*wext*abs(wext) - rout*(wmix-wext)*abs(wmix-wext) -&
                dpexhout
        srexh = sign(rexh,dpexh)
! **** Determine necessity of solving quadratic - Check for special cases
        if (srout==0) then
! **** Outside flow zero
            wout = 0.0
            wrec = wmix
            wexh = wext - wrec
        elseif (srrec==0) then
! **** Recirculation flow zero
            wrec = 0.0
            wout = wmix
            wexh = wext
        elseif (srexh==0) then
! **** Exhaust flow zero
            wexh = 0.0
            wrec = wext
            wout = wmix -wrec
        else
! **** All three unknown flows non-zero:
! **** Calc wrec by solving a quadratic of the form a*wrec**2 + b*wrec + c = 0
! **** Note that aa=a, bb=-b/2, cc=c,  bbb=-b/2a,  ccc=c/a
            aa   = srout+srexh-srrec
            bb   = srout*wmix+srexh*wext
            cc   = srout*wmix*wmix+srexh*wext*wext+dpexhout
            if (aa==0.) then
! **** One root
                wrec = cc/(2.0*bb)
            else
! **** Two roots
                bbb = bb/aa
                ccc = cc/aa
                bbbsq = bbb*bbb
                diff = bbbsq-ccc
                if (diff>0) then
! **** Roots unequal
                    z = sqrt(diff)
                    wrecp = bbb + z
                    wrecn = bbb - z
! **** Check each root for consistency with directions of air flow
! **** already determined
                    wexhp = wext - wrecp
                    woutp = wmix - wrecp
                    if (wexhp*srexh<0.0 .or.&
                        wrecp*srrec<0.0 .or.&
                        woutp*srout<0.0) then
                        rootpos = .false.
                    else
                        rootpos = .true.
                    endif
                    wexhn = wext-wrecn
                    woutn = wmix-wrecn
                    if (wexhn*srexh<0.0 .or.&
                        wrecn*srrec<0.0 .or.&
                        woutn*srout<0.0) then
                        rootneg = .false.
                    else
                        rootneg = .true.
                    endif
                    if (.not.rootpos .and. .not.rootneg) then
                        stop 'type324: no consistent solution'
                    elseif (rootpos .and. rootneg) then
                           stop 'type324: ambiguous solution'
                       elseif (rootpos .and. .not.rootneg) then
                           wrec = wrecp
                       elseif (rootneg .and. .not.rootpos) then
                           wrec = wrecn
                    else
                        stop 'type324: logic fault 1'
                    endif
                elseif (diff==0) then
! **** Equal roots
                    wrec = bbb
                else
! **** Complex roots (coding error, negative resistances ...?)
                    stop  'type 324: complex roots'
                endif
            endif
! **** Calculate remaining unknown flows
            wout = wmix - wrec
            wexh = wext - wrec
! **** Calculate mixed and extract pressures
            pmix = pout - srout*wout*wout
            pext = pexh + srexh*wexh*wexh
        endif
! **** Having calculated flows, calculate mixed temperature and humidity
! **** Model assumes wmix and wext both non-negative. warn if not true.
        if (lastit==1) then
            if (wmix<-wsmall) write(*,*) 'type 324: wmix < ',-wsmall
            if (wext<-wsmall) write(*,*) 'type 324: wext < ',-wsmall
        endif
! **** Check for zero or reverse mixed flow
        if (wmix<=wsmall) then
! **** Reverse mixed flow, take mixed temperature and humidity ratio to be
! **** an average of the outside and extract conditions so as to provide a
! **** reasonable guess for subsequent iterations
            tmix = (tout+text)/2.0
            gmix = (gout+gext)/2.0
        elseif (wmix>wsmall .and. wext<=wsmall) then
! **** Zero/reverse extract flow - Mixed condition is the same as outside
! **** condition
            tmix = tout
            gmix = gout
        elseif (wmix>wsmall .and. wext>wsmall ) then
! **** Normal, forward mixed and extract flows - Calculate mixed conditions
! **** first calculate recirc air conditions if recirc flow positive
            if (wrec>wsmall) then
! **** Test for zero/reverse exhaust flow
                if (wexh<wsmall) then
! **** Reverse exhaust flow - Recirc flow is a mixture of extract and outside
! **** conditions
                    wexhrev = -wexh
                    t(1) = text
                    gx(1) = gext
                    w(1) = wext
                    t(2) = tout
                    gx(2) = gout
                    w(2) = wexhrev
                    call moistmix(t,gx,w,2,trec,grec,wdummy)
                else
! **** Normal, forward exhaust flow - recirc condition is extract condition
                    trec = text
                    grec = gext
                endif
            endif
! **** Calculate mixed air conditions
! **** mixed condition is a combination of the recirc and outside conditions
! **** moistmix also tests for reverse flows and sets the mixed condition
! **** to the forward flow condition, or to an average of the flow conditions
! **** if both are reverse
            t(1) = trec
            gx(1) = grec
            w(1) = wrec
            t(2) = tout
            gx(2) = gout
            w(2) = wout
            call moistmix(t,gx,w,2,tmix,gmix,wdummy)
        else
            stop 'type324: logic fault 2'
        endif
! **** Assign outputs
        yout(1) = tmix
        yout(2) = gmix
        yout(3) = pmix
        yout(4) = pext
        yout(5) = wout
        yout(6) = wrec
        yout(7) = wexh
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo

        return
        end subroutine type324

! *********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Mixing box
! *
! * PURPOSE:    Calculates the mixed air flow rate, temperature and
! *             humidity and the extract air pressure. Uses Legg's
! *             correlations for the resistance of parallel and opposed
! *             blade dampers
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. tri     : Outside air dry bulb temperature                   (C)
! *  2. gout    : Outside air humidity ratio                     (kg/kg)
! *  3. text    : Extract air dry bulb temperature                   (C)
! *  4. gext    : Extract air humidity ratio                     (kg/kg)
! *  5. pout    : Outside air intake pressure                      (kPa)
! *  6. pexh    : Exhaust air outlet pressure                      (kPa)
! *  7. pmix    : Mixed air outlet pressure                        (kPa)
! *  8. wext    : Extract air mass flow rate                      (kg/s)
! *  9. cout    : Outside air damper position (0=closed, 1=open)     (-)
! * 10. crec    : Recirc air damper position (0=open if PAR(16)=0)   (-)
! * 11. cexh    : Extract air damper position (0=closed, 1=open)     (-)
! *
! * OUTPUTS
! * =======
! *  1. tmix    : Mixed air dry bulb temperature                     (C)
! *  2. gmix    : Mixed air humidity ratio                       (kg/kg)
! *  3. wmix    : Mixed air mass flow rate                        (kg/s)
! *  4. pext    : Extract air pressure                             (kPa)
! *  5. wout    : Outside air mass flow rate                      (kg/s)
! *  6. wrec    : Recirc air mass flow rate                       (kg/s)
! *  7. wexh    : Exhaust air mass flow rate                      (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. iparout : Outside air damper: opposed (0) or parallel (1)    (-)
! *  2. iparrec : Recirc air damper: opposed (0) or parallel (1)     (-)
! *  3. iparexh : Exhaust air damper: opposed (0) or parallel (1)    (-)
! *  4. ropenout: Open resist. of outside air damper         (0.001/k.m)
! *  5. ropenrec: Open resist. of recirc air damper          (0.001/k.m)
! *  6. ropenexh: Open resist. of exhaust air damper         (0.001/k.m)
! *  7. fareaout: Face area of outside air damper                   (m2)
! *  8. farearec: Face area of recirc air damper                    (m2)
! *  9. fareaexh: Face area of exhaust air damper                   (m2)
! * 10. fleakout: Leakage for outside air damper  fraction of full flow)
! * 11. fleakrec: Leakage for recirc air damper  (fraction of full flow)
! * 12. fleakexh: Leakage for exhaust air damper (fraction of full flow)
! * 13. rfixout : Fixed resistance in outside air branch     (0.001/k.m)
! * 14. rfixrec : Fixed resistance in recirc air branch      (0.001/k.m)
! * 15. rfixexh : Fixed resistance in exhaust air branch     (0.001/k.m)
! * 16. noninver: 0=invert recirc air damper, 1=not inverted         (-)
! *
! **********************************************************************
! 
!   MAJOR RESTRICTIONS:  Assumes positive mixed and extract flow rates
! 
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
! 
!   DATE:                November 8, 1995
! 
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  FLOWSPLT, MOISTMIX
!   FUNCTIONS  CALLED:   None
! 
!   REVISION HISTORY:    None
! 
!   REFERENCE:
! 
! **********************************************************************
! *
! * INTERNAL VARIABLES
! * ==================
! * rout    : total resistance of outside air branch
! * rrec    : total resistance of recirculation branch
! * rexh    : total resistance of exhaust branch
! * dp      : pressure drop across outside air branch
! * small   : outsideold for significant flow rate
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type325(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=11,no=7,np=16,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real,dimension(20)    :: t,gx,w
        real,dimension(0:1)   :: alegg=(/-1.51,-1.51/),&
                                 blegg=(/0.105,0.0842/)
        real         :: small=1.0e-10,wcritl=3.e-2,wcritu=6.e-2
        real         :: ropenout,ropenrec,ropenexh,fareaout,farearec,&
                        fareaexh,fleakout,fleakrec,fleakexh,rfixout,&
                        rfixrec ,rfixexh,tout,gout,text,gext,pout,&
                        pexh,pmix,wext,cout,crec,cexh,rout,rdamper,rrec,&
                        rexh,wexh,wrec,pext,dp,wout,wqudlin,wmix,tmix,&
                        gmix,wexhrev,wdummy,grec,trec
        integer      :: i,iparout,iparrec,iparexh,ifail,noninver

! **** Read and check parameters
        iparout = nint(par_v(1))
        iparrec = nint(par_v(2))
        iparexh = nint(par_v(3))
        ropenout= par_v(4)
        if (ropenout<=0.0) stop 'type325: ropenout not greater than 0'
        ropenrec= par_v(5)
        if (ropenrec<=0.0) stop 'type325: ropenrec not greater than 0'
        ropenexh= par_v(6)
        if (ropenexh<=0.0) stop 'type325: ropenexh not greater than 0'
        fareaout= par_v(7)
        farearec= par_v(8)
        fareaexh= par_v(9)
        fleakout= par_v(10)
        fleakrec= par_v(11)
        fleakexh= par_v(12)
        rfixout= par_v(13)
        if (rfixout<=0.0) stop 'type325: rfixout not greater than 0'
        rfixrec= par_v(14)
        if (rfixrec<=0.0) stop 'type325: rfixrec not greater than 0'
        rfixexh= par_v(15)
        if (rfixexh<=0.0) stop 'type325: rfixexh not greater than 0'
        noninver = nint(par_v(16))

! **** Read in inputs
        tout = xin(1)                ! Moved on Feb. 6, 2007
        gout = xin(2)
        text = xin(3)
        gext = xin(4)
        pout = xin(5)
        pexh = xin(6)
        pmix = xin(7)
        wext = xin(8)
        cout = xin(9)
        if (noninver/=1) then
            crec = 1.0-xin(10)
        else
            crec = xin(10)
        endif
        cexh = xin(11)
! **** Calculate resistance of each branch
! **** Resistance is sum of damper resistance and fixed resistance
        rout = rdamper(cout,ropenout,fleakout,fareaout,alegg(iparout),&
                       blegg(iparout),iparout) + rfixout
        rrec = rdamper(crec,ropenrec,fleakrec,farearec,alegg(iparrec),&
                       blegg(iparrec),iparrec) + rfixrec
        rexh = rdamper(cexh,ropenexh,fleakexh,fareaexh,alegg(iparexh),&
                       blegg(iparexh),iparexh) + rfixexh
! **** Calculate recirc and exhaust flow rates
! **** Call flow split routine
        call flowsplt(wext,pmix,pexh,0.0,rrec,rexh,wcritl,wcritu,rtolx,&
                      pext,wrec,wexh,ifail)
! **** Check for unsuccesful completion
        if (ifail==1) then
! **** One or more resistances negative
            stop 'type 325: resistances must not be negative'
        elseif (ifail==2) then
! **** Zero resistance for both outlet branches
            stop 'type 325: either rrec or rexh must be non-zero'
        endif
! **** Calculate outside air mass flow rate
        dp=pout-pmix
        wout=wqudlin(rout,dp)
! **** Calculate mixed air flow rate as sum of outside and recirc flow rates
        wmix=wout+wrec
! **** Having calculated flows, calculate mixed temperature and humidity
! **** Model assumes wmix and wext both non-negative. Warn if not true.
        if (wmix<-1.0e-4 .and. lastit==1)&
            write(*,*) 'type 325: wmix < -1.0e-4'
        if (wext<-1.0e-4 .and. lastit==1)&
            write(*,*) 'type 325: wext < -1.0e-4'
! **** Check for reverse mixed flow
        if (wmix<=small) then
! **** Reverse mixed flow, take mixed temperature and humidity ratio to be
! **** an average of the outside and extract conditions so as to provide a
! **** reasonable guess for subsequent iterations
            tmix = (tout+text)/2.0
            gmix = (gout+gext)/2.0
        elseif (wmix>small .and. wext<=small) then
! **** Reverse extract flow - Mixed condition is the same as outside condition
            tmix = tout
            gmix = gout
        elseif (wmix>small .and. wext>small ) then
! **** Normal, forward mixed and extract flows - Calculate mixed conditions
! **** first calculate recirc air conditions if recirc flow positive
            if (wrec>small) then
! **** Test for reverse exhaust flow
                if (wexh<small) then
! **** Reverse exhaust flow - Recirc flow is a mixture of extract and outside
! **** conditions
                    wexhrev=-wexh
                    t(1) = text
                    gx(1) = gext
                    w(1) = wext
                    t(2) = tout
                    gx(2) = gout
                    w(2) = wexhrev
                    call moistmix(t,gx,w,2,trec,grec,wdummy)
                else
! **** Normal, forward exhaust flow - Recirc condition is extract condition
                    trec = text
                    grec = gext
                endif
            endif
! **** Calculate mixed air conditions
! **** Mixed condition is a combination of the recirc and outside conditions
! **** moistmix also tests for reverse flows and sets the mixed condition
! **** to the forward flow condition, or to an average of the flow conditions
! **** if both are reverse
            t(1) = trec
            gx(1) = grec
            w(1) = wrec
            t(2) = tout
            gx(2) = gout
            w(2) = wout
            call moistmix(t,gx,w,2,tmix,gmix,wdummy)
        else
            write(*,*) 'itime = ', itime
            write(*,*) 'type 325: wmix < -1.0e-4'
            write(*,*) 'type 325: wext < -1.0e-4'
            stop 'type325: logic fault'
        endif
! **** Assign outputs
        yout(1) = tmix
        yout(2) = gmix
        yout(3) = wmix
        yout(4) = pext
        yout(5) = wout
        yout(6) = wrec
        yout(7) = wexh
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo

        return
        end subroutine type325

! *********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Mixing box with minimum OA damper - implicit flow
! *
! * PURPOSE:    Calculates the mixed air temperature and humidity and
! *             the mixed air and extract air pressures. Uses Legg's
! *             correlations for the resistance of parallel and opposed
! *             blade dampers. Includes an additional damper in parallel
! *             with the main outside air damper.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. tout    : Outside air dry bulb temperature                   (C)
! *  2. gout    : Outside air humidity ratio                     (kg/kg)
! *  3. text    : Extract air dry bulb temperature                   (C)
! *  4. gext    : Extract air humidity ratio                     (kg/kg)
! *  5. pout    : Outside air intake pressure                      (kPa)
! *  6. pexh    : Exhaust air outlet pressure                      (kPa)
! *  7. wmix    : Mixed air mass flow rate                        (kg/s)
! *  8. wext    : Extract air mass flow rate                      (kg/s)
! *  9. cout    : Outside air damper position (0=closed, 1=open)     (-)
! * 10. crec    : Recirc air damper position (0=open if PAR(16)=0)   (-)
! * 11. cexh    : Extract air damper position (0=closed, 1=open)     (-)
! * 12. cmou    : Minimum outside air damper pos. (0=closed, 1=open) (-)
! *
! * OUTPUTS
! * =======
! *  1. tmix    : Mixed air dry bulb temperature                     (C)
! *  2. gmix    : Mixed air humidity ratio                       (kg/kg)
! *  3. pmix    : Mixed air pressure                               (kPa)
! *  4. pext    : Extract air pressure                             (kPa)
! *  5. wout    : Outside air mass flow rate                      (kg/s)
! *  6. wrec    : Recirc air mass flow rate                       (kg/s)
! *  7. wexh    : Exhaust air mass flow rate                      (kg/s)
! *  8. wmou    : Mass flow rate thru minimum outside air damper  (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. iparout : Outside air damper: opposed (0) or parallel (1)    (-)
! *  2. iparrec : Recirc air damper: opposed (0) or parallel (1)     (-)
! *  3. iparexh : Exhaust air damper: opposed (0) or parallel (1)    (-)
! *  4. iparmou : Min outside air damper: opposed (0) or parallel (1)(-)
! *  5. ropenout: Open resist. of outside air damper         (0.001/k.m)
! *  6. ropenrec: Open resist. of recirc air damper          (0.001/k.m)
! *  7. ropenexh: Open resist. of exhaust air damper         (0.001/k.m)
! *  8. ropenmou: Open resist. of minimum outside air damper (0.001/k.m)
! *  9. fareaout: Face area of outside air damper                   (m2)
! * 10. farearec: Face area of recirc air damper                    (m2)
! * 11. fareaexh: Face area of exhaust air damper                   (m2)
! * 12. fareamou: Face area of minimum outside air damper           (m2)
! * 13. fleakout: Leakage for outside air damper (fraction of full flow)
! * 14. fleakrec: Leakage for recirc air damper  (fraction of full flow)
! * 15. fleakexh: Leakage for exhaust air damper (fraction of full flow)
! * 16. fleakmou: Leakage for minimum outside air damper      (fraction)
! * 17. rfixout : Fixed resistance in outside air branch     (0.001/k.m)
! * 18. rfixrec : Fixed resistance in recirc air branch      (0.001/k.m)
! * 19. rfixexh : Fixed resistance in exhaust air branch     (0.001/k.m)
! * 20. noninver: 0=invert recirc air damper, 1=not inverted         (-)
! *
! **********************************************************************
! 
!   MAJOR RESTRICTIONS:  Assumes positive mixed and extract flow rates
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
! *
! * INTERNAL VARIABLES
! * ==================
! * routd   : resistance of main outside air damper
! * rmoud   : resistance of minimum outside air damper
! * routdp  : resistance of main and minimum dampers in parallel
! * rout    : total resistance of outside air branch
! * rrec    : total resistance of recirculation branch
! * rexh    : total resistance of exhaust branch
! * dpexhout: difference between exhaust and outside pressures
! * dpout   : pressure drop across outside branch
! * dprec   : pressure drop across recirculation branch
! * dpexh   : pressure drop across exhaust branch
! * srout   : "signed" resistance of outside branch
! * srrec   : "signed" resistance of recirculation branch
! * srexh   : "signed" resistance of exhaust branch
! * diff    : (b^2 - 4*a*c)/(2*a)
! * wxxxp   : flow in xxx branch corresponding to positive root
! * wxxxn   : flow in xxx branch corresponding to negative root
! * rootpos : true if positive root consistent with known directions
! * rootneg : true if negative root consistent with known directions
! * wsmall  : threshold for significant flow rate
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type326(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=12,no=8,np=20,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real,dimension(20)    :: t,gx,w
        real,dimension(0:1)   :: alegg=(/-1.51,-1.51/),&
                                 blegg=(/0.105,0.0842/)
        real         :: wsmall=1.0e-6
        logical      :: rootpos,rootneg
        real         :: fleakexh,fleakmou,rfixout ,rfixrec ,rfixexh ,&
                        crec,routd,rmoud,sqrroutd,sqrrmoud,routdp,&
                        rout,rrec,rexh,dpexhout,dpout,srout,dprec,&
                        srrec,dpexh,srexh,wout,wrec,wexh,aa,bb,cc,&
                        bbb,ccc,bbbsq,diff,z,wrecp,wrecn,wexhp,woutp,&
                        wexhn,woutn,pmix,pext,dpd,wmou,tmix,gmix,&
                        wexhrev,wdummy,grec,trec,tout,gout,text,gext,&
                        pout,pexh,wmix,wext,cout,cexh,cmou,ropenout,&
                        ropenrec,ropenexh,ropenmou,fareaout,farearec,&
                        fareaexh,fareamou,fleakout,fleakrec,rdamper,&
                        dpqudlin,wqudlin
        integer      :: i,noninver,iparout,iparrec,iparexh,iparmou

! **** Read in inputs
        tout = xin(1)
        gout = xin(2)
        text = xin(3)
        gext = xin(4)
        pout = xin(5)
        pexh = xin(6)
        wmix = xin(7)
        wext = xin(8)
        cout = xin(9)
        cexh = xin(11)
        cmou = xin(12)
! **** Read and check parameters
        iparout = nint(par_v(1))
        iparrec = nint(par_v(2))
        iparexh = nint(par_v(3))
        iparmou = nint(par_v(4))
        ropenout= par_v(5)
        ropenrec= par_v(6)
        ropenexh= par_v(7)
        ropenmou= par_v(8)
        fareaout= par_v(9)
        farearec= par_v(10)
        fareaexh= par_v(11)
        fareamou= par_v(12)
        fleakout= par_v(13)
        fleakrec= par_v(14)
        fleakexh= par_v(15)
        fleakmou= par_v(16)
        rfixout = par_v(17)
        rfixrec = par_v(18)
        rfixexh = par_v(19)
        noninver= nint(par_v(20))

        if (noninver/=1) then        !  moved 12/22/99
            crec = 1.0-xin(10)
        else
            crec = xin(10)
        endif

! **** Calculate resistance of each branch as sum of damper resistance and
! **** fixed resistance
! **** Outside air branch - Resistance of main and minimum dampers in
! **** parallel is in series with the fixed resistance
        routd = rdamper(cout,ropenout,fleakout,fareaout,alegg(iparout),&
                       blegg(iparout),iparout)
        rmoud = rdamper(cmou,ropenmou,fleakmou,fareamou,alegg(iparmou),&
                       blegg(iparmou),iparmou)
        sqrroutd = sqrt(routd)
        sqrrmoud = sqrt(rmoud)
        routdp =1./(1./routd + 1./rmoud + 2./(sqrroutd*sqrrmoud))
        rout  = routdp + rfixout
        if (rout<=0.0) stop &
            'type326: resistance in outside branch not greater than 0'
! **** Recirc and exhaust branches - one damper in series with fixed
! **** resistance
        rrec = rdamper(crec,ropenrec,fleakrec,farearec,alegg(iparrec),&
                       blegg(iparrec),iparrec) + rfixrec
        if (rrec<=0.0) stop &
            'type326: resistance in recirc branch not greater than 0'
        rexh = rdamper(cexh,ropenexh,fleakexh,fareaexh,alegg(iparexh),&
                       blegg(iparexh),iparexh) + rfixexh
        if (rexh<=0.0) stop &
            'type326: resistance in exhaust branch not greater than 0'
! **** Calculate "signed" resistances (see printed documentation)
! **** test for direction of the outside, return and exhaust flows by
! **** determining the sign of the pressure difference across the
! **** corresponding resistance assuming the flow in question is zero.
        dpexhout = pexh - pout
! **** Outside:
        dpout = rrec*wmix*abs(wmix) - rexh*(wext-wmix)*abs(wext-wmix) -&
                dpexhout
        srout = sign(rout,dpout)
! **** Return:
        dprec = rexh*wext*abs(wext) + rout*wmix*abs(wmix) + dpexhout
        srrec = sign(rrec,dprec)
! **** Exhaust:
        dpexh = rrec*wext*abs(wext) - rout*(wmix-wext)*abs(wmix-wext) -&
                dpexhout
        srexh = sign(rexh,dpexh)
! **** Determine necessity of solving quadratic - Check for special cases
        if (srout==0) then
! **** Outside flow zero
            wout = 0.0
            wrec = wmix
            wexh = wext - wrec
        elseif (srrec==0) then
! **** Recirculation flow zero
            wrec = 0.0
            wout = wmix
            wexh = wext
        elseif (srexh==0) then
! **** Exhaust flow zero
            wexh = 0.0
            wrec = wext
            wout = wmix -wrec
        else
! **** All three unknown flows non-zero:
! **** Calc wrec by solving a quadratic of the form a*wrec**2 + b*wrec + c = 0
! **** Note that aa=a, bb=-b/2, cc=c,  bbb=-b/2a,  ccc=c/a
            aa   = srout+srexh-srrec
            bb   = srout*wmix+srexh*wext
            cc   = srout*wmix*wmix+srexh*wext*wext+dpexhout
            if (aa==0.) then
! **** One root
                wrec = cc/(2.0*bb)
            else
! **** Two roots
                bbb = bb/aa
                ccc = cc/aa
                bbbsq = bbb*bbb
                diff = bbbsq-ccc
                if (diff>0) then
! **** Roots unequal
                    z = sqrt(diff)
                    wrecp = bbb + z
                    wrecn = bbb - z
! **** Check each root for consistency with directions of air flow
! **** already determined
                    wexhp = wext - wrecp
                    woutp = wmix - wrecp
                    if (wexhp*srexh<0.0 .or.&
                        wrecp*srrec<0.0 .or.&
                        woutp*srout<0.0) then
                        rootpos = .false.
                    else
                        rootpos = .true.
                    endif
                    wexhn = wext-wrecn
                    woutn = wmix-wrecn
                    if (wexhn*srexh<0.0 .or.&
                        wrecn*srrec<0.0 .or.&
                        woutn*srout<0.0) then
                        rootneg = .false.
                    else
                        rootneg = .true.
                    endif
                    if (.not.rootpos .and. .not.rootneg) then
                        stop 'type326: no consistent solution'
                    elseif (rootpos .and. rootneg) then
                           stop 'type326: ambiguous solution'
                       elseif (rootpos .and. .not.rootneg) then
                           wrec = wrecp
                       elseif (rootneg .and. .not.rootpos) then
                           wrec = wrecn
                    else
                        stop 'type326: logic fault 1'
                    endif
                elseif (diff==0) then
! **** Equal roots
                    wrec = bbb
                else
! **** Complex roots (coding error, negative resistances ...?)
!         write(*,*) (xin(iii), iii=1,ni)
!         write(*,*) (par(iii), iii=1,np)
!         write(*,*) rout, rrec, rexh
!         write(*,*) srout, srrec, srexh
!         write(*,*) bbbsq, ccc, diff
                    stop  'type 326: complex roots'
                endif
            endif
! **** Calculate remaining unknown flows
            wout = wmix - wrec
            wexh = wext - wrec
! **** Calculate mixed and extract pressures
            pmix = pout - srout*wout*wout
            pext = pexh + srexh*wexh*wexh
        endif
! **** Flow rate through minimum outside air damper
        dpd  = dpqudlin(routdp,wout)
        wmou = wqudlin(rmoud,dpd)
! **** Having calculated flows, calculate mixed temperature and humidity
! **** Model assumes wmix and wext both non-negative. Warn if not true.
        if (lastit==1) then
            if (wmix<-wsmall) write(*,*) 'type 326: wmix < ',-wsmall
            if (wext<-wsmall) write(*,*) 'type 326: wext < ',-wsmall
        endif
! **** Check for zero or reverse mixed flow
        if (wmix<=wsmall) then
! **** Reverse mixed flow, take mixed temperature and humidity ratio to be
! **** an average of the outside and extract conditions so as to provide a
! **** reasonable guess for subsequent iterations
            tmix = (tout+text)/2.0
            gmix = (gout+gext)/2.0
        elseif (wmix>wsmall .and. wext<=wsmall) then
! **** Zero/reverse extract flow - Mixed condition is the same as outside
! **** condition
            tmix = tout
            gmix = gout
        elseif (wmix>wsmall .and. wext>wsmall ) then
! **** Normal, forward mixed and extract flows - Calculate mixed conditions
! **** first calculate recirc air conditions if recirc flow positive
            if (wrec>wsmall) then
! **** Test for zero/reverse exhaust flow
                if (wexh<wsmall) then
! **** Reverse exhaust flow - Recirc flow is a mixture of extract and outside
! **** conditions
                    wexhrev=-wexh
                    t(1) = text
                    gx(1) = gext
                    w(1) = wext
                    t(2) = tout
                    gx(2) = gout
                    w(2) = wexhrev
                    call moistmix(t,gx,w,2,trec,grec,wdummy)
                else
! **** Normal, forward exhaust flow - Recirc condition is extract condition
                    trec = text
                    grec = gext
                endif
            endif
! **** Calculate mixed air conditions
! **** Mixed condition is a combination of the recirc and outside conditions
! **** moistmix also tests for reverse flows and sets the mixed condition
! **** to the forward flow condition, or to an average of the flow conditions
! **** if both are reverse
            t(1) = trec
            gx(1) = grec
            w(1) = wrec
            t(2) = tout
            gx(2) = gout
            w(2) = wout
            call moistmix(t,gx,w,2,tmix,gmix,wdummy)
        else
            stop 'type326: logic fault 2'
        endif
! **** Assign outputs
        yout(1) = tmix
        yout(2) = gmix
        yout(3) = pmix
        yout(4) = pext
        yout(5) = wout
        yout(6) = wrec
        yout(7) = wexh
        yout(8) = wmou
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo

        return
        end subroutine type326

! *********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Mixing box with minimum outside air damper
! *
! * PURPOSE:    Calculates the mixed air flow rate, temperature and
! *             humidity and the extract air pressure. Uses Legg's
! *             correlations for the resistance of parallel and opposed
! *             blade dampers. Includes an additional damper in parallel
! *             with the main outside air damper.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. tri     : Outside air dry bulb temperature                   (C)
! *  2. gout    : Outside air humidity ratio                     (kg/kg)
! *  3. text    : Extract air dry bulb temperature                   (C)
! *  4. gext    : Extract air humidity ratio                     (kg/kg)
! *  5. pout    : Outside air intake gauge pressure                (kPa)
! *  6. pexh    : Exhaust air outlet gauge pressure                (kPa)
! *  7. pmix    : Mixed air gauge pressure                         (kPa)
! *  8. wext    : Extract dry air mass flow rate                  (kg/s)
! *  9. cout    : Outside air damper position (0=closed, 1=open)     (-)
! * 10. cre     : Recirc air damper position (0=open if PAR(16)=0)   (-)
! * 11. cexh    : Extract air damper position (0=closed, 1=open)     (-)
! * 12. cmou    : Minimum outside air damper pos. (0=closed,1=open)  (-)
! *
! * OUTPUTS
! * =======
! *  1. tmix    : Mixed air dry bulb temperature                     (C)
! *  2. gmix    : Mixed air humidity ratio                       (kg/kg)
! *  3. wmix    : Mixed air mass flow rate                        (kg/s)
! *  4. pext    : Extract air gauge pressure                       (kPa)
! *  5. wout    : Total outside air mass flow rate                (kg/s)
! *  6. wrec    : Recirc air mass flow rate                       (kg/s)
! *  7. wexh    : Exhaust air mass flow rate                      (kg/s)
! *  8. wmou    : Mass flow rate thru minimum outside air damper  (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. iparout : Outside air damper: opposed (0) or parallel (1)    (-)
! *  2. iparrec : Recirc air damper: opposed (0) or parallel (1)     (-)
! *  3. iparexh : Exhaust air damper: opposed (0) or parallel (1)    (-)
! *  4. iparmou : Min outside air damper: opposed (0) or parallel (1)(-)
! *  5. ropenout: Open resist. of outside air damper         (0.001/k.m)
! *  6. ropenrec: Open resist. of recirc air damper          (0.001/k.m)
! *  7. ropenexh: Open resist. of exhaust air damper         (0.001/k.m)
! *  8. ropenmou: Open resist. of minimum outside air damper (0.001/k.m)
! *  9. fareaout: Face area of outside air damper                   (m2)
! * 10. farearec: Face area of recirc air damper                    (m2)
! * 11. fareaexh: Face area of exhaust air damper                   (m2)
! * 12. fareamou: Face area of minimum outside air damper           (m2)
! * 13. fleakout: Leakage for outside air damper (fraction of full flow)
! * 14. fleakrec: Leakage for recirc air damper  (fraction of full flow)
! * 15. fleakexh: Leakage for exhaust air damper (fraction of full flow)
! * 16. fleakmou: Leakage for minimum outside air damper      (fraction)
! * 17. rfixout :  Fixed resistance in outside air branch    (0.001/k.m)
! * 18. rfixrec :  Fixed resistance in recirc air branch     (0.001/k.m)
! * 19. rfixexh :  Fixed resistance in exhaust air branch    (0.001/k.m)
! * 20. noninver: 0=invert recirc air damper, 1=not inverted         (-)
! *
! **********************************************************************
! 
!   MAJOR RESTRICTIONS:  Assumes positive mixed and extract flow rates
! 
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
! 
!   DATE:                November 29, 1995
! 
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  FLOWSPLT, MOISTMIX
!   FUNCTIONS  CALLED:   None
! 
!   REVISION HISTORY:    None
! 
!   REFERENCE:
! 
! **********************************************************************
! *
! * INTERNAL VARIABLES
! * ==================
! * routd   : resistance of main outside air damper
! * rmoud   : resistance of minimum outside air damper
! * routdp  : resistance of main and minimum dampers in parallel
! * rout    : total resistance of outside air branch
! * rrec    : total resistance of recirculation branch
! * rexh    : total resistance of exhaust branch
! * dp      : pressure drop across outside air branch
! * dpd     : pressure drop across outside air dampers
! * small   : threshold for significant flow rate
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type327(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=12,no=8,np=20,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real,dimension(20)    :: t,gx,w
        real,dimension(0:1)   :: alegg=(/-1.51,-1.51/),&
                                 blegg=(/0.105,0.0842/)
        real         :: small=1.0e-10,wcritl=3.e-2,wcritu=6.e-2
        real         :: ropenexh,ropenmou,fareaout,farearec,fareaexh,&
                        fareamou,fleakout,fleakrec,fleakexh,fleakmou,&
                        rfixout,rfixrec,rfixexh,tout,gout,text,&
                        gext,pout,pexh,pmix,wext,cout,crec,cexh,cmou,&
                        routd,rmoud,sqrroutd,sqrrmoud,routdp,rout,rrec,&
                        rexh,wexh,wrec,pext,dp,wout,dpd,wmou,wmix,tmix,&
                        gmix,wexhrev,wdummy,grec,trec,ropenout,ropenrec,&
                        rdamper,wqudlin,dpqudlin
        integer      :: i,noninver,ifail,iparout,iparrec,iparexh,iparmou

! **** Read and check parameters
        iparout = nint(par_v(1))
        iparrec = nint(par_v(2))
        iparexh = nint(par_v(3))
        iparmou = nint(par_v(4))
        ropenout= par_v(5)
        ropenrec= par_v(6)
        ropenexh= par_v(7)
        ropenmou= par_v(8)
        fareaout= par_v(9)
        farearec= par_v(10)
        fareaexh= par_v(11)
        fareamou= par_v(12)
        fleakout= par_v(13)
        fleakrec= par_v(14)
        fleakexh= par_v(15)
        fleakmou= par_v(16)
        rfixout = par_v(17)
        rfixrec = par_v(18)
        rfixexh = par_v(19)
        noninver= nint(par_v(20))
! **** Read in inputs
        tout = xin(1)
        gout = xin(2)
        text = xin(3)
        gext = xin(4)
        pout = xin(5)
        pexh = xin(6)
        pmix = xin(7)
        wext = xin(8)
        cout = xin(9)
        if (noninver/=1) then
            crec = 1.0-xin(10)
        else
            crec = xin(10)
        endif
        cexh = xin(11)
        cmou = xin(12)
! **** Calculate resistance of each branch as sum of damper resistance and
! **** fixed resistance
! **** Outside air branch - Resistance of main and minimum dampers in
! **** parallel is in series with the fixed resistance
        routd = rdamper(cout,ropenout,fleakout,fareaout,alegg(iparout),&
                       blegg(iparout),iparout)
        rmoud  = rdamper(cmou,ropenmou,fleakmou,fareamou,alegg(iparmou),&
                       blegg(iparmou),iparmou)
        sqrroutd = sqrt(routd)
        sqrrmoud = sqrt(rmoud)
        routdp =1./(1./routd + 1./rmoud + 2./(sqrroutd*sqrrmoud))
        rout  = routdp + rfixout
        if (rout<=0.0) stop &
            'type327: resistance in outside branch not greater than 0'
! **** Recirc and exhaust branches - One damper in series with fixed
! **** resistance
        rrec  = rdamper(crec,ropenrec,fleakrec,farearec,alegg(iparrec),&
                        blegg(iparrec),iparrec) + rfixrec
        if (rrec<=0.0) stop &
            'type327: resistance in recirc branch not greater than 0'
        rexh  = rdamper(cexh,ropenexh,fleakexh,fareaexh,alegg(iparexh),&
                        blegg(iparexh),iparexh) + rfixexh
        if (rexh<=0.0) stop &
            'type327: resistance in exhaust branch not greater than 0'
! **** Calculate recirc and exhaust flow rates
! **** Call flow split routine
        call flowsplt(wext,pmix,pexh,0.0,rrec,rexh,wcritl,wcritu,rtolx,&
                      pext,wrec,wexh,ifail)
! **** Check for unsuccesful completion
        if (ifail==1) then
! **** One or more resistances negative
            stop 'type 327: resistances must not be negative'
        elseif (ifail==2) then
! **** Zero resistance for both outlet branches
            stop 'type 327: either rret or rexh must be non-zero'
        endif
! **** Total outside air mass flow rate
        dp   = pout-pmix
        wout = wqudlin(rout,dp)
! **** Flow rate through minimum outside air damper
        dpd  = dpqudlin(routdp,wout)
        wmou = wqudlin(rmoud,dpd)
! **** Calculate mixed air flow rate as sum of outside and recirc flow rates
        wmix = wout+wrec
! **** Having calculated flows, calculate mixed temperature and humidity
! **** Model assumes wmix and wext both non-negative. Warn if not true.
        if (wmix<-1.0e-4 .and. lastit==1)&
            write(*,*) 'type 327: wmix < -1.0e-4'
        if (wext<-1.0e-4 .and. lastit==1)&
            write(*,*) 'type 327: wext < -1.0e-4'
! **** Check for reverse mixed flow
        if (wmix<=small) then
! **** Reverse mixed flow, take mixed temperature and humidity ratio to be
! **** an average of the outside and extract conditions so as to provide a
! **** reasonable guess for subsequent iterations
            tmix = (tout+text)/2.0
            gmix = (gout+gext)/2.0
        elseif (wmix>small .and. wext<=small) then
! **** Reverse extract flow - Mixed condition is the same as outside condition
            tmix = tout
            gmix = gout
        elseif (wmix>small .and. wext>small ) then
! **** Normal, forward mixed and extract flows - Calculate mixed conditions
! **** first calculate recirc air conditions if recirc flow positive
            if (wrec>small) then
! **** Test for reverse exhaust flow
                if (wexh<small) then
! **** Reverse exhaust flow - Recirc flow is a mixture of extract and outside
! **** conditions
                    wexhrev = -wexh
                    wexhrev = -wexh
                    t(1) = text
                    gx(1) = gext
                    w(1) = wext
                    t(2) = tout
                    gx(2) = gout
                    w(2) = wexhrev
                    call moistmix(t,gx,w,2,trec,grec,wdummy)
                else
! **** Normal, forward exhaust flow - Recirc condition is extract condition
                    trec = text
                    grec = gext
                endif
            endif
! **** Calculate mixed air conditions
! **** Mixed condition is a combination of the recirc and outside conditions
! **** moistmix also tests for reverse flows and sets the mixed condition
! **** to the forward flow condition, or to an average of the flow conditions
! **** if both are reverse
            t(1) = trec
            gx(1) = grec
            w(1) = wrec
            t(2) = tout
            gx(2) = gout
            w(2) = wout
            call moistmix(t,gx,w,2,tmix,gmix,wdummy)
        else
            write(*,*) 'type 327: wmix =',wmix
            write(*,*) 'type 327: wext =', wext

            stop 'type327: logic fault'
        endif
! **** Assign outputs
        yout(1) = tmix
        yout(2) = gmix
        yout(3) = wmix
        yout(4) = pext
        yout(5) = wout
        yout(6) = wrec
        yout(7) = wexh
        yout(8) = wmou
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo

        return
        end subroutine type327

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:  Two port control valve
! *
! * PURPOSE:     Calculates the inlet pressure from the flow rate and
! *              outlet pressure
! **********************************************************************
! * INPUTS
! * ======
! *  1. ww      : water mass flow rate                            (kg/s)
! *  2. pwo     : outlet water pressure                            (kPa)
! *  3. vstem   : valve stem position                                (-)
! *
! * OUTPUTS
! * =======
! *  1. pwi     : inlet water pressure                             (kPa)
! *
! * PARAMETERS
! * ==========
! *  1. kv      : valve capacity (Kv)                    (m3/hr @ 1 bar)
! *  2. eqpchar : valve curvature parameter (0=linear)               (-)
! *  3. xeqp    : equal percentage/linear break point                (-)
! *  4. sv      : valve rangability                                  (-)
! *  5. cl      : valve leakage (fractional flow)                    (-)
! *  6. xlin    : linear/close-off break point                       (-)
! *
! **********************************************************************
! 
!   MAJOR RESTRICTIONS:  Assumes characteristic has three segments -
!                        close-off, linear and exponential. Linear
!                        characteristic obtained by setting linear/
!                        exponential break point to 1. Close-off
!                        segment is also linear.
! 
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
! 
!   DATE:                November 14, 1995
! 
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   RLINPORT, REQPPORT
! 
!   REVISION HISTORY:    None
! 
!   REFERENCE:           (Based on confidential information from a major
!                        manufacturer.)
! 
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type328(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=6,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: kv
! **** Minimum water-side resistance
        real         :: frsmall=1.0e-8

        real         :: ww,pwo,vstem,eqpchar,xeqp,sv,cl,xlin,frwvalve,&
                        rlinport,reqpport,pwi,dpqudlin

! **** Read in inputs
        ww      = xin(1)
        pwo     = xin(2)
        vstem   = xin(3)
! **** Read in parameters
        kv      = par_v(1)
        eqpchar = par_v(2)
        xeqp    = par_v(3)
        sv      = par_v(4)
        cl      = par_v(5)
        xlin    = par_v(6)
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
        if (frwvalve<frsmall) then
! **** Total resistance almost zero
            write(*,*) 'type328: resistance must not be < ',frsmall
            stop
        endif
! **** Calculate water inlet pressure from flow resistance
        pwi     = pwo+dpqudlin(frwvalve,ww)
! **** Assign output
        yout(1)  = pwi
! **** Allow freezing of algebraic variable
        iostat(1)=1

        return
        end subroutine type328
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:  Two port control valve
! *
! * PURPOSE:     Calculates the flow rate from the inlet pressure and
! *              outlet pressure
! **********************************************************************
! * INPUTS
! * ======
! *  1. pwi     : inlet water gauge pressure                       (kPa)
! *  2. pwo     : outlet water gauge pressure                      (kPa)
! *  3. vstem   : valve stem position                                (-)
! *
! * OUTPUTS
! * =======
! *  1. ww      : water mass flow rate                            (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. kv      : valve capacity (Kv)                    (m3/hr @ 1 bar)
! *  2. eqpchar : valve curvature parameter (0=linear)               (-)
! *  3. xeqp    : equal percentage/linear break point                (-)
! *  4. sv      : valve rangability                                  (-)
! *  5. cl      : valve leakage (fractional flow)                    (-)
! *  6. xlin    : linear/close-off break point                       (-)
! *
! **********************************************************************
! 
!   MAJOR RESTRICTIONS:  Assumes characteristic has three segments -
!                        close-off, linear and exponential. Linear
!                        characteristic obtained by setting linear/
!                        exponential break point to 1. Close-off
!                        segment is also linear.
! 
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
! 
!   DATE:                November 14, 1995
! 
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   RLINPORT, REQPPORT
! 
!   REVISION HISTORY:    None
! 
!   REFERENCE:           (Based on confidential information from a major
!                        manufacturer.)
! 
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type329(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=6,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: kv
! **** Minimum water-side resistance,
        real         :: frsmall=1.0e-8
        real         :: pwi,pwo,vstem,eqpchar,xeqp,sv,cl,xlin,&
                        frwvalve,rlinport,reqpport,dp,ww,wqudlin

! **** Read in inputs
        pwi     = xin(1)
        pwo     = xin(2)
        vstem   = xin(3)
! **** Read in parameters
        kv      = par_v(1)
        eqpchar = par_v(2)
        xeqp    = par_v(3)
        sv      = par_v(4)
        cl      = par_v(5)
        xlin    = par_v(6)
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
        if (frwvalve<frsmall) then
! **** Resistance almost zero
            write(*,*) 'type329: resistance must not be < ',frsmall
            stop
        endif
! **** Calculate water inlet pressure from flow resistance
        dp      = pwi-pwo
        ww      = wqudlin(frwvalve,dp)
! **** Assign output
        yout(1)  = ww
! **** Allow freezing of algebraic variable
        iostat(1)=1
        return
        end subroutine type329
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:  Three port mixing valve
! *
! * PURPOSE:     Calculates the pressure at each inlet port and the flow
! *              at the common port from the flows at each inlet port
! *              and the outlet pressure.
! **********************************************************************
! * INPUTS
! * ======
! *  1. wwflow  : flow rate through flow port                     (kg/s)
! *  2. wwbypas : flow rate through by-pass port                  (kg/s)
! *  3. pwo     : outlet water pressure                            (kPa)
! *  4. vstem   : valve stem position                                (-)
! *
! * OUTPUTS
! * =======
! *  1. pwflow  : inlet pressure at flow port                      (kPa)
! *  2. pwbypas : inlet pressure at by-pass port                   (kPa)
! *  3. wwo     : outlet water flow rate                          (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. ivaltype: valve type: 0=lin/lin, 1=eq%(flow)/lin(byp), 2=lin/eq%
! *  2. kv      : valve capacity (Kv)                    (m3/hr @ 1 bar)
! *  3. eqpchar : valve curvature parameter (equal percentage port)  (-)
! *  4. xeqp    : equal percentage/linear break point                (-)
! *  5. sv      : valve rangability                                  (-)
! *  6. cl      : valve leakage (fractional flow)                    (-)
! *  7. xlin    : linear/close-off break point                       (-)
! *
! **********************************************************************
! 
!   MAJOR RESTRICTIONS:
! 
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
! 
!   DATE:                November 14, 1995
! 
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   RLINPORT,REQPPORT
! 
!   REVISION HISTORY:    None
! 
!   REFERENCE:           (Based on confidential information from a
!                        major manufacturer.)
! 
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type330(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=4,no=3,np=7,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: kv
! **** Minimum water-side resistance,
        real         :: frsmall=1.0e-8

        real         :: wwflow,wwbypas,pwo,vstem,eqpchar,xeqp,sv,cl,&
                        xlin,vstembyp,frvalflo,rlinport,reqpport,&
                        frvalbyp,pwflow,dpqudlin,pwbypas,wwo
        integer      :: i,ivaltype

! **** Read in inputs
        wwflow  = xin(1)
        wwbypas = xin(2)
        pwo     = xin(3)
        vstem   = xin(4)
! **** Read in parameters
        ivaltype= nint(par_v(1))
        kv      = par_v(2)
        eqpchar = par_v(3)
        xeqp    = par_v(4)
        sv      = par_v(5)
        cl      = par_v(6)
        xlin    = par_v(7)
! **** Valve types and resistances
! **** limit stem position to range 0-1
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
            stop 'type330: illegal valve type'
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
        if (frvalflo<frsmall .or. frvalbyp<frsmall) then
! **** Total resistance almost zero
            write(*,*) &
            'type330: water-side flow resistance must not be < ',frsmall
            stop
        endif
! **** Pressure at inlet to flow and by-pass ports
        pwflow  = pwo+dpqudlin(frvalflo,wwflow)
        pwbypas = pwo+dpqudlin(frvalbyp,wwbypas)
! **** Outlet flow rate
        wwo     = wwflow+wwbypas
! **** Assign output values
        yout(1)  = pwflow
        yout(2)  = pwbypas
        yout(3)  = wwo
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo
        return
        end subroutine type330
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     Variable speed drive
! *
! * PURPOSE:        Calculate the rotation speed of an ideal motor
! *                 connected to a variable speed drive incorporating a
! *                 rate limit.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. fspddem : demanded fractional speed                          (-)
! *
! * OUTPUTS
! * =======
! *  1. rspd    : actual rotation speed                          (rev/s)
! *
! * PARAMETERS
! * ==========
! *  1. rspdmax : maximum rotation speed                         (rev/s)
! *  2. ttran   : travel time (lim-lim)                              (s)
! *
! *  SAVED
! *  =====
! *  1. time    : time of previous call for rate limit
! *  2. fspdact : actual fractional speed at previous call
! *  3. fspdact : actual fractional speed at previous timestep
! *
! **********************************************************************
! 
!   MAJOR RESTRICTIONS:  Neglects motor dynamics and drive hysteresis
!                        and non-linearity. Does not treat power
!                        consumption or efficiency.
! 
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
! 
!   DATE:                June 16, 1994
! 
!   INCLUDE FILES:
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:
! 
!   REVISION HISTORY:    None
! 
!   REFERENCE:           Haves, P., Component-Based Modelling of VAV
!                        Systems, Proc. System Simulation in Buildings
!                        '94, Liege, Belgium, December 1994
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! * fspdactp: actual fractional speed at the previous step time
! * dfspdmax: maximum change in fractional speed in one time-step
! * dfspd   : actual change in fractional speed during time-step
! * fspdact : current actual fractional speed
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type333(xin,yout,par_v,saved_v,iostat)
 
        use modsim_head
!        implicit none
        integer,parameter                 :: ni=1,no=1,np=2,ns=3
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        logical      :: quick
        real         :: small=1.e-6
        
! **** Read in inputs
        fspddem  = xin(1)
! **** Read in parameters
        rspdmax  = par_v(1)
        ttran    = par_v(2)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
                do is = 2,ns-1,2
                    saved_v(is) = 0.0
                enddo
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - Update previous sample instant values
            do is=2,ns-1,2
                saved_v(is+1) = saved_v(is)
            enddo
        endif
! **** Determine curent position based on situation at previous step time
! **** Update previous values
        fspdactp = saved_v(3)
! **** Limit demanded fractional speed
        fspddem = max(0.,min(fspddem,1.))
! **** determine if demanded position attained
        if (ttran<=(tstep*small)) then
! **** Response is essentially instantaneous
            quick=.true.
        else
            quick=.false.
! **** Travel time is non-zero - Determine maximum change in fractional
! **** speed in one time-step
            dfspdmax=tstep/ttran
        endif
! **** Difference between current demanded speed and previous actual speed
        dfspd=fspddem-fspdactp
        if (quick .or. abs(dfspd)<=dfspdmax) then
! **** Demanded speed attained
            fspdact=fspddem
        else
! **** Demanded speed not attained
            fspdact=fspdactp+sign(dfspdmax,dfspd)
        endif
! **** Actual rotation speed
        rspd=fspdact*rspdmax
! **** Save time of current call
        saved_v(1)=time
! **** Save provisional values to be used in updating at next step time
        saved_v(2)=fspdact
! **** Output
        yout(1)=rspd
! **** Disallow freezing
        do i=1,no
            iostat(i)=0
        enddo
! **** Return
        return
        end subroutine type333
        
