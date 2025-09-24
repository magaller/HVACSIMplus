! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE:          TYPE 301 - Temperature Sensor
! *
! * PURPOSE:             Calculate the response of a first order linear
! *                      temperature sensor
! *
! * MAJOR RESTRICTIONS:  First order, linear
! *
! * DEVELOPER:           Li Mei and Philip Haves
! *                      Loughborough University
! *
! * LAST MODIFIED:       November 22, 1995
! *
! * INCLUDE FILES:       None
! * SUBROUTINES CALLED:  None
! * FUNCTIONS  CALLED:   None
! *
! * REFERENCE:           Based on HVACSIM+ Type7 (NIST)
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. ti      : temperature measured                               (C)
! *
! * OUTPUTS
! * =======
! *  1. co      : sensor output                                      (-)
! *
! * PARAMETERS
! * ==========
! *  1. tzero   : offset: input for zero output                      (C)
! *  2. tgain   : gain: change in input for 0->1 output              (K)
! *  3. tcon    : time constant                                      (s)
! *  4. cmax    : upper limit of output range                        (-)
! *  5. cmin    : lower limit of output range                        (-)
! *
! *  SAVED
! *  =====
! *  1. time    : time of previous call
! *  2. ts      : solution of differential equation from previous call
! *  3. tsp     : solution of differential equation from previous step
! *               time
! *
! * INTERNAL VARIABLES
! * ==================
! * (None)
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type301(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=1,no=1,np=5,&
                                             ndiffeq=1,ns=1+ndiffeq*2
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: ti,tgain,tcon,cmax,cmin,tsp,ts,aa,bb,tsba,ci,&
                        co,tsbar

! **** Read in input
        ti = xin(1)
! **** Read in parameters
        tzero = par_v(1)
        tgain = par_v(2)
        if (tgain==0.0) then
            stop 'type 301: temperature range must be non-zero'
        endif
        tcon = par_v(3)
        cmax = par_v(4)
        cmin = par_v(5)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.0
            endif
            if (init==0) then
! **** Sensor is in equilibrium with environment at start of simulation
                saved_v(2) = ti
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update value of sensor temperature from
! **** previous time-step
            saved_v(3) = saved_v(2)
        endif
! **** Update previous values
        tsp = saved_v(3)
! **** Determine temperature of sensor element
        if (tcon<=0.0) then
! **** Instantaneous response
            ts = ti
        else
! **** Non-zero time constant - integrate analytically
            aa = -1./tcon
            bb = ti/tcon
            call diffeq(time,aa,bb,tsp,ts,tsbar)
        endif
! **** Determine sensor output
! **** Calculate unclipped sensor output from range and offset
        ci = (ts - tzero)/tgain
! **** limit output to user-specified range
        if (ci>cmax) then
            co = cmax
        elseif (ci<cmin) then
            co = cmin
        else
            co = ci
        endif
! **** Save provisional value to be used in updating at next step time
        saved_v(2) = ts
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = co
! **** Determine whether to allow freezing
! **** Freezing of the output is allowed only if the input is a constant
! **** or a boundary variable and the change in the output is small
        if (tcon==0.0 .or. ((iostat(1)<-1).or.(iostat(1)==2))&
              .and. ((abs(ti - ts))<=(rtolx * abs(ts)+atolx))) then
            iostat(1) = 1
        else
            iostat(1) = 0
        endif
! **** Return
        return
        end subroutine type301

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE:          TYPE 302 - Humidity sensor
! *
! * PURPOSE:             Calculate the response of a first order linear
! *                      humidity sensor
! *
! * MAJOR RESTRICTIONS:  First order, linear
! *
! * DEVELOPER:           Li Mei and Philip Haves
! *                      Loughborough University of Technology
! *
! * LAST MODIFIED:       November 22, 1995
! *
! * INCLUDE FILES:       None
! * SUBROUTINES CALLED:  None
! * FUNCTIONS  CALLED:   FPHI,FTDEW,FTWB
! *
! * REFERENCE:
!
! **********************************************************************
! * INPUTS
! * ======
! *  1. gi      : humidity ratio                                 (kg/kg)
! *  2. ti      : temperature                                        (C)
! *
! * OUTPUTS
! * =======
! *
! *  2. co      : sensor output                                      (-)
! *
! * PARAMETERS
! * ==========
! *  1. mode    : humidity ratio=1, relative humidity=2, dew point=3,
! *               wet bulb=4, deg of saturation=5, specific enthalpy=6
! *  2. xzero   : offset: input for zero output        (sensed quantity)
! *  3. xgain   : gain: change in i/p for 0->1 o/p     (sensed quantity)
! *  4. tcon    : time constant                                      (s)
! *  5. cmax    : upper limit of output range                        (-)
! *  6. cmin    : lower limit of output range                        (-)
! *
! *  SAVED
! *  =====
! *  1. time    : time of previous call
! *  2. hs      : solution of differential equation from previous call
! *  3. hsp     : solution of differential equation from previous
! *               step-time
! *
! * INTERNAL VARIABLES
! * ==================
! * (None)
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type302(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=2,no=1,np=6,&
                                             ndiffeq=1,ns=1+ndiffeq*2
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

! **** Thermophysical constants
        real         :: patmpa=101325.0
        real         :: gi,ti,xzero,xgain,tcon,cmax,cmin,xi,fphi,ftdew,&
                        ftwb,gs,fwphi,fhair,xsp,xs,aa,bb,xsbar,ci,co
        integer      :: mode

! **** Read in input
        gi      = xin(1)
        ti      = xin(2)
! **** Read in parameters
        mode    = nint(par_v(1))
        xzero   = par_v(2)
        xgain   = par_v(3)
        if (xgain==0.0) then
            stop 'type 302: humidity range must be non-zero'
        endif
        tcon    = par_v(4)
        cmax    = par_v(5)
        cmin    = par_v(6)
! **** Calculate required measure of humidity
        select case (mode)
            case (1)
! **** Humidity ratio
               xi = gi
            case (2)
! **** Relative humidity
               xi = fphi(ti,gi,patmpa)
            case (3)
! **** Dew point temperature
               xi = ftdew(gi,patmpa)
            case (4)
! **** Wet bulb temperature
               xi = ftwb(ti,gi,patmpa)
            case (5)
! **** Degree of saturation
               gs  = fwphi(ti,100.,patmpa)
               xi = gi/gs
            case (6)
! **** Specific enthalpy
               xi = fhair(ti,gi)
            case default
        end select

! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Sensor is in equilibrium with environment at start of simulation
                saved_v(2) = xi
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update value of sensor temperature from
! **** previous time-step
            saved_v(3) = saved_v(2)
        endif
! **** Update previous values
        xsp = saved_v(3)
! **** Determine humidity of sensor element
        if (tcon<=0.0) then
! **** Instantaneous response
            xs = xi
        else
! **** Non-zero time constant - integrate analytically
            aa = -1./tcon
            bb = xi/tcon
            call diffeq(time,aa,bb,xsp,xs,xsbar)
        endif
! **** Determine sensor output
! **** Calculate unclipped sensor output from range and offset
        ci = (xs - xzero)/xgain
! **** Limit output to user-specified range
        if (ci>cmax) then
            co = cmax
        elseif (ci<cmin) then
            co = cmin
        else
            co = ci
        endif
! **** Save provisional value to be used in updating at next step time
        saved_v(2) = xs
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = co
! **** Determine whether to allow freezing
! **** Freezing of the output is allowed only if the input is a constant
! **** or a boundary variable and the change in the output is small
        if (tcon==0.0 .or. ((iostat(1)<-1).or.(iostat(1)==2))&
!            .and. ((iostat(2)<-1).or.(iostat(2)==2)) .and.& 
		!MAG comment out iostat(2) code, add .and. below
            .and.((abs(xi - xs))<=(rtolx * abs(xs)+atolx))) then
            iostat(1) = 1
        else
            iostat(1) = 0
        endif
! **** Return
        return
        end subroutine type302
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:   Flow rate sensor
! *
! * PURPOSE:      Calculate the response of a first order linear
! *               flow rate sensor
! **********************************************************************
! * INPUTS
! * ======
! *  1. wi      : mass flow rate                                  (kg/s)
! *
! * OUTPUTS
! * =======
! *  1. co      : sensor output                                      (-)
! *
! * PARAMETERS
! * ==========
! *  1. mode    : mass flow rate = 1, volumetric flow rate = 2
! *               velocity = 3, velocity pressure = 4                (-)
! *  2. ifluid  : air = 1, water = 2                                 (-)
! *  3. area    : cross-sectional area of duct or pipe              (m2)
! *  4. xzero   : offset: input when output is zero    (sensed quantity)
! *  5. xgain   : gain: input change for output 0->1   (sensed quantity)
! *  6. tcon    : time constant                                      (s)
! *  7. cmax    : upper limit of output range                        (-)
! *  8. cmin    : lower limit of output range                        (-)
! *
! * SAVED
! * =====
! *  1. time    : time of previous call                              (s)
! *  2. xsp     : solution of differential equation from previous call
! *  3. xsp     : solution of differential equation from previous step
! *               time
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  First order, linear
!
!   DEVELOPER:           Li Mei and Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 22, 1994
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************
!
!   INTERNAL VARIABLES:  None
!
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type303(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=1,no=1,np=8,&
                                             ndiffeq=1,ns=1+ndiffeq*2
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: rhoa=1.2,rhow=1000.
        real         :: wi,area,xzero,xgain,tcon,cmax,cmin,rho,xi,xsp,&
                        xs,aa,bb,xsbar,ci,co
        integer      :: mode,ifluid

! **** Read in input
        wi = xin(1)
! **** Read in parameters
        mode    = nint(par_v(1))
        ifluid  = nint(par_v(2))
        area    = par_v(3)
        xzero   = par_v(4)
        xgain   = par_v(5)
        if (xgain==0.0) then
            stop 'type 303: mass flow rate range must be non-zero'
        endif
        tcon    = par_v(6)
        cmax    = par_v(7)
        cmin    = par_v(8)
! **** Determine different modes for mass flow sensor
        if (ifluid==1) then
            rho = rhoa
        else
            rho = rhow
        endif
        if (area<=0.0) stop 'type303: area must be > 0'

        select case (mode)
            case (1)
               xi = wi

            case (2)
               xi = wi/rho
            case (3)
               xi = wi/(area*rho)
            case (4)
               xi = 0.5*wi/(rho*area*area)
            case default
        end select

! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Sensor is in equilibrium with environment at start of simulation
                saved_v(2) = xi
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update value of sensor temperature from
! **** previous time-step
            saved_v(3) = saved_v(2)
        endif
! **** Update previous values
        xsp = saved_v(3)
! **** Determine state of sensor element
        if (tcon<=0.0) then
! **** Instantaneous response
            xs = xi
        else
! **** Non-zero time constant - integrate analytically
            aa = -1./tcon
            bb = xi/tcon
            call diffeq(time,aa,bb,xsp,xs,xsbar)
        endif
! **** Determine sensor output
! **** Calculate unclipped sensor output from range and offset
        ci = (xs - xzero)/xgain
! **** Limit output to user-specified range
        if (ci>cmax) then
            co = cmax
        elseif (ci<cmin) then
            co = cmin
        else
            co = ci
        endif
! **** Save provisional value to be used in updating at next step time
        saved_v(2) = xs
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = co
! **** Determine whether to allow freezing
! **** Freezing of the output is allowed only if the input is a constant
! **** or a boundary variable and the change in the output is small
        if (tcon==0.0 .or. ((iostat(1)<-1).or.(iostat(1)==2))&
              .and. ((abs(xi-xs))<=(rtolx * abs(xs)+atolx))) then
            iostat(1) = 1
        else
            iostat(1) = 0
        endif
! **** Return
        return
        end subroutine type303
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:  Total pressure sensor
! *
! * PURPOSE:     Calculate the response of a first order linear
! *              total pressure sensor
! **********************************************************************
! * INPUTS
! * ======
! *  1. pi      : pressure input                                   (kPa)
! *
! * OUTPUTS
! * =======
! *  1. co      : sensor output                                      (-)
! *
! * PARAMETERS
! * ==========
! *  1. pzero   : offset: input for zero output                    (kPa)
! *  2. pgain   : gain: change in input for 0->1 output            (kPa)
! *  3. tcon    : time constant                                      (s)
! *  4. cmax    : upper limit of output range                        (-)
! *  5. cmin    : lower limit of output range                        (-)
! *
! * SAVED
! * =====
! *  1. time    : time of previous call
! *  2. ps      : solution of differential equation from previous call
! *  3. psp     : solution of differential equation from previous step
! *               time
! *
! ***********************************************************************
!
!   MAJOR RESTRICTIONS:  First order, linear
!
!   DEVELOPER:           Li Mei and Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 29, 1994
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************
!
!   INTERNAL VARIABLES:  None
!
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type304(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=1,no=1,np=5,&
                                             ndiffeq=1,ns=1+ndiffeq*2
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: pi,pzero,pgain,tcon,cmax,cmin,psp,ps,aa,bb,&
                        psbar,ci,co

! **** Read in input
        pi    = xin(1)
! **** Read in parameters
        pzero = par_v(1)
        pgain = par_v(2)
        if (pgain==0.0) then
            stop 'type 304: pressure range must be non-zero'
        endif
        tcon  = par_v(3)
        cmax  = par_v(4)
        cmin  = par_v(5)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                 saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Sensor is in equilibrium with environment at start of simulation
                saved_v(2) = pi
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update value of sensor temperature from
! **** previous time-step
            saved_v(3) = saved_v(2)
        endif
! **** Update previous values
        psp = saved_v(3)
! **** Determine state of sensor element
        if (tcon<=0.0) then
! **** Instantaneous response
            ps = pi
        else
! **** Non-zero time constant - integrate analytically
            aa = -1./tcon
            bb = pi/tcon
            call diffeq(time,aa,bb,psp,ps,psbar)
        endif
! **** Determine sensor output
! **** Calculate unclipped sensor output from range and offset
        ci = (ps - pzero)/pgain
! **** Limit output to user-specified range
        if (ci>cmax) then
            co = cmax
        elseif (ci<cmin) then
            co = cmin
        else
            co = ci
        endif
! **** Save provisional value to be used in updating at next step time
        saved_v(2) = ps
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = co
! **** Determine whether to allow freezing
! **** Freezing of the output is allowed only if the input is a constant
! **** or a boundary variable and the change in the output is small
        if (tcon==0.0 .or. ((iostat(1)<-1).or.(iostat(1)==2))&
               .and. ((abs(pi - ps))<=(rtolx * abs(ps)+atolx))) then
            iostat(1) = 1
        else
            iostat(1) = 0
        endif
! **** Return
        return
        end subroutine type304


! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:  Static pressure sensor
! *
! * PURPOSE:     Calculate the response of a first order linear
! *              static pressure sensor
! *********************************************************************
! * INPUTS
! * ======
! *  1. ptot    : pressure input (total)                          (kPa)
! *  2. wi      : mass flow rate                                 (kg/s)
! *
! * OUTPUTS
! * =======
! *  1. co      : sensor output                                     (-)
! *
! * PARAMETERS
! * ==========
! *  1. ifluid  : air = 1, water = 2                                (-)
! *  2. area    : cross-sectional area                             (m2)
! *  3. pzero   : offset: input for zero output       (sensed quantity)
! *  4. pgain   : gain: change in input for 0->1 output (sensed quantity)
! *  5. tcon    : time constant                                     (s)
! *  6. cmax    : upper limit of output range                       (-)
! *  7. cmin    : lower limit of output range                       (-)
! *
! * SAVED
! * =====
! *  1. time    : time of previous call
! *  2. ps      : solution of differential equation from previous call
! *  3. psp     : solution of differential equation from previous step
! *               time
! ***********************************************************************
!
!   MAJOR RESTRICTIONS:  First order, linear
!
!   DEVELOPER:           Li Mei and Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 29, 1994
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
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

        subroutine type305(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=2,no=1,np=7,&
                                             ndiffeq=1,ns=1+ndiffeq*2
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: rhoa=1.2, rhow=1000.
        real         :: ptot,wi,area,pzero,pgain,tcon,cmax,cmin,&
                        rho,pi,psp,ps,aa,bb,psbar,ci,co
        integer      :: ifluid

! **** Read in input
        ptot    = xin(1)
        wi      = xin(2)
! **** Read in parameters
        ifluid  = nint(par_v(1))
        area    = par_v(2)
        pzero   = par_v(3)
        pgain   = par_v(4)
        tcon    = par_v(5)
        cmax    = par_v(6)
        cmin    = par_v(7)
! **** Check air or water flow
        if (ifluid==1) then
            rho = rhoa
        elseif (ifluid==2) then
            rho = rhow
        else
            stop 'type 305: ifluid must be 1 or 2'
        endif
! *****Calculate static pressure
        pi = ptot - 0.001*0.5*wi*wi/(rho*area*area)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Sensor is in equilibrium with environment at start of simulation
                saved_v(2) = pi
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update value of sensor temperature from
! **** previous time-step
            saved_v(3) = saved_v(2)
        endif
! **** Update previous values
        psp = saved_v(3)
! **** Determine state of sensor element
        if (tcon<=0.0) then
! **** Instantaneous response
            ps = pi
        else
! **** Non-zero time constant - integrate analytically
            aa = -1./tcon
            bb = pi/tcon
        call diffeq(time,aa,bb,psp,ps,psbar)
        endif
! **** Determine sensor output
! **** Calculate unclipped sensor output from range and offset
        ci = (ps - pzero)/pgain
! **** Limit output to user-specified range
        if (ci>cmax) then
            co = cmax
        elseif (ci<cmin) then
            co = cmin
        else
            co = ci
        endif
! **** Save provisional value to be used in updating at next step time
        saved_v(2) = ps
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = co
! **** Determine whether to allow freezing
! **** Freezing of the output is allowed only if the input is a constant
! **** or a boundary variable and the change in the output is small
        if (tcon==0.0 .or. ((iostat(1)<-1).or.(iostat(1)==2))&
            !.and. ((iostat(2)<-1).or.(iostat(2)==2)) .and.&
		! MAG comment out iostat(2) code, add .and. below
            .and.((abs(pi-ps))<=(rtolx * abs(ps)+atolx))) then
            iostat(1) = 1
        else
            iostat(1) = 0
        endif
! **** Return
        return
    end subroutine type305

    
    
! * SUBROUTINE:  Pressure modification
! *
! * PURPOSE:     Allow modifications to pressure values for Cx
! *********************************************************************
! * INPUTS
! * ======
! *  1. p1      : Pressure input 1                            (kPa)
! *  2. p2      : Pressure input 2                            (kPa)
! *  3. c1      : Control (type) input                           (-)
! *     0- Output p1 (default) pass-thru
! *     1- Output p2           swap
! *     2- Output p1 + p2      add (i.e. for gradient)
! *
! * OUTPUTS
! * =======
! *  1. to      : Pressure output                                (kPa)
! *
! * PARAMETERS
! * ==========
! *  1. dummy
! *
! * SAVED
! * =====
! ***********************************************************************
!
!   MAJOR RESTRICTIONS:  NOne
!
!   DEVELOPER:           Michael A. Galler
!                        NIST
!
!   DATE:                September 15, 2024
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************

        subroutine type391(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=1,&
                                             ndiffeq=0,ns=0
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: p1,p2,co
        integer      :: c1

! **** Read in input
        p1    = xin(1)
        p2    = xin(2)
        c1    = nint(xin(3))

        if(c1 .eq. 1) then
            co = p2
        else if(c1 .eq. 2) then 
            co = p1 + p2
        else 
            co = p1
        end if
        
        iostat(1) = 1
        yout(1) = co
        
        return
    end subroutine type391
    
    
  
! * SUBROUTINE:  Flow modification
! *
! * PURPOSE:     Allow modifications to flow values for Cx
! *********************************************************************
! * INPUTS
! * ======
! *  1. w1      : Flow input 1                                   (kg/s)
! *  2. w2      : Flow input 2                                   (kg/s)
! *  3. c1      : Control (type) input                           (-)
! *     0- Output w1 (default) pass-thru
! *     1- Output w2           swap
! *     2- Output w1 + w2      add (i.e. for gradient)
! *
! * OUTPUTS
! * =======
! *  1. to      : Flow output                                (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. dummy
! *
! * SAVED
! * =====
! ***********************************************************************
!
!   MAJOR RESTRICTIONS:  NOne
!
!   DEVELOPER:           Michael A. Galler
!                        NIST
!
!   DATE:                September 15, 2024
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************

        subroutine type392(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=1,&
                                             ndiffeq=0,ns=0
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: w1,w2,co
        integer      :: c1

! **** Read in input
        w1    = xin(1)
        w2    = xin(2)
        c1    = nint(xin(3))

        if(c1 .eq. 1) then
            co = w2
        else if(c1 .eq. 2) then 
            co = w1 + w2
        else 
            co = w1
        end if
        
        iostat(1) = 1
        yout(1) = co
        
        return
    end subroutine type392
    
    
! * SUBROUTINE:  Temperature modification
! *
! * PURPOSE:     Allow modifications to temperature values for Cx
! *********************************************************************
! * INPUTS
! * ======
! *  1. t1      : Temperature input 1                            (C)
! *  2. t2      : Temperature input 2                            (C)
! *  3. c1      : Control (type) input                           (-)
! *     0- Output t1 (default) pass-thru
! *     1- Output t2           swap
! *     2- Output t1 + t2      add (i.e. for gradient)
! *
! * OUTPUTS
! * =======
! *  1. to      : temperature output                                (C)
! *
! * PARAMETERS
! * ==========
! *  1. dummy
! *
! * SAVED
! * =====
! ***********************************************************************
!
!   MAJOR RESTRICTIONS:  NOne
!
!   DEVELOPER:           Michael A. Galler
!                        NIST
!
!   DATE:                September 15, 2024
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************

        subroutine type393(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=1,&
                                             ndiffeq=0,ns=0
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: t1,t2,co
        integer      :: c1

! **** Read in input
        t1    = xin(1)
        t2    = xin(2)
        c1    = nint(xin(3))

        if(c1 .eq. 1) then
            co = t2
        else if(c1 .eq. 2) then 
            co = t1 + t2
        else 
            co = t1
        end if
        
        iostat(1) = 1
        yout(1) = co
        
        return
    end subroutine type393

        
! * SUBROUTINE:  Control modification
! *
! * PURPOSE:     Allow modifications to control values for Cx
! *********************************************************************
! * INPUTS
! * ======
! *  1. c1      : Control input 1                            (-)
! *  2. c2      : Control input 2                            (-)
! *  3. c3      : Control (type) input                       (-)
! *     0- Output c1 (default) pass-thru
! *     1- Output c2           swap
! *     2- Output c1 + c2      add (i.e. for gradient)
! *
! * OUTPUTS
! * =======
! *  1. co      : control output                             (-)
! *
! * PARAMETERS
! * ==========
! *  1. dummy
! *
! * SAVED
! * =====
! ***********************************************************************
!
!   MAJOR RESTRICTIONS:  NOne
!
!   DEVELOPER:           Michael A. Galler
!                        NIST
!
!   DATE:                September 15, 2024
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************

        subroutine type394(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=1,&
                                             ndiffeq=0,ns=0
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: c1,c2,co
        integer      :: c3

! **** Read in input
        c1    = xin(1)
        c2    = xin(2)
        c3    = nint(xin(3))

        if(c3 .eq. 1) then
            co = c2
        else if(c3 .eq. 2) then 
            co = c1 + c2
        else 
            co = c1
        end if
        
        iostat(1) = 1
        yout(1) = co
        
        return
    end subroutine type394
    

! * SUBROUTINE:  Other type modification
! *
! * PURPOSE:     Allow modifications to 5-other values for Cx
! *********************************************************************
! * INPUTS
! * ======
! *  1. m1      : Other type input 1                            (-)
! *  2. m2      : Other type input 2                            (-)
! *  3. c3      : Control (type) input                       (-)
! *     0- Output m1 (default) pass-thru
! *     1- Output m2           swap
! *     2- Output m1 + m2      add (i.e. for gradient)
! *
! * OUTPUTS
! * =======
! *  1. co      : control output                             (-)
! *
! * PARAMETERS
! * ==========
! *  1. dummy
! *
! * SAVED
! * =====
! ***********************************************************************
!
!   MAJOR RESTRICTIONS:  NOne
!
!   DEVELOPER:           Michael A. Galler
!                        NIST
!
!   DATE:                September 15, 2024
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************

        subroutine type395(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=1,&
                                             ndiffeq=0,ns=0
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         ::m1,m2,co
        integer      :: c1

! **** Read in input
        m1    = xin(1)
        m2    = xin(2)
        c1    = nint(xin(3))

        if(c1 .eq. 1) then
            co = m2
        else if(c1 .eq. 2) then 
            co = m1 + m2
        else 
            co = m1
        end if
        
        iostat(1) = 1
        yout(1) = co
        
        return
    end subroutine type395

    
! * SUBROUTINE:  Energy type modification
! *
! * PURPOSE:     Allow modifications to 6-Energy values for Cx
! *********************************************************************
! * INPUTS
! * ======
! *  1. e1      : Energy type input 1                            (kJ)
! *  2. e2      : Energy type input 2                            (kJ)
! *  3. c3      : Control (type) input                       (-)
! *     0- Output e1 (default) pass-thru
! *     1- Output e2           swap
! *     2- Output e1 + e2      add (i.e. for gradient)
! *
! * OUTPUTS
! * =======
! *  1. eo      : energy output                             (kJ)
! *
! * PARAMETERS
! * ==========
! *  1. dummy
! *
! * SAVED
! * =====
! ***********************************************************************
!
!   MAJOR RESTRICTIONS:  None
!
!   DEVELOPER:           Michael A. Galler
!                        NIST
!
!   DATE:                September 15, 2024
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************

        subroutine type396(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=1,&
                                             ndiffeq=0,ns=0
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         ::e1,e2,eo
        integer      ::c1

! **** Read in input
        e1    = xin(1)
        e2    = xin(2)
        c1    = nint(xin(3))

        if(c1 .eq. 1) then
            eo = e2
        else if(c1 .eq. 2) then 
            eo = e1 + e2
        else 
            eo = e1
        end if
        
        iostat(1) = 1
        yout(1) = eo
        
        return
    end subroutine type396
    
    
    ! * SUBROUTINE:  Power type modification
! *
! * PURPOSE:     Allow modifications to 7-power values for Cx
! *********************************************************************
! * INPUTS
! * ======
! *  1. q1      : Power type input 1                            (kW)
! *  2. q2      : Power type input 2                            (kW)
! *  3. c3      : Control (type) input                       (-)
! *     0- Output q1 (default) pass-thru
! *     1- Output q2           swap
! *     2- Output q1 + q2      add (i.e. for gradient)
! *
! * OUTPUTS
! * =======
! *  1. qo      : Power output                             (kW)
! *
! * PARAMETERS
! * ==========
! *  1. dummy
! *
! * SAVED
! * =====
! ***********************************************************************
!
!   MAJOR RESTRICTIONS:  None
!
!   DEVELOPER:           Michael A. Galler
!                        NIST
!
!   DATE:                September 15, 2024
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************

        subroutine type397(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=1,&
                                             ndiffeq=0,ns=0
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         ::q1,q2,qo
        integer      ::c1

! **** Read in input
        q1    = xin(1)
        q2    = xin(2)
        c1    = nint(xin(3))

        if(c1 .eq. 1) then
            qo = q2
        else if(c1 .eq. 2) then 
            qo = q1 + q2
        else 
            qo = q1
        end if
        
        iostat(1) = 1
        yout(1) = qo
        
        return
    end subroutine type397
    
! * SUBROUTINE:  Humidity type modification
! *
! * PURPOSE:     Allow modifications to 8-humidity values for Cx
! *********************************************************************
! * INPUTS
! * ======
! *  1. g1      : Other type input 1                            (kg/kg)
! *  2. g2      : Other type input 2                            (kg/kg)
! *  3. c3      : Control (type) input                       (-)
! *     0- Output g1 (default) pass-thru
! *     1- Output g2           swap
! *     2- Output g1 + g2      add (i.e. for gradient)
! *
! * OUTPUTS
! * =======
! *  1. go      : Humidity output                             (kg/kg)
! *
! * PARAMETERS
! * ==========
! *  1. dummy
! *
! * SAVED
! * =====
! ***********************************************************************
!
!   MAJOR RESTRICTIONS:  NOne
!
!   DEVELOPER:           Michael A. Galler
!                        NIST
!
!   DATE:                September 15, 2024
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************

        subroutine type398(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=1,np=1,&
                                             ndiffeq=0,ns=0
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         ::g1,g2,go
        integer      :: c1

! **** Read in input
        g1    = xin(1)
        g2    = xin(2)
        c1    = nint(xin(3))

        if(c1 .eq. 1) then
            go = g2
        else if(c1 .eq. 2) then 
            go = g1 + g2
        else 
            go = g1
        end if
        
        iostat(1) = 1
        yout(1) = go
        
        return
    end subroutine type398
    
    
! * SUBROUTINE:  Variable value modification
! *
! * PURPOSE:     Allow modifications to variable values for Cx
! *********************************************************************
! * INPUTS
! * ======
! *  1. p1      : Pressure type input                          (kPa)
! *  2. f1      : Flow type input                              (kg/s)
! *  3. t1      : Temperature type input                       (C)
! *  4. c1      : Control type input                           (-)
! *  5. o1      : Other type input                             (-)
! *  6. e1      : Energy type input                            (kJ)
! *  7. q1      : Power type input                             (kW)
! *  8. g1      : Humidity type input                          (kg/kg)
! *  9. c2      : control multiplier
! *     All input values are multiplied by #9. Unused outputs should be directed to category index 0.
! *     1- Output g2           swap
! *     2- Output g1 + g2      add (i.e. for gradient)
! *
! * OUTPUTS
! * =======
! *  1. p2      : Pressure type output                          (kPa)
! *  2. f2      : Flow type output                              (kg/s)
! *  3. t2      : Temperature type output                       (C)
! *  4. c2      : Control type output                           (-)
! *  5. o2      : Other type output                             (-)
! *  6. e2      : Energy type output                            (kJ)
! *  7. q2      : Power type output                             (kW)
! *  8. g2      : Humidity type output                          (kg/kg)
! *
! * PARAMETERS
! * ==========
! *  1. dummy
! *
! * SAVED
! * =====
! ***********************************************************************
!
!   MAJOR RESTRICTIONS:  NOne
!
!   DEVELOPER:           Michael A. Galler
!                        NIST
!
!   DATE:                September 15, 2024
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:
!
! **********************************************************************

        subroutine type399(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=9,no=1,np=1,&
                                             ndiffeq=0,ns=0
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         ::p1,f1,t1,c1,o1,e1,q1,g1,c2
        real         ::pout,fout,tout,cout,oout,eout,qout,gout
        
! **** Read in input
        p1 = xin(1)
        f1 = xin(2)
        t1 = xin(3)
        c1 = xin(4)
        o1 = xin(5)
        e1 = xin(6)
        q1 = xin(7)
        g1 = xin(8)
        c2 = xin(9) ! multiplier value, not a switch value

        pout = p1 * c2
        fout = f1 * c2
        tout = t1 * c2
        cout = c1 * c2
        oout = o1 * c2
        eout = e1 * c2
        qout = q1 * c2
        gout = g1 * c2

        iostat(1) = 1
        yout(1) = pout
        yout(2) = fout
        yout(3) = tout
        yout(4) = cout
        yout(5) = oout
        yout(6) = eout
        yout(7) = qout
        yout(8) = gout
        
        return
    end subroutine type399

    
    !MAG TODO- add type396, type397, type398, add all to select()
    
    !MAG TODO- add fault controller/aggregator, one control input multiple control outputs
