! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Fanger PMV and PPD
! *
! * PURPOSE:    Calculate the Predicted Mean Vote (PMV) and the
! *             Predicted Percentage Dissatisfied using the ISO version
! *             of Fanger's algorithm
! *
! **********************************************************************
! * INPUTS
! * ======
! * ta      : Dry bulb temperature                                   (C)
! * tr      : Mean radiant temperature                               (C)
! * vel     : Air speed                                            (m/s)
! * w       : Humidity ratio                                     (kg/kg)
! * patm    : Atmospheric pressure (relative to Standard Presssure)(kPa)
! *
! * OUTPUT
! * ======
! * pmv     : Predicted Mean Vote                                    (-)
! * ppd     : Predicted Percentage Dissatisfied                      (%)
! *
! * PARAMETERS
! * ==========
! * met     : Metabolic rate                           (met = 58.2 W/m2)
! * clo     : Clothing resistance                   (clo = 0.155 m2.C/W)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  None
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University
!
!   DATE:                August 13, 1996
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  None
!   FUNCTIONS  CALLED:   None
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ISO Standard 7730
!
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type461(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=5,no=2,np=2,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: met,mw,icl
        real         :: ta,tr,vel,w,patm ,clo,pa,fcl,hcf,tra,taa,tcla ,&
                        p1,p2,p3,p4,p5,xn,xf,eps,hcn,hc,tcl,hl1,hl2,hl3,&
                        hl4,hl5,hl6,ts,pmv,ppd
        integer      :: m,n

! **** Read in inputs
        ta   = xin(1)
        tr   = xin(2)
        vel  = xin(3)
        w    = xin(4)
        patm = xin(5)*1000. + 101325.0
! **** Read in parameters
        met  = par_v(1)
        m    = met*58.2
        clo  = par_v(2)
        icl  = clo*0.155
! **** Metabolic rate minus external work (external work assumed = 0)
        mw   = m
! **** Water vapor pressure
	pa   = patm*w/(0.62198+w)
! **** Clothing area factor
        if (icl<0.078) then
	    fcl  = 1.0+1.29*icl
        else
	    fcl  = 1.05+0.645*icl
        endif
! **** Convective heat transfer coefficient
	hcf  = 12.1*sqrt(vel)
! **** Absolute temperatures
	tra  = tr+273.
	taa  = ta+273.
! **** Calculate surface temperature of clothing by iteration
! **** First guess for surface temperature
	tcla = taa+(35.5-ta)/(3.5*(6.45*icl+0.1))
! **** Calculation terms
	p1   = icl*fcl
	p2   = p1*3.96
	p3   = p1*100
	p4   = p1*taa
	p5   = 308.7-0.028*mw+p2*(tra/100.)**4
	xn   = tcla/100.
	xf   = xn
	n    = 0
! **** Convergence tolerance (hundredths of Kelvins)
	eps  = 0.00015
! **** Iterate to find surface temperature
 350    xf   = (xf+xn)/2.
! **** Heat transfer coefficient by natural convection
	hcn  = 2.38*abs(100.*xf-taa)**0.25
        if (hcf>hcn) then
            hc   = hcf
        else
            hc   = hcn
        endif
	xn   = (p5+p4*hc-p2*xf**4)/(100.+p3*hc)
	n    = n+1
	if (n>150) then
            write(*,*) 'type 461: iteration for tcl failed to converge'
            goto 550
        endif
	if (abs(xn-xf)>eps) goto 350
	tcl  =  100.*xn-273.
! **** Calculate heat loss components
! **** Heat loss by diffusion through skin
        hl1  = 3.05*0.001*(5733.-6.99*mw-pa)
! **** Heat loss by sweating (comfort)
	if (mw>58.15) then
            hl2  = 0.42*(mw-58.15)
        else
	    hl2  = 0.0
        endif
! **** Latent respiration heat loss
        hl3  = 1.7*0.00001*m*(5867.-pa)
! **** Dry respiration heat loss
        hl4  = 0.0014*m*(34.-ta)
! **** Heat loss by radiation
	hl5  = 3.96*fcl*(xn**4-(tra/100.)**4)
! **** Heat loss by convection
	hl6  = fcl*hc*(tcl-ta)
! **** Calculate pmv and ppd
! **** Thermal sensation tran coeff
	ts   = 0.303*exp(-0.036*m)+0.028
! **** Predicted mean vote
	pmv  = ts*(mw-hl1-hl2-hl3-hl4-hl5-hl6)
	if (abs(pmv)>3.) then
            write(*,*) 'type 461: pmv out of range'
            goto 550
        endif
! **** Predicted percent dissatisfied
        ppd  = 100.-95.*exp(-0.03353*pmv**4-0.2179*pmv**2)
	goto 570
! **** Iteration failed to converge or pmv out of range
 550    pmv  = 999999.
	ppd  = 100.
 570    continue
! **** Assign output
        yout(1) = pmv
        yout(2) = ppd
! **** Allow freezing of algebraic variables
        iostat(1) = 1
        iostat(2) = 1
        return
        end subroutine type461


! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Heat Meter
! *
! * PURPOSE:    Calculate the thermal power delivered by a fluid stream
! *             from the flow rate and the temperature difference, and
! *             integrate the power to obtain the energy transfered
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. w       : Mass flow rate                                  (kg/s)
! *  2. ti      : Inlet temperature                                  (C)
! *  3. to      : Outlet temperature                                 (C)
! *
! * OUTPUTS
! * =======
! *  1. q       : Power                                             (kW)
! *  2. e       : Energy                                            (kJ)
! *
! * PARAMETERS
! * ==========
! *  1. ifluid  : air = 1, water = 2                                 (-)
! *
! * SAVED
! * =====
! *  1. time    : time of previous call                              (s)
! *  2. ep      : solution of differential equation from previous call
! *  3. ep      : solution of differential equation from previous step
! *               time
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  None
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                August 14, 1996
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  DIFFEQ
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

        subroutine type462(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=2,np=1,&
                                             ndiffeq=1,ns=1+ndiffeq*2
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: cpa=1.0, cpw=4.18
        real         :: w,ti,to,ep,cp,q,aa,bb,ebar,e
        integer      :: ifluid

! **** Read in inputs
        w  = xin(1)
        ti = xin(2)
        to = xin(3)
! **** Read in parameters
        ifluid = nint(par_v(1))
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Reset "heat meter"
                saved_v(2) = 0.0
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update value of "meter reading" from
! **** previous time-step
            saved_v(3) = saved_v(2)
        endif
! **** Update previous values
        ep = saved_v(3)
! **** Calculate heat rate
        if (ifluid==1) then
            cp = cpa
        elseif (ifluid==2) then
            cp = cpw
        else
            stop 'type 462: parameter 1 out of range'
        endif
        q = w*cp*(ti-to)
! **** Integrate heat rate analytically
        aa = 0.0
        bb = q
        call diffeq(time,aa,bb,ep,e,ebar)
! **** Save provisional value to be used in updating at next step time
        saved_v(2) = e
! **** Save time of current call
        saved_v(1) = time
! **** Outputs
        yout(1) = q
        yout(2) = e
! **** Disallow freezing
        iostat(1) = 0
! **** Return
        return
        end subroutine type462

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Energy Meter
! *
! * PURPOSE:    Calculate energy consumption by integrating power
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. q       : Power                                             (kW)
! *
! * OUTPUT
! * ======
! *  1. e       : Energy                                            (kJ)
! *
! * PARAMETERS
! * ==========
! *  1.         : (dummy)                                            (-)
! *
! * SAVED
! * =====
! *  1. time    : time of previous call                              (s)
! *  2. ep      : solution of differential equation from previous call
! *  3. ep      : solution of differential equation from previous step
! *               time
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  None
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                August 14, 1996
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  DIFFEQ
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

        subroutine type463(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=1,no=1,np=1,&
                                             ndiffeq=1,ns=1+ndiffeq*2
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: q,ep,aa,bb,ebar,e

! **** Read in inputs
        q  = xin(1)
! **** (no parameters)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
! **** Reset "energy meter"
                saved_v(2) = 0.0
            endif
        endif
        if (time>saved_v(1)) then
! **** First call of timestep - update value of "meter reading" from
! **** previous time-step
            saved_v(3) = saved_v(2)
        endif
! **** Update previous values
        ep = saved_v(3)
! **** Integrate power analytically
        aa = 0.0
        bb = q
        call diffeq(time,aa,bb,ep,e,ebar)
! **** Save provisional value to be used in updating at next step time
        saved_v(2) = e
! **** Save time of current call
        saved_v(1) = time
! **** Outputs
        yout(1) = e
! **** Disallow freezing
        iostat(1) = 0
! **** Return
        return
        end subroutine type463
        
