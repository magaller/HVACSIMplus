! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Fluid resistance
! *
! * PURPOSE:    Calculate the inlet pressure for a fluid resistance
! *             given the outlet pressure and the inlet flow rate.
! *             Uses a linear relationship between pressure drop and
! *             flow rate at low flow rates to avoid numerical problems
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. w       : Mass flow rate                                  (kg/s)
! *  2. po      : Outlet pressure                                  (kPa)
! *
! * OUTPUT
! * ======
! *  1. pi      : Inlet pressure                                   (kPa)
! *
! * PARAMETERS
! * ==========
! *  1. r       : Resistance                                 (0.001/k.m)
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
!   FUNCTIONS  CALLED:   DPQUDLIN
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

        subroutine type341(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=2,no=1,np=1,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: w,po,r,pi,dpqudlin

! **** Read in inputs
        w=xin(1)
        po=xin(2)
! **** Read in parameters
        r=par_v(1)
! **** Call pressure drop function
        pi=po+dpqudlin(r,w)
! **** Assign output
        yout(1)=pi
! **** Allow freezing of algebraic variable
        iostat(1)=1
        return
        end subroutine type341

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Fluid resistance - calculates flow rate
! *
! * PURPOSE:    Calculates the flow rate for a fluid resistance
! *             given the inlet pressure and the outlet pressure.
! *             Uses a linear relationship between pressure drop and
! *             flow rate at low flow rates to avoid numerical problems
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. pi      : Inlet pressure                                   (kPa)
! *  2. po      : Outlet pressure                                  (kPa)
! *
! * OUTPUT
! * ======
! *  1. w       : Mass flow rate                                  (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. r       : Resistance                                 (0.001/k.m)
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
!   FUNCTIONS  CALLED:   DPQUDLIN
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

        subroutine type342(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=2,no=1,np=1,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: pi,po,r,dp,w,wqudlin

! **** Read in inputs
        pi=xin(1)
        po=xin(2)
! **** Read in parameters
        r=par_v(1)
! **** Calculate pressure drop and call flow rate function
        dp=pi-po
        w=wqudlin(r,dp)
! **** Assign output
        yout(1)=w
! **** Allow freezing of algebraic variable
        iostat(1)=1
        return
        end subroutine type342

! *********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Flow split
! *
! * PURPOSE:    Calculate the outlet flow rates and the inlet pressure
! *             for a flow split, given the outlet pressures and the
! *             inlet flow rate
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. wi      : Inlet flow rate                                 (kg/s)
! *  2. po1     : Pressure at outlet 1                             (kPa)
! *  3. po2     : Pressure at outlet 2                             (kPa)
! *
! * OUTPUTS
! * =======
! *  1. wo1     : Mass flow rate at outlet 1                      (kg/s)
! *  2. wo2     : Mass flow rate at outlet 2                      (kg/s)
! *  3. pi      : Inlet pressure                                   (kPa)
! *
! * PARAMETERS
! * ==========
! *  1. ri      : Inlet resistance                           (0.001/k.m)
! *  2. ro1     : Resistance of outlet 1                     (0.001/k.m)
! *  3. ro2     : Resistance of outlet 2                     (0.001/k.m)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  Resistances must not be negative. At least one
!                        outlet resistance must be non-zero.
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 8, 1995
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  FLOWSPLT
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
! * wcritu  : flow threshold, lower limit of pure quadratic solution
! * wcritl  : flow threshold, upper limit of pure linear solution
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type345(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=3,np=3,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: wcritl=3.e-2, wcritu=6.e-2
        real         :: wi,po1,po2,ri,ro1,ro2,wo2,wo1,pi
        integer      :: i,ifail

! **** Read in and check inputs
        wi=xin(1)
        po1=xin(2)
        po2=xin(3)
! **** Check for reverse inlet flow at last iteration
        if (wi<-0.01 .and. lastit==1)&
            write(*,*) 'warning:', ' type 345: ',&
                       'reverse inlet flow'
! **** Read in parameters
        ri=par_v(1)
        ro1=par_v(2)
        ro2=par_v(3)
! **** Call flow split routine
        call flowsplt(wi,po1,po2,ri,ro1,ro2,wcritl,wcritu,rtolx,&
                      pi,wo1,wo2,ifail)
! **** Check for unsuccesful completion
        if (ifail==1) then
! **** One or more resistances negative
            write(*,*) 'Warning:',' type 345: ',&
                       'resistances must not be negative'
            stop
        elseif (ifail==2) then
! **** Zero resistance for both outlet branches
            write(*,*) 'Warning:',' type 345: ',&
                       'both outlet resistances cannot be zero'
            stop
        endif
! **** Assign outputs
        yout(1)=wo1
        yout(2)=wo2
        yout(3)=pi
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo

        return
        end subroutine type345

! *********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Asymmetric flow split
! *
! * PURPOSE:    Calculate the flow rate in the main outlet and the
! *             pressure at the branch outlet of a flow split, given the
! *             inlet and branch flow rates and the main outlet pressure.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. wi      : Inlet flow rate                                 (kg/s)
! *  2. pomain  : Pressure at main outlet                          (kPa)
! *  3. wobranch: Mass flow rate through branch outlet             (kPa)
! *
! * OUTPUTS
! * =======
! *  1. pi      : Inlet pressure                                   (kPa)
! *  2. womain  : Mass flow rate through main outlet              (kg/s)
! *  3. pobranch: Pressure at branch outlet                       (kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. ri      : Inlet resistance                           (0.001/k.m)
! *  2. romain  : Resistance of main outlet                  (0.001/k.m)
! *  3. robranch: Resistance of branch outlet                (0.001/k.m)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  Resistances must not be negative. At least one
!                        outlet resistance must be non-zero.
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 8, 1995
!
!   INCLUDE FILES:       None
!   SUBROUTINES CALLED:  FLOWSPLT
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
! * wcritu  : flow threshold, lower limit of pure quadratic solution
! * wcritl  : flow threshold, upper limit of pure linear solution
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type346(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=3,np=3,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: wi,pomain,wobranch,ri,romain,robranch,womain,&
                        pc,dpqudlin,pi,pobranch
        integer      :: i

! **** Read in and check inputs
        wi       = xin(1)
        pomain   = xin(2)
        wobranch = xin(3)
! **** Check for reverse inlet flow at last iteration
        if (wi<-0.01 .and. lastit==1) then
            write(*,*) 'warning: ',' type 346: ',&
                       'reverse inlet flow'
        endif
! **** Read in parameters
        ri       = par_v(1)
        romain   = par_v(2)
        robranch = par_v(3)
! **** Calculate main outlet flow and inlet and branch outlet pressures
! **** Calculate center node pressure
        womain   = wi-wobranch
        pc       = pomain+dpqudlin(romain,womain)
! **** Calculate inlet and branch outlet pressures
        pi       = pc+dpqudlin(ri,wi)
        pobranch = pc-dpqudlin(robranch,wobranch)
! **** Assign outputs
        yout(1) = pi
        yout(2) = womain
        yout(3) = pobranch
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i) = 1
        enddo

        return
        end subroutine type346

! *********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Flow merge
! *
! * PURPOSE:     Calculate the outlet flow rate and the inlet pressures
! *              for a flow merge, given the outlet pressure and the
! *              inlet flow rates
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. wi1     : Inlet flow rate 1                               (kg/s)
! *  2. wi2     : Inlet flow rate 2                               (kg/s)
! *  3. po      : Pressure at outlet                               (kPa)
! *
! * OUTPUTS
! * =======
! *  1. wo      : Mass flow rate at outlet                        (kg/s)
! *  2. pi1     : Inlet pressure 1                                 (kPa)
! *  3. pi2     : Inlet pressure 2                                 (kPa)
! *
! * PARAMETERS
! * ==========
! *  1. ri1     : Inlet resistance 1                         (0.001/k.m)
! *  2. ri2     : Inlet resistance 2                         (0.001/k.m)
! *  3. ro      : Resistance of outlet                       (0.001/k.m)
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
!   SUBROUTINES CALLED:  FLOWMERG
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

        subroutine type348(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=3,np=3,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: wi1,wi2,po,ri1,ri2,ro,wo,pi2,pi1
        integer      :: i

! **** Read in inputs
        wi1=xin(1)
        wi2=xin(2)
        po=xin(3)
! **** Read in parameters
        ri1=par_v(1)
        ri2=par_v(2)
        ro=par_v(3)
! **** Call flow merge routine
        call flowmerg(wi1,wi2,po,ri1,ri2,ro,pi1,pi2,wo)
! **** Assign outputs
        yout(1)=wo
        yout(2)=pi1
        yout(3)=pi2
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo
        return
        end subroutine type348

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! *
! * SUBROUTINE: Room air mass balance
! *
! * PURPOSE:    Performs a mass balance on six air streams, the supply
! *             air from the HVAC system, the return air to the HVAC
! *             system, leakage to ambient, local mechanical extract
! *             to ambient and flow to two adjacent zones. The flow from
! *             the first adjacent zone is calculated from the pressure
! *             difference and the resistance between that zone and the
! *             zone on which the mass balance is performed, and the
! *             leakage is calculated similarly.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. wsup    : supply air mass flow rate (positive in)         (kg/s)
! *  2. padjzon1: pressure of 1st adjacent zone                    (kPa)
! *  3. wadjzon2: flow rate to 2nd adjacent zone (positive out)   (kg/s)
! *  4. proom   : room pressure                                    (kPa)
! *  5. pamb    : ambient pressure                                 (kPa)
! *
! * OUTPUT
! * ======
! *  1. wret    : return air mass flow rate (positive out)        (kg/s)
! *  2. wadjzon1: flow rate from 1st adjacent zone (positive in)  (kg/s)
! *  3. wleak   : leakage mass flow rate to ambient (positive out)(kg/s)
! *
! * PARAMETERS
! * ==========
! *  1. radjzon1: resistance to 1st adjacent zone            (0.001/k.m)
! *  2. rloss   : leakage resistance                         (0.001/k.m)
! *  3. wfan    : local extract fan mass flow rate (positive out) (kg/s)
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

        subroutine type349(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=5,no=3,np=3,ns=0
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real         :: wsup,padjzon1,wadjzon2,proom,pamb,radjzon1,&
                        rloss,wfan,dp,wadjzon1,wqudlin,wleak,wret
        integer      :: i

! **** Read in inputs
        wsup = xin(1)
        padjzon1 = xin(2)
        wadjzon2 = -xin(3)
        proom = xin(4)
        pamb  = xin(5)
! **** Read in parameters
        radjzon1 = par_v(1)
        rloss = par_v(2)
        wfan  = par_v(3)
! **** Calculate flow rate from 1st adjacent zone
        dp = padjzon1-proom
        wadjzon1=wqudlin(radjzon1,dp)
! **** Calculate leakage flow rate to ambient
        dp = proom-pamb
        wleak=wqudlin(rloss,dp)
! **** Mass balance on room to find return air flow rate
        wret=wsup+wadjzon1+wadjzon2-wfan-wleak
! **** Assign outputs
        yout(1) = wret
        yout(2) = wadjzon1
        yout(3) = wleak
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo
        return
        end subroutine type349
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     Fan or pump
! *
! * PURPOSE:        Calculates inlet pressure and fluid heating rate from
! *                 outlet pressure, mass flow rate, rotation speed and
! *                 inlet temperature
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. w       : mass flow rate                                  (kg/s)
! *  2. po      : outlet pressure                                  (kPa)
! *  3. n       : rotation speed                                 (rev/s)
! *
! * OUTPUTS
! * =======
! *  1. pi      : inlet pressure                                   (kPa)
! *  2. qa      : fluid stream heat addition rate                   (kW)
! *  3. power   : power consumption                                 (kW)
! *
! * PARAMETERS
! * ==========
! *  1.         : coefficient of PN**0 in normalized head curve
! *  2.         : coefficient of PN**1 in normalized head curve
! *  ...
! *  5.         : coefficient of PN**4 in normalized head curve
! *  6.         : coefficient of PN**0 in normalized head curve
! *  7.         : coefficient of PN**1 in normalized head curve
! * ...
! * 10.         : coefficient of PN**4 in normalized head curve
! * 11. d       : diameter                                           (-)
! * 12. mode    : fluid: 1 = air, any other value = water            (-)
! * 13. wncritlo: normalised flow at lower bound of valid region     (-)
! * 14. wncritup: normalised flow at upper bound of valid region     (-)
! *
! * SAVED
! * =====
! *  1.         : flag used to indicate initialization has been performed
! *               (used in HVACSIM+ only)
! *  2. alo/n**2: A coefficient in lower region, divided by N squared
! *  3. rlo     : equivalent resistance in lower region
! *  4. aup/n**2: A coefficient in upper region, divided by N squared
! *  5. rup     : equivalent resistance in upper region
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  Fixed blade angle, no inlet guide vanes, all
!                        inefficiencies appear as heat in the fluid
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                September 30, 1992
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
! * funp()  : normalized head curve function
! * dfunp() : derivative of normalized head curve function wrt to flow
! * funn()  : normalized efficiency function
! * cp      : specific heat of working fluid
! * rho     : density of working fluid
! * wn      : normalized mass flow rate
! * dwn     : fraction of valid range of normalized flow rate
! * anlo    : normalized A coefficient in lower region
! * rnlo    : normalized effective resistance in lower region
! * anup    : normalized A coefficient in upper region
! * rnup    : normalized effective resistance in upper region
! * roff    : effective resistance when fan/pump is off
! * pdrop   : pressure drop when fan/pump switched off
! * alo     : (un-normalized) A coefficient in lower region
! * rlo     : (un-normalized) effective resistance in lower region
! * aup     : (un-normalized) A coefficient in upper region
! * rup     : (un-normalized) effective resistance in upper region
! * prise   : pressure rise when fan/pump switched on
! * pn      : normalized pressure rise
! * eff     : efficiency
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type350(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=3,np=14,ns=5
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: n
        real         :: rhoa=1.2, rhow=1000., cpa=1.0, cpw=4.18
        real         :: z,funp,dfunp,funn,w,po,d,mode,wncritlo,wncritup,&
                        cp,rho,dwn,wn,anlo,rnlo,anup,rnup,roff,pdrop,&
                        dpqudlin,pi,power,alo,rlo,prise,eff,aup,rup,pn,qa
        integer      :: i

! **** Internal functions
! **** "Head curve" - normalized pressure rise as a function of normalized flow
        funp(z)  = par_v(1)+z*(par_v(2)+z*(par_v(3)+z*(par_v(4)+z*par_v(5))))
! **** Derivative of head curve
        dfunp(z) = par_v(2)+z*(2.0*par_v(3)+z*(3.0*par_v(4)+z*4.0*par_v(5)))
! **** Efficiency curve
        funn(z)  = par_v(6)+z*(par_v(7)+z*(par_v(8)+z*(par_v(9)+z*par_v(10))))
! **** Read in inputs
        w        = xin(1)
        po       = xin(2)
        n        = xin(3)
! **** Read in parameters not used in head curve and efficiency functions
        d        = par_v(11)
        mode     = nint(par_v(12))
        wncritlo = par_v(13)
        wncritup = par_v(14)
! **** Outside the range of validity of the head curve polynomial, an ideal
! **** fan/pump model is used (pn=an-rn*wn*|wn|, where pn is the normalised
! **** pressure rise, wn is the normalised flow rate, an is +ve constant and
! **** rn is a +ve constant). an and rn are calculated so that the head curve
! **** and its gradient are continuous.

! **** Set specific heat and density for appropriate fluid
        if (mode==1) then
            cp  = cpa
            rho = rhoa
        else
            cp  = cpw
            rho = rhow
        endif
! *** Limit rotation speed to positive values
        n = max(0.0,n)
! **** "One-off" calculations
! **** Test for first call of simulation, assuming saved array is
! **** initialized to zero
        if (init==0 .and. saved_v(1)==0.0) then
            saved_v(1) = 1.
! **** Check head curve is monotonically decreasing
! **** Test gradient at 50 points across user-defined range
            dwn = (wncritup-wncritlo)/50.
            wn  = wncritlo
            do i=0,50
                if (dfunp(wn)>0.0) then
                    write(*,*)&
                    'type 350: head curve slope positive at wn =',wn
                endif
                wn = wn+dwn
            enddo
! **** Determine coefficents of quadratic fits outside valid region
            anlo = funp(wncritlo)-wncritlo*dfunp(wncritlo)/2.0
            if (anlo<=0.0) stop 'type350: anlo<=0'
            rnlo = -dfunp(wncritlo)/(2.0*abs(wncritlo))
            if (rnlo<=0.0) stop 'type350: rnlo<=0'
            anup = funp(wncritup)-wncritup*dfunp(wncritup)/2.0
            if (anup<=0.0) stop 'type350: anup<=0'
            rnup = -dfunp(wncritup)/(2.0*abs(wncritup))
            if (rnup<=0.0) stop 'type350: rnup<=0'
! **** Store partly or completely un-normalized coefficients
            saved_v(2) = anlo*0.001*rho*d*d
            saved_v(3) = rnlo*0.001/(rho*d*d*d*d)
            saved_v(4) = anup*0.001*rho*d*d
            saved_v(5) = rnup*0.001/(rho*d*d*d*d)
        endif
! **** Calculate pressure rise, fluid heating rate and power consumption
        if (n<=0.0) then
! **** Fan/pump off
! **** Calculate the pressure drop assuming the resistance is equal to
! **** the equivalent resistance in the upper region
            roff  = saved_v(5)
            pdrop = dpqudlin(roff,w)
            pi    = po+pdrop
            power = 0.0
        else
! **** Fan/pump on
! **** Normalized flow rate
            wn = w/(rho*d*d*d*n)
            if (wn<wncritlo) then
! **** Below valid region
! **** Calculate pressure rise from quadratic extrapolation of head curve
                alo   = saved_v(2)*n*n
                rlo   = saved_v(3)
                prise = alo-dpqudlin(rlo,w)
                pi    = po-prise
                eff   = funn(wncritlo)
            elseif (wn>wncritup) then
! **** Above valid region
! **** Calculate pressure rise from quadratic extrapolation of head curve
                aup   = saved_v(4)*n*n
                rup   = saved_v(5)
                prise = aup-dpqudlin(rup,w)
                pi    = po-prise
                eff   = funn(wncritup)
            else
! **** In valid range
! **** Calculate pressure rise from head curve
                pn    = funp(wn)
                prise = 0.001*pn*rho*d*d*n*n
                pi    = po-prise
                eff   = funn(wn)
            endif
! **** Oower consumption - only correct within range of validity of
! **** efficiency polynomial
            if (w>0.0.and.(po-pi)>0.0) then
                power = (po-pi)*w/(rho*eff)
            else
                power = 0.0
            endif
        endif
! **** Rate at which heat is added to fluid
        if (mode==1) then
! **** Air - include effect of fluid work
            qa = power
        else
! **** Water - exclude effect of fluid work
            qa = (po-pi)*w*(1./eff-1.)/rho
        endif
! **** Output
        yout(1) = pi
        yout(2) = qa
        yout(3) = power
! **** Allow freezing
        iostat(1) = 1
        iostat(2) = 1
        iostat(3) = 1
! **** Return
        return
        end subroutine type350
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     Fan or pump (implicit flow)
! *
! * PURPOSE:        Calculates pressure, rotation speed and fluid
! *                 heating rate from mass flow rate using a quadratic
! *                 load line. Uses external iteration on rotation speed.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. ifstatus: fan status (0 = off, 1 = on)                       (-)
! *  2. w       : mass flow rate                                  (kg/s)
! *  3. n       : rotation speed                                 (rev/s)
! *
! * OUTPUTS
! * =======
! *  1. dp      : pressure rise                                    (kPa)
! *  2. qa      : fluid stream heat addition rate                   (kW)
! *  3. power   : power consumption                                 (kW)
! *  4. n       : rotation speed                                 (rev/s)
! *
! * PARAMETERS
! * ==========
! *  1.         : coefficient of PN**0 in normalized head curve
! *  2.         : coefficient of PN**1 in normalized head curve
! *  ...
! *  5.         : coefficient of PN**4 in normalized head curve
! *  6.         : coefficient of PN**0 in normalized head curve
! *  7.         : coefficient of PN**1 in normalized head curve
! * ...
! * 10.         : coefficient of PN**4 in normalized head curve
! * 11. d       : diameter                                           (-)
! * 12. mode    : fluid: 1 = air, any other value + water            (-)
! * 13. wncritlo: normalised flow at lower bound of valid region     (-)
! * 14. wncritup: normalised flow at upper bound of valid region     (-)
! * 15. pstatsp : setpoint for the fan static pressure control loop(kPa)
! * 16. dsarea  : cross sect area of duct at static pressure sensor (m2)
! * 17. rduct   : resistance of duct sys. upstream of sensor (0.001/k.m)
! *
! * SAVED
! * =====
! *  1.         : flag used to indicate initialization has been performed
! *               (used in HVACSIM+ only)
! *  2. alo/n**2: A coefficient in lower region, divided by N squared
! *  3. rlo     : equivalent resistance in lower region
! *  4. aup/n**2: A coefficient in upper region, divided by N squared
! *  5. rup     : equivalent resistance in upper region
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  Fixed blade angle, no inlet guide vanes, all
!                        inefficiencies appear as heat in the fluid
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                September 30, 1992
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
!
!  Calculates actual pressure rise from simple model of external system:
!
!  DP = PTOTSP + RDUCT * W * ABS(W)
!
!  where PTOTSP is the total pressure corresponding to the set-point for
!  the fan static pressure control loop and RDUCT is the resistance of
!  the duct/pipe system upstream of the static pressure sensor.
!
!  The actual pressure rise is compared with the normalised pressure rise
!  calculated from the head curve in order to calculate the rotation speed
!  of the fan, hence its efficiency and hence the power consumption and
!  temperature rise. The rotation speed is both an input and an output.
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! * funp()  : normalized head curve function
! * dfunp() : derivative of normalized head curve function wrt to flow
! * funn()  : normalized efficiency function
! * cp      : specific heat of working fluid
! * rho     : density of working fluid
! * wn      : normalized mass flow rate
! * dwn     : increment of valid range of normalized flow rate
! * anlo    : normalized A coefficient in lower region
! * rnlo    : normalized effective resistance in lower region
! * anup    : normalized A coefficient in upper region
! * rnup    : normalized effective resistance in upper region
! * roff    : effective resistance when fan/pump is off
! * pdrop   : pressure drop when fan/pump switched off
! * alo     : (un-normalized) A coefficient in lower region
! * alodnn  : A coefficient in lower region normalized by fan speed only
! * rlo     : (un-normalized) effective resistance in lower region
! * aup     : (un-normalized) A coefficient in upper region
! * aupdnn  : A coefficient in upper region normalized by fan speed only
! * rup     : (un-normalized) effective resistance in upper region
! * prise   : pressure rise when fan/pump switched on
! * pn      : normalized pressure rise
! * eff     : efficiency
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type351(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=3,no=4,np=17,ns=5
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: n
        real         :: rhoa=1.2, rhow=1000., cpa=1.0, cpw=4.18
        real         :: z,funp,dfunp,funn,w,d,wncritlo,wncritup,pstatsp,&
                        dsarea,rduct,cp,rho,dwn,wn,anlo,rnlo,anup,rnup,&
                        roff,pdrop,dpqudlin,prise,power,ptotsp,alodnn,&
                        rlo,alo,eff,aupdnn,rup,aup,pn,qa
        integer      :: i,ifstatus,mode

! **** Internal functions
! **** "Head curve" - normalized pressure rise as a function of normalized flow
        funp(z) = par_v(1)+z*(par_v(2)+z*(par_v(3)+z*(par_v(4)+z*par_v(5))))
! **** Derivative of head curve
        dfunp(z) = par_v(2)+z*(2.0*par_v(3)+z*(3.0*par_v(4)+z*4.0*par_v(5)))
! **** Efficiency curve
        funn(z) = par_v(6)+z*(par_v(7)+z*(par_v(8)+z*(par_v(9)+z*par_v(10))))
! **** Read in inputs
        ifstatus = nint(xin(1))
        w        = xin(2)
        n        = xin(3)
! **** Read in parameters not used in head curve and efficiency functions
        d        = par_v(11)
        mode     = nint(par_v(12))
        wncritlo = par_v(13)
        wncritup = par_v(14)
        pstatsp  = par_v(15)
        dsarea   = par_v(16)
        rduct    = par_v(17)
! **** Outside the range of validity of the head curve polynomial, an ideal
! **** fan/pump model is used (pn=an-rn*wn*|wn|, where pn is the normalised
! **** pressure rise, wn is the normalised flow rate, an is +ve constant and
! **** rn is a +ve constant). an and rn are calculated so that the head curve
! **** and its gradient are continuous.

! **** Set specific heat and density for appropriate fluid
        if (mode==1) then
            cp=cpa
            rho=rhoa
        else
            cp=cpw
            rho=rhow
        endif
! *** Limit rotation speed to positive values
        n=max(0.0,n)
! **** "One-off" calculations
! **** Test for first call of simulation, assuming saved array is
! **** initialized to zero
        if (init==0 .and. saved_v(1)==0.0) then
            saved_v(1) = 1.0
! **** Check head curve is monotonically decreasing
! **** Test gradient at 50 points across user-defined range
            dwn = (wncritup-wncritlo)/50.
            wn  = wncritlo
            do i=0,50
                if (dfunp(wn)>0.0) then
                    write(*,*)&
                    'type 351: head curve slope positive at wn =',wn
                endif
                wn = wn+dwn
            enddo
! **** Determine coefficents of quadratic fits outside valid region
            anlo = funp(wncritlo)-wncritlo*dfunp(wncritlo)/2.0
            if (anlo<=0.0) stop 'type351: anlo<=0'
            rnlo = -dfunp(wncritlo)/(2.0*abs(wncritlo))
            if (rnlo<=0.0) stop 'type351: rnlo<=0'
            anup = funp(wncritup)-wncritup*dfunp(wncritup)/2.0
            if (anup<=0.0) stop 'type351: anup<=0'
            rnup = -dfunp(wncritup)/(2.0*abs(wncritup))
            if (rnup<=0.0) stop 'type351: rnup<=0'
! **** Store partly or completely un-normalized coefficients
            saved_v(2) = anlo*0.001*rho*d*d
            saved_v(3) = rnlo*0.001/(rho*d*d*d*d)
            saved_v(4) = anup*0.001*rho*d*d
            saved_v(5) = rnup*0.001/(rho*d*d*d*d)
        endif
! **** Calculate pressure rise and temperature rise
        if (ifstatus==0) then
! **** Fan/pump off
! **** Calculate the pressure drop assuming the resistance is equal to
! **** the equivalent resistance in the upper region
            roff  = saved_v(5)
            pdrop = dpqudlin(roff,w)
            prise = -pdrop
            power = 0.0
            n     = 0.0
        else
! **** Fan/pump on
! **** Normalized flow rate
            if (w<=0.0) then
! **** Flow is zero, so normalized flow is zero
                wn = 0.0
            elseif (n==0.0) then
! **** Rotation speed guess is zero - normalized flow rate is infinite
! **** or indeterminate, use 1.0 as a better starting guess
                wn = 1.0
            else
                wn = w/(rho*d*d*d*n)
           endif
! **** Total pressure rise from load line is sum of total pressure at
! **** static pressure sensor and pressure drop across fixed resistance
! **** upstream of sensor
            ptotsp = pstatsp+0.001*0.5*w*w/(rho*dsarea*dsarea)
            prise  = ptotsp+dpqudlin(rduct,w)
            if (prise<0.0) then
                stop 'type 351: negative pressure rise across fan'
            endif
            if (wn<wncritlo) then
! **** Below valid region
! **** Calculate rotation speed by comparing un-normalized and normalized
! **** estimates of the constant term in the quadratic extrapolation of
! **** the head curve
                alodnn = saved_v(2)
                rlo    = saved_v(3)
                alo    = prise+dpqudlin(rlo,w)
                n      = sqrt(alo/alodnn)
                eff    = funn(wncritlo)
            elseif (wn>wncritup) then
! **** Above valid region
! **** Calculate rotation speed by comparing un-normalized and normalized
! **** estimates of the constant term in the quadratic extrapolation of
! **** the head curve
                aupdnn = saved_v(4)
                rup    = saved_v(5)
                aup    = prise+dpqudlin(rup,w)
                n      = sqrt(aup/aupdnn)
                eff    = funn(wncritup)
            else
! **** In valid range
! **** Calculate rotation speed by comparing un-normalized and normalized
! **** estimates of the pressure rise
                pn     = funp(wn)
                n      = sqrt(prise/(0.001*pn*rho*d*d))
                eff    = funn(wn)
            endif
! **** Power consumption
            if (w>0.0.and.(prise)>0.0) then
                power  = (prise)*w/(rho*eff)
            else
                power  = 0.0
            endif
        endif
! **** Rate at which heat is added to fluid
        if (mode==1) then
! **** Air - include effect of fluid work
            qa = power
        else
! **** Water - exclude effect of fluid work
            qa = max(0.0,(prise)*w*(1./eff-1.)/rho)
        endif
! **** Output
        yout(1) = prise
        yout(2) = qa
        yout(3) = power
        yout(4) = n
! **** Allow freezing
        do i=1,no
            iostat(i) = 1
        enddo
! **** Return
        return
        end subroutine type351
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     Fan or pump - temperature rise
! *
! * PURPOSE:        Calculates inlet pressure, fluid heating rate and
! *                 outlet temperature from outlet pressure, mass flow
! *                 rate, rotation speed and inlet temperature
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. w       : mass flow rate                                  (kg/s)
! *  2. po      : outlet pressure                                  (kPa)
! *  3. n       : rotation speed                                 (rev/s)
! *  4. ti      : inlet temperature                                  (C)
! *
! * OUTPUTS
! * =======
! *  1. pi      : inlet pressure                                   (kPa)
! *  2. qa      : fluid stream heat addition rate                   (kW)
! *  3. power   : power consumption                                 (kW)
! *  4. to      : outlet temperature                                 (C)
! *
! * PARAMETERS
! * ==========
! *  1.         : coefficient of PN**0 in normalized head curve
! *  2.         : coefficient of PN**1 in normalized head curve
! *  ...
! *  5.         : coefficient of PN**4 in normalized head curve
! *  6.         : coefficient of PN**0 in normalized head curve
! *  7.         : coefficient of PN**1 in normalized head curve
! * ...
! * 10.         : coefficient of PN**4 in normalized head curve
! * 11. d       : diameter                                           (-)
! * 12. mode    : fluid: 1 = air, any other value = water            (-)
! * 13. wncritlo: normalised flow at lower bound of valid region     (-)
! * 14. wncritup: normalised flow at upper bound of valid region     (-)
! *
! * SAVED
! * =====
! *  1.         : flag used to indicate initialization has been performed
! *               (used in HVACSIM+ only)
! *  2. alo/n**2: A coefficient in lower region, divided by N squared
! *  3. rlo     : equivalent resistance in lower region
! *  4. aup/n**2: A coefficient in upper region, divided by N squared
! *  5. rup     : equivalent resistance in upper region
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  Fixed blade angle, no inlet guide vanes, all
!                        inefficiencies appear as heat in the fluid
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                September 30, 1992
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
! * funp()  : normalized head curve function
! * dfunp() : derivative of normalized head curve function wrt to flow
! * funn()  : normalized efficiency function
! * cp      : specific heat of working fluid
! * rho     : density of working fluid
! * wn      : normalized mass flow rate
! * dwn     : fraction of valid range of normalized flow rate
! * anlo    : normalized A coefficient in lower region
! * rnlo    : normalized effective resistance in lower region
! * anup    : normalized A coefficient in upper region
! * rnup    : normalized effective resistance in upper region
! * roff    : effective resistance when fan/pump is off
! * pdrop   : pressure drop when fan/pump switched off
! * alo     : (un-normalized) A coefficient in lower region
! * rlo     : (un-normalized) effective resistance in lower region
! * aup     : (un-normalized) A coefficient in upper region
! * rup     : (un-normalized) effective resistance in upper region
! * prise   : pressure rise when fan/pump switched on
! * pn      : normalized pressure rise
! * eff     : efficiency
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type352(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=4,no=4,np=14,ns=5
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: n
        real         :: rhoa=1.2, rhow=1000., cpa=1.0, cpw=4.18
        real         :: z,funp,dfunp,funn,w,po,ti,d,wncritlo,wncritup,&
                        cp,rho,dwn,wn,anlo,rnlo,anup,rnup,roff,pdrop,&
                        dpqudlin,pi,power,alo,rlo,prise,eff,aup,rup,&
                        pn,qa,to
        integer      :: i,mode

! **** Internal functions
! **** "Head curve" - normalized pressure rise as a function of normalized flow
        funp(z)  = par_v(1)+z*(par_v(2)+z*(par_v(3)+z*(par_v(4)+z*par_v(5))))
! **** Derivative of head curve
        dfunp(z) = par_v(2)+z*(2.0*par_v(3)+z*(3.0*par_v(4)+z*4.0*par_v(5)))
! **** Efficiency curve
        funn(z)  = par_v(6)+z*(par_v(7)+z*(par_v(8)+z*(par_v(9)+z*par_v(10))))
! **** Read in inputs
        w        = xin(1)
        po       = xin(2)
        n        = xin(3)
        ti       = xin(4)
! **** Read in parameters not used in head curve and efficiency functions
        d        = par_v(11)
        mode     = nint(par_v(12))
        wncritlo = par_v(13)
        wncritup = par_v(14)
! **** Outside the range of validity of the head curve polynomial, an ideal
! **** fan/pump model is used (pn=an-rn*wn*|wn|, where pn is the normalised
! **** pressure rise, wn is the normalised flow rate, an is +ve constant and
! **** rn is a +ve constant). an and rn are calculated so that the head curve
! **** and its gradient are continuous.

! **** Set specific heat and density for appropriate fluid
        if (mode==1) then
            cp  = cpa
            rho = rhoa
        else
            cp  = cpw
            rho = rhow
        endif
! *** Limit rotation speed to positive values
        n = max(0.0,n)
! **** "One-off" calculations
! **** Test for first call of simulation, assuming saved array is
! **** initialized to zero
        if (init==0 .and. saved_v(1)==0.0) then
            saved_v(1) = 1.
! **** Check head curve is monotonically decreasing
! **** Test gradient at 50 points across user-defined range
            dwn = (wncritup-wncritlo)/50.
            wn  = wncritlo
            do i=0,50
                if (dfunp(wn)>0.0) then
                    write(*,*) &
                    'type 352: head curve slope positive at wn =',wn
                endif
                wn = wn+dwn
            enddo
! **** Determine coefficents of quadratic fits outside valid region
            anlo = funp(wncritlo)-wncritlo*dfunp(wncritlo)/2.0
            if (anlo<=0.0) stop 'type352: anlo<=0'
            rnlo = -dfunp(wncritlo)/(2.0*abs(wncritlo))
            if (rnlo<=0.0) stop 'type352: rnlo<=0'
            anup = funp(wncritup)-wncritup*dfunp(wncritup)/2.0
            if (anup<=0.0) stop 'type352: anup<=0'
            rnup = -dfunp(wncritup)/(2.0*abs(wncritup))
            if (rnup<=0.0) stop 'type352: rnup<=0'
! **** Store partly or completely un-normalized coefficients
            saved_v(2) = anlo*0.001*rho*d*d
            saved_v(3) = rnlo*0.001/(rho*d*d*d*d)
            saved_v(4) = anup*0.001*rho*d*d
            saved_v(5) = rnup*0.001/(rho*d*d*d*d)
        endif
! **** Calculate pressure rise, fluid heating rate and power consumption
        if (n<=0.0) then
! **** Fan/pump off
! **** Calculate the pressure drop assuming the resistance is equal to
! **** the equivalent resistance in the upper region
            roff  = saved_v(5)
            pdrop = dpqudlin(roff,w)
            pi    = po+pdrop
            power = 0.0
        else
! **** Fan/pump on
! **** Normalized flow rate
            wn = w/(rho*d*d*d*n)
            if (wn<wncritlo) then
! **** Below valid region
! **** Calculate pressure rise from quadratic extrapolation of head curve
                alo   = saved_v(2)*n*n
                rlo   = saved_v(3)
                prise = alo-dpqudlin(rlo,w)
                pi    = po-prise
                eff   = funn(wncritlo)
            elseif (wn>wncritup) then
! **** Above valid region
! **** Calculate pressure rise from quadratic extrapolation of head curve
                aup   = saved_v(4)*n*n
                rup   = saved_v(5)
                prise = aup-dpqudlin(rup,w)
                pi    = po-prise
                eff   = funn(wncritup)
            else
! **** In valid range
! **** Calculate pressure rise from head curve
                pn    = funp(wn)
                prise = 0.001*pn*rho*d*d*n*n
                pi    = po-prise
                eff   = funn(wn)
            endif
! **** Power consumption - only correct within range of validity of
! **** Efficiency polynomial
            if (w>0.0.and.(po-pi)>0.0) then
                power = (po-pi)*w/(rho*eff)
            else
                power = 0.0
            endif
        endif
! **** Rate at which heat is added to fluid
        if (mode==1) then
! **** Air - include effect of fluid work
            qa = power
        else
! **** Water - exclude effect of fluid work
            qa = (po-pi)*w*(1./eff-1.)/rho
        endif
! **** Temperature rise (no dynamics)
        if (w>0.0) then
            to = ti + qa/(w*cp)
        else
            to = ti
        endif
! **** Output
        yout(1) = pi
        yout(2) = qa
        yout(3) = power
        yout(4) = to
! **** Allow freezing
        iostat(1) = 1
        iostat(2) = 1
        iostat(3) = 1
! **** Return
        return
        end subroutine type352
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     Fan or pump (implicit flow) - temperature rise
! *
! * PURPOSE:        Calculates pressure, rotation speed, fluid heating
! *                 rate and outlet temperature from mass flow rate and
! *                 using a quadratic load line. Uses external iteration
! *                 on rotation speed.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. ifstatus: fan status (0 = off, 1 = on)                       (-)
! *  2. w       : mass flow rate                                  (kg/s)
! *  3. n       : rotation speed                                 (rev/s)
! *  4. tin     : inlet temperature                                  (C)
! *
! * OUTPUTS
! * =======
! *  1. dp      : pressure rise                                    (kPa)
! *  2. qa      : fluid stream heat addition rate                   (kW)
! *  3. power   : power consumption                                 (kW)
! *  4. n       : rotation speed                                 (rev/s)
! *  5. tout    : outlet temperature                                 (C)
! *
! * PARAMETERS
! * ==========
! *  1.         : coefficient of PN**0 in normalized head curve
! *  2.         : coefficient of PN**1 in normalized head curve
! *  ...
! *  5.         : coefficient of PN**4 in normalized head curve
! *  6.         : coefficient of PN**0 in normalized head curve
! *  7.         : coefficient of PN**1 in normalized head curve
! * ...
! * 10.         : coefficient of PN**4 in normalized head curve
! * 11. d       : diameter                                           (-)
! * 12. mode    : fluid: 1 = air, any other value + water            (-)
! * 13. wncritlo: normalised flow at lower bound of valid region     (-)
! * 14. wncritup: normalised flow at upper bound of valid region     (-)
! * 15. pstatsp : setpoint for the fan static pressure control loop(kPa)
! * 16. dsarea  : cross sect area of duct at static pressure sensor (m2)
! * 17. rduct   : resistance of duct sys. upstream of sensor (0.001/k.m)
! *
! * SAVED
! * =====
! *  1.         : flag used to indicate initialization has been performed
! *               (used in HVACSIM+ only)
! *  2. alo/n**2: A coefficient in lower region, divided by N squared
! *  3. rlo     : equivalent resistance in lower region
! *  4. aup/n**2: A coefficient in upper region, divided by N squared
! *  5. rup     : equivalent resistance in upper region
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  Fixed blade angle, no inlet guide vanes, all
!                        inefficiencies appear as heat in the fluid
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                September 30, 1992
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
!
!  Calculates actual pressure rise from simple model of external system:
!
!  DP = PTOTSP + RDUCT * W * ABS(W)
!
!  where PTOTSP is the total pressure corresponding to the set-point for
!  the fan static pressure control loop and RDUCT is the resistance of
!  the duct/pipe system upstream of the static pressure sensor.
!
!  The actual pressure rise is compared with the normalised pressure rise
!  calculated from the head curve in order to calculate the rotation speed
!  of the fan, hence its efficiency and hence the power consumption and
!  temperature rise. The rotation speed is both an input and an output.
!
! **********************************************************************
! * INTERNAL VARIABLES
! * ==================
! * funp()  : normalized head curve function
! * dfunp() : derivative of normalized head curve function wrt to flow
! * funn()  : normalized efficiency function
! * cp      : specific heat of working fluid
! * rho     : density of working fluid
! * wn      : normalized mass flow rate
! * dwn     : increment of valid range of normalized flow rate
! * anlo    : normalized A coefficient in lower region
! * rnlo    : normalized effective resistance in lower region
! * anup    : normalized A coefficient in upper region
! * rnup    : normalized effective resistance in upper region
! * roff    : effective resistance when fan/pump is off
! * pdrop   : pressure drop when fan/pump switched off
! * alo     : (un-normalized) A coefficient in lower region
! * alodnn  : A coefficient in lower region normalized by fan speed only
! * rlo     : (un-normalized) effective resistance in lower region
! * aup     : (un-normalized) A coefficient in upper region
! * aupdnn  : A coefficient in upper region normalized by fan speed only
! * rup     : (un-normalized) effective resistance in upper region
! * prise   : pressure rise when fan/pump switched on
! * pn      : normalized pressure rise
! * eff     : efficiency
! *
! **********************************************************************
!
!   Updated to Fortran 90 March 2, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type353(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=4,no=5,np=17,ns=5
        real,dimension(ni),intent(in)     :: xin
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real         :: n
        real         :: rhoa=1.2, rhow=1000., cpa=1.0, cpw=4.18
        real         :: z,funp,dfunp,funn,w,ti,d,wncritlo,wncritup,&
                        pstatsp,dsarea,rduct,cp,rho,dwn,wn,anlo,rnlo,&
                        anup,rnup,roff,pdrop,dpqudlin,prise,power,ptotsp,&
                        alodnn,rlo,alo,eff,aupdnn,rup,aup,pn,qa,to
        integer      :: i,mode,ifstatus

! **** Internal functions
! **** "Head curve" - normalized pressure rise as a function of normalized flow
        funp(z) = par_v(1)+z*(par_v(2)+z*(par_v(3)+z*(par_v(4)+z*par_v(5))))
! **** Derivative of head curve
        dfunp(z) = par_v(2)+z*(2.0*par_v(3)+z*(3.0*par_v(4)+z*4.0*par_v(5)))
! **** Efficiency curve
        funn(z) = par_v(6)+z*(par_v(7)+z*(par_v(8)+z*(par_v(9)+z*par_v(10))))
! **** Read in inputs
        ifstatus = nint(xin(1))
        w        = xin(2)
        n        = xin(3)
        ti       = xin(4)
! **** Read in parameters not used in head curve and efficiency functions
        d        = par_v(11)
        mode     = nint(par_v(12))
        wncritlo = par_v(13)
        wncritup = par_v(14)
        pstatsp  = par_v(15)
        dsarea   = par_v(16)
        rduct    = par_v(17)
! **** Outside the range of validity of the head curve polynomial, an ideal
! **** fan/pump model is used (pn=an-rn*wn*|wn|, where pn is the normalised
! **** pressure rise, wn is the normalised flow rate, an is +ve constant and
! **** rn is a +ve constant). an and rn are calculated so that the head curve
! **** and its gradient are continuous.

! **** Set specific heat and density for appropriate fluid
        if (mode==1) then
            cp=cpa
            rho=rhoa
        else
            cp=cpw
            rho=rhow
        endif
! *** Limit rotation speed to positive values
        n=max(0.0,n)
! **** "One-off" calculations
! **** Test for first call of simulation, assuming saved array is
! **** initialized to zero
        if (init==0 .and. saved_v(1)==0.0) then
            saved_v(1) = 1.0
! **** Check head curve is monotonically decreasing
! **** Test gradient at 50 points across user-defined range
            dwn = (wncritup-wncritlo)/50.
            wn  = wncritlo
            do i=0,50
                if (dfunp(wn)>0.0) then
                    write(*,*)&
                    'type 353: head curve slope positive at wn =',wn
                endif
                wn = wn+dwn
            enddo
! **** Determine coefficents of quadratic fits outside valid region
            anlo = funp(wncritlo)-wncritlo*dfunp(wncritlo)/2.0
            if (anlo<=0.0) stop 'type353: anlo<=0'
            rnlo = -dfunp(wncritlo)/(2.0*abs(wncritlo))
            if (rnlo<=0.0) stop 'type353: rnlo<=0'
            anup = funp(wncritup)-wncritup*dfunp(wncritup)/2.0
            if (anup<=0.0) stop 'type353: anup<=0'
            rnup = -dfunp(wncritup)/(2.0*abs(wncritup))
            if (rnup<=0.0) stop 'type353: rnup<=0'
! **** Store partly or completely un-normalized coefficients
            saved_v(2) = anlo*0.001*rho*d*d
            saved_v(3) = rnlo*0.001/(rho*d*d*d*d)
            saved_v(4) = anup*0.001*rho*d*d
            saved_v(5) = rnup*0.001/(rho*d*d*d*d)
        endif
! **** Calculate pressure rise and temperature rise
        if (ifstatus==0) then
! **** Fan/pump off
! **** Calculate the pressure drop assuming the resistance is equal to
! **** the equivalent resistance in the upper region
            roff  = saved_v(5)
            pdrop = dpqudlin(roff,w)
            prise = -pdrop
            power = 0.0
            n     = 0.0
        else
! **** Fan/pump on
! **** Normalized flow rate
            if (w<=0.0) then
! **** Flow is zero, so normalized flow is zero
                wn = 0.0
            elseif (n==0.0) then
! **** Rotation speed guess is zero - normalized flow rate is infinite
! **** or indeterminate, use 1.0 as a better starting guess
                wn = 1.0
            else
                wn = w/(rho*d*d*d*n)
           endif
! **** Total pressure rise from load line is sum of total pressure at
! **** static pressure sensor and pressure drop across fixed resistance
! **** upstream of sensor
            ptotsp = pstatsp+0.001*0.5*w*w/(rho*dsarea*dsarea)
            prise  = ptotsp+dpqudlin(rduct,w)
            if (prise<0.0) then
                stop 'type 353: negative pressure rise across fan'
            endif
            if (wn<wncritlo) then
! **** Below valid region
! **** Calculate rotation speed by comparing un-normalized and normalized
! **** estimates of the constant term in the quadratic extrapolation of
! **** the head curve
                alodnn = saved_v(2)
                rlo    = saved_v(3)
                alo    = prise+dpqudlin(rlo,w)
                n      = sqrt(alo/alodnn)
                eff    = funn(wncritlo)
            elseif (wn>wncritup) then
! **** Above valid region
! **** Calculate rotation speed by comparing un-normalized and normalized
! **** estimates of the constant term in the quadratic extrapolation of
! **** the head curve
                aupdnn = saved_v(4)
                rup    = saved_v(5)
                aup    = prise+dpqudlin(rup,w)
                n      = sqrt(aup/aupdnn)
                eff    = funn(wncritup)
            else
! **** In valid range
! **** Calculate rotation speed by comparing un-normalized and normalized
! **** estimates of the pressure rise
                pn     = funp(wn)
                n      = sqrt(prise/(0.001*pn*rho*d*d))
                eff    = funn(wn)
            endif
! **** Power consumption
            if (w>0.0.and.(prise)>0.0) then
                power  = (prise)*w/(rho*eff)
            else
                power  = 0.0
            endif
        endif
! **** Rate at which heat is added to fluid
        if (mode==1) then
! **** Air - include effect of fluid work
            qa = power
        else
! **** Water - exclude effect of fluid work
            qa = max(0.0,(prise)*w*(1./eff-1.)/rho)
        endif
! **** Temperature rise (no dynamics)
        if (w>0.0) then
            to = ti + qa/(w*cp)
        else
            to = ti
        endif
! **** Output
        yout(1) = prise
        yout(2) = qa
        yout(3) = power
        yout(4) = n
        yout(5) = to
! **** Allow freezing
        do i=1,no
            iostat(i) = 1
        enddo
! **** Return
        return
        end subroutine type353

