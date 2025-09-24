! **********************************************************************
     subroutine type107(xin,yout,par_v,saved_v,iostat)
!-----------------------------------------------------------------------
!
!     TYPE 107: Holder for data storage allocation

!               This holder keeps the numbers of state variables
!               same for all three different cases of the 9-zone
!               VCBT model.
!
!     Created : June 12, 2000     Cheol Park
!     Modified: July 17, 2001
!
! INPUTS:
!    xin    1  - 34  inputs
!    xin(1)  sa_tmp
!    xin(2)  ma_tmp
!    xin(3)  oa_tmp
!    xin(4)  ra_tmp
!    xin(5)  za_tmp_1
!    xin(6)  za_tmp_2
!    xin(7)  za_tmp_3
!    xin(8)  oa_hum
!    xin(9)  ra_hum
!    xin(10) sa_press
!    xin(11) zsa_tmp_1
!    xin(12) zsa_tmp_2
!    xin(13) zsa_tmp_3
!    xin(14) zflow_1
!    xin(15) zflow_2
!    xin(16) zflow_3
!    xin(17) sa_flow
!    xin(18) ra_flow
!    xin(19) sa_fan_stat
!    xin(20) ra_fan_stat
!    xin(21)  dmp_dem
!    xin(22)  cc_dem
!    xin(23)  zdmp_ctl_1
!    xin(24)  zdmp_ctl_2
!    xin(25)  zdmp_ctl_3
!    xin(26)  rhc_ctl_1
!    xin(27)  rhc_ctl_2
!    xin(28)  rhc_ctl_3
!    xin(29)  oa_dmp_ctl
!    xin(30)  ra_dmp_ctl
!    xin(31)  ea_dmp_ctl
!    xin(32)  sa_fan_speed
!    xin(33)  ra_fan_speed
!    xin(34)  hc_dem
!
! ------------------------------------------------------------------------

      use modsim_head
      implicit none
      integer,parameter                :: ndiffeq=0,ni=34,no=0,np=1,ns=0
      real,dimension(ni),intent(in)    :: xin
      real,dimension(no)               :: yout
      real,dimension(np),intent(in)    :: par_v
      real,dimension(ns)               :: saved_v
      integer,dimension(no)            :: iostat
      
      integer                          :: j
      real                             :: dummy

      do j= 1, ni
         dummy = xin(j)
      end do

      return
      end subroutine type107

! **********************************************************************
     subroutine type108(xin,yout,par_v,saved_v,iostat)
!-----------------------------------------------------------------------
!
!     TYPE 108: Holder for data storage allocation for the 6-zone model of
!               VCBT controllers
!
!     Created : June 6, 2007     Cheol Park
!
! INPUTS:
!    xin    1  - 46  inputs
!    xin(1)   sa_tmp
!    xin(2)   ma_tmp              
!    xin(3)   oa_tmp              
!    xin(4)   ra_tmp              
!    xin(5)   za_tmp_1            
!    xin(6)   za_tmp_2            
!    xin(7)   za_tmp_3            
!    xin(8)   za_tmp_4            
!    xin(9)   za_tmp_5            
!    xin(10)  za_tmp_6            
!    xin(11)  oa_hum              
!    xin(12)  ra_hum              
!    xin(13)  sa_press            
!    xin(14)  zsa_tmp_1           
!    xin(15)  zsa_tmp_2           
!    xin(16)  zsa_tmp_3           
!    xin(17)  zsa_tmp_4           
!    xin(18)  zsa_tmp_5           
!    xin(19)  zsa_tmp_6           
!    xin(20)  zflow_1             
!    xin(21)  zflow_2             
!    xin(22)  zflow_3             
!    xin(23)  zflow_4             
!    xin(24)  zflow_5             
!    xin(25)  zflow_6             
!    xin(26)  sa_flow             
!    xin(27)  ra_flow             
!    xin(28)  sa_fan_stat         
!    xin(29)  ra_fan_stat         
!    xin(30)  dmp_dem
!    xin(31)  cc_dem
!    xin(32)  zdmp_ctl_1
!    xin(33)  zdmp_ctl_2
!    xin(34)  zdmp_ctl_3
!    xin(35)  zdmp_ctl_4
!    xin(36)  zdmp_ctl_5
!    xin(37)  zdmp_ctl_6
!    xin(38)  rhc_ctl_1
!    xin(39)  rhc_ctl_2
!    xin(40)  rhc_ctl_3
!    xin(41)  rhc_ctl_4
!    xin(42)  rhc_ctl_5
!    xin(43)  rhc_ctl_6
!    xin(44)  sa_fan_speed
!    xin(45)  ra_fan_speed
!    xin(46)  hc_dem
!
! ------------------------------------------------------------------------

      use modsim_head
      implicit none
      integer,parameter                :: ndiffeq=0,ni=46,no=0,np=1,ns=0
      real,dimension(ni),intent(in)    :: xin
      real,dimension(no)               :: yout
      real,dimension(np),intent(in)    :: par_v
      real,dimension(ns)               :: saved_v
      integer,dimension(no)            :: iostat
      
      integer                          :: j
      real                             :: dummy

      do j= 1, ni
         dummy = xin(j)
      end do

      return
    end subroutine type108

    ! **********************************************************************
     subroutine type109(xin,yout,par_v,saved_v,iostat)
!-----------------------------------------------------------------------
!
!     TYPE 109: Holder for data storage allocation for input for HVAC-Cx AHU model (TODO other equipment models)
!     based on TYPE 108
!     Created : July 21, 2023     Mike Galler
!
! INPUTS:
!   1	Occ
!   2	SA sp
!   3	SA Temp
!   4	RA Temp
!   5	MA Temp
!   6	OA Temp
!   7	CW Valve
!   8	HW Valve
!   9	MA/RA damper
!   10	Recirc damper/return recovery
!   11	Humid valve
!   12	OA RH
!   13	Mixing Air Temp- recovery
!   14	RA sp
!   15	OA damper
!   16	RA humidity
     
! TYPAR entry
!     ******************************************************************************
!109  'Patch for data storage allocation                                       '
!  0  0 16  0 1              ! numbers of saved, diff. eq., xin, out, par
!  4  ' 1 Occupancy                                                            '
!  3  ' 2 sasp_tmp                                                             '
!  3  ' 3 sa_tmp                                                               '
!  3  ' 4 ra_tmp                                                               '
!  3  ' 5 ma_tmp                                                               '
!  3  ' 6 oa_tmp                                                               '
!  4  ' 7 chw_vlv                                                              '
!  4  ' 8 hw_vlv                                                               '
!  4  ' 9 ma/ra_damper                                                         '
!  4  '10 recirc damper                                                        '
!  4  '11 humid_vlv                                                            '
!  4  '12 oa_rh                                                                '
!  3  '13 ma_tmp_recovery                                                      '
!  3  '14 rasp_tmp                                                             '
!  4  '15 oa_damper                                                         '
!  4  '16 ra_rh                                                                '
!#
!#
!  1  'dummy '

 
     !MAG TODO: determine value for ni, fill in INPUTS: section above, edit TYPAR.dat file
      use modsim_head
      implicit none
      integer,parameter                :: ndiffeq=0,ni=16,no=0,np=1,ns=0
      real,dimension(ni),intent(in)    :: xin
      real,dimension(no)               :: yout
      real,dimension(np),intent(in)    :: par_v
      real,dimension(ns)               :: saved_v
      integer,dimension(no)            :: iostat
      
      integer                          :: j
      real                             :: dummy

      do j= 1, ni
         dummy = xin(j)
      end do

      return
    end subroutine type109

     
    
! **********************************************************************
!  
!   SUBROUTINE: Type365: Mixing of three moist air streams
!  
!   PURPOSE:    Calculates the mixed temperature, humidity and mass
!               flow rate resulting from mixing three air streams
!  
! **********************************************************************
!   INPUTS
!   ======
!    1. t(1)    : Stream 1 dry bulb temperature                      (C)
!    2. g(1)    : Stream 1 humidity ratio                        (kg/kg)
!    3. w(1)    : Stream 1 mass flow rate                         (kg/s)
!    4. t(2)    : Stream 2 dry bulb temperature                      (C)
!    5. g(2)    : Stream 2 humidity ratio                        (kg/kg)
!    6. w(2)    : Stream 2 mass flow rate                         (kg/s)
!    7. t(3)    : Stream 3 dry bulb temperature                      (C)
!    8. g(3)    : Stream 3 humidity ratio                        (kg/kg)
!    9. w(3)    : Stream 3 mass flow rate                         (kg/s)
!  
!   OUTPUTS
!   =======
!    1. tmix    : Mixed air dry bulb temperature                     (C)
!    2. gmix    : Mixed air humidity ratio                       (kg/kg)
!    3. wmix    : Mixed air mass flow rate                        (kg/s)
!  
!   PARAMETERS
!   ==========
!   (none)
!  
! **********************************************************************
!
!   Major restrictions:
!
!   Developer:           Philip Haves
!                        Loughborough University of Technology
!
!   Date:                November 8, 1995
!
!   Include files:       None
!   Subroutines called:  MOISTMIX
!   Functions  called:   None
!
!   Revision history:    None
!
!   Reference:
!
! **********************************************************************
!
!  Updated : June 12, 2007  Cheol Park, NIST
!            Converted Fortran 77 code into Fortran 90
!            with some changes.
!
! **********************************************************************
!
        subroutine type365(xin,yout,par_v,saved_v,iostat)   
!
        use modsim_head
        implicit none
        integer,parameter                 :: ni=9, no=3, np=1, ns=0
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        integer, parameter                :: nstream = 3
        real, dimension(nstream)          :: tx, gx, wx
        real                              :: tmix, gmix, wmix
        integer                           :: i

!                 Read in inputs
        tx(1) = xin(1)
        gx(1) = xin(2)
        wx(1) = xin(3)
        tx(2) = xin(4)
        gx(2) = xin(5)
        wx(2) = xin(6)
        tx(3) = xin(7)
        gx(3) = xin(8)
        wx(3) = xin(9)
!                 Calculate mixed air conditions
!                 moistmix also tests for reverse flows and sets the mixed
!                 condition to the forward flow condition, or to an average of
!                 the flow conditions if both are reverse

        call moistmix(tx,gx,wx,nstream,tmix,gmix,wmix)

!                 Assign outputs
        yout(1) = tmix
        yout(2) = gmix
        yout(3) = wmix

        do i=1,no
            iostat(i)=1                       ! 4/30/07
        enddo

        return
        end subroutine type365

! ***********************************************************************
! *
! * SUBROUTINE: Mixing box
! *
! * PURPOSE:    Calculates the mixed air flow rate, temperature and
! *             humidity and the extract air pressure. Uses Legg's
! *             correlations for the resistance of parallel and opposed
! *             blade dampers.
! *
! ***********************************************************************
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
! * 13. rfixout :  Fixed resistance in outside air branch    (0.001/k.m)
! * 14. rfixrec :  Fixed resistance in recirc air branch     (0.001/k.m)
! * 15. rfixexh :  Fixed resistance in exhaust air branch    (0.001/k.m)
! * 16. noninver: 0=invert recirc air damper, 1=not inverted         (-)
! *
! ***********************************************************************
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
!   REVISION HISTORY:    Minimum air damper is removed.  Cheol Park, NIST  8/31/01
!   Updated : June 12, 2007  Cheol Park, NIST
!             Converted Fortran 77 code into Fortran 90
!
!   REFERENCE:
!
! ***********************************************************************
! *
! * INTERNAL VARIABLES
! * ==================
! * routd   : resistance of outside air damper
! * rout    : total resistance of outside air branch
! * rrec    : total resistance of recirculation branch
! * rexh    : total resistance of exhaust branch
! * dp      : pressure drop across outside air branch
! * dpd     : pressure drop across outside air dampers
! * small   : threshold for significant flow rate
! *
! ***********************************************************************
!
        subroutine type371(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=11, no=7, np=16, ns=0
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real,dimension(0:1)               :: alegg=(/-1.51,-1.51/),&
                                             blegg=(/0.105,0.0842/)
        real,dimension(2)                 :: t,gx,w
        real                              :: small=1.0e-10,wcritl=3.e-2,&
                                             wcritu=6.e-2
        real          :: ropenout,ropenrec,ropenexh,fareaout,farearec,&
                         fareaexh,fleakout,fleakrec,fleakexh,rfixout,&
                         rfixrec,rfixexh,tout,gout,text,&
                         gext,pout,pexh,pmix,wext,cout,crec,cexh,&
                         routd,rdamper,rout,rrec,rexh,wexh,wrec,pext,&
                         dp,wout,wqudlin,wmix,tmix,gmix,wexhrev,&
                         wdummy,grec,trec
        integer       :: iparout,iparrec,iparexh,ifail,i,noninver

! *** Read and check parameters
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
        rfixout = par_v(13)
        rfixrec = par_v(14)
        rfixexh = par_v(15)
        noninver= nint(par_v(16))
! *** Read in inputs
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
! *** Calculate resistance of each branch as sum of damper resistance and
! *** fixed resistance
        routd = rdamper(cout,ropenout,fleakout,fareaout,alegg(iparout),&
                       blegg(iparout),iparout)
        rout  = routd + rfixout
        if (rout<=0.0) stop&
            'type371: resistance in outside branch not greater than 0'
! *** Recirc and exhaust branches - one damper in series with fixed
! *** resistance
        rrec  = rdamper(crec,ropenrec,fleakrec,farearec,alegg(iparrec),&
                        blegg(iparrec),iparrec) + rfixrec
        if (rrec<=0.0) stop&
            'type371: resistance in recirc branch not greater than 0'
        rexh  = rdamper(cexh,ropenexh,fleakexh,fareaexh,alegg(iparexh),&
                        blegg(iparexh),iparexh) + rfixexh
        if (rexh<=0.0) stop&
            'type371: resistance in exhaust branch not greater than 0'
! *** Calculate recirc and exhaust flow rates
! *** Call flow split routine
        call flowsplt(wext,pmix,pexh,0.0,rrec,rexh,wcritl,wcritu,rtolx,&
                      pext,wrec,wexh,ifail)
! *** Check for unsuccesful completion
        if (ifail==1) then
! *** One or more resistances negative
            stop 'type 327: resistances must not be negative'
        elseif (ifail==2) then
! *** Zero resistance for both outlet branches
            stop 'type 327: either rret or rexh must be non-zero'
        endif
! *** Outside air mass flow rate
        dp   = pout-pmix
        wout = wqudlin(rout,dp)
! *** Calculate mixed air flow rate as sum of outside and recirc flow rates
        wmix = wout+wrec
! *** Having calculated flows, calculate mixed temperature and humidity
! *** model assumes wmix and wext both non-negative. warn if not true.
        if (wmix<-1.0e-4 .and. lastit==1)&
            write(*,*) 'type 327: wmix < -1.0e-4'
        if (wext<-1.0e-4 .and. lastit==1)&
            write(*,*) 'type 327: wext < -1.0e-4'
! *** Check for reverse mixed flow
        if (wmix<=small) then
! *** Reverse mixed flow, take mixed temperature and humidity ratio to be
! *** an average of the outside and extract conditions so as to provide a
! *** reasonable guess for subsequent iterations
            tmix = (tout+text)/2.0
            gmix = (gout+gext)/2.0
        elseif (wmix>small .and. wext<=small) then
! *** Reverse extract flow - mixed condition is the same as outside condition
            tmix = tout
            gmix = gout
        elseif (wmix>small .and. wext>small ) then
! *** Normal, forward mixed and extract flows - calculate mixed conditions
! *** first calculate recirc air conditions if recirc flow positive
            if (wrec>small) then
! *** Test for reverse exhaust flow
                if (wexh<small) then
! *** Reverse exhaust flow - recirc flow is a mixture of extract and outside
! *** conditions
                    wexhrev = -wexh
                    wexhrev = -wexh
                    t(1) = text
                    gx(1) = gext
                    w(1) = wext
                    t(2) = tout
                    gx(2) = gout
                    w(2) = wexhrev
                    call moistmix(t,g,w,2,trec,grec,wdummy)
                else
! *** Normal, forward exhaust flow - recirc condition is extract condition
                    trec = text
                    grec = gext
                endif
            endif
! *** Calculate mixed air conditions
! *** Mixed condition is a combination of the recirc and outside conditions
! *** moistmix also tests for reverse flows and sets the mixed condition
! *** to the forward flow condition, or to an average of the flow conditions
! *** if both are reverse
            t(1) = trec
            gx(1) = grec
            w(1) = wrec
            t(2) = tout
            gx(2) = gout
            w(2) = wout
            call moistmix(t,g,w,2,tmix,gmix,wdummy)
        else
            stop 'type371: logic fault'
        endif
! *** Assign outputs
         
        yout(1) = tmix
        yout(2) = gmix
        yout(3) = wmix
        yout(4) = pext
        yout(5) = wout
        yout(6) = wrec
        yout(7) = wexh
! *** Allow freezing of algebraic variables
        do i=1,no
            iostat(i)=1
        enddo

        return
        end subroutine type371

! ***********************************************************************
! *
! * SUBROUTINE:     Room  (analytical integration)
! *
! * PURPOSE:        Calculate the temperature and humidity ratio in the
! *                 room by performing heat balances and moisture balances.
! *                 Include the effects of air flows from adjacent zones
! *                 and any (inward) air leakage and reverse return flow in addition to the HVAC supply
! *                 air flow. Use DIFFEQ to peform analytical integration.
! *                 Heat gain through the wall surfaces, whcih is computed
! *                 by the shell model, is included.
! *
! **********************************************************************
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
!
!  Revised for the shell model by Cheol Park, NIST  October 27, 1999
!
!  Plenum is not considered.  Nodal point model is changed.
!  The heat gain through the walls is entered as an input.
!  The heat transfer rates are obtained from the shell model, which implements
!  the conduction transfer functions for multilayer structs.
!
!  Updated:  November 12, 2002  Cheol Park
!            Electric use by lighting & equipment is reported.
!  Updated : June 12, 2007  Cheol Park, NIST
!            Converted Fortran 77 code into Fortran 90
!
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
! * 14. fracocc : Fractional occupancy                               (-)
! * 15. fraclig : Fractional lighting heat gain                      (-)
! * 16. fracpow : Fractional equipment heat gain                     (-)
! * 17. troom   : Room temperature                                   (C)
! * 18. groom   : Room humidity ratio                            (kg/kg)
! * 19. qwall   : Heat gain from all surfaces (walls, doors, windows,
! *               floor, and ceiling)                               (kW)
! *
! * OUTPUTS
! * =======
! *  1. troomn  : Room temperature                                   (C)
! *  2. groomn  : Room humidity ratio                            (kg/kg)
! *  3. tret    : Return temperature                                 (C)
! *  4. qsensr  : Sensible heat gains of room                       (kW)
! *  5. wvapr   : Water vapour gains of room                      (kg/s)
! *  6. qrlig + qpow: Electric power by lighting & equipment        (kW)
! *
! * PARAMETERS
! * ==========
! *  1. xcap    : Room air capacity multiplier                       (-)
! *  2. crair   : Capacitance of room air node (unmodified)       (kJ/K)
! *  3. vroom   : Volume of room                                    (m3)
! *  4. noccup  : Number of occupants                                (-)
! *  5. qlite   : Lighting heat gain                                (kW)
! *  6. qequp   : Equipment heat gain                               (kW)
! *  7. nfile   : Zone number (parameter file='zoneN.par', N > 0)    (-)
! *
! * SAVEd
! * =====
! *  1. time    : Time of previous call of TYPE
! *  2. troomn  : Room temperature from previous call
! *  3. troomp  : Room temperature from previous step time
! * 10. groom   : Room humidity ratio from previous call
! * 11. groomp  : Room humidity ratio from previous step time
! *
! **********************************************************************
!
! * INTERNAL VARIABLES
! * ==================
! * mroom   : mass of air in room
! * cpa     : specific heat of dry air
! * cpg     : specific heat of water vapor
! * hfg     : latent heat of water vapor
! * rhoair  : density of air
! * qsperson: sensible heat gain from a person
! * qlperson: latent heat gain from a person
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
! * qpow    : heat gain to room from equipment
! * aa      : A coefficent in dT/dt = A*T + B
! * bb      : B coefficent in dT/dt = A*T + B
! *
! **********************************************************************
!
        subroutine type428(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndiffeq=2,nsv=1+ndiffeq*2
        integer,parameter                 :: ni=19, no=6, np=7, ns=nsv
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

        real        :: mroom,noccup
        real        :: cpa=1.0,cpg=1.805,hfg=2501.,rhoair=1.2,&
                       qsperson=0.095,qlperson=0.045
        real        :: qpow,qrlig,qsensr,tsup,gsup,wsup,tinz1,ginz1,&
                       winz1,tinz2,ginz2,winz2,wret,tamb,gamb,wleak,&
                       fracoc,fracli,fracpo,troom,groom,xcap,crair,&
                       vroom,qlite,qequp,troomp,groomp,wfan,diwsup,&
                       diwinz,diwlea,diwret,diwfan,csup,cinz1,&
                       cinz2,cret,cleak,cfan,cr,qoccs,wvapr,aa,bb,&
                       troomb,groomb,qwall,tret,groomn,troomn,&
                       fracocc,fraclig,fracpow,diwinz1,diwinz2,diwleak
        integer     :: is,i,nfile

        namelist  /inputs/nfile,time,qwall
        namelist  /outputs/time,troomn,groomn,tret,qsensr,qrlig,qpow

! Read in inputs

         tsup    = xin(1)
         gsup    = xin(2)
         wsup    = xin(3)
         tinz1   = xin(4)
         ginz1   = xin(5)
         winz1   = xin(6)
         tinz2   = xin(7)
         ginz2   = xin(8)
         winz2   = xin(9)
         wret    = xin(10)
         tamb    = xin(11)
         gamb    = xin(12)
         wleak   = xin(13)
         fracocc = xin(14)
         fraclig = xin(15)
         fracpow = xin(16)
         troom   = xin(17)    
         groom   = xin(18)
         qwall   = xin(19)

!     Read in parameters

        xcap     = par_v(1)
        crair    = par_v(2)
        vroom    = par_v(3)
        noccup   = par_v(4)
        qlite    = par_v(5)
        qequp    = par_v(6)
        nfile    = nint(par_v(7))

!      if (mod(time, 3600.0)==0) then
!         print inputs
!       end if

!     Initialize at beginning of simulation

        if (itime <= 1) then
            if (init==0 .or. saved_v(1) > time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then

!     Initialize room conditions to 20 deg C, 50% rh

                do is = 2,3
                    saved_v(is) = 20.0
                end do
                do is = 4,5
                    saved_v(is) = 0.0074
                end do
            endif
        endif
        if (time > saved_v(1)) then

!     First call of timestep - update previous sample instant values

            do is = 2,nsv-1,2
                saved_v(is+1) = saved_v(is)
            end do
        endif

!     Update previous values

        troomp  = saved_v(3)
        groomp  = saved_v(5)

!     Set up and solve heat and moisture balances
!     Calculate local extract from mass balance

        wfan = wsup + winz1 - winz2 - wleak - wret

!     Calculate capacity rates entering room.
!     nb wsup and winz1 are positive entering the room, winz2, wleak, wret and
!     wfan are positive leaving the room. Capacity rates are non-zero if the
!     flow is into the zone.

        diwsup  = max(0.0,wsup)
        diwinz1 = max(0.0,winz1)
        diwinz2 = max(0.0,-winz2)
        diwleak = max(0.0,-wleak)
        diwret  = max(0.0,-wret)
        diwfan  = max(0.0,-wfan)
        csup    = diwsup*(cpa+gsup*cpg)
        cinz1   = diwinz1*(cpa+ginz1*cpg)
        cinz2   = diwinz2*(cpa+ginz2*cpg)
        cret    = diwret*(cpa+groom*cpg)    ! 11/24/99
        cleak   = diwleak*(cpa+groom*cpg)   !
        cfan    = diwfan*(cpa+groom*cpg)    !
!        cret    = diwret*(cpa+gamb*cpg)
!        cleak   = diwleak*(cpa+gamb*cpg)
!        cfan    = diwfan*(cpa+gamb*cpg)
        cr      = crair*xcap

!     Calculate internal gains
!             Sensible gains from occupants

        qoccs = fracocc*noccup*qsperson

!             Sensible gains from lighting

        qrlig = fraclig*qlite

!             Sensible gains from equipment

        qpow   = fracpow*qequp

!       Sensible gains to room

        qsensr = qoccs+qrlig+qpow+qwall      ! qwall added  10/27/99

!       Latent gains from occupants

        wvapr = (fracocc*noccup*qlperson)/hfg

!     Heat balance on room air

        aa = -(1./cr)*(csup+cinz1+cinz2+cleak+cret+cfan)
        bb =  (1./cr)*(tsup*csup+tinz1*cinz1+tinz2*cinz2+&
                       tamb*(cleak+cret+cfan)+qsensr)

        call diffeq(time,aa,bb,troomp,troomn,troomb)  ! analytical integration

!       if (mod(time, 1800.0)==0) then
!         print *," time, troomp, troomn :", time, troomp, troomn
!       end if

!     Moisture balance on room air

        mroom = vroom*rhoair
        aa = -(1./mroom)*(diwsup+diwinz1+diwinz2+diwleak+diwret+diwfan)
        bb =  (1./mroom)*(gsup*diwsup+ginz1*diwinz1+ginz2*diwinz2+&
                          gamb*(diwleak+diwret+diwfan)+wvapr)

        call diffeq(time,aa,bb,groomp,groomn,groomb)  ! analytical integration

!     Extracted air temperature

        tret = troomn

!     Save provisional values to be used in updating at next step time

        saved_v(2) = troomn
        saved_v(4) = groomn

!     Save time of current call

        saved_v(1) = time

!     Output

        yout(1) = troomn
        yout(2) = groomn
        yout(3) = tret
        yout(4) = qsensr
        yout(5) = wvapr
        yout(6) = qrlig + qpow                  ! 11/12/02

!       if (mod(time, 3600.0)==0) then
!         print outputs
!       end if

        do i=1,no
            iostat(i) = 0
        end do

        return
        end subroutine type428
! ***********************************************************************
! *
! * SUBROUTINE:     E51 supply fan control
! *
! * PURPOSE:        Control supply static pressure by varying supply fan
! *                 speed. Override fan speed if static pressure exceeds
! *                 high limit.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. ps      : static pressure sensor                           (kPa)
! *  2. psset   : static pressure setpoint                         (kPa)
! *  3. sfstatus: supply fan status (1 = on, 0 = off)                (-)
! *  4. rfstatus: return fan status (1 = on, 0 = off)                (-)
! *  5. occupancy: building occupancy (1 = on, 0 = off)              (-)      ! 11/05/02
! *
! * OUTPUTS
! * =======
! *  1. sfspd   : supply fan speed (0-1)                             (-)
! *
! * PARAMETERS
! * ==========
! *  1. propb   : proportional band                                (kPa)
! *  2. tint    : integral time                                      (s)
! *  3. tdt     : derivative time                                    (s)
! *  4. dband   : deadband                                         (kPa)
! *  5. propbhl : high limit override propotional band             (kPa)
! *  6. pshlimit: high limit override setpoint                     (kPa)
! *  7. closloop: control mode (0=open loop, 1=closed loop)          (-)
! *  8. sfspdman: open loop supply fan speed (0-1)                   (-)
! *  9. tsamp   : sampling interval                                  (s)
! * 10. nfile   : controller number (parameters in file contN.par)   (-)
! *
! * SAVED
! * =====
! *  1. time    : time of previous call
! *  2. nsample : sample number of previous controller execution
! *  3. nsample : sample number of previous controller sample
! *  4. intp    : integral term from previous call
! *  5. intp    : integral term from previous sample
! *  6. difp    : derivative term from previous call
! *  7. difp    : derivative term from previous sample
! *  8. errp    : error from previvous call
! *  9. errp    : error from previvous sample
! * 10. pidoutp : output from previous call
! * 11. pidoutp : output from previous sample
! * 12-19.      : (PAR(1)-PAR(8) read from file
! * 20. sfspd   : controller output from previous sample
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Li Mei and Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 22, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   RETOLOG,DEADBANDT,PIDCONT,SWITCH,SPAN,
!                        LOGICAND
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
!   Revised:  Nov. 5, 2002 C.P.
!   Updated : June 12, 2007  Cheol Park, NIST
!             Converted Fortran 77 code into Fortran 90
!
! **********************************************************************
!
        subroutine type471(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=11,nfp=8
        integer,parameter                 :: ni=5,no=1,np=11,ns=nsv+nfp+1
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

        real         :: ps,psset,propb,tint,tdt,dband,propbhl,pshlimit,&
                        sfspdman,tsamp,sfspdmin,difp,errp,pidoutp,&
                        pssetdb,pidout,pidcont,contout,switch,dummy4,&
                        dummy3,dummy2,dummy1,sfspdhl,sfspd
        integer      :: nfile,is,nsample,i

! *** Declaration of local variables
        real         :: intp
        logical      :: retolog,sfstatus,rfstatus,&
                        closloop,auxdis,logicand,occupancy

! *** Read in inputs
        ps        = xin(1)
        psset     = xin(2)
        sfstatus  = retolog(xin(3))
        rfstatus  = retolog(xin(4))
        occupancy = retolog(xin(5))                                ! 11/5/02

! *** Read in parameters
        propb     = par_v(1)
        tint      = par_v(2)
        tdt       = par_v(3)
        dband     = par_v(4)
        propbhl   = par_v(5)
        pshlimit  = par_v(6)
        closloop  = retolog(par_v(7))
        sfspdman  = par_v(8)
        tsamp     = par_v(9)
        sfspdmin  = par_v(10)     ! minimun fan speed during unoccupied period 11/5/02
        nfile     = par_v(11)
! *** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
                saved_v(2) = 0.
            endif
            if (init==0) then
                do 100 is = 4,(nsv-1),2
                    saved_v(is) = 0.0
100             continue
            endif
        endif
! *** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! *** Run controller if a sample instant
        nsample=nint(saved_v(3))+1
        if (time>=(nsample*tsamp) .or. time==0.0) then
! *** Use parameters from file if controller number not zero
!            if (nfile>0) then
!  *** First call of time-step - read parameters from file and store
!                if (time>saved_v(1)) then
!                    call rfile(nfile,'cont',nfp,481,fpar)
!                    do 150 i=1,nfp
!                        saved_v(nsv+i)=fpar(i)
!150                continue
!                endif
!  *** Get parameters that were read form file and then stored
!                propb    = saved_v(nsv+1)
!                tint     = saved_v(nsv+2)
!                tdt      = saved_v(nsv+3)
!                dband    = saved_v(nsv+4)
!                propbhl  = saved_v(nsv+5)
!                pshlimit = saved_v(nsv+6)
!                closloop = retolog(saved_v(nsv+7))
!                sfspdman = saved_v(nsv+8)
!            endif
! *** Run controller if mode is closed loop else use "manual" value
            if (closloop) then
! *** Closed loop mode
                if (time>saved_v(1)) then
! *** First call of timestep
! *** Update previous sample instant values
                    do is=4,(nsv-1),2
                        saved_v(is+1)=saved_v(is)
                    enddo
                endif
! *** Update previous values
                intp      = saved_v(5)
                difp      = saved_v(7)
                errp      = saved_v(9)
                pidoutp   = saved_v(11)
! *** Determine output of main control loop then apply over-ride if
! *** static pressure too high
! *** Apply deadband around set-point
                call deadband(ps,psset,dband,pssetdb)
! *** pid controller
                pidout    = pidcont(ps,pssetdb,propb,tint,tdt,intp,&
                                    difp,pidoutp,errp,tsamp,1.0,0.0)
! *** Determine whether auxiliary input switch is enabled or disabled
! *** Disable auxiliary input if supply fan and return fan ok
                auxdis    = logicand(sfstatus,rfstatus)
! *** Determine controller output signal - Select pidout if auxdis is
! *** true else select 0.0
                contout   = switch(auxdis,pidout,0.0)
! *** Determine high limit over-ride value - proportional control
                sfspdhl   = pidcont(ps,pshlimit,propbhl,0.0,0.0,dummy1,&
                                  dummy2,dummy3,dummy4,tsamp,1.0,0.0)
! *** Select minimun controller output signal
                if(occupancy) then
                   sfspd  = min(contout,sfspdhl)
                else
                   sfspd  = sfspdmin                                ! 11/5/02
                endif

! *** Save provisional values to be used in updating at next sample instant
                saved_v(4)  = intp
                saved_v(6)  = difp
                saved_v(8)  = errp
                saved_v(10) = pidout
! *** Open loop mode
            else
                sfspd     = sfspdman
            endif

! *** Save current sample number and output
            saved_v(2)  = float(nsample)
            saved_v(20) = sfspd
        else
! *** Not a sample instant, set output to value from previous sample
! *** instant
            sfspd = saved_v(20)
        endif
! *** Save time of current call
        saved_v(1) = time
! *** Output
        yout(1) = sfspd
! *** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo
! *** Return
        return
        end subroutine type471
! ***********************************************************************
! * SUBROUTINE:     E51 return fan volume matching control
! *
! * PURPOSE:        Calculate return fan control signal from difference
! *                 between supply and return air volumetric flow rate.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. supflow : supply air volume flow rate sensor              (m3/s)
! *  2. retflow : return air volume flow rate sensor              (m3/s)
! *  3. sfstatus: supply fan status (1 = on, 0 = off)                (-)
! *  4. rfstatus: return fan status (1 = on, 0 = off)                (-)
! *  5. dflowset: setpoint for flow rate difference               (m3/s)
! *  6. occupancy
! *
! * OUTPUTS
! * =======
! *  1. rfspd   : return fan speed (0-1)                             (-)
! *
! * PARAMETERS
! * ==========
! *  1. propb   : proportional band                               (m3/s)
! *  2. tint    : integral time                                      (s)
! *  3. tdt     : derivative time                                    (s)
! *  4. deadb   : deadband                                        (m3/s)
! *  5. closloop: control mode (0 = open loop, 1 = closed loop)      (-)
! *  6. rfspdman: open loop return fan speed (0-1)                   (-)
! *  7. tsamp   : sampling interval                                  (s)
! *  8. nfile   : controller number (parameters in file contN.par)   (-)
! *  9. rfspdmin : minimum return fan speed
! *
! * SAVED
! * ======
! *  1. time    : time of previous call
! *  2. nsample : sample number of previous controller execution
! *  3. nsample : sample number of previous controller sample
! *  4. intp    : integral term from previous call
! *  5. intp    : integral term from previous sample
! *  6. difp    : derivative term from previous call
! *  7. difp    : derivative term from previous sample
! *  8. errp    : error from previvous call
! *  9. errp    : error from previvous sample
! * 10. pidoutp : PID controller output from previous call
! * 11. pidoutp : PID controller output from previous sample
! * 12-17.      : (PAR(1)-PAR(6) read from file
! * 18. rfspd   : return fan speed from previous sample
! *********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Li Mei and Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 22, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  DEADBAND
!   FUNCTIONS  CALLED:   RETOLOG,PIDCONT,LOGICAND,SWITCH,SPAN,
!                        LOGICAND
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
!   Updated : June 12, 2007  Cheol Park, NIST
!             Converted Fortran 77 code into Fortran 90
!
! **********************************************************************
!
        subroutine type480(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=11,nfp=6
        integer,parameter                 :: ni=6,no=1,np=9,ns=nsv+nfp+1
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

        real      :: supflow,retflow,dflowset,propb,tint,tdt,deadb,&
                     rfspdman,tsamp,rfspdmin,difp,errp,pidoutp,dflow,&
                     setdb,pidout,pidcont,rfspd,switch
        integer   :: nfile,i,is,nsample

! *** Declaration of local variables
        real         :: intp
        logical      :: retolog,sfstatus,rfstatus,&
                        closloop,auxdis,logicand,occupancy

! *** Read in inputs
        supflow  = xin(1)
        retflow  = xin(2)
        sfstatus = retolog(xin(3))
        rfstatus = retolog(xin(4))
        dflowset = xin(5)
        occupancy = retolog(xin(6))                                ! 2/20/03

! *** Read in parameters
        propb    = par_v(1)
        tint     = par_v(2)
        tdt      = par_v(3)
        deadb    = par_v(4)
        closloop = retolog(par_v(5))
        rfspdman = par_v(6)
        tsamp    = par_v(7)
        nfile    = par_v(8)
        rfspdmin = par_v(9)                                           ! 2/20/03
! *** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
                saved_v(2) = 0.
            endif
            if (init==0) then
                do is = 4,(nsv-1),2
                    saved_v(is) = 0.0
                enddo
            endif
        endif
! *** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! *** Run controller if a sample instant
        nsample=nint(saved_v(3))+1
        if (time>=(nsample*tsamp) .or. time==0.0) then
! *** Use parameters from file if controller number not zero
!            if (nfile>0) then
!  *** First call of time-step - read parameters from file and store
!               if (time>saved_v(1)) then
!                    call rfile(nfile,'cont',nfp,482,fpar)
!                    do 150 i=1,nfp
!                        saved_v(nsv+i)=fpar(i)
! 150                continue
!               endif
!  *** Get parameters that were read form file and then stored
!               propb    = saved_v(nsv+1)
!               tint     = saved_v(nsv+2)
!               tdt      = saved_v(nsv+3)
!               deadb    = saved_v(nsv+4)
!               closloop = retolog(saved_v(nsv+5))
!               rfspdman = saved_v(nsv+6)
!            endif
! *** Run controller if mode is closed loop else use "manual" value
            if (closloop) then
! *** Closed loop mode
                if (time>saved_v(1)) then
! *** First call of timestep
! *** Update previous sample instant values
                    do is=4,(nsv-1),2
                        saved_v(is+1) = saved_v(is)
                    enddo
                endif
! *** Update previous values
                intp    = saved_v(5)
                difp    = saved_v(7)
                errp    = saved_v(9)
                pidoutp = saved_v(11)
! *** Difference of supply and return air volume
                dflow   = supflow - retflow
! *** Apply deadband around set-point
                call deadband(dflow,dflowset,deadb,setdb)
! *** Pid controller
                pidout  = pidcont(dflow,setdb,propb,tint,tdt,intp,&
                                 difp,pidoutp,errp,tsamp,1.0,0.0)
! *** Determine whether auxiliary input switch is enabled or disabled
! *** Disable auxiliary input if supply fan and return fan ok
                auxdis = logicand(sfstatus,rfstatus)
! *** Determine controller output signal - - select pidout if auxdis is
! *** true else select 0.0
! *** Select minimun controller output signal
                if(occupancy) then
                   rfspd = switch(auxdis,pidout,0.0)
                else
                   rfspd  = rfspdmin                                ! 2/20/03
                endif

! *** Save provisional values to be used in updating at next sample instant
                saved_v(4)  = intp
                saved_v(6)  = difp
                saved_v(8)  = errp
                saved_v(10) = pidout
! *** Open loop mode
            else
                rfspd     = rfspdman
            endif
! *** Save current sample number and output value
            saved_v(2)  = float(nsample)
            saved_v(18) = rfspd
        else
! *** Not a sample instant, set output(s) to value(s) from previous sample
! *** instant
            rfspd = saved_v(18)
        endif
! *** Save time of current call
        saved_v(1) = time
! *** Output
        yout(1) = rfspd
! *** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo
! *** Return
        return
        end subroutine type480
! ***********************************************************************
!
!     Type492 : Supply air temperature control with both heating and
!               cooling coil valves
!               Type486.for was modified to include the heating coil
!               valve control.
!
!     Created: June 23, 2000   Cheol Park, NIST
!     Modified: August 30, 2000
!     Revised:  January 12, 2005, June 11, 2007
!     Updated : June 12, 2007  Cheol Park, NIST
!               Converted Fortran 77 code into Fortran 90
!
! ***********************************************************************
! * SUBROUTINE:    E51 supply air temperature control
! *
! * PURPOSE:       Calculate demanded position of cooling coil valve and
! *                mixing box dampers. Override at low temperatures or
! *                if supply fan status is not OK. At low temperatures,
! *                set cooling coil valve fully open (to prevent freezing).
! *                (The mixed air damper control (TYPE485) will set the
! *                mixing box to full recirc and the minimum outside air
! *                damper controller (TYPE484) will close the minimum
! *                outside air damper at low temperatures.)
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. tsup    : supply air temperature sensor                      (C)
! *  2. sfstatus: supply fan status (1 = on, 0 = off)                (-)
! *  3. lowtovr : low temperature override signal (1=TRUE, 0=FALSE)  (-)
! *  4. tset    : supply air temperature setpoint                    (C)
! *
! * OUTPUTS
! * =======
! *  1. ddem    : damper cooling demand (0-1)                        (-)
! *  2. cdem    : cooling coil valve demand (0-1)                    (-)
! *  3. hdem    : heating coil valve demand (0-1)                    (-)
! *
! * PARAMETERS
! * ==========
! *  1. propb   : proportional band                                  (K)
! *  2. tint    : integral time                                      (s)
! *  3. tdt     : derivative time
! *  4. u_min   : minimum control value
! *  5. hbreak  : breakpoint between damper and heating coil demand(0-2)
! *  6. cbreak  : breakpoint between damper and cooling coil demand(0-2)
! *  7. u_max   : maximum control value
! *  8. deadb   : deadband                                           (K)
! *  9. closloop: control mode (0=open loop, 1=closed loop)          (-)
! * 10. cdemman : open loop cooling coil demand (0-1)                (-)
! * 11. hdemman : open loop heating coil demand (0-1)                (-)
! * 12. tsamp   : sampling interval                                  (s)
! * 13. nfile   : controller number (parameters in file contN.par)   (-)
! *
! * SAVED
! * =====
! *  1. time    : time of previous call
! *  2. nsample : sample number of previous controller execution
! *  3. nsample : sample number of previous controller sample
! *  4. intp    : integral term from previous call
! *  5. intp    : integral term from previous sample
! *  6. difp    : derivative term from previous call
! *  7. difp    : derivative term from previous sample
! *  8. errp    : error from previous call
! *  9. errp    : error from previous sample
! * 10. pidoutp : output from previous call
! * 11. pidoutp : output from previous sample
! * 12-18.      : (PAR(1)-PAR(7) read from file
! * 19. ddem    : damper demand from previous sample
! * 20. cdem    : cooling coil demand from previous sample
! * 21. hdem    : heating coil demand from previous sample
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  Based on the strategy used in Building E-51
!                        at MIT
!
!   DEVELOPER:           Li Mei and Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 22, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  DEADBAND
!   FUNCTIONS  CALLED:   RETOLOG,PIDCONT,LOGICNOT,LOGICAND,SWITCH,SPAN
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
!
        subroutine type492(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=11,nfp=11
        integer,parameter                 :: ni=4,no=3,np=13,ns=nsv+nfp+3
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

! *** Declaration of local variables
        real        :: tsup,tset,propb,tint,tdt,u_min,hbreak,cbreak,&
                       u_max,deadb,cdemman,hdemman,tsamp,difp,errp,&
                       pidoutp,tsetdb,pidout,pidcont,contout,hdem,&
                       span,cdem,ddem
        integer     :: nfile,is,nsample,i

        real        :: intp
        logical     :: retolog,logicnot,logicand,sfstatus,lowtovr,&
                       nlowtovr,auxdis,closloop
        integer     :: control_stat
        integer     :: itype=492

! *** Read in inputs
        tsup       = xin(1)
        sfstatus   = retolog(xin(2))
        lowtovr    = retolog(xin(3))
        tset       = xin(4)

! *** Read in parameters
        propb      = par_v(1)
        tint       = par_v(2)
        tdt        = par_v(3)
        u_min      = par_v(4)
        hbreak     = par_v(5)
        cbreak     = par_v(6)
        u_max      = par_v(7)
        deadb      = par_v(8)
        closloop   = retolog(par_v(9))
        cdemman    = par_v(10)
        hdemman    = par_v(11)
        tsamp      = par_v(12)
        nfile      = nint(par_v(13))

! *** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
                saved_v(2) = 0.
            endif
            if (init==0) then
                do is = 4,nsv-1,2
                    saved_v(is) = 0.0
                enddo
            endif
        endif

! *** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif

! *** Run controller if a sample instant

        nsample=nint(saved_v(3))+1

        if (time>=(nsample*tsamp) .or. time==0.0) then

! **** Use parameters from file if controller number not zero
           if (nfile>0) then
! **** First call of time-step - read parameters from file and store
                if (time>saved_v(1)) then
                    call rfile(nfile,'cont',nfp,itype,fpar)
                    do i=1,nfp
                        saved_v(nsv+i)=fpar(i)
                    enddo
                endif
                propb      = saved_v(nsv+1)
                tint       = saved_v(nsv+2)
                tdt        = saved_v(nsv+3)
                u_min      = saved_v(nsv+4)
                hbreak     = saved_v(nsv+5)
                cbreak     = saved_v(nsv+6)
                u_max      = saved_v(nsv+7)
                deadb      = saved_v(nsv+8)
                closloop   = retolog(saved_v(nsv+9))
                cdemman    = saved_v(nsv+10)
                hdemman    = saved_v(nsv+11)
            endif

! *** Run controller if mode is closed loop else use "manual" value

            if (closloop) then

! *** Closed loop mode
                if (time>saved_v(1)) then

! *** First call of timestep
! *** Update previous sample instant values

                    do is=4,nsv-1,2
                        saved_v(is+1)=saved_v(is)
                    enddo
                endif

! *** Update previous values
                intp = saved_v(5)
                difp = saved_v(7)
                errp = saved_v(9)
                pidoutp = saved_v(11)

! *** Apply deadband around set-point


                call deadband(tsup,tset,deadb,tsetdb)
! *** Pid controller

                pidout = pidcont(tsup,tsetdb,propb,tint,tdt,intp,difp,&
                                 pidoutp,errp,tsamp,u_max,u_min)

! *** Determine whether auxiliary input switch is enabled or disabled
! *** compliment of low temperature over-ride

!!                nlowtovr = logicnot(lowtovr)

! *** Disable auxiliary input if fan ok and low temperature over-ride
! *** not set

!!                auxdis = logicand(sfstatus,nlowtovr)

! *** Determine auxiliary input value - cooling coil valve fully open
! *** if low temperature over-ride set else fully closed if fan status
! *** not ok and low temperature over-ride not set

!!                auxinp = switch(lowtovr,2.0,0.0)

! *** Determine controller output signal - select pidout if auxdis is
! *** true else select auxinp

!!                contout = switch(auxdis,pidout,auxinp)

! *** Sequence demands

          contout = pidout             ! 8/31/00

          if ( contout >= u_min .and. contout <  hbreak ) then
               control_stat = 1
          else if ( contout >= hbreak .and. contout <= cbreak ) then
               control_stat = 2
          else if ( contout > cbreak .and. contout <= u_max ) then
               control_stat = 3
          end if


          temp_control: select case ( control_stat )

! *** Heating coil demand

             case(1)
                hdem = span(contout,u_min,hbreak,1.0, 0.0)
                cdem = 0.0
                ddem = 0.0

! *** Mixing damper demand

             case(2)
                hdem = 0.0
                cdem = 0.0
                ddem = span(contout,hbreak,cbreak,0.0,1.0)

! *** Cooling coil demand

             case(3)
                hdem = 0.0
                cdem = span(contout,cbreak,u_max,0.0,1.0)
                ddem = 0.0                                ! added 1/12/05

              case default
          end select temp_control

! *** Save provisional values to be used in updating at next sample instant

                saved_v(4) =  intp
                saved_v(6) =  difp
                saved_v(8) =  errp
                saved_v(10) = pidout

! *** Open loop mode (dampers set manually in type265)

            else
                ddem = 0.0
                cdem = cdemman
                hdem = hdemman
            endif

! *** Save current sample number
            saved_v(2) = real(nsample)

! *** Save outputs for use when not a sample instant
            saved_v(23) = ddem
            saved_v(24) = cdem
            saved_v(25) = hdem
        else

! *** Not a sample instant, set output to value from prev sample instant
            ddem = saved_v(23)
            cdem = saved_v(24)
            hdem = saved_v(25)

        endif

! *** Save time of current call

        saved_v(1) = time

! *** Output
        yout(1) = ddem
        yout(2) = cdem
        yout(3) = hdem

!        print *, 'ddem, cdem, hdem =',  ddem, cdem, hdem

! *** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo

! *** Return
        return
        end subroutine type492

! ***********************************************************************
! * SUBROUTINE:     E51 supply air temperature reset
! *
! * PURPOSE:        Reset supply air temperature setpoint so as to keep
! *                 the maximum of several room temperatures at a
! *                 specified value for cooling mode.  The minimum of
! *                 several room temperatures is used for heating mode.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. tz(1)   : temperature in zone 1                              (C)
! *  2. tz(2)   : temperature in zone 2                              (C)
! *  3. tz(3)   : temperature in zone 3                              (C)
! *  4. tz(4)   : temperature in zone 4                              (C)
! *  5. tz(5)   : temperature in zone 5                              (C)
! *  6. tz(6)   : temperature in zone 6                              (C)
! *  7. cdem    : cooling demand
! *  8. hdem    : heating demand
! *
! * OUTPUT
! * ======
! *  1. tsset   : supply air temperature setpoint                    (C)
! *
! * PARAMETERS
! * ==========
! *  1. tzset   : zone temperature setpoint                          (C)
! *  2. propb   : proportional band                                (K/K)
! *  3. tint    : integral time                                      (s)
! *  4. tssetmax: upper limit of the output                          (C)
! *  5. tssetmin: lower limit of the output                          (C)
! *  6. tsset0  : output at zero error (P), initial output (PI)      (C)
! *  7. deadb   : deadband                                           (K)
! *  8. ninputs : number of inputs                                   (-)
! *  9. tsamp   : sampling interval                                  (s)
! * 10. nfile   : controller number (parameters in file contN.par)   (s)
! *
! * SAVED
! * =====
! *  1. time    : time of previous call
! *  2. nsample : sample number of previous controller execution
! *  3. nsample : sample number of previous controller sample
! *  4. intp    : integral term from previous call
! *  5. intp    : integral term from previous sample
! *  6. tssetp  : output from previous call
! *  7. tssetp  : output from previous sample
! *  8-14.      : (PAR(1)-PAR(7) read from file
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:  No derivative action
!
!   DEVELOPER:           Li Mei and Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 22, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  DEADBAND
!   FUNCTIONS  CALLED:   BIGGEST, PIDCONT
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
!
!   Revised on August 31, 2000  Cheol Park, NIST
!    to include AHU heating coil.
!
!   Minor revision on September 21, 2001
!    to reduce the switching over between cooling and heating.
!
!   Updated : June 12, 2007  Cheol Park, NIST
!             Converted Fortran 77 code into Fortran 90
!
! **********************************************************************
!
        subroutine type493(xin,yout,par_v,saved_v,iostat)
!
        use modsim_head
        implicit none
        integer,parameter                 :: nsv=7,nfp=7
        integer,parameter                 :: ni=8,no=1,np=10,ns=nsv+nfp
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

! *** Declaration of local variables
        real,dimension(ni) :: tz
        real               :: intp
        real               :: cdem_min=0.001,hdem_min=0.001 ! 9/21/01
        real               :: tsset,hdem,cdem,tzset,propb,tint,tssetmax,&
                              tssetmin,tsset0,deadb,tssetp,tzmax,biggest,&
                              tzsetdb,dtsmax,dtsmin,dummy,pidcont,tzmin,&
                              smallest,tsamp
        integer            :: ipath,ninputs,nfile,i,nsample,is



        namelist /name493/time, ipath, cdem, hdem, tsset

! *** Read in parameters
        tzset      = par_v(1)
        propb      = par_v(2)
        tint       = par_v(3)
        tssetmax   = par_v(4)
        tssetmin   = par_v(5)
        tsset0     = par_v(6)
        deadb      = par_v(7)
        ninputs    = nint(par_v(8))
        tsamp      = par_v(9)
        nfile      = par_v(10)
! *** Read in inputs
        do i=1,ninputs
            tz(i) = xin(i)
        enddo
        cdem = xin(7)                    ! 8/31/00
        hdem = xin(8)

! *** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
                saved_v(2) = 0.
            endif
            if (init==0) then
! *** Initialize integral term and controller output
                saved_v(4) = 0.0
                saved_v(6) = tsset0
            endif
        endif
! *** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! *** Run controller if a sample instant
        nsample = nint(saved_v(3))+1
        if (time>=(nsample*tsamp) .or. time==0.0) then
! *** Use parameters from file if controller number not zero
            if (nfile>0) then
!  *** First call of time-step - read parameters from file and store
                if (time>saved_v(1)) then
                    call rfile(nfile,'cont',nfp,493,fpar)
                    do 150 i=1,nfp
                        saved_v(nsv+i)=fpar(i)
 150                continue
                endif
!  *** Get parameters that were read form file and then stored
                tzset      = saved_v(nsv+1)
                propb      = saved_v(nsv+2)
                tint       = saved_v(nsv+3)
                tssetmax   = saved_v(nsv+4)
                tssetmin   = saved_v(nsv+5)
                tsset0     = saved_v(nsv+6)
                deadb      = saved_v(nsv+7)
            endif
            if (time>saved_v(1)) then
! *** first call of timestep
! *** update previous sample instant values
                do is=4,(nsv-1),2
                    saved_v(is+1)=saved_v(is)
                enddo
            endif
! *** Update previous values
            intp   = saved_v(5)
            tssetp = saved_v(7)
            if(cdem > cdem_min) then                 ! 8/30/00
               ipath = 1
! *** Delect the highest room temperature
               tzmax = biggest(tz,ninputs)
! *** Apply dead-band
               call deadband(tzmax,tzset,deadb,tzsetdb)
! *** Pi control
! *** Setpoint can vary between tsset0+dtsmax and tsset0+dtsmin
               dtsmax = tssetmax-tsset0
               dtsmin = tssetmin-tsset0
! *** Pid output is offset of setpoint from tsset0
               tsset  = tsset0+pidcont(tzmax,tzsetdb,propb,tint,&
                  0.0,intp,dummy,tssetp,dummy,tsamp,dtsmax,dtsmin)
            else if (hdem > hdem_min) then          !  9/21/01
               ipath = 2
! *** Select the lowest room temperature
               tzmin = smallest(tz,ninputs)
! *** Apply dead-band
               call deadband(tzmin,tzset,deadb,tzsetdb)
! *** Pi control
! *** Setpoint can vary between tsset0+dtsmax and tsset0+dtsmin
               dtsmax = tssetmax-tsset0
               dtsmin = tssetmin-tsset0
! *** Pid output is offset of setpoint from tsset0
               tsset  = tsset0+pidcont(tzmin,tzsetdb,propb,tint,&
                  0.0,intp,dummy,tssetp,dummy,tsamp,dtsmax,dtsmin)
            else
               ipath = 3
               tsset = tssetp                    ! 9/21/01
            endif

! *** Save current sample number
            saved_v(2) = float(nsample)
! *** Save provisional values to be used in updating at next sample instant
            saved_v(4) = intp
            saved_v(6) = tsset
        else
! *** Not a sample instant, set output(s) to value(s) from previous sample
! *** instant
            tsset = saved_v(6)
        endif
! *** Save time of current call
        saved_v(1) = time
! *** Output
        yout(1) = tsset
! *** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo

!        print name493
! *** Return
        return
        end subroutine type493
! ***********************************************************************
! * SUBROUTINE:     E51 economizer control
! *
! * PURPOSE:        Determine the economizer control mode by comparing
! *                 outside and return air enthalpy
! *
!
!   Revised:       November 19, 2001   Cheol Park, NIST
!                  A deadband of enthalpy is included to minimizing
!                  switching of operation.
!   Revised:       Jan. 12, 2005
!                  Set .TRUE., if OA enthalpy is less than RA enthalpy.
!                      .FALSE., if OA enthalpy is greater than RA enthalpy.
!
! **********************************************************************
! * INPUTS
! * ======
! *  1. tout    : outdoor air temperature                            (C)
! *  2. hout    : outdoor air humidity ratio (0-1)                   (-)
! *  3. tret    : return air temperature                             (C)
! *  4. hret    : return air humidity ratio  (0-1)                   (-)
! *
! * OUTPUTS
! * =======
! *  1. econ    : economizer status (1 = OA enthalpy < RA enthalpy)  (-)! 01/12/05
! *
! * PARAMETERS
! * ==========
! *  1. tsamp   : sampling interval                                  (s)
! *  2. deadenth: deadband of enthapy                            (kJ/kg)   11/19/01
! *
! * SAVED
! * =====
! * 1.  time    : time of previous call of TYPE
! * 2.  nsample : sample number of previous controller execution
! * 3.  nsample : sample number of previous controller sample
! * 4.  econ    : controller output from previous sample
! *
! **********************************************************************
!
        subroutine type495(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=3,nfp=0
        integer,parameter                 :: ni=4,no=1,np=2,ns=nsv+nfp+1
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

        real            :: raenth,oaenth,hret,tret,hout,tout,tsamp,&
                           deadenth,fhair
        integer         :: nsample,i

! *** Declaration of local variables
        real         :: logtore
        logical      :: econ,retolog,compare

        namelist /econo/ time, tout,hout,tret,hret,oaenth, raenth,&
                         econ
! *** Read in inputs
        tout = xin(1)
        hout = xin(2)
        tret = xin(3)
        hret = xin(4)
! *** Read in parameter
        tsamp = par_v(1)
        deadenth = par_v(2)
! *** Initialize at beginning of simulation
        if (itime<=1) then
          if (init==0 .or. saved_v(1)>time) then
		saved_v(1) = -999999.
		saved_v(2) = 0.
	    endif
        endif
! *** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! *** Run controller if a sample instant
        nsample=nint(saved_v(3))+1
        if (time>=(nsample*tsamp) .or. time==0.0) then
! *** Run controller
! *** Calculate enthalpy of outside air
            oaenth = fhair(tout,hout)           ! 11/08/01
! *** Calculate enthalpy of return air
            raenth = fhair(tret,hret)           ! 11/08/01
! *** Determine economizer control signal
!            econ = compare(oaenth,raenth)

            if(oaenth >= raenth + deadenth) then              ! 01/12/05
               econ = .false.
            elseif(oaenth < raenth - deadenth) then
               econ = .true.
            else
               econ = retolog(saved_v(4))
            endif

! *** Save current sample number and output
	      saved_v(2) = float(nsample)
	      saved_v(4) = logtore(econ)
        else
! *** Not a sample instant, set output(s) to value(s) from previous sample
! *** instant
            econ = retolog(saved_v(4))
        endif
! *** Save time of current call
        saved_v(1) = time
! *** Output

!        if (econ) print econo

        yout(1) = logtore(econ)
! *** Allow freezing of algebraic variables
        do i=1,no
            iostat(i) = 1
        enddo
! *** Return
	return
	end subroutine type495

! ***********************************************************************
! * SUBROUTINE:    VAV room temperature control with reheat
! *
! * PURPOSE:       PI control of room temperature.
! *                Calculates the normalized demanded flow rate for a
! *                pressure-independent VAV box and the valve position
! *                for an associated reheat coil.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. tzon    : space temperature sensor                           (C)
! *  2. tsup    : supply air temperature sensor                      (C)
! *  3. sfstatus: supply fan status (1 = on, 0 = off)                (-)
! *  4. occupancy: building occupancy (1 = on, 0 = off)              (-)
! *
! * OUTPUTS
! * =======
! *  1. vdem    : normalised volumetric flow demand (0-1)            (-)
! *  2. rdem    : reheat coil demand (0-1)                           (-)
! *  3. contout : room demand (heating and cooling) (-1 - +1)        (-)
! *
! * PARAMETERS
! * ==========
! *  1. tsetheat: heating setpoint for occupied zone                 (C)
! *  2. tsetcool: cooling setpoint for occupied zone                 (C)
! *  3. tsetheatu: heating setpoint for unoccupied zone              (C)
! *  4. tsetcoolu: cooling setpoint for unoccupied zone              (C)
! *  5. propb   : proportional band                                  (K)
! *  6. tint    : integral time                                      (s)
! *  7. tdt     : derivative time                                    (s)
! *  8. rbreak  : breakpoint btwn damper and reheat coil signals (-1-+1)
! *  9. tdmin   : minimum turndown ratio                             (-)
! * 10. closloop: control mode (0=open loop, 1=closed loop)          (-)
! * 11. vdemman : open loop demanded normalized volumetric flow rate (-)
! * 12. rdemman : open loop reheat coil demand (0-1)                 (-)
! * 13. tsamp   : sampling interval                                  (s)
! * 14. nfile   : controller number (parameters in file contN.par)   (-)
! *
! * SAVED
! * =====
! *  1. time    : time of previous call
! *  2. nsample : sample number of previous controller execution
! *  3. nsample : sample number of previous controller sample
! *  4. intp    : integral term from previous call
! *  5. intp    : integral term from previous sample
! *  6. difp    : derivative term from previous call
! *  7. difp    : derivative term from previous sample
! *  8. errp    : error from previous call
! *  9. errp    : error from previous sample
! * 10. pidoutp : output from previous call
! * 11. pidoutp : output from previous sample
! * 12. suphot  : hysteresis output from previous call
! * 13. suphotp : hysteresis output from previous sample
! * 14-23.      : (PAR(1)-PAR(10) read from file
! * 24. vdem    : normalised volumetric flow set-point from previous call
! * 25. rdem    : reheat coil demand from previous call
! * 26. contout : room demand (heating and cooling) (-1 - +1)
! *
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 22, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  DEADBAND
!   FUNCTIONS  CALLED:   COMPHYS,RETOLOG,PIDCONT,LOGTORE,SWITCH,SPAN
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
!   Revicsed:            Cheol Park, NIST   Oct. 30, 2002
!                        To include setpoints of unoccupied periods.
!   Updated : June 12, 2007  Cheol Park, NIST
!             Converted Fortran 77 code into Fortran 90
!
! **********************************************************************
!
        subroutine type496(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=13,nfp=10
        integer,parameter                 :: ni=4,no=3,np=14,ns=nsv+nfp+3
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar
        real         :: tzon,tsup,tsetheat,tsetcool,tsetheatu,tsetcoolu,&
                        propb,tint,tdt,rbreak,tdmin,vdemman,rdemman,&
                        tsamp,difp,errp,pidoutp,tset,deadb,tsetdb,&
                        pidout,pidcont,switch,hour,vdem,span,rdem,contout
        integer      :: nfile,is,nsample,i

! *** Declaration of local variables
        real         :: logtore,intp
        logical      :: retolog,sfstatus,suphotp,suphot,comphys,closloop
        logical      :: occupancy

! *** Read in inputs
        tzon      = xin(1)
        tsup      = xin(2)
        sfstatus  = retolog(xin(3))
        occupancy = retolog(xin(4))

        if (occupancy) then
           sfstatus = .true.
        else
           sfstatus = .false.
        end if

! *** Read in parameters
        tsetheat  = par_v(1)
        tsetcool  = par_v(2)
        tsetheatu = par_v(3)
        tsetcoolu = par_v(4)
        propb     = par_v(5)
        tint      = par_v(6)
        tdt       = par_v(7)
        rbreak    = par_v(8)
        tdmin     = par_v(9)
        closloop  = retolog(par_v(10))
        vdemman   = par_v(11)
        rdemman   = par_v(12)
        tsamp     = par_v(13)
        nfile     = nint(par_v(14))
! *** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
                saved_v(2) = 0.
            endif
            if (init==0) then
                contout = 0.0          ! 12/21/99
                do is = 4,nsv-1,2
                    saved_v(is) = 0.0
                enddo
            endif
        endif
! *** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! *** Run controller if a sample instant
        nsample=nint(saved_v(3))+1
        if (time>=(nsample*tsamp) .or. time==0.0) then
! *** Run controller if mode is closed loop else use "manual" value
            if (closloop) then
! *** Closed loop mode
                if (time>saved_v(1)) then
! *** First call of timestep
! *** Update previous sample instant values
                    do is=4,nsv-1,2
                        saved_v(is+1)=saved_v(is)
                    enddo
                endif
! *** Update previous values
                intp    = saved_v(5)
                difp    = saved_v(7)
                errp    = saved_v(9)
                pidoutp = saved_v(11)
                suphotp = retolog(saved_v(13))

                if(occupancy) then
                   tset  = (tsetcool+tsetheat)/2.   ! occupied period  10/30/02
                   deadb = (tsetcool-tsetheat)/2.
                else
                   tset  = (tsetcoolu+tsetheatu)/2. ! unoccupied period 10/30/02
                   deadb = (tsetcoolu-tsetheatu)/2.
                endif

! *** Apply deadband around set-point
                call deadband(tzon,tset,deadb,tsetdb)
! *** Pid controller - output in the range -1 - +1
                pidout = pidcont(tzon,tsetdb,propb,tint,tdt,intp,difp,&
                                 pidoutp,errp,tsamp,1.0,-1.0)
! *** Set output to zero if fan status not ok
                contout = switch(sfstatus,pidout,0.0)

!   Revised code         11/16/01
                hour = time/3600

                if(contout <= rbreak) then
                     vdem = span(contout,rbreak,-1.0,tdmin,1.0)
                     rdem = 0.0
                else
                     vdem = tdmin
                     rdem = span(contout,rbreak,1.0,0.0,1.0)
                endif

! *** Save provisional values to be used in updating at next sample instant
                saved_v(4)  = intp
                saved_v(6)  = difp
                saved_v(8)  = errp
                saved_v(10) = pidout
                saved_v(12) = logtore(suphot)
! *** Open loop mode - set manually
            else
                vdem = vdemman
                rdem = rdemman
            endif
! *** Save current sample number
            saved_v(2) = float(nsample)
! *** Save outputs for use when not a sample instant
            saved_v(24) = vdem
            saved_v(25) = rdem
            saved_v(26) = contout
        else
! *** Not a sample instant, set output to value from prev sample instant
            vdem    = saved_v(24)
            rdem    = saved_v(25)
            contout = saved_v(26)
        endif
! *** Save time of current call
        saved_v(1) = time
! *** Output
        yout(1) = vdem
        yout(2) = rdem
        yout(3) = contout

! *** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo
! *** Return
        return
        end subroutine type496
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     E51 modulated mixed air damper control
! *
! * PURPOSE:        Determines the mixed air damper control signal.
! *                 The damper positions are determined by a combination
! *                 of the supply air temperature controller (TYPE486)
! *                 and proportional control of the mixed air
! *                 temperature. If the supply fan status is not OK or
! *                 the outside enthalpy is higher than the return
! *                 enthalpy or low temperature over-ride is set,
! *                 the mixing box is set to full recirculation.
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. tmix    : mixed air temperature sensor                       (C)
! *  2. mbcdem  : mixing box cooling demand (0-1)                    (-)
! *  3. sfstatus: supply fan status (1 = on, 0 = off)                (-)
! *  4. lowtovr : low temperature over-ride (1 = low temperature)    (-)
! *  5. econ    : economizer status (1 = OA enthalpy > RA enthalpy)  (-)
! *  6. occupancy: building occupancy (1 = on, 0 = off)              (-)      ! 11/01/02
! *
! * OUTPUTS
! * =======
! *  1. ddemo   : outside air damper demanded position (0-1)         (-)
! *  2. ddemr   : return air damper demanded position (0-1)          (-)
! *  3. ddeme   : exhaust air damper demanded position (0-1)         (-)
! *
! * PARAMETERS
! * ==========
! *  1. tset    : mixed air temperature setpoint                     (C)
! *  2. propb   : proportional band                                  (K)
! *  3. tint    : integral time                                      (s)
! *  4. tdt     : derivative time                                    (s)
! *  5. deadb   : deadband                                           (K)
! *  6. closloop: control mode (0 = open loop, 1 = closed loop)      (-)
! *  7. ddemmano: open loop outside air damper position (0-1)        (-)
! *  8. ddemmanr: open loop return air damper position (0-1)         (-)
! *  9. ddemmane: open loop exhaust air damper position (0-1)        (-)
! * 10. tsamp   : sampling interval                                  (s)
! * 11. nfile   : controller number (parameters in file contN.par)   (-)
! *
! * SAVED
! * =====
! *  1. time    : time of previous call
! *  2. nsample : sample number of previous controller execution
! *  3. nsample : sample number of previous controller sample
! *  4. intp    : integral term from previous call
! *  5. intp    : integral term from previous sample
! *  6. difp    : derivative term from previous call
! *  7. difp    : derivative term from previous sample
! *  8. errp    : error from previous call
! *  9. errp    : error from previous sample
! * 10. pidoutp : output from previous call
! * 11. pidoutp : output from previous sample
! * 12-20.      : (PAR(1)-PAR(9) read from file
! * 21. ddemo   : Outside air damper demanded from previous call
! * 22. ddemr   : Return air damper demanded from previous call
! * 23. ddeme   : Exhaust air damper demanded from previous call
! **********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Li Mei and Philip Haves
!                        Loughborough University of Technology
!
!   DATE:                November 22, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   RETOLOG,LOGICOR,LOGICNOT,LOGICAND,SWITCH,SPAN
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
!
!   revised:             November 1, 2002 Cheol Park
!                        To inculde unoccupied operation.
!
! **********************************************************************

        subroutine type497(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=11,ni=6,no=3,np=12,&
                                             nfp=9,ns=nsv+nfp+3
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

! **** Declaration of local variables
        real         :: mbcdem, intp           ! changed 12/6/1999
        logical      :: retolog,logicor,logicnot,logicand
        logical      :: sfstatus,lowtovr,econ,closloop,orloeo,norloeo,&
                        auxdis
        real         :: tmix,tset,propb,tint,tdt,deadb,ddemmano,&
                        ddemmanr,ddemmane,tsamp,difp,errp,pidoutp,&
                        tsetdb,pidout,pidcont,passthru,switch,ddemo,&
                        ddemr,ddeme,ddemo_temp,ddemomin
        integer      :: nfile,is,nsample,i
        logical      :: occupancy
        integer      :: itype=497

! **** Read in inputs
        tmix        = xin(1)
        mbcdem      = xin(2)
        sfstatus    = retolog(xin(3))
        lowtovr     = retolog(xin(4))
        econ        = retolog(xin(5))
        occupancy   = retolog(xin(6))

! **** Read in parameters
        tset        = par_v(1)
        propb       = par_v(2)
        tint        = par_v(3)
        tdt         = par_v(4)
        deadb       = par_v(5)
        closloop    = retolog(par_v(6))
        ddemmano    = par_v(7)
        ddemmanr    = par_v(8)
        ddemmane    = par_v(9)
        tsamp       = par_v(10)
        nfile       = par_v(11)
	ddemomin    = par_v(12)

! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
                saved_v(2) = 0.
            endif
            if (init==0) then
                do is = 4,(nsv-1),2
                    saved_v(is) = 0.0
                enddo
            endif
        endif
! **** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! **** Run controller if a sample instant
        nsample=nint(saved_v(3))+1
        if (time>=(nsample*tsamp) .or. time==0.0) then
! **** Read parameters from file if controller number not zero
            if (nfile>0) then
! **** First call of time-step - read parameters from file and store
                if (time>saved_v(1)) then
                    call rfile(nfile,'cont',nfp,itype,fpar)
                    do i=1,nfp
                        saved_v(nsv+i)=fpar(i)
                    enddo
                endif
! **** Get parameters that were read form file and then stored
                tset     = saved_v(nsv+1)
                propb    = saved_v(nsv+2)
                tint     = saved_v(nsv+3)
                tdt      = saved_v(nsv+4)
                deadb    = saved_v(nsv+5)
                closloop = retolog(saved_v(nsv+6))
                ddemmano = saved_v(nsv+7)
                ddemmanr = saved_v(nsv+8)
                ddemmane = saved_v(nsv+9)
            endif
            if (closloop) then
! **** Closed loop mode
! **** First call of timestep - update previous sample instant values
                if (time>saved_v(1)) then
                    do is = 4,nsv-1,2
                        saved_v(is+1) = saved_v(is)
                    enddo
                endif
! **** Update previous values
                intp    = saved_v(5)
                difp    = saved_v(7)
                errp    = saved_v(9)
                pidoutp = saved_v(11)
! **** Apply deadband around set-point
                call deadband(tmix,tset,deadb,tsetdb)
! **** Pid controller
                pidout = pidcont(tmix,tsetdb,propb,tint,tdt,intp,difp,&
                                 pidoutp,errp,tsamp,2.0,0.0)
! **** Determine whether auxiliary input switch is enabled or disabled
! **** not low temperature over-ride or economizer
                orloeo = logicor(lowtovr,econ)
                norloeo = logicnot(orloeo)
! **** Auxiliary input disabled if supply fan status ok and neither
! **** low temperature over-ride or economizer is set
                auxdis = logicand(sfstatus,norloeo)
! **** Determine mixing box position demanded by supply air temperature
! **** controller - pass input through if
! **** auxilliary input disabled else set to zero (full recirculation)
                passthru = switch(auxdis,mbcdem,0.0)
! **** Determine demanded mixing box position - lesser of the mixed air
! **** temperature controller output and the (modoified) signal from the
! **** supply air temperature controller
!                ddemo = min(pidout,passthru)
! **** Return and exhaust air dampers demanded position
!                ddemr = 1.0 - ddemo
!                ddeme = ddemo

! ****** MODIFICATION to eliminate P control of mixed air temperature ******
! 
! pidout is set to a large number so that ddemo_temp will always use the
! value of passthru
!
                pidout = 999.99      ! changed 7/18/01 (jh)

                ddemo_temp = min(pidout,passthru)

! if passthru = 1.0, ddemo = 1.0
!                    ddemr = 0.0
!                    ddeme = 1.0

! if passtrhu = 0.0, ddemo = ddemomin = 0.15   ! assuming 15% oa
!                    ddemr = 0.85
!                    ddeme = 0.15

                if (occupancy) then            ! building occupied       11/1/02
		     ddemo = ddemo_temp * (1.0 - ddemomin) + ddemomin ! ddemomin = minimum oa damper signal
                else
		     ddemo = ddemomin          ! building unoccupied      11/1/02
                end if
		ddemr = 1.0 - ddemo
		ddeme = ddemo

! **** save provisional values to be used in updating at next sample instant
                saved_v(4) = intp
                saved_v(6) = difp
                saved_v(8) = errp
                saved_v(10) = pidout
! **** open loop mode
            else
                ddemo = ddemmano
                ddemr = ddemmanr
                ddeme = ddemmane
            endif

! **** Save current sample number
            saved_v(2)  = float(nsample)
! **** Save provisional values to be used in updating at next sample instant
            saved_v(21) = ddemo
            saved_v(22) = ddemr
            saved_v(23) = ddeme
        else
! ****Not a sample instant, set output to value from prev sample instant
            ddemo = saved_v(21)
            ddemr = saved_v(22)
            ddeme = saved_v(23)
        endif
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = ddemo
        yout(2) = ddemr
        yout(3) = ddeme
! **** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo
! **** Return
        return
        end subroutine type497

! ***********************************************************************
!
! TYPE499.for:     Air damper control
!
!   Created:  January 13, 2005  C.P., NIST
!   Updated : June 12, 2007  Cheol Park, NIST
!             Converted Fortran 77 code into Fortran 90
!
! *********************************************************************
!  INPUTS
!  ======
!   1. ddem    : damper demand (0-1)                                (-)
!   2. econ    : economizer status(TRUE: OA enthalpy < RA enthalpy) (-)
!   3. occupancy: building occupancy (1 = on, 0 = off)              (-)      ! 11/01/02
!   4. cdem    : cooling coil valve demand (0-1)                    (-)
!
!  OUTPUTS
!  =======
!   1. ddemo   : outside air damper demanded position (0-1)         (-)
!   2. ddemr   : return air damper demanded position (0-1)          (-)
!   3. ddeme   : exhaust air damper demanded position (0-1)         (-)
!
!  PARAMETERS
!  ==========
!   1. closloop: control mode (0 = open loop, 1 = closed loop)      (-)
!   2. ddemmano: open loop outside air damper position (0-1)        (-)
!   3. ddemmanr: open loop return air damper position (0-1)         (-)
!   4. ddemmane: open loop exhaust air damper position (0-1)        (-)
!   5. ddemomin: minimum OA damper ventilation position             (-)
!
!  SAVED
!  =====
!   1. time    : time of previous call
!   2. ddemo   : Outside air damper demanded from previous call
!   3. ddemr   : Return air damper demanded from previous call
!   4. ddeme   : Exhaust air damper demanded from previous call
! *********************************************************************
!
        subroutine type499(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ni=4,no=3,np=5,ns=4
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

! *** Declaration of local variables
        logical      :: retolog,econ,closloop,occupancy
        real         :: ddemo,ddemomin,ddem,cdem,ddemmano,ddemmanr,&
                        ddemmane,ddemr,ddeme
        integer      :: is,i

        namelist /ty498/time, econ, ddemomin, ddemo
                                        ! Read in inputs
        ddem        = xin(1)
        econ        = retolog(xin(2))
        occupancy   = retolog(xin(3))
        cdem        = xin(4)
                                        ! Read in parameters
        closloop    = retolog(par_v(1))
        ddemmano    = par_v(2)
        ddemmanr    = par_v(3)
        ddemmane    = par_v(4)
	  ddemomin    = par_v(5)
                                        ! Initialize at beginning of simulation
        if (itime <= 1) then
            if (init == 0 .or. saved_v(1) > time) then
                saved_v(1) = -999999.
            endif
            if (init==0) then
                do is = 2,ns
                    saved_v(is) = 0.0
                enddo
            endif
        endif
        if (time >= saved_v(1) .or. time == 0.0) then
                                        ! Closed loop mode
            if (closloop) then
                if (occupancy) then     ! Building occupied
		     ddemo = ddem * (1.0 - ddemomin) + ddemomin
                     if(cdem > 0.0) then
                        if(econ) then
                           ddemo = 1.0
                        else
                           ddemo = ddemomin
                        endif
                     endif
                else
		     ddemo = ddemomin   ! Building unoccupied
                end if
 		ddemr = 1.0 - ddemo
 		ddeme = ddemo
                                        ! Open loop mode
           else
               ddemo = ddemmano
               ddemr = ddemmanr
               ddeme = ddemmane
           endif
                                        ! Save provisional values
           saved_v(2) = ddemo
           saved_v(3) = ddemr
           saved_v(4) = ddeme
        else
                                        ! If not a sample instant,
                                        ! set output to previously saed  value
           ddemo = saved_v(2)
           ddemr = saved_v(3)
           ddeme = saved_v(4)
        endif
                                        ! Save time of current call
        saved_v(1) = time
                                        ! Output

        yout(1) = ddemo
        yout(2) = ddemr
        yout(3) = ddeme
                                        ! Disallow freezing
        do i=1,no
           iostat(i) = 0
        enddo

        return
        end subroutine type499
! ***********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! ***********************************************************************
! * SUBROUTINE:  Read inputs from a file
! *
! * PURPOSE:     Get control signals from another program (e.g. a text
! *              editor) by reading from an intermediate file every step time
! *
! ***********************************************************************
! *
! * INPUTS
! * ======
! *  1. (dummy)
! *
! * OUTPUTS
! * =======
! *  1.         : simulation input  1
! *  2.         : simulation input  2
! *  3.         : simulation input  3
! *  4.         : simulation input  4
! *  5.         : simulation input  5
! *  6.         : simulation input  6
! *  7.         : simulation input  7
! *  8.         : simulation input  8
! *  9.         : simulation input  9
! * 10.         : simulation input 10
! * 11.         : simulation input 11
! * 12.         : simulation input 12
! * 13.         : simulation input 13
! * 14.         : simulation input 14
! * 15.         : simulation input 15
! * 16.         : simulation input 16
! *
! * PARAMETERS
! * ==========
! *  1. nfile   : file number (FILE='inputN.par', N > 0)             (-)
! *  2. tsamp   : sample time (interval between reads)               (s)
! *  3. realt   : real time scaling factor (0=no wait, 1=real time)  (-)
! *  4. itext   : text output to screen (0=no, 1=yes)                (-)
! *  5. noutx   : number of values to read                           (-)
! *
! * SAVED
! * =====
! *  1. time     : time of previous call of TYPE
! *  2. nsample  : sample number of previous read
! *  3-18        : file data from previous sample
! *
! ***********************************************************************
!
!   MAJOR RESTRICTIONS:  First order, linear
!
!   DEVELOPER:           Philip Haves
!                        Loughborough University of Technology
!
!   REVISION DATE:       May 19, 1996
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:  RFILE, DATETIME, TCORRECT
!   FUNCTIONS  CALLED:   ELAPSED
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! ***********************************************************************

       subroutine type504(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: ndiffeq=0,ni=1,no=16,np=5,&
                                             ns=2+no
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(no)                :: store
        character(len=30)                 :: timstamp

        real         :: tsamp,realt,rtime,stime,etime,rlate,elapsed
        integer      :: irtimes0,irtimem0,nfile,itext,nsample,irtimes,&
                        irtimem,i,is,noutx
        integer      :: itype=504

! **** Save starting time in seconds and milliseconds (must be stored as
! **** integers)
        save irtimes0,irtimem0
! **** Read in inputs
! **** Input is dummy - no need to read
! **** Read in parameters
        nfile = nint(par_v(1))
        tsamp = par_v(2)
        realt = par_v(3)
        itext = par_v(4)
        noutx = nint(par_v(5))

        stime = time                           ! ??????????
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>stime) then
                saved_v(1) = -999999.
! **** Initialize whole saved_v array to avoid spurious values for
! **** unused outputs
                do is=2,ns
                    saved_v(is) = 0.
                enddo
                do i=1,no
                   store(i) = 0.0            ! 12/21/99
                enddo
            endif
        endif
! **** Accelerated simulation time
        rtime=realt*stime
! **** Read file if a sample instant
        nsample=nint(saved_v(2))
        if (stime>=(nsample*tsamp)) then
! **** Sample instant
            if (stime>saved_v(1)) then
! **** First call of timestep - read file
                if(realt>0.0) then
! **** Real time - read clock and wait if necessary
! **** check for first step time of simulation
                    if (itime<=1) then
! **** First step time
! **** Read real-time clock and get starting time (seconds and milliseconds)
                        call datetime(irtimes0,irtimem0,timstamp)
! **** Correct for first timestep
                        call tcorrect(irtimes0,irtimem0,rtime)
                        if (itext==1) write(*,900) timstamp
                    else
! **** Not first timestep - wait for real time if synchronized, then read file
! **** update previous values
                        call datetime(irtimes,irtimem,timstamp)
                        etime=elapsed(irtimes,irtimem,irtimes0,irtimem0)
                        rlate=etime-rtime
                        if (rlate>0.0) then
! **** Late - write time to screen and procede
                            write(*,901) stime,rlate
                        else
! **** Early, wait in loop
                            if (itext==1) then
                                write(*,902) stime,-rlate
                            endif
 100                        call datetime(irtimes,irtimem,timstamp)
                            etime=elapsed(irtimes,irtimem,irtimes0,irtimem0)
                            if (rtime>etime) goto 100
                        endif
                    endif
                endif
! **** Get input from file
                call rfile(nfile,'input',noutx,itype,store)
                do i=1,no
                    yout(i) = store(i)
                enddo
! **** Write values to screen if reqired
                if (itext==1) then
                    write(*,903) timstamp
                    write(*,910) (yout(i),i=1,8)
                    if (noutx>=9) then
                        write(*,911) (yout(i),i=9,16)
                    endif
                endif
! **** Update number of next sample
                saved_v(2)=real(nint(stime/tsamp) + 1)
! **** Update previous outputs
                do i=1,noutx
                    saved_v(i+2)=yout(i)
                enddo
            else
! **** Not first iteration - use previous values
                do i=1,no
                    yout(i)=saved_v(2+i)
                enddo
            endif
        else
! **** Not a sample instant - use previous values
            do i=1,no
                yout(i)=saved_v(2+i)
            enddo
        endif
! **** Save time of current call
        saved_v(1) = stime
! **** Disallow freezing (outputs read from file may change)
        do i=1,no
            iostat(i) = 1
        enddo
! **** Return
        return
! **** Format statements
900     format('T504: starting time = ',a30)
901     format('T504: simulation time =',f10.2,&
               ' simulation is ',f6.2,' late')
902     format('T504: simulation time =',f10.2,&
               ', spare time =',f6.2)
903     format(a30)
910     format('T504: vars read (1-8) :',8f7.3)
911     format('T504: vars read (9-16):',8f7.3)
        end subroutine type504

! ***********************************************************************
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
! * 10. occupancy: Occupnacy (1/0)                                   (-)
! *
! * OUTPuts
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
! * PARAmeters
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
! * 24. coil_id
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
! ***********************************************************************
!
!   Revised:  Dec. 10, 2002  Cheol Park, NIST
!   Updated : June 12, 2007  Cheol Park, NIST
!             Converted Fortran 77 code into Fortran 90
!
! **********************************************************************
!
        subroutine type530(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        use cool_plant
        implicit none
        integer,parameter                 :: ndiffeq=1
        integer,parameter                 :: ni=10,no=15,np=24,&
                                             ns=1+ndiffeq*6
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real          :: axo,ap,asx,ditube,ai,afree,aext,aratio,bratio,&
                         rifin,rofin,heifin,aflow,hydiam,confin,contube,&
                         vcpfin,vcptube,vstembyp,frvalflo,frvalbyp,&
                         frtotflo,frtotbyp,dp,ww,wwbypas,frprim,&
                         frtotcoi,frcoibyp,frtot,dpprim,dpcoibyp,pai,&
                         ta_fire,tsp,aap,bbp,twbo,rho,ho,bf,shr,bb,&
                         aa,aabar,bbbar,tsbar,effect,&
                         qa,wwprim,tret,two,twi,vstem,rlinport,fwphi,&
                         tai,gi,gsat,pao,wa,pwi,pwo,tsdyn,hcoil,wcoil,&
                         dotube,thitube,spafin,thifin,fra,frwcoil,&
                         frwbypas,eqpchar,sv,cl,ts,tao,go,&
                         fracfin,acgross,ttubelen,reqpport,wqud,dpqud,&
                         dpqudlin
        real          :: lcoil,kv
        integer       :: is,i,ifault,nrow,ntpr,ncir,icountfl,&
                         ntube,matube,mafin,coil_id,ivaltype
        logical       :: dynamic,psychro

                  ! thermophysical constants
        real,dimension(4)   :: volcpm = (/2455.,3475.,3915.,2250./),&
                               conma  = (/0.221,0.393,0.0453,0.0013/)

                  ! miscellaneous constants (pi, atmospheric pressure, minimum
                  ! water-side resistance
        real                :: pi = 3.14159, patmpa = 101325.,&
                               frsmall = 1.0e-8

                  ! break-points in valve characteristics (see documentation for
                  ! type329 two port valve)
        real                :: xlin = 0.01, xeqp = 0.33

        logical             :: occupancy, retolog
        
        namelist /coil_data1/time, coil_id, vstem, twi, two, tret,&
          wwprim, qa
        namelist /coil_data2/xin, yout

        if(time > start_chiller .and. time < stop_chiller) then
           cooling = .true.
        else
           cooling = .false.
        endif
                         ! Read in inputs
        tai     = xin(1)
        gi      = xin(2)
                         ! Limit humidity ratio to saturated value
        gsat    = fwphi(tai,100.0,patmpa)
        gi      = min(gi,gsat)
        pao     = xin(3)
        wa      = xin(4)
        twi     = xin(5)
        pwi     = xin(6)
        pwo     = xin(7)
        vstem   = xin(8)
        tsdyn   = xin(9)
        occupancy = retolog(xin(10))         !  12/4/02

                          ! Read in parameters
                          ! Dynamic or steady state
        if (nint(par_v(1))==1) then
            dynamic = .true.
        else
            dynamic = .false.
        endif
                          ! Fault type
        ifault = nint(par_v(2))
                          ! Peform additional pyschrometric calculations?
        if (nint(par_v(3))/=0) then
            psychro = .true.
        else
            psychro = .false.
        endif
                          ! Coil geometry and materials
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
                          ! Air and water flow resistances
        fra     = par_v(16)
        frwcoil = par_v(17)
        frwbypas= par_v(18)
                          ! Valve characteristics
        ivaltype= nint(par_v(19))                            ! 4/24/07
        kv      = par_v(20)
        eqpchar = par_v(21)
        sv      = par_v(22)
        cl      = par_v(23)
        coil_id = nint(par_v(24))                             ! added 11/18/02

        if(coil_id == 2 .or. coil_id == 7 .or. coil_id == 12) then  ! <<<<< hard code <<<<<
           if (.not. cooling) then
              ts  = tai
              tao = tai
              go  = gi
              two = tai                        ! 12/17/02
!              two = twi
              tret= twi
              go to 10
           else
              if(twi > tai) then               ! 12/18/02
                 twi = tai
                 two = tai
              endif
           end if

        endif

        if(.not. occupancy) then                ! occupancy 12/4/02
           ts  = tai
           tao = tai
           go  = gi                             !  12/17/02
           two = tai
!           two = twi
           tret= twi
           goto 10
        endif

                          ! Heat exchanger geometry
        if (nrow<=2) then
                          ! One or two rows -> cross flow
            icountfl=0
        elseif (nrow>2) then
            if (ifault/=1) then
                          ! Three or more rows and no connection fault ->
                          ! counterflow
                icountfl=1
            else
                          ! Three or more rows and connection fault ->
                          ! parallel flow
                icountfl=-1
            endif
        endif

                          ! Heat transfer areas
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
                          ! Valve types and resistances
                          ! Limit stem position to range 0-1
        vstem   = max( 0.0, min(1.0,vstem) )
        vstembyp= 1.-vstem
                          ! Water flow rates through coil and bypass
        if (ivaltype>=0 .and. ivaltype<=2) then
                          ! Common port carries common flow
                          ! Flow port - Select valve type and calculate
                          ! resistance
            if (ivaltype==0 .or. ivaltype==2) then
                          ! Linear characteristic - two segments
                          ! (linear and close-off)
                frvalflo=rlinport(vstem,kv,sv,xlin,cl)
            elseif (ivaltype==1) then
                          ! Equal percentage characteristic -
                          ! three segments (exponential, linear and close-off)
                frvalflo=reqpport(vstem,kv,xeqp,eqpchar,xlin,sv,cl)
            endif
                          ! By-pass port - Select valve type and calculate
                          ! resistance
            if (ivaltype==0 .or. ivaltype==1) then
                          ! Linear characteristic - two segments
                          ! (linear and close-off)
                frvalbyp=rlinport(vstembyp,kv,sv,xlin,cl)
            else
                          ! Equal percentage characteristic -
                          ! three segments (exponential, linear and close-off)
                frvalbyp=reqpport(vstembyp,kv,xeqp,eqpchar,xlin,sv,cl)
            endif
                          ! Total water-side resistances - coil + valve,
                          ! by-pass + valve
            frtotflo  = frwcoil+frvalflo
            frtotbyp  = frwbypas+frvalbyp
            if (frtotflo<frsmall .or. frtotbyp<frsmall) then
                          ! total resistance almost zero
                write(*,*)&
                'type530: water-side flow resistance must not be < ',&
                frsmall
                stop
            endif
                          ! Pressure drop
            dp      = pwi-pwo
                          ! Flow rate through coil
            ww      = wqud(frtotflo,dp)
                          ! Flow rate through by-pass
            wwbypas = wqud(frtotbyp,dp)
                          ! Flow rate in primary circuit
            wwprim  = ww+wwbypas
        elseif (ivaltype==3 .or. ivaltype==4) then
                          ! Faulty installation
            if (ivaltype==3) then
                          ! Flow and common ports swopped
                          ! Resistance of flow port - equal percentage
                          ! characteristic
                frprim   = reqpport(vstem,kv,xeqp,eqpchar,xlin,sv,cl)
                          ! Resistance of bypass port - linear characteristic
                frvalbyp = rlinport(vstembyp,kv,sv,xlin,cl)
                          ! Resistance of coil branch (common port has zero
                          ! resistance)
                frtotcoi = frwcoil
                          ! Resistance of bypass branch (balancing valve +
                          ! bypass port)
                frtotbyp = frwbypas + frvalbyp
                          ! Parallel resistance of coil and bypass branches
                frcoibyp = (frtotcoi**(-0.5)+frtotbyp**(-0.5))**2
                          ! Total resistance of subsystem
                frtot    = frcoibyp + frprim
                          ! Pressure drop across subsystem
                dp       = pwi-pwo
                          ! Flow rate in primary circuit
                wwprim   = wqud(frtot,dp)
                dpprim   = dpqud(frprim,wwprim)
                dpcoibyp = dp - dpprim
                          ! Flow rate through coil
                ww       = wqud(frtotcoi,dpcoibyp)
                          ! Flow rate through by-pass
                wwbypas  = wqud(frtotbyp,dpcoibyp)
            elseif (ivaltype==4) then
                          ! Bypass and common ports swopped
                          ! resistance of flow port - equal percentage
                          ! characteristic
                frvalflo = reqpport(vstem,kv,xeqp,eqpchar,xlin,sv,cl)
                          ! Resistance of bypass port - linear characteristic
                frprim   = rlinport(vstembyp,kv,sv,xlin,cl)
                          ! Resistance of coil branch (flow port + coil)
                frtotcoi = frvalflo + frwcoil
                          ! Resistance of bypass branch (common port has zero
                          ! resistance)
                frtotbyp = frwbypas
                          ! Parallel resistance of coil and bypass branches
                frcoibyp = (frtotcoi**(-0.5)+frtotbyp**(-0.5))**2
                          ! Total resistance of subsystem
                frtot    = frcoibyp + frprim
                          ! Pressure drop across subsystem
                dp       = pwi-pwo
                          ! Flow rate in primary circuit
                wwprim   = wqud(frtot,dp)
                dpprim   = dpqud(frprim,wwprim)
                dpcoibyp = dp - dpprim
                          ! Flow rate through coil
                ww       = wqud(frtotcoi,dpcoibyp)
                          ! Flow rate through by-pass
                wwbypas  = wqud(frtotbyp,dpcoibyp)
            endif
        else
            stop 'valve type out of range'
        endif
                          ! Calculate air inlet pressure from flow resistance
        pai     = pao+dpqudlin(fra,wa)

                          ! Bypass the calculation of coil duty and outlet
                          ! conditions, if the inlet air temperature is out
                          ! of range due to fire.  9/30/98

     ta_fire = 50.0                     ! for CFAST
     if (tai < ta_fire) then

                          ! Calculate coil duty and outlet conditions
        if (dynamic) then
                          ! Initialize at beginning of simulation
            if (itime<=1) then
                if (init==0 .or. saved_v(1)>time) then
                    saved_v(1) = -999999.
               endif
                if (init==0) then
                          ! Capacitance of coil initially at inlet air
                          ! temperature (coil "off")
                    saved_v(2) = tai
                    saved_v(4) = 0.0
                    saved_v(6) = 0.0
                endif
            endif
                          ! Every time-step
            if (time>saved_v(1)) then
                          ! First call of time-step - update value of
                          ! temperature rise, aa and bb from previous time-step
                do is = 2,ns-1,2
                    saved_v(is+1) = saved_v(is)
                enddo
            endif
                          ! Update previous values
            tsp = saved_v(3)
            aap = saved_v(5)
            bbp = saved_v(7)
            call licoilab(tai,gi,wa,twi,ww,tsdyn,&
                          psychro,icountfl,lcoil,wcoil,dotube,thitube,&
                          thifin,ttubelen,asx,ditube,afree,aext,&  ! MAG change as to asx
                          aratio,bratio,rifin,rofin,heifin,aflow,&
                          hydiam,confin,contube,vcpfin,vcptube,&
                          aa,bb,tao,go,two,qa,shr,bf,ho,rho,twbo)
                          ! Integrate differential equation for surface
                          ! temperature
            aabar = (aa + aap)/2.0
            bbbar = (bb + bbp)/2.0
            call diffeq(time,aabar,bbbar,tsp,ts,tsbar)
                          ! Save time of current call
            saved_v(1) = time
                          ! Save provisional value to be used in updating at
                          ! next step time
            saved_v(2) = ts
            saved_v(4) = aa
            saved_v(6) = bb
        else
                          ! Steady state model - output is value of surface
                          ! temperature
            call licoilss(tai,gi,wa,twi,ww,&
                          psychro,icountfl,lcoil,wcoil,thitube,&
                          thifin,ditube,afree,aext,&
                          aratio,bratio,rifin,rofin,heifin,aflow,&
                          hydiam,confin,contube,&
                          ts,tao,go,two,qa,shr,effect,bf,ho,rho,twbo)
        endif
                          ! Return water temperature
        if (wwprim>1.0e-10) then
            tret    = (ww*two + (wwprim-ww)*twi)/wwprim
        else
            tret    = two
        endif
     else

       ts  = tai                   !  tai > 50.0
       tao = tai
       go  = gi
       two = twi
       tret= twi

     end if
                          ! Assign output values
   10 continue

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
                        ! Effectiveness not defined for dynamic case
          yout(11)  = -99.9
      else
          yout(11)  = effect
      endif

      yout(12) = bf

      if (psychro) then
                        ! Auxilliary outputs
          yout(13) = ho
          yout(14) = rho/100.
          yout(15) = twbo
      else
                        ! Auxilliary outputs not calculated - set values to
                        ! flag this
          yout(13) = -99.9
          yout(14) = -99.9
          yout(15) = -99.9
      endif
                        ! Allow freezing of algebraic variables
      if (.not.dynamic) then
          iostat(1)=1
      else
          iostat(1)=0
      endif
      do i=2,no
          iostat(i)=1
      enddo

      if (coil_id == 2 .and. time > 234000.0 ) then                ! <<<<< hard code <<<<<<<
!      if (coil_id == 2 .or. coil_id == 7 .or. coil_id == 12) then
!         print coil_data1
!         print coil_data2
      end if

        return
        end subroutine type530

