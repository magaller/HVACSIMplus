! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     PID controller
! *
! * PURPOSE:        PID controller
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. y       : controlled variable (sensor signal)                (-)
! *  2. w       : setpoint                                           (-)
! *
! * OUTPUTS
! * =======
! *  1. u       : control signal                                     (-)
! *
! * PARAMETERS
! * ==========
! *  1. propb   : proportional band                                  (-)
! *  2. tint    : integral time                                      (s)
! *  3. tdt     : derivative time                                    (s)
! *  4. dband   : deadband                                           (-)
! *  5. closloop: control mode (0=open loop, 1=closed loop)          (-)
! *  6. uman    : open loop control signal                           (-)
! *  7. umin    : lower limit on control signal                      (-)
! *  8. umax    : upper limit on control signal                      (-)
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
! * 20. up      : controller output from previous sample
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
!   FUNCTIONS  CALLED:   RETOLOG, PIDCONT
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
!
!   Updated to Fortran 90 April 19, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type441(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=11,ni=2,no=1,np=10,&
                                             nfp=8,ns=nsv+nfp+1
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

! **** Declaration of local variables
        real         :: intp
        logical      :: retolog,closloop
        real         :: y,w,propb,tint,tdt,dband,uman,umin,umax,tsamp,&
                        difp,errp,up,wdb,pidcont,ux
        integer      :: nfile,is,nsample,i
        integer      :: itype=441

! **** Read in inputs
        y         = xin(1)
        w         = xin(2)
! **** Read in parameters
        propb     = par_v(1)
        tint      = par_v(2)
        tdt       = par_v(3)
        dband     = par_v(4)
        closloop  = retolog(par_v(5))
        uman      = par_v(6)
        umin      = par_v(7)
        umax      = par_v(8)
        tsamp     = par_v(9)
        nfile     = par_v(10)
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
! **** Use parameters from file if controller number not zero
            if (nfile>0) then
! **** First call of time-step - read parameters from file and store
                if (time>saved_v(1)) then
                    call rfile(nfile,'cont',nfp,itype,fpar)
                    do i=1,nfp
                        saved_v(nsv+i)=fpar(i)
                    enddo
                endif
! **** Get parameters that were read from file and then stored
                propb    = saved_v(nsv+1)
                tint     = saved_v(nsv+2)
                tdt      = saved_v(nsv+3)
                dband    = saved_v(nsv+4)
                closloop = retolog(saved_v(nsv+5))
                uman     = saved_v(nsv+6)
                umin     = saved_v(nsv+7)
                umax     = saved_v(nsv+8)
            endif
! **** Run controller if mode is closed loop else use "manual" value
            if (closloop) then
! **** Closed loop mode
                if (time>saved_v(1)) then
! **** First call of timestep
! **** Update previous sample instant values
                    do is=4,(nsv-1),2
                        saved_v(is+1)=saved_v(is)
                    enddo
                endif
! **** Update previous values
                intp = saved_v(5)
                difp = saved_v(7)
                errp = saved_v(9)
                up   = saved_v(11)
! **** Apply deadband around set-point
                call deadband(y,w,dband,wdb)
! **** Pid controller
                ux = pidcont(y,wdb,propb,tint,tdt,intp,difp,&
                               up,errp,tsamp,umax,umin)
! **** Save provisional values to be used in updating at next sample instant
                saved_v(4)  = intp
                saved_v(6)  = difp
                saved_v(8)  = errp
                saved_v(10) = ux
! **** Open loop mode
            else
                ux = uman
            endif
! **** Save current sample number and output
            saved_v(2)  = float(nsample)
            saved_v(20) = ux
        else
! **** Not a sample instant, set output to value from previous sample
! **** instant
            ux = saved_v(20)
        endif
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = ux
! **** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo
! **** Return
        return
        end subroutine type441
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
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
! **********************************************************************
!
!   Updated to Fortran 90 April 19, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type481(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=11,ni=4,no=1,np=10,&
                                             nfp=8,ns=nsv+nfp+1
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

 ! **** Declaration of local variables
        real         :: intp
        logical      :: retolog,sfstatus,rfstatus,closloop,auxdis,&
                        logicand
        real         :: ps,psset,propb,tint,tdt,dband,propbhl,pshlimit,&
                        sfspdman,tsamp,difp,errp,pidoutp,pssetdb,pidout,&
                        pidcont,contout,switch,dummy4,dummy3,dummy2,&
                        dummy1,sfspdhl,sfspd
        integer      :: nfile,is,nsample,i
        integer      :: itype=481

! **** Read in inputs
        ps        = xin(1)
        psset     = xin(2)
        sfstatus  = retolog(xin(3))
        rfstatus  = retolog(xin(4))
! **** Read in parameters
        propb     = par_v(1)
        tint      = par_v(2)
        tdt       = par_v(3)
        dband     = par_v(4)
        propbhl   = par_v(5)
        pshlimit  = par_v(6)
        closloop  = retolog(par_v(7))
        sfspdman  = par_v(8)
        tsamp     = par_v(9)
        nfile     = par_v(10)
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
! **** Use parameters from file if controller number not zero
            if (nfile>0) then
! **** First call of time-step - read parameters from file and store
                if (time>saved_v(1)) then
                    call rfile(nfile,'cont',nfp,itype,fpar)
                    do i=1,nfp
                        saved_v(nsv+i)=fpar(i)
                    enddo
                endif
! **** Get parameters that were read form file and then stored
                propb    = saved_v(nsv+1)
                tint     = saved_v(nsv+2)
                tdt      = saved_v(nsv+3)
                dband    = saved_v(nsv+4)
                propbhl  = saved_v(nsv+5)
                pshlimit = saved_v(nsv+6)
                closloop = retolog(saved_v(nsv+7))
                sfspdman = saved_v(nsv+8)
            endif
! **** Run controller if mode is closed loop else use "manual" value
            if (closloop) then
! **** Closed loop mode
                if (time>saved_v(1)) then
! **** First call of timestep
! **** Update previous sample instant values
                    do is=4,(nsv-1),2
                        saved_v(is+1)=saved_v(is)
                    enddo
                endif
! **** Update previous values
                intp      = saved_v(5)
                difp      = saved_v(7)
                errp      = saved_v(9)
                pidoutp   = saved_v(11)
! **** Determine output of main control loop then apply over-ride if
! **** static pressure too high
! **** Apply deadband around set-point
                call deadband(ps,psset,dband,pssetdb)
! **** Pid controller
                pidout    = pidcont(ps,pssetdb,propb,tint,tdt,intp,&
                                    difp,pidoutp,errp,tsamp,1.0,0.0)
! **** Determine whether auxiliary input switch is enabled or disabled
! **** disable auxiliary input if supply fan and return fan ok
                auxdis    = logicand(sfstatus,rfstatus)
! **** Determine controller output signal - select pidout if auxdis is
! **** true else select 0.0
                contout   = switch(auxdis,pidout,0.0)
! **** Determine high limit over-ride value - proportional control
                sfspdhl   = pidcont(ps,pshlimit,propbhl,0.0,0.0,dummy1,&
                                  dummy2,dummy3,dummy4,tsamp,1.0,0.0)
! **** Select minimun controller output signal
                sfspd     = min(contout,sfspdhl)
! **** Save provisional values to be used in updating at next sample instant
                saved_v(4)  = intp
                saved_v(6)  = difp
                saved_v(8)  = errp
                saved_v(10) = pidout
! **** Open loop mode
            else
                sfspd     = sfspdman
            endif
! **** Save current sample number and output
            saved_v(2)  = float(nsample)
            saved_v(20) = sfspd
        else
! **** Not a sample instant, set output to value from previous sample
! **** instant
            sfspd = saved_v(20)
        endif
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = sfspd
! **** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo
! **** Return
        return
        end subroutine type481

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
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
! **********************************************************************
!
!   Updated to Fortran 90 April 19, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type482(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=11,ni=5,no=1,np=8,&
                                             nfp=6,ns=nsv+nfp+1
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

 ! **** Declaration of local variables
        real         :: intp
        logical      :: retolog,sfstatus,rfstatus,closloop,auxdis,&
                        logicand
        real         :: supflow,retflow,dflowset,propb,tint,&
                        tdt,deadb,rfspdman,tsamp,difp,errp,pidoutp,&
                        dflow,setdb,pidout,pidcont,rfspd,switch
        integer      :: nfile,is,nsample,i
        integer      :: itype=482

! **** Read in inputs
        supflow  = xin(1)
        retflow  = xin(2)
        sfstatus = retolog(xin(3))
        rfstatus = retolog(xin(4))
        dflowset = xin(5)
! **** Read in parameters
        propb    = par_v(1)
        tint     = par_v(2)
        tdt      = par_v(3)
        deadb    = par_v(4)
        closloop = retolog(par_v(5))
        rfspdman = par_v(6)
        tsamp    = par_v(7)
        nfile    = par_v(8)
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
! **** Use parameters from file if controller number not zero
            if (nfile>0) then
! **** First call of time-step - read parameters from file and store
               if (time>saved_v(1)) then
                    call rfile(nfile,'cont',nfp,itype,fpar)
                    do i=1,nfp
                        saved_v(nsv+i)=fpar(i)
                    enddo
               endif
! **** Get parameters that were read form file and then stored
               propb    = saved_v(nsv+1)
               tint     = saved_v(nsv+2)
               tdt      = saved_v(nsv+3)
               deadb    = saved_v(nsv+4)
               closloop = retolog(saved_v(nsv+5))
               rfspdman = saved_v(nsv+6)
            endif
! **** Run controller if mode is closed loop else use "manual" value
            if (closloop) then
! **** Closed loop mode
                if (time>saved_v(1)) then
! **** First call of timestep
! **** Update previous sample instant values
                    do is=4,(nsv-1),2
                        saved_v(is+1) = saved_v(is)
                    enddo
                endif
! **** Update previous values
                intp    = saved_v(5)
                difp    = saved_v(7)
                errp    = saved_v(9)
                pidoutp = saved_v(11)
! **** Difference of supply and return air volume
                dflow   = supflow - retflow
! **** Apply deadband around set-point
                call deadband(dflow,dflowset,deadb,setdb)
! **** Pid controller
                pidout  = pidcont(dflow,setdb,propb,tint,tdt,intp,&
                                 difp,pidoutp,errp,tsamp,1.0,0.0)
! **** Determine whether auxiliary input switch is enabled or disabled
! **** disable auxiliary input if supply fan and return fan ok
                auxdis = logicand(sfstatus,rfstatus)
! **** Determine controller output signal - - select pidout if auxdis is
! **** true else select 0.0
                rfspd = switch(auxdis,pidout,0.0)
! **** Save provisional values to be used in updating at next sample instant
                saved_v(4)  = intp
                saved_v(6)  = difp
                saved_v(8)  = errp
                saved_v(10) = pidout
! **** Open loop mode
            else
                rfspd     = rfspdman
            endif
! **** Save current sample number and output value
            saved_v(2)  = float(nsample)
            saved_v(18) = rfspd
        else
! **** Not a sample instant, set output(s) to value(s) from previous sample
! **** instant
            rfspd = saved_v(18)
        endif
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = rfspd
! **** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo
! **** Return
        return
        end subroutine type482

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     E-51 return fan reset control
! *
! * PURPOSE:        Determine set-point for return fan controller
! *
! ***********************************************************************
! * INPUTS
! * ======
! *  1. mindpr  : minimun outside air damper position (open=TRUE)     (-)
! *
! * OUTPUTS
! * =======
! *  1. dflowset: return fan controller set-point                  (m3/s)
! *
! * PARAMETERS
! * ==========
! *  1. dflowsop: return fan set-point when min OA damper open     (m3/s)
! *  2. tsamp   : sample time                                         (s)
! *  3. nfile   : controller number (parameters in file contN.par)    (-)
! *
! * SAVED
! * ======
! *  1. time    : time of previous call
! *  2. nsample : sample number of previous controller execution
! *  3. nsample : sample number of previous controller sample
! *  4.         : PAR(1) read from file
! *  5. dflowset: return fan controller set-point from previous call
! *
! **********************************************************************
!
!   Updated to Fortran 90 April 19, 2007 Cheol Park, NIST
!
! ************************************************************************

        subroutine type483(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=3,ni=1,no=1,np=3,&
                                             nfp=1,ns=nsv+nfp+1
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar
        integer                           :: itype=483

! **** Declaration of local variables
        logical      :: mindpr,retolog
        real         :: dflowsop,tsamp,dflowset,switch
        integer      :: nfile,nsample,i

! **** Read in inputs
        mindpr = retolog(xin(1))
! **** Read in parameters
        dflowsop = par_v(1)
        tsamp    = par_v(2)
        nfile    = par_v(3)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
                saved_v(2) = 0.
            endif
        endif
! **** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! **** Run controller if a sample instant
        nsample=nint(saved_v(3))+1
        if (time>=(nsample*tsamp) .or. time==0.0) then
! **** Use parameters from file if controller number not zero
            if (nfile>0) then
                if (time>saved_v(1)) then
! **** First call of timestep - read parameters from file and store
                    call rfile(nfile,'cont',nfp,itype,fpar)
                    saved_v(nsv+1)=fpar(1)
                endif
                dflowsop  = saved_v(nsv+1)
             endif
! **** Run controller
! **** Set-point for return fan controller is zero if minimum outside air
! **** Damper is closed
            dflowset = switch(mindpr,dflowsop,0.0)
! **** Save current sample number and output value
            saved_v(2) = float(nsample)
            saved_v(5) = dflowset
        else
! **** Not a sample instant, set output(s) to value(s) from previous sample
! **** instant
            dflowset = saved_v(5)
        endif
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = dflowset
! **** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo
! **** Return
	    return
	    end subroutine type483
! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     E51 minimum outside air damper control
! *
! * PURPOSE:        Determine the demanded position of the minimun
! *                 outside air damper (fully open or fully closed)
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. tsup    : supply air temperature sensor                      (C)
! *  2. tout    : outside air temperature sensor                     (C)
! *  3. sfstatus: supply fan status (1 = on, 0 = off)                (-)
! *
! * OUTPUTS
! * =======
! *  1. damprpos: minimun outside air damper position (0 or 1)       (-)
! *
! * PARAMETERS
! * ==========
! *  1. tsuplim : supply air temperature limit                       (C)
! *  2. toutlim : outside air temperature limit                      (C)
! *  3. tsamp   : sampling interval                                  (s)
! *  4. nfile   : controller number (parameters in file contN.par)   (-)
! *
! * SAVED
! * =====
! * 1.  time    : time of previous call of TYPE
! * 2.  nsample : sample number of previous controller execution
! * 3.  nsample : sample number of previous controller sample
! * 4-5.        : PAR(1)-PAR(2) read from file
! * 6.  damprpos: (controller output from previous sample)
! *
!**************************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Li Mei and Philip Haves
!                        Loughborough University of Technology
!
!   TSUPE:               November 22, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   RETOLOG,LOGICNOT,LOGICAND,LOGICOR,COMPARE
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
!
!   Updated to Fortran 90 April 19, 2007 Cheol Park, NIST
!
! *************************************************************************

        subroutine type484(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=3,ni=3,no=1,np=4,&
                                             nfp=2,ns=nsv+nfp+1
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

! **** Declaration of local variables
        logical      :: retolog,sfstatus,l1,l2,l3,l4,l5,compare
        logical      :: logicand,logicor,logicnot
        real         :: tsup,tout,tsuplim,toutlim,tsamp,damprpos
        integer      :: nfile,nsample,i
        integer      :: itype=484

! **** Read in inputs
        tsup      = xin(1)
        tout      = xin(2)
        sfstatus = retolog(xin(3))
! **** Read in parameters
        tsuplim    = par_v(1)
        toutlim    = par_v(2)
        tsamp      = par_v(3)
        nfile      = par_v(4)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
		saved_v(1) = -999999.
		saved_v(2) = 0.
	    endif
        endif
! **** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! **** Run controller if a sample instant
        nsample=nint(saved_v(3))+1
        if (time>=(nsample*tsamp) .or. time==0.0) then
! **** Use parameters from file if controller number not zero
            if (nfile>0) then
                if (time>saved_v(1)) then
! **** First call of timestep - read parameters from file and store
                    call rfile(nfile,'cont',nfp,itype,fpar)
                    saved_v(nsv+1) = fpar(1)
                    saved_v(nsv+2) = fpar(2)
                endif
                tsuplim  = saved_v(nsv+1)
                toutlim  = saved_v(nsv+2)
             endif
! **** Run controller
! **** Compare outside air temperature with outside air temperature
! **** limit - true if toutlim > tout
            l1 = compare(toutlim,tout)
! **** Compare discharge air temperature with discharge air temperature
! **** limit - true if tsuplim > tsup
            l2 = compare(tsuplim,tsup)
! **** l1 nor l2
            l3 = logicor(l1,l2)
            l4 = logicnot(l3)
! **** Determine minimum outside air damper control signal
            l5 = logicand(sfstatus,l4)
            if (l5) then
! **** Damper is open if supply fan ok and tout > toutlim and tsup > tsuplim
                damprpos = 1.0
            else
                damprpos = 0.0
            endif
! **** Save current sample number and output
            saved_v(2) = float(nsample)
            saved_v(6) = damprpos
        else
! **** Not a sample instant, set output(s) to value(s) from previous sample
! **** instant
            damprpos = saved_v(6)
        endif
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = damprpos
! **** Disallow freezing (parameters read from file may change)
        do i=1,no
            iostat(i) = 0
        enddo
! **** Return
        return
        end subroutine type484
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
!   Updated to Fortran 90 April 19, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type485(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=11,ni=5,no=3,np=11,&
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
                        ddemr,ddeme
        integer      :: nfile,is,nsample,i
        integer      :: itype=485

! **** Read in inputs
        tmix        = xin(1)
        mbcdem      = xin(2)
        sfstatus    = retolog(xin(3))
        lowtovr     = retolog(xin(4))
        econ        = retolog(xin(5))
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
                ddemo = min(pidout,passthru)
! **** Return and exhaust air dampers demanded position
                ddemr = 1.0 - ddemo
                ddeme = ddemo
! **** Save provisional values to be used in updating at next sample instant
                saved_v(4) = intp
                saved_v(6) = difp
                saved_v(8) = errp
                saved_v(10) = pidout
! **** Open loop mode
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
        end subroutine type485

!***********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
!***********************************************************************
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
! *
! * PARAMETERS
! * ==========
! *  1. propb   : proportional band                                  (K)
! *  2. tint    : integral time                                      (s)
! *  3. tdt     : derivative time                                    (s)
! *  4. cbreak  : breakpoint between damper and cooling coil demand(0-2)
! *  5. deadb   : deadband                                           (K)
! *  6. closloop: control mode (0=open loop, 1=closed loop)          (-)
! *  7. cdemman : open loop cooling coil demand (0-1)                (-)
! *  8. tsamp   : sampling interval                                  (s)
! *  9. nfile   : controller number (parameters in file contN.par)   (-)
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
!   Updated to Fortran 90 April 19, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type486(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=11,ni=4,no=2,np=9,&
                                             nfp=7,ns=nsv+nfp+2
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

! **** Declaration of local variables
        real         :: intp
        logical      :: retolog,logicnot,logicand
        logical      :: sfstatus,lowtovr,nlowtovr,auxdis,closloop
        real         :: tsup,tset,propb,tint,tdt,cbreak,deadb,cdemman,&
                        tsamp,difp,errp,pidoutp,tsetdb,pidout,pidcont,&
                        auxinp,switch,contout,ddem,span,cdem
        integer      :: nfile,is,nsample,i
        integer      :: itype=486

! **** Read in inputs
        tsup       = xin(1)
        sfstatus   = retolog(xin(2))
        lowtovr    = retolog(xin(3))
        tset       = xin(4)
! **** Read in parameters
        propb      = par_v(1)
        tint       = par_v(2)
        tdt        = par_v(3)
        cbreak     = par_v(4)
        deadb      = par_v(5)
        closloop   = retolog(par_v(6))
        cdemman    = par_v(7)
        tsamp      = par_v(8)
        nfile      = nint(par_v(9))
! **** Initialize at beginning of simulation
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
! **** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! **** Run controller if a sample instant
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
                cbreak     = saved_v(nsv+4)
                deadb      = saved_v(nsv+5)
                closloop   = retolog(saved_v(nsv+6))
                cdemman    = saved_v(nsv+7)
            endif
! **** Run controller if mode is closed loop else use "manual" value
            if (closloop) then
! **** Closed loop mode
                if (time>saved_v(1)) then
! **** First call of timestep
! **** Update previous sample instant values
                    do is=4,nsv-1,2
                        saved_v(is+1)=saved_v(is)
                    enddo
                endif
! **** Update previous values
                intp = saved_v(5)
                difp = saved_v(7)
                errp = saved_v(9)
                pidoutp = saved_v(11)
! **** Apply deadband around set-point
                call deadband(tsup,tset,deadb,tsetdb)
! **** Pid controller
                pidout = pidcont(tsup,tsetdb,propb,tint,tdt,intp,difp,&
                                 pidoutp,errp,tsamp,2.0,0.0)
! **** Determine whether auxiliary input switch is enabled or disabled
! **** compliment of low temperature over-ride
                nlowtovr = logicnot(lowtovr)
! **** Disable auxiliary input if fan ok and low temperature over-ride
! **** not set
                auxdis = logicand(sfstatus,nlowtovr)
! **** Determine auxiliary input value - cooling coil valve fully open
! **** if low temperature over-ride set else fully closed if fan status
! **** not ok and low temperature over-ride not set
                auxinp = switch(lowtovr,2.0,0.0)
! **** Determine controller output signal - select pidout if auxdis is
! **** true else select auxinp
                contout = switch(auxdis,pidout,auxinp)
! **** Sequence demands
! **** mixing damper demand
                if (cbreak/=0.0) then
                    ddem = span(contout,0.0,cbreak,0.0,1.0)
                else
                    stop 'type 486: cbreak = 0'
                endif
! **** Cooling coil demand
                if (cbreak/=2.0) then
                    cdem = span(contout,cbreak,2.0,0.0,1.0)
                else
                    stop 'type 486: cbreak = 2'
                endif
! **** Save provisional values to be used in updating at next sample instant
                saved_v(4) = intp
                saved_v(6) = difp
                saved_v(8) = errp
                saved_v(10) = pidout
! **** Open loop mode (dampers set manually in type265)
            else
                ddem = 0.0
                cdem = cdemman
            endif
! **** Save current sample number
            saved_v(2) = float(nsample)
! **** Save outputs for use when not a sample instant
            saved_v(19) = ddem
            saved_v(20) = cdem
        else
! **** Not a sample instant, set output to value from prev sample instant
            ddem = saved_v(19)
            cdem = saved_v(20)
        endif
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = ddem
        yout(2) = cdem
! **** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo
! **** Return
        return
        end subroutine type486

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     E51 economizer control
! *
! * PURPOSE:        Determine the economizer control mode by comparing
! *                 outside and return air enthalpy
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. tout    : outdoor air temperature sensor                     (C)
! *  2. rhout   : outdoor air relative humidity sensor (0-1)         (-)
! *  3. tret    : return air temperature sensor                      (C)
! *  4. rhret   : return air relative humidity sensor (0-1)          (-)
! *
! * OUTPUTS
! * =======
! *  1. econ    : economizer status (1 = OA enthalpy > RA enthalpy)  (-)
! *
! * PARAMETERS
! * ==========
! *  1. tsamp   : sampling interval                                  (s)
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
!   Updated to Fortran 90 April 19, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type487(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=3,ni=4,no=1,np=1,&
                                             nfp=0,ns=nsv+nfp+1
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat

! **** Declaration of local variables
        real         :: logtore
        logical      :: econ,retolog,compare
        real         :: tout,rhout,tret,rhret,tsamp,oaenth,enthalpy,raenth
        integer      :: nsample,i

! **** Read in inputs
        tout = xin(1)
        rhout = xin(2)    ! changed 12/6/1999
        tret = xin(3)
        rhret = xin(4)    ! changed 12/6/1999
! **** Read in parameter
        tsamp = par_v(1)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
		saved_v(1) = -999999.
		saved_v(2) = 0.
	    endif
        endif
! **** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! **** Run controller if a sample instant
        nsample=nint(saved_v(3))+1
        if (time>=(nsample*tsamp) .or. time==0.0) then
! **** Run controller
! **** Calculate enthalpy of outside air
            oaenth = enthalpy(tout,rhout)
! **** Calculate enthalpy of return air
            raenth = enthalpy(tret,rhret)
! **** Determine economizer control signal
            econ = compare(oaenth,raenth)
! **** Save current sample number and output
	    saved_v(2) = float(nsample)
	    saved_v(4) = logtore(econ)
        else
! **** Not a sample instant, set output(s) to value(s) from previous sample
! **** instant
            econ = retolog(saved_v(4))
        endif
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = logtore(econ)
! **** Allow freezing of algebraic variables
        do i=1,no
            iostat(i) = 1
        enddo
! **** Return
	    return
	    end subroutine type487

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     E51 low temperature override control
! *
! * PURPOSE:        Calculate low temperature override signal for
! *                 cooling coil and maximum outside air damper control
! *
! **********************************************************************
! * INPUTS
! * ======
! *  1. tsup    : discharge air temperature                          (C)
! *  2. tout    : outdoor air temperature                            (C)
! *  3. sfstatus: supply fan status (1 = on, 0 = off)                (-)
! *
! * OUTPUTS
! * =======
! *  1. lowtovr : low temperature override signal                    (-)
! *
! * PARAMETERS
! * ==========
! *  1. tsuplim : supply air temperature limit                       (C)
! *  2. toutlim : outdoor air temperature limit                      (C)
! *  3. tsamp   : sample time                                        (s)
! *  4. nfile   : controller number (parameters in file contN.par)   (-)
! *
! * SAVED
! * =====
! *  1. time    : time of previous call of TYPE
! *  2. nsample : sample number of previous controller execution
! *  3. nsample : sample number of previous controller sample
! *  4-5.       : PAR(1)-PAR(2) read from file
! *  6. lowtovr : controller output from previous sample
! *
!***********************************************************************
!
!   MAJOR RESTRICTIONS:
!
!   DEVELOPER:           Li Mei and Philip Haves
!                        Loughborough University of Technology
!
!   TSUPE:               November 22, 1995
!
!   INCLUDE FILES:
!   SUBROUTINES CALLED:
!   FUNCTIONS  CALLED:   RETOLOG,LOGICNOT,LOGICAND,LOGICOR,COMPARE,LOGTORE
!
!   REVISION HISTORY:    None
!
!   REFERENCE:           ASHRAE 825-RP Final Report
!
! **********************************************************************
!
!   Updated to Fortran 90 April 19, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type488(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=3,ni=3,no=1,np=4,&
                                             nfp=2,ns=nsv+nfp+2
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar
        integer                           :: itype=488

! **** Declaration of local variables
        real         :: logtore
        logical      :: retolog,sfstatus,nsfstat,l1,l2,l3,l4,lowtovr,&
                        compare
        logical      :: logicnot,logicand,logicor
        real         :: tsup,tout,tsuplim,toutlim,tsamp
        integer      :: nfile,nsample,i

! **** Read in inputs
        tsup      = xin(1)
        tout      = xin(2)
        sfstatus  = retolog(xin(3))
! **** Read in parameters
        tsuplim   = par_v(1)
        toutlim   = par_v(2)
        tsamp     = par_v(3)
        nfile     = par_v(4)
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
                saved_v(2) = 0.
            endif
        endif
! **** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! **** Run controller if a sample instant
        nsample=nint(saved_v(3))+1
        if (time>=(nsample*tsamp) .or. time==0.0) then
! **** Use parameters from file if controller number not zero
            if (nfile>0) then
                if (time>saved_v(1)) then
! **** First call of timestep - read parameters from file and store
                    call rfile(nfile,'cont',nfp,itype,fpar)
                    saved_v(nsv+1)=fpar(1)
                    saved_v(nsv+2)=fpar(2)
                endif
                tsuplim   = saved_v(nsv+1)
                toutlim   = saved_v(nsv+2)
             endif
! **** Run controller
! **** Not supply fan status
            nsfstat = logicnot(sfstatus)
! **** Compare outside air temperature with outside air temperature
! **** limit - true if tout < toutlim
            l1 = compare(toutlim,tout)
! **** True if supply fan is not on and outside air temperature is lower
! **** than limit
            l2 = logicand(nsfstat,l1)
! **** Compare discharge air temperature with discharge air temperature
! **** limit - true if tsup < tsuplim
            l3 = compare(tsuplim,tsup)
! **** If supply fan is on and discharge air temperature is lower
! **** than limit
            l4 = logicand(sfstatus,l3)
! **** Determine low temperature over-ride - set if supply fan not ok and
! **** tout < toutlim, or if supply fan ok and tsup < tsuplim
            lowtovr = logicor(l2,l4)
! **** Save current sample number and output
            saved_v(2) = float(nsample)
            saved_v(6) = logtore(lowtovr)
        else
! **** Not a sample instant, set output(s) to value(s) from previous sample
! **** instant
            lowtovr = retolog(saved_v(6))
        endif
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = logtore(lowtovr)
! **** Disallow freezing (parameters read from file may change)
        do i=1,no
            iostat(i) = 0
        enddo
! **** Return
        return
        end subroutine type488

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
! * SUBROUTINE:     E51 supply air temperature reset
! *
! * PURPOSE:        Reset supply air temperature setpoint so as to keep
! *                 the maximum of several room temperatures at a
! *                 specified value.
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
!   Updated to Fortran 90 April 19, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type489(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=7,ni=6,no=1,np=10,&
                                             nfp=7,ns=nsv+nfp
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

! **** Declaration of local variables
        real,dimension(ni)                :: tz
        real         :: intp
        real         :: tzset,propb,tint,tssetmax,tssetmin,tsset0,&
                        deadb,tssetp,tzmax,biggest,tzsetdb,dtsmax,&
                        dtsmin,dummy,tsset,pidcont,tsamp
        integer      :: ninputs,nfile,i,nsample,is
        integer      :: itype=489

! **** Read in parameters
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
! **** Read in inputs
        do i=1,ninputs
            tz(i) = xin(i)
        enddo
! **** Initialize at beginning of simulation
        if (itime<=1) then
            if (init==0 .or. saved_v(1)>time) then
                saved_v(1) = -999999.
                saved_v(2) = 0.
            endif
            if (init==0) then
! **** Initialize integral term and controller output
                saved_v(4) = 0.0
                saved_v(6) = tsset0
            endif
        endif
! **** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! **** Run controller if a sample instant
        nsample = nint(saved_v(3))+1
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
! **** Get parameters that were read form file and then stored
                tzset      = saved_v(nsv+1)
                propb      = saved_v(nsv+2)
                tint       = saved_v(nsv+3)
                tssetmax   = saved_v(nsv+4)
                tssetmin   = saved_v(nsv+5)
                tsset0     = saved_v(nsv+6)
                deadb      = saved_v(nsv+7)
            endif
            if (time>saved_v(1)) then
! **** First call of timestep
! **** Update previous sample instant values
                do is=4,(nsv-1),2
                    saved_v(is+1)=saved_v(is)
                enddo
            endif
! **** Update previous values
            intp   = saved_v(5)
            tssetp = saved_v(7)
! **** Select the highest room temperature
            tzmax = biggest(tz,ninputs)
! **** Apply dead-band
            call deadband(tzmax,tzset,deadb,tzsetdb)
! **** Pi control
! **** Setpoint can vary between tsset0+dtsmax and tsset0+dtsmin
            dtsmax = tssetmax-tsset0
            dtsmin = tssetmin-tsset0
! **** Pid output is offset of setpoint from tsset0
            tsset  = tsset0+pidcont(tzmax,tzsetdb,propb,tint,0.0,intp,&
                               dummy,tssetp,dummy,tsamp,dtsmax,dtsmin)
! **** Save current sample number
            saved_v(2) = float(nsample)
! **** Save provisional values to be used in updating at next sample instant
            saved_v(4) = intp
            saved_v(6) = tsset
        else
! **** Not a sample instant, set output(s) to value(s) from previous sample
! **** instant
            tsset = saved_v(6)
        endif
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = tsset
! **** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo
! **** Return
        return
        end subroutine type489

! **********************************************************************
! * Copyright ASHRAE.   Control Simulation Testbed
! **********************************************************************
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
! *
! * OUTPUTS
! * =======
! *  1. vdem    : normalised volumetric flow demand (0-1)            (-)
! *  2. rdem    : reheat coil demand (0-1)                           (-)
! *  3. contout : room demand (heating and cooling) (-1 - +1)        (-)
! *
! * PARAMETERS
! * ==========
! *  1. tsetheat: heating setpoint for zone                          (C)
! *  2. tsetcool: cooling setpoint for zone                          (C)
! *  3. propb   : proportional band                                  (K)
! *  4. tint    : integral time                                      (s)
! *  5. tdt     : derivative time                                    (s)
! *  6. rbreak  : breakpoint btwn damper and reheat coil signals (-1-+1)
! *  7. tdmin   : minimum turndown ratio                             (-)
! *  8. closloop: control mode (0=open loop, 1=closed loop)          (-)
! *  9. vdemman : open loop demanded normalized volumetric flow rate (-)
! * 10. rdemman : open loop reheat coil demand (0-1)                 (-)
! * 11. tsamp   : sampling interval                                  (s)
! * 12. nfile   : controller number (parameters in file contN.par)   (-)
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
! **********************************************************************
!
!   Updated to Fortran 90 April 19, 2007 Cheol Park, NIST
!
! **********************************************************************

        subroutine type490(xin,yout,par_v,saved_v,iostat)

        use modsim_head
        implicit none
        integer,parameter                 :: nsv=13,ni=3,no=3,np=12,&
                                             nfp=10,ns=nsv+nfp+3
        real,dimension(ni),intent(in)     :: xin                      
        real,dimension(no),intent(out)    :: yout                     
        real,dimension(np),intent(in)     :: par_v                    
        real,dimension(ns),intent(in out) :: saved_v                  
        integer,dimension(no),intent(out) :: iostat
        real,dimension(nfp)               :: fpar

! **** Declaration of local variables
        real         :: logtore, intp
        logical      :: retolog,sfstatus,suphotp,suphot,comphys,&
                        closloop
        real         :: tzon,tsup,tsetheat,tsetcool,propb,tint,tdt,&
                        rbreak,tdmin,vdemman,rdemman,tsamp,contout,&
                        difp,errp,pidoutp,tset,deadb,tsetdb,pidout,&
                        pidcont,switch,vdem,span,rdem
        integer      :: nfile,is,nsample,i
        integer      :: itype=490

! **** Read in inputs
        tzon     = xin(1)           
        tsup     = xin(2)           
        sfstatus = retolog(xin(3))
! **** Read in parameters
        tsetheat = par_v(1)
        tsetcool = par_v(2)
        propb    = par_v(3)
        tint     = par_v(4)
        tdt      = par_v(5)
        rbreak   = par_v(6)
        tdmin    = par_v(7)
        closloop = retolog(par_v(8))
        vdemman  = par_v(9)
        rdemman  = par_v(10)
        tsamp    = par_v(11)
        nfile    = nint(par_v(12))
! **** Initialize at beginning of simulation
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
! **** Update number of previous sample if new step-time
        if (time>saved_v(1)) then
            saved_v(3) = saved_v(2)
        endif
! **** Run controller if a sample instant
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
                tsetheat   = saved_v(nsv+1)
                tsetcool   = saved_v(nsv+2)
                propb      = saved_v(nsv+3)
                tint       = saved_v(nsv+4)
                tdt        = saved_v(nsv+5)
                rbreak     = saved_v(nsv+6)
                tdmin      = saved_v(nsv+7)
                closloop   = retolog(saved_v(nsv+8))
                vdemman    = saved_v(nsv+9)
                rdemman    = saved_v(nsv+10)
            endif
! **** Run controller if mode is closed loop else use "manual" value
            if (closloop) then
! **** Closed loop mode
                if (time>saved_v(1)) then
! **** First call of timestep
! **** Update previous sample instant values
                    do is=4,nsv-1,2
                        saved_v(is+1)=saved_v(is)
                    enddo
                endif
! **** Update previous values
                intp    = saved_v(5)
                difp    = saved_v(7)
                errp    = saved_v(9)
                pidoutp = saved_v(11)
                suphotp = retolog(saved_v(13))
! **** Apply deadband around set-point
                tset  = (tsetcool+tsetheat)/2.
                deadb = (tsetcool-tsetheat)/2.
                call deadband(tzon,tset,deadb,tsetdb)
! **** Pid controller - output in the range -1 - +1
                pidout = pidcont(tzon,tsetdb,propb,tint,tdt,intp,difp,&
                                 pidoutp,errp,tsamp,1.0,-1.0)
! **** Set output to zero if fan status not ok
                contout = switch(sfstatus,pidout,0.0)
! **** Sequence demands
! **** Vav damper demand - depends on whether supply air is hotter than
! **** Zone air - comparison uses 1 k hysteresis
                suphot = comphys(tsup,tzon,1.0,suphotp)
! **** Minimum turn-down ratio
                if (suphot) then
! **** Increased heating demand produces higher demanded flow rate
                    vdem = span(contout,-1.0,rbreak,tdmin,1.0) ! revised 12/6/99
                else
! **** Increased heating demand produces lower demanded flow rate
                    vdem = span(contout,rbreak,-1.0,tdmin,1.0) ! revised 12/6/99
                endif
! **** Reheat coil demand
                rdem = span(contout,rbreak,1.0,0.0,1.0)        ! revised 12/6/99
! **** Save provisional values to be used in updating at next sample instant
                saved_v(4)  = intp
                saved_v(6)  = difp
                saved_v(8)  = errp
                saved_v(10) = pidout
                saved_v(12) = logtore(suphot)
! **** Open loop mode - set manually
            else
                vdem = vdemman
                rdem = rdemman
            endif
! **** Save current sample number
            saved_v(2) = float(nsample)
! **** Save outputs for use when not a sample instant
            saved_v(24) = vdem
            saved_v(25) = rdem
            saved_v(26) = contout
        else
! **** Not a sample instant, set output to value from prev sample instant
            vdem    = saved_v(24)
            rdem    = saved_v(25)
            contout = saved_v(26)
        endif
! **** Save time of current call
        saved_v(1) = time
! **** Output
        yout(1) = vdem
        yout(2) = rdem
        yout(3) = contout
! **** Disallow freezing
        do i=1,no
            iostat(i) = 0
        enddo   
! **** Return
        return
        end  subroutine type490
                
