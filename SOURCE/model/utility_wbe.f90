! ********************************************************************

	real function pid(p,o,u,l,ms,ss,g,ti,td,integ,dif,ppid,&
                      pp,rsec,nseq)                                   ! 11/19/02

! ********************************************************************
!
!   PID controller module - engineering units - formatted for use with
!   TREND macros
!
!   P Haves, University of Oxford, August 1989
!
!   Uses feedback form of integral and derivative action
!   ("PID Algorithms and their Computer Implementation",
!   D.W. Clarke, Oxford University Engineering Laboratory
!   Report No. 1482/83)
!
!       pid  =  controller output (0-100)
!    	p    =  process variable
!  	o    =  occupation setpoint
!   	u    =  unoccupied setpoint
!	l    =  manual level
!	ms   =  manual select (0=manual override disabled, 1=enabled)
!	ss   =  setpoint select (0=unoccupied, 1=occupied)
!	g    =  proportional gain (% o/p for unit error)
!	ti   =  integral time (sec)
!	td   =  derivative time (sec)
!   	integ=  integral term from previous sample instant (entry),
!		current value of integral term (return)
!   	dif  =  derivative term from previous sample instant (entry),
!		current value of derivative term (return)
!       ppid =  controller output from previous sample instant
!	pp   =  process variable from previous sample instant
! 	rsec =  reschedule time (sec) - sequence time
!	nseq =  number of times in sequence table
!
! ********************************************************************

        implicit none
        real                       :: l,iterm,integ                                       ! 11/19/02
        integer                    :: ss
        logical                    :: zeroi,zerod
        real                       :: tsamp,rsec,p,o,u,g,ti,td,dif,ppid,pp,&
                                      e,pterm,beta,dterm,d1,d2
        integer                    :: ms,nseq

!   Check for manual override
        if (ms==1) then
	    pid=l
	    pid=min(100.0,max(0.0,pid))
	    return
	elseif (ms/=0) then
	    stop 'ms error in function pid'
	endif
!  Sample time
	tsamp=rsec/nseq
!  Error
	if (ss==1) then
	    e=o-p
	elseif (ss==0) then
	    e=u-p
	else
	    stop 'ss error in function pid'
	endif
!
!  Controller
!  ----------
	zeroi=(ti<=0.0).or.((tsamp/ti)<1.0e-3)
	if (zeroi) then
	    pterm=g*e
	    iterm=0.0
	else
	    pterm=g*e*(1.0+0.5*tsamp/ti)
	    beta=(2.0*ti-tsamp)/(2.0*ti+tsamp)
	    if (beta<0.0) write(*,*) 'func pid: ti < tsamp/2'
 	    iterm=beta*integ+(1.0-beta)*ppid
	endif
!  Derivative term  - incorrectly implemented !!!!!
	zerod=(td/tsamp)<1.0e-3
	if (zerod) then
	    dterm=0.0
	else
	    d1=(0.2*td-tsamp)/(0.2*td+tsamp)
	    if (d1<0.0) write(*,*) 'func pid: td < 5*tsamp'
	    d2=2.0*td/(0.2*td+tsamp)
	    dterm=d1*dif + d2*(p-pp)
	endif
!  Controller output
	pid=pterm+iterm
	pid=min(100.0,max(0.0,pid))
!  Update previous integral and derivative terms
	integ=iterm
	dif=dterm

	return
	end function pid

! ********************************************************************

	real function rsclfr(e,f,g)

! ********************************************************************

        implicit none
        real                       :: e,f,g                                       ! 11/19/02

!  Rescale from 0-100 to e-f
	rsclfr=e+g*(f-e)/100.0
	if (e<f) then
	    rsclfr=min(f,max(e,rsclfr))
	else
	    rsclfr=min(e,max(f,rsclfr))
	endif

	return
	end function rsclfr

! ********************************************************************

	real function rsclto(e,f,g)

! ********************************************************************

        implicit none
        real                       :: e,f,g                                       ! 11/19/02

!  Rescale to 0-100 from e-f
	rsclto=100.0*(g-e)/(f-e)
	rsclto=min(100.0,max(0.0,rsclto))

	return
	end function rsclto

! ********************************************************************

	real function xlimat(e,f,g)

! ********************************************************************

        implicit none
        real                       :: e,f,g                                       ! 11/19/02

!  Limit g to range e-f
	if (g>e) then
	    xlimat=e
	elseif (g<f) then
	    xlimat=f
	else
	    xlimat=g
	endif

	return
	end function xlimat

! ********************************************************************

	real function add(e,f,g,h)

! ********************************************************************
        implicit none
        real                       :: e,f,g,h                                       ! 11/19/02

	add=e*g+f*h

	return
	end function add

! ********************************************************************

	integer function icomp(e,f)

! ********************************************************************
        implicit none
        real                       :: e,f                                       ! 11/19/02

	if (f>e) then
	    icomp=1
	else
	    icomp=0
	endif

	return
	end function icomp

! ********************************************************************

	integer function ihyst(g,e,f,d)

! ********************************************************************
        implicit none
        real                       :: e,f,g                                       ! 11/19/02
	integer                    :: d

	if (g<(e-f/2.0)) then
	    ihyst=0
	elseif (g>(e+f/2.0)) then
	    ihyst=1
	else
	    ihyst=d
	endif
	d=ihyst

	return
	end function ihyst

! ********************************************************************

	integer function itimer(inx,ondelay,minon,offdelay,inpx,&
                        tinon,tinoff,itimerp,touton)

! ********************************************************************

        use modsim_head
        implicit none
        real                       :: ondelay,offdelay,tinon,tinoff,touton                                        ! 11/19/02
        integer                    :: inx,inpx,itimerp
!  Logic timer module
	    real                       :: minon

!  Check for change in input
	if (inx==1 .and. inpx==0) tinon=time
	if (inx==0 .and. inpx==1) tinoff=time
!  Output on if input on for longer than on-delay
	if ( (inx==1 .and. (time-tinon)>=ondelay)&
!  or on previously and minimum on-time not exceeded
        .or. (itimerp==1 .and. (time-touton)<minon)&
!  or input off time less than off-delay
        .or. ((time-tinoff)<offdelay) )  then
	    itimer=1
	else
	    itimer=0
	endif
!  Update previous states
	if (itimer==1 .and. itimerp==0) touton=time
	inpx=inx
	itimerp=itimer
	return
	end function itimer

