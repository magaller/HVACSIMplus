! **********************************************************************

       subroutine rdenv(isshel)

! ----------------------------------------------------------------------
!
!     WEATHER AND CONDUCTION TRANSFER FUNCTIONS:
!        Read the weather data file and the file of conduction transfer
!        functions of building constructs.
!        Interpolate hourly weather data.
!        Calculate sines of solar altitude and azimuth angles.
!
!     December 6, 1984  Cheol Park
!     Last Updated: August 1, 1985  D. Clark
!     Updated:  October 18, 2002  Cheol Park & Mike Galler, NIST
!
!     Updated:  May 2, 2007  C. Park
!               Converted to Fortran 90.
!
!     INPUTS FROM FILE  ifile7 (ctfdata.dat):
!       idstr:    ID number of structures, [1,10] (-)
!       nctfa:    Number of CTF terms, [1,20] (-)
!       norda:    Number of order, [1,5] (-)
!       u:        Overall conductance (W/m**2-K)
!       x:        Conduction transfer function (W/m**2-K)
!       y:        Conduction transfer function (W/m**2-K)
!       z:        Conduction transfer function (W/m**2-K)
!       rj:       Flux referring CTF           (W/m**2)
!
!     INPUT FROM FILE ifile8 (weather.dat):
!       month:    Month of year (-)
!       alat:     Latitude (degree)
!       along:    Longitude (degree)
!       tzn:      Time zone number (-)
!                 =4 --- Atlantic
!                 =5 --- Eastern
!                 =6 --- Central
!                 =7 --- Mountain
!                 =8 --- Pacific
!       isflag:   Flag identifying source of solar radiation data in
!                 weather file.
!
!       hr:       Time of day (h)
!       t:        Outside air dry-bulb temperature (C)      -- toa
!       w:        Outside air humidity ratio (-)            -- woa
!       p:        Outside air barometric pressure (kPa)     -- poa
!       v:        Outside wind speed (m/s)                  -- vwind
!       soldir:   Direct normal solar radiation (W/m**2)    -- hdn
!       solsky:   Sky diffusive radiation (W/m**2)          -- hsky
!       solhoz:   Total horizontal radiation (W/m**2)       -- hhor
!
!     SUBPROGRAMS CALLED:
!       spline, speval
!     REFERENCES:
!
!
! **********************************************************************

      use modsim_head
      implicit none

      integer,parameter              :: nwdata=25
      logical                        :: sw1=.false.,eof=.false.
      character(len=3),dimension(12) :: mname= &
                                       (/'Jan','Feb','Mar','Apr','May','Jun',&
                                         'Jul','Aug','Sep','Oct','Nov','Dec'/)
      character(len=28),dimension(3) :: source= &
                                       (/'Weather tape                ',&
                                         'Clear sky design day method ',&
                                         'Cloudy sky design day method'/)

      integer                        :: isshel,i,idstr,n,isflag,iday1,&
                                        nday,iday,nw,ii,month,day,ios

      real                           :: dtr=0.0174533,time_old,&
                                        tzilch,hour,tzn,along,alat,b,&
                                        decl,hss,hsr,h1,h2,h,eqtm,&
                                        cosh,coshs,cosw,cosz,sinl,cosl,&
                                        sind,cosd
      real,dimension(nwdata)         :: hr,t,w,p,v,soldir,solsky,solhoz,&
                                        fdpt,fdpw,fdpp,fdpv,fdpd,fdps,fdph
      integer,dimension(12)          :: mdays=(/0,31,59,90,120,151,181,212,&
                                                243,273,304,334/)

!     At initial time, open the weather, and CTF files.
!     read CTF data.

      if(.not.sw1) then
        sw1=.true.
        tzero=tzilch
        hour=time/3600.+tzero
        do i=1,maxstr
          read(ifile7,*,end=20) idstr,nctfa(i),norda(i),u(i)
          read(ifile7,*) (xc(n,i),n=0,nctfa(i))
          read(ifile7,*) (yc(n,i),n=0,nctfa(i))
          read(ifile7,*) (zc(n,i),n=0,nctfa(i))
          read(ifile7,*) (rj(n,i),n=1,norda(i))
        enddo

!     Write weather-related information to simulation summary file

20      write(ifile3,1000) isshel,tshell
        read(ifile8,*) month,iday1,alat,along,tzn,isflag
        write(ifile3,2000) alat,along,iday1,mname(month),source(isflag)

        sinl=sin(alat*dtr)
        cosl=cos(alat*dtr)
        nday=mdays(month)+iday1
        b=dtr*360.*(nday-81)/364.
        eqtm=(9.87*sin(2*b)-7.53*cos(b)-1.5*sin(b))/60.
        decl=dtr*23.45*sin(dtr*360.*(284+nday)/365.)
        sind=sin(decl)
        cosd=cos(decl)

!     Read  first set of weather data.

        do i=1,nwdata
          ii=i
          read(ifile8,*,end=999) month,day,hr(i),t(i),w(i),p(i),v(i),&
                             soldir(i),solsky(i),solhoz(i)
        enddo

        hr(nwdata) = 24.0         ! added on July 31, 2000

        print *,'--  first weather data set has been read '
        goto 55
      endif

!     Compute solar declination, equation of time, and sunset hour
!     angle, once a day.

      hour=time/3600.
      if(hour+tzilch>=24.) then
        tzilch=tzilch-24.
        nday=nday+1
        b=dtr*360.*(nday-81)/364.
        eqtm=(9.87*sin(2*b)-7.53*cos(b)-1.5*sin(b))/60.
        decl=dtr*23.45*sin(dtr*360.*(284+nday)/365.)
        sind=sin(decl)
        cosd=cos(decl)
        coshs=-sinl*sind/(cosl*cosd)
        hss=acos(coshs)
        hsr=-hss
      endif
      hour=hour+tzero
      if(eof .and. hour>hr(nwdata)+0.001) goto 998

!     Read another set of weather data if necessary

      if(eof .or. (hour < hr(nwdata-1))) goto 70
      iday=0
      if(hr(nwdata-3)>24.) iday=24.
      do i=1,4
        n=nwdata-4+i
        hr(i) =hr(n)-iday
        t(i)  =t(n)
        w(i)  =w(n)
        p(i)  =p(n)
        v(i)  =v(n)
        soldir(i)=soldir(n)
        solsky(i)=solsky(n)
        solhoz(i)=solhoz(n)
      enddo
      do i=5,nwdata
        ii=i
        read(ifile8,*,end=999) month,day,hr(i),t(i),w(i),p(i),v(i),&
                             soldir(i),solsky(i),solhoz(i)
!        if(hr(i)<hr(i-1)) hr(i)=hr(i)+24.  !  MAG comment out this line
      enddo

!     8 lines are added  10/18/02 MAG
      if(hr(1) > 24) then
         hr(1) = hr(1) - 24
      endif
      do i =2,nwdata
         if(hr(i)<hr(i-1)) then
            hr(i)=hr(i)+24.
         endif
      enddo

55    nw=nwdata-1                           ! changed 10/18/02
60    call spline(nw,hr,t,fdpt)
      call spline(nw,hr,w,fdpw)
      call spline(nw,hr,p,fdpp)
      call spline(nw,hr,v,fdpv)
      call spline(nw,hr,soldir,fdpd)
      call spline(nw,hr,solsky,fdps)
      call spline(nw,hr,solhoz,fdph)

      if(hour>=hr(nw)) then
        hour=hour-24.
        tzero=tzero-24.
      endif

!     Interpolate the weather data.

70    call speval(nw,hr,t,fdpt,hour,toa)
      call speval(nw,hr,w,fdpw,hour,woa)
      call speval(nw,hr,p,fdpp,hour,poa)
      call speval(nw,hr,v,fdpv,hour,vwind)
      call speval(nw,hr,soldir,fdpd,hour,hdn)
      call speval(nw,hr,solsky,fdps,hour,hsky)
      call speval(nw,hr,solhoz,fdph,hour,hhor)

      if(woa<0.)   woa   =0.
      if(poa<0.)   poa   =0.
      if(vwind<0.) vwind =0.
      if(hdn<0.)   hdn   =0.
      if(hsky<0.)  hsky  =0.
      if(hhor<0.)  hhor  =0.
      if(hour>=24.) hour=hour-24.

!     Solar geometry is calculated half a time step ago, to represent
!     the average solar position over the time step.  if sunrise or
!     sunset occurs during the time step, exclude the night-time
!     portion of the time step.

      h1=(15.*(hour-12.+eqtm+tzn)-along)*dtr
      h2=h1-dtr*tstep/240.
      if(h2<hsr .and. h1>hsr) h2=hsr
      if(h1>hss .and. h2<hss) h1=hss
      h=0.5*(h1+h2)

!     Calculate the directional cosines

      cosh=cos(h)
      cosz=sind*sinl+cosd*cosl*cosh
      cosw=cosd*sin(h)

!     Determine the sines of solar altitude and azimuth angles.

      ssalt=cosz
      ssazm=cosw/sqrt(1.-ssalt*ssalt)

!     Make sure the sun is turned off at night.

      if(cosz<0.) then
        hdn=0.
        hsky=0.
        hhor=0.
      endif

!     Copy interpolated weather data into state vector

      if(ntoa>0 .and. ndent(3)>=0) then
        n=ndent(3)+ntoa
        state(n)=toa
        tstate(n)=toa
      endif
      if(nwoa>0 .and. ndent(8)>=0) then
        n=ndent(8)+nwoa
        state(n)=woa
        tstate(n)=woa
      endif
      if(npoa>0 .and. ndent(1)>=0) then
        n=ndent(1)+npoa
        state(n)=poa
        tstate(n)=poa
      endif
      if(ndent(7)>=0) then
        if(ndn>0) then
          n=ndent(7)+ndn
          state(n)=hdn
          tstate(n)=hdn
        endif
        if(nsky>0) then
          n=ndent(7)+nsky
          state(n)=hsky
          tstate(n)=hsky
        endif
        if(nhor>0) then
          n=ndent(7)+nhor
          state(n)=hhor
          tstate(n)=hhor
        endif
      endif

      return

1000  format(' Building shell model in superblock ',i2,':'/&
             '   Constant time step    tshell = ',f8.2/)
2000  format(' Weather data:  Latitude =',f8.3,'  Longitude =',f8.3/&
             ' Starting date: ',i4,1x,a3,'.'/&
             ' Source:  ',a28/)

999   if(ii<=4) then
        print *,' *****  Error:  Insufficient weather data  *****'
        print *,'        At least four records are required'
        goto 9999
      endif
      eof=.true.
      nw=ii-1
      goto 60

998   print *,' *****  Error:  Insufficient weather data  *****'
      print *,'        Simulation time exceeds time of last record in'
      print *,'        weather data file'

9999  stop ' ********  Simulation terminated  ********'
      end subroutine rdenv
      
! **********************************************************************

      SUBROUTINE SPLINE(N,X,Y,FDP)

! ----------------------------------------------------------------------
!
!     SPLINE : computes the second derivatives needed in cubic
!             spline interpolation. The original program was written by
!              J. H. Ferziger Ref.[1]. A little modification is made.
!
!     July 18, 1984 C.P.
!
!     n:         Number of data points
!     x:         Array containing the values of the independent variable
!                (Assume to be in ascending order)
!     y:         Array containing the values of the function at the data
!                points given in the X array
!
!     fdp:       Output array which contains the second derivatives of
!                the interpolating cubic spline.
!
!     REFERENCE:
!       [1]  Joel H. Ferziger, "Numerical Methods for Engineering
!                    Application," John Wiley & Sons, 1981, pp.17-18.
!
! **********************************************************************

     implicit none

     integer                    :: i,n
     integer,parameter          :: nmax=25
     real                       :: lamda,t
     real,dimension(nmax)       :: x,y,a,b,c,r,fdp

!     Compute the coefficients and the RHS of the equations.
!     This routine uses the cantilever condition. The parameter
!     LAMDA is set to 1. But this can be user-modified.

!     A,B,C are the three diagonals of the tridiagonal system,
!     and R is the right hand side.

      lamda = 1.
      c(1)=x(2)-x(1)
      do i=2,n-1
        c(i)=x(i+1)-x(i)
        a(i)=c(i-1)
        b(i)=2.*(a(i)+c(i))
        r(i)=6.*((y(i+1)-y(i))/c(i)-(y(i)-y(i-1))/c(i-1))
      enddo
      b(2)=b(2)+lamda*c(1)
      b(n-1)=b(n-1)+lamda*c(n-1)

!     Tridiagonal solver subroutine.
!     but the notaion is clumsy so we will solve directly.

      do i=3,n-1
        t=a(i)/b(i-1)
        b(i)=b(i)-t*c(i-1)
        r(i)=r(i)-t*r(i-1)
      enddo
      fdp(n-1)=r(n-1)/b(n-1)
      do i=2,n-2
        fdp(n-i)=(r(n-i)-c(n-i)*fdp(n-i+1))/b(n-i)
      enddo
      fdp(1)=lamda*fdp(2)
      fdp(n)=lamda*fdp(n-1)

      return
      end
! **********************************************************************

      subroutine speval(n,x,y,fdp,xx,f)

! ----------------------------------------------------------------------
!
!     speval : evaluates the cubic spline for given the derivatives
!              computed by subroutine SPLINE.
!
!     xx:        Value of independent variable for which an interpolated
!                value is reqested
!     f:         The interpolated result
!
! **********************************************************************

      implicit none

      integer                    :: i,n
      integer,parameter          :: nmax=25
      real                       :: dxm,dxp,del,f,xx
      real,dimension(nmax)       :: x,y,fdp

!     Range checking by Mike
	if( xx<=x(1) ) xx = xx + 24
	if( xx>x(nmax-1) ) xx = xx - 24

	if(( xx<=x(1) ).or.(xx>x(nmax-1))) then
		print *,"shell-speval  ### range error ###  xx = ",xx,&
     	" x(1) = ",x(1)," x(max-1) = ",x(nmax-1)
	endif

!     Find the proper interval.

      do i=1,n-1
        if(xx<=x(i+1)) goto 20
      enddo

!     Evaluate the cubic.

20    dxm=xx-x(i)
      dxp=x(i+1)-xx
      del=x(i+1)-x(i)
      f=fdp(i)*dxp*(dxp*dxp/del-del)/6.&
        +fdp(i+1)*dxm*(dxm*dxm/del-del)/6.&
        +y(i)*dxp/del+y(i+1)*dxm/del

      return
      end

