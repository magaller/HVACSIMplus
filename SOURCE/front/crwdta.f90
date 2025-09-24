!  **********************************************************************
!
!      crwdta : Create the weather data file called weather.dat.
!
!  ----------------------------------------------------------------------
!
!      Version 1.0:  Jan. 23, 1985 C.P.
!      Version 2.0:  Aug. 13, 1985 D.R.C. and J.L.
!      Version 3.0:  Dec.  1, 1986 C.P.
!      Version 4.0:  April 27, 1989 C.P.
!      Version 5.0:  July 28, 2000 Cheol Park, NIST

!      Version 6.0:  October 17, 2006 C.P.
!                    Converted Fortran 77 code into Fortran 90.
!
!      variables:
!        month:    Month of year
!        latd:     Latitude (degrees)
!        long:     Longitude (degrees)
!        tzn:      Time zone number (-)
!                  =4 --- Atlantic
!                  =5 --- Eastern
!                  =6 --- Central
!                  =7 --- Mountain
!                  =8 --- Pacific
!
!        hr:       Time of day (h)
!        t:        Dry-bulb outdoor temperature (C)
!        w:        Humidity ratio (-)
!        p:        Barometri!  pressure (kPa)
!        v:        Wind speed (m/s)
!        c:        Sky clearness number [0,1], (-)
!        soldir:   Direct solar radiation (W/m**2)
!        solsky:   Sky diffusive radiation (W/m**2)
!        solhoz:   Total horizontal radiation (W/m**2)
!
!      subprograms called:
!        wtpinp, solar
!
!  **********************************************************************

      module crwdta_comm
      implicit none
      real                       :: latd,long,tzn
      integer                    :: isflag,day,month
      integer                    :: inp=5,ifile8=8,ifile9=9
      integer,dimension(12)      :: mday=(/0,31,59,90,120,151,181,&
                                           212,243,273,304,334/)
      real,dimension(24)         :: hr,t,w,p,v,soldir,solsky,solhoz
      real,dimension(24,14)      :: wdta
      real                       :: twopi=6.2832,dtr=0.0174533,&
                                    gsc=1367.0,smallz=0.1

      end module crwdta_comm

!  **********************************************************************

      program   crwdta
      use crwdta_comm
      implicit none

      character(len=20)           :: fname ='            ',&
                                     fname1='hvacsim.met ',&
                                     fname2='            '
!      Read input information.

      print *, '   ************************************************* '
      print *, '   *                                               * '
      print *, '   *         Creating a weather data file          * '
      print *, '   *                                               * '
      print *, '   ************************************************* '
      print *, '  '

      print *, ' Enter latitude, longitude, and time zone:'
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,*) latd,long,tzn

5     print *,' Enter one of the following:'
      print *,' 1 - to process the weather data in file wtpout.dat'
      print *,'   (previously read from weather file by program rdwdf)'
      print *,' 2 - to generate clear sky design data'
      print *,' 3 - to generate cloudy sky design data'
      print *
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,*) isflag
      if(isflag<1.or.isflag>3) goto 5

!      Give name to output file.

      print *,' Enter output file name (up to 40 characters)'
      print *,' or carriage return for default name: ',fname1
      write(*,fmt='(a4)',advance='no') ' => '
      read(unit=inp,fmt='(a20)') fname

      if(fname==fname2) then
        fname=fname1
      endif

25    open(unit=ifile8,file=fname,status='replace')

      if(isflag==1) then
        call wtpinp
      else
        call solar
      endif

      stop  '---- End of creating weather file ------'
      end program crwdta
!  **********************************************************************
!
!      wtpinp : Read the output file of rdwdf program, wtpout.dat.
!
!  ----------------------------------------------------------------------
!      Jan. 23, 1985 C.P.
!      Updated: July 26, 1985 D.C.R. and J.L.
!      Updated: Dec. 1, 1986 C.P.
!      Updated: July 28, 2000  Cheol Park
!      Updated: October 17, 2006 C.P.
!               Converted Fortran 77 code into Fortran 90.
!
!      inputs:
!       See descriptions in the program, rdwdf.
!
!      local variables:
!        b:      Intermediate result for E
!        e:      Equation of time (h)
!        decl:   Declination (rads)
!        gzeron: Extraterrestial normal radiation (W/m**2)
!        coshs:  Cosine of the sunset angle
!        hrss:   Hour of sunset (h)
!        hrsr:   Hour of sunrise (h)
!        tss:    Time of sunset for a particular location (h)
!        tsr:    Time of sunrise for a particular location (h)
!        h:      Current hour angle (rads)
!        hlast:  Previous hour angle (rads)
!        hbar:   Average of the current and previous hour angles (rads)
!        cosz:   Cosine of the solar zenith angle
!        izero:  Intergrated extraterrestial radiation (W/m**2)
!        akt:    Hourly clearness index
!
!      subprogram called:
!        wf, copyfl
!
!      reference:
!      [1] Erbs, D.G., Klein S.A. and Duffie J.A., "Estimation of the
!                  Diffuse Radiation Fraction for Hourly, Daily, and
!                  Monthly-Average Global Radiation", Solar Energy
!                  Laboratory, University of Wisconsin, Madison, 1981
!
!  **********************************************************************

      subroutine wtpinp

      use crwdta_comm
      implicit none
      logical                         :: first=.true.,code=.true.
      integer                         :: sttn,wyr
      integer,dimension(24)           :: whr
      integer,dimension(24,8)         :: idta
      real                            :: izero
      character(len=132),dimension(2) :: head
      real                            :: sinl=0.0,cosl=1.0
      integer                         :: jflag=1,ioerr
      integer                         :: j,nday,ihrss,ihrsr
      real                            :: wf,sumdir,sumsky,sumhoz,b,e,&
                                         decl,gzeron,sind,cosd,coshs,&
                                         hrss,hrsr,adjust,tsr,tss,h,hlast,&
                                         hbar,cosz,akt


      open(unit=ifile9,file='wtpout.dat',status='old',iostat=ioerr)
      if(ioerr/=0) then
        print *, ' File open error'
        return
      end if
      read(unit=ifile9,fmt='(a132/a132)') head(1),head(2)
10    do j=1,24
        read(unit=ifile9,fmt=200,iostat=ioerr)&
          sttn,wyr,month,day,whr(j),t(j),wdta(j,11),p(j),&
          v(j),idta(j,1),wdta(j,8),wdta(j,9),wdta(j,10),wdta(j,12)
        200 format(1x,i5,4i5,f10.2,3f10.1,i5,4f10.2)
        if(ioerr>0) then
          print *,' Error in reading input data'
          return
        elseif(ioerr<0) then
          print *,' End of input data file '
          return
        endif

        w(j)     =wf(wdta(j,11),p(j))
        hr(j)    =whr(j)

        if(wdta(j,9)<=0) then
          soldir(j)=0.0
        else
          soldir(j)=wdta(j,9)
        endif

        if(wdta(j,10)<=0) then
          solsky(j)=0.0
        else
          solsky(j)=wdta(j,10)
        endif

        if(wdta(j,12)<=0) then
          solhoz(j)=0.0
        else
          solhoz(j)=wdta(j,12)
        endif

      end do

!      Examine which solar components are available from wtpout.dat file
!      and set jflag.

      if(first) then

        cosl=cos(latd*dtr)
        sinl=sin(latd*dtr)

        sumdir=0.0
        sumsky=0.0
        sumhoz=0.0

        do j=1,24
          sumdir=sumdir+soldir(j)
          sumsky=sumsky+solsky(j)
          sumhoz=sumhoz+solhoz(j)
        end do

              jflag=1
        if(sumhoz>1.0) then
              jflag=2
          if(sumsky>1.0) then
              jflag=3
            if(sumdir>1.0) then
              jflag=4
            endif
          endif
        elseif(sumdir>1.0 .or. sumsky>1.0) then
              jflag=5
        endif
        first=.false.
      endif

      hr(24)=24.0                                !   10/16/2006
      nday=mday(month)+day
      b=twopi*(nday-81)/364.
      e=(9.87*sin(2*b)-7.53*cos(b)-1.5*sin(b))/60.

      decl=23.45*sin(twopi*(284+nday)/365.0)*dtr
      gzeron=gsc*(1+0.033*cos(twopi*nday/365.0))
      sind=sin(decl)

      cosd=cos(decl)
      coshs=-sinl*sind/(cosl*cosd)

      hrss=acos(coshs)/(15.0*dtr)
      hrsr=-hrss
      adjust=(long/15.0)+12-e-tzn
      tss=hrss+adjust
      tsr=hrsr+adjust
      ihrsr=tsr
      ihrss=tss+1

!      If sky or direct radiation information is missing, then calculate
!      their value using the total horizontal radiation.

      if(jflag==2 .or. jflag==3) then
        do j=1,24

!          Check if hour angle is greater than sunset angle.

          if(j<=ihrsr.or.j>ihrss) then
            soldir(j)=0.
            solsky(j)=solhoz(j)
          else
            h=(15.0*(hr(j)-12.+e+tzn)-long)*dtr
            hlast=h-15.0*dtr

            if(cos(hlast)<coshs) hlast=-acos(coshs)
            if(cos(h)<coshs) h=acos(coshs)

            hbar=0.5*(h+hlast)
            cosz=sind*sinl+cosd*cosl*cos(hbar)

            if(cosz <= smallz) then                             ! 11/4/02
               cosz = smallz
            endif

            izero=(24./twopi)*gzeron*(cosl*cosd*(sin(h)-sin(hlast))+&
              (h-hlast)*sinl*sind)

            if(izero>0)then
              akt=solhoz(j)/izero
            else
              akt=0.
            endif

            if(solhoz(j)==0) then
              solsky(j)=0.0
              soldir(j)=0.0

            else
              if(soldir(j)>0.and.solsky(j)==0) then
                solsky(j)=solhoz(j)-soldir(j)*cosz

              elseif(solsky(j)>0.and.soldir(j)==0) then
                soldir(j)=(solhoz(j)-solsky(j))/cosz

              elseif(soldir(j)==0.and.solsky(j)==0) then
 
!              Hourly diffuse correlation of erbs ref.[1].

                if(akt<=0.22) then
                  solsky(j)=solhoz(j)*(1.0 -0.09*akt)

                else if(akt>0.22.and.akt<0.80) then
                  solsky(j)=solhoz(j)*(0.9511+akt*(-0.1604+akt*(4.388+&
                    akt*(-16.638+akt*(12.336)))))

                else if(akt>=0.80) then
                  solsky(j)=solhoz(j)*0.165
                endif

                soldir(j)=(solhoz(j)-solsky(j))/cosz

              endif
            endif
          endif
       end do

!      If the total horizontal radiation is missing, then calculate
!      the value using the direct and diffuse radiation data.

      elseif(jflag==5) then
        do j=1,24
          h=(15.0*(hr(j)-12.+e+tzn)-long)*dtr
          hlast=h-15.0*dtr
          if(cos(hlast)<coshs) hlast=-acos(coshs)
          if(cos(h)<coshs) h=acos(coshs)
          hbar=0.5*(h+hlast)
          cosz=sind*sinl+cosd*cosl*cos(hbar)

          if(cosz <= smallz) then                             ! 11/4/02
             cosz = smallz
          endif

          solhoz(j)=soldir(j)*cosz+solsky(j)

        end do
      endif

!      Write 24 hours of data.

      call copyfl(code)

      goto 10

      return
      end subroutine wtpinp

!  **********************************************************************

!     wf: Compute humidity ratio using dew point temperature.
!         This rouitne is based on the similar routine in Ref.[2].

!  ----------------------------------------------------------------------
!
!      reference:
!      [1] Brokaw, R.S., "Calculation of Flue Losses for High-efficiency
!                        Furnaces and Appliances," ASHRAE J., Jan. 1979,
!                        pp. 49-51.
!      [2] Park, C.,Kelly,G.E., and Kao, J.Y.,"Economizer Algorithms for
!                        Energy Management and Control Systerms,"
!                        NBSIR 84-2832, Feb. 1984, p.B-3.
!
!  **********************************************************************

      real function wf(dp,p)
      implicit none
      real                :: dp,p,psw
      psw=3.376*exp(15.463-7284./(1.8*dp+424.0))
      wf=0.62198*psw/(p-psw)

      return
      end function wf
!  **********************************************************************
!
!      solar : Create clear sky or cloudy sky weather data.
!
!  ----------------------------------------------------------------------

!      July 31, 1985 D.R.C. and J.L.
!      Updated: October 17, 2006 C.P.
!               Converted Fortran 77 code into Fortran 90.
!
!      local variables:
!        ac:   Apparent solar constant (W/m**2)
!        bc:   Exponential attenuation (decay) coefficient
!        cc:   Diffuse fraction factor
!        taua: Transmittance
!        bw:   Correction factor
!        For other variables see subroutine wtpinp
!
!      subprograms called:
!        db,copyfl
!
!      reference:
!      [1] Machler, M.A. & Iqbal, M. , "A Modification of the ASHRAE
!                  Clear Sky Irradiation Model", ASHRAE Trans. Vol.91,
!                  Part 1, 1985
!      [2] Duffie & Beckman ,"Solar Engineering of Thermal Processes"
!                  p.79
!      [3] Erbs, D.G., Klein S.A. and Duffie J.A., "Estimation of the
!                  Diffuse Radiation Fraction for Hourly, Daily, and
!                  Monthly-Average Global Radiation", Solar Energy
!                  Laboratory, University of Wisconsin, Madison, 1981
!      [4] Erbs D.G., Klein S.A. and Beckman W.A., "A Simple Distribution
!                  Method for Two-Dimensonal Temperature/Humidity Bin
!                  Data", ASHRAE Transactions 1985, Vol.91, Pt.2
!
!  **********************************************************************

      subroutine solar

      use crwdta_comm
      implicit none
      logical               :: code=.true.
      integer               :: n=1
      real,dimension(12)    :: &
               ac=(/ 1202.0, 1187.0, 1164.0, 1130.0, 1106.0, 1092.0,&
                     1093.0, 1107.0, 1136.0, 1166.0, 1190.0, 1204.0/),&
               bc=(/ 0.141, 0.142, 0.149, 0.164, 0.177, 0.185,&
                     0.186, 0.182, 0.165, 0.152, 0.144, 0.141/),&
               cc=(/ 0.103, 0.104, 0.109, 0.120, 0.130, 0.137,&
                     0.138, 0.134, 0.121, 0.111, 0.106, 0.103/)

      integer,dimension(12) :: moday=(/31,28,31,30,31,30,31,31,30,31,30,31/)
      real                  :: cor=1.0,vis=0.0,sw=0.0,ppo=1.0,kbar=0.5,izero

      integer                         :: j,nday,ihrss,ihrsr,noday,iday,&
                                         i,lday,kday,k,jday
      real                            :: b,e,decl,gzeron,sind,cosd,coshs,&
                                         hrss,hrsr,adjust,tsr,tss,h,hlast,&
                                         hbar,cosz,akt,rh,vc,pc,tmax,tmin,&
                                         cosl,sinl,sa,sb,frac,ams,expo,&
                                         taux,taua,bwx,bw,tauax,rh_percent

      print *,' Enter initial day and month, and number of days'
      print *,' for which weather calculations will be made'
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,*) iday,month,noday

      print *,' Enter pressure (kpa), wind speed (m/s), and'
      print *,' relative humidity (%)'
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,*) pc,vc,rh_percent
      rh = rh_percent/100.0

      print *,' Enter minimum and maximum temperatures (C):'
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,*) tmin,tmax

      do i=1,24
        p(i)=pc
        v(i)=vc
      end do

      if(isflag==2) then
        print *,' Enter visibility (km); if value unknown, use 0:'
        write(*,fmt='(a4)',advance='no') ' => '
        read(inp,*) vis
        if(vis==0) then
          n=1
          print *,' Enter geographic correction factor'
          print *,' [ASHRAE Fund. 1981, p.27.8]; if value unknown, ',&
                  ' use 1:'
          write(*,fmt='(a4)',advance='no') ' => '
          read(inp,*) cor
        else
          n=2
          print *,' Enter precipitable water (mm); if value unknown, ',&
                  ' use 0:'
          write(*,fmt='(a4)',advance='no') ' => '
          read(inp,*) sw
          if(sw/=0.) n=3
        endif
        ppo=pc/101.3

      else
        print *,'Enter daily clearness index:'
        write(*,fmt='(a4)',advance='no') ' => '
        read(inp,*) kbar
      endif

!      Determine if the number of days is equal to number of days in
!      initial month.

      lday=iday+noday-1
2     if(lday<=moday(month)) then
        kday=lday
        k=1
      else
        kday=moday(month)
        k=2
      endif

      day_loop : do i=iday,kday

        day=i
        cosl=cos(latd*dtr)
        sinl=sin(latd*dtr)
        nday=mday(month)+day
        b=twopi*(nday-81)/364.0

        e=(9.87*sin(2.0*b)-7.53*cos(b)-1.5*sin(b))/60.0
        decl=23.45*sin(twopi*(284+nday)/365.0)*dtr
        gzeron=gsc*(1.0+0.033*cos(twopi*nday/365.0))
        sind=sin(decl)
        cosd=cos(decl)

        coshs=-sinl*sind/(cosl*cosd)
        sa=0.409+0.5016*sin((acos(coshs))-1.0472)
        sb=0.6609- 0.4767*sin((acos(coshs))-1.0472)
        hrss=acos(coshs)/(15.0*dtr)
        hrsr=-hrss

        adjust=(long/15.0)+12.0-e-tzn
        tss=hrss+adjust
        tsr=hrsr+adjust
        ihrsr=tsr
        ihrss=tss+1
  
        hour_loop : do j=1,24
          hr(j)=j                                           ! 10/16/2006
          h=(15.0*(hr(j)-12.+e+tzn)-long)*dtr
          hlast=h-15.0*dtr
          if(cos(h)<coshs) h=acos(coshs)
          if(cos(hlast)<coshs) hlast=-acos(coshs)
          hbar=0.5*(h+hlast)
          cosz=sind*sinl+cosd*cosl*cos(hbar)

          if(cosz <= smallz) then                             ! 11/4/02
             cosz = smallz
          endif

          izero=(24.0/twopi)*gzeron*(cosl*cosd*(sin(h)-sin(hlast))+&
            (h-hlast)*sinl*sind)

!      Determine the fraction of the extraterrestial radiation for
!      sunrise and sunset.

          if(hr(j)==(ihrsr+1)) then
            frac=1.0-(tsr-ihrsr)
            izero=izero*frac
          else if(hr(j)==ihrss) then
            frac=tss-(ihrss-1)
            izero=izero*frac
          endif

          if(hr(j)<=ihrsr.or.hr(j)>ihrss) then
            solhoz(j)=0.0
            solsky(j)=0.0
            soldir(j)=0.0
          else
            if(isflag==2) then
              ams=ppo/cosz

!      Clear sky design day calculations ref.[1].

              if(n==1) then
                if(cosz<=1.0e-3) then
                  expo=0.0
                  soldir(j)=0.0
                  solsky(j)=0.0
                else
                  expo=exp(-bc(month)*ams)
                  soldir(j)=ac(month)*expo*cor
                  solsky(j)=soldir(j)*cc(month)
                endif

              else
                tauax=1.0 -1.13*vis**(-0.57)
                if(tauax <= 0.0) then
                  stop ' -- input data not reasonable ---'
                endif
                taua=tauax**(ams**0.85)
                soldir(j)=gzeron*taua*0.775**(ams**0.50)
                if(n==3) then
                  bwx=1.0223 -0.00149*sw
                  if(bwx <= 0.0) then
                    stop ' -- input data not reasonable ---'
                  endif                                  
                  bw=bwx**(ams**0.27)
                  soldir(j)=soldir(j)*bw
                endif
                solsky(j)=soldir(j)*(0.1+3.0/vis)
              endif
              solhoz(j)=soldir(j)*cosz+solsky(j)

!      Cloudy sky design day calculations.

            else

!      Equation for hourly clearness index akt, is derived from
!      relationship given in ref.[2].

              akt=kbar*(sa+sb*cos(hbar))
              solhoz(j)=akt*izero

!      Hourly diffuse correlation of erbs ref.[3].

              if(akt<=0.22) then
                solsky(j)=solhoz(j)*(1.0 -0.09*akt)

              else if(akt>0.22.and.akt<0.80) then
                solsky(j)=solhoz(j)*(0.9511+akt*(-0.1604+akt*(4.388+&
                    akt*(-16.638+akt*(12.336)))))

              else if(akt>=0.80) then
                solsky(j)=solhoz(j)*0.165
              endif
              soldir(j)=(solhoz(j)-solsky(j))/cosz
            endif

          endif
        end do hour_loop

!      Calculate dry-bulb temperature distribution.

        call db(tmax,tmin)

!      Calculate the humidity ratio, assuming constant relative
!      humidity (as suggested by ref. [4]).

        call humidy(rh,pc)

!      Write weather information in output file.

        call copyfl(code)

      end do day_loop

!      If number of days exceed days contained in the initial month,
!      then proceed with the following month.

      if(k==1) then
        return
      else
        jday=moday(month)-(iday-1)
        lday=noday-jday                
        month=month+1                  
        iday=1                         
        go to 2
      endif

      return
      end subroutine solar
!  **********************************************************************
!
!      db : Approximated dry-bulb temperature calculation
!           using empirical values shown in the table, Ref[1].
!           This result is valid for summer time only.
!
!  ----------------------------------------------------------------------
!      November 14, 1984 C.P.
!      Updated July 30, 1985 J.L.
!      Updated: October 17, 2006 C.P.
!               Converted Fortran 77 code into Fortran 90.
!
!      inputs:
!        range:  Daily temperature range (C)
!        percen: Percentage of the daily range (-)
!                from Table 3, p.26.6, ASHRAE handbook 1981
!
!      reference:
!      [1] ASHRAE Fundamentals Handbook 1981, p.26.6
!
!  **********************************************************************

      subroutine db(tmax,tmin)

      use crwdta_comm
      implicit none
      real,dimension(24)  :: percent=(/87.,92.,96.,99.,100.,98.,93.,84.,&
                                       71.,56.,39.,23.,11., 3., 0., 3.,&
                                       10.,21.,34.,47.,58.,68.,76.,82./)
      real                :: tmax,tmin,range
      integer             :: i

      range=tmax-tmin
      do i=1,24
        t(i)=tmax-range*percent(i)/100.0
      end do

      return
      end subroutine db

!  **********************************************************************
!
!      humidy : Calculate the humidity ratio for 24 hours,
!                 assuming constant relative humidity.
!
!  ----------------------------------------------------------------------
!      Aug. 12, 1985 D.R.C. and J.L.
!      Updated: October 17, 2006 C.P.
!               Converted Fortran 77 code into Fortran 90.
!
!      inputs:
!        rh: Relative humidity
!        pc: Atmospheric  presure (kPa)
!        t:  Dry-bulb temperature (C)
!
!      local variables:
!        at:  Absolute temperature (K)
!        pws: Saturation pressure (Pa)
!
!      outputs:
!        w: Humidity ratio
!
!      subprograms called:
!        None
!
!      reference:
!      [1] ASHRAE Fundamenals Handbook 1985, p.6.4
!
!  **********************************************************************

      subroutine humidy(rh,pc)

      use crwdta_comm
      implicit none
      real             :: c1=-5674.5359,    c2=6.3925247,&
                          c3=-0.9677843e-2, c4=0.62215701e-6,&
                          c5=0.20747825e-8, c6=0.9484024e-12,&
                          c7=4.1635019,     c8=-5800.2206,&
                          c9=1.3914993,     c10=-0.04860239,&
                          c11=0.41764768e-4,c12=-0.14452093e-7,&
                          c13=6.5459673
      real             :: at,pws,rh,pc
      integer          :: j

      do j=1,24
        at=t(j)+273.15
        if(t(j)<=0.0) then
          pws=exp(c1/at+c2+at*(c3+at*(c4+at*(c5+at*c6)))+c7*alog(at))
        else
          pws=exp(c8/at+c9+at*(c10+at*(c11+at*c12))+c13*alog(at))
        endif
        w(j)=(0.62198*rh*pws)/((pc*1000.0)-rh*pws)
      end do
      return
      end subroutine humidy

!  **********************************************************************
!
!      copyfl : Write 24 hours of data.
!
!  ----------------------------------------------------------------------
!      Aug. 2, 1985 D.R.C. and J.L.
!      July 28, 2000  C.P.
!      Updated: October 17, 2006 C.P.
!               Converted Fortran 77 code into Fortran 90.
!

!      subprograms called:
!        None
!
!  **********************************************************************

      subroutine copyfl(code)

      use crwdta_comm
      implicit none
      logical             :: code
      integer             :: k
      real                :: hr_zero

!      Write title.

      if(code) then                              !   7/28/00
        write(ifile8,100)month,day,latd,long,tzn,isflag
        100 format(1x,2i5,3f10.2,i5)

!       Write information corresponding to the first hour of
!       first day as initial values for modsim.

        hr_zero = 0.0                          ! 10/17/2006
        write(ifile8,200)month,day,hr_zero,t(1),w(1),p(1),v(1),&
            soldir(1),solsky(1),solhoz(1)
        200 format(1x,2i5,f5.1,7f10.4)                             
          code=.false.                                             
      endif

!     Write information corresponding to one day.
                                                                   
      do k=1,24
        write(ifile8,200)month,day,hr(k),t(k),w(k),p(k),v(k),&
          soldir(k),solsky(k),solhoz(k)
      end do

      return
      end subroutine copyfl
