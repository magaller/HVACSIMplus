!  **********************************************************************
!
!      rdwdf : Read a whole year weather data file.
!
!               This routine was made based on the TARP package writtren
!               by G.N. Walton , Ref[1].
!
!      July 31, 1984  Cheol Park, N.B.S.
!      Updated : Jan. 24, 1985 C.P.
!      Updated : Sept. 11, 2003 C.P., NIST
!                Added a capability to read WYEC2 data
!      Updated : May 27, 2008 C.P.
!                Code converted from Fortran 77 to Fortran 90, and
!                renamed to rdwdf from rdtape.
!
!  ----------------------------------------------------------------------
!
!      inputs:
!        itype:    Type of data format
!                  1 : TRY
!                  2 : TMY
!                  3 : SOLMET
!                  4 : WYEC
!                  5 : WYEC2
!        isttn:    Station number
!        iyear:    Year (e.g. 1953)
!        mostrt:   Month of start date
!        mostop:   Month of stop date
!        dystrt:   Day of start date
!        dystop:   Day of stop date
!
!      subprograms called:
!        rdwtp, jds, wrtfil
!
!      REFERENCES:
!      [1] Walton, G.N., "Thermal Analysis Research Program Reference
!          Manual," Nat'l Bureau of Standards, NBSIR 83-2655, March 1983.
!
!      [2] "SOLMET User's Manual, Volume 1 Hourly Solar Radiation-Surface
!          Meteorological Observations," NOAA National Climati!  Center,
!          Asheville, NC, August 1978.
!      [3] Hall, I.J., Prairie, R.R., Anderson, H.E., and Boes, E.C.,
!          "Generation of typical meteorological years for
!          26 SOLMET stations," ASHRAE Trans. DE-79-2,No.3,
!          1979,pp.507-518.
!      [4] "User's Manual Typical Meteorological Year Data,
!          Tape Deck 9734", Sandia Lab.,Albuquerque, NM.
!      [5] "WYEC2 User's Manual and Toolkit," ASHRAE, 1996.
!
!  **********************************************************************

      module rdwdf_head

      integer                  :: isttn,iyear,jdnext,wflag
      integer                  :: sttn,wdate,wyr,wmo,wdy
      integer,dimension(24)    :: whr
      integer,dimension(24,8)  :: idta
      real,dimension(24,14)    :: wdta
      integer                  :: ifile=7,ofile=8,inp=5
      integer                  :: ioerr,l
      
      end module rdwdf_head

!  **********************************************************************

      program   rdwdf
      use rdwdf_head
      implicit none
      integer               :: dystrt,dystop
      integer               :: nwdays,itype,mostrt,mostop,jdstrt,jdstop,iflag
      character(len=12)     :: wtpin
      data nwdays/0/

!*    namelist  /indata/ itype,isttn,iyear,mostrt,dystrt,mostop,dystop &
!*             /jddata/ jdstrt,jdstop

      print *, ' Enter input file name up to 12 characters ---'
      write(*,fmt='(a4)',advance='no') ' => '
      read (inp,fmt='(a12)') wtpin
      open(unit=ifile,file=wtpin,access='sequential',&
           form='formatted',status='old',iostat=ioerr)
      if(ioerr>0) then
        stop '****  error: file is not available ******'
      endif
      open(unit=ofile,file='wtpout.dat')
      close(unit=ofile,status='delete')
      open(unit=ofile,file='wtpout.dat')
      rewind (unit=ifile)

!      Read input data

      print *,' What is the type of weather data format? '
      print *,' Enter 1 for (TRY), 2 (TMY), 3 (SOLMET), 4 (WYEC),'
      print *,'       5 (WYEC2)'
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,*) itype
      print *,' Where is the weather station?'
      print *,' Enter station ID number '
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,*) isttn
      print *,' Enter the year (4 digits)'
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,*) iyear
      print *,' Type the start date:  Month,Day '
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,*) mostrt,dystrt
      print *,' Type the stop date:   Month,Day '
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,*) mostop,dystop
!*    print indata

!      Convert the conventional date into Julian day. Ref.[1]

      call jds(1,jdstrt,iyear,mostrt,dystrt)
      if(mostrt<=mostop) then
        call jds(1,jdstop,iyear,mostop,dystop)
      else
        call jds(1,jdstop,iyear+1,mostop,dystop)
      endif

!      Read the whole year weather data file.

!*    print jddata
      jdnext=jdstrt
10    call rdwtp(itype,jdstrt)
      if(wflag/=0) then
        stop '---- End of job due to end-of-file or error -----'
      endif
      nwdays=nwdays+1
      if(jdnext<jdstop) then
        if(jdnext==jdstrt) then
          print *,'------- The start day -------'
          print *,' sttn=',sttn,' wyr=',wyr,' wmo=',wmo,' wdy=',wdy
          iflag=1
        else
          iflag=2
        endif
        call wrtfil(iflag)
        jdnext=jdnext+1
        goto 10
      elseif(jdnext==jdstop) then
        print *,' -------  The stop day -------'
        print *,' sttn=',sttn,' wyr=',wyr,' wmo=',wmo,' wdy=',wdy
        call wrtfil(iflag)
        close (unit=ifile,status='keep')
        print 100,nwdays
        100 format(5x,i5,' Days written on the output file'/)
        stop '--------- Normal end of job ------'
      endif
      end program rdwdf

!  **********************************************************************

      subroutine rdwtp(itype,jdstrt)

!  ----------------------------------------------------------------------
!
!      rdwtp : Read weather data file.
!
!      July 31, 1984 C.P.
!      Updated :  Sept. 11, 2003
!      Updated : October 13, 2006 C.P.
!                Code converted from Fortran 77 to Fortran 90.
!
!      subprogams called:
!        rdtry, rdtmy, rdsolm, rdwyec
!
!  **********************************************************************

      use rdwdf_head
      implicit none
      integer         :: jdstrt,itype
      logical         :: first=.true.

!      Read file according to its type.

10    select case(itype)
        case(1)
          call rdtry
        case(2)
          call rdtmy
        case(3)
          call rdsolm
        case(4)
          call rdwyec
        case(5)
          call rdwyec2
        case default
      end select

!      Write the first record.

      if(first) then
        print *,' ---- The first day of the weather data ----'
        print *,' sttn=',sttn,' wyr=',wyr,' wmo=',wmo,' wdy=',wdy
        first=.false.
      endif

!      Test the position.

      if(wflag==1) then
        print *,'***** End of file ******'
        return
      elseif(wflag==-1) then
        print *,'***** Read error *******'
      else
        if(isttn/=sttn .or. jdstrt>wdate) then
          goto 10
        endif
      endif

!      Return to the main program if station and required date
!      are matched with the readout of file.

      if(isttn==sttn .and. jdnext==wdate) then
          return

!        No matched date is found until the end of file reaches.

      else
        wflag=1
        print *,'****** Error: date not found in file *******'
      endif

      return
      end subroutine rdwtp

!  **********************************************************************

      subroutine rdsolm

!  ----------------------------------------------------------------------
!
!      Read 24 hours of data from a NOAA SOLMET format file.
!      This routine is based on the similar BLAST routine by
!      G. Walton and L. Lawrie.
!      TCY (Typical Cooling Year) format data file also uses this format.
!
!      Modifications were made by C. Park.
!      Jan. 17, 1985
!      Updated : October 13, 2006 C.P.
!                Code converted from Fortran 77 to Fortran 90.
!
!      variables:
!      tmiss  - Missing value
!      cnvrj  - Unit conversion factor for radiation ( kJ/m**2 -> W/m**2)
!
!      The TMY format data positions and variables are ---- [Ref.2]
!         position       variable       field description
!         5-9          sttn           WBAN Station number
!         10-11        wyr            Year     (2 digits)
!         12-13        wmo            Month
!         14-15        wdy            Day
!         16-17        whr(hour)      Hour
!         23-27        wdta(hour,8)   Extraterrestrial radiation(kJ/m**2)
!         28-32        wdta(hour,9)   Direct radiation (kJ/m**2)
!         33-37        wdta(hour,10)  Diffuse radiation (kJ/m**2)
!         48-52        wdta(hour,14)  Observed horiz. radiation (kJ/m**2)
!         53-57        wdta(hour,13)  Engineering corrected radiation
!                                                               (kJ/m**2)
!         58-62        wdta(hour,12)  Std. year corr. horiz. radiation
!         91           idta(hour,5)   Rain indicator (If = 2,3,6, OR 8)
!         92           idta(hour,6)   Rain indicator (If = 2 OR 3)
!         103-107      wdta(hour,3)   Station pressure (kPa)
!         108-111      wdta(hour,1)   Dry bulb temperature (C)
!         112-115      wdta(hour,11)  Dew point temperature (C)
!         116-118      wdta(hour,6)   Wind direction (degree)
!         119-122      wdta(hour,5)   Wind speed (m/s)
!         123-124      idta(hour,1)   Total sky cover [0,10]
!         125-128      idta(hour,7)   Type of lowest cloud
!         163          idta(hour,4)   Snow cover indicator
!
!  **********************************************************************

      use rdwdf_head
      implicit none
      real          :: tmiss=9999.0,smiss=-1000.0,cnvrj=0.2777777777778
      integer       :: imiss=-1000

      wflag=0

!      Read 24 hour data.

      do l=1,24
        read(unit=ifile,fmt=100,iostat=ioerr)&
         sttn,wyr,wmo,wdy,whr(l),wdta(l,8),wdta(l,9),wdta(l,10),&
         wdta(l,14),wdta(l,13),wdta(l,12),idta(l,5),&
         idta(l,6),wdta(l,3),wdta(l,1),wdta(l,11),&
         wdta(l,6),wdta(l,5),idta(l,1),idta(l,7),idta(l,4)
        100 format(4x,i5,4i2,5x,3(1x,f4.0),10x,3(1x,f4.0),&
            28x,2i1,10x,f5.2,2f4.1,f3.0,f4.1,i2,2x,i2,34x,i1)
        if(ioerr>0) then                 ! data reading error
          wflag=-1
          return
        elseif(ioerr<0) then             ! end of file
          wflag=1
          return
        endif
      end do
      call jds(1,wdate,iyear,wmo,wdy)
      if(wdate/=jdnext) return

!      Convert units and fill out undefined data

      do l=1,24

!        Solar

        if(wdta(l,12)==tmiss) then
          if(wdta(l,13)/=tmiss) then
            wdta(l,12)=wdta(l,13)
          else
            wdta(l,12)=wdta(l,14)
          end if
        end if
        if(wdta(l,12)<tmiss) then
          wdta(l,12)=wdta(l,12)*cnvrj
        else
          wdta(l,12)=smiss
        end if
        if(wdta(l,8)<tmiss) then
          wdta(l,8)=wdta(l,8)*cnvrj
        else
          wdta(l,8)=smiss
        end if
        if(wdta(l,9)<tmiss) then
          wdta(l,9)=wdta(l,9)*cnvrj
        else
          wdta(l,9)=smiss
        end if
        if(wdta(l,10)<tmiss) then
          wdta(l,10)=wdta(l,10)*cnvrj
        else
          wdta(l,10)=smiss
        end if

!        Cloud amount

        if(idta(l,1)<0 .or. idta(l,1)>10) then
          idta(l,1)=imiss
        endif
      end do
      return
      end subroutine rdsolm
!  **********************************************************************

      subroutine rdtmy

!  ----------------------------------------------------------------------
!
!      Read 24 hours of data from a NOAA Typical Meteorological Year(TMY)
!      format file. This routine is based on the similar BLAST routine by
!      G. Walton and L. Lawrie.
!
!      Modifications were made by C. Park.
!      Jan. 17, 1985
!      Updated : October 13, 2006 C.P.
!                Code converted from Fortran 77 to Fortran 90.
!
!      variables:
!      tmiss  - Missing value
!      cnvrj  - Unit conversion factor for radiation ( kJ/m**2 -> W/m**2)
!
!      The TMY format data positions and variables are ---- [Ref.3 & 4]
!         position       variable       field description
!         1-5          sttn           WBAN Station number
!         6-7          wyr            Year     (2 digits)
!         8-9          wmo            Month
!         10-11        wdy            Day
!         12-13        whr(hour)      Hour
!         19-23        wdta(hour,8)   Extraterrestrial radiation(kJ/m**2)
!         24-28        wdta(hour,9)   Direct radiation (kJ/m**2)
!         29-33        wdta(hour,10)  Diffuse radiation (kJ/m**2)
!         44-48        wdta(hour,14)  Observed horiz. radiation (kJ/m**2)
!         49-53        wdta(hour,13)  Engineering corrected radiation
!                                                               (kJ/m**2)
!         54-58        wdta(hour,12)  Std. year corr. horiz. radiation
!         87           idta(hour,5)   Rain indicator (If = 2,3,6, OR 8)
!         88           idta(hour,6)   Rain indicator (If = 2 OR 3)
!         99-103       wdta(hour,3)   Station pressure (kPa)
!         104-107      wdta(hour,1)   Dry bulb temperature (C)
!         108-111      wdta(hour,11)  Dew point temperature (C)
!         112-114      wdta(hour,6)   Wind direction (degree)
!         115-118      wdta(hour,5)   Wind speed (m/s)
!         119-120      idta(hour,1)   Total sky cover [0,10]
!         121-122      idta(hour,7)   Type of lowest cloud
!         123          idta(hour,4)   Snow cover indicator
!         124-132                     Blank
!
!  **********************************************************************

      use rdwdf_head
      implicit none
      real          :: tmiss=9999.0,smiss=-1000.0,cnvrj=0.2777777777778
      integer       :: imiss=-1000

      wflag=0

!      Read 24 hour data.

      do l=1,24
        read(unit=ifile,fmt=100,iostat=ioerr)&
         sttn,wyr,wmo,wdy,whr(l),wdta(l,8),wdta(l,9),wdta(l,10),&
         wdta(l,14),wdta(l,13),wdta(l,12),idta(l,5),&
         idta(l,6),wdta(l,3),wdta(l,1),wdta(l,11)
        100 format(i5,4i2,5x,3(1x,f4.0),10x,3(1x,f4.0),&
         28x,2i1,10x,f5.2,2f4.1,f3.0,f4.1,i2,i2,i1,9x)
        if(ioerr>0) then                 ! data reading error
          wflag=-1
          return
        elseif(ioerr<0) then             ! end of file
          wflag=1
          return
        endif
      end do
      call jds(1,wdate,iyear,wmo,wdy)
      if(wdate/=jdnext) return

!      Convert units and fill out undefined data.

      do l=1,24

!        Solar

        if(wdta(l,12)==tmiss) then
          if(wdta(l,13)/=tmiss) then
            wdta(l,12)=wdta(l,13)
          else
            wdta(l,12)=wdta(l,14)
          end if
        end if
        if(wdta(l,12)<tmiss) then
          wdta(l,12)=wdta(l,12)*cnvrj
        else
          wdta(l,12)=smiss
        end if
        if(wdta(l,8)<tmiss) then
          wdta(l,8)=wdta(l,8)*cnvrj
        else
          wdta(l,8)=smiss
        end if
        if(wdta(l,9)<tmiss) then
          wdta(l,9)=wdta(l,9)*cnvrj
        else
          wdta(l,9)=smiss
        end if
        if(wdta(l,10)<tmiss) then
          wdta(l,10)=wdta(l,10)*cnvrj
        else
          wdta(l,10)=smiss
        end if

!        Cloud amount

        if(idta(l,1)<0 .or. idta(l,1)>10) then
          idta(l,1)=imiss
        endif
      end do
      return

      end subroutine rdtmy
!  **********************************************************************

      subroutine rdtry

!  ----------------------------------------------------------------------
!
!      Read 24 hours of data from NOAA Test Reference Year (TRY) format file.
!      This routine is based on the similar BLAST routine by G. Walton
!      and L. Lawrie.
!
!      Modification was made by Cheol Park, NBS
!      Jan. 18, 1985
!      Updated : October 13, 2006 C.P.
!                Code converted from Fortran 77 to Fortran 90.
!
!      variables:
!      tmiss  - Missing value
!      cftmp  - Temperature conversion factor (F -> C)
!      cfpre  - Pressure conversion factor (inHg -> kPa)
!      cfspd  - Speed conversion factor (knots -> m/s)
!      cnvlg  - Radiatin conversion factor (Langleys -> W/m**2)
!
!      The  TRY format data positions and variables are -----
!         position     variable       field description
!         1-5          sttn           Station number
!         6-8          wdta(hour,1)   Dry bulb temperature (F)
!         9-11         wdta(hour,2)   Wet bulb temperature (F)
!         12-14        wdta(hour,11)  Dew point temperature (F)
!         15-17        wdta(hour,6)   Wind direction (degree)
!         18-20        wdta(hour,5)   Wind speed (knots)
!         21-24        wdta(hour,3)   Station pressure(inches of mercury)
!         25           idta(hour,5)   Weather (If =7 rain, =8 snowing)
!         26-27        idta(hour,1)   Total sky cover
!         30           idta(hour,7)   Type of lowest cloud
!         56-59        wdta(hour,12)  Solar radiation (Langleys)
!         70-73        wyr            Year (4 digits)
!         74-75        wmo            Month
!         76-77        wdy            Day
!         78-80                       Blank
!
!  **********************************************************************

      use rdwdf_head
      implicit none
      real          :: tmiss=999.0,smiss=-1000.0
      integer       :: imiss=-1000
      real          :: cftmp=0.5555555556,deltmp=-32.0,cfpre=3.386389,&
                       cfspd=0.5144444,cnvlg=11.62

      wflag=0

!      Read 24 hour data.

      do l=1,24
        read(unit=ifile,fmt=100,iostat=ioerr)&
         sttn,wdta(l,1),wdta(l,2),wdta(l,11),wdta(l,6),&
         wdta(l,5),wdta(l,3),idta(l,5),idta(l,1),&
         idta(l,7),wdta(l,12),wyr,wmo,wdy
        100 format(i5,5f3.0,f4.2,i1,i2,2x,i1,25x,f4.1,10x,i4,2i2,3x)
        if(ioerr>0) then                 ! data reading error
          wflag=-1
          return
        elseif(ioerr<0) then             ! end of file
          wflag=1
          return
        endif
      end do
      call jds(1,wdate,iyear,wmo,wdy)
      if(wdate/=jdnext) return

!      Convert units and fill out undefined data

      do l=1,24
        whr(l)=l-1
        wdta(l,1)=(wdta(l,1)+deltmp)*cftmp
        wdta(l,2)=(wdta(l,2)+deltmp)*cftmp
        wdta(l,11)=(wdta(l,11)+deltmp)*cftmp
        wdta(l,3)=wdta(l,3)*cfpre
        wdta(l,5)=wdta(l,5)*cfspd

!        Solar

        wdta(l,8)=smiss
        wdta(l,9)=smiss
        wdta(l,10)=smiss
        if(wdta(l,12)<tmiss) then
          wdta(l,12)=wdta(l,12)*cnvlg
        else
          wdta(l,12)=smiss
        end if

!        Cloud amount

        if(idta(l,1)<0 .or. idta(l,1)>10) then
          idta(l,1)=imiss
        endif
      end do
      return
      end subroutine rdtry
!  **********************************************************************

      subroutine rdwyec

!  ----------------------------------------------------------------------
!
!      Read 24 hours of data from Weather Year for Energy Calculations
!      (WYEC) format file.  This routine is based on the RDTRY subroutine.
!
!      Modification was made by Cheol Park, NBS
!      Jan. 18, 1985
!      Revised:  July 28, 2000   C.P. NIST
!      Updated : October 13, 2006 C.P.
!                Code converted from Fortran 77 to Fortran 90.
!
!      variables:
!      tmiss  - Missing value
!      cftmp  - Temperature conversion factor (F -> C)
!      cfpre  - Pressure conversion factor (inHg -> kPa)
!      cfspd  - Speed conversion factor (knots -> m/s)
!      cnvlg  - Radiatin conversion factor (Langleys -> W/m**2)
!
!      The  WYEC  format data positions and variables are  -- [Ref. 1]
!         position     variable       field description
!         1-5          sttn           Station number
!         6-8          wdta(hour,1)   Dry bulb temperature (F)
!         9-11         wdta(hour,2)   Wet bulb temperature (F)
!         12-14        wdta(hour,11)  Dew point temperature (F)
!         15-17        wdta(hour,6)   Wind direction (degree)
!         18-20        wdta(hour,5)   Wind speed (knots)
!         21-24        wdta(hour,3)   Station pressure(inches of mercury)
!         25           idta(hour,5)   Weather (If =7 rain, =8 snowing)
!         26-27        idta(hour,1)   Total sky cover
!         30           idta(hour,7)   Type of lowest cloud
!         56-59        wdta(hour,12)  Solar radiation (Langleys)
!         70-73        wyr            Year (4 digits)
!         74-75        wmo            Month
!         76-77        wdy            Day
!         78-79        whr(hour)      Hour (0 TO 23)
!         80           idta(hour,4)   Snow indicator
!
!  **********************************************************************

      use rdwdf_head
      implicit none
      real          :: tmiss=999.0,smiss=-1000.0
      integer       :: imiss=-1000
      real          :: cftmp=0.5555555556,deltmp=-32.0,cfpre=3.386389,&
                       cfspd=0.5144444,cnvlg=11.62

      wflag=0

      do l=1, 24
        read(unit=ifile,fmt=100,iostat=ioerr)&
         sttn,wdta(l,1),wdta(l,2),wdta(l,11),wdta(l,6),&
         wdta(l,5),wdta(l,3),idta(l,5),idta(l,1),&
         idta(l,7),wdta(l,12),wyr,wmo,wdy,whr(l),idta(l,4)
         whr(l) = l - 1
        100 format(i5,5f3.0,f4.2,i1,i2,2x,i1,25x,f4.1,10x,i4,3i2,i1)
        if(ioerr>0) then                 ! data reading error
          wflag=-1
          return
        elseif(ioerr<0) then             ! end of file
          wflag=1
          return
        endif
      end do

      call jds(1,wdate,iyear,wmo,wdy)
!      print *, 'rdwyec2: ',l,wmo,wdy,whr(l),wdate,jdnext
      if(wdate/=jdnext) return

!      Convert units and fill out undefined data

      do l=1,24
        wdta(l,1)=(wdta(l,1)+deltmp)*cftmp
        wdta(l,2)=(wdta(l,2)+deltmp)*cftmp
        wdta(l,11)=(wdta(l,11)+deltmp)*cftmp
        wdta(l,3)=wdta(l,3)*cfpre
        wdta(l,5)=wdta(l,5)*cfspd

!        Solar

        wdta(l,8)=smiss
        wdta(l,9)=smiss
        wdta(l,10)=smiss
        if(wdta(l,12)<tmiss) then
          wdta(l,12)=wdta(l,12)*cnvlg
        else
          wdta(l,12)=smiss
        end if

!        Cloud amount

        if(idta(l,1)<0 .or. idta(l,1)>10) then
          idta(l,1)=imiss
        endif
      end do
      return
      end subroutine rdwyec
!  **********************************************************************

      subroutine rdwyec2

!  ----------------------------------------------------------------------
!
!      Read 24 hours of data from Weather Year for Energy Calculations
!      (WYEC2) ASHRAE file.
!
!      Modifications were made by C. Park.
!      Sept. 11, 2003
!      Updated : October 13, 2006 C.P.
!                Code converted from Fortran 77 to Fortran 90.
!
!      variables:
!      tmiss  - Missing value
!      cnvrj  - Unit conversion factor for radiation
!               Convert hourly solar energy on unit area into
!               solar power per unit area ( kJ/m**2 -> W/m**2)
!
!      The WYEC2 format data positions and variables are ---- [Ref.5]
!         position       variable       field description
!         1-5          sttn           WBAN Station number
!         7-8          wyr            Year     (2 digits)
!         9-10         wmo            Month
!         11-12        wdy            Day
!         13-14        whr(hour)      Hour
!         15-18        wdta(hour,8)   Extraterrestrial radiation(kJ/m**2)
!         19-22        wdta(hour,12)  Global horiz. radiation (kJ/m**2)
!         25-28        wdta(hour,9)   Direct radiation (kJ/m**2)
!         31-34        wdta(hour,10)  Diffuse radiation (kJ/m**2)
!         76           idta(hour,5)   Rain indicator (0 - 8)
!         84-88        wdta(hour,3)   Station pressure (kPa*100)
!         90-93        wdta(hour,1)   Dry bulb temperature (C*10)
!         95-98        wdta(hour,11)  Dew point temperature (C*10)
!         100-102      wdta(hour,6)   Wind direction (degree)
!         104-107      wdta(hour,5)   Wind speed (m/s*10)
!         109-110      idta(hour,1)   Total sky cover [0,10]
!         112-113      idta(hour,7)   Opaque sky cover [0, 10]
!         115          idta(hour,4)   Snow cover
!         116-118                     Blank
!
!  **********************************************************************

      use rdwdf_head
      implicit none
      real          :: tmiss=9999.0,smiss=-1000.0,cnvrj=0.2777777777778
      integer       :: imiss=-1000

      wflag=0

!      Read 24 hour data.

      do l=1,24
        read(unit=ifile,fmt=100,iostat=ioerr)&
         sttn,wyr,wmo,wdy,whr(l),wdta(l,8),wdta(l,12),wdta(l,9),&
         wdta(l,10),idta(l,5),wdta(l,3),wdta(l,1),wdta(l,11),&
         wdta(l,6),wdta(l,5),idta(l,1),idta(l,7),idta(l,4)
        100 format(i5,1x,4i2,2f4.0,2(2x,f4.0),41x,i1,7x,f5.0,&
         2(1x,f4.0),1x,f3.0,1x,f4.0,2(1x,i2),1x,i1,1x)

        if(ioerr>0) then                 ! data reading error
          wflag=-1
          return
        elseif(ioerr<0) then             ! end of file
          wflag=1
          return
        endif
      end do
      call jds(1,wdate,iyear,wmo,wdy)
      if(wdate/=jdnext) return

!      Convert units and fill out undefined data.

      do l=1,24

!        Solar

        if(wdta(l,12)<tmiss) then
          wdta(l,12)=wdta(l,12)*cnvrj
        else
          wdta(l,12)=smiss
        end if
        if(wdta(l,8)<tmiss) then
          wdta(l,8)=wdta(l,8)*cnvrj
        else
          wdta(l,8)=smiss
        end if
        if(wdta(l,9)<tmiss) then
          wdta(l,9)=wdta(l,9)*cnvrj
        else
          wdta(l,9)=smiss
        end if
        if(wdta(l,10)<tmiss) then
          wdta(l,10)=wdta(l,10)*cnvrj
        else
          wdta(l,10)=smiss
        end if

        wdta(l,3) = wdta(l,3)/100           ! pressure
        wdta(l,1) = wdta(l,1)/10            ! dry bulb temperature
        wdta(l,11)= wdta(l,11)/10           ! dew point temperature

!        Cloud amount

        if(idta(l,1)<0 .or. idta(l,1)>10) then
          idta(l,1)=imiss
        endif
      end do
      return
      end subroutine rdwyec2
!  **********************************************************************

      subroutine wrtfil(iflag)

!  ----------------------------------------------------------------------
!
!      wrtfil : Write data on the output file
!
!      OUTPUTS:
!
!      Jan. 14, 1985 C.P.
!      Updated : October 13, 2006 C.P.
!                Code converted from Fortran 77 to Fortran 90.
!
!  **********************************************************************

      use rdwdf_head
      implicit none
      integer         :: iflag

      if(iflag==1) then
        write(unit=ofile,fmt=200)
         200    format(t2,'sttn',t10,'yr',t16,'mo',t20,'day',t26,'hr',&
          t33,'db',t43,'dp',t54,'p',t64,'ws',t71,'cc',t77,'izero',t87,&
          'ibeam',t97,'isky',t106,'ithorz'/t33,'(c)',t43,'(c)',t52,&
          '(kpa)',t63,'(m/s)',t77,'(w/m**2)',t87,'(w/m**2)',t97,&
          '(w/m**2)',t106,'(w/m**2)')
      endif

      do l = 1, 24
        write(unit=ofile,fmt=100)&
         sttn,wyr,wmo,wdy,whr(l),wdta(l,1),wdta(l,11),wdta(l,3),&
         wdta(l,5),idta(l,1),wdta(l,8),wdta(l,9),wdta(l,10),wdta(l,12)
        100 format(1x,i5,4i5,f10.2,3f10.1,i5,4f10.2)
      end do
      return
      end subroutine wrtfil
!  **********************************************************************

      subroutine jds(jflag,jd,year,month,day)

!  ----------------------------------------------------------------------
!
!      Evaluate Julian date. [Ref. 1]
!      input or output:
!          jd    - Julian date
!          year  - Year (4 digits)
!          month - Month
!          day   - Day of the month
!      input:
!          jflag - Control flag
!                  Compute year, month, and day if jflag = 0;
!                  Compute jd from year, month, and day if jflag = 1.
!
!  **********************************************************************

      implicit none
      integer              :: year,month,day,jflag,jd,l,n

      if(jflag/=0) then
        l=(month-14)/12
        jd =day-32075+1461*(year+4800+l)/4 &
         +367*(month-2-l*12)/12-3*((year+4900+l)/100)/4

      elseif(jflag==0) then
        l=jd+68569
        n=4*l/146097
        l=l-(146097*n+3)/4
        year=4000*(l+1)/1461001
        l=l-1461*year/4+31
        month=80*l/2447
        day=l-2447*month/80
        l=month/11
        month=month+2-12*l
        year=100*(n-49)+year+l
      endif
      return
      end subroutine jds
