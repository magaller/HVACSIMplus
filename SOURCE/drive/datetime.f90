! ----------------------------------------------------------------------

      subroutine datetime(is1970,ims,timstamp)

! *** Read system clock and return time as integers and date and time as
! *** a character string. Zero of time is Jan 1, 1970
!
!       Created: Phil Haves, LUT
!
!       Revised: May 8, 1997  Cheol Park, NIST
!
! -------------------------------------------------------------------------

      implicit none
      integer, dimension(8)  ::  dt
      character  (len=10)    ::  date, time, zone
      character  (len=30)    ::  timstamp
      integer                ::  iyr,imon,iday,ihr,imin,isec,i100th,&
                                 ims,is1970,nday,iy

      integer, dimension(12) ::  ndaynly=(/0,31,59,90,120,151,&
                                         181,212,243,273,304,334/)

      integer, dimension(12) ::  ndayly =(/0,31,60,91,121,152,&
                                          182,213,244,274,305,335/)

      namelist /name_time_stamp/ TIMSTAMP

! *** Read system clock/calendar and calculate integer seconds (IS) and
! *** milliseconds (IMS) since Jan 1, 1970

      call date_and_time(date, time, zone, dt)

      iyr  = dt(1)
      imon = dt(2)
      iday = dt(3)
      ihr  = dt(5)
      imin = dt(6)
      isec = dt(7)
      ims  = dt(8)
      i100th = nint(ims/10.0)
      timstamp = ' date: '//date//' time: '//time

!   Calculate number of days from Jan 1 1970 to Dec 31 of previous year
      nday = 0
      do iy=1970,(iyr-1)
!   Test for leap year and calculate number of days since start of year

  	if (mod(iy,4)==0 .and. iy/=2000) then

!   Leap year
           nday = nday+ndayly(12)
  	else
!   Not a leap year
           nday = nday+ndaynly(12)
  	endif
      end do

!   Add in completed days in current year
!   Test for leap year and calculate number of days since start of year
      if (mod(iyr,4)==0 .and. iyr/=2000) then
!   Leap year
         nday = nday+ndayly(imon)
      else
!   Not a leap year
         nday = nday+ndaynly(imon)
      endif
!   Number of seconds since start of 1970
      is1970=isec+60*(imin+60*(ihr+24*(iday+nday)))
!          print *, ' is1970 =  ', is1970

!   Milliseconds - rounded to nearest 100th of a second
      ims=10*i100th

!   Build a string of the form 'date: yymmdd  time: hh:mm:ss'

!          read(timstamp,900) iyr,imon,iday,ihr,imin,isec,i100th

!          print name_time_stamp
!          print 900,  iyr,imon,iday,ihr,imin,isec,i100th

 900    format('   date: ',i5,2i3,' time: ',3(i2,':'),'.',i2)
      end  subroutine datetime

