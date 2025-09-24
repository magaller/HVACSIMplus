! **********************************************************************
!
!     Compute specific heat, desity, vapor pressure, air exchange rate,
!     and humidity ratio at saturation state.
!
!      December 26, 1984 Cheol Park
!      Upadted to Fortran 90  Nov. 16, 2006 Cheol Park, NIST
!
!       cp:       Specific heat of moist air (kJ/kg-K)
!       densit:   Density of moist air (kg/m**3)
!       airex:    Air exchange rate (1/h)
!       sairex:   Standard air exchange rate (1/h)
!                  Living space --  1.5 for leaky building
!                                   1.0 for standard building
!                                   0.5 for modern-type building
!                  Attic space  -- 20.0 for mechanical ventilation
!                                   6.0 for natural ventilation
!                  Crawl space  --  3.0
!       p:        Barometric pressure (kPa)
!       psw:      Saturated vapor pressoure (kPa)
!       pwf:      Vapor pressure (kPa)
!       t:        Dry-bulb temperture (C)
!       tia:      Inside air dry-bulb temperature (C)
!       toa:      Outside air dry-bulb temperature (C)
!       vwind:    Wind speed (m/s)
!       w:        Humidity ratio (-)
!       wsatf:    Saturated humidity ratio (-)
!
!     References:
!     [1] Threlekld, J.L., Thermal Environmental Engineering,
!         2nd ed., Prentice-Hall, 1970.
!     [2] Kusuda, T. and Saitoh, T., "Simplified Heating and Cooling
!         Energy   Analysis Calculations for Residential Applications,"
!         National Bureau of Standards, NBSIR 80-1961, July 1980.
!     [3] Braokaw, R.S., "Calculation for Flue Losses for
!         High-Efficiency  Furnaces and Appliances,"
!         ASHRAE J., Jan. 1979, pp. 49-51.
!     [4] Park,C., Kelly, G.E., and Kao, J.Y., "Economizer Algorithms
!         for Energy Management and Control Systems," Nat'l Bur. of
!         Standards, NBSIR 84-2832, Feb. 1984.
!     [5] Walton, G.N.,"TARP REf. Manual," NBSIR-83-2655,
!         NBS, March 1983.
!
! **********************************************************************

      real function cp(w)

      implicit none
      real             :: w

!     Compute the specific heat of air for the given humidity ratio,
!     ref. [1].

      cp=1.0+1.805*w

      return
      end function cp
! ----------------------------------------------------------------------
      real function densit(p,t,w,pw)

      implicit none
      real             :: p,t,w,pw

!     Compute the density of moist air, ref. [1].

      densit=(p-pw)*(1.0+w)/(0.287055*(t+273.0))

      return
      end function densit
! ----------------------------------------------------------------------
      real function pwf(w,p)

      implicit none
      real             :: w,p

!     Compute vapor pressure of moist air, ref. [1].

      pwf=w*p/(w+0.62198)

      return
      end function pwf
! ----------------------------------------------------------------------
      real function airex(sairex,vwind,toa,tia)

      implicit none
      real             :: sairex,vwind,toa,tia

!     Calculate the air exchange rate for the wind speed and the
!     temperature difference between indoor and outdoor air ref.[2].

      airex=sairex*(0.15+0.013*2.2369*vwind+0.009*abs(toa-tia))/0.695

      return
      end function airex
! ----------------------------------------------------------------------
      real function wsatf(t,p)

      implicit none
      real             :: t,p,psw

!     Compute humidity ratio at saturation state
!     ref.[3,4].

      psw=3.376*exp(15.463-7284.0/(1.8*t+424.0))
      wsatf=0.62198*psw/(p-psw)

      return
      end function wsatf
! **********************************************************************

      real function hiscf(tis,tia,tilt)

! ----------------------------------------------------------------------
!
!     HISCF : Convective heat transfer coefficient of the inner surface
!            to the zone air (W/m**2 K). See p.79, eq. J.2.2 of Ref.[5].
!
!     June 4, 1984 C.P.
!
!       tis:      Inner surface temperature (C)
!       tia:      Indoor air dry-bulb temperature (C)
!       tilt:     Tilt angle of the surface from the horizontal(Degree)
!       vwind:    Reported wind speed (m/s)
!       irofs:    Surface roughness index (-)
!                 = 1 --- Stucco
!                 = 2 --- Brick, rough plaster
!                 = 3 --- Concrete
!                 = 4 --- Clear pine
!                 = 5 --- Smooth plaster
!                 = 6 --- Glass, paint on pine
!       fview:    Surface view angle factor (-)
!       emitis:   Emissivity of the inner surface  (-)
!
! ----------------------------------------------------------------------

      implicit none
      logical                :: upward
      real                   :: tis,tia,tilt
      real                   :: dtr= 0.0174533

      if(tilt>=0.0 .and. tilt<=90.0) then
        if(tis>=tia) then
          upward=.true.
        else
          upward=.false.
        endif
      else
        if(tis>=tia) then
          upward=.false.
        else
          upward=.true.
        endif
      endif

      if(upward) then
        hiscf=9.482*(abs(tia-tis))**0.333/(7.238-abs(cos(tilt*dtr)))
      else
        hiscf=1.810*(abs(tia-tis))**0.333/(1.382+abs(cos(tilt*dtr)))
      endif

      return
      end function hiscf
! ----------------------------------------------------------------------

      real function hosf(vwind,irofs)

!     HOSF : Convective plus radiant heat transfer coefficient of the
!            outer surface (W/m**2 K). See p.71, eq. I.1.2 of Ref.[5].

      implicit none
      real                :: vwind
      integer             :: irofs
      real,dimension(3,6) :: a
      data a/11.58,5.894,0.0,12.49,4.065, 0.028,&
             10.79,4.192,0.0, 8.23,4.000,-0.057,&
             10.22,3.100,0.0, 8.23,3.330,-0.036/

      hosf=a(1,irofs)+a(2,irofs)*vwind+a(3,irofs)*vwind*vwind

      return
      end function hosf
! ----------------------------------------------------------------------

      real function hisrf(tis,fview,emitis)

!     HISRF : Radiant heat transfer coefficient of the inner surface
!             to the zone air (W/m**2 K). See p.69, eq. H.2.3 of Ref.[5].

      implicit none
      real                   :: tis,fview,emitis
      real                   :: sigma= 5.6697e-8

      if(fview<=0.0 .or.emitis<=0.0) then
        hisrf=4.0*sigma*(tis+273.0)**3
      else
        hisrf=4.0*sigma*(tis+273.0)**3/(1.0/fview+(1.0-emitis)/emitis)
      endif

      return
      end function hisrf
! **********************************************************************

      subroutine view(area,f,maxns,ns,maxzn,izn)

! ----------------------------------------------------------------------
!
!     VIEW : Mean Radiant Tempetature Network(MRTN) Method
!            This program was made based on the TARP package
!            written by G.N. Walton.
!
!      July 13, 1984 Cheol Park
!
!       area:     Surface area (m**2)
!       f:        View factor (-)
!       maxns:    Maximum number of surfaces
!       ns:       Number of surfaces
!       maxzn:    Maximum number of zones
!       izn:      Zone ID
!       maxitr:   Maximum number of iterations
!
!     Refernces:
!     [1] CARROL, J.A.,"A Comparision of Radiant Interchange
!                  Algorithms," ASME/SSEA Conference,April 1981.
!     [2] Walton, G.N.,"TARP REf. Manual," NBSIR-83-2655,
!                  NBS, March 1983.
!
! **********************************************************************

      implicit none
      integer,parameter            :: pp=selected_real_kind(15)
      integer                      :: i,j,maxns,ns,maxzn,izn
      integer                      :: maxitr=100
      real,dimension(maxns,maxzn)  :: area,f
      real(kind=pp)                :: sum,suma,x

!     initial guess of f(i,izn) and suma

      sum=0.0
      do i=1 ,ns
        f(i,izn)=1.0
        sum=sum+area(i,izn)
      end do
      suma=sum
      if(suma<=0.0) return

!     Calculate new value of f(i,izn) and use it for computing
!     new suma.

      do j=1,maxitr
        x=0.0
        do i=1 ,ns
          f(i,izn)=1.0/(1.0-area(i,izn)*f(i,izn)/suma)
          x=x+f(i,izn)*area(i,izn)
        end do

!       As a convergence test, compare the new suma, which is
!       represented by x, with the old suma.

        if(abs(1.0-x/suma)<0.1e-5) return

!       if convergence is not satisfied, replace suma with x
!       iterate again.

        suma=x
      end do
      stop  '******   MRT network did not converge ***********'

      end subroutine view
