!=======================================================================
!
!                  AIR PROPERITES
!
!=======================================================================
!
!  Ref:  Irvine, T.F.Jr., and Liley, P.E., "Steam and Gas Tables
!        with Computer Equations," Academic Press, Inc., 1984.
!
!***********************************************************************

      subroutine cpcva(tca,cpa,cva,gamma,sonic)

! ----------------------------------------------------------------------
!
!     This subroutine takes Celsius air temperature, tca, and computes:
!
!     cpa:   specific heat of air at constant pressure [kJ/(kg*K)]
!     cva:   specific heat of air at constant volume [kJ/(kg*K)]
!     gamma: the ratio cpa/cva [dimensionless]
!     sonic: speed of sound in air (m/s)
!
!***********************************************************************

      implicit none
      real    :: tca,cpa,cva,gamma,sonic,t
      real    :: a0=1.03409,a1=-0.284887e-3,a2=0.7816818e-6,&
                 a3=-0.4970786e-9,a4=0.1077024e-12,tconv=273.15,&
                 r=0.287040

      t=tca+tconv
      cpa=a0+t*(a1+t*(a2+t*(a3+t*a4)))
      cva=cpa-r
      gamma=cpa/cva
      sonic=sqrt(gamma*r*t)

      return
      end subroutine cpcva
!***********************************************************************

      real function ha(tc)

! ----------------------------------------------------------------------
!
!     Enthalpy of air [kJ/kg] as a function of temperature [C]
!
!***********************************************************************

      implicit none
      real    :: tc,t
      real    :: a0=12.074,a1=0.924502,a2=0.115984e-3,a3=-0.563568e-8,&
                 tconv=273.15

      t=tc+tconv
      ha=a0+t*(a1+t*(a2+t*a3))

!     Note that internal energy,u, is equal to ha-r*t, where r=0.287040

      return
      end function ha
!***********************************************************************

      real function phia(tca)

! ----------------------------------------------------------------------
!
!     Entropy function, phi(t) = [ cpa(t)/t dt ] integrated from t0 to t
!      where t0 is a reference temperature at which entropy = 0.
!      phia has same units as entropy [kJ/(kg*K)]
!     Note that s2 - s1 = phia(tc2) - phia(tc1) - r*log(p2/p1)
!      where r = 0.287040 [kJ/(kg-K)]
!     Other useful relationships:
!      isentropic pressure function, log(p/p0) = phia(tca) / r
!      isentropic volume function,   log(v/v0) = log(r*t) - log(p/p0)
!
!***********************************************************************

      implicit none
      real    :: tca,t
      real    :: a0=1.386989,a1=0.184930e-3,a2=0.95,tconv=273.15

      t=tca+tconv
      phia=a0+a1*t+a2*log(t)

      return
      end function phia
!***********************************************************************

      real function tphia(phia)

! ----------------------------------------------------------------------
!
!     Temperature [C] of air as a function of the entropy function, phia
!
!***********************************************************************

      implicit none
      real    :: phia,pr
      real    :: a0=-8800.92,a1=1269.74,a2=-61.9391,a3=1.03530
      real    :: tconv=273.15,r=0.287040

      pr=phia/r
      tphia=a0+pr*(a1+pr*(a2+pr*a3))-tconv

      return
      end function tphia
!***********************************************************************

      real function visca(tc)

! ----------------------------------------------------------------------
!
!     Dynamic viscosity [(N*s)/(m*m)] of air given Celsius temperature
!
!***********************************************************************


      implicit none
      real    :: tc,t
      real    :: a0=-0.98601,a1=9.080125e-2,a2=-1.17635575e-4,&
                 a3=1.2349703e-7,a4=-5.7971299e-11,&
                 b0=4.8856745,b1=5.43232e-2,b2=-2.4261775e-5,&
                 b3=7.9306e-9,b4=-1.10398e-12,tconv=273.15

      t=tc+tconv

      if(t >= 600.0) then
        visca=b0+t*(b1+t*(b2+t*(b3+t*b4)))
      else
        visca=a0+t*(a1+t*(a2+t*(a3+t*a4)))
      end if

      visca=visca*1.0e-6

      return
      end function visca
!***********************************************************************

      real function aka(tc)

! ----------------------------------------------------------------------
!
!     Thermal conductivity of air [kW/(m*K)] given Celsius temperature
!
!***********************************************************************

      implicit none
      real    :: tc,t
      real    :: c0=-2.276501e-3,c1=1.2598485e-4,&
                 c2=-1.4815235e-7,c3=1.73550646e-10,c4=-1.066657e-13,&
                 c5=2.47663035e-17,tconv=273.15

      t=tc+tconv
      aka=0.001*(c0+t*(c1+t*(c2+t*(c3+t*(c4+t*c5)))))

      return
      end function aka
