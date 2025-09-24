!=======================================================================
!
!                  STEAM AND LIQUID WATER PROPERTIES
!
!=======================================================================
!
!  References:
!
!    1)  Irvine, T.F.Jr., and Liley, P.E., "Steam and Gas Tables
!        with Computer Equations," Academic Press, Inc., 1984.
!    2)  Van Wylen, G.J., and Sonntag, R.E., Fundamentals of Classical
!        Thermodynamics (SI Version).  New York: John Wiley and Sons.
!    3)  Chapman, A.J.  Heat Transfer (3rd Edition).  New York:
!        Macmillan  Publishing Co. (1974).
!    4)  Chi, J.
!    5)  Karlekar, B.V., and Desmond, R.M., Engineering Heat Transfer.
!        St. Paul: West Publishing Co., 1977.
!    6)  CRC Handbook of Chemistry and Physics, 61st Ed. (1980-1981).
!
!  This file contains the following functions:
!
!  Superheated Steam:
!       vs(pkpa,tc)     eq. from [1]
!       hs(pkpa,tc)     eq. from [1]
!       ss(pkpa,tc)     eq. from [1]
!       cps(tc)         eq. from [2]
!       tpss(pkpa,s)    derived using functions SS and CPS
!       vissph(tc)      curve fit to data in [3]
!       steamk(t)       curve fit to data in [3]
!       cvs(v,t)        adapted from subroutine in [4]
!
!  Saturated Steam:
!       tsats(pkpa)     eq. from [1]
!       psats(tc)       eq. from [1]
!       vsats(pkpa,tc)  eq. from [1]
!       hsats(tc)       eq. from [1]
!       ssats(tc)       eq. from [1]
!       vissv(p)        curve fit to data in [3]
!
!  Two Phase:
!       hfg(tc)         eq. from [1]
!
!  Saturated Liquid Water:
!       vsatw(tc)       eq. from [1]
!       hsatw(tc)       eq. from [1]
!       ssatw(tc)       eq. from [1]
!       wmu(t)          curve fit to data in [5] for T above 100 C.
!       wcp(t)          curve fit to data in [5] for T above 100 C.
!       wk(t)           curve fit to data in [5] for T above 100 C.
!
!  Liquid Water at 1 Atmosphere:
!       wmu(t)          eq. from [6]
!       wcp(t)          eq. from [6]
!       wk(t)           curve fit to data in [6]
!       wrho(t)         eq. from [6]
!
!***********************************************************************

      real function tsats(pkpa)

! ----------------------------------------------------------------------
!
!  Saturation temp. of steam (C) as a function of pressure (kPa)
!
!***********************************************************************

      implicit none
      real :: pkpa,p
      real :: a1=42.6776,b1=-3892.70,c1=-9.48654,tconv=-273.15,&
              a2=-387.592,b2=-12587.5,c2=-15.2578,pconv=0.001

      p=pkpa*pconv

      if(p < 12.33) then
        tsats=tconv+a1+b1/(log(p)+c1)
      else
        tsats=tconv+a2+b2/(log(p)+c2)
      endif

      return
      end function tsats
!***********************************************************************

      real function psats(tc)

! ----------------------------------------------------------------------
!
!  Saturation pressure of steam (kPa) as a function of temperature (C)
!
!***********************************************************************

      implicit none
      real :: tc,t,plog
      real :: a0=10.4592,a1=-0.404897e-2,a2=-0.417520e-4,a3=0.368510e-6,&
              a4=-0.101520e-8,a5=0.865310e-12,a6=0.903668e-15,&
              a7=-0.199690e-17,a8=0.779287e-21,a9=0.191482e-24,&
              a10=-3968.06,a11=39.5735,tconv=273.15,pconv=1000.

      t=tc+tconv
      plog=a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*a9&
       )))))))) + a10/(t-a11)
      psats=pconv*exp(plog)

      return
      end function psats
!***********************************************************************

      real function vsats(pkpa,tc)

! ----------------------------------------------------------------------
!
!  Sat. specific volume of steam (m3/kg) given sat. t (C) and p (kPa)
!
!***********************************************************************

      implicit none
      real :: pkpa,tc,tr,y
      real :: a=1.0,b=1.6351057,c=52.584599,d=-44.694653,e1=-8.9751114,&
              e2=-0.43845530,e3=-19.179576,e4=36.765319,e5=-19.462437,&
              tcr=647.3,pcr=22.089,vcr=3.155e-3,tcnv=273.15,pcnv=0.001

      tr=(tcr-tc-tcnv)/tcr
      y=a+b*tr**(1./3.)+c*tr**(5./6.)+d*tr**0.875
      y=y+tr*(e1+tr*(e2+tr*(e3+tr*(e4+tr*e5))))
      vsats=y*pcr*vcr/(pkpa*pcnv)

      return
      end function vsats
!***********************************************************************

      real function vsatw(tc)

! ----------------------------------------------------------------------
!
!  Sat. specific volume of water (m3/kg) given sat. t (C)
!
!***********************************************************************

      implicit none
      real :: tc,tr,y
      real :: a=1.0,b=-1.9153882,c=12.015186,d=-7.8464025,e1=-3.8886414,&
              e2=2.0582238,e3=-2.0829991,e4=0.82180004,e5=0.47549742,&
              tcr=647.3,vcr=3.155e-3,tcnv=273.15

      tr=(tcr-tc-tcnv)/tcr
      y=a+b*tr**(1./3.)+c*tr**(5./6.)+d*tr**0.875
      y=y+tr*(e1+tr*(e2+tr*(e3+tr*(e4+tr*e5))))
      vsatw=y*vcr

      return
      end function vsatw
!***********************************************************************

      real function hsatw(tc)

! ----------------------------------------------------------------------
!
!  Sat. enthalpy of liquid water (kJ/kg) given sat. t (C)
!
!***********************************************************************

      implicit none
      real :: tc,tk,tr,y
      real :: e11=624.698837,e21=-2343.85369,e31=-9508.12101,hfcr=2099.3,&
              e41=71628.7928,e51=-163535.221,e61=166531.093,tcnv=273.15,&
              e71=-64785.4585,a2=0.8839230108,e12=-2.67172935,&
              e22=6.22640035,e32=-13.1789573,e42=-1.91322436,&
              e52=68.7937653,e62=-124.819906,e72=72.1435404,a3=1.0,&
              b3=-0.441057805,c3=-5.52255517,d3=6.43994847,&
              e13=-1.64578795,e23=-1.30574143,tcr=647.3
!
      tk=tc+tcnv
      tr=(tcr-tk)/tcr

      if(tk<300.0)then
        y=tr*(e11+tr*(e21+tr*(e31+tr*(e41+tr*(e51+tr*(e61+tr*e71))))))
      else if(tk<600.0)then
        y=tr*(e12+tr*(e22+tr*(e32+tr*(e42+tr*(e52+tr*(e62+tr*e72))))))
        y=y+a2
      else
       y=a3+b3*tr**(1./3.)+c3*tr**(5./6.)+d3*tr**0.875+tr*(e13+tr*e23)
      endif

      hsatw=y*hfcr

      return
      end function hsatw
!***********************************************************************

      real function hfg(tc)

! ----------------------------------------------------------------------
!
!  Latent heat of vaporization of water (kJ/kg) given sat. t (C)
!
!***********************************************************************

      implicit none
      real :: tc,tr,y
      real :: e1=-3.87446,e2=2.94553,e3=-8.06395,e4=11.5633,e5=-6.02884,&
              b=0.779221,c=4.62668,d=-1.07931,hfgtp=2500.9,tcr=647.3,&
              tcnv=273.15
!
      if(tc<0.) then
        tc=0.
      endif

      if(tr<0.) then
        hfg=0.
        return
      else
        tr=(tcr-tc-tcnv)/tcr
      endif

      y=b*tr**(1./3.)+c*tr**(5./6.)+d*tr**0.875
      y=y+tr*(e1+tr*(e2+tr*(e3+tr*(e4+tr*e5))))
      hfg=y*hfgtp

      return
      end function hfg
!***********************************************************************

      real function hsats(tc)

! ----------------------------------------------------------------------
!
!  Enthalpy of saturated steam (kJ/kg) given sat. t (C)
!
!***********************************************************************

      implicit none
      real :: tc,tr,y
      real :: e1=-4.81351884,e2=2.69411792,e3=-7.39064542,e4=10.4961689,&
              e5=-5.46840036,b=0.457874342,c=5.08441288,d=-1.48513244,&
              a=1.0,tcr=647.3,hcr=2099.3,tcnv=273.15

      tr=(tcr-tc-tcnv)/tcr
      y=a+b*tr**(1./3.)+c*tr**(5./6.)+d*tr**0.875
      y=y+tr*(e1+tr*(e2+tr*(e3+tr*(e4+tr*e5))))
      hsats=y*hcr

      return
      end function hsats

!***********************************************************************

      real function ssatw(tc)

! ----------------------------------------------------------------------
!
!  Sat. entropy of liquid water [kJ/(kg K)] given sat. t (C)
!
!***********************************************************************

      implicit none
      real :: tc,tk,tr,y
      real :: e11=-1836.92956,e21=14706.6352,e31=-43146.6046,scr=4.4289,&
              e41=48606.6733,e51=7997.5096,e61=-58333.9887,tcnv=273.15,&
              e71=33140.0718,a2=0.912762917,e12=-1.75702956,tcr=647.3,&
              e22=1.68754095,e32=5.82215341,e42=-63.3354786,&
              e52=188.076546,e62=-252.344531,e72=128.058531,a3=1.0,&
              b3=-0.324817650,c3=-2.990556709,d3=3.2341900,&
              e13=-0.678067859,e23=-1.91910364

      tk=tc+tcnv
      tr=(tcr-tk)/tcr

      if(tk<300.0)then
        y=tr*(e11+tr*(e21+tr*(e31+tr*(e41+tr*(e51+tr*(e61+tr*e71))))))
      else if(tk<600.0)then
        y=tr*(e12+tr*(e22+tr*(e32+tr*(e42+tr*(e52+tr*(e62+tr*e72))))))
        y=y+a2
      else
        y=a3+b3*tr**(1./3.)+c3*tr**(5./6.)+d3*tr**0.875+tr*(e13+tr*e23)
      endif

      ssatw=y*scr

      return
      end function ssatw

!***********************************************************************

      real function ssats(tc)

! ----------------------------------------------------------------------
!
!  Entropy of saturated steam [kJ/(kg K)] given sat. t (C)
!
!***********************************************************************

      implicit none
      real :: tc,tr,y
      real :: e1=-4.34839,e2=1.34672,e3=1.75261,e4=-6.22295,e5=9.99004,&
              a=1.0,b=0.377391,c=-2.78368,d=6.93135,tcr=647.3,scr=4.4289,&
              tcnv=273.15

      tr=(tcr-tc-tcnv)/tcr
      y=a+b*tr**(1./3.)+c*tr**(5./6.)+d*tr**0.875
      y=y+tr*(e1+tr*(e2+tr*(e3+tr*(e4+tr*e5))))
      ssats=y*scr

      return
      end function ssats

!***********************************************************************

      real function vs(pkpa,tc)

! ----------------------------------------------------------------------
!
!  Specific volume of superheated steam (m3/kg) given p (kPa) and t (C)
!
!***********************************************************************


      implicit none
      real :: pkpa,tc,p,t,ts
      real :: r=4.61631e-4,b1=5.27993e-2,b2=3.75928e-3,b3=0.022,em=40.0,&
              a0=-3.741378,a1=-4.7838281e-3,a2=1.5923434e-5,tcnv=273.15,&
              a3=10.0,c1=42.6776,c2=-3892.70,c3=-9.48654,pcnv=0.001,&
              c4=-387.592,c5=-12587.5,c6=-15.2578

      p=pkpa*pcnv
      t=tc+tcnv

      if(p >= 12.33) then
        ts=c4+c5/(log(p)+c6)
      else
        ts=c1+c2/(log(p)+c3)
      endif

      vs=r*t/p-b1*exp(-b2*t)+(b3-exp(a0+ts*(a1+ts*a2)))/(a3*p)&
         *exp((ts-t)/em)

      return
      end function vs
!***********************************************************************

      real function hs(pkpa,tc)

! ----------------------------------------------------------------------
!
!  Enthalpy of superheated steam (kJ/kg) given p (kPa) and t (C)
!
!***********************************************************************

      implicit none
      real :: pkpa,tc,p,t,ts,a0,a1,a2,a3
      real :: b11=2041.21,b12=-40.4002,b13=-0.48095,b21=1.610693,&
              b22=5.472051e-2,b23=7.517537e-4,b31=3.383117e-4,&
              b32=-1.975736e-5,b33=-2.87409e-7,b41=1707.82,b42=-16.99419,&
              b43=6.2746295e-2,b44=-1.0284259e-4,b45=6.4561298e-8,em=45.0,&
              c1=42.6776,c2=-3892.70,c3=-9.48654,pcnv=0.001,&
              c4=-387.592,c5=-12587.5,c6=-15.2578,tcnv=273.15

      p=pkpa*pcnv
      t=tc+tcnv

      if(p >= 12.33) then
        ts=c4+c5/(log(p)+c6)
      else
        ts=c1+c2/(log(p)+c3)
      endif

      a0=b11+p*(b12+p*b13)
      a1=b21+p*(b22+p*b23)
      a2=b31+p*(b32+p*b33)
      a3=b41+ts*(b42+ts*(b43+ts*(b44+ts*b45)))

      hs=a0+t*(a1+t*a2)-a3*exp((ts-t)/em)

      return
      end function hs
!***********************************************************************

      real function ss(pkpa,tc)

! ----------------------------------------------------------------------
!
!  Entropy of superheated steam [kJ/(kg K)] given p (kPa) and t (C)
!
!***********************************************************************

      implicit none
      real :: pkpa,tc,p,t,ts
      real :: a0=4.6162961,a1=1.039008e-2,a2=-9.873085e-6,a3=5.43411e-9,&
              a4=-1.170465e-12,b1=-0.4650306,b2=0.001,b3=10.0,c0=1.777804,&
              c1=-1.802468e-2,c2=6.854459e-5,c3=-1.184434e-7,em=85.0,&
              c4=8.142201e-11,e1=42.6776,e2=-3892.70,e3=-9.48654,&
              e4=-387.592,e5=-12587.5,e6=-15.2578,tcnv=273.15

      p=pkpa*b2
      t=tc+tcnv

      if(p >= 12.33) then
        ts=e4+e5/(log(p)+e6)
      else
        ts=e1+e2/(log(p)+e3)
      endif

      ss=a0+t*(a1+t*(a2+t*(a3+t*a4)))+b1*log(b2+p*b3)-exp((ts-t)/em)*&
         (c0+ts*(c1+ts*(c2+ts*(c3+ts*c4))))

      return
      end function ss

!***********************************************************************

      real function tpss(p,s)

! ----------------------------------------------------------------------
!
!  Temperature (C) of steam,  given p (kPa) and s [kJ/(kg K)]
!
!***********************************************************************

      implicit none
      integer :: i
      real    :: p,s,to,so,ta,sa,t,ssats,cps,ss
      real    :: e1=42.6776,e2=-3892.70,e3=-9.48654,pcnv=0.001,&
                 e4=-387.592,e5=-12587.5,e6=-15.2578,tabs=273.15

!  Compare input entropy with saturation value

      if(p >= 12330.0) then
        to=e4-tabs+e5/(log(p*pcnv)+e6)
      else
        to=e1-tabs+e2/(log(p*pcnv)+e3)
      endif

      so=ssats(to)

      if(so>=s)then
        tpss=to
        return
      endif

!  Initial guess ta is based on assumption of constant specific heat.
!  Subsequent approximations made by interpolation.

      ta=(to+tabs)*(1.+(s-so)/cps(to))-tabs
      sa=ss(p,ta)

      do i=1,10
        t=ta+(to-ta)*(s-sa)/(so-sa)
        if(abs(t-ta)<0.05) goto 900
        to=ta
        so=sa
        ta=t
        sa=ss(p,ta)
      enddo

      write(1,100)
 100  format(' Warning: function tpss fails to converge')
 900  tpss=t

      return
      end function tpss
!***********************************************************************

      real function cps(t)

! ----------------------------------------------------------------------
!
!       Determine specific heat of steam, cp, (kJ/kg/K) given temp. (C)
!
!       Specific heat equation from "Fundamentals of Classical
!       Thermodynamics-si version" by van Wylen and Sonntag
!       Table a.9, pg. 683.
!
!       Valid for t between 300-3500 K   max error = .43%
!
!***********************************************************************

      implicit none
      real :: t,t1,tk
      real :: c1=143.05,c2=-183.54,c3=82.751,c4=-3.6989,e1=0.25,e2=0.5

      tk=t+273.15
      t1=tk/100.
      cps=(c1+c2*t1**e1+c3*t1**e2+c4*t1)/18.015

      if(tk<300.0.or.tk>3500.0)then
        write(1,100)
 100    format(' ',' warning: function cps: t out of range')
      endif

      return
      end function cps
!***********************************************************************

      real function cvs(v,t)

! ----------------------------------------------------------------------
!
!       Calculate cv (kJ/kg/K) given v (m3/kg) and t (C)
!
!***********************************************************************

      implicit none
      real :: v,t,tr,ve
      real :: tc=1165.11,tfr=459.67,b1=0.0063101
      real,dimension(7)  :: a=(/0.99204818,-33.137211,416.29663,0.185053,&
                                5.475,-2590.5815,113.95968/)

      tr=9./5.*t+32.+tfr
      ve=(v-b1)/0.062428
      cvs=a(1)+a(2)/sqrt(tr)+a(3)/tr-a(4)*a(5)**2*tr/tc**2*&
          exp(-a(5)*tr/tc)*(a(6)/ve+a(7)/ve**2)
      cvs=cvs*4.1868

      return
      end function cvs

!***********************************************************************

      real function vissv(p)

! ----------------------------------------------------------------------
!
!  Calculates the dynamic viscosity (kg/m-s) of saturated
!  vapor given the pressure (kPa).  Correlation obtained from
!  a curve fit of data from 'Heat Transfer' by Alan J. Chapman, 1974.
!
!***********************************************************************

      implicit none
      real :: p,psi
      real :: c1=0.0314,c2=2.9675e-05,c3=-1.60583e-08,c4=3.768986e-12

!     Covert pressure from kPa to psi

      psi=p/6.894757
      vissv=c1+c2*psi+c3*psi**2+c4*psi**3

!     Convert viscosity from lbm/ft-hr to kg/m-s

      vissv=vissv*4.1338e-04

      return
      end function vissv
!***********************************************************************

      real function vissph(t)

! ----------------------------------------------------------------------
!
!   Calculates the dynamic viscosity (kg/m-s) of superheated
!   steam given the temperature (C).  The correlation is obtained
!   from a curve fit at atmospheric pressure from 'Heat Transfer'
!   by Alan J. Chapman, 1974. (Note: There is little  variation in
!   viscosity at higher pressures.)
!
!***********************************************************************
!
      implicit none
      real :: t,tf
      real :: c1=0.0183161,c2=5.7067e-05,c3=-1.42253e-08,c4=7.241555e-12

!     Convert temperature from C to F

      tf=t*1.8+32.
      vissph=c1+c2*tf+c3*tf**2+c4*tf**3

!     Convert viscosity from lbm/ft-hr to kg/m-s

      vissph=vissph*4.1338e-04

      return
      end function vissph
!***********************************************************************

      real function steamk(t)

! ----------------------------------------------------------------------
!
!     Calculates thermal conductivity of superheated steam (kW/m-C)
!     given the temperature (C).  Curve fit from data in 'Heat Transfer'
!     by Alan J. Chapman, 1974.
!
!***********************************************************************
!
      implicit none
      real :: t,tf
      real :: c1=0.824272,c2=0.00254627,c3=9.848539e-08

!     Convert temperature from C to F

      tf=t*1.8+32.
      steamk=(c1+c2*tf+c3*tf**2)*.01

!     Convert k from btu/hr-ft-F to kW/m-C

      steamk=steamk*0.0017308

      return
      end function steamk
!***********************************************************************

      real function wrho(tw)

! ----------------------------------------------------------------------
!
!  Density eq. for water at 1 atm., from CRC Handbook of Chem. & Phys.,
!   61st edition (1980-1981), p. f-6.  Density (kg/m3) given temp. (C).
!
!***********************************************************************
!
      implicit none
      real :: tw
      real :: ar0=999.83952,ar1=16.945176,ar2=-0.0079870401,ar6=0.01687985,&
              ar3=-46.170461e-06,ar4=105.56302e-09,ar5=-280.54253e-12

      wrho=(ar0+tw*(ar1+tw*(ar2+tw*(ar3+tw*(ar4+tw*ar5)))))/&
           (1.+ar6*tw)

      return
      end function wrho

!***********************************************************************

      real function wmu(tw)

! ----------------------------------------------------------------------
!
!  Viscosity equations for water at 1 atm., from CRC Handbook (op.cit.),
!    page f-51.  wmu in kg/meter-sec; for centipoise, multiply by 1000.
!    for temps > 100 C, fit to data from Karlekar & Desmond (saturated).
!
!***********************************************************************

      implicit none
      real :: tw
      real :: am0=-3.30233,am1=1301.0,am2=998.333,am3=8.1855,am4=0.00585,&
              am5=1.002,am6=-1.3272,am7=-0.001053,am8=105.0,&
              a10=0.68714,a11=-0.0059231,a12=2.1249e-05,a13=-2.69575e-08

      if(tw<20.0) then
        wmu=10.**(am0+am1/(am2+(tw-20.)*(am3+am4*(tw-20.))))*100.
      elseif(tw>100.0) then
        wmu=a10+tw*(a11+tw*(a12+tw*a13))
      else
        wmu=am5*10.**((tw-20.)*(am6+(tw-20.)*am7)/(tw+am8))
      endif

      wmu=0.001*wmu

      return
      end function wmu
!***********************************************************************

      real function wk(tw)

! ----------------------------------------------------------------------
!
!  Thermal conductivity equation from linear least-squares fit to data
!   in CRC Handbook (op.cit.), page e-11; temps. from 270 K to 620 K.
!    temperature in celsius, wk in [kW/(m K)].  Values at one atmosphere
!    for t from 0 to 100 C, at saturation for t above 100.
!
!***********************************************************************

      implicit none
      real :: tw
      real :: ak0=0.560101,ak1=0.00211703,ak2=-1.05172e-05,&
              ak3=1.497323e-08,ak4=-1.48553e-11

      wk=0.001*(ak0+tw*(ak1+tw*(ak2+tw*(ak3+tw*ak4))))

      return
      end function wk
!**********************************************************************

      real function wcp(tw)

! ----------------------------------------------------------------------
!
!  Specific heat of water at 1 atmosphere, 0 to 100 C.  equation from
!    linear least-squares regression of data from CRC Handbook (op.cit.)
!    page d-174; in J/g-C (or kJ/kg-C).
!    For temps > 100, fit to data from Karlekar & Desmond (saturated).
!
!***********************************************************************

      implicit none
      real :: tw
      real :: acp0=4.21534,acp1=-0.00287819,acp2=7.4729e-05,&
              acp3=-7.79624e-07,acp4=3.220424e-09,acp5=2.9735,&
              acp6=0.023049,acp7=-0.00013953,acp8=3.092474e-07

      if(tw>100.0) then
        wcp=acp5+tw*(acp6+tw*(acp7+tw*acp8))
      else
        wcp=acp0+tw*(acp1+tw*(acp2+tw*(acp3+tw*acp4)))
      endif

      return
      end function wcp

