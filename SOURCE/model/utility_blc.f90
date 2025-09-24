! *********************************************************************

      real function cpf(tc1,tc2)

! ---------------------------------------------------------------------
!
!     Combustion product specific heat
!
!     INPUTS:
!        tc1    Air temperature (C)
!        tc2    Combustion product temperature (C)
!        pco2   Number of moles of CO2
!        ph2o   Number of moles of H2O
!        po2    Number of moles of O2
!        pso2   Number of moles of SO2
!        pn2    Number of moles of N2
!
!     OUTPUT:
!        cpf    Specific heat of combustion product (kJ/(kg-C))
!
! *********************************************************************

      use blc_head
      implicit none
      real,dimension(6,2)   :: c
      real,dimension(30)    :: a,b
      real                  :: tc2,tc1,t1,t2,tm
      integer               :: j,i,m,jj

      data a / 4.2497678,  -6.912652e-3,  3.1602134e-5,  -2.9715432e-8,&
            9.510358e-12,        2.1701,  1.0378115e-2,  -1.0733938e-5,&
            6.3459175e-9,-1.6280701e-12,     4.1565016,  -1.7244334e-3,&
            5.6982316e-6, -4.5930044e-9, 1.4233654e-12,      3.6916148,&
           -1.3332552e-3,    2.65031e-6, -9.768834e-10, -9.9772234e-14,&
               3.7189946, -2.5167288e-3,  8.5837353e-6,  -8.2998716e-9,&
            2.708218e-12,     3.2257132,  5.6551207e-3,  -2.4970208e-7,&
           -4.2206766e-9, 2.1392733e-12/
      data b / 1.1795744,  1.0950594e-2,  -4.062213e-6,  7.1370281e-10,&
          -4.7490353e-14,     4.4129266,  3.1922896e-3,   -1.297823e-6,&
           2.4147446e-10, 1.6742986e-14,     2.6707532,   3.0317115e-3,&
            -8.535157e-7, 1.1790853e-10,-6.1973568e-15,      2.8545761,&
            1.5976316e-3, -6.2566254e-7, 1.1315849e-10,  -7.689707e-15,&
               3.5976129,  7.8145603e-4,  -2.238667e-7,  4.2490159e-11,&
          -3.3460204e-15,     5.1982451,  2.0595095e-3,  -8.6254450e-7,&
           1.6636523e-10,-1.1847837e-14/

!     Set t1 and t2 in degree K

      t1=273.15+tc1
      t2=273.15+tc2
      tm=1000.
      c=0.

      if(t1<1000.) then
        if(t2<1000.) then
          do i=1,6
            m=5*(i-1)
            do j=1,5
              jj=m+j
              c(i,1)=c(i,1)+a(jj)*(t2**j-t1**j)/j
            enddo
            c(i,1)=c(i,1)*8.32066/(t2-t1)
          enddo
          cpf=(pco2*c(2,1)+ph2o*c(3,1)+pn2*c(4,1)&
              +po2*c(5,1)+pso2*c(6,1))/(pco2*44.&
              +ph2o*18.+pn2*28.+po2*32.+pso2*64.)
        else
          do i=1,6
            m=5*(i-1)
            do j=1,5
              jj=m+j
              c(i,1)=c(i,1)+a(jj)*(tm**j-t1**j)/j
              c(i,2)=c(i,2)+b(jj)*(t2**j-tm**j)/j
            enddo
            c(i,1)=(c(i,1)+c(i,2))*8.32066/(t2-t1)
          enddo
          cpf=(pco2*c(2,1)+ph2o*c(3,1)+pn2*c(4,1)&
              +po2*c(5,1)+pso2*c(6,1))/(pco2*44.&
              +ph2o*18.+pn2*28.+po2*32.+pso2*64.)
        endif
      else
        if(t2>1000.) then
          do i=1,6
            m=5*(i-1)
            do j=1,5
              jj=m+j
              c(i,2)=c(i,2)+b(jj)*(t2**j-t1**j)/j
            enddo
            c(i,2)=c(i,2)*8.32066/(t2-t1)
          enddo
          cpf=(pco2*c(2,2)+ph2o*c(3,2)+pn2*c(4,2)&
              +po2*c(5,2)+pso2*c(6,2))/(pco2*44.&
              +ph2o*18.+pn2*28.+po2*32.+pso2*64.)
        else
          do i=1,6
            m=5*(i-1)
            do j=1,5
              jj=m+j
              c(i,1)=c(i,1)+a(jj)*(tm**j-t2**j)/j
              c(i,2)=c(i,2)+b(jj)*(t1**j-tm**j)/j
            enddo
            c(i,1)=(c(i,1)+c(i,2))*8.32066/(t1-t2)
          enddo
          cpf=(pco2*c(2,1)+ph2o*c(3,1)+pn2*c(4,1)&
              +po2*c(5,1)+pso2*c(6,1))/(pco2*44.&
              +ph2o*18.+pn2*28.+po2*32.+pso2*64.)
        endif
      endif
      return
      end function cpf

! *********************************************************************

      real function gef(ppco2,pph2o,xl,tg)

! ---------------------------------------------------------------------
!
!     Gas emissivity
!
!     INPUTS:
!        ppco2   Partial pressure of CO2 (atm)
!        pph2o   Partial pressure of H2O (atm)
!        xl      Length of combustion (ft)
!        tg      Radiating gas temperature (C)
!
!     OUTPUT:
!        gef     Gas emissivity
!
! *********************************************************************

      implicit none
      real          :: tg,xl,pph2o,ppco2,xx,c1,c2,trg,gex

!     Select gas emissivity constants C1 & C2

      xx=3.2808*xl*(ppco2+pph2o)
      if(xx>0.2) then
        c1=287.
        c2=0.4
      else
        c1=406.9
        c2=0.62
      endif

!     Compute gas emissivity

      trg=1.8*(tg+273.15)
      gex=xx**c2*c1/trg
      gef=gex
      return
      end function gef

! ---------------------------------------------------------------------

      real function gs(at,cs,ge,se)

! ---------------------------------------------------------------------
!
!     Gas/Sink exchange area
!
!     INPUTS:
!        at    Total area (m2)
!        cs    Ratio of active furnace surface area to total area
!        ge    Gas emissivity
!        se    Active surface emissivity
!
!     OUTPUT:
!        gs    Radiation exchange area (m2)
!
! ---------------------------------------------------------------------

      implicit none
      real         :: se,ge,cs,at,gsx

      gsx=1./ge+1./(cs*se)-1.
      gsx=at/gsx
      gs=gsx
      return
      end function gs

! ----------------------------------------------------------------------

      real function hcovf(re,pr,k,d)

! ----------------------------------------------------------------------
!
!     Convective heat transfer coefficient
!
!     INPUTS:
!       re      Reynolds number
!       pr      Prandtl number
!       k       Thermal conductivity (kW/(m-C))
!       d       Hydaulic diameter (m)
!
!     OUTPUT:
!       hcovf   Convective heat transfer coefficient (kW/(m2-C))
!
! ----------------------------------------------------------------------

      implicit none
      real         :: k,nsn
      real         :: d,re,pr

!     Test if the flow is laminar or turbulent based on Reynolds number.

      if (re>=2000.0) then
        nsn=0.023*re**0.8*pr**0.4
      else
        nsn=3.66
      endif

!     Calculate the convective heat transfer coefficient using the
!     nusselt number.

      hcovf=nsn*k/d

      return
      end function hcovf

! *********************************************************************

      subroutine prdpp(carb,hydr,oxyg,xntr,sulf,xair)

! -------------------------------------------------------
!
!     Combustion products
!
!     INPUTS:
!       pt     Pressure (atm)
!       xair   Excess air fraction
!       carb   Number of moles of carbon
!       hydr   Number of moles of hydrogen
!       oxyg   Number of moles of oxygen
!       xntr   Number of moles of nitrogen
!       sulf   Number of moles of sulfur
!
!     OUTPUTS:
!       pco2   Partial pressure of CO2 (atm)
!       ph2o   Partial pressure of H2O (atm)
!       pso2   Partial pressure of SO2 (atm)
!       pn2    Partial pressure of N2  (atm)
!       po2    Partial pressure of O2  (atm)
!       ratf   Mass ratio of air to fuel (-)
!       rptf   Mass ratio of product to fuel (-)
!       rwtf   Mass ratio of water in product to fuel (-)
!
! *********************************************************************
!
      use blc_head
      implicit none
      real           :: xair,sulf,xntr,oxyg,hydr,carb,r00,rn0,&
                        wtf,wta,wtp,wtw,amult

!     Moles of oxygen and nitrogen to be added

      r00=2.*carb+0.5*hydr+2.*sulf-oxyg
      rn0=3.76*r00

!     Moles of co2, h2o, o2, so2 and n2 in products

      pco2=carb
      ph2o=0.5*hydr
      po2=0.5*r00*xair
      pso2=sulf
      pn2=0.5*xntr+0.5*rn0*(1.+xair)

!     Weight of fuel, air, and products

      wtf=12.*carb+hydr+16.*oxyg+14.*xntr+32.*sulf
      wta=(16.*r00+14.*rn0)*(1.+xair)
      wtp=(44.*pco2+18.*ph2o+32.*po2+64.*pso2+28.*pn2)
      wtw=18.*ph2o

!     Compute partial pressures and mass ratios

      amult=pt/(pco2+ph2o+po2+pso2+pn2)
      pco2=pco2*amult
      ph2o=ph2o*amult
      pso2=pso2*amult
      po2=po2*amult
      pn2=pn2*amult
      ratf=wta/wtf
      rptf=wtp/wtf
      rwtf=wtw/wtf
      return
      end subroutine prdpp

! *********************************************************************

      subroutine prdpr(t,amm,akk)

! ----------------------------------------------------------
!
!     Viscosity and conductivity of combustion product
!
!     INPUTS:
!       pco2   Partial pressure of CO2 (atm)
!       ph2o   Partial pressure of H2O (atm)
!       pso2   Partial pressure of SO2 (atm)
!       pn2    Partial pressure of N2  (atm)
!       po2    Partial pressure of O2  (atm)
!
!     OUTPUTS:
!       amm    Product dynamic viscosity (kg/(m-s))
!       akk    Product thermal conductivity (kW.(m-C))
!
! *********************************************************************

      use blc_head
      implicit none
      real,dimension(11)  :: x,a,b
      real,dimension(5)   :: am,ak
      real                :: akk,amm,t,tf,tab1
      integer             :: np

      data x / -10.,    260.,    440.,    620.,    800.,    1070.,&
               1520.,   2420.,   3140.,   3500.,   4160./
      data a / 0.0387,  0.0553,  0.0646,  0.073,   0.0806,  0.0911,&
               0.1062,  0.1320,  0.1500,  0.1583,  0.1710/
      data b / 0.01287, 0.01944, 0.02333, 0.02692, 0.03022, 0.03483,&
               0.04178, 0.05348, 0.0612,  0.0646,  0.0709/
      data am/ 0.889,   0.685,   1.123,   0.964,   0.889/
      data ak/ 0.925,   0.906,   1.037,   0.983,   0.925/

      np=11
      amm=pco2*am(1)+ph2o*am(2)+po2*am(3)+pn2*am(4)+pso2*am(5)
      akk=pco2*ak(1)+ph2o*ak(2)+po2*ak(3)+pn2*ak(4)+pso2*ak(5)
      tf=32.+1.8*t
      amm=(4.1333e-04)*amm*tab1(np,x,a,tf)
      akk=(1.731e-03)*akk*tab1(np,x,b,tf)
      return
      end subroutine prdpr

! *********************************************************************

      real function taff(hhv,tbs)

! --------------------------------------------------------------------
!
!     Adiabatic flame temperature
!
!     INPUTS:
!        hhv   Higher heating value of fuel (kJ/kg)
!        rwtf  Weight of H2O per kg of fuel (kg/kg)
!        rptf  Weigth of product per kg of fuel (kg/kg)
!        pco2  Number of moles of CO2
!        ph2o  Number of moles of H2O
!        po2   Number of moles of O2
!        pso2  Number of moles of SO2
!        pn2   Number of moles of N2
!        tbs   Base temperature (C)
!
!     OUTPUT:
!        taf   Adiabatic flame temperature (C)
!
! *********************************************************************

      use blc_head
      implicit none
      real       :: tbs,hhv,hvv,taf1,cp,cpf,taf2,t
      integer    :: itr

      itr=0
      hvv=hhv-2442.*rwtf
      taf1=0.995*hvv/rptf+tbs
100   cp=cpf(taf1,tbs)
      taf2=hvv/rptf/cp+tbs
      t=abs(taf1-taf2)
      taf1=taf2
      itr=itr+1
      if(itr<10 .and. t>0.1) then
        goto 100
      endif
      taff =taf1
      return
      end function taff

! *********************************************************************

      real function tab1(np,x,y,x1)

! ---------------------------------------------------------------------
!
!     Determine Y(x1) value from tabulated y vs x values
!
!     INPUTS/OUTPUTS:
!        np   Number of data points
!        x    Data points of independent variables
!        y    Data points of dependent variables
!        tab1 Output dependent variable
!        x1   Input independent variable
!
! *********************************************************************

      implicit none
      integer              :: np,i
      real,dimension(np)   :: x,y
      real                 :: x1,xmin,xmax,y1,dxt,dyt,rdxt,dyodx,dx

!     Set out the range y1(x1) values

      xmin=x(1)
      xmax=x(np)
      if(x1<=xmin) then
        y1=y(1)
      elseif(x1>xmax) then
        y1=y(np)
      else

!     Interpolate for y1(x1) value

        i=1
10      if (x1<=x(i)) then
          dxt=x(i)-x(i-1)
          dyt=y(i)-y(i-1)
          rdxt=1./dxt
          dyodx=dyt*rdxt
          dx=x1-x(i-1)
          y1=y(i-1)+dx*dyodx
        else
          i=i+1
          goto 10
        endif
      endif
      tab1=y1
      return
      end function tab1

