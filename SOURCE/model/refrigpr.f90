!=======================================================================
!
!               REFRIGERANT PROPERTIES
!
!=======================================================================
!
!       Coefficients for refrigerant equations from "Refrigerant
!       Equations", R.C. Downing, Trans. ASHRAE 1974 vol 80 part 2 p158
!
!***********************************************************************
!
     module refrig_head

!     common/refrid/ident(12)
!     common/lde/al(12),bl(12),cl(12),dl(12),el(12),fl(12),gl(12),
!  &  tcr(12),atf(12)
!     common/vpe/a(12),b(12),c(12),d(12),e(12),f(12),atfp(12)
!     common/eee/x(12),y(12)
!     common/hce/aa(12),bb(12),cc(12),dd(12),ff(12),atfppp(12)

!     common/eos/r(12),bv(12),a2(12),b2(12),c2(12),a3(12),
!  &  b3(12),c3(12),a4(12),b4(12),c4(12),a5(12),b5(12),
!  &  c5(12),a6(12),b6(12),c6(12),  ak(12),av(12),cp(12),  tcrp(12),
!  &  atfpp(12)

!     common/eos/r(12),bv(12),a2(12),b2(12),c2(12),
!  &  a3(12),b3(12),c3(12),a4(12),b4(12),c4(12),a5(12),b5(12),
!  &  c5(12),a6(12),b6(12),c6(12),  ck(12),av(12),cp(12),  tc(12),
!  &  atfpp(12)

  integer,parameter   :: nr=12
  real,dimension(nr)  :: ident=(/11,12,13,14,21,22,23,113,114,&
                                500,502,318/),&
!  ident(1):   R-11
!  ident(2):   R-12
!  ident(3):   R-13
!  ident(4):   R-14
!  ident(5):   R-21
!  ident(6):   R-22
!  ident(7):   R-23
!  ident(8):   R-113
!  ident(9):   R-114
!  ident(10):  R-500
!  ident(11):  R-502
!  ident(12):  R-318

        al=(/34.57,34.84,36.06996,39.06,116.37962,&
             32.76,32.7758,122.872,36.32,31.00,35.0,38.70/),&

        bl=(/57.63811,53.341187,54.395124,69.568489,&
             -0.03106808,54.634409,63.37784,-0.0128,61.146414,&
             43.562,53.48437,70.858318/),&

        cl=(/43.6322,0.,0.,4.5866114,-0.0000501,36.74892,&
             -25.30533,0.0000636,0.,74.709,63.86417,23.609759/),&

        dl=(/-42.82356,18.69137,8.512776,36.1716662,0.,&
             -22.2925657,144.16182,0.,16.418015,-87.583,-70.08066,&
             15.989182/),&

        el=(/36.70663,0.,0.,-8.058986,0.,20.4732886,-106.1328,&
             0.,0.,56.483,48.47901,-8.9243856/),&

        fl=(/0.,21.98396,25.879906,0.,0.,0.,&
             0.,0.,17.476838,0.,0.,0./),&

        gl=(/0.,-3.150994,9.589006,0.,0.,0.,0.,&
             0.,1.119828,0.,0.,0./)

  real,dimension(nr)  :: &
        tcr=(/848.07,693.3,543.60,409.50,812.9,664.50,&
              538.33,877.0,753.95,681.59,639.56,699.27/),&

        atf=(/459.67,459.7,459.69,459.69,459.6,459.69,&
              459.69,459.6,459.69,459.69,459.67,459.69/)


  real,dimension(nr)  :: &
        a=(/42.14702865,39.88381727,25.967975,20.71545389,&
            42.7908,29.35754453,328.90853,33.0655,27.071306,&
            17.780935,10.644955,15.63242/),&

        b=(/-4344.343807,-3436.632228,-2709.538217,&
            -2467.505285,-4261.34,-3845.193152,-7952.76913,&
            -4330.98,-5113.7021,-3422.69717,-3671.153813,&
            -4301.063/),&

        c=(/-12.84596753,-12.47152228,-7.17234391,&
            -4.69017025,-13.0295,-7.86103122,-144.5142304,&
            -9.2635,-6.3086761,-3.63691,-0.369835,-2.128401/),&

        d=(/4.0083725e-03,4.73044244e-03,2.545154e-03,&
            6.4798076e-04,3.9851e-03,2.1909390e-03,0.24211502,&
            2.0539e-03,6.913003e-04,5.0272207e-04,&
            -1.746352e-03,-1.19759e-03/),&

        e=(/0.0313605356,0.,0.280301091,0.770707795,&
            0.,0.445746703,-2.1280665e-04,0.,0.78142111,&
            0.4629401,0.8161139,0.6625898/),&

        f=(/862.07,0.,546.00,424.,0.,686.1,9.434955e-08,&
            0.,768.35,695.57,654.,714./)

  real,dimension(nr)  :: &
        atfp=(/459.67,459.7,459.67,459.69,459.6,459.69,&
               459.69,459.6,459.69,459.67,459.67,459.69/),&

        r=(/0.078117,0.088734,0.102728,0.1219336,0.10427,&
            0.124098,0.15327,0.05728,0.062780807,0.10805000,&
            0.096125,0.053645698/),&

        bv=(/0.00190,0.0065093886,0.0048,0.0015,0.,0.002,&
             0.00125,0.,0.005914907,0.006034229,0.00167,0.0060114165/),&

        a2=(/-3.126759,-3.40972713,-3.083417,-2.162959,&
             -7.316,-4.353547,-4.679499,-4.035,-2.3856704,-4.549888,&
             -3.2613344,-1.8947274/),&

        b2=(/1.318523e-03,1.59434848e-03,2.341695e-03,&
             2.135114e-03,4.6421e-03,2.407252e-03,3.472778e-03,&
             2.618e-03,1.0801207e-03,2.308415e-03,2.0576287e-03,&
             9.8484745e-04/),&

        c2=(/-35.76999,-56.7627671,-18.212643,-18.941131,&
             0.,-44.066868,-159.775232,0.,-6.5643648,-92.90748,&
             -24.24879,-28.542156/)

  real,dimension(nr)  :: &
        a3=(/-0.025341,0.0602394465,0.058854,4.404057e-03,&
             -0.20382376,-0.017464,0.012475,-0.0214,0.034055687,&
              0.08660634,0.034866748,0.026479892/),&

        b3=(/4.875121e-05,-1.87961843e-05,-5.671268e-05,&
             1.282818e-05,3.593e-04,7.62789e-05,7.733388e-05,&
             5.00e-05,-5.3336494e-06,-3.141665e-05,&
             -8.6791313e-06,-6.862101e-06/),&

        c3=(/1.220367,1.31139908,0.571958,0.539776,0.,&
             1.483763,5.941212,0.,0.16366057,2.742282,&
             0.33274779,0.66384636/)

  real,dimension(nr)  :: &
        a4=(/1.687277e-03,-5.4873701e-04,-1.026061e-03,&
             1.921072e-04,0.,2.310142e-03,2.068042e-03,0.,&
             -3.857481e-04,-8.726016e-04,-8.5765677e-04,-2.4565234e-04/),&

        b4=(/-1.805062e-06,0.,1.338679e-06,-3.918263e-07,&
             0.,-3.605723e-06,-3.684238e-06,0.,0.,0.,7.0240549e-07,0./),&

        c4=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.022412368,0./)

  real,dimension(nr)  :: &
        a5=(/-2.358930e-05,0.,5.290649e-06,-4.481049e-06,0.,&
             -3.724044e-05,-3.868546e-05,0.,1.6017659e-06,&
             -1.375958e-06,8.8368967e-06,6.0887086e-07/),&

        b5=(/2.448303e-08,3.468834e-09,-7.395111e-09,&
             9.062318e-09,0.,5.355465e-08,6.455643e-08,0.,&
             6.2632341e-10,9.149570e-09,-7.9168095e-09,8.269634e-10/),&

        c5=(/-1.478379e-04,-2.54390678e-05,-3.874233e-05,&
             -4.836678e-05,0.,-1.845051e-04,-7.394214e-04,0.,&
             -1.0165314e-05,-2.102661e-04,-3.7167231e-04,-3.849145e-05/)

  real,dimension(nr)  :: &
        a6=(/1.057504e+08,0.,7.378601e+07,5.838823e+07,0.,&
             1.363387e+08,7.502357e+07,0.,0.,0.,-3.8257766e+07,0./),&

        b6=(/-9.472103e+04,0.,-7.435565e+04,-9.263923e+04,0.,&
             -1.672612e+05,-1.114202e+05,0.,0.,0.,5.5816094e+04,0./),&

        c6=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.5378377e+09,0./)

  real,dimension(nr)  :: &
        ak=(/4.50,5.475,4.00,4.00,0.,4.2,5.50,0.,3.0,5.475,4.2,5./),&
        ck=(/4.50,5.475,4.00,4.00,0.,4.2,5.50,0.,3.0,5.475,4.2,5./),&
        av=(/580.,0.,625.,661.199997,0.,548.2,520.0,0.,0.,0.,609.,0./),&
        cp=(/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,7.e-07,0./),&
        tc=(/848.07,693.3,543.60,409.50,812.9,664.50,538.33,&       ! <<<<
             877.0,753.95,681.59,639.56,699.27/),&
        atfpp=(/459.67,459.7,459.67,459.69,459.69,459.69,&
                459.69,459.69,459.69,459.69,459.67,459.69/)


!     "a" and "b" in equation of state are called "av" and "bv"
!     "c'" is called "cp".
!     "k" is called "ak".

  real,dimension(nr)  :: &
        aa=(/0.023815,8.0945e-03,0.01602,0.0300559282,0.0427,&
             0.02812836,0.07628087,0.07963,0.0175,0.026803537,&
             0.020419,0.0225178157/),&

        bb=(/2.798823e-04,3.32662e-04,2.823e-04,2.3704335e-04,&
             1.40e-04,2.255408e-04,-7.561805e-06,1.159e-04,3.49e-04,&
             2.8373408e-04,2.996802e-04,3.69907814e-04/),&

        cc=(/-2.123734e-07,-2.413896e-07,-1.159e-07,&
             -2.85660077e-08,0.,-6.509607e-08,3.9065696e-07,0.,-1.67e-07,&
             -9.7167893e-08,-1.409043e-07,-1.64842522e-07/),&

        dd=(/5.999018e-11,6.72363e-11,0.,-2.95338805e-11,0.,0.,&
             -2.454905e-10,0.,0.,0.,2.210861e-11,2.152780846e-11/),&

        ff=(/-336.80703,0.,0.,0.,0.,257.341,0.,0.,0.,0.,0.,0./)

  real,dimension(nr)  :: &
        atfppp=(/459.67,459.7,459.69,459.69,459.6,459.69,459.69,&
                 459.6,459.69,459.69,459.67,459.69/)

! lowercase letters "a,b,c,d,f" in heat capacity equation
! are called "aa,bb,cc,dd,ff".

  real,dimension(nr)  :: &
        x=(/50.5418,39.556551,20.911,86.102162,0.,62.4009,0.,&
             25.198,25.3396621,46.4734,35.308,12.19214242/),&

        y=(/-0.0918395,-0.016537936,-0.05676,0.36172528,0.,&
            -0.0453335,0.,-.40552,-0.11513718,-0.09012707564,&
            -0.07444,-0.16828871/)

  real  :: aj=0.185053,tconv1=273.15,tconv2=1.8,tconv3=459.67,&
           tcnv1=273.15,tcnv2=1.8,scnv=4.1868,vcnv=16.01846326,&
           pcnv=0.14503768,hcnv=2.326

  end module refrig_head

!
!***********************************************************************
!
    real function psat(tcels,ir)
!
! ----------------------------------------------------------------------
!
!  Given tsat (Celsius), calculate psat (kPa)
!
!  Conversion factors from psia to kPa, and from Celsius to Kelvin: ???????????
!       p (kPa) = 6.89492 * p (psia)
!       t (Celsius) = ( t (Rankine) / 1.8 ) - 273.15
!  Constants for R-11 (ASHRAE Trans. 1974,vol.80,part 2, page 158)
!
!***********************************************************************

!--        common/vpe/a(12),b(12),c(12),d(12),e(12),f(12),atfp(12)
     use refrig_head
     implicit none
     integer  :: ir
     real     :: tcels,t,clog
     real     :: pconv=6.89476

     t=(tcels+tconv1)*tconv2
     clog=log10(abs(f(ir)-t))
     psat=pconv*10.0**(a(ir)+b(ir)/t+c(ir)*log10(t)+d(ir)*t+&
         e(ir)*(f(ir)-t)/t*clog)

     return
     end function psat
!
!***********************************************************************
!
     real function tsat(pkpa,rtol,ir)
!
! ----------------------------------------------------------------------
!
!   Given vapor psat (kPa), calculate tsat (Celsius)
!
!***********************************************************************

!--        dimension aa(3),bb(3),cc(3),dd(3),ee(3),ff(3)
!--        common/vpe/a(12),b(12),c(12),d(12),e(12),f(12),atfp(12)
     use refrig_head
     implicit none
     integer  :: ir,icor,i
     real     :: pkpa,rtol,p,plog,t0,tg,clog,fn,ft
     real     :: cnvlog=0.43429448,pconv=0.145036
     real,dimension(3)  :: &
                 aa1=(/428.544,350.466,348.394/),&
                 bb1=(/0.084447,.079999,.119359/),&                
                 cc1=(/-3.7483e-05,-3.9864e-05,-6.0199e-05/),&     
                 dd1=(/69.5313,57.5202,47.9408/),&                 
                 ee1=(/12.896,10.4648,18.4989/),&                  
                 ff1=(/3.79718,3.18191,0.0/)                       

! ---The following correlation for tsat(p) is used as initial guess.
! ---the correlation for R-11 was derived for p between 0.2 and 600.
! ---psia, and t between -76.F and 383.F;  s.d.=0.123 in range (45dF)
! ---the correlation for R-12 was derived for t between -140.F and
! ---210.F; s.d.=0.063 in range (170dF)
! ---the correlation for R-500 was derived for t between -80.F and
! ---200.F; s.d.=0.094 in range (136dF)

     icor=ir
     if(ir>2) then
       icor=3
     else
       icor=ir
     end if

     p=pkpa*pconv
     plog=log10(p)
     t0=aa1(icor)+p*(bb1(icor)+p*cc1(icor))+plog*(dd1(icor)+plog*&
        (ee1(icor)+plog*ff1(icor)))

     do i=1,30
       tg=t0

       if(f(ir)>0.0) then
         clog=log10(f(ir)-tg)
       else
         clog=0.0
       end if

       fn=a(ir)+b(ir)/tg+c(ir)*log10(tg)+d(ir)*tg+e(ir)*&
          (f(ir)-tg)/tg*clog-plog                                   
       ft=-b(ir)/tg**2+cnvlog*c(ir)/tg+d(ir)-e(ir)*(cnvlog/tg+&     
          f(ir)*clog/tg**2)                                         
       t0=tg-fn/ft                                                  
                                                                    
!      Note rel. error tolerance is applied to fahrenheit temp.     
                                                                    
       if(abs(t0-tg)<(rtol*abs(t0-tconv3)+0.001)) go to 20
     end do

     write(*,1)p,t0
  1  format(' warning--tsat(p) fails to converge in 30 iterations'/&
            '   p =',g11.4,'  t=?',g11.4)
20   tsat=(t0/tconv2)-tconv1

     return
     end function tsat
!
!***********************************************************************
!
        function tvsat(vm3kg,rtol,ir)
!
! ----------------------------------------------------------------------
!
!  Given saturation specific volume, calculate saturation temp.
!
!***********************************************************************
!
!--        dimension x1(3),x2(3),x3(3),x4(3),x5(3),x6(3),x7(3)
!--        common/vpe/a(12),b(12),c(12),d(12),e(12),f(12),atfp(12)
!        common/eos/r(12),bv(12),a2(12),b2(12),c2(12),
!     &  a3(12),b3(12),c3(12),a4(12),b4(12),c4(12),a5(12),b5(12),
!     &  c5(12),a6(12),b6(12),c6(12),ck(12),av(12),cp(12),tc(12),
!     &  atfpp(12)

     use refrig_head
     implicit none
     integer  :: ir,icor,i
     real     :: rtol,vm3kg,tvsat,vl,t0,vinv,cktc,tg,clog,psl,&  ! <<< 11/28/06
                 dpsldt,ps,dpsdt,ex1,ex2,p,dpdt,fnc,dfncdt
     real     :: cnvlog=0.43429448,vconv=16.01846326
     real,dimension(3)  :: &
                 x1=(/596.107,484.805,486.375/),&
                 x2=(/-158.41188,-128.7714,-128.07559/),&
                 x3=(/42.7895,34.7526,30.9753/),&
                 x4=(/-4.3556,-4.3363,-4.32086/),&
                 x5=(/-7.36803,-5.44105,-2.08663/),&
                 x6=(/5.10203,4.40104,5.02922/),&
                 x7=(/-1.01335,-0.99385,-2.95400/)

! ---Initial guess:  the correlation for R-11 from -80 to 365F;
! ---s.d.=0.179 with 216 d.F.
! ---                the correlation for R-12 from -140 to 210F;
! ---s.d.=0.1043 with 169 d.F.
! ---                the correlation for R-500 from -80 to 200F;
! ---s.d.=0.0215 with 134 d.F.

     if(ir>2) then
       icor=3
     else
       icor=ir
     end if

     vl=log10(vm3kg*vconv)
     t0=x1(icor)+vl*(x2(icor)+vl*(x3(icor)+vl*(x4(icor)+vl*&
        (x5(icor)+vl*(x6(icor)+vl*x7(icor))))))

     vinv=1.0/(vm3kg*vconv-bv(ir))
     cktc=ck(ir)/tc(ir)
     do i=1,50
       tg=t0

       if(tg<=f(ir)) then
         clog=log10(f(ir)-tg)
       else
         clog=0.0
       endif

!      Vapor pressure equation

       psl=a(ir)+b(ir)/tg+c(ir)*log10(tg)+d(ir)*tg+e(ir)*&
          (f(ir)-tg)/tg*clog                                                    
       dpsldt=-b(ir)/tg**2+cnvlog*c(ir)/tg+d(ir)-e(ir)*(cnvlog/&                
              tg+f(ir)*clog/tg**2)                                              
       ps=10.0**psl                                                             
       dpsdt=dpsldt*ps/cnvlog                                                   
                                                                                
!      Equation of state                                                        
                                                                                
       ex1=exp(-cktc*tg)                                                        
       ex2=ex1*cktc                                                             
       p=vinv*(r(ir)*tg+vinv*((a2(ir)+b2(ir)*tg+c2(ir)*ex1)+&                   
         vinv*((a3(ir)+b3(ir)*tg+c3(ir)*ex1)+vinv*((a4(ir)+b4(ir)*&             
         tg+c4(ir)*ex1)+vinv*(a5(ir)+b5(ir)*tg+c5(ir)*ex1)))))                  
       dpdt=vinv*(r(ir)+vinv*((b2(ir)-c2(ir)*ex2)+vinv*((b3(ir)-&               
            c3(ir)*ex2)+vinv*((b4(ir)-c4(ir)*ex2)+vinv*(b5(ir)-c5(ir)*&         
            ex2)))))                                                            
       fnc=ps-p                                                                 
       dfncdt=dpsdt-dpdt                                                        
       t0=tg-fnc/dfncdt                                                         
                                                                                
!      Note rel. error tolerance is applied to Fahrenheit temp.                 
                                                                                
       if(abs(t0-tg).lt.(rtol*abs(t0-tconv3)+0.001))goto 20                     
     end do

     write(*,1)
  1  format(' Warning--tvsat(p) fails to converge in 50 iterations')
 20  tvsat=(t0/tconv2)-tconv1

     return
     end function tvsat
!
!***********************************************************************
!
     real function pgas(vm3kg,tcels,ir)
!
! ----------------------------------------------------------------------
!
!       Given v (m3/kg) and t (Celsius), calculate p (kPa)
!
!       Conversion factors:
!       t (Celsius) = t (Rankine) / 1.8 - 273.15
!       p (kPa) = 6.89492 * p (psia)
!       v (m3/kg) = v (cu.ft/lbm) / 16.018
!
!***********************************************************************

!        common/eos/r(12),bv(12),a2(12),b2(12),c2(12),
!     &  a3(12),b3(12),c3(12),a4(12),b4(12),c4(12),a5(12),b5(12),
!     &  c5(12),a6(12),b6(12),c6(12),ck(12),av(12),cp(12),tc(12),
!     &  atfpp(12)

     use refrig_head
     implicit none
     integer  :: ir
     real     :: vm3kg,tcels,v,t,t1,t2,t3,t4,t5,vinv,p,e1
     real     :: vconv=16.01846326,pconv=6.89492

     v=vm3kg*vconv
     t=(tcels+tconv1)*tconv2
     e1=exp(-ck(ir)*t/tc(ir))
     t1=r(ir)*t
     t2=a2(ir)+b2(ir)*t+c2(ir)*e1
     t3=a3(ir)+b3(ir)*t+c3(ir)*e1
     t4=a4(ir)+b4(ir)*t+c4(ir)*e1
     t5=a5(ir)+b5(ir)*t+c5(ir)*e1

     vinv=1.0/(v-bv(ir))
     p=((((t5*vinv+t4)*vinv+t3)*vinv+t2)*vinv+t1)*vinv
     pgas=p*pconv

     return
     end function pgas

!***********************************************************************
!
     real function vgas(pkpa,tcels,rtol,ir)
!
! ----------------------------------------------------------------------
!
!  Given p (kPa) and t (Celsius), calculate v (m3/kg)
!
!  Conversion factors:
!       t (Celsius) = t (Rankine) / 1.8 - 273.15
!       p (kPa) = 6.89492 * p (psia)
!       v (m3/kg) = v (cu.ft/lbm) / 16.018
!
!***********************************************************************

!        dimension x11(3),x22(3),x33(3),x44(3)
!        common/eos/r(12),bv(12),a2(12),b2(12),c2(12),
!     &  a3(12),b3(12),c3(12),a4(12),b4(12),c4(12),a5(12),b5(12),
!     &  c5(12),a6(12),b6(12),c6(12),ck(12),av(12),cp(12),tc(12),
!     &  atfpp(12)

     use refrig_head
     implicit none
     integer  :: ir,icor,i
     real     :: pkpa,tcels,rtol,p,t,t1,t2,t3,t4,t5,v,vinv0,vinv,&
                 fn,fv,e1
     real     :: vconv=0.062427961,pconv=0.145036
     real,dimension(3)  :: &
                 x11=(/-0.016739,-0.0227961,-0.038263/),&
                 x22=(/-1.48616e-05,-9.68957e-06,-9.922e-06/),&
                 x33=(/3.3399,4.0473,7.7941/),&
                 x44=(/-299.608,-294.012,-579.584/)


     p=pkpa*pconv
     t=(tcels+tconv1)*tconv2
     e1=exp(-ck(ir)*t/tc(ir))
     t1=r(ir)*t
     t2=a2(ir)+b2(ir)*t+c2(ir)*e1
     t3=a3(ir)+b3(ir)*t+c3(ir)*e1
     t4=a4(ir)+b4(ir)*t+c4(ir)*e1
     t5=a5(ir)+b5(ir)*t+c5(ir)*e1

! ---Initial guess:  the correlation for R-11 from 0.4 to 575 psia,
! ---and from tsat to 1059.67R;  s.d.=0.001659 in range 647 d.F.
! ---                the correlation for R-12 from 0.4 to 575 psia,
! ---and from tsat to 1059.67R;  s.d.=0.002143 in range 647 d.F.
! ---                the correlation for R-500 from 0.4 to 575 psia,
! ---and from tsat to 1059.67R;  s.d.=0.002562 in range 647 d.F.


     if(ir>2) then
       icor=3
     else
       icor=ir
     end if

     v=t1/p+x22(icor)*p+x11(icor)+(x33(icor)+x44(icor)/t1)/t1
     vinv0=1.0/(v-bv(ir))

     do i=1,30
       vinv=vinv0
       fn=((((t5*vinv+t4)*vinv+t3)*vinv+t2)*vinv+t1)*vinv-p             
!      fv=((((-5.*t5*vinv-4.*t4)*vinv-3.*t3)*vinv-2.*t2)*vinv-t1)*&     
!         vinv*vinv                                                     
       fv=(((5.0*t5*vinv+4.0*t4)*vinv+3.0*t3)*vinv+2.0*t2)*vinv+t1      
       vinv0=vinv-fn/fv                                                 
       if(abs(vinv0-vinv)<=rtol*vinv0)goto 20
     end do
     write(*,1)p,t
  1  format(' warning--vgas(p,t) unconverged after 30 iterations'/&
       '   p=',g11.4,'  t=',g11.4)
 20  vgas=(1.0/vinv0+bv(ir))*vconv

     return
     end function vgas
!
!***********************************************************************
!
     real function hgas(pkpa,vm3kg,tcels,ir)
!
! ----------------------------------------------------------------------
!
!    Units are kPa, m3/kg, Celsius; hgas in kJ/kg.
!
!    Caution: in converting hgas from btu/lbm to kJ/kg, it's left with
!          english reference point, not converted to si reference point.
!
!***********************************************************************
!
!        common/eos/r(12),bv(12),a2(12),b2(12),c2(12),a3(12),
!     &  b3(12),c3(12),a4(12),b4(12),c4(12),a5(12),b5(12),c5(12),
!     &  a6(12),b6(12),c6(12),ak(12),av(12),cp(12),tc(12),atfpp(12)
!        common/hce/aa(12),bb(12),cc(12),dd(12),ff(12),atfppp(12)
!        common/eee/x(12),y(12)

     use refrig_head
     implicit none
     integer  :: ir
     real     :: pkpa,vm3kg,tcels,fa3,fa4,fa5,fc3,fc5,fb,fc,fd,&
                 p,t,v,vi,ck1

     fa3=0.5*a3(ir)
     fa4=(1.0/3.0)*a4(ir)
     fa5=0.25*a5(ir)
     fc3=0.5*c3(ir)
     fc5=0.25*c5(ir)
     fb=0.5*bb(ir)
     fc=(1.0/3.0)*cc(ir)
     fd=0.25*dd(ir)

     p=pkpa*pcnv
     t=(tcels+tcnv1)*tcnv2
     v=vm3kg*vcnv

!----Eq. from ASHRAE Trans vol 80 part 2 1974 pp 158-169;
!--  Consts have been divided by 2,3, or 4 to simplefy eq.

     vi=1.0/(v-bv(ir))
     ck1=ak(ir)*t/tc(ir)
     hgas=x(ir)-ff(ir)/t+(((fd*t+fc)*t+fb)*t+aa(ir))*t+aj*&
         (p*v+(((fa5*vi+fa4)*vi+fa3)*vi+a2(ir))*vi+(1.0+&
         ck1)*exp(-ck1)*(((fc5*vi+c4(ir))*vi+fc3)*vi+c2(ir))*vi)
     hgas=hgas*hcnv

     return
     end function hgas
!
!***********************************************************************
!
     real function sgas(vg,tcels,ir)
!
! ----------------------------------------------------------------------
!
!   Units are m3/kg, Celsius, and kJ/(kg-K).
!
!   Caution: in converting sgas from btu/(lbm-r), it's left with
!            english reference point, not si reference point.
!
!***********************************************************************
!
!        common/eos/r(12),bv(12),a2(12),b2(12),c2(12),a3(12),
!     &  b3(12),c3(12),a4(12),b4(12),c4(12),a5(12),b5(12),c5(12),
!     &  a6(12),b6(12),c6(12),ak(12),av(12),cp(12),tc(12),atfpp(12)
!        common/hce/aa(12),bb(12),cc(12),dd(12),ff(12),atfppp(12)
!        common/eee/x(12),y(12)

     use refrig_head
     implicit none
     integer  :: ir
     real     :: tcels,vg,fkc,fkd,fkf,fkb5,fkb3,fkb4,fkc3,&
                 fkc5,t,v,vi

     fkc=0.5*cc(ir)
     fkd=(1.0/3.0)*dd(ir)
     fkf=0.5*ff(ir)
     fkb5=0.25*b5(ir)
     fkb3=0.5*b3(ir)
     fkb4=(1.0/3.0)*b4(ir)
     fkc3=0.5*c3(ir)
     fkc5=0.25*c5(ir)

     t=(tcels+tcnv1)*tcnv2
     v=(vg*vcnv)-bv(ir)
     vi=1.0/v
     sgas=y(ir)+aa(ir)*log(t)-fkf/t**2+((fkd*t+fkc)*&
          t+bb(ir))*t+aj*(r(ir)*log(v)-(((fkb5*vi+fkb4)*vi+&
          fkb3)*vi+b2(ir))*vi+ak(ir)/tc(ir)*exp(-ak(ir)*t/tc(ir))*&
          (((fkc5*vi+c4(ir))*vi+fkc3)*vi+c2(ir))*vi)
     sgas=sgas*scnv

     return
     end function sgas
!
!***********************************************************************
!
     real function hps(pkpa,skj,rtol,ir)
!
!***********************************************************************
!
!        dimension a1(3),a2(3),a3(3),a4(3),a5(3),a6(3),a7(3)
!        common/eos/r(12),bv(12),fa2(12),b2(12),c2(12),fa3(12),
!     &  b3(12),c3(12),fa4(12),b4(12),c4(12),fa5(12),b5(12),c5(12),
!     &  fa6(12),b6(12),c6(12),ak(12),av(12),cp(12),tc(12),atfpp(12)
!        common/hce/aa(12),bb(12),cc(12),dd(12),ff(12),atfppp(12)

     use refrig_head
     implicit none
     integer  :: ir,icor,i
     real     :: pkpa,skj,rtol,fc3,fc5,p,s,plog,t0,t,ts,tsat,tcel,&
                 vm3kg,vgas,v,vi,fn,sgas,dfndt,hgas
     real,dimension(3)  ::  &
                 a11=(/208.874,91.129,-9.883/),&
                 a21=(/-1129.03,-110.28,535.055/),&
                 a31=(/10713.9,7848.37,4361.40/),&
                 a41=(/639.463,588.596,444.853/),&
                 a51=(/-4.3971,-0.5176,3.4711/),&
                 a61=(/3.1021,3.0040,3.1513/),&
                 a71=(/-0.0532,-0.1525,-0.1817/)

     fc3=0.5*c3(ir)
     fc5=0.25*c5(ir)

     p=pkpa*pcnv
     s=skj/scnv
     plog=log10(p)

! ---Initial guess:  the correlation for R-11 from 0.4 to 575 psia,
! ---and from tsat to 1059.67R;  s.d.=1.049 in range 329 d.F.
! ---                the correlation for R-12 from 0.4 to 575 psia,
! ---and from tsat to 1059.67R;  s.d.=1.170 in range 329 d.F.
! ---                the correlation for R-500 from 0.4 to 575 psia,
! ---and from tsat to 1059.67R;  s.d.=0.5654 in range 329 d.F.

     if(ir>2) then
       icor=3
     else
       icor=ir
     end if

     t0=a11(icor)+s*(a21(icor)+s*a31(icor)+plog*a41(icor))+plog*&
        (a51(icor)+plog*(a61(icor)+plog*a71(icor)))/s

     do i=1,30
       t=t0

       if(t<ts) then
         t=ts
       else
         ts=(tsat(pkpa,rtol,ir)+tcnv1)*tcnv2
       end if

       tcel=t/tcnv2-tcnv1                                                    
       vm3kg=vgas(pkpa,tcel,rtol,ir)                                         
       v=vm3kg*vcnv                                                          
       vi=1.0/(v-bv(ir))
       fn=sgas(vm3kg,tcel,ir)/scnv-s                                         
       dfndt=(dd(ir)*t+cc(ir))*t+bb(ir)+aa(ir)/t+ff(ir)/t**3-aj*(ak(ir)/&
             tc(ir))**2*exp(-ak(ir)*t/tc(ir))*(((fc5*vi+c4(ir))*vi+&
             fc3)*vi+c2(ir))*vi
       t0=t-fn/dfndt                                                         
                                                                             
!      Note rel. error tolerance is applied to Fahrenheit temp.              
                                                                             
       if(abs(t-t0)<(rtol*abs(t0-tconv3)+0.001))goto 20
     end do
     write(*,1) p,s
  1  format(' warning--hps(p,s) not converged in 30 iterations'/&
            '   p= ',g11.4,' s= ',g11.4)
  20 tcel=t0/tcnv2-tcnv1
     vm3kg=vgas(pkpa,tcel,rtol,ir)
     hps=hgas(pkpa,vm3kg,tcel,ir)

     return
     end function hps
!
!***********************************************************************
!
     real function tph(pkpa,hkjkg,rtol,ir)
!
!***********************************************************************
!
!        common/eos/r(12),bv(12),fa2(12),b2(12),c2(12),fa3(12),
!     &  b3(12),c3(12),fa4(12),b4(12),c4(12),fa5(12),b5(12),c5(12),
!     &  fa6(12),b6(12),c6(12),ak(12),av(12),cp(12),tc(12),atfpp(12)
!        common/hce/aa(12),bb(12),cc(12),dd(12),ff(12),atfppp(12)

     use refrig_head
     implicit none
     integer  :: ir,icor,i
     real     :: pkpa,hkjkg,rtol,fkec3,fkec5,p,h,plog,ts,vs,hs,tsat,&
                 vgas,hgas,dh,tseu,t0,t,tcel,vm3kg,v,vi,fn,dfndt
     real,dimension(3)  ::  &
                 x1=(/7.93703,8.2267,7.06948/),&
                 x2=(/-0.0207321,-0.027269,-0.0217131/),&
                 x3=(/7.4643e-05,9.0386e-05,5.67e-05/),&
                 x4=(/-0.265717,-0.268452,-0.24733/),&
                 x5=(/-0.36136,-0.42599,-0.37881/),&
                 x6=(/0.0043632,0.0043449,0.0031129/)

     fkec3=0.5*c3(ir)
     fkec5=0.25*c5(ir)

     p=pkpa*pcnv
     h=hkjkg/hcnv
     plog=log10(p)

! ---Initial guess:  the correlation for R-11 from 0.4 to 575 psia,
! ---and from tsat to 1059.67R;  s.d.=0.8744 in range 250 d.F.
! ---                the correlation for R-12 from 0.4 to 575 psia,
! ---and from tsat to 1059.67R;  s.d.=1.436 in range 250 d.F.
! ---                the correlation for R-500 from 0.4 to 575 psia,
! ---and from tsat to 1059.67R;  s.d.=1.542 in range 250 d.F.

     if(ir>2) then
       icor=3
     else
       icor=ir
     endif

     ts=tsat(pkpa,rtol,ir)
     vs=vgas(pkpa,ts,rtol,ir)
     hs=hgas(pkpa,vs,ts,ir)/hcnv

     if(h<hs) then
       tph=ts/tcnv2-tcnv1
       return
     end if

     dh=h-hs
     tseu=(ts+tcnv1)*tcnv2
     t0=tseu+dh*((x1(icor)+dh*(x2(icor)+dh*x3(icor)))+plog*&
        (x4(icor)+plog*(x5(icor)+dh*x6(icor))))

     do i=1,30
       t=t0

       if(t<ts) then
         t=ts
       else
         ts=(tsat(pkpa,rtol,ir)+tcnv1)*tcnv2
       endif

       tcel=t/tcnv2-tcnv1                                                     
       vm3kg=vgas(pkpa,tcel,rtol,ir)                                          
       v=vm3kg*vcnv                                                           
       vi=1.0/(v-bv(ir))
       fn=hgas(pkpa,vm3kg,tcel,ir)/hcnv-h                                     
       dfndt=ff(ir)/t**2+((dd(ir)*t+cc(ir))*t+bb(ir))*t+aa(ir)-aj*&
             (ak(ir)/tc(ir))**2*t*exp(-ak(ir)*t/tc(ir))*(((fkec5*vi+&
              c4(ir))*vi+fkec3)*vi+c2(ir))*vi
       t0=t-fn/dfndt                                                          
                                                                              
!      Note rel. error tolerance is applied to Fahrenheit temp.               
                                                                              
       if(abs(t-t0)<(rtol*abs(t0-tconv3)+0.001))goto 20
     end do
     write(*,1)p,h
  1  format(' Warning--tph(p,h) fails to converge in 30 iterations'/&
            '   p= ',g11.4,' h= ',g11.4)
  20 tph=t0/tcnv2-tcnv1

     return
     end function tph
!
!***********************************************************************
!
     real function tvh(vm3kg,hkjkg,rtol,ir)
!
!***********************************************************************
!
!        dimension x1(3),x2(3),x3(3),x4(3),x5(3),x6(3),x7(3),
!     &  x8(3),x9(3),x10(3)
!        common/eos/r(12),bv(12),fa2(12),b2(12),c2(12),fa3(12),
!     &  b3(12),c3(12),fa4(12),b4(12),c4(12),fa5(12),b5(12),c5(12),
!     &  fa6(12),b6(12),c6(12),ak(12),av(12),cp(12),tc(12),atfpp(12)
!        common/hce/aa(12),bb(12),cc(12),dd(12),ff(12),atfppp(12)

     use refrig_head
     implicit none
     integer  :: ir,icor,i
     real     :: vm3kg,hkjkg,rtol,fakec3,fakec5,h,v,vlog,t0,t,tcel,&
                 pkpa,pgas,vi,fn,hgas,dfndt
     real,dimension(3)  ::  &
                 x1=(/-443.043,-281.582,-315.191/),&
                 x2=(/13.0283,12.3854,11.0026/),&
                 x3=(/-0.0405824,-0.0430709,-0.0336587/),&
                 x4=(/7.87e-05,9.0069e-05,5.66285e-05/),&
                 x5=(/-40.165,-34.544,-37.647/),&
                 x6=(/10.437,9.153,9.886/),&
                 x7=(/-1.2046,-1.0396,-1.0856/),&
                 x8=(/0.29215,0.28148,0.24745/),&
                 x9=(/-0.026045,-0.02607,-0.02337/),&
                 x10=(/6.7734e-04,-7.153e-04,-5.0261e-04/)

     fakec3=0.5*c3(ir)
     fakec5=0.25*c5(ir)

     if(ir>2) then
       icor=3
     else
       icor=ir
     end if

     h=hkjkg/hcnv
     v=vm3kg*vcnv
     vlog=log10(v)

! ---Initial guess:  the correlation for R-11 from 0.4 to 575 psia,
! ---and from tsat to 1059.67R;  s.d.=0.3627 in range 246 d.F.
! ---                the correlation for R-12 from 0.4 to 575 psia,
! ---and from tsat to 1059.67R;  s.d.=0.7989 in range 246 d.F.
! ---                the correlation for R-500 from 0.4 to 575 psia,
! ---and from tsat to 1059.67R;  s.d.=0.8470 in range 246 d.F.

     t0=x1(icor)+h*(x2(icor)+h*(x3(icor)+h*x4(icor)))+&
        vlog*((x5(icor)+vlog*(x6(icor)+vlog*x7(icor)))+&
        h*(x8(icor)+vlog*x9(icor)+h*x10(icor)))

     do i=1,30
       t=t0
       tcel=t/tcnv2-tcnv1                                                        
       pkpa=pgas(vm3kg,tcel,ir)                                                  
!      ps=psat(tcel,ir)                                                          
!      if(pkpa.gt.ps)pkpa=ps                                                     
       vi=1.0/(v-bv(ir))
       fn=hgas(pkpa,vm3kg,tcel,ir)/hcnv-h                                        
       dfndt=ff(ir)/t**2+((dd(ir)*t+cc(ir))*t+bb(ir))*t+aa(ir)-aj*&
             (ak(ir)/tc(ir))**2*t*exp(-ak(ir)*t/tc(ir))*(((fakec5*&
             vi+c4(ir))*vi+fakec3)*vi+c2(ir))*vi
       t0=t-fn/dfndt                                                             
                                                                                 
!      Note rel. error tolerance is applied to Fahrenheit temp.                  
                                                                                 
       if(abs(t-t0)<(rtol*abs(t0-tconv3)+.001))goto 20
     end do
     write(*,1)v,h
  1  format(' Warning--tvh(v,h) fails to converge in 30 iterations'/&
            '   v= ',g11.4,' h= ',g11.4)
  20 tvh=t0/tcnv2-tcnv1

     return
     end function tvh
!
!***********************************************************************
!
     real function dhlat(psatsi,vsatsi,tsatsi,ir)
!
! ----------------------------------------------------------------------
!
!  Latent heat of vaporization.  constants are for R-11.
!
!***********************************************************************
!
!        common/vpe/a(12),b(12),c(12),d(12),e(12),f(12),atfp(12)

     use refrig_head
     implicit none
     integer  :: ir,icor,i
     real     :: psatsi,vsatsi,tsatsi,psat,tsat,vsat,tinv,flog,rholiq
     real     :: ln10=2.302585093,log10e=0.4342944819

     psat=psatsi*pcnv
     tsat=(tsatsi+tcnv1)*tcnv2
     vsat=vsatsi*vcnv
     tinv=1.0/tsat

     if(f(ir)>0.0) then
       flog=log10(f(ir)-tsat)
     else
       flog=0.0
     end if

     dhlat=hcnv*(vsat-1.0/rholiq(tsatsi,ir))*aj*tsat*(psat*ln10*&
           (d(ir)+tinv*(-b(ir)*tinv+c(ir)/ln10-e(ir)*(log10e+tinv*&
           f(ir)*flog))))

     return
     end function dhlat

!***********************************************************************
!
     real function rholiq(tcel,ir)
!
!***********************************************************************

!        common/lde/al(12),bl(12),cl(12),dl(12),el(12),fl(12),
!     &  gl(12),tc(12),atf(12)

     use refrig_head
     implicit none
     integer  :: ir
     real     :: tcel,t,tf
     real     :: rhocnv=16.01846326
     real     :: e11=0.3333333333,e21=0.6666666667,e31=1.333333333

     t=(tcel+tcnv1)*tcnv2
     tf=1.0-t/tc(ir)
     rholiq=al(ir)+bl(ir)*tf**e11+cl(ir)*tf**e21+dl(ir)*tf+el(ir)*&
     tf**e31+fl(ir)*sqrt(tf)+gl(ir)*tf*tf
     rholiq=rholiq*rhocnv

     return
     end function rholiq
!
!***********************************************************************
!
     real function cv(vm3kg,tcel,ir)
!
! ----------------------------------------------------------------------
!
!  Specific heat of vapor at constant volume.
!
!***********************************************************************
!
!        common/eos/r(12),bv(12),a2(12),b2(12),c2(12),
!     &  a3(12),b3(12),c3(12),a4(12),b4(12),c4(12),a5(12),b5(12),
!     &  c5(12),a6(12),b6(12),c6(12),ak(12),av(12),cp(12),tc(12),
!     &  atfpp(12)
!        common/hce/aa(12),bb(12),cc(12),dd(12),ff(12),atfppp(12)
     use refrig_head
     implicit none
     integer  :: ir
     real     :: vm3kg,tcel,t,v,vinv
     real     :: cvcnv=4.1868

     t=(tcel+tcnv1)*tcnv2
     v=vm3kg*vcnv
     vinv=1.0/(v-bv(ir))
     cv=aa(ir)+t*(bb(ir)+t*(cc(ir)+dd(ir)*t))+ff(ir)/(t*t)-aj*&
        ak(ir)*ak(ir)*t*exp(-ak(ir)*t/tc(ir))/(tc(ir)*tc(ir))*&
        vinv*(c2(ir)+vinv*(c3(ir)/2.0+vinv*(c4(ir)/3.0+vinv*c5(ir)/&
        4.0)))
     cv=cv*cvcnv

     return
     end function cv
!
!***********************************************************************
!
     subroutine cpcv(vm3kg,tcel,ir,cp1,cv,gamma,sonic)
!
! ----------------------------------------------------------------------
!
!  Given specific volume and temperature, calculate cp, cv, cp/cv,
!   and sonic velocity.
!
!***********************************************************************
!
!        common/eos/r(12),bv(12),a2(12),b2(12),c2(12),
!     &  a3(12),b3(12),c3(12),a4(12),b4(12),c4(12),a5(12),b5(12),
!     &  c5(12),a6(12),b6(12),c6(12),ak(12),av(12),cpr(12),tc(12),
!     &  atfpp(12)
!        common/hce/aa(12),bb(12),cc(12),dd(12),ff(12),atfppp(12)
     use refrig_head
     implicit none
     integer  :: ir
     real     :: vm3kg,tcel,cp1,cv,gamma,sonic,t,v,vinv,aktc,ex1,ex2,&
                 dpdv,dpdt
     real     :: cvcnv=4.1868

     t=(tcel+tcnv1)*tcnv2
     v=vm3kg*vcnv
     vinv=1.0/(v-bv(ir))
     aktc=ak(ir)/tc(ir)
     ex1=exp(-aktc*t)
     ex2=aktc*ex1
     cv=aa(ir)+t*(bb(ir)+t*(cc(ir)+dd(ir)*t))+ff(ir)/(t*t)-aj*aktc*ex2*&
        vinv*(c2(ir)+vinv*(c3(ir)/2.+vinv*(c4(ir)/3.+vinv*c5(ir)/4.0)))
     dpdv=-vinv*vinv*(r(ir)*t+vinv*(2.0*(a2(ir)+b2(ir)*t+c2(ir)*ex1)+&
           vinv*(3.0*(a3(ir)+b3(ir)*t+c3(ir)*ex1)+vinv*&
           (4.0*(a4(ir)+b4(ir)*t+c4(ir)*ex1)+vinv*5.0*(a5(ir)+b5(ir)*t+&
           c5(ir)*ex1)))))
     dpdt=vinv*(r(ir)+vinv*((b2(ir)-c2(ir)*ex2)+vinv*((b3(ir)-&
          c3(ir)*ex2)+vinv*((b4(ir)-c4(ir)*ex2)+vinv*(b5(ir)-c5(ir)*&
          ex2)))))
     cp1=cv-aj*t*dpdt*dpdt/dpdv
     gamma=cp1/cv
     sonic=0.3048*v*sqrt(857.36091*t*dpdt*dpdt/cv-4633.056*dpdv)

! ***Includes conversion to meters/sec from ft/sec.

     cv=cv*cvcnv
     cp1=cp1*cvcnv

     return
     end subroutine cpcv
