! **********************************************************************

      real function hyster(ca,cr,cro,told,hys)

! ---------------------------------------------------------------------
!
!       Models hysteresis of actuators.
!         ca  : rel. position of actuator (0<=ca<=1)
!         cr  : rel. position of valve or damper [0<=cr<=(1.-hys)]
!         cro : value of cr at time told
!         told: previous value of time
!         hyster: cr scaled such that 0<=hyster<=1
!
! **********************************************************************

      use modsim_head
      implicit none
      real :: ca,cr,cro,told,hys,c,slope,slack

      if(itime<=1 .and. init==0) then
        told=time
                                  
        if(cro<0.0) then          
          cro=0.0                 
        else                      
          cro=ca-hys              
        endif                     
                                  
        cr=cro                    
        hyster=cr/(1.0-hys)
        return                    
      endif

      if(time>told) then
        cro=cr
      endif

      told=time
      c=max(ca,0.0)
      c=min(c,1.0)
      slope=1.0/(1.0-hys)
      slack=c-cro

      if(slack>hys) then
        cr=c-hys
      elseif(slack<0.0) then
        cr=c
      else
        cr=cro
      endif

      hyster=slope*cr

      return
      end function hyster
! ***********************************************************************

      real function delay(tinew,tau,a1,tin,dzold,time0)

! ---------------------------------------------------------------------
!
!     Models transport delays in ducting components.
!     the temperature distribution in the duct is modeled by a fifth
!     order time dependent polynomial. The coefficients of the
!     polynomial are calculated at each time step.
!
!     Modified:  March 5, 1987 c.p.
!
! ***********************************************************************

      use modsim_head
      implicit none
      integer                   :: i,j,k,l,i2
      real                      :: tinew,tau,tin,dzold,time0,x,dz,xx,r,z
      integer,parameter         :: m=5,n=m+1
      real,dimension(n)         :: a1
      real,dimension(m,n)       :: am1

!     Initialize coefficients if this is the first call.
!     If the transport delay is zero, return the input value.

      if(itime==1 .and. init==0) then
        tin=tinew
        a1(1)=tin
        time0=time
        delay=tin

        if(tau>0.0) then
          dzold=tstep/tau
        else
          dzold=2.0
        endif

        do i=2,n
          a1(i)=0.
        enddo

        return
      endif

!     Set the time step equal to minimum time step if tstep < tmin

      if(tstep<tmin) then
        tstep=tmin
      endif

      if(tau>0.0) then
        if(tau>1.0e+10) then
           tau=1.0e+10
        endif
        dz=tstep/tau
      else
        dz=2.0
      endif

!     If first call of the timestep, evaluate new coefficients.

      if(time>time0) then

!     Set up matrix and constant vector to solve for new set of
!     coefficients.

        if(dzold <= 0.5) then
          do i=1,m
            am1(i,n)=0.0
            x=1.0
            xx=1.0

            if(i==1) then
              if(dzold>=0.2) then
                z=0.5*dzold
              else
                z=dzold
              endif
            elseif(i==2 .and. dzold>=0.2) then
                z=dzold
            else
              i2=i
              z=0.2*real(i2)
            endif

            do j=1,m
              xx=xx*z
              am1(i,j)=xx
              if(z>dzold) then
                x=x*(z-dzold)
                am1(i,n)=am1(i,n)+a1(j+1)*x
              endif
            enddo

            if(z<=dzold) then
              am1(i,n)=(a1(1)-tin)*z/dzold
            else
              am1(i,n)=am1(i,n)+a1(1)-tin
            endif
          enddo

!     Solve for coefficients by gaussian elimination.

          do i=2,m
            do j=i,m
              r=am1(j,i-1)/am1(i-1,i-1)
              do k=i,n
                am1(j,k)=am1(j,k)-r*am1(i-1,k)
              enddo
            enddo
          enddo

          do i=2,m
            k=n+1-i
            r=am1(k,n)/am1(k,k)
            do j=i,m
              l=n-j
              am1(l,n)=am1(l,n)-r*am1(l,k)
            enddo
          enddo

          a1(1)=tin
          do i=1,m
            a1(i+1)=am1(i,n)/am1(i,i)
          enddo

!     Last dz between 0.5 and 1: reset to linear temperature
!     distribution.

        elseif(dzold > 0.5 .and. dzold < 1.0) then
          do i=2,n
            a1(i)=0.0
          enddo
          a1(2)=a1(1)-tin
          a1(1)=tin

!     Last dz beyond inlet: reset to uniform temperature.

        else
          do i=2,n
            a1(i)=0.0
          enddo
          a1(1)=tin
        endif
      endif

!     Calculate output.

      time0=time
      tin=tinew
      dzold=dz
      if(dz>=1.0) then
        delay=tinew
      else
        delay=a1(1)
        x=1.0
        do i=2,n
          x=x*(1.0-dz)
          delay=delay+a1(i)*x
        enddo
      endif

      return
      end function delay

! =======================================================================
!     The following lines of code consist of 4 subroutines:
!       sufed, besi, besk, and polfit.
!     These subroutines are used once per cooling coil per simulation,
!     at time=0, to obtain coeff's for dry fin efficiency equation.
! ***********************************************************************

      subroutine sufed(roh,cof)

! ***********************************************************************
      implicit none
      integer                   :: i,ie1,ie2,ie3,ie4,ie5,ie6
      real                      :: roh,fai,r1,r2,ro,r1i1,r2k1,r2i1,&
                                   k1,r2i0,r2k0,fed,r1k1
      real,dimension(60)        :: x,y
      real,dimension(10)        :: cof

!     roh = (outside tube diam.)/(effective fin diam.)

      fai=0.02
      do i=1,60
        fai=fai+0.035
        r1=fai/(1.0-roh)
        r2=r1*roh
        ro=2.0*roh/(fai*(1.+roh))
        call besi(r1,1,r1i1,ie1)
        call besk(r2,1,r2k1,ie2)
        call besi(r2,1,r2i1,ie3)
        call besk(r1,1,r1k1,ie4)
        call besi(r2,0,r2i0,ie5)
        call besk(r2,0,r2k0,ie6)
        fed=ro*(r1i1*r2k1-r2i1*r1k1)/(r2i0*r1k1+r1i1*r2k0)
        x(i)=fai
        y(i)=fed
      enddo

      call polfit(x,y,60,1.e-05,4,cof)

      return
      end subroutine sufed
! ***********************************************************************

      subroutine besi(x,n,bi,ier)

! ----------------------------------------------------------------------
!
!     Subroutine besi to calculate the modified Bessel function
!      of order 0 to n
!
!     call besi(x,n,bi,ier)
!
!      x    argument of Bessel function
!      n    order of Bessel function, greater than or equal to zero
!      bi   resultant value of i Bessel function
!      ier  resultant error code:
!          ier = 0   no error
!          ier = 1   n < 0
!          ier = 2   x < 0
!          ier = 3   bi < 10**(-30),     bi is set to 0
!          ier = 4   x > n & x > 90,     bi is set to 10**38
!
! ***********************************************************************

      implicit none
      integer                   :: n,ier,i,k
      real                      :: x,bi,fn,xx,term,fk,fi
      real                      :: pi=3.141592653,tol=1.0e-06

      ier=0
      bi=1.
      if(x == 0.0 .and. n == 0) return
      if(n < 0) then
        ier=1
        return
      else if (x < 0.0) then
        ier = 2
        return
      else if (x > 12.0 .and. x > n) then
        if (x > 90.0) then
          ier = 4
          bi=10.0**30
          return
        endif
        fn=4*n*n
        xx=0.125/x
        term=1.0
        bi=1.0
        do k=1,30
          if(abs(term)<=abs(tol*bi)) then
              bi=bi*exp(x)/sqrt(2.0*pi*x)
              return
          endif
          fk=(2*k-1)**2
          term=term*xx*(fk-fn)/real(k)
          bi=bi+term
        enddo
      endif
      xx=x/2.0
      term=1.0
      if(n>0) then
        do i=1,n
          fi=i
          if(abs(term)<1.e-30*fi/xx) then
            ier=3
            bi=0.0
            return
          endif
          term=term*xx/fi
        enddo
      endif

      bi=term
      xx=xx*xx
      do k=1,1000
        if(abs(term)<=abs(bi*tol)) return
        fk=k*(n+k)
        term=term*(xx/fk)
        bi=bi+term
      enddo

      return
      end subroutine besi
! ***********************************************************************

      subroutine besk(x,n,bk,ier)

! ----------------------------------------------------------------------
!
!     Subroutine besk  to compute the k Bessel function for a given
!                      argument and order
!
!      call besk(x,n,bk,ier)
!
!       x    the argument of the k Bessel function desired
!       n    the order of the k Bessel function desired
!       bk   the resultant k Bessel function
!       ier  resultant error code:
!             ier=0  no error
!             ier=1  n is negative
!             ier=2  x is zero or negative
!             ier=3  x > 85, bk < 10**-38; bk set to 0.
!             ier=4  bk > 10**38; bk set to 10**38
!
!      Note: n must be greater than or equal to zero
!
!      Method:
!       Computes zero order and first order Bessel functions using
!       series approximations and then computes n th order function
!       using recurrence relation.
!       Recurrence relation and polynomial approximation technique
!       as described by A.J.M. Hitchcock, 'Polynomial approximations
!       to Bessel functions of order zero and one and to related
!       functions,' M.T.A.C., v.11, 1957, pp. 86-88, and G.N. Watson,
!       'A treatise on the theory of Bessel functions,' Cambridge
!       University Press, 1958, p.62
!
! ***********************************************************************

      implicit none
      integer                   :: n,ier,j,l
      real                      :: x,bk,g0,gj,a,b,c,g1,x2j,fact,hj,rj
      real,dimension(12)        :: t
      real                      :: gjmax=1.0e38

      bk=0.0
      g0=0.0
      gj=0.0

      if(n<0) then
        ier=1
        return
      else if(x<=0.0) then
        ier=2
        return
      else if(x>85.0) then
        ier=3
        bk=0.0
        return
      endif

      ier=0

!     Use polynomial approximation if x > 1.

      if(x>1.0) then
        a=exp(-x)
        b=1.0/x
        c=sqrt(b)
        t(1)=b

        do l=2,12
          t(l)=t(l-1)*b
        enddo

        if(n/=1) then

!     Compute k0 using polynomial approximation

          g0=a*(1.2533141-.1566642*t(1)+0.08811128*t(2)-0.09139095&
           *t(3)+0.1344596*t(4)-0.2299850*t(5)+0.3792410*t(6)-0.5247277&
           *t(7)+0.5575368*t(8)-0.4262633*t(9)+0.2184518*t(10)&
           -0.06680977*t(11)+0.009189383*t(12))*c
          if(n == 0) then
            bk=g0
            return
          endif
        endif

!   Compute k1 using polynomial approximation

        g1=a*(1.2533141+0.4699927*t(1)-0.1468583*t(2)+0.1280427*t(3)&
         -0.1736432*t(4)+0.2847618*t(5)-0.4594342*t(6)+0.6283381*t(7)&
         -0.6632295*t(8)+0.5050239*t(9)-0.2581304*t(10)+0.07880001*t(11)&
         -0.01082418*t(12))*c
        if(n == 1) then
          bk=g1
          return
        endif
      else

!   Use series expansion if x <= 1.

        b=x/2.0
        a=0.5772157+log(b)
        c=b*b
        if(n/=1) then

!   Compute k0 using series expansion

          g0=-a
          x2j=1.0
          fact=1.0
          hj=0.0
          do j=1,6
            rj=1.0/real(j)
            x2j=x2j*c
            fact=fact*rj*rj
            hj=hj+rj
            g0=g0+x2j*fact*(hj-a)
          enddo
          if(n==0) then
              bk=g0
              return
          endif
        endif

!   Compute k1 using series expansion

        x2j=b
        fact=1.0
        hj=1.0
        g1=1.0/x+x2j*(0.5+a-hj)
        do j=2,8
          x2j=x2j*c
          rj=1.0/real(j)
          fact=fact*rj*rj
          hj=hj+rj
          g1=g1+x2j*fact*(0.5+(a-hj)*real(j))
        enddo
        if(n==1) then
            bk=g1
            return
        endif
      endif

!   From k0 and k1 compute kn using recurrence relation

      do j=2,n
        gj=2.0*(real(j)-1.0)*g1/x+g0

        if(gj>=gjmax) then
          ier=4
          gj=gjmax
          exit
        else
          g0=g1
          g1=gj
        endif
      enddo

      bk=gj

      return
      end subroutine besk
! ***********************************************************************

      subroutine polfit(x,y,n,tol,last,cof)

! ---------------------------------------------------------------------
!
!   Subroutine polfit fits polynomial of order from 1 to last to the
!   ordered pairs of data points x,y
!
! ***********************************************************************

      implicit none
      integer                   :: n,last,i,j,l,nord,kk,ik,kkk,ii
      real                      :: tol,b,c,s1,s2
      real,dimension(10)        :: sumx,sumy,cof
      real,dimension(60)        :: x,y
      real,dimension(10,10)     :: a

      sumx(1)=n
      sumx(2)=0.0
      sumx(3)=0.0
      sumy(1)=0.0
      sumy(2)=0.0

      l=last+1
      do i=1,l
        cof(i)=0.0
      enddo

      do i=1,n
        sumx(2)=sumx(2)+x(i)
        sumx(3)=sumx(3)+x(i)*x(i)
        sumy(1)=sumy(1)+y(i)
        sumy(2)=sumy(2)+x(i)*y(i)
      enddo

      nord=1
 91   l=nord+1
      kk=l+1
      do i=1,l
        do j=1,l
          ik=j-1+i
          a(i,j)=sumx(ik)
        enddo
        a(i,kk)=sumy(i)
      enddo

      do i=1,l
        a(kk,i)=-1.
        kkk=i+1
        do j=kkk,kk
          a(kk,j)=0.0
        enddo

        c=1.0/a(1,i)
        do ii=2,kk
          do j=kkk,kk
            a(ii,j)=a(ii,j)-a(1,j)*a(ii,i)*c
          enddo
        enddo

        do ii=1,l
          do j=kkk,kk
            a(ii,j)=a(ii+1,j)
          enddo
        enddo
      enddo

      s2=0.0
      do j=1,n
        s1=a(1,kk)
        do i=1,nord
          s1=s1+a(i+1,kk)*x(j)**i
        enddo
        s2=s2+(s1-y(j))*(s1-y(j))
      enddo

      b=n-l
      if(s2 > 0.0001) then
        s2=sqrt(s2/b)
      endif

      do i=1,l
        j=i-1
        cof(i)=a(i,kk)
      enddo

      if(nord<last) then
        if(s2>tol) then
          nord=nord+1
          j=2*nord
          sumx(j)=0.0
          sumx(j+1)=0.0
          sumy(nord+1)=0.0
                                                                
          do i=1,n                                              
            sumx(j)=sumx(j)+x(i)**(j-1)                         
            sumx(j+1)=sumx(j+1)+x(i)**j                         
            sumy(nord+1)=sumy(nord+1)+y(i)*x(i)**nord           
          enddo                                                 
                                                                
          goto 91
        endif
      endif

      return
      end subroutine polfit
