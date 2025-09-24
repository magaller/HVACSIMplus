!***********************************************************************

      real function bakdif(s,id,isblk)

!----------------------------------------------------------------------
!
!     BAKDIF : Calculates derivatives using backward difference formulas.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated  : Jan. 30, 1987 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     subroutines called: none
!
!***********************************************************************

      use modsim_head
      implicit none
      integer               :: isblk,id,j,kp1,i
      real                  :: s

      j=ip(isblk,1)
      spast(id,j)=s-stold(idevar(id))
      bakdif=0.
      kp1=korder(isblk)+1
      j=ip(isblk,kp1)
      ah(isblk,j)=0.

!     Ref. 2, eqn. 2.

      do i=1,kp1
        bakdif=bakdif-ah(isblk,i)*spast(id,i)
      enddo
      if(step(isblk)<tmin) then
        step(isblk)=tmin
      endif
      bakdif=bakdif/step(isblk)

      return
      end function bakdif

!***********************************************************************

      real function caln(isblk,irej)

! ----------------------------------------------------------------------
!
!     CALN : Calculate minimum time step ratio.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 7, 1984
!
!       Updated on :  June 29, 1989 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     subroutines called: predik
!
!***********************************************************************

      use modsim_head
      implicit none
      integer              :: irej,isblk,i,id,ivar
      real                 :: e,ec,fn,predik

      irej=0
      caln=-1.0
      do i=1,ndsb(isblk)
        id=indesb(isblk,i)
        ivar=idevar(id)
        if(icheck(ivar)/=2) then

!       Ref. 2, eqn. 5b.

          e=dnk(isblk)*abs(state(ivar)-predik(id,isblk))

!       Ref. 2, Appendix II, step BDF-4(b).

          if(e> 1.0e-20) then
            ec=(rtolx*abs(state(ivar))+atolx)*step(isblk)/ttime
            if(e>3.*ec) irej=1
            fn=ec/e
            if(caln<=0.) caln=fn
            caln=min(caln,fn)
          endif
        endif
      enddo
      if(caln<=0.) caln=2.

      return
      end function caln

!***********************************************************************

      integer function ecntrl(isblk)

! ----------------------------------------------------------------------
!
!     ECNTRL : Calculate truncation errors, and determines
!              time step and integration order changes.
!              See ref. 2, appendix II, step BDF-4.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated on: Oct. 25, 1984 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     subroutines called: caln,norder,saveco
!
!***********************************************************************

      use modsim_head
      implicit none
      integer            :: isblk,k,irej,is,ia0,ia1,ika,i,j
      real               :: nk,nkm1,nkp1,nmax
      real               :: fk,caln

      ecntrl=0
      k=korder(isblk)
      fk=1./real(k)

!     Calculate the minimum time step ratio for the present order.
!     If korder+1 steps have been taken since the last attempt to
!     change order, calculate the ratio for orders korder+1 and
!     korder-1. Otherwise increment the step counter, change the time
!     step, and return.

      nk=caln(isblk,irej)
      if(irej/=0) then
        ecntrl=1
        step(isblk)=step(isblk)*nk
        if(step(isblk)<tmin) then
          step(isblk)=tmin
        endif
        if(kstep(isblk)<k+1) then
          kstep(isblk)=kstep(isblk)+1
        endif
        return
      endif

      nk=caln(isblk,irej)**fk
      if(kstep(isblk)<k+1) then
        if((step(isblk)>2.*tmin).and.(nk>2.)) nk=2.
        step(isblk)=step(isblk)*nk
        if(step(isblk)<tmin) then
          step(isblk)=tmin
        endif
        kstep(isblk)=kstep(isblk)+1
        return
      endif

!     Save the present coefficient vectors, a and g.

      nkm1=0.
      nkp1=0.
      is=0
      call saveco(is,k,isblk)

!     If the present order is one, the order cannot be decreased.

      if(korder(isblk)/=1) then

!     Decrease the order by one.

        is=-1
        call norder(is,isblk)
        fk=1./real(korder(isblk))

!     Calculate the time step ratio for order korder-1.

        nkm1=caln(isblk,irej)**fk

!     If the present order is six, the order cannot be increased.

        if(k<6) then

!     Replace the a and g vectors with those saved above.

          is=1
          call saveco(is,k,isblk)
        endif
        korder(isblk)=k
      endif
      if(k/=6) then

!     Calculate the time step ratio for order korder+1.

        is=1
        call norder(is,isblk)
        fk=1./real(korder(isblk))
        nkp1=caln(isblk,irej)**fk
      endif

!     Determine the order for which the time step ratio is largest and
!     replace the present order.

      nmax=max(nk,nkm1,nkp1)
      korder(isblk)=k
      if(nk>=nmax) then
        is=1
        call saveco(is,k,isblk)
      elseif(nkm1>=nmax) then
        is=1
        call saveco(is,k,isblk)
        is=-1
        call norder(is,isblk)
      else
        korder(isblk)=korder(isblk)+1
      endif

      if((step(isblk)>2.*tmin).and.(nmax>2.)) nmax=2.
      step(isblk)=step(isblk)*nmax
      if(step(isblk)<tmin) then
        step(isblk)=tmin
      endif
      kstep(isblk)=0
      ia0=ipa(isblk,1)
      ia1=ipa(isblk,2)
      ika=ipa(isblk,korder(isblk)+2)

!     Calculate new values for f(n,k) and d(n,k).

      dnk(isblk)=am(isblk,ia0,ia1)/am(isblk,ia0,ika)
      fnk(isblk)=1.
      do i=0,korder(isblk)
        j=ipa(isblk,i+2)
        fnk(isblk)=fnk(isblk)*am(isblk,ia0,j)
      enddo

      return
      end function ecntrl

!***********************************************************************

      subroutine iperm(isblk)

! ----------------------------------------------------------------------
!
!     IPERM : Permute vectors IP and IPA.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated on : Oct. 25, 1984 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     subroutines called: none
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                  :: isblk,kord,ik,i,j

      kord=korder(isblk)
      ik=ip(isblk,kord+1)
      do i=1,kord
        j=kord+2-i
        ip(isblk,j)=ip(isblk,j-1)
      enddo
      ip(isblk,1)=ik
      ik=ipa(isblk,8)
      do i=1,7
        j=9-i
        ipa(isblk,j)=ipa(isblk,j-1)
      enddo
      ipa(isblk,1)=ik

      return
      end subroutine iperm

!***********************************************************************

      subroutine norder(is,isblk)

! ----------------------------------------------------------------------
!
!     NORDER : Increases or decreases the order of integration.
!       NOTE: when the order is decreased, the permutation vector
!       IP is changed as well as the ordering of arrays SPAST, A and G.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 8, 1984
!
!       Updated on : Oct. 25, 1984 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     subroutines called: none
!
!***********************************************************************

      use modsim_head
      implicit none
      integer             :: isblk,is,kord,ia0,kp1,ika,i,j,ja,l,ia1,kk,id
      real                :: gg,s

      kord=korder(isblk)
      ia0=ipa(isblk,1)
      kp1=kord+1

!     IS is a switch. is>0 - increase order by one; IS<0 - decrease
!     order by one.

      if(is<0) goto 60

!     Increase order by one.

      ip(isblk,kp1+1)=kp1+1
      ika=ipa(isblk,kp1+2)

!     Ref. 2, eqn. 26a.

      do i=1,kp1
        j=ip(isblk,i)
        ja=ipa(isblk,i+1)
        g(isblk,j)=g(isblk,j)*am(isblk,ia0,ika)/am(isblk,ja,ika)
      enddo
      g(isblk,kp1+1)=1.
      do i=1,kp1
        j=ip(isblk,i)
        g(isblk,kp1+1)=g(isblk,kp1+1)-g(isblk,j)
      enddo

!     The gh vector is required for the truncation error calculation.

      gg=0.
      do i=1,kp1
        l=kp1+1-i
        j=ip(isblk,l+1)
        gg=gg-g(isblk,j)
        gh(isblk,j)=gg
      enddo
      ika=ipa(isblk,kp1+1)

!     Ref. 2, eqn. 27a.

      do i=1,kp1
        j=ip(isblk,i)
        ja=ipa(isblk,i)
        a(isblk,j)=a(isblk,j)*am(isblk,ia0,ika)/am(isblk,ja,ika)
      enddo
      ia1=ipa(isblk,2)
      j=ip(isblk,1)
      a(isblk,j)=a(isblk,j)-am(isblk,ia0,ia1)/am(isblk,ia0,ika)
      a(isblk,kp1+1)=0.
      do i=1,kp1
        l=kp1+1-i
        j=ip(isblk,l)
        a(isblk,kp1+1)=a(isblk,kp1+1)-a(isblk,j)
      enddo
      kord=kord+1
      korder(isblk)=kord
      ika=ipa(isblk,kord+2)

!     d(n,k+1) is required for the truncation error calculation.

      dnk(isblk)=am(isblk,ia0,ia1)/am(isblk,ia0,ika)
      return

!     Decrease order by one.

60    ika=ipa(isblk,kp1+1)
      ia1=ipa(isblk,2)

!     Ref. 2, eqn. 26b.

      do i=1,kord
        j=ip(isblk,i)
        ja=ipa(isblk,i+1)
        g(isblk,j)=g(isblk,j)*am(isblk,ja,ika)/am(isblk,ia0,ika)
      enddo
      ika=ipa(isblk,kp1)

!     Ref. 2, eqn. 27b.

      do i=2,kord
        j=ip(isblk,i)
        ja=ipa(isblk,i)
        a(isblk,j)=a(isblk,j)*am(isblk,ja,ika)/am(isblk,ia0,ika)
      enddo
      j=ip(isblk,1)
      a(isblk,j)=a(isblk,j)+am(isblk,ia0,ia1)/am(isblk,ia0,ika)

!     Search the ip vector. the first norder positions must contain the
!     first norder values of a, g, and spast. This is not necessarily
!     the case when the order is decreased.

      do i=1,kp1
        j=i
        if(ip(isblk,j)==kp1) exit
      enddo

!     Switch the appropriate values of a, g, and spast.

      kk=ip(isblk,kp1)
      do i=1,ndsb(isblk)
        id=indesb(isblk,i)
        s=spast(id,kp1)
        spast(id,kp1)=spast(id,kk)
        spast(id,kk)=s
      enddo
      s=g(isblk,kp1)
      g(isblk,kp1)=g(isblk,kk)
      g(isblk,kk)=s
      s=a(isblk,kp1)
      a(isblk,kp1)=a(isblk,kk)
      a(isblk,kk)=s
      ip(isblk,j)=kk
      ip(isblk,kp1)=kp1
      kord=kord-1
      korder(isblk)=kord

!     The gh vector is required for the truncation error calculation.

      gg=0.
      do i=1,kord
        l=kord+1-i
        j=ip(isblk,l+1)
        gg=gg-g(isblk,j)
        gh(isblk,j)=gg
      enddo
      ika=ipa(isblk,kord+2)

!     d(n,k-1) is required for the truncation error calcualtion.

      dnk(isblk)=am(isblk,ia0,ia1)/am(isblk,ia0,ika)

      return
      end subroutine norder

!***********************************************************************

      real function predik(id,isblk)

! ----------------------------------------------------------------------
!
!     PREDIK : Calculate predicted values for the next time step.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated on: Oct. 25, 1984 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     subroutines called: none
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                :: isblk,id,i

      spast(id,ip(isblk,1))=0.
      predik=spred(id)

!     Ref. 2, eqn. 3.

      do i=1,korder(isblk)+1
        predik=predik+gh(isblk,i)*spast(id,i)
      enddo

      return
      end function predik

!***********************************************************************

      subroutine reset(is,isblk)

! ----------------------------------------------------------------------
!
!     RESET : Reset the differential equation integrator.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated on : Oct. 25, 1984 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     subroutines called : iperm
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                :: isblk,is,i,id,i1,i2,j

      if(is==0) then
        do i=1,ndsb(isblk)
          id=indesb(isblk,i)
          spast(id,1)=0.
        enddo
      else
        i1=ipa(isblk,2)
        i2=ipa(isblk,3)
        j=ip(isblk,2)
        do i=1,ndsb(isblk)
          id=indesb(isblk,i)
          spast(id,1)=spast(id,j)*step(isblk)/am(isblk,i1,i2)
        enddo
      endif
      korder(isblk)=1
      kstep(isblk)=0
      do i=1,8
        do j=i,8
        am(isblk,i,j)=step(isblk)*real(j-i)
        am(isblk,j,i)=-am(isblk,i,j)
        enddo
      enddo
      g(isblk,1)=2.
      g(isblk,2)=-1.
      a(isblk,1)=-1.
      a(isblk,2)=1.
      dnk(isblk)=0.5
      fnk(isblk)=2.*step(isblk)**2
      do i=1,7
        ip(isblk,i)=i
        ipa(isblk,i)=i
      enddo
      ipa(isblk,8)=8
      call iperm(isblk)

      return
      end subroutine reset

!***********************************************************************

      subroutine saveco(is,k,isblk)

! ----------------------------------------------------------------------
!
!     SAVECO : Save or replace coefficients.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated on : May 27, 1992 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     subroutines called: none
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                 :: isblk,k,is,i
      real,dimension(7)       :: gs,as1
      integer,dimension(7)    :: ips
      save ips

      if(is<=0) then
        do i=1,k+1
          gs(i)=g(isblk,i)
          as1(i)=a(isblk,i)
          ips(i)=ip(isblk,i)
        enddo
      else
        do i=1,k+1
          g(isblk,ip(isblk,i))=gs(ips(i))
          a(isblk,ip(isblk,i))=as1(ips(i))
        enddo
      endif

      return
      end subroutine saveco

!***********************************************************************

      subroutine update(isblk)

! ----------------------------------------------------------------------
!
!     UPDATE :
!       Update coefficients for functions BAKDIF and PREDIK.
!       Assume that permutation vectors IP and IPA are set for new time
!       step.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated : Jan. 30, 1987 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     subroutines called: none
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                    :: isblk,kp1,i1,ia0,ia1,i,l,j,&
                                    ja,k,jnu,ikp1,jp1,jap1,iap1
      real                       :: s,ss,g1,gg,t1
      integer,dimension(maxsbk)  :: istep
      data istep/maxsbk*0/

!      DATA ISTEP/MAXSBK*0/

      istep(isblk)=istep(isblk)+1
      kp1=korder(isblk)+1
      i1=ip(isblk,1)
      ia0=ipa(isblk,1)
      ia1=ipa(isblk,2)
      if(step(isblk)<tmin) then
        step(isblk)=tmin
      endif

!     Update the a matrix. See Ref. 2, Appendix II, step BDF-1.

      do i=1,6
        l=9-i
        j=ipa(isblk,l)
        am(isblk,ia0,j)=step(isblk)+am(isblk,ia1,j)
      enddo
      am(isblk,ia0,ia1)=step(isblk)
      do i=2,8
        j=ipa(isblk,i)
        am(isblk,j,ia0)=-am(isblk,ia0,j)
      enddo

!     Calculate t1, f(n+1,k), g1, and d(n+1,k). See Ref. 2, Appendix II,
!     step BDF-2.

      t1=dnk(isblk)*fnk(isblk)
      fnk(isblk)=1.
      do i=1,kp1
        j=ipa(isblk,i+1)
        fnk(isblk)=fnk(isblk)*am(isblk,ia0,j)
      enddo
      g1=-fnk(isblk)/t1
      j=ipa(isblk,kp1+1)
      dnk(isblk)=am(isblk,ia0,ia1)/am(isblk,ia0,j)

!     Every ten steps the g vector is calculated exactly using eqn. 21
!     of Ref. 2. In this case it is not necessary to permute vectors
!     a and g or use step BDF-2(d).

      if(mod(istep(isblk),10)==0) then
        do i=2,kp1
          j=ip(isblk,i)
          ja=ipa(isblk,i+1)
          s=am(isblk,ia0,ja)
          do k=1,kp1
            if(i==k) cycle
            jnu=ipa(isblk,k+1)
            s=s*am(isblk,ja,jnu)
          enddo
          g(isblk,j)=fnk(isblk)/s
        enddo
      else

!       Permute vectors a and g in inverse order of the permutation
!       vector ip so that they are in proper order for the next time
!       step.

        s=a(isblk,i1)
        ss=g(isblk,i1)
        do i=2,kp1
          k=ip(isblk,i)
          l=ip(isblk,i-1)
          a(isblk,l)=a(isblk,k)
          g(isblk,l)=g(isblk,k)
        enddo
        ikp1=ip(isblk,kp1)
        a(isblk,ikp1)=s
        g(isblk,ikp1)=ss

!       Update vector g. see Ref. 2, Appendix II, step BDF-2(d,e). The
!       gh vector is calculated at the same time.

        do i=1,korder(isblk)
          jp1=ip(isblk,i+1)
          jap1=ipa(isblk,i+2)
          g(isblk,jp1)=g1*a(isblk,jp1)*am(isblk,ia1,jap1)/&
                       am(isblk,ia0,jap1)
        enddo
      endif

      gg=0.
      do i=1,korder(isblk)
        j=korder(isblk)+1-i
        l=ip(isblk,j+1)
        gg=gg-g(isblk,l)
        gh(isblk,l)=gg
      enddo
      g(isblk,i1)=1.+gg
      iap1=ipa(isblk,kp1+1)

!     Update vector a. See Ref. 2, Appendix II, step BDF-3. The vector
!     ah is calculated at the same time.

      do i=1,korder(isblk)
        j=ip(isblk,i)
        jp1=ip(isblk,i+1)
        jap1=ipa(isblk,i+1)
        a(isblk,jp1)=dnk(isblk)*g(isblk,j)*am(isblk,jap1,iap1)/&
                     am(isblk,ia0,jap1)
      enddo
      gg=0.
      do i=1,korder(isblk)
        k=korder(isblk)+1-i
        j=ip(isblk,k)
        jp1=ip(isblk,k+1)
        gg=gg-a(isblk,jp1)
        ah(isblk,j)=gg
      enddo
      a(isblk,i1)=gg

      return
      end subroutine update

