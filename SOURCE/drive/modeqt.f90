!***********************************************************************

      subroutine block(isblk,iblk,iset)

! ----------------------------------------------------------------------
!
!     BLOCK : Calculate a part of a new state vector.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Calculate part of a new state vector by calling a block. If a
!       block  contains more than one unit, the variables specified
!       by ISOLVE  are determined by solving the set of simultaneous
!       equations defined by the chosen variables.
!       Called by subroutines SUPERB and SUPFNC.
!       When a Jacobian matrix for a superblock is being calculated,
!       input  parameter ISET is set equal to -1 and equation solver
!       SNSQ2 is used.
!       SNSQ2 is an abreviated version of SNSQ which requires much less
!       execution time.
!
!       Modified into FORTRAN77 by Cheol Park, September 6, 1984
!
!       Updated on : December 21, 1984 Dan Clark
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     subroutines called: inputs,output,select,predik,snsq1,snsq2,fnc
!
!***********************************************************************

      use modsim_head
	  use nist_solver
      implicit none
      external fnc
      integer               :: iset,iblk,isblk,neqs,i,k,nprint,&
                               ivar,kk,id,maxfev,njev,nfev,info,j,l,&
                               iu,mode
      real                  :: jac,predik,xtol
      integer,parameter     :: neqib=mseqib
      integer,parameter     :: mlb=neqib-1,mub=mlb,lrb=neqib*(neqib+1)/2

      character(len=4),dimension(mseqib):: label
      integer,dimension(mseqib)         :: number
      real,dimension(minoiu)            :: xin,out
      real,dimension(mseqib)            :: s,fvec,diag,qtf,wa1,wa2,wa3,wa4
      real,dimension(maxdeq)            :: dout
      real,dimension(mseqib,mseqib)     :: fjac
      real,dimension(lrb)               :: r
!	  NITSOL
      integer nitinfo(6), ipar(1), iterm 
      real nitftol, stptol, rwork(100000), rpar(1) 
      external jacv, ddot, dnrm2

!     If no variables are solved simultaneously, the unit(s) are called
!     and snsq is not called.

      neqs=nsolve(iblk)
      if(neqs<=0) then
        do i=1,nunits(iblk)
          iu=iblock(iblk,i)
          call inputs(iu,xin)
          call select(iu,xin,out)
          call output(iu,out)
        enddo
        return
      endif

!     More than one unit in the block. Initial guess input to snsq is
!     vector s. If a variable is determined by a de, the initial guess
!     is calculated by function predik.

      loop_30: do k=1,neqs
        ivar=isolve(iblk,k)
        if(ndsb(isblk)>0) then
          do kk=1,ndsb(isblk)
            id=indesb(isblk,kk)
            if(idevar(id)==ivar) then
              s(k)=predik(id,isblk)
              if(s(k)<state(ivar)) then
                isign(ivar)=-1
              else
                isign(ivar)=1
              endif
              cycle loop_30
            endif
          enddo
        endif
        s(k)=state(ivar)
      enddo loop_30

!     Write diagnostic information to ifile3, if appropriate

      if((iprint/=0).and.(time>=tprton).and.(time<=tprtof)) then
        nprint=iprint
        write(ifile3,1000) iblk,time,neqs
        call vlab(neqs,maxblk,mseqib,isolve,iblk,label,number)
        write(ifile3,2000) (label(i),number(i),i=1,neqs)
      else
        nprint=0
      endif

!     Call snsq.

      maxfev=200*(neqs+1)
      mode=2
      do k=1,neqs
        diag(k)=1./(abs(s(k))+atolx)
      enddo

!     If a jacobian matrix is being calculated, call snsq2. Otherwise,
!     call snsq1.

   
      xtol = xtola                              ! <<< 11/28/06
      if(iset==-1) then
        call snsq2(fnc,fnc,2,neqs,s,fvec,fjac,neqib,xtol,maxfev,mlb,&
                   mub,0.,diag,mode,1.,nprint,info,nfev,njev,r,lrb,&
                   qtf,wa1,wa2,wa3,wa4,iblk,isblk)
      else
	  solution_method: select case (solver)
!--------------------------------------------------------------------------
		case default
		call snsq1(fnc,fnc,2,neqs,s,fvec,fjac,neqib,xtol,maxfev,mlb,&
                   mub,0.,diag,mode,1.,nprint,info,nfev,njev,r,lrb,&
                   qtf,wa1,wa2,wa3,wa4,iblk,isblk)
!--------------------------------------------------------------------------
		case(1)
			nitftol=user_ftol	
			stptol=xtol		
			call nitsol(neqs,s,fnc,jacv,nitftol,stptol,nitinput,&
     		   nitinfo,rwork,rpar,ipar,iterm,ddot,dnrm2,&
     		   iblk,isblk)
			info=1	!ZC
!--------------------------------------------------------------------------
	  end select solution_method			
      endif
!--------------------------------------------------------------------------

      if(info<1) stop

!     Set flag for "last iteration" to enable error messages from types
      lastit=1

!     Calculate outputs from the block using the inputs determined by
!     snsq.

      big_loop: do j=1,nunits(iblk)
        iu=iblock(iblk,j)
        call inputs(iu,xin)
        do kk=1,neqs
          do l=1,nin(iu)
            if(in(iu,l)==isolve(iblk,kk)) xin(l)=s(kk)
          enddo

!         If a unit contains differential equations, save the solved
!         variable for later use.

          if(nde(iu)>0) then
            do i=1,nde(iu)
              id=inde(iu,i)
              if(idevar(id)==isolve(iblk,kk)) dout(id)=s(kk)
            enddo
          endif
        enddo

!       Call the type subroutine and store the outputs.

        call select(iu,xin,out)
        call output(iu,out)

!       If the unit contains differential equations, replace the
!       derivative (which is returned by the type subroutine) with the
!       solved variable which was saved above.

        if(nde(iu)>0) then
          do i=1,nde(iu)
            id=inde(iu,i)
            tstate(idevar(id))=dout(id)
          enddo
        endif
      enddo big_loop

!     Reset "last iteration" flag

      lastit=0

1000  format(/1x,'block ',i3,', time =',f10.3,':  number of active ',&
             'simultaneous eqs. =',i3)
2000  format(5(5x,a4,i2,5x))

      return
      end subroutine block

!***********************************************************************

      subroutine fnc(n,x,fvec,iflag,iblk,isblk)

! ----------------------------------------------------------------------
!
!     FNC : Calculate the residual functions for block IBLK.
!       Called by nonlinear equation solver SNSQ.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated on : December 21, 1984 D.C.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     subroutines called: inputs,select,bakdif
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                 :: isblk,iblk,iflag,ified,i,j,&
                                 nlout,l,k,kk,id,n,iu
      real                    :: fmax,bakdif
      real,dimension(n)       :: x,fvec
      real,dimension(minoiu)  :: xin,out,lout,kout

!     Write intermediate solutions.

      if(iflag==0) then
        fmax=0.
        do i=1,n
          fmax=max(fmax,abs(fvec(i)))
        enddo
        write(ifile3,1000) iblk,fmax
        write(ifile3,2000) (fvec(i),i=1,n)
        write(ifile3,3000) 
        write(ifile3,2000) (x(i),i=1,n)
        return
      endif

!     Begin loop to call units.

      big_loop: do j=1,nunits(iblk)
      iu=iblock(iblk,j)
      nlout=0

!     Check each unit for outputs which are solved simultaneously.

      loop_30: do l=1,nout(iu)
        do k=1,n
          if(iout(iu,l)==isolve(iblk,k)) then
            nlout=nlout+1
            lout(nlout)=l
            kout(nlout)=k
            cycle loop_30
          endif
        enddo
      enddo loop_30

!     Skip unit if no outputs are solved simultaneously.

      if(nlout>0) then

!       Set up the input vector for the unit containing equation k.

        call inputs(iu,xin)
        do kk=1,nsolve(iblk)
          do l=1,nin(iu)
          if(in(iu,l)==isolve(iblk,kk)) xin(l)=x(kk)
          enddo
        enddo

!       Call the type subroutine and calculate the residual.

        call select(iu,xin,out)

!       If the equation labeled by k is a differential equation, the
!       residual function is the difference between the derivative
!       calculated by the type subroutine and the backward difference
!       calculated by function bakdif.

        loop_60: do kk=1,nlout
          l=lout(kk)
          k=kout(kk)
          if(nde(iu)>0) then
            do i=1,nde(iu)
              id=inde(iu,i)
              if(idevar(id)==iout(iu,l)) then
                fvec(k)=(out(l)-bakdif(x(k),id,isblk))*tstep
                cycle loop_60
              endif
            enddo
          endif
          fvec(k)=out(l)-x(k)
        enddo loop_60
      endif
      enddo big_loop

1000  format(1x,'snsq1: intermediate solution, block',i4/1x,'maximum',&
        ' residual=',e9.2,'     residuals:')
2000  format(1x,5e15.6)
3000  format(1x,'intermediate solution:')

      return
      end subroutine fnc

!***********************************************************************

      subroutine superb(isblk)

! ----------------------------------------------------------------------
!
!     SUPERB : Solve a superblock.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, September 6, 1984
!
!       Updated on : December 6, 1984 Dan Clark
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     SUBROUTINES CALLED: SNSQ,BLOCK,SUPFNC
!
!***********************************************************************

      use modsim_head
      implicit none
      external  supfnc
      integer               :: iset,iblk,isblk,neqs,i,k,&
                               ivar,kk,id,maxfev,njev,nfev,info,j,l,&
                               nprint,mode
      real                  :: jac,epsfcn,xtol
      integer,parameter     :: neqis=mseqis
      integer,parameter     :: mls=neqis-1,mus=mls,lrs=neqis*(neqis+1)/2
      character(len=4),dimension(mseqis):: label
      integer,dimension(mseqis)         :: number
      real,dimension(mseqis)            :: s,fvec,diag,qtf,wa1,wa2,wa3,wa4
      real,dimension(mseqis,mseqis)     :: fjac
      real,dimension(lrs)               :: r

!     If no variables are solved simultaneously, call all blocks and
!     return.

      neqs=nssolv(isblk)
      if(neqs>0) then

!       Set up initial guess for snsq.

        do i=1,neqs
          ivar=issolv(isblk,i)
          s(i)=state(ivar)
        enddo

!     Write diagnostic information to ifile3, if appropriate

        if((iprint/=0).and.(time>=tprton).and.(time<=tprtof)) then
          nprint=iprint
          write(ifile3,1000) isblk,time,neqs
          call vlab(neqs,maxsbk,mseqis,issolv,isblk,label,number)
          write(ifile3,2000) (label(i),number(i),i=1,neqs)
        else
          nprint=0
        endif

!       Call snsq.

        xtol = xtola                     ! <<< 11/28/06
        maxfev=200*(neqs+1)
        epsfcn=10.*xtol
        mode=2
        do k=1,neqs
          diag(k)=1./(abs(s(k))+atolx)
        enddo
        
        call snsq(supfnc,supfnc,2,neqs,s,fvec,fjac,neqis,xtol,maxfev,mls,&
              mus,epsfcn,diag,mode,1.,nprint,info,nfev,njev,r,lrs,qtf,&
              wa1,wa2,wa3,wa4,isblk)
        if(info<1) stop

!       Save output from snsq.

        do i=1,neqs
          ivar=issolv(isblk,i)
          tstate(ivar)=s(i)
          state(ivar)=s(i)
        enddo
      endif

!     Call all blocks.

      iset=0
      do i=1,nsuper(isblk)
        iblk=isuper(isblk,i)
        if(ibstat(iblk)/=0) then
          call block(isblk,iblk,iset)
        endif
      enddo

1000  format(/1x,'superblock ',i3,', time =',f10.3,':  # of active ',&
             'simultaneous eqs. =',i3)
2000  format(5(5x,a4,i2,5x))

      return
      end subroutine superb

!***********************************************************************

      subroutine supfnc(n,x,fvec,iflag,isblk,iset)

! ----------------------------------------------------------------------
!
!     SUPFNC : Calculate the residual functions for superblock ISBLK.
!              Called by subroutine SNSQ.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated on : Oct. 25, 1984 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     subroutines called: block
!
!***********************************************************************

      use modsim_head
      implicit none
      integer               :: iset,isblk,i,iblk,j,k,iflag,n
      real                  :: fmax
      real,dimension(n)     :: x,fvec

!     Write intermediate solution.

      if(iflag==0) then
        fmax=0.
        do i=1,n
          fmax=max(fmax,abs(fvec(i)))
        enddo
        write(ifile3,1000) isblk,fmax
        write(ifile3,2000) (fvec(i),i=1,n)
        write(ifile3,3000) 
        write(ifile3,2000) (x(i),i=1,n)
        return
      endif

!     Set state variables passed from snsq.

      do i=1,n
        tstate(issolv(isblk,i))=x(i)
        state(issolv(isblk,i))=x(i)
      enddo

!     Call a block only if an output is solved simultaneously.

      loop_40: do i=1,nsuper(isblk)
        iblk=isuper(isblk,i)
        if(ibstat(iblk)/=0) then
          do j=1,noutb(iblk)
            do k=1,n
              if(ioutb(iblk,j)==issolv(isblk,k)) then
                call block(isblk,iblk,iset)
                cycle loop_40     
              endif
            enddo
          enddo
        endif
      enddo loop_40

!     Calculate residual functions.

      do i=1,n
        fvec(i)=tstate(issolv(isblk,i))-x(i)
      enddo

!     If iset=1, replace unfrozen state variables with block outputs.

      if(iset==1) then
        do i=1,nsuper(isblk)
          iblk=isuper(isblk,i)
          do j=1,noutb(iblk)
            if(ioutb(iblk,j)>0) then
              if(abs(icheck(ioutb(iblk,j)))/=2) then
                state(ioutb(iblk,j))=tstate(ioutb(iblk,j))
              endif
            endif
          enddo
        enddo
      endif

1000  format(1x,'snsq:  intermediate solution, superblock',i4/&
             1x,'maximum residual =',e9.2,'     residuals:')
2000  format(1x,5e15.6)
3000  format(1x,'intermediate solution:')

      return
      end subroutine supfnc

