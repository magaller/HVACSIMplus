!***********************************************************************

      subroutine asembl

!-----------------------------------------------------------------------
!
!     ASEMBL : Assemble block input and output information.
!
!       C. R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated : April 13, 1989  C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                  :: i,j,iblk,i1,k,i2,i3,isblk,iu

!     Assemble block input and output arrays.

      do i=1,nsblok
        do j=1,nsuper(i)
          iblk=isuper(i,j)                               
          noutb(iblk)=0                                  
          i1=0                                           
          do k=1,nunits(iblk)                            
            iu=iblock(iblk,k)                            
            do i2=1,nout(iu)                             
              if(iout(iu,i2)>0) then
                noutb(iblk)=noutb(iblk)+1                
                i1=i1+1                                  
                ioutb(iblk,i1)=iout(iu,i2)               
              endif                                      
            enddo                                        
          enddo                                          
          i1=0                                           
          ninb(iblk)=0                                   
          do k=1,nunits(iblk)                            
            iu=iblock(iblk,k)                            
            loop_50: do i2=1,nin(iu)
              if(i1>0) then
                do i3=1,i1
                  if(inb(iblk,i3)==in(iu,i2)) cycle loop_50
                enddo
              endif
              i1=i1+1                                    
              ninb(iblk)=ninb(iblk)+1                    
              inb(iblk,i1)=in(iu,i2)                     
            enddo loop_50
          enddo
        enddo
      enddo

!     Assemble vectors identifying DE's solved in each superblock.

      do isblk=1,nsblok
        ndsb(isblk)=0
        i1=0
        do i=1,nsuper(isblk)
          iblk=isuper(isblk,i)
          do j=1,nunits(iblk)
            iu=iblock(iblk,j)
            if(nde(iu)>0) then
              ndsb(isblk)=ndsb(isblk)+nde(iu)
              do k=1,nde(iu)
                i1=i1+1
                indesb(isblk,i1)=inde(iu,k)
              enddo
            endif
          enddo
        enddo
      enddo

      return
      end subroutine asembl

!***********************************************************************

      subroutine bactiv(isblk)

!-----------------------------------------------------------------------
!
!     BACTIV : Determine block activity status.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated : April 13, 1989 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                    :: isblk,ibk,iblk,j,ivar,k,i
      logical                    :: act

!     Begin loops to check for superblock equations.

      loop_30: do ibk=1,nsuper(isblk)
        iblk=isuper(isblk,ibk)
        ibstat(iblk)=-1

!       Check the block for superblock equations.
!       If superblock isblk has no simultaneous equations, continue.

        if(nssolv(isblk)>0) then

!         Check the outputs of block iblk.
                                                                                
          do j=1,noutb(iblk)                                                    
            ivar=ioutb(iblk,j)                                                  
            if(ivar>0) then
            do k=1,nssolv(isblk)                                                
              if(ivar==issolv(isblk,k)) then
                                                                                
!             If any output of the block is solved by a superblock,
!             mark the block.
                                                                                
                ibstat(iblk)=2
                cycle loop_30
              endif
            enddo                                                               
            endif                                                               
          enddo                                                                 
        endif

!       Mark any blocks not previously marked.

        if(ibstat(iblk)<0) then
          if(nsolve(iblk)<=0) then
            ibstat(iblk)=1
          else
            ibstat(iblk)=2
          endif
        endif
      enddo loop_30

!     Begin loop to mark inactive blocks.

      loop_50: do ibk=1,nsuper(isblk)
        iblk=isuper(isblk,ibk)

!       If block iblk has been marked active with some variables unfrozen,
!       skip over it.

        if(ibstat(iblk)/=2) then

!       Check block connections to block iblk. if all inputs to block
!       iblk are frozen or are time independent boundary conditions,
!       the block is marked inactive.

          do i=1,ninb(iblk)
            ivar=inb(iblk,i)
            act= (icheck(ivar)==-3).or.&
                 (icheck(ivar)==-1).or.&
                 (icheck(ivar)== 0).or.&
                 (icheck(ivar)== 1).or.&
                 (icheck(ivar)== 3)
            if(act) cycle loop_50
          enddo
          ibstat(iblk)=0
        endif
      enddo loop_50

      return
      end subroutine bactiv

!***********************************************************************

      subroutine frzvar(isblk)

! ----------------------------------------------------------------------
!
!     FRZVAR : Freeze variables.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol park, August 24, 1984
!
!       Updated on : Oct. 25, 1984 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!       A variable is added or removed from arrays ISOLVE and ISSOLV
!       according to the change in the variable from one time step
!       to the next. If the change in the variable does not exceed
!       the error tolerance, the variable is frozen.
!
!       SUBROUTINES CALLED: NONE
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                    :: isblk,ibk,iblk,ndeb,nsolv,iv,ivar,&
                                    i,j,id,iu
      integer,dimension(muntib)  :: idrv

!     Begin loops to check variables solved in blocks.

      loop_60: do ibk=1,nsuper(isblk)
        iblk=isuper(isblk,ibk)
        ndeb=0
        nsolv=0
        isolve(iblk,1)=0
        loop_30: do iv=1,njsolv(iblk)
          ivar=jsolve(iblk,iv)

!         Count and mark the equations which are not frozen, excluding DE's.

          do i=1,nunits(iblk)
            iu=iblock(iblk,i)
            if(nde(iu)>0) then
              do j=1,nde(iu)
                id=inde(iu,j)
                if(idevar(id)==ivar) then
                  ndeb=ndeb+1
                  idrv(ndeb)=id
                  cycle loop_30
                endif
              enddo
            endif
          enddo

!         If a variable was unfrozen or is marked unfreezable, mark the
!         variable.

          if(icheck(ivar)==3) then
            icheck(ivar)=1
            nsolv=nsolv+1
          elseif(icheck(ivar)==0) then
            nsolv=nsolv+1

!         If the change in the variable does not exceed the error tolerance,
!         the variable is frozen.

          elseif(abs(state(ivar)-stold(ivar))<=&
                  (rtolx*abs(state(ivar))+atolx)) then
            icheck(ivar)=2
          else
            icheck(ivar)=1
            nsolv=nsolv+1
          endif
        enddo loop_30

!       Mark variables solved by de's. de's may not be frozen unless all
!       other  variables are frozen.

        if(ndeb>0) then
          do i=1,ndeb
            id=idrv(i)
            ivar=idevar(id)
            if(icheck(ivar)==3) icheck(ivar)=0
            if(nsolv>0) icheck(ivar)=0
            if(icheck(ivar)==0) cycle
            if(idechk(id)==1) icheck(ivar)=2
          enddo
        endif

!       Assemble the equation definition array (isolve).

        isolve(iblk,1)=0
        nsolv=0
        do iv=1,njsolv(iblk)
          ivar=jsolve(iblk,iv)
          if(icheck(ivar)==2) cycle
          nsolv=nsolv+1
          isolve(iblk,nsolv)=ivar
        enddo
        nsolve(iblk)=nsolv
      enddo loop_60

!     Begin the loops which check variables solved in superblocks

      nsolv=0
      if(njsslv(isblk)>0) then
        issolv(isblk,1)=0
        do iv=1,njsslv(isblk)
          ivar=jssolv(isblk,iv)

!     If the variable was unfrozen, it may not be refrozen until the
!     next time step.

          if(icheck(ivar)==3) then
            icheck(ivar)=1
            nsolv=nsolv+1
            issolv(isblk,nsolv)=ivar
          elseif(icheck(ivar)==0) then
            nsolv=nsolv+1
            issolv(isblk,nsolv)=ivar
          else
            icheck(ivar)=2

!           If the error tolerance is not exceeded, freeze the variable.

            if(abs(state(ivar)-stold(ivar))>=&
                rtolx*abs(state(ivar))+atolx) then
              icheck(ivar)=1
              nsolv=nsolv+1
              issolv(isblk,nsolv)=ivar
            endif
          endif
        enddo
      endif

      nssolv(isblk)=nsolv
      do ibk=1,nsuper(isblk)
        iblk=isuper(isblk,ibk)
        do iv=1,noutb(iblk)
          ivar=ioutb(iblk,iv)
          if(ivar>0.and.icheck(ivar)==-1) then
            icheck(ivar)=-2
            if(abs(state(ivar)-stold(ivar))>=&
                 rtolx*abs(state(ivar))+atolx) then
              icheck(ivar)=-1
            endif
          endif
        enddo
      enddo

      return
      end subroutine frzvar

!***********************************************************************

      subroutine inputs(iunit,x)

! ----------------------------------------------------------------------
!
!     INPUTS : Return inputs for unit IUNIT.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated on : Oct. 25, 1984 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!     SUBROUTINES CALLED: NONE
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                    :: iunit,i
      real,dimension(minoiu)     :: x

      do i=1,nin(iunit)
        x(i)=state(in(iunit,i))
      enddo

      return
      end subroutine inputs

!***********************************************************************

     integer function inscan(isblk)

! ----------------------------------------------------------------------
!
!     INSCAN :
!       Scan superblock inputs to detect changes larger than the error
!       tolerance.
!
!       C.R. Hill, National Bureau of Standards
!       April, 1983
!
!       Modified into FORTRAN77 by Cheol Park, September 5, 1984
!
!       Updated on : May 31, 1986 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                 :: isblk,i,ivar
      real                    :: err1

      inscan=0
      do i=1,ninsb(isblk)
        ivar=insb(isblk,i)
        err1=rtolx*abs(state(ivar))+atolx
        if(abs(state(ivar)-stold(ivar))>err1) then
          inscan=1
          return
        endif
      enddo

      return
      end function inscan

!***********************************************************************

      subroutine intliz(treset,ireset)

! ----------------------------------------------------------------------
!
!     INTLIZ : Initialize the simulation.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, September 20, 1984
!
!       Updatd on : May 31, 1986 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!       SUBROUTINES CALLED: RESET
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                  :: ireset,isblk,i,j,iblk,ivar,i1,k,id
      real                     :: treset            ! <<< 11/28/06

!     Initialize time counters.

      time=0.
      itime=0
      mtstep=1
      mtime=0
      ireset=0
      treset=0.
      tstep=tmin
      do isblk=1,nsblok
        ideflg(isblk)=0
        mstep(isblk)=1
        msbt(isblk)=0
        step(isblk)=tmin
        tnext(isblk)=treprt(isblk)
        sbtime(isblk)=0.
      enddo

!     Initialize arrays which define the equation sets.

      do i=1,nsblok
        nssolv(i)=njsslv(i)
        do j=1,njsslv(i)
          issolv(i,j)=jssolv(i,j)
        enddo
      enddo
      do i=1,nblock
        nsolve(i)=njsolv(i)
        do j=1,njsolv(i)
          isolve(i,j)=jsolve(i,j)
        enddo
      enddo

!     Initialize the variable status array.

      do i=1,nblock
        ibstat(i)=2
      enddo

      do i=1,nstate
        icheck(i)=-4
      enddo

      do i=1,nsblok
        do j=1,njsslv(i)
          icheck(jssolv(i,j))=0
        enddo
      enddo

      do i=1,nblock
        do j=1,njsolv(i)
          icheck(jsolve(i,j))=0
        enddo
      enddo

      do i=1,nbound
        icheck(ibound(i))=-3
      enddo

      do iblk=1,nblock
        do i=1,noutb(iblk)
          ivar=ioutb(iblk,i)
          if(ivar>0.and.icheck(ivar)==-4) then
            icheck(ivar)=-1
          endif
        enddo
      enddo

!     Assemble the superblock input array for input scanning.

      do isblk=1,nsblok
        i1=0
        do i=1,nsuper(isblk)
          iblk=isuper(isblk,i)
          loop_130: do j=1,ninb(iblk)
            ivar=inb(iblk,j)
            if(icheck(ivar)==-4) cycle loop_130
            do k=1,njsslv(isblk)
              if(ivar==jssolv(isblk,k)) cycle loop_130
            enddo
            do k=1,njsolv(iblk)
              if(ivar==jsolve(iblk,k)) cycle loop_130
            enddo
            i1=i1+1
            insb(isblk,i1)=ivar
          enddo loop_130
        enddo
        ninsb(isblk)=i1
      enddo

!     Initialize the de status vector and reset the integration
!     algorithm.

      if(nd<=0) return
      do id=1,nd
        idechk(id)=0
      enddo

      do isblk=1,nsblok
        call reset(0,isblk)
      enddo

      return
      end subroutine intliz

!***********************************************************************

      subroutine output(iunit,x)

! ----------------------------------------------------------------------
!
!     OUTPUT : Store outputs from unit IUNIT.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, August 24, 1984
!
!       Updated on : Oct. 25, 1984 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!       SUBROUTINES CALLED: NONE
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                     :: iunit,i
      real,dimension(minoiu)      :: x

      do i=1,nout(iunit)
        if(iout(iunit,i)<=0) cycle
        tstate(iout(iunit,i))=x(i)
      enddo

      return
      end subroutine output

!***********************************************************************

      subroutine restat(iset,iblk,isblk)

! ----------------------------------------------------------------------
!
!     RESTAT :
!       Reset the outputs of block IBLK when variables are unfrozen.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Modified into FORTRAN77 by Cheol Park, Oct. 2, 1984
!
!       Updated on : Oct. 25, 1984 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                      :: isblk,i,ivar,iblk,iset

      if(nsolve(iblk)<=0 ) return
      do i=1,nsolve(iblk)
        ivar=isolve(iblk,i)
        state(ivar)=stold(ivar)
      enddo
      if(iset==0 .or. nssolv(isblk)<=0) return
      do i=1,nssolv(isblk)
        ivar=issolv(isblk,i)
        state(ivar)=stold(ivar)
      enddo

      return
      end subroutine restat

!***********************************************************************

      integer function unfrez(isblk,iblk)

! ----------------------------------------------------------------------
!
!     UNFREZ : Unfreeze variables.
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Checks for frozen variables. If a variable is frozen and changes
!       by more than the error tolerance, it is unfrozen.
!
!       Modified into FORTRAN77 by Cheol Park, September 5, 1984
!
!       Updated on : December 6, 1984 C.P.
!       Updated to Fortran 90  December 27, 2006 C. Park
!
!       SUBROUTINES CALLED: NONE
!
!***********************************************************************

      use modsim_head
      implicit none
      integer                  :: iblk,isblk,ndeb,k,i,iv,id,ideflg1,iu 
      real                     :: hys=0.5,err1

      unfrez=0
      ndeb=0
      do k=1,nunits(iblk)
        iu=iblock(iblk,k)
        ndeb=ndeb+nde(iu)
      enddo
      if(nsolve(iblk)<ndeb) then
        ideflg1=1
      else
        ideflg1=0
      endif

!     Begin loop to check each block input.

20    loop_70: do i=1,noutb(iblk)
        iv=ioutb(iblk,i)
        if(iv>0.and.icheck(iv)==2) then
          err1=hys*(rtolx*abs(state(iv))+atolx)
          do k=1,njsolv(iblk)
            if(iv==jsolve(iblk,k)) goto 40
          enddo
          cycle loop_70

!         Check if the variable is solved by a de.

40        if(ndeb>0) then
            do id=1,nd
              if(iv==idevar(id)) then
                if(ideflg1==2.or.(ideflg1/=2.and.idechk(id)/=1)) then
                  nsolve(iblk)=nsolve(iblk)+1
                  isolve(iblk,nsolve(iblk))=iv
                  icheck(iv)=3
                  unfrez=1
                endif
                cycle loop_70
              endif
            enddo
          endif

!       If error tolerance is exceeded, unfreeze the variable.
!       If an algebraic variable unfreezes, unfreeze all de's.

          if(abs(tstate(iv)-stold(iv))>=err1) then
            nsolve(iblk)=nsolve(iblk)+1
            isolve(iblk,nsolve(iblk))=iv
            icheck(iv)=3
            unfrez=1
            if(ideflg1==1) then
              ideflg1=2
              goto 20
            endif
          endif
        endif
      enddo loop_70

!     Begin loops to determine variables solved in superblocks.

      if(njsslv(isblk)>0) then
        loop_100: do i=1,noutb(iblk)
          iv=ioutb(iblk,i)
          if(iv>0.and.icheck(iv)==2) then
            do k=1,njsslv(isblk)
              if(iv==jssolv(isblk,k)) goto 90
            enddo
            cycle loop_100

!           If error tolerance is exceeded, unfreeze the variable.

90          err1=hys*(rtolx*abs(state(iv))+atolx)
            if(abs(tstate(iv)-stold(iv))>=err1) then
              nssolv(isblk)=nssolv(isblk)+1
              issolv(isblk,nssolv(isblk))=iv
              icheck(iv)=3
              if(unfrez==0) then
                unfrez=-2
              elseif(unfrez==1) then
                unfrez=-1
              endif
            endif
          endif
        enddo loop_100
      endif

      do i=1,noutb(iblk)
        iv=ioutb(iblk,i)
        if(iv>0.and.icheck(iv)==-2) then
          err1=hys*(rtolx*abs(state(iv))+atolx)
          if(abs(tstate(iv)-stold(iv))>=err1) then
            icheck(iv)=-1
          endif
        endif
      enddo

      return
      end function unfrez

