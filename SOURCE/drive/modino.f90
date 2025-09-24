!***********************************************************************

      subroutine bounds(treset,ireset)

!   ----------------------------------------------------------------------
!
!       BOUNDS : Read time-dependent boundary conditions.
!
!         C.R. Hill, National Bureau of Standards
!         February, 1983
!
!         Modified into FORTRAN77 by Cheol Park, September 11, 1984
!
!         Modifications were made by Dan Clark, Dec. 1984.
!
!         Updated to Fortran 90  April 27, 2007 C. Park
!
!     Data is read from the file open to unit IFILE4. The data consists of
!     TIME and NBOUND data for at least four times. The state vector is
!     changed at every call to subroutine BOUNDS. Lagrangian interpolation
!     is used to calculate the new state variables (3rd order). If the
!     simulation time exceeds the maximum time in the input file, the
!     state vector is returned unchanged and a warning message is issued.
!
!     If a large change in a boundary variable occurs, the DE integration
!     algorithm must be reset.
!
!       subroutines called: none
!
!***********************************************************************

      use modsim_head
      implicit none

      integer                   :: ireset,ieof,nrold,i,it,nrec,&
                                   nbak,ib,nn,j
      integer,parameter         :: pd=selected_real_kind(15) ! double precision
      real(kind=pd)             :: yout,term
      real                      :: treset,eps,tset,tt
      real,dimension(maxbnd)    :: xx
      real,dimension(4)         :: t
      real,dimension(maxbnd,4)  :: x

      save t,tt,x,xx,nrec,nrold,ieof,eps
      data ieof,eps/0,0./

!       On the first call, read the first four variables and times, and
!       initialize the end of file flag ieof, and precision limit eps.

      if(itime==1) then
        ieof=0
        nrold=0
        eps=-0.5*tmin
        read(ifile4,*) t(1),(x(i,1),i=1,nbound)
        if(t(1)>time) then
          write(ifile3,2000)t(1),time,time
2000      format(' ***Warning: first value of time in boundary file,',&
           f10.3/'  exceeds first simulation time,',f10.3/'  first',&
           ' boundary time has been changed to ',f10.3/)
          t(1)=time
        endif
        do it=2,4
          read(ifile4,*) t(it),(x(i,it),i=1,nbound)
        enddo
        nrec=4
      endif

!       Test for simulation time less than t(1) due to time step rejection

      if(time<t(1)) then
        if(treset>0.) then
          nrold=nrold-1
          ireset=0
          treset=0.
        endif
        nbak=nrec-nrold+4
        do i=1,nbak
          backspace ifile4
        enddo
        do it=1,4
          read(ifile4,*) t(it),(x(i,it),i=1,nbound)
        enddo
        nrec=nrold
      else
        nrold=nrec
      endif

!       If an end of file is detected, write message and set flag.

      if(ieof==0) then

!         Test for reset condition and need to read another time
!         record from the boundary condition file.

        if((time-t(3)>=eps).and.(treset.le.0.)) then
20        read(ifile4,*,end=110) tt,(xx(i),i=1,nbound)
          nrec=nrec+1
          if(abs(tt-t(4)) < 1.0e-37) then
!          if(tt==t(4)) then
            treset=abs(tt)
          else
            do it=1,3
              t(it)=t(it+1)
              do i=1,nbound
                x(i,it)=x(i,it+1)
              enddo
            enddo
            t(4)=tt
            do i=1,nbound
              x(i,4)=xx(i)
            enddo

!             Test for need to read another time record.

            if((time-t(3)>=eps).and.(treset.le.0.)) goto 20
          endif
        endif

!       Test for simulation time greater than time of the last record.

      elseif(time>t(4)) then
        return
      endif

!       Test for reset condition. increment the reset counter if necessary

50    if((time-treset>=eps).and.(treset>0.)) ireset=ireset+1
      tset=time
      if(ireset==2) then
        do i=1,nbound
          x(i,1)=xx(i)
        enddo
        t(1)=tt
        do it=2,4
          read(ifile4,*) t(it),(x(i,it),i=1,nbound)
        enddo
        nrec=nrec+3
        tset=treset+tmin
      elseif(ireset==1) then
        tset=treset
      endif

!       Calculate interpolated values of the boundary variables.

      do ib=1,nbound
        nn=ibound(ib)
        yout=0.
        do i=1,4
          term=x(ib,i)
          do j=1,4
            if(i/=j .and. abs(t(i)-t(j))>1.0e-10) then
              term=term*dble(tset-t(j))/dble(t(i)-t(j))
            endif
          enddo
          yout=yout+term
        enddo
        tstate(nn)=yout
        state(nn)=yout
      enddo
      return

110   ieof=1
      print 1000, time
      write(ifile3,1000) time
      if(time.le.t(4)) goto 50

1000  format(1x,'***** End of file encountered on boundary condition'&
             ,' file at time',f10.2)

      return
      end subroutine bounds
!***********************************************************************

      subroutine indata(bldshl,ispos,isshel,isview,nstvec,tmax,&
        tstop,view,intshl,mout)

!   ----------------------------------------------------------------------
!
!       INDATA : Enter input data interactively before a simulation.
!
!       March 15, 1985 Cheol Park & Dan Clark
!
!       Added MOUT to select time base for writing output. Jan. 16, 1986
!       Added PRINT *,'=>' to indicate user's input.       April 25, 1989
!       Added FRZ to able or disable the freezing variables Jan. 16, 1997 C.P.
!       Updated to Fortran 90  April 27, 2007 C. Park
!
!       subprograms called: opnfil
!
!***********************************************************************

      use modsim_head
      implicit none

      logical          :: bldshl,intshl,view         ! changed 1/16/97
      character(len=1) :: answer
!      character(len=4),dimension(maxlbl) :: lbls
      integer          :: mout,nstvec,isview,isshel,idump,&
                          monitr,ians,i
      real             :: tmax,tstop
      integer,dimension(10) :: ispos

 !    Initialize default values

      bldshl=.false.
      intshl=.false.
      view=.false.
      frz = .true.
      init=0
      idump=0
      monitr=0
      iprint=0
      tprton=0.
      tprtof=0.
      ispos=0        ! added  12/3/99
      isshel=0
      isview=0
      nstvec=0

      print *
      print *,' Enter minimum time step, maximum time step,',&
              ' and simulation stopping time:'
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,*) tmin,tmax,tstop

      print *,'  Is the building shell model used? <n>   '
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,fmt='(a1)') answer
      if(answer=='y' .or. answer=='Y') bldshl=.true.
      print *,'  Will the initialization file be called? <n>   '
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,fmt='(a1)') answer
      if(answer=='y' .or. answer=='Y') then
        init=1
      endif

!    ask questions pertaining to building shell

      if(bldshl) then
        if(init==0) then
          print *,' Simulate building shell only? <n>   '
          write(*,fmt='(a4)',advance='no') ' => '
          read(inp,fmt='(a1)') answer
          if(answer=='y' .or. answer=='Y') then
            intshl=.true.
	        print *,'   Only the superblock containing the building',&
                    ' model will be called.'
          else
	        print *,'   All superblocks will be called.'
          endif
        endif

        print *
	    print *,'  What is the index number of the superblock for ',&
                'the building shell?    '
        write(*,fmt='(a4)',advance='no') ' => '
        read(inp,*) isshel
        print *,'  isshel =', isshel

        print *
	    print *,' Enter the time of day (in hours after midnight) '
        print *,' at which the simulation is to begin :  '
        write(*,fmt='(a4)',advance='no') ' => '
        read(inp,*) tzero

      endif

!    ask questions pertaining to all simulations

!    open input and output files

      call opnfil(bldshl)

      print *,'  -- The outputs can be written to the output file'
      print *,'  -- based on either simulation time or reported time.'
      print *
      print *,' Do you want to use reported time for outputs <n>?  '
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,fmt='(a1)') answer
      if(answer=='y' .or. answer=='Y') then
        mout=2
      else
        mout=1
      endif

      print *                                     ! added jan. 16, 1997
      print *,' Do you wish to disable freezing variable feature <n>?  '
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,fmt='(a1)') answer
      if(answer=='y' .or. answer=='Y') then
        frz = .false.
      else
        frz = .true.
      endif

      print *
      print *,' Do you want diagnostic information to be written <n>? '
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,fmt='(a1)') answer
      if(answer=='y' .or. answer=='Y') then
        print *,'  Five levels of detail are available:'
        print *
        print *,' 0:  No diagnostics'
        print *,' 1:  Trace which state variables are solved ',&
                'simultaneously'
        print *,' 2:  Trace simultaneously solved state variables and',&
                ' their initial'
        print *,'       and final values on each call to the equation',&
                ' solver'
        print *,' 3:  Trace simultaneously solved state variables and',&
                ' their estimated'
        print *,'       values on every iteration within the equation',&
                ' solver'
        print *,' 4:  Perform the previous trace and dump jacobian',&
                ' matrices'
        print *
        print *,'  Enter 0, 1, 2,3, or 4:'
        write(*,fmt='(a4)',advance='no') ' => '
        read(inp,*) ians
        if(ians==1) then
          iprint=-1
        elseif(ians==2) then
          iprint=10000
        elseif(ians==3) then
          iprint=1
        elseif(ians==4) then
          iprint=2
        endif
        if(ians>=1 .and. ians.le.4) then
	      print *,' Enter starting and stopping times for diagnostic',&
                  ' output:  '
          write(*,fmt='(a4)',advance='no') ' => '
          read(inp,*) tprton,tprtof
        endif
      endif

      print *,' Would you like to monitor simulation on screen? <n> '
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,fmt='(a1)') answer
      if(answer=='y' .or. answer=='Y') then
        view=.true.
    	print *,' Enter the index number of the superblock to monitor '
    	print *,' or zero (0) to monitor all superblocks.'
        write(*,fmt='(a4)',advance='no') ' => '
        read(inp,*) isview
100     print *,' --During the simulation, up to five state variables',&
                ' can be viewed-- '
    	print *,' Enter the number of state variables to be viewed.'
        write(*,fmt='(a4)',advance='no') ' => '
        read(inp,*) nstvec
        if(nstvec > 0) then
          if(nstvec>5) goto 100
          do i=1,maxlbl
            print *,'  ', i,' = ',lbls(i)
          enddo
          print *
    	  print *,' Enter the category number (above) and index number'
          print *,' for each of the ',nstvec,' variables to be viewed.'
          write(*,fmt='(a4)',advance='no') ' => '
          read(inp,*) (ispos(i),ispos(i+5),i=1,nstvec)
        endif
      endif

      return
      end subroutine indata
!***********************************************************************

      subroutine opnfil(bldshl)

!   ----------------------------------------------------------------------
!
!       OPNFIL:  Opens input and output files.
!
!         D.R. Clark,  National Bureau of Standards, December 21, 1984
!         Updated: March 8, 1985 C.P.
!
!	Modified: July 11, 1996 by P Haves, Loughborough University, UK
!	  - default filename and extensions changed to use unique extensions
!	  - file names may be up 46 characters.
!           - single filename may be entered in place of default filenames
!
!       Updated to Fortran 90  April 27, 2007 C. Park
!
!       subprograms called: none
!
!***********************************************************************

!   (NB modify character and format statements below when changing file
!    name lengths)
      use modsim_head
      implicit none

      logical                        :: bldshl
      character(len=11),dimension(8) :: deflt
      character(len=46)              :: fname,cr
      character(len=1)               :: ans,delim,charct,cr1
      character(len=7)               :: stat
      character(len=24),dimension(8) :: filtyp
      integer                        :: lname=46,ldeflt=11,lext=3
      integer                        :: numfil,i,icount,ideflt,j,ios
      integer,dimension(8)           :: ifile

      data deflt(1) /'hvacsim.dfn'/     ! old name : modeldef.dat 1/23/97
      data deflt(2) /'hvacsim.bnd'/     ! old name : boundary.dat
      data deflt(3) /'hvacsim.ini'/     ! old name : initin.dat
      data deflt(4) /'hvacsim.fin'/     ! old name : initout.dat
      data deflt(5) /'hvacsim.out'/     ! old name : modout.dat
      data deflt(6) /'hvacsim.sum'/     ! old name : modsum.dat
      data deflt(7) /'hvacsim.ctf'/     ! old name : ctfdata.dat
      data deflt(8) /'hvacsim.met'/     ! old name : weather.dat

      data filtyp(1) /'model definition file,  '/
      data filtyp(2) /'boundary variable file, '/
      data filtyp(3) /'initial state file,     '/
      data filtyp(4) /'final state file,       '/
      data filtyp(5) /'output data file,       '/
      data filtyp(6) /'simulation summary file,'/
      data filtyp(7) /'ctf file,               '/
      data filtyp(8) /'weather data file,      '/

      data cr /'                                              '/
      data stat /'old    '/
      data delim /'.'/
      data cr1/' '/
      ifile(1)=ifile2
      ifile(2)=ifile4
      ifile(3)=ifile5
      ifile(4)=ifile6
      ifile(5)=ifile1
      ifile(6)=ifile3
      ifile(7)=ifile7
      ifile(8)=ifile8

      if(bldshl) then
        numfil=8
      else
        numfil=6
      endif

      print *,' Use same file names for all files? (y/n) <y>  '
      write(*,fmt='(a4)',advance='no') ' => '
      read(inp,1)ans
   1  format(a1)
      if(ans=='n' .or. ans=='N')then
        do i=1,numfil
          if(i==3 .and. init/=1) cycle
          icount=0

!           request input

 100      print 2, filtyp(i),deflt(i)
   2      format(' Enter the name of the ',a24/' or carriage return',&
                 ' for default name:  ',a11)
          write(*,fmt='(a4)',advance='no') ' => '
          fname=cr
          read(inp,3)fname
   3      format(a46)
          if(fname==cr) then
            fname(1:ldeflt)=deflt(i)
          else

!           Add delimiter and extension to file name,
!           unless delimiter is already present

            do j=1,lname
              charct=fname(j:j)
              if(charct==delim) goto 400
              if(charct==' ') goto 300
            enddo
            if(fname(j:j)==delim)goto 400
 300        fname(j:j)=delim
            fname(j+1:j+lext)=deflt(i)(ldeflt-2:ldeflt)
          endif
! MAG moved 400 from endif line to below

!           Open file

400      if(i>=4 .and. i<=6) then
            open(unit=ifile(i),file=fname)
            close(unit=ifile(i),status='delete')
            stat='unknown'
          else
            stat='old    '
          endif
          open(unit=ifile(i),file=fname,iostat=ios,status=stat)
          if(ios/=0) then
            print 4, fname,ios
   4        format(' Error opening file ',a50,'     status =',i6)
            icount=icount+1
            if(icount>=3) then
              stop ' Task terminated'
            else
              goto 100
            endif
          endif
        enddo
      else

!         Use same names for all files

!         Get filename
520     write(*,530)
530     format(' Enter the name for all files to open or'/,&
               ' hit carriage return for default filename <hvacsim>')
        write(*,fmt='(a4)',advance='no') ' => '
        fname=cr
        read(inp,3) fname
        if (fname==cr1) then		 ! 7/11/97
          j = len('hvacsim.')
          fname(1:j-1)=deflt(1)(1:j-1)
        else
          do j=1,lname
             charct=fname(j:j)
             if (charct==delim) then
               print *,' Enter filename without extension'
               go to 520
             endif
             if (charct==' ') exit
          enddo
        end if

        do i=1,numfil
          if(i==3 .and. init/=1) cycle
          fname=fname(1:j-1)//delim//deflt(i)(ldeflt-2:ldeflt)
          write(*,*)
          write(*,fmt='(5x,"  File name : ",a46)') fname
          if(i>=4 .and. i<=6) then
            open(unit=ifile(i),file=fname,iostat=ios)
            if(ios/=0) then
              print *, fname,ios
   5          format(' Error opening file ',a50/' status =',i6)
              stop ' Task terminated'
            endif
            close(unit=ifile(i),status='delete')
            stat='unknown'
          else
            stat='old    '
          endif
          open(unit=ifile(i),file=fname,iostat=ios,status=stat)
          if(ios/=0) then
            print *, fname,ios
            stop ' Task terminated'
          endif
        enddo
      endif

      return
      end subroutine opnfil
!***********************************************************************

      subroutine rconf                        ! <<< 11/28/06
!      subroutine rconf(treprt)

!   ----------------------------------------------------------------------
!
!       RCONF : Read the model definition file.
!
!         C.R. Hill, National Bureau of Standards
!         February, 1983
!
!       Modifications by DRC, May 1984: read NSAVED from model definition
!       file; create COMMON/XINIT/INIT,NSAVED; if INIT = 1, read
!       STATE, TSTATE, STOLD, and SAVED from logical IFILE5.
!
!       Modified into FORTRAN77 by Cheol Park,  September 21, 1984
!
!       Updated:  Feb. 10, 1995 Cheol Park
!                 To read the new format of model definition file
!
!       Updated to Fortran 90  April 27, 2007 C. Park
!
!       subprograms called: none
!
!***********************************************************************

      use modsim_head
      implicit none

      character(len=80)         :: descript
      integer                   :: i,j,nn,isblk
      real                      :: time1,tstop1,tmax1,tmin1
!      real,dimension(maxsbk)    :: treprt           ! <<< 11/28/06

      read(ifile2,1000) descript                     ! added 2/10/95
      read(ifile2,1000) title
      read(ifile2,1000) descript
      read(ifile2,2000) nstate,nsblok
      read(ifile2,1000) descript
      read(ifile2,2000) (nsuper(i),i=1,nsblok)
      nblock=0
      do i=1,nsblok
        nblock=nblock+nsuper(i)
      enddo
      read(ifile2,1000) descript
      read(ifile2,3000) (state(i),i=1,nstate)
      read(ifile2,1000) descript
      read(ifile2,2000) (ndent(i),i=1,maxlbl)
      do i=1,nstate
        isign(i)=1
        stold(i)=state(i)
        tstate(i)=state(i)
      enddo
      read(ifile2,1000) descript
      read(ifile2,2000) (nunits(i),i=1,nblock)
      read(ifile2,1000) descript
      read(ifile2,2000) (njsslv(i),i=1,nsblok)
      read(ifile2,1000) descript
      read(ifile2,2000) (njsolv(i),i=1,nblock)
      read(ifile2,1000) descript
      do i=1,nsblok
        read(ifile2,2000) (isuper(i,j),j=1,nsuper(i))
      enddo
      read(ifile2,1000) descript
      do i=1,nblock
        read(ifile2,2000) (iblock(i,j),j=1,nunits(i))
      enddo
      nu=0
      do i=1,nblock
        nu=nu+nunits(i)
      enddo
      read(ifile2,1000) descript
      read(ifile2,2000) (iunits(i),i=1,nu)
      read(ifile2,1000) descript
      read(ifile2,2000) (nin(i),i=1,nu)
      read(ifile2,1000) descript
      do i=1,nu
        read(ifile2,2000) (in(i,j),j=1,nin(i))
      enddo
      read(ifile2,1000) descript
      read(ifile2,2000) (nout(i),i=1,nu)
      read(ifile2,1000) descript
      do i=1,nu
        read(ifile2,2000) (iout(i,j),j=1,nout(i))
      enddo
      read(ifile2,1000) descript
      read(ifile2,1000) descript
      do i=1,nsblok
        nn=njsslv(i)
        if(nn.le.0) nn=1
        read(ifile2,2000) (jssolv(i,j),j=1,nn)
      enddo
      read(ifile2,1000) descript
      read(ifile2,1000) descript
      do i=1,nblock
        nn=njsolv(i)
        if(nn.le.0) nn=1
        read(ifile2,2000) (jsolve(i,j),j=1,nn)
      enddo
      read(ifile2,1000) descript
      read(ifile2,2000) (nde(i),i=1,nu)
      read(ifile2,1000) descript
      nd=0
      do i=1,nu
        nn=nde(i)
        if(nn.le.0) nn=1
        read(ifile2,2000) (inde(i,j),j=1,nn)
        nd=nd+nde(i)
      enddo
      nn=nd
      if(nn.le.0) nn=1
      read(ifile2,1000) descript
      read(ifile2,2000) (idevar(i),i=1,nn)
      read(ifile2,1000) descript
      read(ifile2,2000) (isaved(i),i=1,nu)
      read(ifile2,1000) descript
      read(ifile2,2000) (jpar(i),i=1,nu)
      read(ifile2,1000) descript
      read(ifile2,2000) npar,nsaved
      read(ifile2,1000) descript
      read(ifile2,3000) (par(i),i=1,npar)
      read(ifile2,1000) descript
      read(ifile2,2000) nbound
      read(ifile2,1000) descript
      read(ifile2,2000) (ibound(i),i=1,nbound)
      read(ifile2,1000) descript
      read(ifile2,2000) (nreprt(i),i=1,nsblok)
      read(ifile2,1000) descript
      read(ifile2,3000) (treprt(i),i=1,nsblok)
      read(ifile2,1000) descript
      read(ifile2,1000) descript
      read(ifile2,1000) descript
      do isblk=1,nsblok
        read(ifile2,2000) (ireprt(isblk,i),i=1,nreprt(isblk))
        read(ifile2,2000) (ident(isblk,i,1),i=1,nreprt(isblk))
        read(ifile2,2000) (ident(isblk,i,2),i=1,nreprt(isblk))
      enddo
      read(ifile2,1000) descript
      read(ifile2,3000) rtolx,atolx,xtola,ttime
      read(ifile2,1000) descript
      read(ifile2,2000) (ifzopt(i),i=1,nsblok)
      read(ifile2,1000) descript
      read(ifile2,2000) (insopt(i),i=1,nsblok)
      close(ifile2)                                   ! added 2/10/95

!       If init is 1, read the initialization file.

      if(init==1) then
        read(ifile5,4000) tmin1,tmax1,tstop1,time1,nstate,nsaved
        read(ifile5,3000) (state(i), i=1,nstate)
        read(ifile5,5000) (isign(i), i=1,nstate)
        do i=1,nstate
          tstate(i)=state(i)
          stold(i)=state(i)
        enddo
        if(nsaved>0) then
          read(ifile5,3000) (saved(i), i=1,nsaved)
        endif
      else
        do i=1,nsaved
          saved(i)=0.0
        enddo
      endif

1000  format(a80)
2000  format(16i5)               !  03/07/2007
! 2000  format(20i4)
3000  format(5g15.6)
4000  format(4f10.2,2i5)
5000  format(25i3)

      return
      end subroutine rconf
!***********************************************************************

      subroutine report(tstop,ireset,isblk,mout,label)

!   ----------------------------------------------------------------------
!
!       REPORT :
!         Write reports to the printer. Reports are made at equally spaced
!         times, based on interpolated values of the reported variables.
!
!         C.R. Hill, National Bureau of Standards
!         February, 1983
!
!         Modification was made by Dan Clark, April 1984.
!
!         Modified into FORTRAN77 by Cheol Park, September 14, 1984
!
!         Updated on : Nov. 20, 1984 D.C.
!
!         Add MOUT to select time base for writing output. Jan. 7, 1987
!
!         June 27, 1997  Cheol Park, NIST
!         To adjust the size of the simulation summary file, hvacsim.sum,
!         the maximum time for writing data,WRTMAX, is defined in the file,
!         hvacsim.par.
!
!         Updated to Fortran 90  April 27, 2007 C. Park
!
!         subroutines called: none
!
!***********************************************************************

      use modsim_head
      implicit none

      character(len=4),dimension(maxlbl) :: label,line
      integer                            :: mout,isblk,ireset,nn,i,j,nline,&
                                            l,ir,nvar,id
      integer,dimension(8)               :: number
      real                               :: tstop
      real,dimension(mrptis)             :: var
      real,dimension(maxsbk,2)           :: rtime
      real,dimension(maxsbk,mrptis,2)    :: rvar

      save rtime,rvar

!      data wrtmax/ 3600.0 / ! wrtmax is defined in hvacsim.par

      nn=nreprt(isblk)

!       Make report to the data file at each simulation time.

      if(mout==1) then
        write(ifile1,1000) isblk,time
        do i=1,nn
          var(i)=state(ireprt(isblk,i))
        enddo
        write(ifile1,2000) (var(i),i=1,nn)
      endif

!       On the first call, store the reported variables and times and
!       return.

!-ph  Previous values should not be used for interpolation immediately
!-ph  following reset.
!-ph
!-ph  if((itime<=1).and.(ireset==2)) then
      if((itime<=1).or.(ireset==2)) then
        rtime(isblk,2)=time
        do i=1,nn
          rvar(isblk,i,2)=state(ireprt(isblk,i))
        enddo
        return
      endif
      rtime(isblk,1)=rtime(isblk,2)
      do j=1,nn
        rvar(isblk,j,1)=rvar(isblk,j,2)
      enddo
      rtime(isblk,2)=time
      do i=1,nn
        rvar(isblk,i,2)=state(ireprt(isblk,i))
      enddo

!       Test if a report is due.

      if(time<tnext(isblk)) return

!       Make report to the list file.

50    if(time <= wrtmax+tstep ) then
        write(ifile3,3000) isblk
        write(ifile3,4000) tnext(isblk)
      endif
      nline=nn/8
      if(8*nline.ne.nn) nline=nline+1

!       Begin loop to write reports.

      do l=1,nline
        ir=8*(l-1)
        nvar=nn-ir
        if(nvar>8) nvar=8

!       Begin loop to write a single line.

        do i=1,nvar
          id=ident(isblk,ir+i,1)

!         Calculate interpolated values.
 
          var(i)=rvar(isblk,ir+i,1)+&
                  (rvar(isblk,ir+i,2)-rvar(isblk,ir+i,1))*&
                  (tnext(isblk)-rtime(isblk,1))/&
                  (rtime(isblk,2)-rtime(isblk,1))

!       Assemble variable labels.

          number(i)=ident(isblk,ir+i,2)

          line(i)=lbls(id)                       ! <<< 11/30/06
        enddo
 
        if(time <= wrtmax+tstep) then
          write(ifile3,5000) (line(i),number(i),i=1,nvar)
          write(ifile3,6000) (var(i),i=1,nvar)
        endif
      enddo

!       Make report to the data file at the reported time.

      if(mout==2) then
        write(ifile1,1000) isblk,tnext(isblk)
        do i=1,nn
          var(i)=state(ireprt(isblk,i))
        enddo
        write(ifile1,2000) (var(i),i=1,nn)
      endif

!       Write the reported time on a terminal screen or equivalent to
!       monitor the simulation progress.

      if(isblk==1) then
        print 4000, tnext(isblk)
      endif

!       Set time for next report. if next report is not due, return.

      tnext(isblk)=tnext(isblk)+treprt(isblk)         
      if((time>=tnext(isblk)).and.(tnext(isblk)<=tstop)) goto 50

1000  format('SUPERBLOCK',i2,f10.2)
2000  format(1p5g15.6)
3000  format(/1x,'********** SUPERBLOCK',i3,' **********')
4000  format(1x,'time=',f10.2)
5000  format(1x,8(a4,i3,3x))
6000  format(1x,1p8g10.3)

      return
      end subroutine report
!***********************************************************************

      subroutine sumary(tmax,tstop)

!   ----------------------------------------------------------------------
!
!       SUMARY : Write a summary of configuration of simulation
!
!         C.R. Hill, National Bureau of Standards
!         February, 1983
!
!         Modified into FORTRAN77 by Cheol Park, August 16, 1984
!         Updated to Fortran 90  April 27, 2007 C. Park
!
!         subprograms called: none
!
!***********************************************************************

      use modsim_head
      implicit none

      character(len=4),dimension(msumary)  :: label
      integer                              :: nlbl,i,ns1,j,isblk,nn,nnn,i1,k,&
                                              itype,nnnn,kk,iu       ! 11/28/06
      integer,dimension(msumary)           :: number
      integer,dimension(1,maxlbl)          :: lbl
      real                                 :: tmax,tstop

      write(ifile3,1000)
      write(ifile3,2000) title
      write(ifile3,3000) nsblok,nblock,nu
      nlbl=0
      do i=2,maxlbl
        if(ndent(i)<=0) cycle
        nlbl=nlbl+1
        lbl(1,nlbl)=ndent(i)
      enddo
      nlbl=nlbl+1
      lbl(1,nlbl)=nstate
      call vlab(nlbl,1,maxlbl,lbl,1,label,number)
      write(ifile3,3200)nstate
      write(ifile3,3400) (number(i),label(i),i=1,nlbl)
      write(ifile3,3600)
      ns1=1
      do i=1,nlbl
        write(ifile3,3800) label(i),(state(j),j=ns1,lbl(1,i))
        ns1=lbl(1,i)+1
      enddo
      write(ifile3,3900) nbound
      if(nbound>0) then
        call vlab(nbound,1,maxbnd,ibound,1,label,number)
        write(ifile3,6000) (label(j),number(j),j=1,nbound)
      endif
      write(ifile3,13000)
      write(ifile3,12000) rtolx,atolx,xtola,ttime        
      do isblk=1,nsblok
        write(ifile3,4000) isblk
        write(ifile3,4200) ifzopt(isblk)
        write(ifile3,4400) insopt(isblk)
        nn=nreprt(isblk)
        write(ifile3,4600) nn
        if(nn>0) then
          call vlab(nn,maxsbk,mrptis,ireprt,isblk,label,number)
          write(ifile3,6000) (label(j),number(j),j=1,nn)
        endif
        nnn=njsslv(isblk)
        write(ifile3,5000) nnn
        if(nnn>0) then
          call vlab(nnn,maxsbk,mseqis,jssolv,isblk,label,number)
          write(ifile3,6000) (label(j),number(j),j=1,nnn)
        endif
        do i1=1,nsuper(isblk)
          i=isuper(isblk,i1)
          write(ifile3,7000) i
          nnn=njsolv(i)
          write(ifile3,5000) nnn
          if(nnn>0) then
            call vlab(nnn,maxblk,mseqib,jsolve,i,label,number)
            write(ifile3,6000) (label(k),number(k),k=1,nnn)
          endif
          do j=1,nunits(i)
            iu=iblock(i,j)
            itype=iunits(iu)
            write(ifile3,8000) iu,itype
            nnn=nin(iu)
            write(ifile3,9000) nnn
            if(nnn>0) then
              call vlab(nnn,maxunt,minoiu,in,iu,label,number)
              write(ifile3,6000) (label(k),number(k),k=1,nnn)
            endif
            nnn=nout(iu)
            write(ifile3,10000) nnn
            if(nnn>0) then
              call vlab(nnn,maxunt,minoiu,iout,iu,label,number)
              write(ifile3,6000) (label(k),number(k),k=1,nnn)
            endif
            write(ifile3,11000)
            nnn=jpar(iu)
            nnnn=npar
            do kk=1,nu
              if(kk/=iu .and. jpar(kk)>=nnn)&
                                  nnnn=min(nnnn,jpar(kk)-1)
            enddo
            if(nnnn<nnn) then
              write(ifile3,11500)
            else
              write(ifile3,12000) (par(k),k=nnn,nnnn)
            endif
          enddo
        enddo
      enddo
      write(ifile3,14000) tmin,tmax,tstop

1000  format(10x,'***** Program MODSIM *****'/8x,&
             'A MODular SIMulation program'/)
2000  format(a80)
3000  format(1x,/1x,i2,' superblocks',5x,i3,' blocks',5x,i4,&
     &       ' units')
3200  format(/1x,i4,' state variables:')
3400  format(8(2x,i3,1x,a4))
3600  format(1x/1x,'initial state vector:')
3800  format(/2x,a4,':'/(2x,1p5g15.6))
3900  format(/1x,i2,' time dependent boundary variables:')
4000  format(1x,/1x,'***** superblock',i3,' *****')
4200  format(/'  superblock simultaneous equation unfreezing option,',&
              ' ifzopt =',i3)
4400  format('  superblock input scan option, insopt =',i3)
4600  format(/1x,i2,' reported variables:')
5000  format(1x/1x,i2,' simultaneous equations; variables:')
6000  format(8(1x,a4,i3,2x))
7000  format(1x,/1x,'***** block',i3,' *****')
8000  format(1x/1x,'unit',i4,5x,'type',i3)
9000  format(1x/1x,i3,' inputs:')
10000 format(1x/1x,i3,' outputs:')
11000 format(1x/1x,'parameters:')
11500 format(1x,'  none')
12000 format(1x,1p5g15.5)
13000 format(/1x,'error tolerances:  rtolx, atolx, xtol, ttime:')
14000 format(/1x,'--------------------------------------------------',&
             '---------'/1x,'tmin =',f12.3,'  tmax =',f12.3,'  tstop =',&
             f12.3/1x,'-----------------------------------------------',&
             '------------'/)
!
      return
      end subroutine sumary

!***********************************************************************

      subroutine vlab(n,m1,m2,ivar,indx,label,number)

!   ----------------------------------------------------------------------
!
!       VLAB : Return label and number for variables in IVAR
!
!       C.R. Hill, National Bureau of Standards
!       February, 1983
!
!       Updated : Nov. 20, 1984 D.C.
!       Updated to Fortran 90  April 27, 2007 C. Park
!
!       subroutines called: none
!
!***********************************************************************
      use modsim_head
      implicit none

      integer                              :: n  
      character(len=4)                     :: null='null'
      character(len=4),dimension(n)        :: label
      integer                              :: indx,m1,m2,i,j,k
      integer,dimension(n)                 :: number
      integer,dimension(msumary)           :: id
      integer,dimension(m1,m2)             :: ivar

!       Note: dimension of id is max(mseqis,mseqib,minoiu,mrptis,maxlbl)

      ndnt=ndent                          ! vector operation
      do i=1,n
        do j=1,maxlbl
          k=j
          if(ndnt(j)>=0) exit
        enddo
        do j=2,maxlbl
          if(ndnt(j)>0) then
            if(ivar(indx,i)<=ndnt(j)) exit
            k=j
          endif
        enddo
        number(i)=ivar(indx,i)-ndnt(k)
        id(i)=k
      enddo
      do i=1,n
        label(i)=lbls(id(i))
        if(number(i)<=0) label(i)=null
      enddo

      return
      end subroutine vlab

