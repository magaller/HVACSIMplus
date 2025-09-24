!=======================================================================
! SLIMCON SOURCE  VERSION 6.0
!=======================================================================
!
!  Created : Nov. 28, 1984   C.R. Hill and D. Clark (NBS)      ver. 2.1
!  Revised : Sep. 11, 1990   C. Park (NIST) and P. Haves (UOX) ver. 4.0
!  Revised : February 8, 1995  C. Park (NIST)                  ver. 5.2
!
! - version 5.3  July 9, 1996:
!                Changes made by P Haves, Loughborough University, UK to
!                support the "simulation test-bed" produced by ASHRAE 825-RP
!                (1) Long file names supported (46 chars + '.' + 3 char
!                    extension = 50 chars)
!                (2) SUBROUTINE FILEOP: uses ".sim" instead of ".SIM" and
!                    ".dfn" instead of ".DAT". Lower case is used for
!                    convenience when file names are case sensitive and ".dfn"
!                    is used as a unique file extension for model definition
!                    files
!                (3) SUBROUTINE FILEOP: Tries to open TYPES data file TYPARF in
!                    current directory or file TYPARFF specified with full path
!                    - allows custom version in current directory to be used if
!                    present, otherwise default is used.
!                (4) Allow algebraic variables as well as differential
!                    variables to be inputs and oututs of the same unit -
!                    enables equation solver to be used to solve algebriac
!                    equations within TYPEs.
!
! - version 5.4  Oct. 30, 1997:   Cheol Park, NIST
!                Made it an option to allow algebraic variables to be
!                inputs and outputs of the same unit for the equation solver
!
! - version 6.0  April 10, 2007  C. Park, NIST
!                Updated code from Fortran 77 to Fortran 90.
!
! ----------------------------------------------------------------------
!
! This source code, when compiled with a FORTRAN 90/95 compiler, and linked
! will produce the program SLIMCON. SLIMCON is used to produce a model
! definition file which describes the structure and characteristics
! of the system to be simulated by the main simulation program, MODSIM.
! SLIMCON is part of the HVACSIM+ dynamic simulation package for
! building heating, ventilation, and air conditions systems and controls,
! plus other building systems. Simulation work files produced by the
! program HVACGEN are processed by the program SLIMCON to produce a
! model definition file.
!
!=======================================================================
! The model definition file is written to logical unit LUD. Input to
! SLIMCON is read from logical unit LUF.
! Numerical values are assigned in subroutine FILEOP.
!
!=======================================================================
!
!------- The following parameters are used to define the array sizes
!        in the program. To change a parameter, this change must be
!        made in ALL subroutines in the program. The meanings of the
!        parameters are:
!
!           maxsbk = maximum superblocks in the simulation
!           maxblk = maximum blocks in the simulation
!           maxunt = maximum units in the simulation
!           maxstv = maximum state variables in the simulation
!           maxdeq = maximum differential equations in the simulation
!           maxpar = maximum unit parameters in the simulation
!           maxsav = maximum saved variables in the simulation
!           maxbnd = maximum boundary variables in the simulation
!           mblkis = maximum blocks in a superblock
!           muntib = maximum units in a block
!           mdeqis = maximum differential equations in a superblock
!           mdeqiu = maximum differential equations in a unit
!           mseqis = maximum simultaneous equations in a superblock
!           mseqib = maximum simultaneous equations in a block
!           mbndis = maximum boundary conditions in a superblock
!           minoib = maximum inputs or outputs in a block
!           minoiu = maximum inputs or outputs in a unit
!           mrptis = maximum reported variables in a superblock
!           mpariu = maximum parameters in a unit
!           maxlbl = maximum number of variable categories (labels)
!
! *********************************************************************

      program   slimcon

      use hvacsim_par
      use slim
      implicit none
      integer,dimension(maxunt)   :: iunits
      integer,dimension(maxstv)   :: iouts,ioutu,ioutn
      integer                     :: isav,nstout,isblk,ib0,nbk,nu,id,jjpar,&
                                     nn,k,idnt,ndeq,nsave,nub,&
                                     nio,np,nd,nstate,iub,jb,npar,nnsolv,&
                                     nnn,iid,iu1,iu2,nnnn,ivar,nse,ib1,&
                                     ib2,ni,no,nrv,nr
      integer                     :: nsblok,itype,iude,nnin,nnout,&
                                     nnpar,nsaved,nblock,nbound,iu 
      integer                     :: i,ii,iblk,j,l,jj,kk                               
                                            

!     Open work and model definition files

10      call fileop

        write(*,*) '----------------------------------------'
        write(*,*) ' iu  itype nsaved iude nnin nnout nnpar '
        write(*,*) '----------------------------------------'

!     Title read in from work file.

        read(luf,fmt='(a80)',iostat=ioerr) title
         if(ioerr /= 0) then
          go to 10
        end if

!     Read number of superblocks from the work file.
 
        read(luf,*) nsblok

!     Error tolerances read from the work file

        read(luf,*) rtolx,atolx,xtol,ttime

!     Initialize for main loop to define superblocks, blocks and unit
!     connections.

        id     = 1
        jjpar  = 1
        isav   = 1
        nblock = 0
        nstout = 0
        do i=1,maxlbl
         ndent(i) = 0
        end do

!     Begin main loop to define superblocks, blocks and unit
!     connections.

        do isblk=1,nsblok
         ndeqsb(isblk) = 0

!     Input number of blocks in superblock isblk.

         read(luf,*) nsuper(isblk)
         ib0 = nblock
         nblock = nblock + nsuper(isblk)
         nbk = nsuper(isblk)

!     Begin subloop to define blocks.

         do ii=1,nbk
          iblk = ib0 + ii
          isuper(isblk,ii) = iblk

!     Input number of units in block iblk.

          read(luf,*) nunits(iblk)
          nu = nunits(iblk)

!     Begin subloop to define unit connections.

          do i=1,nu

!     Input unit number and type number.

           read(luf,*) iu,itype

!     Call typar to obtain parameters for the type module assigned to
!     unit iu.

           call typar(itype,nsaved,iude,nnin,nnout,nnpar,0)
           isaved(iu) = isav
           isav = isav + nsaved
           iblock(iblk,i) = iu
           iunits(iu) = itype
           nin(iu) = nnin
           nout(iu) = nnout
           maxpun(iu) = nnpar
           write(*,5000) iu,itype,nsaved,iude,nnin,nnout,nnpar
           5000 format(2i5,i7,4i5)

!           5000 format('  unit,type,nsaved,nde,nin,nout,npar=', &
!                        2i5,i7,4i5)

!     Input the input connections for unit iu.

           nn = nin(iu)
           read(luf,*) (in(iu,k),k=1,nn)
           do k=1,nn
            idnt = iident(k)
            ndent(idnt) = max(in(iu,k),ndent(idnt))
           end do

!     Input the output connections for unit iu.

           nn = nout(iu)
           read(luf,*) (iout(iu,k),k=1,nn)
           do k=1,nn
            idnt = iodent(k)
            ndent(idnt) = max(iout(iu,k),ndent(idnt))
           end do

!     Input parameters for unit iu.

           nn = jjpar + nnpar - 1
           read(luf,*) (par(k),k=jjpar,nn)
           jpar(iu) = 0
           jpar(iu) = jjpar
           jjpar= nn + 1

!     Define the connections for differential equations in unit iu.

           nn = jpar(iu)
           call typar(itype,nsaved,iude,nnin,nnout,nnpar,nn)
           nde(iu) = iude
           inde(iu,1) = 0
           idevar(id) = 0

           nn = nde(iu)
           ndeqsb(isblk) = iude + ndeqsb(isblk)
           if(nn>0) then
            do k=1,nn
             inde(iu,k) = id
             id = id + 1
            end do
           end if
          end do
         end do
        end do

        write(*,*)
        write(*,*) '                   ========= SLIMCON SUMMARY =========='
        write(*,*)
        call report(nsblok,maxsbk,6)
        call report(nblock,maxblk,1)
        ndeq = id - 1
        if(ndeq<=0) id = 1
        call report(ndeq,maxdeq,3)
        nsave = isav - 1
        call report(nsave,maxsav,5)
        nu = 0
        nub = 0
        do i=1,nblock
          nn = nunits(i)
          if(nn>nub) nub = nn
          nu = nu + nn
        end do
        call report(nu,maxunt,8)
        call report(nub,muntib,19)
        nn  = 0
        nio = 0
        np  = 0
        do i = 1,nu
          if(nde(i)>nn)    nn = nde(i)
          if(nin(i)>nio)   nio = nin(i)
          if(nout(i)>nio)  nio = nout(i)
          if(maxpun(i)>np) np = maxpun(i)
        end do
        call report(nn,mdeqiu,11)
        call report(nio,minoiu,14)
        call report(np,mpariu,15)
        nu = 0
        nd = 0
        do i=1,nsblok
          nn = nsuper(i)
          if(nn>nu) nu = nn
          if(ndeqsb(i)>nd) nd = ndeqsb(i)
        end do
        call report(nu,mblkis,9)
        call report(nd,mdeqis,10)

!     Define the structure of the state vector.

        nstate = 0
        do i=1,maxlbl
          ndnt(i) = nstate
          nstate= nstate + ndent(i)
        end do
        call report(nstate,maxstv,7)

!     Set elements of the state vector structure definition array to
!     -1  which correspond to null categories of variables.

        do i=1,maxlbl
          if(ndent(i)<=0) ndnt(i) = -1
        end do

!     Begin loop to assign storage for the state vector.

        do iblk=1,nblock
         maxibl(iblk) = 0
         maxobl(iblk) = 0
         nu = nunits(iblk)
         do i=1,nu
          iu = iblock(iblk,i)
          call typar(iunits(iu),nsaved,iude,nnin,nnout,nnpar,0)
          maxibl(iblk) = maxibl(iblk) + nnin
          maxobl(iblk) = maxobl(iblk) + nnout

!     Assign indices for the input connections.

          nn = nin(iu)
          loop_160: do j=1,nn
           k = iident(j)
           in(iu,j) = ndnt(k) + in(iu,j)
           do iub = iu-i+1,iu
	       if ( iub <= 0 ) cycle loop_160
            do jb = 1,nin(iub)
              if(.not.(jb>=j.and.iu==iub)) then
                if(in(iu,j)==in(iub,jb)) then
                  maxibl(iblk) = maxibl(iblk) - 1
                  cycle loop_160
                end if
              end if
            end do
           end do
          end do loop_160

!     Assign indices for the output connections.

          nn = nout(iu)
          do j=1,nn
           if(iout(iu,j)==0)  then
            maxobl(iblk) = maxobl(iblk) - 1
           elseif(iout(iu,j)>0) then
            k = iodent(j)
            iout(iu,j) = ndnt(k) + iout(iu,j)
            nstout = nstout + 1
            iouts(nstout) = iout(iu,j)
            ioutu(nstout) = iu
            ioutn(nstout) = j
           end if
          end do

!     Assign indices for the differential equations.

          nn = nde(iu)
          if(nn>0) then
           do j=1,nn
            k = inde(iu,j)
            idevar(k) = iout(iu,j)
           end do
          end if
         end do
        end do

        nio = 0
        do i=1,nblock
          nn = nunits(i)
          if(maxibl(i)>nio)  nio = maxibl(i)
          if(maxobl(i)>nio)  nio = maxobl(i)
        end do
        call report(nio,minoib,12)

        npar = jjpar - 1
        call report(npar,maxpar,4)

!     Begin loop to define variables which are solved simultaneously
!     within blocks.

        do iblk=1,nblock
         nn = nunits(iblk)
         nnsolv = 0
         isolve(iblk,1) = 0
         do i=1,nn
          iu = iblock(iblk,i)
          nnn = nde(iu)
          if(nnn>0) then
           do ii=1,nnn
            iid = inde(iu,ii)
            nnsolv = nnsolv + 1
            isolve(iblk,nnsolv) = idevar(iid)
           end do
          end if
         end do
         do i=1,nn
          loop_270: do ii=1,nn
!-ph  Don't exclude variables which are inputs and outputs of same unit,
!     if algbraic is true.

           if((i==ii) .and. (.not.algbraic)) cycle loop_270      ! 10/30/97 cp
           iu1 = iblock(iblk,i)
           iu2 = iblock(iblk,ii)
           nnn = nin(iu1)
           nnnn = nout(iu2)
           loop_260: do j=(nde(iu2)+1),nnnn
            ivar = iout(iu2,j)
            if(ivar>0) then
             do k=1,nnn
              if(ivar==in(iu1,k)) then
               do l=1,nnsolv
                if(isolve(iblk,l)==ivar) cycle loop_260
               end do
               nnsolv = nnsolv + 1
               isolve(iblk,nnsolv) = ivar
              end if
             end do
            end if
           end do loop_260
          end do loop_270
         end do
         nsolve(iblk) = nnsolv
        end do

!     End loop to define block equations.

        nse = 0
        do i=1,nblock
          if(nsolve(i)>nse) nse = nsolve(i)
        end do
        call report(nse,mseqib,17)

!     Begin loop to define variables solved simultaneously within
!     superblocks.

        do isblk=1,nsblok
         nbk = nsuper(isblk)
         nnsolv = 0
         issolv(isblk,1) = 0
         do i=1,nbk
          ib1 = isuper(isblk,i)
          nn = nunits(ib1)
          do ii=1,nbk
           if(i/=ii) then
            ib2 = isuper(isblk,ii)
            nnn = nunits(ib2)
            do j=1,nn
             iu1 = iblock(ib1,j)
             ni = nin(iu1)
             do jj=1,nnn
              iu2 = iblock(ib2,jj)
              no = nout(iu2)
              loop_340: do k=1,ni
               do kk=1,no
                ivar = iout(iu2,kk)
                if(in(iu1,k)==ivar .and. ivar/=0) then
                 do l=1,nnsolv
                  if(issolv(isblk,l)==ivar) cycle loop_340
                 end do
                 nnnn = nsolve(ib2)
                 do l=1,nnnn
                  if(isolve(ib2,l)==ivar) cycle loop_340
                 end do
                 nnsolv = nnsolv + 1
                 issolv(isblk,nnsolv) = ivar
                end if
               end do
              end do loop_340
             end do
            end do
           end if
          end do
         end do
         nssolv(isblk) = nnsolv
        end do

!     End loop to define superblock equations.

        nse = 0
        do i=1,nsblok
          if(nssolv(i)>nse) nse = nssolv(i)
        end do
        call report(nse,mseqis,18)

!     Input initial state vector.

        read(luf,*) (state(i),i=1,nstate)

!     Input variables which are time dependent boundary conditions.

        read(luf,*) nbound
        call report(nbound,maxbnd,2)
        read(luf,*) (ibound(i),i=1,nbound)
        call tdbvis(nsblok,nbound)

!     Input variables which are to be reported during the simulation.

        nrv = 0
        do isblk=1,nsblok
         read(luf,*) nreprt(isblk),treprt(isblk)
         if(nreprt(isblk)>nrv) nrv = nreprt(isblk)
         nr = nreprt(isblk)
         read(luf,*) (ireprt(isblk,i), i=1,nr)
         read(luf,*) (ident(isblk,i,1),i=1,nr)
         read(luf,*) (ident(isblk,i,2),i=1,nr)
        end do
        call report(nrv,mrptis,16)
        write(*,*) 
 
!     Input superblock variable freezing options.

        read(luf,*) (ifzopt(i),i=1,nsblok)

!     Input superblock input scan options.

        read(luf,*) (insopt(i),i=1,nsblok)
 
!     Write data to the model definition file.

        write(lud,*) 'title: simulation title'
        write(lud,fmt='(1x,a79)') title
        write(lud,*) 'nstate,nsblok: # of state variables, # of SBS'
        write(lud,3000) nstate,nsblok
        3000 format(16i5)
!        3000 format(20i4)
        write(lud,*) 'nsuper(s): # of blocks in each SB'
        write(lud,3000) (nsuper(i),i=1,nsblok)
        write(lud,*) 'state(i): vector of state variable initial values'
        write(lud,4000) (state(i),i=1,nstate)
        4000 format(5e15.6)
        write(lud,*) 'ndent(i): state variable identification vector'
        write(lud,3000) (ndnt(i),i=1,8)
        write(lud,*) 'nunits(b): # of units in each block B'
        write(lud,3000) (nunits(i),i=1,nblock)
        write(lud,*) 'njsslv(s): # of simultaneous eqs in each SB'
        write(lud,3000) (nssolv(i),i=1,nsblok)
        write(lud,*) 'njsolv(b): # of simultaneqou eqs in each block'
        write(lud,3000) (nsolve(i),i=1,nblock)
        write(lud,*) 'isuper(s,i): array of block numbers in each SB'
        do i=1,nsblok
         nn = nsuper(i)
         write(lud,3000) (isuper(i,j),j=1,nn)
        end do
        write(lud,*) 'iblock(b,i): array of unit numbers in each block'
        nu = 0
        do i=1,nblock
         nn = nunits(i)
         write(lud,3000) (iblock(i,j),j=1,nn)
         nu =nu + nn
        end do
        write(lud,*) 'iunits(u): array of type #s for each unit'
        write(lud,3000) (iunits(i),i=1,nu)
        write(lud,*) 'nin(u): number of inputs to unit u'
        write(lud,3000) (nin(i),i=1,nu)
        write(lud,*) 'in(u,i): array of input connections for unit U'
        do i=1,nu
         nn = nin(i)
         write(lud,3000) (in(i,j),j=1,nn)
        end do
        write(lud,*) 'nout(u): number of outputs from unit U'
        write(lud,3000) (nout(i),i=1,nu)
        write(lud,*) 'iout(u,i): array of output connections for unit U'
        do i=1,nu
         nn = nout(i)
         write(lud,3000) (iout(i,j),j=1,nn)
        end do
        write(lud,*) 'jssolv(s,i): array of variables solved'
        write(lud,*) ' simultaneously within each sb (between blocks)'
        do i=1,nsblok
         nn = nssolv(i)
         if(nn<=0) nn = 1
         write(lud,3000) (issolv(i,j),j=1,nn)
        end do
        write(lud,*) 'jsolve(b,i): array of variables solved'
        write(lud,*) ' simultaneously within each block'
        do i=1,nblock
         nn = nsolve(i)
         if(nn<=0) nn = 1
         write(lud,3000) (isolve(i,j),j=1,nn)
        end do
        write(lud,*) 'nde(u): # of differential eqs in unit U'
        write(lud,3000) (nde(i),i=1,nu)
        write(lud,*) 'inde(u,i): de index for the ith de in unit U'
        do i=1,nu
         nn = nde(i)
         if(nn<=0) nn = 1
         write(lud,3000) (inde(i,j),j=1,nn)
        end do
        id = id - 1
        if(id<=0) id = 1
        write(lud,*) 'idevar(d): variable index for de #d'
        write(lud,3000) (idevar(i),i=1,id)
        write(lud,*) 'isaved(u): index of first saved var. for unit U'
        write(lud,3000) (isaved(i),i=1,nu)
        write(lud,*) 'jpar(u): index of first parameter for unit U'
        write(lud,3000) (jpar(i),i=1,nu)
        write(lud,*) 'npar,nsaved: # of parameters & saved variables'
        isav = isav - 1
        write(lud,3000) npar,isav
        write(lud,*) 'par(i): array of parameters for all units'
        write(lud,4000) (par(i),i=1,npar)
        write(lud,*) 'nbound: # of time-dependent boundary variables'
        write(lud,3000) nbound
        write(lud,*) 'ibound(i) state variable indices of boundary variables.'
        write(lud,3000) (ibound(i),i=1,nbound)
        write(lud,*) 'nreprt(s): # of reported variables in each SB'
        write(lud,3000) (nreprt(i),i=1,nsblok)
        write(lud,*) 'treprt(s): reporting interval for each SB'
        write(lud,4000) (treprt(i),i=1,nsblok)
        write(lud,*) 'ireprt(i): indices of reported variables'
        write(lud,*) 'ident1: category # of reported variables'
        write(lud,*) 'ident2: position in category of reported variables'
        do isblk=1,nsblok
         nr = nreprt(isblk)
         write(lud,3000) (ireprt(isblk,i),i=1,nr)
         write(lud,3000) (ident(isblk,i,1),i=1,nr)
         write(lud,3000) (ident(isblk,i,2),i=1,nr)
        end do
        write(lud,*) 'rtolx,atolx, xtol,ttime: error tolerances'
        write(lud,4000) rtolx,atolx,xtol,ttime
        write(lud,*) 'ifzopt(s): sb variable unfreezing option vector'
        write(lud,3000) (ifzopt(i),i=1,nsblok)
        write(lud,*) 'insopt(s): sb input scan option vector'
        write(lud,3000) (insopt(i),i=1,nsblok)

        call varchk(nblock,nsblok,nbound)
        call outchk(nstout,iouts,ioutu,ioutn)

        stop 'Model definition file completed'
        end program slimcon

!**********************************************************************

      subroutine fileop

! ---------------------------------------------------------------------
!
!     This routine opens the input and output files
!
!**********************************************************************

      use slim
      implicit none
      character(len=1)          :: charct,answer,delim = '.'
      integer                   :: lname = 46, lext = 3
      character(len=50)         :: fullnm
      character(len=12)         :: typarf='typar.dat   '
      integer                   :: i,option

      write(*,*) ' ----------------------------------------------------------'
      write(*,*) ' *                                                        *'
      write(*,*) ' *                    SLIMCON                             *'
      write(*,*) ' *                                                        *' 
      write(*,*) ' * Converts simulation work file to model definition file *'
      write(*,*) ' *             Version 6.0 ( April 10, 2007)              *'
      write(*,*) ' *                                                        *'
      write(*,*) ' ----------------------------------------------------------'

10    write(*,*)' Enter the simulation file name (up to 46 characters)'
      write(*,*)' without any extension, or carriage return to end.'
      write(*,fmt='(a4)',advance='no')' => '
      fullnm='                                                  '       

!    Next lines are the only read

      read(*,fmt='(a46)') fullnm
      if(fullnm(1:12)=='           ')&
         stop 'No model definition file created.'

      print *,' Do you want to allow algebraic variables to be inputs'
      print *,'  and outputs of the same unit? (y/n) (default= no) '
      write(*,fmt='(a4)',advance='no')' => '
      read(*,fmt='(a1)') answer
      if (answer .eq. 'y' .or. answer .eq. 'Y') then       !10/30/97  cp
         algbraic = .true.
      else
         algbraic = .false.
      end if

!---------- Detect illegal filename -------------------------------------
!  This area may require customization for specific systems which allow
!  more or less different characters in a filename.

      do i=1,lname
       charct = ' '
       charct = fullnm(i:i)
       if(charct==' ') goto 30
       if((charct>='a'.and.charct<='z').or. &
         (charct>='A' .and. charct<='Z') .or. &
         (charct>='0'.and.charct<='9'.and.i/=1)) then
         cycle
       else
         write(*,4000) i
         4000 format(' Character #',i2,' in filename is illegal')
         goto 10
       end if
      end do

      i = lname + 1
30    fullnm(i:i) = delim
      if(delim==' ') i = i - 1

      fullnm(i+1:i+lext) = 'sim'

      open(luf,file=fullnm,status='old',iostat=ioerr)
      if(ioerr /= 0) then
        write(*,*) ' Cannot open the work file.'
        go to 10
      end if

      fullnm(i+1:i+lext) = 'dfn'
      open(lud,file=fullnm,iostat=ioerr)
      if(ioerr /= 0) then
        write(*,*) ' Cannot open the the model definition file.'
        go to 10
      end if

!     Attempt to open types data file in current directory. if not present,
!      open default types data file

      open(luty,file=typarf,status='old',iostat=ioerr)
      if(ioerr == 0) then
        return
      else
        open(luty,file=typarff,status='old',iostat=ioerr)
        if(ioerr /= 0) then
          write(*,7000) typarf, typarff
          7000 format(1x,'Cannot open ',a12,' or ',a50)
          stop ' No model definition file created.'
        end if
      end if
       
      return
      end subroutine fileop

!**********************************************************************

        subroutine typar(itype,nsaved,iude,nnin,nnout,nnpar,inpar)
!
! ---------------------------------------------------------------------
!
!   This routine provides information about the types.
!   If inpar is less than or equal to zero, the routine determines
!   nsaved, iude, nnin, nnout, nnpar, iodent, iident. If inpar is
!   greater than zero, the routine continues reading the types data
!   file at the same location as for the previous typar call and
!   determines the number of differential equations using the
!   differential equation modifiers.
!
!**********************************************************************
!
      use slim
      implicit none
      character(len=1)         :: star,pound
      character(len=6)         :: compar
      character(len=80)        :: line,blank
      logical                  :: logand,lsum
      integer                  :: itype,iude,nnin,nnout,nnpar,nsaved,&
                                  inpar,icheck,iargue,ipchek,i
      real                     :: ckvalu                             

      data  blank(1:40)  /'                                        '/
      data  blank(41:80) /'                                        '/
!
!     Star is used as a flag to locate the beginning of type info.
!     icheck is used to select the correct type from typar.dat
!
      if(inpar<=0) then
        rewind luty

10      read(luty,fmt='(a1)',iostat=ioerr) star
        if(ioerr > 0) then
          write(*,4000) itype
          4000 format(1x,'Read error on type number ',i6)
          return
        else if(ioerr < 0) then
          write(*,5000) itype
          5000 format(1x,'Cannot find entry for type',i6, &
          ' in type information file.')
          return
        end if
        if(star /= '*') then
          goto 10
        end if

        read(luty,*,iostat=ioerr) icheck
        if(ioerr > 0) then
          write(*,4000) itype
          return
        else if(ioerr < 0) then
          write(*,5000) itype
          return
        end if
        if(icheck /= itype) then
          goto 10
        end if

        read(luty,*,iostat=ioerr) nsaved,iude,nnin,nnout,nnpar
        if(ioerr > 0) then
          write(*,6000) itype
          6000 format(1x,'Read error for type parameters for type',i6)
          return
        else if(ioerr < 0) then
          write(*,7000) itype
          7000 format(1x,'Cannot find information on type ',i6, &
          ' at end of file')
          return
        end if

        do i=1,nnin
          read(luty,*,iostat=ioerr) iident(i)
          if(ioerr > 0) then
            write(*,6000) itype
            return
          else if(ioerr < 0) then
            write(*,7000) itype
            return
          end if
        end do

        read(luty,fmt='(a1)',iostat=ioerr) pound
          if(ioerr > 0) then
            write(*,6000) itype
            return
          else if(ioerr < 0) then
            write(*,7000) itype
            return
          end if

        do i=1,nnout
          read(luty,*,iostat=ioerr) iodent(i)
          if(ioerr > 0) then
            write(*,6000) itype
            return
          else if(ioerr < 0) then
            write(*,7000) itype
            return
          end if
        end do

!     ipchek is the variable used to choose the parameter
!     that must be checked to determine the proper number
!     of differential equations; iude
!     compar is the equality or inequality used to select
!     the correct if statement
!     ckvalu is the variable that defines the value that
!     must be checked with par(pcheck)
!     iargue is the variable used in the calculation of iude
!
!     Initialize logical 'and' check variables

      else

        logand = .false.

!     Read in differential equation number modifiers

        do i=1,10

!-ph  Initialise lsum inside do loop.

         lsum = .true.

!-ph  check for the end of de number modifiers
!-ph  read line into buffer and check for termination cahracter '#'

         line = blank
         read(luty,fmt='(a80)',iostat=ioerr) line
         if(ioerr /= 0) then
           exit ! go to 50
         end if

         if (line(1:1)=='#') then
           exit ! goto 50
         else
            read(line,*) ipchek,iargue,compar,ckvalu
         end if

         ipchek = ipchek + inpar - 1

!   Determine if there is a logical "and" of modifiers

         if(compar(4:6)=='and' .or. compar(4:6)=='AND')then
           compar(4:6) = '   '
           logand = .true.
         else
           logand = .false.
         end if

!   Determine if differential equation modifer is true

         if(compar(1:2) == 'lt' .or. compar(1:2) == 'LT' ) then
           if(.not.(par(ipchek) .lt. ckvalu)) lsum = .false.
         else if(compar(1:2) == 'gt' .or. compar(1:2) == 'GT') then
           if(.not.(par(ipchek) .gt. ckvalu)) lsum = .false.
         else if(compar(1:2) == 'le' .or. compar(1:2) == 'LE') then
           if(.not.(par(ipchek) .le. ckvalu)) lsum = .false.
         else if(compar(1:2) == 'ge' .or. compar(1:2) == 'GE') then
           if(.not.(par(ipchek) .ge. ckvalu)) lsum = .false.
         else if(compar(1:2) == 'eq' .or. compar(1:2) == 'EQ') then
           if(.not.(par(ipchek) == ckvalu)) lsum = .false.
         else if(compar(1:2) == 'ne' .or. compar(1:2) == 'NE') then
           if(.not.(par(ipchek) .ne. ckvalu)) lsum = .false.
         end if
         if(.not.(logand).and.lsum) iude = iude - iargue
        end do

50      if(iude<=0) iude = 0
        if(logand) write(*,3000) itype
        3000 format(1x,'Warning ->  too many ands in', &
        ' differential equation modifers, type=',i5)

      end if
      return
      end subroutine typar
!***********************************************************************

      subroutine report(number,max,index)

!***********************************************************************

      implicit none
      integer                              :: index,max,number,left
      real                                 :: pct 
      character(len=40),dimension(20)      :: descrp
      character(len=8),dimension(19)       :: name
      
      data descrp(1)  / 'Blocks in the simulation ...............'/
      data descrp(2)  / 'Time dependent boundary variables ......'/
      data descrp(3)  / 'Differential equations in the simulation'/
      data descrp(4)  / 'Unit parameters in the simulation ......'/
      data descrp(5)  / 'Saved variables in the simulation ......'/
      data descrp(6)  / 'Superblocks in the simulation ..........'/
      data descrp(7)  / 'State variables in the simulation ......'/
      data descrp(8)  / 'Units in the simulation ................'/
      data descrp(9)  / 'Blocks in the largest superblock .......'/
      data descrp(10) / 'Differential equations in one superblock'/
      data descrp(11) / 'Differential equations in one unit .....'/
      data descrp(12) / 'Inputs or outputs in a single block ....'/
      data descrp(13) / 'Boundary conditions in one superblock ..'/
      data descrp(14) / 'Inputs or outputs in a single unit .....'/
      data descrp(15) / 'Parameters in a single unit ............'/
      data descrp(16) / 'Reported variables in one superblock ...'/
      data descrp(17) / 'Simultaneous equations in a single block'/
      data descrp(18) / 'Simultaneous equations in one superblock'/
      data descrp(19) / 'Units in a single block ................'/
      data descrp(20) / '........................................'/

      data name(1)  /'maxblk ='/
      data name(2)  /'maxbnd ='/
      data name(3)  /'maxdeq ='/
      data name(4)  /'maxpar ='/
      data name(5)  /'maxsav ='/
      data name(6)  /'maxsbk ='/
      data name(7)  /'maxstv ='/
      data name(8)  /'maxunt ='/
      data name(9)  /'mblkis ='/
      data name(10) /'mdeqis ='/
      data name(11) /'mdeqiu ='/
      data name(12) /'minoib ='/
      data name(13) /'mbndis ='/
      data name(14) /'minoiu ='/
      data name(15) /'mpariu ='/
      data name(16) /'mrptis ='/
      data name(17) /'mseqib ='/
      data name(18) /'mseqis ='/
      data name(19) /'muntib ='/

      pct = (float(number) / float(max)) * 100.
      left = max - number
      if(left>2) then
        write(*,1000) number,descrp(index),name(index),max,pct
        1000  format(1x,10x,i4,1x,a40, 1x, a8,i4,' (',f5.1,'%)')
      else if(left<0) then
        write(*,2000) number,descrp(index),name(index),max,pct
        2000 format(1x,'! Error ->',i4,1x,a40, 1x, a8,i4,' (',f5.1,'%)')
      else
        write(*,3000) number,descrp(index),name(index),max,pct
        3000 format(1x,'Warning ->',i4,1x,a40, 1x, a8,i4,' (',f5.1,'%)')
      end if

      return
      end subroutine report
!***********************************************************************

      subroutine varchk(nblock,nsblok,nbound)

! ---------------------------------------------------------------------
!
!  This routine checks to see if any of the time-dependent boundary
!  variables are also solved simultaneously.
!
!***********************************************************************

      use hvacsim_par
      use slim
      implicit none
      integer                      :: nblock,nbound,nsblok,npnt,i,j,k,l

!     Check variables solved simultaneously in blocks

      do i = 1,nblock
       do j = 1,nsolve(i)
        do k = 1,nbound
         if(ibound(k)==isolve(i,j)) then
           do l = 1,8
             if(ibound(k)<=ndnt(l).and.ndnt(l)/=-1) exit
           end do
20         l = l - 1

!     The next line added by P.Haves to handle correctly case
!     where next variable cataegory is empty ( and 20 lines below )

           if ( ndnt(l) == -1) l = l - 1

           npnt = ibound(k) - ndnt(l)
           print 1000, npnt,l
         end if
        end do
       end do
      end do

!     Check variables solved simultaneously in superblocks

      do i = 1,nsblok
       do j = 1,nssolv(i)
        do k = 1,nbound
         if(ibound(k)==issolv(i,j)) then
           do l = 1,8
            if(ibound(k)<=ndnt(l).and.ndnt(l)/=-1) exit
           end do
70         l = l - 1

           if ( ndnt(l) == -1) l = l - 1

           npnt = ibound(k) - ndnt(l)
           write(*,1000) npnt,l
           1000 format(1x,'! Error ->',' variable',i5,' in category', &
           &i5,' is',' solved simultaneously and is ',/,11x, &
           &'a boundary variable')
         end if
        end do
       end do
      end do


      return
      end subroutine varchk 
!**********************************************************************

       subroutine outchk(ns,is,iu,inx)

! ---------------------------------------------------------------------
!
! This routine checks to see if any two or more outputs are identified
! as a single state variable.
!
!**********************************************************************

       use hvacsim_par
       use slim
       implicit none
       integer,dimension(maxstv)     :: is,iu,inx
       integer                       :: ierr,n,ns,iplus1,i,j
       
       ierr = 0
       n = ns-1
       do i=1,n
        iplus1 = i + 1
        do j=iplus1,ns
         if(is(i)==is(j)) then
          write(*,1000) inx(i),iu(i),inx(j),iu(j)
          1000 format(' Error:  output',i4,' of unit',i4, &
              &' is identical to output',i4,' of unit',i4)
          ierr = 1
         end if
        end do
       end do
       if(ierr/=0) then
         write(*,*) ' Set of equations may be overdetermined.'
       end if

       return
       end subroutine outchk
!**********************************************************************

       subroutine tdbvis(nsblok,nbound)

!**********************************************************************

      use hvacsim_par
      use slim
      implicit none
      integer               :: nsblok,nbound,ibnd,max,isblk,nbk,i,&
                               ib1,nn,j,iu1,ni,k


      do isblk=1,nsblok
       nbk = nsuper(isblk)
       ntdbvs(isblk) = 0
       do i=1,nbk
        ib1 = isuper(isblk,i)
        nn = nunits(ib1)

        do j=1,nn
         iu1 = iblock(ib1,j)
         ni = nin(iu1)
         do k=1,ni
          do ibnd = 1,nbound
           if(in(iu1,k)==ibound(ibnd)) then
            ntdbvs(isblk) = ntdbvs(isblk) + 1
           end if
          end do
         end do
        end do
       end do
      end do

      max = 0
      do isblk=1,nsblok
       if(ntdbvs(isblk)>max) then
        max = ntdbvs(isblk)
       end if
      end do

      call report(max,mbndis,13)

      return
      end subroutine tdbvis

