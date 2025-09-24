!***********************************************************************
!
!     CTF : Creation of a file for conduction transfer functions.
!
!         Main calculation routines were programmed
!            by G. N. Walton ,Ref. [1].
!         Front-end routines were written by C. Park, S. Lopez,
!            & D. Clark.
!
!     July 25, 1984. Cheol Park, N.B.S.
!
!     Last updated :   April 26, 1989 C.P.
!
!     Revised: Sept. 8, 1999 Cheol Park
!              Changed the format for reading thermal property data
!
!     Update: January 11, 2007 Cheol Park, NIST
!              Fortran 77 code converted into Fortran 90.
!
! ----------------------------------------------------------------------
!
!     INPUTS:
!       list:     Output control
!                 =0 --- minimal output
!                 =1 --- simple output
!                 =2 --- detailed output
!                 =3 --- root search
!       metric:   Units of material properties
!                 =1 --- Metric
!                 =0 --- English
!       nsec:     Time step (s)
!       nl:       Number of material layers in the wall construct
!                 (max(nl)=maxnl)
!       idlayr:   Identification index of the material layer in FILE 8
!                 (THERM2.DAT  ), which is the input file to BANKTP
!
!     INPUTS FROM FILE 10 (direct access file) = output file of BANKTP
!       mnames:   Mnemonic name of the material (max(words)=lenwrd)
!       length:   Thickness of the layer     (m)        (ft)
!       cond:     Conductivity of material   (W/m**2-K/m)
!                                                    (Btu/h-ft**2-F/ft)
!       dens:     Density of material        (kg/m**3)  (lb/ft**3)
!       spht:     Specific heat of material  (kJ/kg-K)  (Btu/lb-F)
!       rval:     Thermal resistance         (m**2-K/W  (ft**2-F-h/Btu)
!
!     OUTPUTS:
!       sec:      Time step                    (seconds)
!       istr:     ID number of construct
!       nctf:     Number of CTF terms          (max(nctf)=maxnrf)
!       nord:     Number of order              (max(NORD)=5)
!       uval:     Overall conductance          (W/m**2-K)
!       ctfx:     Conduction transfer function (W/m**2-K)
!       ctfy:     Conduction transfer function (W/m**2-K)
!       ctfz:     Conduction transfer function (W/m**2-K)
!       ctfq:     Flux referring CTF           (W/m**2)
!
!     OUTPUTS IN FILE 9 (CTFDATA.DAT):
!       istr,nctf,nord,uval,ctfx,ctfy,ctfz,ctfq
!
!-----------------------------------------------------------------------
!
!     NOMENCLATURE for the main routine:
!
!       aa           Lower bound of the interval for root search
!       acode        ASHRAE code for structural material
!       ad           Upper left element of the layer derivative matrix
!       al           Upper left element of the layer matrix
!       as           Upper left element of the total construct matrix
!       at           Upper left element of the total derivative matrix
!       bb           Upper bound of the interval for root search
!       bd           Upper right element of the layer derivative matrix
!       beta         Roots
!       bl           Upper right element of the layer matrix
!       bs           Upper right element of the total construct matrix
!       bt           Upper right element of the total derivative matrix
!       cc           Thermal capacitance (layer thickness * density
!                    * specific heat)
!       cd           Lower left element of the layer derivative matrix
!       cetype       Integer error type
!       cl           Lower left element of the layer matrix
!       cnd          Overall conductance of the construct
!       code         Code number of each layer of construction
!       cond         Conductivity of each layer of construct
!                    (from outside)
!       conv         Convergence criterion (=1.0D-4)
!       count        Number of error occurance of ETYPE=2
!       cs           Lower left element of the total construct matrix
!       ct           Lower left element of the total derivative matrix
!       ctfq         Flux conduction transfer functions
!       ctfx         Adjusted X (external) CTF
!       ctfy         Adjusted Y (cross) CTF
!       ctfz         Adjusted Z (internal) CTF
!       ctq          Flux conduction transfer functions
!       ctx          X (external) conduction transfer functions
!       cty          Y (cross) conduction transfer functions
!       ctz          Z (internal) conduction transfer functions
!       d            1.0-RATIO
!       dd           Lower right element of the layer derivative matrix
!       dens         Material density
!       dl           Lower right element of the total construct matrix
!       dt           Lower right element of total derivative matrix
!       ebbt         Common ratio
!
!       etype        Error type
!                    = 3 causes stop
!                    = 2 causes stop, if this error type occurs 31 times.
!
!       itr          Number of iteration of root finding
!
!       key          Flag integer
!                    = 1 for computation of roots
!                    = 2 for computation of exact response
!
!       length       Thickness of the layer
!       lenwrd       Maximum number of characters in a word
!
!       list         Output control
!                    = 0 for minimal output
!                    = 1 for simple output
!                    = 2 for detailed output
!                    = 3 for root search
!
!       maxnl        Maximum number of layers in a construct
!       maxnrf       Maximum number of response factors
!       maxnrt       Maximum number of roots for CTF calculation
!       maxstr       Maximum number of structures
!       maxtfh       Maximum number of temperature and flux histories
!       mnames       Names of materials
!
!       mode         Dump instruction
!                    = 1 for describing conductive layers
!                    = 2 for printing roots and ratios
!                    = 3 for printing all orders of CTF
!
!       msg          Error message
!       ncl          Number of conductive layers in a construct
!       nctf         Number of conduction transfer function coefficients
!       nl           Number of layers in a construct
!       nord         Number of order
!       nratio       Number of RATIO
!       nrf          Number of temperature CTF calculated
!       nrt          Number of roots computed
!       nsec         Sampling time in second (integer)
!       order        Order ( number of fulx coefficients) of CTF
!       pfc          Components of pulse functions
!       r            Component of heat flux
!       rc           Square-root value of (RR*CC)
!       repeat       If true, repeat root search with smaller SINC
!       rfc          Response factor components
!       rr           Thermal resistance of single layer slab
!       rtres        Root-related residue elements
!       rval         Thermal resistance of material
!       r1           Sum of root-related residues, RTRES(j,1)
!       r2           Sum of root-realted residues, RTRES(j,2)
!       r3           Sum of root related residues, RTRES(j,3)
!       sevr         Warning message
!       sinc         Root search increment
!       spht         Specific heat of material
!       sumx         Sum of X CTF to test convergence
!       sumy         Sum of Y CTF to test convergence
!       sumz         Sum of Z CTF to test convergence
!       tc           Temperature of center
!       tinc         Time increment for CTF calculation
!       tol          Tolerance in searching roots (=1.0D-12)
!
!       units        Units of material properties
!                    = 1 for metric units
!                    = 2 for English units
!
!       uval         Overall conductance of construct
!       worst        Indicator for worst case
!       wt           Weight of material per unit area
!       w1           Root
!       x1           Initial value for root finding
!       x2           Root in exact solution
!       zrres        Zero-related residue elements
!
! ----------------------------------------------------------------------
!
!     SUBPROGRAMS CALLED:
!       Front-end routines : BANKTP, THERMP
!       Main routines:       DER,DUMPRF,ERROR,ILLINI,INITRF,MATRIX,
!                            RFCOMP,SEARCH,ZERORE
!     REFERENCES:
!     [1]  Walton, G.N., "Thermal Analysis Research Program Reference
!                        Manual," Nat'l Bureau of Standards,
!                        NBSIR 83-2655, March 1983.
!     [2]  Kusuda, T., "NBSLD, the Computer Program for Heating and
!                        Cooling Loads in Buildings,"
!                        Nat'l Bureau of Standards, BSS 69, July 1976.
!
! **********************************************************************
    module precision1
                                                                                   
      implicit none                                                                
      save                                                                         
      integer,parameter  :: ps=selected_real_kind(6)  ! single precision           
      integer,parameter  :: pd=selected_real_kind(15) ! double precision           
                                                                                   
    end module precision1                                                          
                                                                                   
!    **********************************************************************        
                                                                                   
    module ctf_comm                                                                
                                                                                   
    use precision1                                                                 
    save                                                                           
                                                                                   
    integer,parameter    :: maxnl=10,maxnrt=60,maxnrf=20,&                         
                            maxord=5,maxstr=10,lenwrd=43                           
                                                                                   
    character(len=1)     :: answer                                                 
    character(len=3)     :: acodex
    character(len=126)   :: title                                                  
    character(len=lenwrd):: mnamesx                                                
    real(kind=ps)        :: lengthx,condx,densx,sphtx,rvalx,wtx,wt                 
    integer              :: idata,nrow_title=12
    integer              :: list
    integer              :: i,j,k,n,jmax,jj,ios                                    
    integer,dimension(maxnl)               :: idlayr                               
                                                                                   
    logical              :: repeat
    integer              :: order
    integer              :: ncl,nrt,nrf
    real(kind=pd)        :: cnd,tinc, sumrc, sinc,zero,&                           
                            half, one, two, three, six
    character(len=lenwrd),dimension(maxnl) :: mnames
    real(kind=ps),dimension(maxnl)         :: length,cond,spht,dens,rval
    real(kind=pd),dimension(0:5,-1:maxnrf) :: ctx,cty,ctz
    real(kind=pd),dimension(5)             :: ctq
    real(kind=pd),dimension(maxnl)         :: rr,rc
    real(kind=pd),dimension(maxnrt)        :: beta,zrres
    real(kind=pd),dimension(maxnrt,3)      :: rtres

    end module ctf_comm                                                            

! **********************************************************************

    program ctf
    use ctf_comm
    implicit none
    character (len=12)           :: fname,fname1='therm2.dat  ',&
                                    cr='            '

    open(10,access='direct',recl=126)

    print *,' Enter the name of the thermal property data file,'
    print *,' or carriage return for default name: ',fname1
    write(*,fmt='(a4)',advance='no') ' => '
    read (*,fmt='(a12)') fname
    if(fname==cr) fname=fname1
    open(8,file=fname,iostat=ios)

    if(ios/=0) then
      print *,' Error: thermal property data file'
      print *,' does not exist.'
      stop ' *** Execution terminated ***'
    endif

    big_loop: do
      print *,'***************************************************'
      print *,'                select your choice:'                         
      print *,'***************************************************'         
      print *,'Enter A, C, D, or E'                                         
      print *                                                               
      print *,'A : Add thermal property data to the construction'
      print *,'     materials database'
      print *,'C : Create a ctf data file'
      print *,'D : Display the contents of the database file'
      print *,'E : End'
      write(*,fmt='(a4)',advance='no') ' => '
      read (*,fmt='(a1)') answer
                                                                            
      if(answer=='A' .or. answer=='a' ) then                                
        call thermp                                                         
      elseif(answer=='C' .or. answer=='c') then
        call make_ctf                                                       
      elseif(answer=='D' .or. answer=='d') then
        rewind 8                                                            
        do j=1,nrow_title                                                   
           read(8,fmt='(a80)') title
           write(*,fmt='(a80)') title
        end do                                                              
        do
          read(8,1000,iostat=ios) idata,mnamesx,&                           
                  lengthx,condx,densx,sphtx,rvalx,wt,acodex                 
          1000 format(i3,2x,a43,f7.4,f7.3,f6.0,f9.3,f7.3,f8.1,4x,a3)
          if(ios/=0) then
            exit                                                            
          endif                                                             
          write(*,1000) idata,mnamesx,lengthx,condx,densx,&                 
                    sphtx,rvalx,wt,acodex                                   
        enddo                                                               
      elseif(answer=='E' .or. answer=='e' ) then
        close(10,status='delete')
        stop ' ********** End of ctf run **********'                        
      else
        print *,'Invalid option !'                                          
      endif                                                                 
    enddo big_loop
    close(10,status='delete')

    stop '----------- end of ctf run -------------------'
    end program ctf

! **********************************************************************

      subroutine make_ctf

! ----------------------------------------------------------------------
!
!     make_ctf :  To create a ctf data file for modsim
!
!     January 11, 2007 Cheol Park, NIST
!
! **********************************************************************

      use ctf_comm
      implicit none                      
      character (len=12)           :: fname,&
                                      fname2='ctfinput.dat',&
                                      fname3='hvacsim.ctf ',&
                                      cr='            '
      real,dimension(0:maxnrf)     :: ctfx,ctfy,ctfz
      real,dimension(5)            :: ctfq
      real                         :: sec,uval
      integer                      :: metric,nsec,istr,nctf,nord,nl

!     Open input and output files

      print *,' Enter the name of the ctf definition file,'
      print *,' or carriage return for default name: ',fname2
      write(*,fmt='(a4)',advance='no') ' => '
      read (*,500) fname
      if(fname==cr) fname=fname2
      open(12,file=fname,status='replace')

      print *,' Enter the name of the ctf output file,'
      print *,' or carriage return for default name: ',fname3
      write(*,fmt='(a4)',advance='no') ' => '
      read (*,500) fname
      if(fname==cr) fname=fname3
      open(9,file=fname,status='replace')

!     Read material layers from the direct access file.
!     THERM2.DAT supplied with HVACSIM+ is in metric units.

      print *,'What kind units are used for material properties?'
      print *
      print *,'Enter  1 for metric units'
      print *,'       2 for standard(English)units'
      write(*,fmt='(a4)',advance='no') ' => '
      read *, metric

!     Create a direct access temporary file.

      call banktp

!     Read input information interactively.

      print *,'What kind of output do you want?'
      print *,' 0 :   for a very simple output,'
      print *,' 1 :   for a less simple output,'
      print *,' 2 :   for detailed output   or'
      print *,' 3 :   for root search'
      write(*,fmt='(a4)',advance='no') ' => '
      read *, list

      print *,'What is the time interval for ctf calculation in s ?'
      write(*,fmt='(a4)',advance='no') ' => '
      read *, sec
      nsec=sec
      write(9,300) sec
      write(12,300) sec

      istr=1
10    uval=0
      nctf=0
      nord=0
      do i=1,maxnl
        idlayr(i)=0
      enddo
      print *,'This construct id number (istr) is  ',istr
      print *,'How many layers in this construct? (max. = ',maxnl,' )'
      write(*,fmt='(a4)',advance='no') ' => '
      read *, nl
      print *,'Enter the layer id numbers with most outer layer first'
      write(*,fmt='(a4)',advance='no') ' => '
      read *, (idlayr(i),i=1,nl)

!     Store input information to input file

      write(12,400) nl,(idlayr(i),i=1,nl)

      tinc=dble(real(nsec))

      do  i=1,nl
        jj=idlayr(i)                                                       
        read (10,rec=jj) jj,mnamesx,lengthx,condx,densx, sphtx,rvalx

        mnames(i) = mnamesx
        select case (metric)
            case (1)
              length(i)=lengthx
              cond(i)=condx        
              spht(i)=sphtx        
              dens(i)=densx        
              rval(i)=rvalx
            case (2)
              length(i)=lengthx*0.3048
              cond(i)=condx*1.729577206
              spht(i)=sphtx*4.184
              dens(i)=densx*16.01846337
              rval(i)=rvalx*0.1762280394
            case default
              length(i)=lengthx
              cond(i)=condx        
              spht(i)=sphtx        
              dens(i)=densx        
              rval(i)=rvalx
        end select
      end do

!     Compute conduction transfer functions.

      call initrf(nl,nctf,nord,uval,ctfx,ctfy,ctfz,ctfq)

!     Print results.

      print *, 'istr: ',istr
      print *, 'nctf: ',nctf,'  nord: ',nord
      print *, 'uval: ',uval
      print *, 'ctfx: ',(ctfx(i),i=0,nctf)
      print *, 'ctfy: ',(ctfy(i),i=0,nctf)
      print *, 'ctfz: ',(ctfz(i),i=0,nctf)
      print *, 'ctfq: ',(ctfq(i),i=1,nord)

!     Store ctf for each structure

      if(nord==0) then
        nord=1
        ctfq(1)=0.0
      endif
      write(9,200) istr,nctf,nord,uval
      write(9,300) (ctfx(i),i=0,nctf)
      write(9,300) (ctfy(i),i=0,nctf)
      write(9,300) (ctfz(i),i=0,nctf)
      write(9,300) (ctfq(i),i=1,nord)

      istr=istr+1
      if(istr<=maxstr) then
        print *,'Do you want to continue? <y/n>  '
        write(*,fmt='(a4)',advance='no') ' => '
        read 100, answer
        if(answer=='y'.or.answer=='Y') then
          go to 10
        else
          endfile 9
        endif
      endif

100   format(a1)
200   format(1x,3i4,g15.6)
300   format(5g15.6)
400   format(11i4)
500   format(a12)
600   format(i3,2x,a3,2x,a60,f8.4,f8.3,f7.0,f6.3,f7.3,f7.1)

      return
      end subroutine make_ctf

! **********************************************************************

      subroutine thermp

! ----------------------------------------------------------------------
!
!     THERMP :  To create a data bank for thermal properties of each
!               layer of construct.
!
!     July 25, 1984, Cheol Park, N.B.S.
!     Updated:   April 27, 1989 C.P.
!     Update: January 11, 2007 Cheol Park, NIST
!              Fortran 77 code converted into Fortran 90.
!
! **********************************************************************

      use ctf_comm
      implicit none                      
      integer            :: ndata,id,idata_end
      character(len=3)   :: acode = "xxx"

      print *,'************************************************'
      print *,'*                                              *'
      print *,'*   You are now adding thermal property data   *'
      print *,'*   into the existing file.                    *'
      print *,'*   Please enter your data with care.          *'
      print *,'*                                              *'
      print *,'************************************************'
      print *,' '

!     Read existing records in the file 8.

      rewind 8
      do j=1,nrow_title
         read(8,fmt='(a80)') title
      end do

      idata = 1
      do
        read(8,100,iostat=ios) id,mnamesx,lengthx,condx,densx,&
                           sphtx,rvalx,wt,acodex                      
        100 format(i3,2x,a43,f7.4,f7.3,f6.0,f9.3,f7.3,f8.1,4x,a3)     
        if(ios==0) then
          if(id/=0) then
            write(10,rec=idata) id,mnamesx,lengthx,condx,&
                                densx,sphtx,rvalx,wt,acodex
            print *,' Reading existing data number = ', idata
            idata=idata+1
          else
            exit
          endif
        elseif(ios<0) then    ! end of file
          exit
        elseif(ios>0) then    ! file reading error
          stop 'Error in reading thermal properties data'
        endif
      enddo
      idata=idata-1

!     Add new data.

20    ndata=idata+1
      print *,' Enter the name of the construction material'
      print *,' (up to 40 characters)'
      write(*,fmt='(a4)',advance='no') ' => '
      print *,'=> '
      read (*,fmt='(a43)') mnamesx
      print *,'   Thickness (m) =  '
      read *, lengthx
      print *,'   conductivity (w/m-K) =  '
      read *, condx
      print *,'   Density (kg/m**3) =  '
      read *, densx
      print *,'   Specific heat (kj/kg-K) =  '
      read *, sphtx
      print *,'   Thermal resistance (m**2-K/W) =  '
      read *, rvalx

      wt=lengthx*densx

!     Print entered values.

      print *, mnamesx
      print *,' idata=',ndata,' length=',lengthx,' cond=',condx
      print *,' dens=',densx,' spht=',sphtx,' rval=',rvalx
      print *,' wt =',wt,' acode=',acode

!     Write new data on the file

      write(10,rec=ndata) ndata,mnamesx,lengthx,condx,densx,&
                          sphtx,rvalx,wt,acode

      rewind 8
      do j=1,nrow_title
        read(8,fmt='(a80)') title
      enddo
      do  idata=1,ndata
        read(10,rec=idata) id,mnamesx,lengthx,condx,densx,&
                           sphtx,rvalx,wt,acode
        write(8,100) id,mnamesx,lengthx,condx,densx,sphtx,&
                           rvalx,wt,acode
      end do

      print *,'Do you want to continue?  <y/n>  '
      write(*,fmt='(a4)',advance='no') ' => '
      read (*,fmt='(a1)') answer
      if(answer=='y' .or. answer=='Y') then
        idata=ndata
        goto 20
      else
        print *,'-- End of addition ---'
      endif

      return
      end subroutine thermp
! **********************************************************************

      subroutine banktp

! ----------------------------------------------------------------------
!
!     BANKTP :  To create a data bank for thermal properties of each
!               layer of construct.
!
!     February 1979 Cheol Park, N.B.S.
!
!     Updated : April 27, 1989 C.P.
!     Update: January 11, 2007 Cheol Park, NIST
!              Fortran 77 code converted into Fortran 90.
!
! **********************************************************************

      use ctf_comm
      implicit none

      rewind 8
      do j=1,nrow_title
         read(8,fmt='(a80)') title
      end do

      do
        read(8,100,iostat=ios) idata,mnamesx,lengthx,condx,densx,&
                               sphtx,rvalx,wt,acodex
        100 format(i3,2x,a43,f7.4,f7.3,f6.0,f9.3,f7.3,f8.1,4x,a3)
        if(ios<0) then
          exit
        elseif(ios>0) then
          stop 'Error in reading thermal properties data'
        endif
        if(idata/=0) then
 !          write(*,100) idata,mnamesx,lengthx,condx,densx,&
 !                       sphtx,rvalx,wt,acodex
           write(10,rec=idata) idata,mnamesx,lengthx,condx,densx,&
                        sphtx,rvalx,wt,acodex
        endif
      enddo

      return
      end subroutine banktp
! **********************************************************************
!
      subroutine der
!
! ----------------------------------------------------------------------
!
!     DER : Calculate the total construct and total derivative
!           matrices for a single value of the Laplace parameter
!
!     Feb. 1984, George N. Walton, N.B.S.
!     Update: January 11, 2007 Cheol Park, NIST
!              Fortran 77 code converted into Fortran 90.
!
! **********************************************************************

      use ctf_comm
      implicit none                                                      
      real(kind=pd)    :: ad,bd,cd,dd,al,bl,cl,dl,as,bs,cs,ds,at,bt,ct,dt
      real(kind=pd)    :: es, et, c, s, p                                
      common /rf4/ ad,bd,cd,al,bl,cl,as,bs,cs,ds,at,bt,ct,dt             
      equivalence (ad,dd), (al,dl)                                       

!     Compute rtres for each root.

      do j=1,nrt
        as=one
        ds=one
        bs=zero
        cs=zero
        at=zero
        bt=zero

        do i=1,ncl

!       Define layer and layer derivative matrix elements

          if (rc(i)==zero) then

!       Case 1:  rc(i) is zero.

            al=one
            bl=rr(i)
            cl=zero
            ad=zero
            bd=zero
            cd=zero

!       Case 2:  rc(i) is non-zero.

          else
            p=beta(j)*rc(i)
            c=cos(p)
            s=sin(p)
            al=c
            bl=s*rr(i)/p
            cl=-s*p/rr(i)
            ad=s*rc(i)*half/beta(j)
            bd=(s/p-c)*rr(i)*half/(beta(j)*beta(j))
            cd=(rc(i)*rc(i)*c*half+ad)/rr(i)
          end if

!     Total derivative matrix calculation

          et=at*al+bt*cl + as*ad+bs*cd
          bt=at*bl+bt*dl + as*bd+bs*dd
          at=et

!     Total construct matrix calculation

          es=as*al+bs*cl
          bs=as*bl+bs*dl
          as=es
          es=cs*al+ds*cl
          ds=cs*bl+ds*dl
          cs=es
        enddo

!     Combine residue elements for each root.

        rtres(j,2)=one/(bt*beta(j)**4)
        rtres(j,1)=rtres(j,2)*ds
        rtres(j,3)=rtres(j,2)*as

      enddo

      return
      end subroutine der
! **********************************************************************
!
      subroutine dumprf(mode)
!
! ----------------------------------------------------------------------
!
!     DUMPRF : Dump conduction transfer calculation.
!
!       MODE=1:  Describe conductive layers.
!       MODE=2:  Print roots and ratios.
!       MODE=3:  Print all orders of conduction transfer functions.
!
!        DUMPRF initiated by RNOS=3.
!
!     Feb. 1984, George N. Walton, N.B.S.
!     Update: January 11, 2007 Cheol Park, NIST
!              Fortran 77 code converted into Fortran 90.
!
! **********************************************************************

      use ctf_comm
      implicit none
      integer, intent(in)                  :: mode
      character(len=10)                        :: mnam

      select case (mode)
        case(1)
          write(*,112)
          do n=1,ncl
            mnam = mnames(n)(1:10)
            write(*,113) mnam,length(n),cond(n),spht(n),&
               dens(n),rval(n),rc(n)
          enddo
          write(*,114) sumrc,cnd,tinc                          
        case(2)
          write(*,121) nrt,sinc
          do n=1,nrt                        
            write(*,123) n,beta(n)          
          enddo                             
        case(3)
          write(*,131) nrf,order
          write(*,134)                         
          do n=0,nrf                           
            write(*,133) n,(ctx(k,n),k=0,5)    
          enddo                                
          write(*,134)                         
          do n=0,nrf                           
            write(*,133) n,(cty(k,n),k=0,5)    
          enddo                                
          write(*,134)                         
          do n=0,nrf                           
            write(*,133) n,(ctz(k,n),k=0,5)    
          enddo                                
          write(*,134)                         
          write(*,133) order,(ctq(k),k=1,5)
        case default
      end select

  112 format(/,'  layer',t18,'l',8x,'k',7x,'cp',8x,'d',8x,'r',7x,'rc')
  113 format(1x,a10,8f9.3)
  114 format(/,' sumrc =',f10.5,'  cnd =',f10.5,'  tinc =',f10.2)
  121 format(/,'  number of roots =',i3,';  search increment =',f8.6,//,&
       '    n        root')
  123 format(i5,1x,2f14.9)
  131 format(/,'  number of response factors =',i3,';  order =',i3,//,&
        4x,'n',8x,'o 0',9x,'o 1',9x,'o 2',9x,'o 3',9x,'o 4',9x,'o 5')
  133 format(i5,2x,6f12.6)
  134 format(' ')

      return
      end subroutine dumprf
! **********************************************************************
!
      subroutine errorx(msg,cetype)
!
! ----------------------------------------------------------------------
!
!     ERRORX : Print error message and take action.
!
!     ETYPE=3 causes stop;  31st ETYPE=2 causes stop.
!     The error type is usually treated as a constant in the calling
!     program.
!
!     Feb. 1984, George N. Walton, N.B.S.
!     Update: January 11, 2007 Cheol Park, NIST
!              Fortran 77 code converted into Fortran 90.
!
! **********************************************************************

      implicit none
                                                                 
      integer           :: jlo,jhi,j,ibound                      
      integer           :: count=0, etype, worst=-1, cetype      
      character (len=*) :: msg                                   
      character (len=7),dimension(0:3) :: sevr                   
      data sevr /'       ','warning','severe ','fatal  '/        


!     Statement function: bounded integer

      ibound(jlo,j,jhi)=max(jlo,min(j,jhi))

      if(cetype>=0) then
        etype=ibound(0,cetype,3)
        worst=max(worst,etype)
        write(*,fmt='(a4,a7,2x,a)') ' ** ',sevr(etype),msg
        if(etype==3) stop 'fatal input error'
        if(etype==2) count=count+1
        if(count>30) stop 'too many severe errors'
      else
        cetype=worst
      end if

      return
      end subroutine errorx
! **********************************************************************
!
      subroutine illini(aa,ff,bb,gg,tol)
!
! ----------------------------------------------------------------------
!
!     ILLINI : Computation of root in the interval [AA,BB] using
!              modified false position method.
!     Convergent solution is found when the upper bound minus the
!     lower bound is less than TOL.  If a root is not found, Set
!     REPEAT to TRUE.
!
!     Feb. 1984, George N. Walton, N.B.S.
!     Update: January 11, 2007 Cheol Park, NIST
!              Fortran 77 code converted into Fortran 90.
!
! **********************************************************************

      use ctf_comm
      implicit none
      integer                       :: itr
      real(kind=pd),intent(in)      :: aa,ff,bb,gg                       
      real(kind=pd),intent(in)  :: tol                               
      real(kind=pd) :: ad,bd,cd,dd,al,bl,cl,dl,as,bs,cs,ds,at,bt,ct,dt   
      real(kind=pd) :: a, b, f, g, fw0, fw1, w1                          
      common /rf4/ ad,bd,cd,al,bl,cl,as,bs,cs,ds,at,bt,ct,dt             
      equivalence (ad,dd), (al,dl)                                       
      equivalence (fw1,bs)                                               

      a=aa
      b=bb
      f=ff
      g=gg
      fw1=ff
      itr=0

!     Root finder iteration

   10 continue
      fw0=fw1
      itr=itr+1

      if(itr>25) then

!     Repeat root search.

        repeat=.true.
        return
      endif

!     Approximate root.

      w1=b-g*(b-a)/(g-f)

!     Calculate new root approximation.

      call matrix(w1)

      if(list>=3) print *, 'i matrix: ',w1,bs

!     Separate the cases.

      if(f*fw1 < 0) then

!     Case 1: approximation on right side of root.

        b=w1
        g=fw1

!     Last two approximations on right side of root?

        if(fw0*fw1>zero) then
          f=half*f
          go to 40
        endif

      elseif(f*fw1 > 0) then

!     Case 2: approximation on left side of root.

        a=w1
        f=fw1

      elseif(f*fw1 == 0) then

        goto 50

      endif

!     Last two approximations on left side of root?

      if(fw0*fw1>zero) then
        g=half*g
      endif

!     Root found if tolerance met.

   40 continue
      if(abs(b-a)>tol) go to 10
      w1=half*(a+b)

!     Case 3: exact root found.

   50 continue
      nrt=nrt+1
      if(nrt<=maxnrt) then
        beta(nrt)=w1
      endif

      return
      end subroutine illini
! **********************************************************************
!
      subroutine initrf(nl,nctf,nord,uval,ctfx,ctfy,ctfz,ctfq)
!
! ----------------------------------------------------------------------
!
!     INITRF : Direct the calculation of conduction transfer functions
!
!     This is the root subroutine of the CTF segment.
!
!     Feb. 1984, George N. Walton, N.B.S.
!     Update: January 11, 2007 Cheol Park, NIST
!              Fortran 77 code converted into Fortran 90.
!
! **********************************************************************

      use ctf_comm
      implicit none
      integer, intent(in)                      :: nl
      integer, intent(out)                     :: nctf,nord
      real, intent(out)                        :: uval
      real, dimension(0:maxnrf),intent(out)    :: ctfx,ctfy,ctfz
      real, dimension(5)                       :: ctfq
      real(kind=pd)                            :: rsum
      real(kind=pd),dimension(maxnl)           :: cc
      integer                                  :: nqrtr
!
      zero=0.0d0
      half=0.5d0
      one=1.0d0
      two=2.0d0
      three=3.0d0
      six=6.0d0
      sumrc=zero
      rsum=zero
      ncl=0
      do i=1,nl
        ncl=ncl+1
        if(cond(i)<=0.1e-06) go to 10
        rr(ncl)=dble(length(i))/dble(cond(i))
        cc(ncl)=dble(length(i)*dens(i)*spht(i)*1000.)
        if(cc(ncl)>=0.1d-06) go to 20
   10   rr(ncl)=dble(rval(i))
        cc(ncl)=zero
        if(rr(ncl)>=0.1e-06) go to 20
        call errorx('insufficient data for response factors',2)
        call dumprf(1)
        return
   20   rsum=rsum+rr(ncl)
      enddo
      cnd=one/rsum
      if(uval<=0.0) uval=cnd

!     Adjust to specified u-value.

      rsum=zero
      do i=1,ncl
        rr(i)=rr(i)*cnd/uval
        rsum=rsum+rr(i)
        rc(i)=sqrt(rr(i)*cc(i))
        sumrc=sumrc+rc(i)
      enddo
      cnd=one/rsum

!     Compute conduction transfer functions (ctf) for multi-layer
!     construct.

      nqrtr=0

      if(list>=1) call dumprf(1)

      if(sumrc<=0.1d-06) go to 60
      sinc=min(six*six/tinc,two/sumrc)

!     Find roots for residue calculations.

   40 continue
        repeat=.false.
        call search

        if(list>=2) call dumprf(2)

        if(repeat) go to 50
        if(nrt<=0) go to 60

!     Calculation derivative & total construct matrices and
!     combine residue elements for each root.

        call der

!     Calculate zero residue elements.

        call zerore

!     Calculate response factor components, response factors, and
!     higher order conduction transfer functions.

        call rfcomp(nctf,nord,ctfx,ctfy,ctfz,ctfq)

        if(list>=1) call dumprf(3)

!     Check response factors.

        if(.not.repeat) return
   50   continue
        nqrtr=nqrtr+1
        if(nqrtr>3) go to 70

!     Recompute roots using smaller non-repetitive search increment.

        sinc=sinc*0.247
        go to 40

!     No roots found. quick responding surface. set ctf to u-value.

   60   continue
      nctf=0
      nord=0
      ctfx(0)=cnd
      ctfy(0)=cnd
      ctfz(0)=cnd

!     Run termination

   70 continue
      if(list<2) then
        call dumprf(1)
        call dumprf(2)
        call dumprf(3)
      endif
      call errorx('response factors not computed (respns)',2)

      return
      end subroutine initrf
! **********************************************************************
!
      subroutine matrix(w)
!
! ----------------------------------------------------------------------
!
!     MATRIX : Evaluate the conduction matrix for a multi-layered slab.
!
!     Feb. 1984, George N. Walton, N.B.S.
!     Update: January 11, 2007 Cheol Park, NIST
!              Fortran 77 code converted into Fortran 90.
!
! **********************************************************************

      use ctf_comm
      implicit none                                                         
                                                                            
      real(kind=pd),intent(in)             :: w                             
      real(kind=pd) :: ad,bd,cd,dd,al,bl,cl,dl,as,bs,cs,ds,at,bt,ct,dt      
      real(kind=pd) :: es, p, s                                             
      common /rf4/ ad,bd,cd,al,bl,cl,as,bs,cs,ds,at,bt,ct,dt                
      equivalence (ad,dd), (al,dl)                                          

      as=one
      bs=zero
      do i=1,ncl

!     Define layer matrix elements.

        if (rc(i)==zero) then

!     Case 1   rc(i) is zero.

          al=one
          bl=rr(i)
          cl=zero

!     Case 2 : rc(i) is non-zero.

        else
          p=w*rc(i)
          s=sin(p)
          al=cos(p)
          bl=s*rr(i)/p
          cl=-s*p/rr(i)
        end if

!     Total construct matrix calculation.

        es=as*al+bs*cl
        bs=as*bl+bs*dl
        as=es
      enddo

      return
      end subroutine matrix
! **********************************************************************
!
      subroutine rfcomp(nctf,nord,ctfx,ctfy,ctfz,ctfq)
!
! ----------------------------------------------------------------------
!
!     RFCOMP : Compute conduction transfer functions for
!              multi-layered slabs from the roots and residues.
!
!     Feb. 1984, George N. Walton, N.B.S.
!     Update: January 11, 2007 Cheol Park, NIST
!              Fortran 77 code converted into Fortran 90.
!
!**********************************************************************

      use ctf_comm
      implicit none                                                           
                                                                              
      integer, intent(out)                     :: nctf,nord                   
      real,dimension(0:maxnrf)                 :: ctfx,ctfy,ctfz              
      real,dimension(5)                        :: ctfq                        
      real(kind=pd),dimension(maxnrt)          :: pfc,ebbt                    
      real(kind=pd),dimension(5)               :: sumx,sumy,sumz,d,r          
      real(kind=pd),dimension(3,3)             :: rfc                         
      real(kind=pd)                            :: conv, r1, r2, r3, ti,&      
                                                  x, y, z                     
      integer                                  :: nratio                      

!     Initialization.

      ti=one/tinc
      conv=1.0d-4
      do j=1,3
        do i=1,3
          rfc(i,j)=zero
        enddo
      enddo
      nratio=min(nrt,5)
      do i=-1,maxnrf
        do k=0,5
          ctx(k,i)=zero
          cty(k,i)=zero
          ctz(k,i)=zero
        enddo
      enddo
      do j=1,nrt
        pfc(j)=one
        ebbt(j)=exp(-beta(j)*beta(j)*tinc)
      enddo
      d(1)=one-ebbt(1)
      do k=2,nratio
        d(k)=d(k-1)*(one-ebbt(k))
      enddo
      do k=1,nratio
        d(k)=one/d(k)
        sumx(k)=zero
        sumy(k)=zero
        sumz(k)=zero
      enddo
      ctx(0,0)=cnd
      cty(0,0)=cnd
      ctz(0,0)=cnd
      jmax=nrt

!     Compute conduction transfer functions until arrays filled
!     or convergence tests met.

      big_loop: do i=0,maxnrf

!     Compute root-related residues, and add root-related
!     and zero-related residues.

        r1=zero
        r2=zero
        r3=zero
        do j=jmax,1,-1
          pfc(j)=ebbt(j)*pfc(j)
          if(pfc(j)<1.e-30) jmax=j
          r1=r1+rtres(j,1)*pfc(j)
          r2=r2+rtres(j,2)*pfc(j)
          r3=r3+rtres(j,3)*pfc(j)
          if(list>=4) print *,'rfc ',jmax,r1,r2,r3
        enddo
        rfc(1,1)=(r1+zrres(1))*ti
        rfc(2,1)=(r2+zrres(2))*ti
        rfc(3,1)=(r3+zrres(3))*ti
        if(list>=3) print *,'zrs ',(zrres(j),j=1,3)

!     Calculate response factors.

        ctx(0,i)=ctx(0,i)+rfc(1,1)-two*rfc(1,2)+rfc(1,3)
        cty(0,i)=cty(0,i)+rfc(2,1)-two*rfc(2,2)+rfc(2,3)
        if(list>=4) print *,'cty ',i,(rfc(2,j),j=1,3),cty(0,i)
        ctz(0,i)=ctz(0,i)+rfc(3,1)-two*rfc(3,2)+rfc(3,3)
        if(cty(0,i)<-1.0e-09) repeat=.true.
        if(list>=3)print *,'rfc: ',i,rfc(1,1),ctx(0,i),rfc(2,1),cty(0,i)

!     Calculate higher order conduction transfer functions.

        do k=1,nratio
          ctx(k,i)=ctx(k-1,i)-ebbt(k)*ctx(k-1,i-1)
          cty(k,i)=cty(k-1,i)-ebbt(k)*cty(k-1,i-1)
          ctz(k,i)=ctz(k-1,i)-ebbt(k)*ctz(k-1,i-1)
          sumx(k)=sumx(k)+ctx(k,i)
          sumy(k)=sumy(k)+cty(k,i)
          sumz(k)=sumz(k)+ctz(k,i)
          if(abs((cnd-sumx(k)*d(k))/cnd)>conv) cycle
          if(abs((cnd-sumy(k)*d(k))/cnd)>conv  .and.&
             cty(k,i)>5.d-06)  cycle
          if(abs((cnd-sumz(k)*d(k))/cnd)>conv) cycle

!     Convergence tests met.  order is first k meeting all tests.

          order=k
          nrf=i
          go to 70
        enddo

!     Convergence tests not met.

        order=nratio
        nrf=i
        if(repeat) return

!     Time shift response factor components.

        do n=1,3
          rfc(n,3)=rfc(n,2)
          rfc(n,2)=rfc(n,1)
          rfc(n,1)=zero
        enddo
      enddo big_loop

!     Report failure to converge.

      call errorx('ctf calculation not convergent',2)
      i=nratio
      write(*,fmt=99) sumx(i)*d(i),sumy(i)*d(i),sumz(i)*d(i),cnd
   99 format(' estimated conductivities:  x=',f10.5,'  y=',f10.5,&
        '  z=',f10.5,/,' calculated conductivity =',f10.5)
      order=nratio
   70 continue

!     Adjust ctf for exact u-value.

      x=cnd/(sumx(order)*d(order))
      y=cnd/(sumy(order)*d(order))
      z=cnd/(sumz(order)*d(order))
      do i=0,nrf
        ctfx(i)=ctx(order,i)*x
        ctfy(i)=cty(order,i)*y
        ctfz(i)=ctz(order,i)*z
      enddo

!     Calculate flux transfer functions.

      do k=1,5
        r(k)=zero
      enddo
      do k=1,order
        r(k)=ebbt(k)
      enddo
      ctq(1)=r(1)+r(2)+r(3)+r(4)+r(5)
      ctq(2)=-(r(1)*(r(2)+r(3)+r(4)+r(5))+r(2)*(r(3)+r(4)+r(5))+&
             r(3)*(r(4)+r(5))+r(4)*r(5) )
      ctq(3)=r(1)*r(2)*(r(3)+r(4)+r(5))+r(1)*r(3)*(r(4)+r(5))+&
             r(1)*r(4)*r(5)+r(2)*r(3)*(r(4)+r(5))+r(4)*r(5)*(r(2)+r(3))
      ctq(4)=-(r(1)*r(2)*(r(3)*r(4)+r(3)*r(5)+r(4)*r(5))+&
             r(3)*r(4)*r(5)*(r(1)+r(2)) )
      ctq(5)=r(1)*r(2)*r(3)*r(4)*r(5)
      do k=1,order
        ctfq(k)=ctq(k)
      enddo
      nord=order
      nctf=nrf

      return
      end subroutine rfcomp
! **********************************************************************
!
      subroutine search
!
! ----------------------------------------------------------------------
!
!     SEARCH : Determine the upper and lower bounds within
!              which a root must exist. Three methods are used to
!              determine root bounding:
!     (1) Successive values of BS have opposite sign,
!     (2) Successive values of AS have opposite sign: 2 roots bound,
!     (3) Middle of 3 values is smallest (absolute value): 2 roots.
!
!     Feb. 1984, George N. Walton, N.B.S.
!     Update: January 11, 2007 Cheol Park, NIST
!              Fortran 77 code converted into Fortran 90.
!
! **********************************************************************

      use ctf_comm
      implicit none
      real(kind=pd) :: ad,bd,cd,dd,al,bl,cl,dl,as,bs,cs,ds,at,bt,ct,dt
      real(kind=pd) :: b0,b1,b2,b3,b4,bh,p0,p1,p2,p3,p4,ph,tol,ul,  &
                       a2,a3,a4,ah,pid
      common /rf4/ ad,bd,cd,al,bl,cl,as,bs,cs,ds,at,bt,ct,dt
      equivalence (ad,dd), (al,dl)

!     initialization.

      nrt=0
      ul=6.6d0/sqrt(tinc)
      if(ncl==1) go to 80
      tol=1.0d-12
      p0=1.0d-3*sinc
      call matrix(p0)
      if(list>=3) print *, 's matrix: ',p0,bs,as
      b0=bs
      p2=p0
      b2=b0
      a2=as
      n=0

!     Begin search loop.

   10 continue

!     Increment upper bound.

      ph=p2+sinc

!     Determine bs at ph.

      call matrix(ph)
      if(list>=3) print *, 's matrix: ',ph,bs,as
      bh=bs
      ah=as
      if(ah*a2<zero) then
        n=n+1
      endif

      if(b2*bh<=zero) then

!     Find root of b in bounded interval.

        call illini(p2,b2,ph,bh,ph*tol)
        n=0
        if(ah*as<zero) n=1
        p0=ph
        p2=ph
        b0=bh
        b2=bh
        a2=ah

      else if (n==2) then

!     Find two roots of b.

        p4=p2
        a4=a2
        p3=half*(p2+ph)

        do i=1,25
          call matrix(p3)
          if(list>=3) print *, 'a matrix: ',p3,bs,as
          b3=bs
          if(b3*b2<zero) exit
          a3=as
          p1=p3
          if(abs(a3-a4)<1.0e-8) then
            repeat=.true.
            return
          endif
          p3=(a3*p4-a4*p3)/(a3-a4)
          if(p3<=p2 .or. p3>=ph) then
            repeat=.true.
            return
          endif
          a4=a3
          p4=p1
        enddo

        call illini(p2,b2,p3,b3,p3*tol)
        call illini(p3,b3,ph,bh,ph*tol)
        n=0
        if(ah*as<zero) n=1
        p2=ph
        b2=bh
        p0=ph
        b0=bh
        a2=ah

      else if (abs(b2)<abs(b0) .and. abs(b2)<abs(bh)) then

!     Find two roots of b.

        p4=ph
        b4=bh
        do i=1,25
          p1=half*(p0+p2)
          call matrix(p1)
          if(list>=3) print *, 'b matrix: ',p1,bs
          b1=bs

          if(b1*b0<=zero) then
            call illini(p0,b0,p1,b1,p1*tol)
            call illini(p1,b1,p2,b2,p2*tol)
            go to 70
          endif

          p3=half*(p2+p4)
          call matrix(p3)
          if(list>=3) print *, 'b matrix: ',p3,bs
          b3=bs

          if(b3*b4<=zero) then
            call illini(p2,b2,p3,b3,p3*tol)
            call illini(p3,b3,p4,b4,p4*tol)
            go to 70
          endif

          if(abs(b1)<abs(b2)) then
            p4=p2
            b4=b2
            p2=p1
            b2=b1
          else if(abs(b3)<abs(b2)) then
            p0=p2
            b0=b2
            p2=p3
            b2=b3
          else
            p0=p1
            b0=b1
            p4=p3
            b4=b3
          end if
        enddo

        repeat=.true.
        return

   70   continue
        n=0
        if(ah*as<zero) n=1
        p2=ph
        b2=bh
        p0=ph
        b0=bh
        a2=ah

!     No root(s) found.

      else
        p0=p2
        b0=b2
        p2=ph
        b2=bh
        a2=ah
      end if

      if(repeat .or. nrt>=maxnrt) then
        return
      endif

      if(p2<ul) then
        go to 10
      else
        return
      endif

!     Roots for single layer construct.

   80 continue
      pid=4.d0*atan(one)
      if(rc(1)>zero) then
        do j=1,maxnrt
          beta(j)=j*pid/rc(1)                          
          if(list>=3) print *, 'direct: ',j,beta(j)    
          if(beta(j)>ul) go to 100                     
        enddo                                          
        call errorx('root search incomplete',2)        
  100   continue                                       
        if(beta(1)<=ul) then                           
          nrt=j-1                                      
        endif                                          
      endif

      return
      end subroutine search
! **********************************************************************
!
      subroutine zerore
!
! ----------------------------------------------------------------------
!
!     ZERORE : Calculate the residues related to the double
!              pole at zero.
!
!     Feb. 1984, George N. Walton, N.B.S.
!     Update: January 11, 2007 Cheol Park, NIST
!              Fortran 77 code converted into Fortran 90.
!
! **********************************************************************

      use ctf_comm
      implicit none                                                         
                                                                            
      real(kind=pd) :: ad,bd,cd,dd,al,bl,cl,dl,as,bs,cs,ds,at,bt,ct,dt      
      real(kind=pd) :: p                                                    
      common /rf4/ ad,bd,cd,al,bl,cl,as,bs,cs,ds,at,bt,ct,dt                
      equivalence (ad,dd), (al,dl)                                          

      bs=zero
      at=zero
      bt=zero
      ct=zero
      dt=zero

      do i=1,ncl

!     Define layer and layer derivative matrix elements.

        p=rc(i)*rc(i)
        bl=rr(i)
        ad=p*half
        bd=p*rr(i)/six
        cd=p/rr(i)

!     Calculation of total derivative matrix elements.

        bt=at*bl+bt+bd+bs*dd
        at=at+ad+bs*cd
        dt=ct*bl+dt+dd
        ct=ct+cd

!     Calculation of total construct matrix elements.

        bs=bs+bl
      enddo

!     Calculation of zero residues.

      zrres(2)=-cnd*cnd*bt
      zrres(1)=zrres(2)+cnd*dt
      zrres(3)=zrres(2)+cnd*at

      return
      end subroutine zerore
