!***********************************************************************
!
!     MODSIM   : A MODular SIMulation program
!                 Main program of HVACSIM+ package
!
!                         Version 20.0
!
!               National Institute of Standards and Technology
!               Building and Fire Research Laboratory
!               Building Environment Division
!               Gaithersburg, Maryland 20899-8631
!               U.S.A.
!
! ----------------------------------------------------------------------
!
!     Original MODSIM program written by  C. Ray Hill, NBS, June 1983
!
!     Modification into FORTRAN77 by Cheol Park, NBS, September 11, 1984
!
!     Modification for reading boundary values and initialization
!     file by  Dan Clark, NBS, Oct. 18, 1983 and May 1, 1984
!
!     Modification for building shell by Cheol Park, NBS, Dec. 27, 1984
!
!     Updated by Cheol Park & Dan Clark, May 6, 1985 & May 31, 1986
!
! ----- Version 5.0.2 February 5, 1987     Cheol Park
!
!       (1) MODSIM5, TYPE51, REPORT, INDATA, BAKDIF, and UPDATE
!           were modified.
!       (2) Parameters of MAXDEQ and MAXPAR were changed into
!           MAXDEQ=90, and MAXPAR=1200 to handle a larger simulation.
!
! ----  Version 5.2.3  August 3, 1988    Cheol Park
!
!       Modified for IBM P! or alike.
! 
! ----  Version 6.2.0  April 25, 1989
!   
!       (1) Parameter values were changed into:
!            MAXSBK=5, MDEQIS=30, MINOIB=70.
!       (2) I/O information code was slightly changed.
!       (3) All DO loops are indented.
! 
! ----  Version 6.2.2  September 18, 1989
!
!       Parameter MAXSBK=10.
!
! ----  Version 6.2.6  February 14, 1992
!
!       The label "RVPS" is changed into "OTHR", and "AHUM" into "HUMT".
!
! ----  Version 7.0  October 7, 1992
!
!       Use of INCLUDE statement, non-ANSI standard, in place of PARAMETERs
!
! ----  Version 8.0  March 7, 1994
!
!       TYPES.FOR was updated according to Y. Zhu's suggestions.
!
! ----  Version 8.1  February 10, 1995
!
!       MODINO.FOR was modified to read new format of model definition file,
!       which is the output of SLIMCON.  The data file has descriptions of
!       each data.
!
! ----  Version 10.0  July 11, 1996 P Haves, Loughborough Univ., UK
!       Supports the "simulation test-bed" produced by ASHRAE 825-RP 
!       (1) COMMON BLOCKS /INFOCB/ and /CHRON2/ added. /INFOCB/ 
!           facilitates warning messages 
!       (2) MAIN PROGRAM: write results to output file only if reporting
!           interval for current superblock is greater than zero. Avoids
!           frequent writes to output files when reporting interval is
!           inadvertantly set to zero
!       (3) SUBROUTINE OPNFIL: Long file names supported (46 chars + '.' + 3 
!           char extension = 50 chars)
!       (4) SUBROUTINE OPNFIL: unique extensions used for each file.
!           Lower case is used for convenience when file names are case sensitive.
!       (5) SUBROUTINE REPORT: avoid interpolation across a reset
!       (6) SUBROUTINE SUMMARY: NUMBER and LABEL dimensioned using PARAMETER
!           MSUMMARY included in modsim.par
!       (7) SUBROUTINE VLAB: ID dimensioned using PARAMETER
!           MSUMMARY included in modsim.par
!       (8) SUBROUTINE (SUPER)BLOCK: list active variables in five columns  
!           rather than eight to match the format of the intermediate solution  
!           and the Jacobian
!       (9) SUBROUTINE (SUPER)FNC: report intermediate residuals in addition  
!           to intermediate variable values
!                    
!       Jan. 16, 1997   C. Park, NIST
!       Added FRZ for selecting either freezing or unfreezing variables.
!
! ----  Version 20.0  July 2, 2007  Cheol Park, NIST
!       Updated to Fortran 90 code. 
!
!-----------------------------------------------------------------------
!  DISCLAIMER
!
!  This program is furnished by the government and used by any
!  recipient  with the express understanding that the United
!  States Government makes no warranty, expressed or implied,
!  concerning the accuracy, completeness, reliability, usability,
!  or suitability for any particular purpose of the information and data
!  contained in this program or furnished in connection therewith,
!  and the United States shall be under no liability whatsoever
!  to any person by reason of any use made thereof.  This program
!  belongs to the government.  Therefore, the recipient
!  further agrees not to assert any proprietary rights therein or
!  to represent this program to anyone as other than a government
!  program.
!
! ----------------------------------------------------------------------------
!
!   NOMENCLATURE:
!
!       a               Alpha vector [Ref.2]
!       ah              Alpha "hat" vector [Ref.2]
!       am              "A" matrix [Ref.2]
!       atolx           Absolute error tolerance (Used to calculate error
!                       criterion for differental equations and to determine
!                       the bounds for variable freezing.)
!       dnk             d-function [Ref.2]
!       fnk             f-function [Ref.2]
!       g               Gamma vector [Ref.2]
!       gh              Gamma "hat" vector [Ref.2]
!       iblock(b,i)     Unit number in BLOCK B
!       ibound(i)       Indices of state variables which are time dependent
!                       boundary conditions
!
!       ibstat(b)       Block status vector
!                       = 0   inactive block
!                       = 1   active block, all variables frozen
!                       = 2   active block, some variables not frozen
!
!       icall(s)        Index number of SUPERBLOCK S called at simulation time
!
!       icheck(i)       State variable status vector
!                       = -4  time independent boundary condition
!                       = -3  time dependent boundary condition
!                       = -2  frozen output  (not solved simultaneously)
!                       = -1  output (not solved simultaneously)
!                       =  0  solved simultaneously, may not be frozen
!                       =  1  solved simultaneously, may be frozen
!                       =  2  frozen
!                       =  3  unfrozen
!
!       icntrl          Indicator for time step rejection
!
!       idechk(d)       Differential equation status vector
!                       =  0  differential equation may not be frozen
!                       =  1  differential equation may be frozen
!
!       ideflg          Flag for controlling activity of differential equation
!       ident(i)        Reported variable identification vector
!       idevar(d)       Variable index for differential equation D
!       ifilen          Logical device number
!
!       ifzopt(s)       Superblock variable unfreezing option vector
!                       = 0  equations not recalculated
!                       = 1  equations recalculated; unfrozen variables added
!                       = 2  equations recalculated; all equations added
!
!       in(u,i)         Array of input connections for UNIT U
!       inb(b,i)        Array of input connections for BLOCK B
!       inde(u,i)       Index of the I-th DE in UNIT U
!       indesb(s,i)     Index of the I-th DE in SUPERBLOCK S
!       init            Flag for use of initialization file
!                       = 0  initialization from the initial state vector
!                       = 1  use of the initialization file
!       inp             Logical number of input device
!       insb(s,i)       Array of inputs to SUPERBLOCK S for input scanning
!
!       insopt(s)       Superblock input scan option vector
!                       = 0  inputs not scanned
!                       = 1  inputs scanned
!
!       iout(u,i)       Array of output connections for UNIT U
!       ioutb(b,i)      Array of output connections for BLOCK B
!       ip(i)           Permutation vector [Ref.2]
!       ipa(i)          Permutation vector for AM(I,J) [Ref.2]
!       iprint          Print-out selector
!       irej            Indicator for time step rejection
!       ireprt(i)       Indices of state variables reported after each
!                       time step
!       ireset          Reset counter
!       isaved(u)       Index of first saved workspace variable for UNIT U
!       isign           Indicator of the direction of state variable change
!       isolve(b,i)     Array of variables solved simultaneously within
!                       BLOCK B
!       ispos           Indices of state variables for viewing on screen
!       issolv(s,i)     Array of variables solved simultaneously within
!                       SUPERBLOCK S
!       isuper(s,i)     Array of BLOCK numbers in SUPERBLOCK S
!       isview          Index of SUPERBLOCK for viewing on screen
!       itime           Time step counter
!       iunits(u)       Vector of component types for UNIT U
!       jpar(u)         Index of first parameter for UNIT U
!       jsolve(b,i)     Array of variables solved simultaneously within
!                       BLOCK B
!       jssolv(s,i)     Array of variables solved simultaneously within
!                       SUPERBLOCK S
!       korder          Order of integration [Ref.2]
!       kstep           Number of steps since last attempt to change order
!       lbls            Labels
!       maxblk          Maximum blocks in the simulation
!       maxbnd          Maximum boundary variables in the simulation
!       maxdeq          Maximum differential equations in the simulation
!       maxlbl          Maximum labels
!       maxpar          Maximum unit parameters in the simulation
!       maxsav          Maximum saved variables in the simulation
!       maxsbk          Maximum superblocks in the simulation
!       maxstv          Maximum state variables in the simulation
!       maxunt          Maximum units in the simulation
!       mblkis          Maximum blocks in a superblock
!       mdeqis          Maximum differential equations in a superblock
!       mdeqiu          Maximum differential equations in a unit
!       minoib          Maximum inputs or outputs in a block
!       minois          Maximum inputs of outputs in a superblock
!       minoiu          Maximum inputs or outputs in a unit
!
!       mode            Indicator of mode of calculation method
!                       = 1  original forward difference method  [Ref.1]
!                       = 2  improved bi-directional method with TSTEP in FCN
!
!       mprtis          Maximum reported variables in a superblock
!       mseqib          Maximum simultaneous equations in a block
!       mseqis          Maximum simultaneous equations in a superblock
!       munitb          Maximum units in a block
!       nblock          Number of blocks in the simulation
!       nbound          Number of boundary conditions in the simulation
!       nd              Number of differential equations in the simulation
!       nde(u)          Number of differential equation in UNIT U
!       ndent(i)        State variable identification vector
!       ndsb(s)         Number of differential equation in SUPERBLOCK S
!       nin(u)          Number of inputs to UNIT U
!       ninb(b)         Number of inputs to BLOCK B
!       ninsb           Number of inputs to SUPERBLOCK S for input scanning
!       njsolv(b)       Number of simultaneous equations in BLOCK B
!       njsslv(s)       Number of simultaneous equations in SUPERBLOCK S
!       nout(u)         Number of outputs from UNIT U
!       noutb(b)        Number of outputs from BLOCK B
!       npar            Number of parameters in the simulation
!       nreprt(s)       Number of reported variables in SUPERBLOCK S
!       nsaved          Number of saved variables in the simulation
!       nsblok          Number of superblocks in the simulation
!       nsolve(b)       Number of simultaneous equations in BLOCK B
!       nssolv(s)       Number of simultaneous equations in SUPERBLOCK S
!       nstate          Number of state variables in the simulation
!       nstve!          Number of state variables for viewing on a screen
!       nsuper(s)       Number of blocks in SUPERBLOCK S
!       ntime           Number of calls of the simulation main routine
!       nu              Number of units in the simulation
!       nunits(b)       Number of units in BLOCK B
!       par(i)          Vector of parameters for all units
!       rtolx           Relative error tolerance (Used to calculate error
!                       criterion for differental equations and to determine
!                       the bounds for variable freezing.)
!       saved(i)        Saved workspace
!       sbtime(s)       Simulation time for SUPERBLOCK S
!       spast(d,i)      Buffer for saved past variables ( The saved variables
!                       are used to calculate backward difference.)
!       spred(d)        Buffer for saved past variables ( The saved variables
!                       are used to calculate predicted values for the next
!                       time step.)
!       state(i)        Vector of state variables
!       step(s)         Time step for SUPERBLOCK S
!       stold(i)        Vector of old state variables
!       time            Simulation time
!       title           Title in the simulation input data file
!       tmax            Maximum time step of the simulation
!       tmin            Minimum time step of the simulation
!       tprtof          Stopping time of writing diagonastic information
!       tprton          Start time of writing diagonastic information
!       treprt(s)       Time interval for writing reported state variables
!       treset          Reset time
!       tstate(i)       Vector of temporary state variables
!       tstep           Time step
!       tstop           Stopping time of the simulation
!       ttime           Time interval defining error control for differential
!                       equations [Ref.2]
!       ufzflg          Flag for unfreezing state variables
!       view            Indicator for viewing on screen
!       xtol            Error tolerance for the equation solver SNSQ [Ref.1]
!
! ----------------------------------------------------------------------------
!
!   DESCRIPTION OF MODULAR PROGRAMMING:
!
!       TYPE            A subroutine which models an individual component.
!
!                       NOTE: type modules which contain differential
!                       equations output the derivative of the corresponding
!                       variable. The variable is determined by solution of
!                       the (nonlinear) equation by subroutine SNSQ [Ref.1].
!
!                       NOTE: the outputs of type subroutines must be ordered
!                       so that derivatives are output in the first elements
!                       of the output vector.
!                       The TYPE subroutines have a call list of the form:
!
!                       SUBROUTINE TYPEn(XIN,OUT,PAR,SAVED,IOSTAT)
!
!                   --  XIN is a vector of inputs.
!                   --  OUT is a vector of outputs.
!                   --  PAR is a vector of parameters.
!                   --  SAVED is a saved workspace.
!                   --  IOSTAT is an input/output status vector. On entry
!                       this vector contains the status of the input
!                       variables (see ICHECK). On exit it contains flags
!                       which enable (IOSTAT=1) or disable (IOSTAT=0)
!                       variable freezing for the output variables.
!
!       UNIT            The implementation of a type in a simulation.
!                       A system may have several units of the same type;
!                       however, the simulation program has only one
!                       subroutine for each component type.
!                       NOTE: The unit numbers may be chosen arbitrarily
!                       with respect to the blocks and types; however,
!                       the unit numbers must be in sequence from 1 to
!                       NU, the total number of units in the simulation.
!
!       BLOCK           A set of units which can be viewed as comprising a
!                       superunit. Certain variables within a block are
!                       determined using a subroutine which solves sets of
!                       simultaneous nonlinear equations.
!
!       SUPERBLOCK      A set of blocks. Variables not solved within the
!                       blocks are solved simultaneously within a superblock.
!                       Each superblock is solved independently, but may be
!                       linked to other superblocks. Differential equations
!                       within each superblock are integrated independently
!                       so that the time step can vary among the superblocks.
!
!       STATE VECTOR    The set of state variables which define the state of
!                       the system being simulated.
!
!       VARIABLE STATUS A state variable can be: (1) a time independent
!                       boundary condition, (2) a time dependent boundary
!                       condition, (3) a superblock output, (4) active - i.e.
!                       solved simultaneously, or (5) frozen - i.e. removed
!                       from the set of active variables.
!
!       BLOCK STATUS    A block can be: (1) active or (2) inactive - all
!                       variables and inputs frozen or time independent
!                       boundary conditions.
!
!       VARIABLE        Information about the status of inputs is passed to
!       INTERLOCK       the type subroutines, and information is passed back
!                       from the types which determines whether or not the
!                       outputs may be frozen. This information flow can be
!                       used to establish interlocks to inhibit variable
!                       freezing, depending on the internal status of the
!                       types and status of the inputs.
!
!   VARIABLES IN MODEL DEFINITION FILE:  ( SB = superblock )
!
!           Variable Name    Format   Description
!           -------------    ------   -----------
!       1.  title            (a80)    Simulation title
!       2.  nstate,nsblok    (20i4)   # of state variables, # of SBs
!       3.  nsuper(s)        (20i4)   # of blocks in each SB
!       4.  state(i)         (5g15.6) Vector of state variable initial values
!       5.  ndent(i)         (20i4)   State variable identification vector
!       6.  nunits(b)        (20i4)   # of units in each block B
!       7.  njsslv(s)        (20i4)   # of simultaneous eqs in each SB
!       8.  njsolv(b)        (20i4)   # of simultaneqou eqs in each block
!       9.  isuper(s,i)      (20i4)   Array of block numbers in each SB
!       10. iblock(b,i)      (20i4)   Array of unit numbers in each block
!       11. iunits(u)        (20i4)   Array of component type #s for each unit
!       12. nin(u)           (20i4)   Number of inputs to unit U
!       13. in(u,i)          (20i4)   Array of input connections for unit U
!       14. nout(u)          (20i4)   Number of outputs from unit U
!       15. iout(u,i)        (20i4)   Array of output connections for unit U
!       16  jssolv(s,i)      (20i4)   Array of variables solved simultaneously
!                                     within each SB (between blocks)
!       17. jsolve(b,i)      (20i4)   Array of variables solved simultaneously
!                                     within each block
!       18. nde(u)           (20i4)   # of differential eqs in unit U
!       19. inde(u,i)        (20i4)   DE index for the Ith DE in unit U
!       20. idevar(d)        (20i4)   Variable index for DE #D
!       21. isaved(u)        (20i4)   Index of first saved variable for unit U
!       22. jpar(u)          (20i4)   Index of first parameter for unit U
!       23. npar,nsaved      (20i4)   Number of parameters & saved variables
!       24. par(i)           (5g15.6) Array of parameters for all units
!       25. nbound           (20i4)   # of time-dependent boundary variables
!       26. ibound(i)        (20i4)   State variable indices of boundary var's
!       27. nreprt(s)        (20i4)   # of reported variables in each SB
!       28  treprt(s)        (5g15.6) Reporting interval for each SB
!       29. ireprt(i)        (20i4)--+for-+Indices of reported variables
!       30. ident(s,i,1)     (20i4)  |each|Category # of reported variables
!       31. ident(s,i,2)     (20i4)--+SB--+Position in category of reported var
!       32. rtolx,atolx,              Error tolerances
!           xtol,ttime       (5g15.6)
!       33. ifzopt(s)        (20i4)   SB variable unfreezing option vector
!       34. insopt(s)        (20i4)   SB input scan option vector
!
!
!   SUBPROGRAMS CALLED:
!
!       asembl,bactiv,block,bounds,ecntrl,frzvar,indata,inscan,intliz,
!       iperm,opnfil,rconf,report,reset,restat,sumary,superb,unfrez,update
!
!      ------ RELATED TO THE BUILDING SHELL MODEL -----
!
!       rdenv
!
!   REFERENCES:
!
!      [1] NBS, "Guide to Avaiable Mathematical Software(GAMS),"
!          Center for Applied Math., National Bureau of Standards,
!          Washington, DC, Oct. 1981.
!          Information on subrotuine SNSQ can be found here.
!
!      [2] Brayton, R.K., Gustavson, F.G., and Hachtel, G.D.,
!          "A New Efficient Algorithms for Solving Differential-Algebraic
!          Systems Using Implicit Backward Differential Formulas,"
!          Proc. IEEE, Vol. 60, No. 1, Jan. 1972.
!          Information on the method used to integrate differential equations
!          can be found here.
!
! ***********************************************************************

    program modsim

    use modsim_head
    implicit none
    logical                   :: bldshl,bscall,intshl,view
    integer                   :: ecntrl,unfrez
    integer                   :: mout,nstvec,isview,isshel,i,j,&
                                 mmax,mtstop,msteps,iset,mreset,&
                                 ntime,ireset,id,ivar,istop,ncall,&
                                 isblk,ic,mfrez,nfrez,iblk,k,isw,&
                                 icntrl,i1,i2,i3,inscan
    integer,dimension(maxsbk) :: icall
    integer,dimension(10)     :: ispos
    real                      :: tmax,tstop,treset
    character(len=4),dimension(maxlbl):: label

      print *
      print *,'       *************************************************'        
      print *,'       *                                               *'        
      print *,'       *   MODSIM   : A MODular SIMulation program     *'
      print *,'       *      Main program of HVACSIM+ package         *'        
      print *,'       *                                               *'        
      print *,'       *                 version 20.8                  *'        
      print *,'       *                                               *'        
      print *,'       * National Institute of Standards & Technology  *'        
      print *,'       *   Gaithersburg, Maryland 20899-8631  U.S.A.   *'        
      print *,'       *                                               *'        
      print *,'       *************************************************'        
      print *                                                                   
                                                                                
!     Enter data interactively as requested.                                    
                                                                                
      call indata(bldshl,ispos,isshel,isview,nstvec,tmax,&
                     tstop,view,intshl,mout)
                                                                                
!     Read model definition data.                                               
                                                                                
      call rconf                         ! <<< 11/28/06
                                                                                
!     Assemble block information vectors.                                       
                                                                                
      call asembl                                                               
                                                                                
!     Write summary of the simulation configuration.                            
                                                                                
      call sumary(tmax,tstop)                                                   
                                                                                
!     Initialize simulation.                                                    
                                                                                
      if(nstvec>0) then                                                         
        do i=1,nstvec                                                           
          j=i+5                                                                 
          ispos(i)=ispos(j)+ndent(ispos(i))                                     
        end do                                                                  
      endif                                                                     
      tshell=0.                                                                 
      if(bldshl) read(ifile7,*) tshell                                          
      mmax=int(0.5+dble(tmax)/dble(tmin))                                       
      mtstop=int(0.5+dble(tstop)/dble(tmin))                                    
      msteps=int(0.5+dble(tshell)/dble(tmin))                                   
      iset=0                                                                    
      mreset=0                                                                  
      ntime=0                                                                   
      call intliz(treset,ireset)                                                
                      
      open(unit=43,file='HVS State.csv', ACTION='write', STATUS='replace')
      ! MAG 190911 print time label and header line
	  write(43, fmt="(a)", advance="no") "sim time, "
	  write(43, '(9999(i4, a))') (i, ",", i=1,nstate)

!     Begin simulation.                                                         
                                                                                
      print *,'   ------- simulation begins ---------'                          
                                                                                
!     If there are de's in the problem, save the variables in spred.            
                                                                                
10    if(nd>0) then                                                             
        do id=1,nd                                                              
          ivar=idevar(id)                                                       
          spred(id)=state(idevar(id))                                           
        end do                                                                  
      endif                                                                     
      itime=itime+1                                                             
      bscall=.false.                                                            
                                                                                
!     If intshl is .true., call only the superblock containing                  
!     the building shell after the first time step                              
                                                                                
      if(intshl .and. itime > 1) then                                           
        mstep(isshel)=int(dble(step(isshel))/dble(tmin))                        
        mtstep=msbt(isshel)+mstep(isshel)-mtime                                 
        mtime=mtime+mtstep                                                      
        if(mtime>=mtstop) istop=1                                               
        time=real(mtime)*tmin                                                   
        ncall=1                                                                 
        icall(1)=isshel                                                         
        bscall=.true.                                                           
        goto 55                                                                 
      endif                                                                     
                                                                                
!     Determine the next time and the superblock to be called.                  
                                                                                
      do isblk=1,nsblok                                                         
        mstep(isblk)=int(dble(step(isblk))/dble(tmin))                          
        if(mstep(isblk)>mmax .and. isblk/=isshel) mstep(isblk)=mmax             
      end do                                                                    
      mtstep=msbt(1)+mstep(1)-mtime                                             
      do isblk=1,nsblok                                                         
        mtstep=min(mtstep,msbt(isblk)+mstep(isblk)-mtime)                       
      end do                                                                    
      mtime=mtime+mtstep                                                        
                                                                                
!     At the end of simulation, call all superblocks.                           
                                                                                
      istop=0                                                                   
      if(mtime>=mtstop) then                                                    
        istop=1                                                                 
        mtime=mtstop                                                            
      endif                                                                     
                                                                                
      time=real(mtime)*tmin                                                     
      ncall=0                                                                   
      do isblk=1,nsblok                                                         
        if(mtime>=msbt(isblk)+mstep(isblk) .or. istop/=0) then                  
          mstep(isblk)=mtime-msbt(isblk)                                        
          step(isblk)=real(mstep(isblk))*tmin                                   
          ncall=ncall+1                                                         
          icall(ncall)=isblk                                                    
        endif                                                                   
      end do                                                                    
                                                                                
!     If there are time dependent boundary conditions, call bounds.             
                                                                                
55    if((.not.intshl) .and. nbound>0) then                                     
        call bounds(treset,ireset)                                              
                                                                                
!       Test for a reset condition.                                             
                                                                                
        if(ireset/=0) then                                                      
                                                                                
!-ph    Reset even if there are no de's in order to prevent invalid             
!-ph    interpolation when calculating outputs                                  
!-ph        if(nd>0) then                                                       
                                                                                
!         If there is a reset condition, integrate up to                        
!         the reset time and then reset the integration algorithm.              
                                                                                
          if(ireset>1) then                                                     
            mtime=mreset+1                                                      
            time=real(mtime)*tmin                                               
            tstep=tmin                                                          
            ireset=0                                                            
            treset=0.                                                           
            ncall=0                                                             
            do isblk=1,nsblok                                                   
              if(isblk==isshel .and. (.not. bscall)) cycle                      
              msbt(isblk)=mtime-1                                               
              mstep(isblk)=1                                                    
              step(isblk)=tmin                                                  
              call reset(0,isblk)                                               
              icall(isblk)=isblk                                                
              ncall=ncall+1                                                     
            end do                                                              
          else                                                                  
            mtime=mtime-mtstep                                                  
!-ph    Rounding down can lead to reset condition being ignored (28.3.89)       
!-ph           mreset=int(dble(treset)/dble(tmin))                              
            mreset=int(0.5+dble(treset)/dble(tmin))                             
            mtstep=mreset-mtime                                                 
            do i=1,nblock                                                       
              if(ibstat(i)==0) ibstat(i)=1                                      
            end do                                                              
            if(mtstep<1) then                                                   
              mtstep=1                                                          
              ireset=0                                                          
              treset=0.                                                         
            endif                                                               
            mtime=mtime+mtstep                                                  
            time=real(mtime)*tmin                                               
            ncall=0                                                             
            do isblk=1,nsblok                                                   
              if(isblk==isshel .and. (.not. bscall)) cycle                      
              mstep(isblk)=mreset-msbt(isblk)                                   
              step(isblk)=real(mstep(isblk))*tmin                               
              icall(isblk)=isblk                                                
              ncall=ncall+1                                                     
            end do                                                              
          endif                                                                 
                                                                                
!        If there are no de's, reset the flags when ireset=2.                   
                                                                                
!-ph    Case of no de's already treated (see above)                             
!-ph           if there are no de's, reset the flags when ireset=2.             
!-ph                                                                            
!-ph        elseif(ireset>1) then                                               
!-ph          treset=0.                                                         
!-ph          ireset=0                                                          
!-ph        endif                                                               
        endif                                                                   
      endif                                                                     
                                                                                
90    if(bldshl) call rdenv(isshel)                                             
                                                                                
!     Call all superblocks scheduled to be called at this time.                 
                                                                                
      big_loop: do ic=1,ncall                                                   
                                                                                
      isblk=icall(ic)                                                           
100   tstep=step(isblk)                                                         
                                                                                
!     If there are de's in the simulation, update the integration               
!     coefficients.                                                             
                                                                                
      if((ndsb(isblk)>0).and.(ideflg(isblk)==0)) then                           
        call update(isblk)                                                      
      endif                                                                     
                                                                                
!     Call superblock isblk.                                                    
                                                                                
      call superb(isblk)                                                        
                                                                                
!     Check for unfrozen variables.                                             
                                                                                
      mfrez=0                                                                   
      nfrez=0                                                                   
      do i=1,nsuper(isblk)                                                      
        iblk=isuper(isblk,i)                                                    
                                                                                
!       Subroutine unfrez checks for unfrozen variables.                        
                                                                                
        j=unfrez(isblk,iblk)                                                    
        ufzflg(i)=j                                                             
        if(abs(j)==1) nfrez=1                                                   
        if(j<0) mfrez=1                                                         
      end do                                                                    
                                                                                
!     If variables were frozen, skip the following block.                       
                                                                                
      if(mfrez+nfrez/=0) then                                                   
                                                                                
!       If no superblock equations were unfrozen or ifzopt is set to 0,         
!       service unfrozen block equations only.                                  
                                                                                
        if(.not.((mfrez==0).or.(ifzopt(isblk)==0))) then                        
                                                                                
!         Reset the state vector.                                               
                                                                                
          do i=1,nsuper(isblk)                                                  
            iblk=isuper(isblk,i)                                                
            call restat(1,iblk,isblk)                                           
          end do                                                                
                                                                                
!         If ifzopt is set to 2, unfreeze all superblock equations.             
                                                                                
          if(ifzopt(isblk)==2) then                                             
                                                                                
!           Set all superblock variables to unfrozen status.                    
                                                                                
            nssolv(isblk)=njsslv(isblk)                                         
            do i=1,njsslv(isblk)                                                
              ivar=jssolv(isblk,i)                                              
              issolv(isblk,i)=ivar                                              
              icheck(ivar)=3                                                    
            end do                                                              
          endif                                                                 
                                                                                
!         Determine block activity status and call the superblock,              
!         again to service the unfrozen superblock variables.                   
                                                                                
          call bactiv(isblk)                                                    
          call superb(isblk)                                                    
                                                                                
!         Check for unfrozen block variables.                                   
                                                                                
          nfrez=0                                                               
          do i=1,nsuper(isblk)                                                  
            iblk=isuper(isblk,i)                                                
            j=unfrez(isblk,iblk)                                                
            ufzflg(i)=j                                                         
            if(abs(j)==1) nfrez=1                                               
          end do                                                                
        endif                                                                   
        if(nfrez/=0) then                                                       
                                                                                
!         Service unfrozen block variables.                                     
                                                                                
          call bactiv(isblk)                                                    
          do i=1,nsuper(isblk)                                                  
            if(abs(ufzflg(i))==1) then                                          
              iblk=isuper(isblk,i)                                              
              call restat(0,iblk,isblk)                                         
              call block(isblk,iblk,iset)                                       
                                                                                
!             If additional block variables unfreeze, set all block             
!             variables to unfrozen status and call the block again.            
                                                                                
              j=unfrez(isblk,iblk)                                              
              if(abs(j)==1) then                                                
                nsolve(iblk)=njsolv(iblk)                                       
                do k=1,njsolv(iblk)                                             
                  ivar=jsolve(iblk,k)                                           
                  isolve(iblk,k)=ivar                                           
                  icheck(ivar)=3                                                
                end do                                                          
                call restat(0,iblk,isblk)                                       
                call block(isblk,iblk,iset)                                     
              endif                                                             
            endif                                                               
          end do                                                                
        endif                                                                   
      endif                                                                     
                                                                                
!     Replace state vector.                                                     
                                                                                
!     Begin loops to check all outputs of the superblock.                       
                                                                                
      do i=1,nsuper(isblk)                                                      
        iblk=isuper(isblk,i)                                                    
        do j=1,noutb(iblk)                                                      
          ivar=ioutb(iblk,j)                                                    
                                                                                
!       If the output is frozen or is time independent, skip it.                
                                                                                
          if(ivar>0.and.icheck(ivar)>=-1.and.icheck(ivar)/=2)then               
            state(ivar)=tstate(ivar)                                            
          endif                                                                 
        end do                                                                  
      end do                                                                    
                                                                                
!     Perform error control and updating for de's.                              
!     if there is no de's in superblock, set time-step to tmax.                 
                                                                                
      isw=1                                                                     
      if(ndsb(isblk)<=0) then                                                   
        step(isblk)=tmax                                                        
      else                                                                      
        ideflg(isblk)=1                                                         
        do i=1,ndsb(isblk)                                                      
          id=indesb(isblk,i)                                                    
          ivar=idevar(id)                                                       
          if(icheck(ivar)/=2) then                                              
            isw=2                                                               
            goto 220                                                            
          endif                                                                 
        end do                                                                  
        call reset(0,isblk)                                                     
      endif                                                                     
                                                                                
!     If all de's in superblock are frozen, double the time-step                
!     but don't exceed tmax.                                                    
                                                                                
      step(isblk)=min(2.*step(isblk),tmax)                                      
                                                                                
!     Calculate new time step and integration order.                            
                                                                                
220   if(isw==2) then                                                           
        ideflg(isblk)=0                                                         
        icntrl=ecntrl(isblk)                                                    
        if((icntrl==0).or.(mstep(isblk)==1).or.isblk==isshel) then              
                                                                                
!         Save backward differences.                                            
                                                                                
          do i=1,ndsb(isblk)                                                    
            id=indesb(isblk,i)                                                  
            i1=ip(isblk,1)                                                      
            ivar=idevar(id)                                                     
            if((kstep(isblk)>=korder(isblk)+1)&
              .and.(korder(isblk)<6)) then
              i2=korder(isblk)+2                                                
              i3=ip(isblk,i2-1)                                                 
              spast(id,i2)=spast(id,i3)                                         
            endif                                                               
            spast(id,i1)=state(ivar)-spred(id)                                  
          end do                                                                
          call iperm(isblk)                                                     
        else                                                                    
                                                                                
!         If the step is rejected, reset the integration algorithm.             
                                                                                
          do i=1,nsuper(isblk)                                                  
            iblk=isuper(isblk,i)                                                
            call restat(1,iblk,isblk)                                           
          end do                                                                
          mstep(isblk)=int(dble(step(isblk))/dble(tmin))                        
          if(mstep(isblk)<=0) mstep(isblk)=1                                    
          step(isblk)=real(mstep(isblk))*tmin                                   
          if(view) print *,' superblock ',isblk,&
            ':  time step rejected:  pstep = ',step(isblk)
          call reset(1,isblk)                                                   
          time=real(msbt(isblk)+mstep(isblk))*tmin                              
          if((.not.intshl) .and. nbound>0) then                                 
            call bounds(treset,ireset)                                          
          endif                                                                 
          if(bldshl) call rdenv(isshel)                                         
          goto 100                                                              
        endif                                                                   
      endif                                                                     
                                                                                
!     Write report.                                                             
                                                                                
      if (treprt(isblk)>0.0) then
        call report(tstop,ireset,isblk,mout,label)
      endif
                                                                                
      if(view .and. ((isblk==isview).or.(isview==0))) then                      
        ntime=ntime+1                                                           
        print *,'sb ',isblk,':','  time =',time,'   ntime =',ntime,&
                '   tstep =',tstep,'   pstep=',step(isblk),&
                '   order = ',korder(isblk)
        if(nstvec>0) then                                                       
          write(*,fmt='(1p5g15.6)')(state(ispos(i)),i=1,nstvec)                 
        endif                                                                   
      endif                                                                     
                                                                                
!     If bldg shell is used, set time step of bldg shell superblock             
                                                                                
      if(isblk==isshel) then                                                    
        step(isblk)=tshell                                                      
        if(itime==1 .and. msteps>1) step(isblk)=tshell-tmin                     
      endif                                                                     
                                                                                
!     Freeze variables if possible.                                             
                                                                                
      if(frz .and. itime>1) then         ! 1/16/97                              
                                                                                
        call frzvar(isblk)                                                      
                                                                                
!       Determine block activity status.                                        
                                                                                
        call bactiv(isblk)                                                      
      endif                                                                     
                                                                                
!     Increment superblock time counters.                                       
                                                                                
      msbt(isblk)=msbt(isblk)+mstep(isblk)                                      
      sbtime(isblk)=real(msbt(isblk))*tmin                                      
      if(msbt(isblk)/=mtime) then                                               
        mstep(isblk)=int(dble(step(isblk))/dble(tmin))                          
        if(mstep(isblk)>mmax) then                                              
          mstep(isblk)=mmax                                                     
        endif                                                                   
        if(msbt(isblk)+mstep(isblk)>mtime) then                                 
          mstep(isblk)=mtime-msbt(isblk)                                        
        endif                                                                   
        step(isblk)=real(mstep(isblk))*tmin                                     
        time=real(msbt(isblk)+mstep(isblk))*tmin                                
        do i=1,ndsb(isblk)                                                      
          id=indesb(isblk,i)                                                    
          ivar=idevar(id)                                                       
          spred(id)=state(ivar)                                                 
        end do                                                                  
        do i=1,nsuper(isblk)                                                    
          iblk=isuper(isblk,i)                                                  
          do j=1,noutb(iblk)                                                    
            ivar=ioutb(iblk,j)                                                  
            if(ivar>0) then                                                     
              stold(ivar)=state(ivar)                                           
              tstate(ivar)=state(ivar)                                          
            endif                                                               
          end do                                                                
        end do                                                                  
        if((.not.intshl) .and. nbound>0) then                                   
          call bounds(treset,ireset)                                            
        endif                                                                   
        if(bldshl) call rdenv(isshel)                                           
        goto 100                                                                
      endif                                                                     
      end do big_loop                                                           
                                                                                
!     Scan all superblock inputs for changes larger than the error              
!     tolerance.  if the tolerance is exceeded, call the superblock.            
                                                                                
      ncall=0                                                                   
      do isblk=1,nsblok                                                         
                                                                                
!       Skip the superblock containing the building shell                       
                                                                                
        if(isblk==isshel) cycle                                                 
                                                                                
!       If the superblock input scan option is set to 0, skip the               
!       superblock.  subroutine inscan scans superblock inputs.                 
                                                                                
        if((msbt(isblk)/=mtime).and.(insopt(isblk)>0)&
          .and.(inscan(isblk)>0)) then
          mstep(isblk)=mtime-msbt(isblk)                                        
          step(isblk)=real(mstep(isblk))*tmin                                   
          ncall=ncall+1                                                         
          icall(ncall)=isblk                                                    
        endif                                                                   
      end do                                                                    
      if((.not.intshl) .and. ncall>0) goto 90                                   
      do i=1,nstate                                                             
        if(state(i)<stold(i)) isign(i)=-1                                       
        if(state(i)>stold(i)) isign(i)=1                                        
        stold(i)=state(i)                                                       
        tstate(i)=state(i)                                                      
      end do    
      
   	  write(43, fmt="(f7.0)", advance="no") time
	  write(43, fmt="(a)", advance="no") ", "
	  write(43, '(9999(e16.8, a))') (state(i), ",", i=1,nstate)

      if(time<(tstop-0.5*tmin)) then                                            
        goto 10                                                                 
      else                                                                      
                                                                                
!       Save simulation results as an initialization file for the next          
!        run.                                                                   
                                                                                
        write(ifile6,fmt='(4f10.2,2i5)')&
                          tmin,tmax,tstop,time,nstate,nsaved
        write(ifile6,fmt='(1p5g15.6)') (state(i),i=1,nstate)                    
        write(ifile6,fmt='(25i3)') (isign(i),i=1,nstate)                        
        if(nsaved>0) then                                                       
          write(ifile6,fmt='(1p5g15.6)') (saved(i),i=1,nsaved)                  
        endif                                                                   
        print *,' --Final state file has been written -----'            
      endif                                                                     
                                                                                
      stop '   ------------End of simulation ---------'                         
      end program modsim                                                        
                                                                                

