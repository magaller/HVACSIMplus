!=======================================================================
! HVACGEN SOURCE FILE #2 OF 6, VERSION 5.0
!=======================================================================

! hvacgen2.f90

! This source code, when compiled with a FORTRAN 90/95 compiler, and linked
! with the compiled source code from the other five source code files,
! will produce the program HVACGEN which is the "front-end" program for
! the HVACSIM+ dynamic simulation package for building heating,
! ventilation, and air conditions systems and controls, plus other
! building systems. This program produces simulation work files which
! must be processed by the program SLIMCON to produce a model definition
! file which describes the structure and characteristics of the system
! to be simulated by the main simulation program, MODSIM.

!=======================================================================

!       CREATE MODULE

!------ This module transfers control to the proper subroutines
!       as indicated by the user's selection of creating a SImulation,
!       BLock,UNit.
!------ CHOICE is the variable which indicates what the user
!       would like to create, i.e., a SImulation, BLock, Unit.
!=======================================================================

SUBROUTINE create(crtflg,edflg,abort)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=2) :: choice

  100   WRITE(*,1)
  1     FORMAT(//1X,'Create a:',//, ' SImulation',/, ' BLock',/, ' UNit')
  messag=' '
  CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    crtflg=0
    blkflg=0
    simflg=0
    supflg=0
    abort = .false.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 100
  END IF

  choice=item
!!--- if(choice == 'SU' .or. choice == 'su') CALL CRSUP(CRTFLG,EDFLG)
IF (choice == 'BL'.OR.choice == 'bl') THEN
  CALL crblk(crtflg,edflg,rtn_crblk )
  if(rtn_crblk == 2) then
    crtflg=0
    blkflg=0
    simflg=0
    supflg=0
    abort = .true.
    return
  endif
ELSE IF (choice == 'UN'.OR.choice == 'un') THEN
  CALL crunit(crtflg,edflg,rtn_crunit,STATUS)
  if(rtn_crunit == 2) then
    crtflg=0
    blkflg=0
    simflg=0
    supflg=0
    abort = .true.
    return
  endif
ELSE IF (choice == 'SI'.OR.choice == 'si') THEN
  CALL crsim(crtflg,edflg,rtn_crsim)
  if(rtn_crsim == 2) then
    crtflg=0
    blkflg=0
    simflg=0
    supflg=0
    abort = .true.
    return
  endif
END IF
crtflg=0
blkflg=0
simflg=0
supflg=0
abort = .false.
RETURN
END SUBROUTINE create
!=======================================================================

!       SUBROUTINE CRUNIT

!------ This subroutine is called by the CREATE module and this
!       subroutine is used to "create" a UNit to be used in a work
!       file.
!------ UNIT is assigned to the unit number as inputted by the user.
!------ TYPE is assigned to the type number that is entered by
!       the user.
!------ STAR is a flag used to indicate the beginning of a type as
!       read from the types data file TYPAR.DAT .
!------ ICHECK is used to define the type that is read in from the data
!       file and is used in comparing whether the type matches that
!       specified by TYPE.
!------ DSCRIB is a character variable assigned to the description
!       of the particular type.
!------ NSAVED indicates the number of saved variables used in a
!       particular type. This is not important to the subroutine.
!------ IUDE is the number of differential equations used in the type.
!------ NNIN is the number of inputs required by the type.
!------ NNOUT is the number of outputs required by the type.
!------ NNPAR is the number of parameters for the type.
!------ IODENT is an array variable used to store the number
!       designating the type of output connection, i.e., Pressure,
!       Temperature, Flow rate, etc.
!------ IIDENT is the array used to store the number designating the
!       type of output connection.
!------ OMSSGE an array variable that contains the description of
!       each of the output connections for the type.
!------ IMSSGE an array that contains the description of the input
!       connections for the type.
!------ PCHECK is the variable used to flag the beginning of the
!       parameter information for the type.
!------ PARNUM is assigned to a parameter number for the type.
!------ PMSSGE contains the description of the parameter as selected
!       by the variable PARNUM and the UNIT number.
!------ PVALUE is doubly subscripted variable that contains the value
!       of the parameter as designated by PARNUM.
!------ INDXIN is the doubly subscripted array variable that contains
!       the input indices for the unit as specified by the UNIT number.
!------ INDXOT is the doubly subscripted array that contains the output
!       indices for the unit.
!------ STATUS is a variable that indicates the type of error
!       that has ocurred, if any in the subroutine: 1 is an
!       end of file error, 2 is a read error.
!------ INTYPE is a doubley subscripted array that contains the type
!       if input - .i.e. pressure is 1, temperature is 2 , etc. This
!       corresponds to the indice assigned in INDXIN.
!------ OTTYPE is a doubley subscripted array that contains the type
!       of output - .i.e. pressure is 1, temperature is 2, etc. This
!       corresponds to the indice assigned in INDXOT.
!------ UNTITL is an array that contains the titles for the sim-
!       ulations by unit number.
!======================================================================

SUBROUTINE crunit(crtflg,edflg,rtn_crunit,STATUS)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=15) :: oprmpt,iprmpt
CHARACTER (LEN=80) :: lntype
character(len=1),dimension(72)           :: double_dash='='

!------ Variable initialization

STATUS=0
rtn_recrt = 0                                                
!------ Flags used in the EDit module are set which indicate where
!       which CReate module made the call.

IF (edflg == 1) GO TO 100
crtflg=1
IF (blkflg == 1) crtflg=2
IF (supflg == 1) crtflg=3
IF (simflg == 1) crtflg=4

!------ Unit number is assigned

UNIT=1
IF (blkflg == 1) UNIT=icnt
IF (supflg == 1) UNIT=supunt
IF (simflg == 1) UNIT=simunt

!------ The type description or type number is entered

100     messag=' Enter the type number (or TYPES for list of types)'
CALL datain(item,anumbr,REAL(maxtyp),1.,icheck,lntype,report, 2,3,0,messag)
CALL recrt(item,crtflg,edflg,rtn_recrt)
IF(rtn_recrt == 2) THEN
  rtn_crunit = 2
  return
elseif(rtn_recrt == 1) then
  go to 100
END IF
IF (report /= 3) icheck = 0

!------ Information is read in from the types data file TYPAR.DAT

CALL typein(4)
IF (istat == 2) THEN
  CALL rite(' That type was not found in the types data file!')
  CALL holdit
  STATUS=1
  RETURN
END IF
IF (istat /= 0) GO TO 112

!------ If the call originated from edit (insert unit) then the
!       necessary assignments are made.

IF (edflg == 1) THEN
  untyp(UNIT)=icheck
  blkin(UNIT)=nnin
  blkout(UNIT)=nnout
  blkpar(UNIT)=nnpar
END IF

!------ Types title is echoed to the terminal

WRITE(*,500) dscrib
500     FORMAT(/1X,a72,/)

!------ A loop is entered to define the input indices.

WRITE(*,502)
502     FORMAT(/,' INPUTS',/)
DO  j=1,nnin
  WRITE(*,fmt='(/72a1)') double_dash             
  WRITE(*,12) imssge(j)
  12      FORMAT(1X,'index for ',a60)
  
!------ The input type is assigned to the array INTYPE, corresponds with
!       the indices as defined by INDXIN(UNIT,J)
  
  intype(UNIT,j)=iident(j)
  
!------ A call is made to PROMPTD to obtain the label for the index
  
  CALL promptd(iident,j,iprmpt)
  130     messag=iprmpt
  CALL datain(item,anumbr,REAL(maxunt*minoiu),0.,indxin(UNIT,j),  &
      line,report,2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    rtn_crunit = 2
    RETURN
  elseif(rtn_recrt == 1) then
    go to 130
  END IF
END DO

!------ The loop to define the output indices is entered.

WRITE(*,501)
501     FORMAT(/,' OUTPUTS',/)
DO  i=1,nnout
  WRITE(*,fmt='(/72a1)') double_dash             
  WRITE(*,11) omssge(i)
  11      FORMAT(1X,'index for ',a60)
  
!------ The output type is assigned to the array OTTYPE, this
!       corresponds with the indice assigned by ONDXOT(UNIT,J)
  
  ottype(UNIT,i)=iodent(i)
  
!------ A call is made to PROMPTD to obtain the label for the index
  
  CALL promptd(iodent,i,oprmpt)
  120     messag=oprmpt
  CALL datain(item,anumbr,REAL(maxunt*minoiu),0.,indxot(UNIT,i),  &
      line,report,2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    rtn_crunit = 2
    RETURN
  elseif(rtn_recrt == 1) then
    go to 120
  END IF
END DO

!------ A loop is entered to define the parameters

DO  k=1,nnpar
  WRITE(*,fmt='(/72a1)') double_dash             
  WRITE(*,14) parnum(k),pmssge(k)
  14        FORMAT(1X,'PARAMETER ',i2,1X,a60)
  140       messag=' ENTER THE VALUE'
  CALL datain(item,pvalue(UNIT,k),1.e35,-1.e35,intnum,line,  &
      report,2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    rtn_crunit = 2
    RETURN
  elseif(rtn_recrt == 1) then
    go to 140
  END IF
END DO

!------ Error flags for the read statements
350     IF (blkflg == 1) RETURN
IF (edflg == 1) THEN
  STATUS=0
  RETURN
END IF
CALL fsave(crtflg,edflg,rtn_fsave)
IF(rtn_fsave == 1) then
  rtn_crunit = 1
  return
elseif(rtn_fsave == 2) THEN
  rtn_crunit = 2
  RETURN
END IF
112     CALL rite(' Read error on TYPES data file')
STATUS=2
rtn_crunit = 1  
300     RETURN
END SUBROUTINE crunit
!=======================================================================

!     SUBROUTINE CRBLK

!------ This subroutine creates a BLock
!------ BLKFLG is used as a flag for the create unit subroutine
!       indicating that it was called by the create block subroutine.
!------ ICNT is used as a counter for the number of units.
!------ UNTNUM is an array that stores the unit numbers after each
!       call to the create unit subroutine. Ths subscripts are
!       dependent on the original calling module - i.e. the superblock
!       or block module.
!------ UNTYP is an array that stores the type number of each unit
!       after the call to the create unit sub.
!------ BLOCK is the variable that sets the default value of the
!       number of blocks to 1. This is reassigned if this module is
!       called by the create SUperblock module.
!------ BLKIN is an array that contains the number of input.
!------ BLKOUT is an array that defines the number of outputs.
!------ BLKPAR is an array that defines the number of parameters.
!=======================================================================

SUBROUTINE crblk(crtflg,edflg,rtn_crblk)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=1) :: answer

!       Flags are set for the EDit module.

crtflg=2
IF (supflg == 1) crtflg=3
IF (simflg == 1) crtflg=4

!------ Block number is assigned

BLOCK=1

!------ Block count is handled by the create SUperblock module if
!       the call was origninally from that module.

IF (supflg == 1) BLOCK=blkcnt

!------ Block count is handled by the create SImulation module it
!       the call originated there.

IF (simflg == 1) BLOCK=simblk

WRITE(*,3) BLOCK
3       FORMAT(' BLock number= ',i2)

!------ Flag for the create UNit module is set so there will not be
!       an automatic save from that subroutine.

blkflg=1

!------ Initialize a variable to count the number of units

icnt=1
120     CALL crunit(crtflg,edflg,rtn_crunit,STATUS)
if(rtn_crunit == 2) THEN
  rtn_crblk = 2
  RETURN
END IF
125     IF (simflg == 1) GO TO 127
IF (supflg == 1) GO TO 126
untnum(icnt)=UNIT
untyp(icnt)=icheck
blkin(untnum(icnt))=nnin
blkout(untnum(icnt))=nnout
blkpar(untnum(icnt))=nnpar
GO TO 128

!------ If original call is from the SUperblock create module
!       the unit count is different than the block unit count.

126     untnum(supunt)=UNIT
untyp(supunt)=icheck
blkin(untnum(supunt))=nnin
blkout(untnum(supunt))=nnout
blkpar(untnum(supunt))=nnpar
blkunt(blkcnt)=icnt
GO TO 128

!------ If original call is from the create SImulation module the
!       unit count is different than the superblock unit count.

127     untnum(simunt)=UNIT
untyp(simunt)=icheck
blkin(untnum(simunt))=nnin
blkout(untnum(simunt))=nnout
blkpar(untnum(simunt))=nnpar
blkunt(simblk)=icnt

!------ Count variable is incremented

128     simunt=simunt+1
supunt=supunt+1
icnt=icnt+1
IF (STATUS == 1)icnt=icnt-1
IF (STATUS == 1) simunt=simunt-1
IF (STATUS == 1) supunt=supunt-1
j=muntib+1-icnt
WRITE(*,1)muntib
1       FORMAT(1X,'Maximum number of UNits per BLock= ',i2)
WRITE(*,2) j
2       FORMAT(1X,'Number of available UNits= ',i2)
IF (j <= 0) GO TO 131

!------ User is prompted for the addition of more units

130     messag=' Do you wish to continue entering UNits, Y/N'
CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
CALL recrt(item,crtflg,edflg,rtn_recrt)
if(rtn_recrt == 1) then
  go to 130
elseif(rtn_recrt == 2) THEN
  rtn_crblk = 2
  RETURN
END IF
answer=item
IF (answer == 'Y' .OR. answer == 'y') GO TO 120

!------ Was BLock call by the create superblock module?

131     IF (supflg == 1) then
           rtn_crblk = 1
           return
        endif   
CALL fsave(crtflg,edflg,rtn_fsave)
if(rtn_fsave == 1) then
  rtn_crblk = 1
elseif(rtn_fsave == 2) then
  rtn_crblk = 2  
endif
350     RETURN
END SUBROUTINE crblk
!======================================================================

!     SUBROUTINE CRSUP

!------ This subroutine creates a superblock
!------ SUPFLG is a variable indicating a call by the superblock module
!------ BLKCNT is a variable used to count the number of blocks per
!       superblock.
!------ BLKNUM is an array that stores the block numbers.
!------ SUPUNT is a variable that counts the number of units used at
!       the superblock level. This count is continuous at this level.
!------ SUPCNT counts the number of superblocks and sets the default
!       to 1 if a call is not made by the create SImulation module.
!------ UNTCNT is unused and will be deleted.
!------ BLKUNT is an array that stores the number of units for each
!       block.
!=======================================================================

SUBROUTINE crsup(crtflg,edflg,rtn_crsup)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=1) :: answer
integer           :: maxbk 

!------ Superblock number is assigned and other variables
!       are initialized.

supcnt=1
supflg=1
blkcnt=1
supunt=1

!------ Flags for the EDit module are set.

crtflg=3
IF (simflg == 1) crtflg=4

IF (simflg == 1) supcnt=simsup
WRITE(*,126) supcnt
126     FORMAT(' Superblock number= ',i2)
120     CALL crblk(crtflg,edflg,rtn_crblk)
IF(rtn_crblk == 2) THEN
  rtn_crsup = 2
  RETURN
END IF

125     blknum(blkcnt)=BLOCK
simcnt(supcnt)=blkcnt
IF (simflg == 1) blknum(simblk)=BLOCK

!------ Count variable is incremented.

simblk=simblk+1
blkcnt=blkcnt+1
maxbk=mblkis+1-blkcnt
WRITE(*,1)mblkis
1       FORMAT(' Maximum number of BLocks per SUperblock= ',i2)
WRITE(*,2)maxbk
2       FORMAT(' Number of available BLocks=',i2)

!------ User is prompted for addition of more BLocks

130     messag=' Do you wish to continue entering BLocks, Y/N'
CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
CALL recrt(item,crtflg,edflg,rtn_recrt)
IF(rtn_recrt == 2) THEN
  rtn_crsup = 2
  RETURN
elseif(rtn_recrt == 1) then
  goto 130
END IF
answer=item
IF (answer == 'Y' .OR. answer == 'y') GO TO 120

!------ Did call originate from the create SImulation module?

IF (simflg == 1) RETURN
CALL fsave(crtflg,edflg,rtn_fsave)
if(rtn_fsave == 2) then
  rtn_crsup = 2
endif

RETURN
END SUBROUTINE crsup
!=======================================================================

!     SUBROUTINE CRSIM

!------ This module creates a simulation.
!------ SIMSUP counts the number of superblocks per simulation.
!------ SIMFLG is a flag indicating call by SImulation module.
!------ SIMBLK counts number of blocks per simulation.
!------ SIMUNT counts number of units per simulation.
!------ SIMCNT array that stores the number of blocks per superblock.
!------ SUPNUM array that stores the number of superblocks.
!------ STATE a vector array that stores the initial values of the
!       state vector.
!------ PRES,FLOW,TEMP,CNTR,OTHR,ENRG,AHUM are all variables that store
!       the maximum variable indice values for Pressure,Flow,Temp, etc.
!       respectively.
!------ BNDRY contains the state vector position of the boundary
!       variables to be used in the simulation.
!------ BVPRES,BVFLOW,BVTEMP,BVCNTR,BVOTHR,BVENRG,BVPOWR,BVAHUM, are
!       used to define state vector position of the boundary variables
!       for each of the following; Pressure,Temperature,Flow, etc.
!------ RPPRES,RPFLOW,RPTEMP,RPCNTR,RPOTHR,RPENRG,RPPOWR,RPAHUM, are
!       used to define the state vector position for the variables
!       that are to be used for reporting during a simulation.
!------ NTYP contains the type description of each variable in the
!       state vector.
!------ NINDX contains the state vector position of the variables
!       corresponding to NTYP and the RP____ variables.
!------ REPINV array containing the reporting interval per superblock.
!------ RTOLX,ATOLX,XTOL,TTIME, are all variables required by MODSIM
!       to run a simulation (See MODSIM for more details)
!------ SCAN contains the scanning option per superblock.
!------ FREEZE contains the freezing option per superblock.
!------ TITLSM contains the tile for the simulation.
!======================================================================

SUBROUTINE crsim(crtflg,edflg,rtn_crsim)

use hvaccomm
implicit none
include 'hvacgen.inc'

CALL crsim1(crtflg,edflg,abort)
IF (abort) then
  rtn_crsim = 2
  RETURN
endif
CALL crsim2(crtflg,edflg,abort)
IF (abort) then
  rtn_crsim = 2
  RETURN
endif
CALL crsim3(crtflg,edflg,abort)
IF (abort) then
  rtn_crsim = 2
  RETURN
endif
CALL crsim4(crtflg,edflg,abort)
if(abort) then
  rtn_crsim = 2
  RETURN
endif
END SUBROUTINE crsim
!=======================================================================

!     SUBROUTINE CRSIM1

!=======================================================================

SUBROUTINE crsim1(crtflg,edflg,abort)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=1) :: answer
integer           :: maxsup

!------ Variables are initialized.

simsup=1
simflg=1
simblk=1
simunt=1

!------ Flags are set for the EDit module.

crtflg=4
120     CALL crsup(crtflg,edflg,rtn_crsup)
IF (rtn_crsup == 2) THEN
  abort = .true.
  RETURN
END IF
125     supnum(supcnt)=simsup

!------ Count variable is incremented.

simsup=simsup+1
maxsup=maxsbk+1-simsup
WRITE(*,1)maxsbk
1       FORMAT(' Maximum number of SUperblocks per SImulation= ',i2)
WRITE(*,2) maxsup
2       FORMAT(' Number of available SUperblocks= ',i2)

!------ User is prompted for the addition of more SUperblocks.

130     messag=' Do you wish to continue entering SUperblocks, Y/N'
CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
CALL recrt(item,crtflg,edflg,rtn_recrt)
IF (rtn_recrt == 2) THEN
  abort = .true.
  RETURN
elseif(rtn_recrt == 1) then
  go to 130
END IF

answer=item
IF (answer == 'Y' .OR. answer == 'y') GO TO 120

!------ Entering the simulation title

messag=' Enter the title for this SImulation'
CALL datain(item,anumbr,0.,0.,intnum,line,report,3,2,1,messag)

!--commentin out the line below allows any title to be used, but
!  does  not  allow reserved words (e.g.ABORT) to have any effect. To
!  replace, put label 499 in front of assignment to variable MESSAG
!  above.

titlsm=line

!------ Entering the error tolerances, RTOLX, ATOLX, XTOL, TTIME for
!       the simulation.

!------ The defaults are:

rtolx=.0001
atolx=.00001
xtol=.0002
ttime=1.0
WRITE(*,10)
10      FORMAT(' The default values for the error tolerances are as foll',  &
    'ows:'/)
WRITE(*,11)rtolx,atolx
11      FORMAT(' RTOLX=',e15.6,2X,'ATOLX=',e15.6)
WRITE(*,12)xtol,ttime
12      FORMAT(' XTOL=',1X,e15.6,2X,'TTIME=',e15.6,/)
599     messag=' Use these default values, Y/N'
CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
CALL recrt(item,crtflg,edflg,rtn_recrt)
IF (rtn_recrt == 2) THEN
  abort = .true.
  RETURN
elseif(rtn_recrt == 1) then
  go to 599
END IF

answer=item
IF (answer == 'Y'.OR.answer == 'y') GO TO 244
598     messag=' Enter RTOLX'
CALL datain(item,rtolx,1.e35,-1.e35,intnum,line,report, 2,3,0,messag)
CALL recrt(item,crtflg,edflg,rtn_recrt)
IF (rtn_recrt == 2) THEN
  abort = .true.
  RETURN
elseif(rtn_recrt == 1) then
  go to 598
END IF
597     messag=' Enter ATOLX'
CALL datain(item,atolx,1.e35,-1.e35,intnum,line,report, 2,3,0,messag)
CALL recrt(item,crtflg,edflg,rtn_recrt)
IF (rtn_recrt == 2) THEN
  abort = .true.
  RETURN
elseif(rtn_recrt == 1) then
  go to 597
END IF
596     messag=' Enter XTOL'
CALL datain(item,xtol,1.e35,-1.e35,intnum,line,report, 2,3,0,messag)
CALL recrt(item,crtflg,edflg,rtn_recrt)
IF (rtn_recrt == 2) THEN
  abort = .true.
  RETURN
elseif(rtn_recrt == 1) then
  go to 596
END IF
595     messag=' Enter TTIME'
CALL datain(item,ttime,1.e35,-1.e35,intnum,line,report, 2,3,0,messag)
CALL recrt(item,crtflg,edflg,rtn_recrt)
IF (rtn_recrt == 2) THEN
  abort = .true.
  RETURN
elseif(rtn_recrt == 1) then
  go to 595
END IF
244     abort = .false.
RETURN
END SUBROUTINE crsim1
!=======================================================================

!     SUBROUTINE CRSIM2

!=======================================================================

SUBROUTINE crsim2(crtflg,edflg,abort)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer,dimension(maxbnd)           :: bvpres,bvflow,bvtemp,bvcntr, &
                                       bvothr,bvenrg,bvpowr,bvahum
integer,dimension(maxsbk,mrptis)    :: rppres,rpflow,rptemp,rpcntr, &
                                       rpothr,rpenrg,rppowr,rpahum

!------ Initilize variables that are used to count the number of
!      occurrences of the various types, ie.,Pressure, Temperature, etc.

pres=0
flow=0
temp=0
cntr=0
othr=0
enrg=0
powr=0
ahum=0

!------ Arrays are zeroed as required

DO  i=1,maxbnd
  bvpres(i)=0
  bvflow(i)=0
  bvtemp(i)=0
  bvcntr(i)=0
  bvothr(i)=0
  bvenrg(i)=0
  bvpowr(i)=0
  bvahum(i)=0
END DO
DO  i=1,maxsbk
  DO  j=1,mrptis
    rppres(i,j)=0
    rpflow(i,j)=0
    rptemp(i,j)=0
    rpcntr(i,j)=0
    rpothr(i,j)=0
    rpenrg(i,j)=0
    rppowr(i,j)=0
    rpahum(i,j)=0
  END DO
END DO
DO  i=1,maxsbk
  DO  j=1,mdmrep
    ntyp(i,j)=0
    nindx(i,j)=0
    repvar(i,j)=0
  END DO
END DO
DO  i=1,mdmbnd
  bndry(i)=0
END DO

DO  ii=1,simunt-1
  DO  jj=1,blkin(ii)
    IF (intype(ii,jj) == 1.AND.indxin(ii,jj) >= pres) pres= indxin(ii,jj)
    IF (intype(ii,jj) == 2.AND.indxin(ii,jj) >= flow) flow= indxin(ii,jj)
    IF (intype(ii,jj) == 3.AND.indxin(ii,jj) >= temp) temp= indxin(ii,jj)
    IF (intype(ii,jj) == 4.AND.indxin(ii,jj) >= cntr) cntr= indxin(ii,jj)
    IF (intype(ii,jj) == 5.AND.indxin(ii,jj) >= othr) othr= indxin(ii,jj)
    IF (intype(ii,jj) == 6.AND.indxin(ii,jj) >= enrg) enrg= indxin(ii,jj)
    IF (intype(ii,jj) == 7.AND.indxin(ii,jj) >= powr) powr= indxin(ii,jj)
    IF (intype(ii,jj) == 8.AND.indxin(ii,jj) >= ahum) ahum= indxin(ii,jj)
  END DO
  DO  jj=1,blkout(ii)
    IF (ottype(ii,jj) == 1.AND.indxot(ii,jj) >= pres) pres= indxot(ii,jj)
    IF (ottype(ii,jj) == 2.AND.indxot(ii,jj) >= flow) flow= indxot(ii,jj)
    IF (ottype(ii,jj) == 3.AND.indxot(ii,jj) >= temp) temp= indxot(ii,jj)
    IF (ottype(ii,jj) == 4.AND.indxot(ii,jj) >= cntr) cntr= indxot(ii,jj)
    IF (ottype(ii,jj) == 5.AND.indxot(ii,jj) >= othr) othr= indxot(ii,jj)
    IF (ottype(ii,jj) == 6.AND.indxot(ii,jj) >= enrg) enrg= indxot(ii,jj)
    IF (ottype(ii,jj) == 7.AND.indxot(ii,jj) >= powr) powr= indxot(ii,jj)
    IF (ottype(ii,jj) == 8.AND.indxot(ii,jj) >= ahum) ahum= indxot(ii,jj)
  END DO
END DO
WRITE(*,1001) pres,flow,temp,cntr,othr,enrg,powr,ahum
1001    FORMAT(' PRES=',i3,' FLOW=',i3,' TEMP=',i3,' CNTR=',i3,' OTHR=',  &
    i3,' ENRG=',i3,' POWR=',i3,' AHUM=',i3)
WRITE(*,89)
89      FORMAT(///,' Entering Variable Initial values:',/)

IF (pres <= 0) GO TO 550
DO  i=1,pres
  WRITE(*,3)i,iounit(1)
  3         FORMAT(/,' ENTER PRESSURE',i3,a20)
  600       messag=' '
  CALL datain(item,state(i),1.e35,-1.e35,intnum,line,report, 2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 600
  END IF
END DO

550     IF (flow <= 0) GO TO 551
DO  j=1,flow
  WRITE(*,4)j,iounit(2)
  4         FORMAT(/,' ENTER FLOWRATE',i3,a20)
  601       messag=' '
  CALL datain(item,state(j+pres),1.e35,-1.e35,intnum,line,  &
      report,2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 601
  END IF
END DO

551     IF (temp <= 0) GO TO 552
DO  l=1,temp
  WRITE(*,5)l,iounit(3)
  5         FORMAT(/,' ENTER TEMPERATURE',i3,a20)
  602       messag=' '
  CALL datain(item,state(l+pres+flow),1.e35,-1.e35,intnum,line,  &
      report,2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 602
  END IF
END DO

552     IF (cntr <= 0) GO TO 553
DO  m=1,cntr
  WRITE(*,6)m,iounit(4)
  6         FORMAT(/,' ENTER CONTROL',i3,a20)
  603       messag=' '
  CALL datain(item,state(m+pres+temp+flow),1.e35,-1.e35,intnum,  &
      line,report,2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 603
  END IF
END DO

553     IF (othr <= 0) GO TO 554
DO  n=1,othr
  WRITE(*,7)n,iounit(5)
  7         FORMAT(/,' ENTER OTHER',i3,a20)
  604       messag=' '
  CALL datain(item,state(n+pres+flow+temp+cntr),1.e35,-1.e35,  &
      intnum,line,report,2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 604
  END IF
END DO

554     IF (enrg <= 0) GO TO 555
DO  ii=1,enrg
  WRITE(*,8)ii,iounit(6)
  8         FORMAT(/,' ENTER ENERGY',i3,a20)
  605       messag=' '
  CALL datain(item,state(ii+pres+temp+flow+cntr+othr),1.e35,  &
      -1.e35,intnum,line,report,2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 605
  END IF
END DO

555     IF (powr <= 0) GO TO 558
DO  jj=1,powr
  WRITE(*,9)jj,iounit(7)
  9         FORMAT(/,' ENTER POWER',i3,a20)
  606       messag=' '
  CALL datain(item,state(jj+pres+temp+flow+cntr+enrg+othr),  &
      1.e35, -1.e35,intnum,line,report,2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 606
  END IF
END DO

558     IF (ahum <= 0) GO TO 557
DO  jj=1,ahum
  WRITE(*,19)jj,iounit(8)
  19        FORMAT(/,' ENTER HUMIDITY',i3,a20)
  607       messag=' '
  CALL datain(item,state(jj+pres+temp+flow+cntr+othr+enrg+powr)  &
      ,1.e35,-1.e35,intnum,line,report,2,3,0,messag)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 607
  END IF
END DO

557     abort = .false.
RETURN
END SUBROUTINE crsim2
!=======================================================================

!     SUBROUTINE CRSIM3

!=======================================================================

SUBROUTINE crsim3(crtflg,edflg,abort)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer,dimension(maxbnd)           :: bvpres,bvflow,bvtemp,bvcntr, &
                                       bvothr,bvenrg,bvpowr,bvahum

!------ Boundary variables are entered

WRITE(*,88)
88      FORMAT(///,' Entering Boundary Variables:',/)

DO  i=1,pres
  750     messag=' Enter a PRESSURE boundary variable or CR to move on'
  CALL datain(item,anumbr,REAL(pres),0.,bvpres(i),line,report, 2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 750
  END IF
  bndry(i)=bvpres(i)
  IF (bvpres(i) <= 0) bndry(i)=0
  IF (bvpres(i) <= 0) GO TO 710
END DO

710     DO  i=1,flow
  751       messag=' Enter a FLOW boundary variable or CR to move on'
  CALL datain(item,anumbr,REAL(flow),0.,bvflow(i),line,report, 2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 751
  END IF
  bndry(i+pres)=bvflow(i)+pres
  IF (bvflow(i) <= 0) bndry(i+pres)=0
  IF (bvflow(i) <= 0) GO TO 711
END DO

711     DO  i=1,temp
  752       messag=' Enter a TEMPERATURE boundary variable or CR to move on'
  CALL datain(item,anumbr,REAL(temp),0.,bvtemp(i),line,report, 2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 752
  END IF
  bndry(i+pres+flow)=bvtemp(i)+pres+flow
  IF (bvtemp(i) <= 0) bndry(i+pres+flow)=0
  IF (bvtemp(i) <= 0) GO TO 712
END DO

712     DO  i=1,cntr
  753       messag=' Enter a CONTROL boundary variable or CR to move on'
  CALL datain(item,anumbr,REAL(cntr),0.,bvcntr(i),line,report, 2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 753
  END IF
  bndry(i+pres+flow+temp)=bvcntr(i)+pres+flow+temp
  IF (bvcntr(i) <= 0) bndry(i+pres+flow+temp)=0
  IF (bvcntr(i) <= 0) GO TO 713
END DO

713     DO  i=1,othr
  754       messag=' Enter an OTHER boundary variable or CR to move on'
  CALL datain(item,anumbr,REAL(othr),0.,bvothr(i),line,report, 2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 754
  END IF
  bndry(i+pres+flow+temp+cntr)=bvothr(i)+pres+flow+temp+cntr
  IF (bvothr(i) <= 0) bndry(i+pres+flow+temp+cntr)=0
  IF (bvothr(i) <= 0) GO TO 714
END DO

714     DO  i=1,enrg
  755       messag=' Enter an ENERGY boundary variable or CR to move on'
  CALL datain(item,anumbr,REAL(enrg),0.,bvenrg(i),line,report, 2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 755
  END IF
  bndry(i+pres+flow+temp+cntr+othr)=bvenrg(i)+pres+flow+temp +cntr+othr
  IF (bvenrg(i) <= 0) bndry(i+pres+flow+temp+cntr+othr)=0
  IF (bvenrg(i) <= 0) GO TO 715
END DO

715     DO  i=1,powr
  756       messag=' Enter a POWER boundary variable or CR to move on'
  CALL datain(item,anumbr,REAL(powr),0.,bvpowr(i),line,report, 2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 756
  END IF
  bndry(i+pres+flow+temp+cntr+othr+enrg)=bvpowr(i)+pres+flow  &
      +temp+cntr+othr+enrg
  IF (bvpowr(i) <= 0) bndry(i+pres+flow+temp+cntr+othr+enrg)=0
  IF (bvpowr(i) <= 0) GO TO 716
END DO

716     DO  i=1,ahum
  757       messag=' Enter a HUMIDITY boundary variable or CR to move on'
  CALL datain(item,anumbr,REAL(ahum),0.,bvahum(i),line,report, 2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 757
  END IF
  bndry(i+pres+flow+temp+cntr+othr+enrg+powr)=bvahum(i)+pres  &
      +flow+temp+cntr+othr+enrg+powr
  IF (bvahum(i) <= 0) bndry(i+pres+flow+temp+cntr+othr+enrg+powr)=0
  IF (bvahum(i) <= 0) GO TO 856
END DO

856     abort = .false.
RETURN
END SUBROUTINE crsim3
!=======================================================================

!     SUBROUTINE CRSIM4

!=======================================================================

SUBROUTINE crsim4(crtflg,edflg,abort)

use hvacsim_par
use hvaccomm
include 'hvacgen.inc'
integer,dimension(maxsbk,mrptis)    :: rppres,rpflow,rptemp,rpcntr, &
                                       rpothr,rpenrg,rppowr,rpahum


!------ Reported variables are entered

!------ Initialize count variables

DO  j=1,simsup-1
  WRITE(*,992)j
  992     FORMAT(///,' Entering Reported Variables for Superblock ',i2,//)
  957     messag=' Enter the reporting interval for this SUperblock in seconds'
  CALL datain(item,repinv(j),1.e35,0.,intnum,line,report, 2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 957
  END IF

  DO  i=1,pres
    950       messag=' Enter a PRESSURE reported variable or CR to move on'
    CALL datain(item,anumbr,REAL(pres),0.,rppres(j,i),line,report,  &
        2,3,0,messag)
    CALL recrt(item,crtflg,edflg,rtn_recrt)
    IF(rtn_recrt == 2) THEN
      abort = .true.
      RETURN
    elseif(rtn_recrt == 1) then
      go to 950
    END IF
    repvar(j,i)=rppres(j,i)
    IF (rppres(j,i) <= 0) repvar(j,i)=0
    IF (rppres(j,i) <= 0) GO TO 910
    ntyp(j,i)=1
    nindx(j,i)=rppres(j,i)
  END DO
  
  910     DO  i=1,flow
    951       messag=' Enter a FLOW reported variable or CR to move on'
    CALL datain(item,anumbr,REAL(flow),0.,rpflow(j,i),line,report,  &
        2,3,0,messag)
    CALL recrt(item,crtflg,edflg,rtn_recrt)
    IF(rtn_recrt == 2) THEN                        
      abort = .true.                               
      RETURN                                       
    elseif(rtn_recrt == 1) then                    
      go to 951                                    
    END IF                                         
    repvar(j,i+pres)=rpflow(j,i)+pres
    IF (rpflow(j,i) <= 0) repvar(j,i+pres)=0
    IF (rpflow(j,i) <= 0) GO TO 911
    ntyp(j,i+pres)=2
    nindx(j,i+pres)=rpflow(j,i)
  END DO
  
  911     DO  i=1,temp
    952       messag=' Enter a TEMPERATURE reported variable or CR to move on'
        
    CALL datain(item,anumbr,REAL(temp),0.,rptemp(j,i),line,report,  &
        2,3,0,messag)
    CALL recrt(item,crtflg,edflg,rtn_recrt)
    IF(rtn_recrt == 2) THEN                        
      abort = .true.                               
      RETURN                                       
    elseif(rtn_recrt == 1) then                    
      go to 952
    END IF                                         
    repvar(j,i+pres+flow)=rptemp(j,i)+pres+flow
    IF (rptemp(j,i) <= 0) repvar(j,i+pres+flow)=0
    IF (rptemp(j,i) <= 0) GO TO 912
    nindx(j,i+pres+flow)=rptemp(j,i)
    ntyp(j,i+pres+flow)=3
  END DO
  
  912     DO  i=1,cntr
    953       messag=' Enter a CONTROL reported variable or CR to move on'
    CALL datain(item,anumbr,REAL(cntr),0.,rpcntr(j,i),line,report,  &
        2,3,0,messag)
    CALL recrt(item,crtflg,edflg,rtn_recrt)
    IF(rtn_recrt == 2) THEN                        
      abort = .true.                               
      RETURN                                       
    elseif(rtn_recrt == 1) then                    
      go to 953
    END IF                                         
    repvar(j,i+pres+flow+temp)=rpcntr(j,i)+pres+flow+temp
    IF (rpcntr(j,i) <= 0) repvar(j,i+pres+flow+temp)=0
    IF (rpcntr(j,i) <= 0) GO TO 913
    ntyp(j,i+pres+flow+temp)=4
    nindx(j,i+pres+flow+temp)=rpcntr(j,i)
  END DO
  
  913     DO  i=1,othr
    954       messag=' Enter an OTHER reported variable or CR to move on'
    CALL datain(item,anumbr,REAL(othr),0.,rpothr(j,i),line,report,  &
        2,3,0,messag)
    CALL recrt(item,crtflg,edflg,rtn_recrt)
    IF(rtn_recrt == 2) THEN                        
      abort = .true.                               
      RETURN                                       
    elseif(rtn_recrt == 1) then                    
      go to 954
    END IF                                         
    repvar(j,i+pres+flow+temp+cntr)=rpothr(j,i)+pres+flow+temp +cntr
    IF (rpothr(j,i) <= 0) repvar(j,i+pres+flow+temp+cntr)=0
    IF (rpothr(j,i) <= 0) GO TO 914
    ntyp(j,i+pres+flow+temp+cntr)=5
    nindx(j,i+pres+flow+temp+cntr)=rpothr(j,i)
  END DO
  
  914     DO  i=1,enrg
    955       messag=' Enter an ENERGY reported variable or CR to move on'
    CALL datain(item,anumbr,REAL(enrg),0.,rpenrg(j,i),line,report,  &
        2,3,0,messag)
    CALL recrt(item,crtflg,edflg,rtn_recrt)
    IF(rtn_recrt == 2) THEN                        
      abort = .true.                               
      RETURN                                       
    elseif(rtn_recrt == 1) then                    
      go to 955
    END IF                                         
    repvar(j,i+pres+flow+temp+cntr+othr)=rpenrg(j,i)+pres+flow+ temp+cntr+othr
    IF (rpenrg(j,i) <= 0) repvar(j,i+pres+flow+temp+cntr+othr)=0
    IF (rpenrg(j,i) <= 0) GO TO 915
    ntyp(j,i+pres+flow+temp+cntr+othr)=6
    nindx(j,i+pres+flow+temp+cntr+othr)=rpenrg(j,i)
  END DO
  
  915     DO  i=1,powr
    956       messag=' Enter a POWER reported variable or CR to move on'
    CALL datain(item,anumbr,REAL(powr),0.,rppowr(j,i),line,report,  &
        2,3,0,messag)
    CALL recrt(item,crtflg,edflg,rtn_recrt)
    IF(rtn_recrt == 2) THEN                        
      abort = .true.                               
      RETURN                                       
    elseif(rtn_recrt == 1) then                    
      go to 956
    END IF                                         
    repvar(j,i+pres+flow+temp+cntr+othr+enrg)=rppowr(j,i)+pres+  &
        flow+temp+cntr+othr+enrg
    IF (rppowr(j,i) <= 0) repvar(j,i+pres+flow+temp+cntr+othr+enrg)=0
    IF (rppowr(j,i) <= 0) GO TO 916
    ntyp(j,i+pres+flow+temp+cntr+othr+enrg)=7
    nindx(j,i+pres+flow+temp+cntr+othr+enrg)=rppowr(j,i)
  END DO
  
  916     DO  i=1,ahum
    960       messag=' Enter a HUMIDITY reported variable or CR to move on'
    CALL datain(item,anumbr,REAL(ahum),0.,rpahum(j,i),line,report,  &
        2,3,0,messag)
    CALL recrt(item,crtflg,edflg,rtn_recrt)
    IF(rtn_recrt == 2) THEN                        
      abort = .true.                               
      RETURN                                       
    elseif(rtn_recrt == 1) then                    
      go to 960
    END IF                                         
    repvar(j,i+pres+flow+temp+cntr+othr+enrg+powr)=rpahum(j,i)  &
        +pres+flow+temp+cntr+othr+enrg+powr
    IF (rpahum(j,i) <= 0) repvar(j,i+pres+flow+temp+cntr+othr+enrg +powr)=0
    IF (rpahum(j,i) <= 0) GO TO 958
    ntyp(j,i+pres+flow+temp+cntr+othr+enrg+powr)=8
    nindx(j,i+pres+flow+temp+cntr+othr+enrg+powr)=rpahum(j,i)
  END DO
  
  958     messag=' Enter the variable freezing option for this SUperblock: &
     & 0,1,or 2.'
  CALL datain(item,anumbr,2.,0.,freeze(j),line,report, 2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 958
  END IF
  959     messag=' Enter the variable input scan option for this SUperblock: &
      & 0 or 1.'
  CALL datain(item,anumbr,1.,0.,scan(j),line,report,2,3,0,messag)
  CALL recrt(item,crtflg,edflg,rtn_recrt)
  IF(rtn_recrt == 2) THEN
    abort = .true.
    RETURN
  elseif(rtn_recrt == 1) then
    go to 959
  END IF
END DO
CALL fsave(crtflg,edflg,rtn_fsave)
350     RETURN
END SUBROUTINE crsim4
!=======================================================================

!       SUBROUTINE RECRT

!------ This subroutine is called by the CREATE module and handles
!       all calls to other selected modules. This subroutine will
!       only allow calls to the FILES and HELP modules from the
!       CREATE module.
!------ POINTR is described in REWORD.

!=======================================================================

SUBROUTINE recrt(item,crtflg,edflg,rtn_recrt)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=1) :: canned
CHARACTER (LEN=2) :: answer

answer=item
IF (answer == 'AB'.OR.answer == 'ab') THEN
  messag=' Are you sure you wish to ABort? NOTHING will be saved . Yes or No?'
  CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
  canned = item
  IF (canned == 'Y'.OR.canned == 'y') then
    rtn_recrt = 2
  else
    rtn_recrt = 1
  endif
  
!------ A call is made to HELP with program control transferred there.
  
ELSE IF (answer == 'HE'.OR.answer == 'he') THEN
  CALL help(crtflg,edflg,abort)
  IF(abort) then
    rtn_recrt = 2
  else
    rtn_recrt = 1
  endif
  
!------ A call is made to the view module with control shifted there.
  
ELSE IF (answer == 'VI'.OR.answer == 'vi') THEN
  CALL view(crtflg,edflg,abort)
  IF(abort) then
    rtn_recrt = 2
  else
    rtn_recrt = 1
  endif
  
!------ A call is made to the edit module with control shifted there.
  
ELSE IF (answer == 'ED'.OR.answer == 'ed') THEN
  CALL edit(crtflg,edflg,abort)
  IF(abort) then
    rtn_recrt = 2
  else
    rtn_recrt = 1
  endif

!-------get contents of TYPAR.DAT file
  
ELSE IF (answer == 'TY'.OR.answer == 'ty') THEN
  CALL types
  IF(abort) then
    rtn_recrt = 2
  else
    rtn_recrt = 1
  endif
END IF
RETURN
END SUBROUTINE recrt
!=======================================================================

!     SUBROUTINE TYPES

!   This routine is used to produce a list of the types available in the
!   types file on the console.

!=======================================================================

SUBROUTINE types

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer            :: ifirst,ilast

ifirst = 1
ilast = 0
lincnt = 1
DO  i = 1,maxtyp
  icheck=i
  CALL typein(1)
  IF (istat == 0) THEN
    IF (ilast > ifirst) THEN
      WRITE(*,1)ifirst,ilast
      1  FORMAT(1X,'TYPES',i3,' THROUGH',i3,' NOT FOUND')
      CALL scroll(lincnt)
    END IF
    IF (ilast == ifirst) THEN
      WRITE(*,2)ilast
      2  FORMAT(1X,'TYPE',i3,' NOT FOUND')
      CALL scroll(lincnt)
    END IF
    ifirst = i + 1
    WRITE(*,3)i,dscrib
      3  FORMAT(1X,'TYPE',i3,' = ',a65)
    CALL scroll(lincnt)
  ELSE
    ilast = i
  END IF
END DO
RETURN
END SUBROUTINE types
! ===========================================================================

