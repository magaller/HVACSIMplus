!=======================================================================
! HVACGEN SOURCE FILE #3 OF 6, VERSION 5.0
!=======================================================================

! hvacgen3.f90

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

!       FILES MODULE

!=======================================================================

!     SUBROUTINE fsave

!------ This subroutine is called by FILES module and selects the
!       type of SAve that should be performed based on the file
!       extension that is passed by the OPNFIL subroutine (See OPENFIL
!       for more details).
!------ EXTNSN is the variable assigned to the file extension. The
!       extension SIM defines a SImulation, BLK defines BLock, and
!       UNT defines UNit.
!=======================================================================

SUBROUTINE fsave(crtflg,edflg,rtn_fsave)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

irdflg=0
250     CALL opnfil(crtflg,edflg,abort)
IF (abort) then
  rtn_fsave = 2
endif  
CALL rite(' Saving to work file....')
IF (extnsn == 'sim') THEN
  CALL savsim(crtflg,edflg,rtn_savsim)
  if(rtn_savsim == 2) then
    CALL rite(' Aborting save...')
    CALL holdit
    rtn_fsave = 2
    RETURN
  elseif(rtn_savsim == 1) then
    rtn_fsave = 1
    return
  endif
ELSE IF (extnsn == 'blk') THEN
  CALL savblk(1,crtflg,edflg,rtn_savblk)
  if(rtn_savblk == 2) then
    CALL rite(' Aborting save...')
    CALL holdit
    rtn_fsave = 2
    RETURN
  elseif(rtn_savblk == 1) then
    rtn_fsave = 1
    return
  endif
ELSE IF (extnsn == 'unt') THEN
  CALL savunt(crtflg,edflg,rtn_savunt)
  if(rtn_savunt == 2) then
    CALL rite(' Aborting save...')
    CALL holdit
    rtn_fsave = 2
    RETURN
  elseif(rtn_savunt == 1) then
    rtn_fsave = 1
    return
  endif
ELSE
  CALL rite(' Illegal file extension!')
  CALL holdit
  GO TO 250
END IF

RETURN
END SUBROUTINE fsave
!=======================================================================

!     SUBROUTINE READIN

!------ This subroutine is called by FILES module and selects the
!       type of REad that should be performed based on the file
!       extension that is passed by the OPNFIL subroutine (See OPENFIL
!       for more details).
!------ EXTNSN is the variable assigned to the file extension. The
!       extension SIM defines a SImulation, SUP defines SUperblock,
!       BLK defines BLock, UNT defines UNit.
!------ IRDFLG a flag that tells the OPNFIL module that the call
!       was from the read module.
!=======================================================================

SUBROUTINE readin(crtflg,edflg,rtn_readin)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

irdflg=1

250     CALL opnfil(crtflg,edflg,abort)
IF(abort) THEN
  rtn_readin = 2
  RETURN
END IF
CALL rite(' Reading from work file....')
IF (extnsn == 'sim') THEN
  CALL rdsim(crtflg,edflg,rtn_rdsim)
  IF(rtn_rdsim == 1) THEN
    rtn_readin = 1
    return
  elseif(rtn_rdsim == 2) then
    rtn_readin = 2
    return
  elseif(rtn_rdsim == 3) then
    go to 250
  END IF
ELSE IF (extnsn == 'blk') THEN
  CALL rdblk(crtflg,edflg,rtn_rdblk)
  IF(rtn_rdblk == 1) THEN
    rtn_readin = 1
    return
  elseif(rtn_rdblk == 2) then
    rtn_readin = 2
    return
  elseif(rtn_rdblk == 3) then
    go to 250
  END IF
ELSE IF (extnsn == 'unt') THEN
  CALL rdunt(crtflg,edflg,rtn_rdunt)
  IF(rtn_rdunt == 1) THEN
    rtn_readin = 1
    return
  elseif(rtn_rdunt == 2) then
    rtn_readin = 2
    return
  elseif(rtn_rdunt == 3) then
    go to 250
  END IF
END IF
CALL rite(' Illegal file extension!')
CALL  holdit
GO TO 250
RETURN
END SUBROUTINE readin
!=======================================================================

!      SUBROUTINE OPNFIL

!----- This subroutine is used to enter the filename and then open
!      the file.
!----- FULLNM is the filename entered by the user.
!----- CHARCT is used to store the filename obtained by stepping
!      through the variable FULLNM one character at a time until
!      a space is detected or 8 characters are reached. This defines
!      the variable name w/o extension. the extension is then tagged
!      on in correspondence to the module originally calling this
!      module.
!----- EXTNSN is the file extension which is limited to a maximum
!      of three.
!----- IRDFLG is an input argument with the value of 1 of a file is
!      to be opened for read and 0 if a file is to be opened for write.
!=======================================================================

SUBROUTINE opnfil(crtflg,edflg,abort)

! The names of files used under different operating systems are
! unfortunately not standard. This routine uses a scheme to encode the
! type of work file into the file name. The assumption is made that
! any computer file name can be characterized as having a name section,
! an extension section, and an extension delimiter. The total file name
! is assummed to be 50 characters or less. Some systems may not have
! extensions and the delimiter may not be present. However, for this
! program the extension is used to encode the work file type. Therefore
! in systems without a formal extension, part of the filename is used
! to code the work file type, and the name section is assummed to be
! shorter. The parameters below are used to customize this routine for
! various systems. LNAME is the length of the name which must be no
! longer than 49.
! LEXT is the length of the extension which should be no smaller than 1
! or larger than 3. The delimiter DELIM, is usually '.'. If DELIM
! is set to ' ', a blank, then no delimiter will be used. If no
! delimiter is used, LFILE + LEXT = the total filename length
! which should not exceed 50. IF a delimiter is used,
! the total filename lenght is = LFILE + LEXT + 1.

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

SAVE

INTEGER, PARAMETER           :: lname=46
INTEGER, PARAMETER           :: lext=3
CHARACTER (LEN=1), PARAMETER :: delim='.'

LOGICAL                      :: mistak
CHARACTER (LEN=1)            :: canned,charct
CHARACTER (LEN=2)            :: answer
CHARACTER (LEN=46)           :: curfil
CHARACTER (LEN=50)           :: fullnm

DATA curfil/'**********************************************'/

!------- Determine if filename has already been entered ----------------

abort = .false.
mistak = .false.

!------- Get a filename from the console -------------------------------

104  fullnm='            '
  filnme='        '
  IF (irdflg == 1.OR.curfil ==  &
        '**********************************************' .OR.mistak) THEN
    WRITE(messag,30303)lname
    30303       FORMAT(' Enter the filename (Maximum of ',i2,' characters)')
    CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1, messag)
    mistak = .false.
  ELSE
    WRITE(*,30301)curfil
    30301       FORMAT(' Current Filename for Save is: ',a46)
    messag=' OK?, Y/N '
    iok=-13
    CALL datain(item,anumbr,1.,0.,iok,line,report,2,0,0,messag)
    CALL reedt(item,crtflg,edflg,rtn_reedt)
    IF(rtn_reedt == 1) THEN
      GO TO 104
    ELSE IF(rtn_reedt == 2) THEN
      GO TO 102
    END IF
    answer=item
    IF (iok == 0.AND.(answer == 'N'.OR.answer == 'n')) THEN
      messag=' Enter the new filename '
      CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1, messag)
    ELSE IF (answer == 'Y'.OR.answer == 'y'.OR.iok == 0) THEN
      item = curfil
    END IF
  END IF
102     answer=item
IF (edflg == 1.OR.edflg == 4)curfil=item

!--------- Check for ABort request when filename is entered ------------

IF (answer == 'AB'.OR.answer == 'ab') THEN
  messag=' Are you sure you wish to ABORT ?, Y/N'
  CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
  canned=item
  IF (canned == 'Y'.OR.canned == 'y') THEN
    abort = .true.
    RETURN
  END IF
  GO TO 104
END IF

!--------- Detect Illegal filename -------------------------------------
!  This area may require customization for specific systems which allow
!  more or less different characters in a filename.

fullnm=item(1:lname)
DO  i=1,lname
  charct=' '
  charct=fullnm(i:i)
  IF (charct == ' ') GO TO 120
  IF ((charct >= 'A'.AND.charct <= 'Z').OR.  &
        (charct >= 'a'.AND.charct <= 'z').OR.  &
        (charct >= '0'.AND.charct <= '9'.AND.i /= 1)) THEN
    CYCLE
  ELSE
    WRITE(*,30302) i
    30302       FORMAT(' Character #',i1,' in filename is illegal')
    mistak = .true.
    item = ' '
    GO TO 104
  END IF
END DO

!------ Determine three letter file extension --------------------------

120     IF (crtflg == 1.OR.edflg == 1.OR.edflg == 5.OR.edflg == 9) THEN
  extnsn='unt'
ELSE IF (crtflg == 2.OR.edflg == 2.OR.edflg == 6.OR.edflg == 10  &
      .OR.edflg == 11) THEN
  extnsn='blk'
ELSE IF (crtflg == 3.OR.edflg == 3.OR.edflg == 7) THEN
  extnsn='sup'
ELSE IF (crtflg == 4.OR.edflg == 4.OR.edflg == 8) THEN
  extnsn='sim'
ELSE
  extnsn = '???'
END IF

!--------- Merge extension name with file name -------------------------
!  This area may require customization for specific systems which have
!  different lengths and specifications for filenames. The extension is
!  used as a unique identifier for the type of work file. Currently
!  configured with a 46 letter filename followed by a period and a
!  three character extension.

IF (i == lname.AND.charct /= ' ') i = i + 1
fullnm(i:i) = delim
IF (delim == ' ') i = i - 1
fullnm(i+1:i+lext) = extnsn

!-------- File OPEN section for reading from files ---------------------
!   Opens a file with name FULLNM.

IF (irdflg == 1) THEN
  OPEN(luf,ERR=300,FILE=fullnm,STATUS='OLD')
  RETURN
ELSE
  OPEN(luf,FILE=fullnm)
  RETURN
END IF

300      CALL rite(' Cannot open a file with the name specified')
item = ' '
GO TO 104
END SUBROUTINE opnfil
!=======================================================================

!     SUBROUTINE RDUNT

!------ This subroutine REads a UNit under the file name given
!------ Variables are defined in the create unit subroutine CRUNT
!------ and open file module OPNFIL.

!=======================================================================

SUBROUTINE rdunt(crtflg,edflg,rtn_rdunt)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=1)    :: answer
INTEGER              :: ijunk

REWIND luf

!------ If EDFLG is 9 then the call was from the module used to insert
!       units in simulations.

IF (edflg == 9) THEN
  READ(luf,*,ERR=400,END=400) ijunk,icheck
  
  110       FORMAT(16I5)
!  110       FORMAT(20I4)
  
!------ Unit number passed by the insert unit subroutine if EDFLG=9.
  
  untyp(UNIT)=icheck
  GO TO 109
END IF
READ(luf,*,ERR=400,END=400) UNIT,icheck

!------ Call made to subroutine TYPEIN  to obtain the the information
!       for NSAVED,IUDE,NNIN,NNOUT,NNPAR,DSCRIB from the types data file
!       TYPAR.DAT.

109     CALL typein(2)
IF (istat /= 0) GO TO 400

blkin(UNIT)=nnin
blkout(UNIT)=nnout
blkpar(UNIT)=nnpar

!------ Read the information from the UNit file.

115     IF (edflg == 9) THEN
  messag=' Set unit variable indices to zero (Y/N) ?'
  CALL datain(item,anumbr,1.,0.,iok,line,report,2,0,0,messag)
  CALL reedt(item,crtflg,edflg,rtn_reedt)
  IF(rtn_reedt == 1) THEN
    GO TO 115
  ELSE IF(rtn_reedt == 2) THEN
    REWIND luf
    CLOSE(luf)
    rtn_rdunt = 3
    RETURN
  END IF
  answer = item
  IF (answer /= 'Y'.AND.answer /= 'y'.AND.answer /= 'N'.AND.  &
      answer /= 'n') GO TO 115
ELSE
  answer = 'N'
END IF
READ(luf,*,ERR=400,END=400) (indxin(UNIT,j),j=1,nnin)
READ(luf,*,ERR=400,END=400) (indxot(UNIT,i),i=1,nnout)
IF (answer == 'Y'.OR.answer == 'y') THEN
  DO  i = 1,nnout
    indxot(UNIT,i) = 0
  END DO
  DO  j = 1,nnin
    indxin(UNIT,j) = 0
  END DO
END IF
READ(luf,*,ERR=400,END=400) (pvalue(UNIT,k),k=1,nnpar)
111     FORMAT(5E15.6)
REWIND luf
CLOSE (luf)
rtn_rdunt = 1
RETURN
400     WRITE(*,30401) icheck
30401   FORMAT(' ERROR reading unit type ',i4,' information from file')
END SUBROUTINE rdunt
!=======================================================================

!     SUBROUTINE RDBLK

!-----  This subroutine reads in information for use by the  INSBLK
!       and VEWBLK (insert and view block subroutines.

!=======================================================================

SUBROUTINE rdblk(crtflg,edflg,rtn_rdblk)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=1)        :: answer
INTEGER                  :: ioff
INTEGER,dimension(8)     :: maxvc
!EQUIVALENCE (pres,maxvc(1))

maxvc = pres                    !  matrix operation    8/16/2006

!------ If EDFLG is set to 10 by the insert block module then
!       the information is read into temporary storage for use
!       by that subroutine.

REWIND luf

ioff = 0
IF (edflg == 10.OR.edflg == 6.OR.edflg == 11) THEN
  READ(luf,*,ERR=400,END=400) nunct
  110       FORMAT(16I5)
!  110       FORMAT(20I4)
  IF (edflg == 10) THEN
    115         messag=' Automatic renumbering of variable indices (Y/N) ? '
    CALL datain(item,anumbr,1.,0.,iok,line,report,2,0,0,messag)
    CALL reedt(item,crtflg,edflg,rtn_reedt)
    IF(rtn_reedt == 1) THEN                     
      GO TO 115                                 
    ELSE IF(rtn_reedt == 2) THEN                
      REWIND luf                                
      CLOSE(luf)                                
      rtn_rdblk = 3
      RETURN                                    
    END IF                                      
    answer = item
    IF (answer /= 'N'.AND.answer /= 'n'.AND.  &
        answer /= 'Y'.AND.answer /= 'y') GO TO 115
  ELSE
    answer='N'
  END IF
  DO  i=1,nunct
    READ(luf,*,ERR=400,END=400)numtem(i),typtem(i)
    
!------ A call is made to the subroutine TYPEIN obtain the
!       information needed from the types data file.
    
    icheck=typtem(i)
    CALL typein(3)
    IF (istat /= 0) GO TO 400
    tin(i)=nnin
    tout(i)=nnout
    tpar(i)=nnpar
    
!------ The information is read in from the file.
    
    READ(luf,*,ERR=400,END=400)(xintem(i,j),j=1,nnin)
    READ(luf,*,ERR=400,END=400)(xottem(i,j),j=1,nnout)
    IF (answer == 'Y'.OR.answer == 'y') THEN
      DO  j = 1,nnin
        IF (xintem(i,j) /= 0) THEN
          xintem(i,j) = xintem(i,j) + maxvc(iident(j))
        END IF
      END DO
      DO  j = 1,nnout
        IF (xottem(i,j) /= 0) THEN
          xottem(i,j) = xottem(i,j) + maxvc(iodent(j))
        END IF
      END DO
    END IF
    READ(luf,*,ERR=400,END=400)(valtem(i,j),j=1,nnpar)
    111         FORMAT(5E15.6)
  END DO
ELSE
  REWIND luf
  CLOSE(luf)
  rtn_rdblk = 2
  RETURN
END IF
REWIND luf
CLOSE(luf)
rtn_rdblk = 1
RETURN
400     WRITE(*,30401) icheck
30401   FORMAT(' ERROR reading unit type ',i4,' information from file')
410     REWIND luf
END SUBROUTINE rdblk
!=======================================================================

!     SUBROUTINE RDSIM

!        This subroutine reads the simulation information that
!        is created with the create SImulation module.
!------  STCNT is the variable that counts the number of variables in
!        the state vector.
!------  TITLSM is used for the simulation title
!------  DECSUP contains the number of superblocks in the simulation.
!------  RTOLX,ATOLX,XTOL, and TTIME are varibles that contain the
!        error tolerances.
!------  SIMCNT counts the number of blocks per superblock.
!------  BLKUNT counts the number of units per block.
!------  UNTNUM array containing the unit numbers.
!------  UNTYP array containing the types numbers.
!------  INDXIN indices for the input connections.
!------  INDXOT indices for the output connections.
!------  PVALUE paramter values for the unit.
!------  STCNT total number of variables in the state vector.
!------  STATE vector containg the initial values for the state
!        variables.
!------  SCNT,SCNT2, and SCNT3 are simple counters.

!=======================================================================

SUBROUTINE rdsim(crtflg,edflg,rtn_rdsim)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=1) :: endfil
INTEGER           :: maxsct,lcnt,mcnt

REWIND luf

250     READ(luf,114,ERR=400,END=400)endfil
114     FORMAT(a1)
IF (endfil /= '*') GO TO 250
READ(luf,*,ERR=400,END=400) decsup,decblk,decunt
READ(luf,*,ERR=400,END=400) (simcnt(i),i=1,decsup)
READ(luf,*,ERR=400,END=400) pres,flow,temp,cntr,othr,enrg, powr,ahum
READ(luf,*,ERR=400,END=400) stcnt
READ(luf,*,ERR=400,END=400) (scnt(k),k=1,decsup)
REWIND luf
l=1
lcnt=1
mcnt=1
maxsct=simcnt(1)
 110     FORMAT(16I5)
! 110     FORMAT(20I4)
111     FORMAT(5E15.6)

!------ Reading the new information.

READ(luf,113)titlsm
113     FORMAT(a80)
READ(luf,*) decsup
READ(luf,*) rtolx,atolx,xtol,ttime

DO  nn=1,decsup
  READ(luf,*) simcnt(nn)
  DO  m=lcnt,maxsct
    READ(luf,*)blkunt(m)
    DO  n=1,blkunt(m)
      READ(luf,*)untnum(l),untyp(untnum(l))
      icheck=untyp(untnum(l))
      CALL typein(2)
      IF (istat /= 0) GO TO 400
      blkin(untnum(l))=nnin
      blkout(untnum(l))=nnout
      blkpar(untnum(l))=nnpar
      
!------- Read information for unit.
      
      READ(luf,*)(indxin(untnum(l),j),j=1,blkin(untnum(l)))
      READ(luf,*)(indxot(untnum(l),k),k=1,blkout(untnum(l)))
      READ(luf,*) (pvalue(untnum(l),mm),mm=1,blkpar(untnum(l)))
      l=l+1
    END DO
    mcnt=mcnt+1
  END DO
  lcnt=mcnt
  maxsct=simcnt(nn+1)+(lcnt-1)
END DO

!------- Reading the state vector.

READ(luf,*) (state(i),i=1,stcnt)

!------- Reading the boundary variables from the work file

READ(luf,*)bcnt
READ(luf,*) (bndry(k),k=1,bcnt)

!------- Reading reported variables to the work file

DO  k=1,decsup
  READ(luf,*)scnt(k),repinv(k)
  112       FORMAT(i4,e15.6)
  READ(luf,*) (repvar(k,m),m=1,scnt(k))
  READ(luf,*) (ntyp(k,m),m=1,scnt(k))
  READ(luf,*) (nindx(k,m),m=1,scnt(k))
END DO

READ(luf,*) (freeze(m),m=1,decsup)
READ(luf,*) (scan(m),m=1,decsup)
REWIND luf
CLOSE(luf)
rtn_rdsim = 1
RETURN
400     CALL rite(' ERROR reading simulation information from file  ')
REWIND luf
CLOSE(luf)
rtn_rdsim = 3
RETURN
END SUBROUTINE rdsim
!=======================================================================

!     SUBROUTINE SAVUNT

!------ This subroutine saves a UNit under the file name given and
!------ makes a backup of the old file if one exists.
!------ Variables are defined in the create unit subroutine CRUNT
!------ and open file module OPNFIL.

!=======================================================================

SUBROUTINE savunt(crtflg,edflg,rtn_savunt)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

110     FORMAT(16I5)
!110     FORMAT(20I4)
111     FORMAT(5E15.6)

!------ Write the new information to the UNit file

REWIND luf
WRITE(luf,110) UNIT,icheck
WRITE(luf,110) (indxin(UNIT,j),j=1,nnin)
WRITE(luf,110) (indxot(UNIT,i),i=1,nnout)
WRITE(luf,111) (pvalue(UNIT,k),k=1,nnpar)
ENDFILE luf
CLOSE(luf)
rtn_savunt = 1
RETURN
END SUBROUTINE savunt
!=======================================================================

!     SUBROUTINE SAVBLK

!------ This subroutine save a BLock under the filename given
!------ Variables are defined in the create block, CRBLK subroutine.
!------ UNCT is used to define the number of units in the block.
!=======================================================================

SUBROUTINE savblk(lstunt,crtflg,edflg,rtn_savblk)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

INTEGER, INTENT(IN)      :: lstunt
INTEGER                  :: ucnt,nucnt,nseq,numunt

!------ ICNT is decremented and renamed to obtain the proper
!       number of units.

ucnt=icnt-1
110     FORMAT(16I5)
! 110     FORMAT(20I4)
111     FORMAT(5E15.6)

!------ Writing the new information to the BLock file

REWIND luf
nucnt = ucnt + lstunt - 1
WRITE(luf,110) ucnt

nseq = 1
DO  m=lstunt,nucnt
  numunt = untnum(m)
  WRITE(luf,110) nseq,untyp(numunt)
  WRITE(luf,110) (indxin(numunt,j),j=1,blkin(numunt))
  WRITE(luf,110) (indxot(numunt,k),k=1,blkout(numunt))
  WRITE(luf,111) (pvalue(numunt,l),l=1,blkpar(numunt))
  nseq = nseq + 1
END DO
ENDFILE luf
CLOSE(luf)
rtn_savblk = 1
RETURN
END SUBROUTINE savblk
!=======================================================================

!      SUBROUTINE SAVSIM

!        This subroutine saves the simulation information that
!        is created with the create SImulation module.
!------- L is a counter that advances the unit number during the
!        save simulation process.
!------- LCNT is used to advance the counter in the loop that
!        controls the variable BLKUNT (number of units per block).
!        This is required so that as the superblock conter is
!        incremented the block counter does not start at 1 again.
!        This is necessary because the block numbers are sequential.
!------- MCNT is a preliminary counter, which LCNT is set equal to
!        because LCNT cannot be incremented within the loop as this
!        would continue to increment the loop counter.
!------- MAXSCT is the variable that sets the maximum limit on the
!        block loop counter, M. MAXSCT is based on SIMCNT (number
!        of blocks per superblock).
!------  STCNT is the variable that counts the number of variables in
!        the state vector.
!------  TITLSM is used for the simulation title
!------  DECSUP contains the number of superblocks in the simulation.
!------  RTOLX,ATOLX,XTOL, and TTIME are varibles that contain the
!        error tolerances.
!------  SIMCNT counts the number of blocks per superblock.
!------  BLKUNT counts the number of units per block.
!------  UNTNUM array containing the unit numbers.
!------  UNTYP array containing the types numbers.
!------  INDXIN indices for the input connections.
!------  INDXOT indices for the output connections.
!------  PVALUE paramter values for the unit.
!------  STCNT total number of variables in the state vector.
!------  STATE vector containg the initial values for the state
!        variables.
!------  SCNT,SCNT2, and SCNT3 are simple counters.
!=======================================================================

SUBROUTINE savsim(crtflg,edflg,rtn_savsim)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

INTEGER                           :: maxsct,lcnt,mcnt
INTEGER,dimension(mdmbnd)         :: bvwrt
INTEGER,dimension(maxsbk,mdmrep)  :: rpwrt

!------- Was the call originally from the edit module or create

IF (edflg == 4) GO TO 109

!------ Counters are decremented to obtain proper count.

decblk=simblk-1
decunt=simunt-1
decsup=simsup-1
109     l=1
lcnt=1
mcnt=1
maxsct=simcnt(1)
110     FORMAT(16I5)
! 110     FORMAT(20I4)
111     FORMAT(5E15.6)

!------ Writing the new information to a superblock file.

WRITE(luf,113)titlsm
113     FORMAT(a80)
WRITE(luf,3110) decsup
3110    FORMAT(i4,t50,'(Superblocks in simulation)')
WRITE(luf,3111) rtolx,atolx,xtol,ttime
3111    FORMAT(4E15.6,t65,'(Error Tol.)')
DO  nn=1,decsup
  WRITE(luf,3112) simcnt(nn), nn
  3112      FORMAT(i4,t50,'(Blocks in SB#',i2,')')
  DO  m=lcnt,maxsct
    WRITE(luf,3113) blkunt(m), m
    3113        FORMAT(i4,t50,'(Units in BLK#',i2,')')
    DO  n=1,blkunt(m)
      WRITE(luf,3114) untnum(l),untyp(untnum(l))
      3114          FORMAT(2I4,t10,40('-'),'(Unit #,Type #)',10('-'))
      WRITE(luf,110)(indxin(untnum(l),j),j=1,blkin(untnum(l)))
      WRITE(luf,110)(indxot(untnum(l),k),k=1,blkout(untnum(l)))
      WRITE(luf,111) (pvalue(untnum(l),mm),mm=1,blkpar(untnum(l)))
      l=l+1
    END DO
    mcnt=mcnt+1
  END DO
  lcnt=mcnt
  maxsct=simcnt(nn+1)+(lcnt-1)
END DO

!------- Writing the state vector to the work file.

stcnt=pres+flow+temp+cntr+othr+enrg+powr+ahum
IF (edflg == 4) GO TO 809
WRITE(luf,111) (state(i),i=1,stcnt)

!------- Writing the boundary variables to the work file

bcnt=0
DO  j=1,stcnt
  IF (bndry(j) > 0) bcnt=bcnt+1
  IF (bndry(j) > 0) bvwrt(bcnt)=bndry(j)
END DO
WRITE(luf,3115) bcnt
3115    FORMAT(i4,t30,'(Boundary Variables in simulation)')
WRITE(luf,110) (bvwrt(k),k=1,bcnt)

!------- Writing reported variables to the work file

DO  k=1,simsup-1
  scnt(k)=0
  DO  l=1,stcnt
    IF (repvar(k,l) > 0) scnt(k)=scnt(k)+1
    IF (repvar(k,l) > 0) rpwrt(k,scnt(k))=repvar(k,l)
  END DO
  WRITE(luf,3116) scnt(k),repinv(k), k
  3116      FORMAT(i4,e15.6,t25,'(Reported Var. & Interval in SB#',i2,')')
  WRITE(luf,110) (rpwrt(k,m),m=1,scnt(k))
  scnt2(k)=0
  DO  l=1,stcnt
    IF (repvar(k,l) > 0) scnt2(k)=scnt2(k)+1
    IF (repvar(k,l) > 0) rpwrt(k,scnt2(k))=ntyp(k,l)
  END DO
  WRITE(luf,110) (rpwrt(k,m),m=1,scnt2(k))
  scnt3(k)=0
  DO  l=1,stcnt
    IF (repvar(k,l) > 0) scnt3(k)=scnt3(k)+1
    IF (repvar(k,l) > 0) rpwrt(k,scnt3(k))=nindx(k,l)
  END DO
  WRITE(luf,110) (rpwrt(k,m),m=1,scnt3(k))
END DO
WRITE(luf,3117) (freeze(m),m=1,decsup)
3117    FORMAT(20i4)                                   ! 4/4/07
! 3117    FORMAT(10I4,t50,'(Freezing Options)')
WRITE(luf,3117) (scan(m),m=1,decsup)
! 3118    FORMAT(10I4,t50,'(Input Scan Options)')

GO TO 802

!------- If call was from Edit module a special save is performed.
!------- Saving state vector.

809     WRITE(luf,111) (state(i),i=1,stcnt)

!------- Saving Boundary Variables.

WRITE(luf,3115)bcnt
WRITE(luf,110) (bndry(k),k=1,bcnt)

!------- Writing reported variables to the work file

DO  k=1,decsup
  WRITE(luf,3116) scnt(k),repinv(k), k
  WRITE(luf,110) (repvar(k,m),m=1,scnt(k))
  WRITE(luf,110) (ntyp(k,m),m=1,scnt(k))
  WRITE(luf,110) (nindx(k,m),m=1,scnt(k))
END DO
WRITE(luf,3117) (freeze(m),m=1,decsup)  ! 11/8/99
WRITE(luf,3117) (scan(m),m=1,decsup)    ! 11/8/99

802     WRITE(luf,810)
810     FORMAT('**************** SUMMARY OF WORK FILE ***************')
!810     FORMAT('****************OLD****WORK****FILE****ENDS****HERE****&
!   &***********')
WRITE(luf,3119) decsup,decblk,decunt
3119    FORMAT(3i4,t50,'(Superblocks, Blocks, Units)')
WRITE(luf,3120) (simcnt(i),i=1,decsup) ! 11/8/99
3120    FORMAT(20i4 /20i4,t85,'(Blocks in Superblocks)')             ! 3/9/07
! 3120    FORMAT(10I4,t50,'(Blocks in Superblocks)')
WRITE(luf,3121) pres,flow,temp,cntr,othr,enrg,powr,ahum
3121    FORMAT(8i4,t50,'(Variables per category)')
WRITE(luf,3122) stcnt
3122    FORMAT(i4,t50,'(State Variables)')
WRITE(luf,3123) (scnt(k),k=1,decsup)     ! 11/8/99
3123    FORMAT(20i4 /20i4,t85,'(Reported Var. in Superblocks)')      ! 3/9/07
! 3123    FORMAT(10I4,t50,'(Reported Var. in Superblocks)')

ENDFILE luf
CLOSE(luf)
rtn_savsim = 1
RETURN
END SUBROUTINE savsim
!=======================================================================

!     SUBROUTINE TYPEIN(IREQ)

!-- This routine provides HVACGEN with a single connection to the TYPES
!   data file. On input IREQ may have the following values:
!   IREQ = 1 - Get type description only in DSCRIB.
!        = 2 - Get above plus number of inputs (NNIN), number of outputs
!              (NNOUT), and number of parameters (NNPAR).
!        = 3 - Get above plus arrays of input and output variable
!              categories (INTYPE,OTTYPE) and character descriptions
!              (IMSSGE,OMSSGE).
!        = 4   Get above plus parameter number (PARNUM), and parameter
!              description (PMSSGE).

!   On output, ISTAT will indicate the success or failure of the
!   information access. Possible values are:
!   ISTAT = 0 - Information successfully obtained.
!         = 1 - A problem was experienced in obtaining the data
!         = 2 - Information on the type was not found in the file

!=======================================================================

SUBROUTINE typein(ireq)

use hvacsim_par
use hvaccomm
implicit none
INTEGER       :: ireq

CHARACTER (LEN=12)            :: typarf='typar.dat   '

INTEGER,dimension(maxtyp)     :: locate,locpar
CHARACTER (LEN=1)             :: star,pcheck
LOGICAL                       :: first
integer                       :: lineno,itype,iude,nsaved,k,locno,locnop, &
                                 nin,nout,npar

SAVE locate,locpar,first
DATA first/.true./,locate /maxtyp*0/, locpar /maxtyp*0/

!-----------------------------------------------------------------------
! On first subroutine call, open file and initialize type pointer table.
! The filename used in the OPEN statement below is system dependent and
! should be changed to the name of the TYPES data file.

IF (first) THEN
  CALL rite(' INITIALIZING TYPES INFORMATION...')
  OPEN(lutypr,FILE=typarf,STATUS='OLD',ERR=3)
  GO TO 7
  3      OPEN(lutypr,FILE=typarff,STATUS='OLD',ERR=5)
  GO TO 7
  5      WRITE(*,6) typarf, typarff
  6      FORMAT(' Cannot open ',a12,' or ',a50)
  STOP ' Fatal error'
  7      CONTINUE
  OPEN(lufb,STATUS='UNKNOWN',FILE='DIRECT.DAT',ACCESS='DIRECT',  &
      RECL=128,FORM='FORMATTED')            ! direct acess file saved
  REWIND lutypr
  
  lineno=0
  10      READ(lutypr,1000,END=60,ERR=100) star
  IF (star == '*') THEN
    lineno=lineno+1
    WRITE(lufb,1000,REC=lineno) star
    lineno=lineno+1
    READ(lutypr,*,ERR=100) itype,dscrib
    WRITE(lufb,2000,REC=lineno) itype,dscrib
    locate(itype)=lineno
    lineno=lineno+1
    READ(lutypr,*,ERR=100) nsaved,iude,nnin,nnout,nnpar
    WRITE(lufb,3000,REC=lineno) nsaved,iude,nnin,nnout,nnpar
    
    DO  k=1,nnin           ! reversed order 9/4/92
      lineno=lineno+1
      READ(lutypr,*,ERR=100) iident(k),imssge(k)
      WRITE(lufb,2000,REC=lineno) iident(k),imssge(k)
    END DO
    
    READ(lutypr,1000) pcheck ! added for clarity  9/4/92
    
    DO  k=1,nnout          ! reversed order 9/4/92
      lineno=lineno+1
      READ(lutypr,*,ERR=100) iodent(k),omssge(k)
      WRITE(lufb,2000,REC=lineno) iodent(k),omssge(k)
    END DO
    
    lineno=lineno+1
    40        READ(lutypr,1000) pcheck
    IF (pcheck == '#') THEN
      WRITE(lufb,1000,REC=lineno) pcheck
      locpar(itype)=lineno
    ELSE
      GO TO 40
    END IF
    
    IF (nnpar >= 1) THEN
      DO  k=1,nnpar
        lineno=lineno+1
        READ(lutypr,*,ERR=100) parnum(k),pmssge(k)
        WRITE(lufb,2000,REC=lineno) parnum(k),pmssge(k)
      END DO
    END IF
  END IF
  GO TO 10
  60      CLOSE (lutypr)
  first=.false.
END IF

!     Position pointer in TYPES data file

IF (icheck > maxtyp) THEN
  CALL rite(' TYPE NUMBER TOO LARGE')
  RETURN
END IF
locno = locate(icheck)
locnop= locpar(icheck)
IF (locno == 0) THEN
  istat=2
  RETURN
END IF

!     Read information from direct access file

READ(lufb,2000,REC=locno) itype,dscrib
IF (itype /= icheck) THEN
  CALL rite('---- MISMATCH ON TYPE -----')
  CALL holdit
  RETURN
END IF
lineno=locno+1

IF (ireq >= 2) THEN
  READ(lufb,3000,REC=lineno) nsaved,iude,nnin,nnout,nnpar
  nin=nnin
  nout=nnout
  npar=nnpar
END IF

IF (ireq >= 3) THEN
  DO  k=1,nnin                    ! resersed order 9/4/92
    lineno=locno+1+k
    READ(lufb,2000,REC=lineno) iident(k),imssge(k)
  END DO
  
  DO  k=1,nnout                   ! resersed order 9/4/92
    lineno=locno+1+nnin+k
    READ(lufb,2000,REC=lineno) iodent(k),omssge(k)
  END DO
END IF

IF (ireq >= 4 .AND. nnpar >= 1) THEN
  READ(lufb,1000,REC=locnop) pcheck
  
  DO  k=1,nnpar
    lineno=locnop+k
    READ(lufb,2000,REC=lineno) parnum(k),pmssge(k)
  END DO
END IF
istat=0
RETURN

100   CALL rite(' --- ERROR in reading TYPES data file -----')
CALL rite('     Check TYPES data file format is correct  ')

1000  FORMAT(a1)
2000  FORMAT(i5,a80)
3000  FORMAT(5I5)
RETURN
END SUBROUTINE typein

