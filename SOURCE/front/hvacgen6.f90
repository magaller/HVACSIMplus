!=======================================================================
! HVACGEN SOURCE FILE #6 OF 6, VERSION 5.0
!=======================================================================

! hvacgen6.f90

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
! CONTINUATION OF EDIT MODULE (REPLACE, INSERT, DELETE, SAVE)
!=======================================================================

!       SUBROUTINE INSSIM

!------ This subroutine allows the user to insert a UNit, or BLock
!       from a simulation.

!=======================================================================

SUBROUTINE inssim(crtflg,edflg,rtn_inssim)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer              :: ioptn

!------ The insert or replace option is selected

101     CALL rite(' Enter the line number of the desired option:')
CALL rite(' ')
CALL rite(' 1 Insert unit')
CALL rite(' 2 Replace unit')
CALL rite(' 3 Insert block')
CALL rite('   or use carriage return to continue')

messag=' '
CALL datain(item,anumbr,3.,0.,ioptn,line,report,2,0,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 101
ELSE IF(rtn_reedt == 2) THEN
  rtn_inssim = 2
  RETURN
END IF

select case (ioptn)
  case (1)
    CALL insunt(crtflg,edflg,rtn_insunt)
    if(rtn_insunt == 1) then
      go to 101
    elseif(rtn_insunt == 2) then
      rtn_inssim = 2
      return
    endif
  case (2)
    CALL repsim(crtflg,edflg,rtn_repsim)
    if(rtn_repsim == 1) then
      rtn_inssim = 1
      return
    elseif(rtn_repsim == 2) then
      rtn_inssim = 2
      return
    endif
  case (3)
    CALL insblk(crtflg,edflg,rtn_insblk)
    if(rtn_insblk == 1) then
      go to 101
    elseif(rtn_insblk == 2) then
      rtn_inssim = 2
      return
    endif
  case (0)
    rtn_inssim = 1
    return
  case default
end select

END SUBROUTINE inssim
!=======================================================================

!     SUBROUTINE INSUNT

!------ This subroutine is used to insert a unit in a simulation.

!=======================================================================

SUBROUTINE insunt(crtflg,edflg,rtn_insunt)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer              :: nblk,nposit,nsup,nubmax,nshift,insert

100     messag=' Enter the Block number where the unit is to be inserted'
CALL datain(item,anumbr,REAL(maxblk),1.,nblk,line,report, 2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_insunt = 2
  RETURN
END IF
102     WRITE(*,1)nblk
1       FORMAT(' Insert a unit in Block',i2)
CALL okay(item,anumbr,10.,1.,intnum,same,abort)
IF (abort) then
  rtn_insunt = 2
  RETURN
endif
IF (.NOT.same) THEN
  nblk=intnum
  GO TO 102
END IF

!----- A call is made to find the superblock location of the block
!      where the Unit is to be inserted.

CALL inschk(nblk,nposit,nsup,l)
IF (nposit == 0) THEN
  CALL rite(' That block is not in the simulation ')
  GO TO 100
END IF
IF (blkunt(nposit) == muntib) THEN
  nubmax=muntib
  WRITE(*,50901) nubmax
  50901     FORMAT(' There are ',i2,' units in this Block',  &
      ' and no more may be added !')
  GO TO 300
END IF

!----- UNTNUM array is shifted to form a hole for the new unit.

nshift=1
DO  i=1,nblk
  nshift=blkunt(i)+nshift
END DO
DO  i=decunt+1,nshift,-1
  untnum(i)=untnum(i-1)
END DO
untnum(nshift)=0

!----- Total unit number in the block is incremented.

blkunt(nblk)=blkunt(nblk)+1

!----- Maximum unit number in the simulation is found.

insert=0
DO  i=1,decunt+1
  IF (untnum(i) > insert) insert=untnum(i)
END DO
insert=insert+1
untnum(nshift)=insert

!----- This value is incremented by one reflecting the inserted unit.

108     CALL rite(' Select one of the following:  ')
CALL rite(' ')
CALL rite(' 1 Read the Unit from a file.  ')
CALL rite(' 2 Create the new Unit at this time. ')

messag=' '
CALL datain(item,anumbr,2.,1.,inswer,line,report,2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 108
ELSE IF(rtn_reedt == 2) THEN
  rtn_insunt = 2
  RETURN
END IF

IF (inswer == 1) THEN
  UNIT=untnum(nshift)
  
!---- EDFLG set to 9 for the call to the read unit module.
  
  edflg=9
  CALL rite(' ')
  CALL rite(' What is the name of the unit work file?')
  CALL readin(crtflg,edflg,rtn_readin)
  IF(rtn_readin == 2) THEN
    rtn_insunt = 2
    RETURN
  END IF
  109       edflg=4
  WRITE(*,4)insert
  4         FORMAT( '  UNIT number ',i4,' is being inserted',//,  &
      '  ---> Note: The inputs and outputs for this unit should',/  &
      ,' be checked to be certain they are correct. This also holds',/  &
      ,' true for the state, boundary and reported variables.')
  GO TO 300
  
ELSE IF (inswer == 2) THEN
  UNIT=untnum(nshift)
  edflg=1
  CALL crunit(crtflg,edflg,rtn_crunit,STATUS)
  IF(rtn_crunit == 2) THEN
    rtn_insunt = 2
    RETURN
  END IF
  110       edflg=4
  IF (STATUS == 1) then
    rtn_insunt = 2
    RETURN
  endif
  WRITE(*,4) insert
  GO TO 300
END IF

GO TO 108

!----- Total units in the simulation is incremented.

300     decunt=decunt+1

!----- The subroutine that computes the new state vector is called.

CALL recalc(rtn_recalc)
IF(rtn_recalc == 1) THEN
  rtn_insunt = 2
  RETURN
endif

CALL fsave(crtflg,edflg,rtn_fsave)
if(rtn_fsave == 1) then
  rtn_insunt = 1
elseif(rtn_fsave == 2) then
  rtn_insunt = 2
endif

RETURN
END SUBROUTINE insunt
!=======================================================================

!      SUBROUTINE INSCHK

!----- This subroutine returns the number of units in the Block
!      where the insert is taking place along with the superblock
!      location.

!=======================================================================

SUBROUTINE inschk(nblk,m,nn,l)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

INTEGER, INTENT(IN OUT)                   :: nblk
integer                                   :: lcnt,mcnt,maxsct

l=1
lcnt=1
mcnt=1
maxsct=simcnt(1)

DO  nn=1,decsup
  DO  m=lcnt,maxsct
    IF (nblk == m) GO TO 200
    DO  n=1,blkunt(m)
      l=l+1
    END DO
    mcnt=mcnt+1
  END DO
  lcnt=mcnt
  maxsct=simcnt(nn+1)+(lcnt-1)
END DO

m=0
200    RETURN
END SUBROUTINE inschk
!=======================================================================

!     SUBROUTINE INSBLK

!------ This subroutine is used to insert blocks into a simulation.

!=======================================================================

SUBROUTINE insblk(crtflg,edflg,rtn_insblk)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer              :: nsup,nposit,numblk,numunt,inunt,mdec,mfinal

100     messag=  &
    ' Enter the Superblock number where the block is to be inserted'
CALL datain(item,anumbr,REAL(maxsbk),1.,nsup,line,report, 2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_insblk = 2
  RETURN
END IF
112     WRITE(*,1)nsup
1       FORMAT(' Insert a block in Superblock',i2)
CALL okay(item,anumbr,10.,1.,intnum,same,abort)
if(abort) then
  rtn_insblk = 2
  RETURN
endif
IF (.NOT.same) THEN
  nsup=intnum
  GO TO 112
END IF

!----- A call is made to find the superblock and unit information.
call inck2(nsup,nposit)
IF (nposit == 0) THEN
  CALL rite(' That superblock is not in the simulation')
  GO TO 100
END IF
IF (simcnt(nposit) >= mblkis) THEN
  WRITE(*,2)mblkis
  2         FORMAT(' There are ',i3,  &
      ' blocks in this Superblock and no more may be added !')
  GO TO 300
END IF

!----- The type of insert is selected by the user.

101     CALL rite(' Select one of the following: ')
CALL rite(' ')
CALL rite(' 1 Read the Block from a file. ')
CALL rite(' 2 Create the new Block at this time.')

messag=' '
CALL datain(item,anumbr,2.,1.,inswer,line,report,2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 101
ELSE IF(rtn_reedt == 2) THEN
  rtn_insblk = 2
  RETURN
END IF
IF (inswer == 2) THEN
  102       messag= ' Enter the number of units in the Block to be inserted'
  CALL datain(item,anumbr,REAL(maxunt),1.,nunct,line,report, 2,3,0,messag)
  CALL reedt(item,crtflg,edflg,rtn_reedt)
  IF(rtn_reedt == 1) THEN                      
    GO TO 102                                  
  ELSE IF(rtn_reedt == 2) THEN                 
    rtn_insblk = 2                             
    RETURN                                     
  END IF                                       
  GO TO 111
ELSE
  
!----- EDFLG is set to 10 indicating a call from the INSBLK module.
  
  edflg=10
  CALL rite(' What is the block that is to be read in ?')
  CALL readin(crtflg,edflg,rtn_readin)
  IF(rtn_readin == 1) THEN
    GO TO 111
  ELSE IF(rtn_readin == 2) THEN
    rtn_insblk = 2                             
    RETURN                                     
  END IF                                       
END IF

!----- Number of blocks up to and including the superblock
!      location of the inserted block is computed.

111     numblk=0
DO  i=1,nposit
  numblk=numblk+simcnt(i)
END DO

!----- Number of units up to and including the insertion point is
!      computed.

numunt=0
DO  i=1,numblk
  numunt=numunt+blkunt(i)
END DO

!----- Total number of blocks in the specified superblock are
!      incremented.

simcnt(nposit)=simcnt(nposit)+1

!----- Maximum unit number in the simulation is determined.

inunt=0
DO  i=1,decunt
  IF (untnum(i) > inunt) inunt=untnum(i)
END DO

!----- Total number of blocks in the simulation is incremented.

decblk=decblk+1

!----- Old unit count is saved before it is incremented.

mdec=decunt

!----- Total number of units in the simulation is incremented
!      depending on the value assigned to NUNCT.

decunt=decunt+nunct

!----- Array containing the number of units in a block is shifted
!      forming a hole for the new block.

DO  i=decblk,numblk+2,-1
  blkunt(i)=blkunt(i-1)
END DO

!----- Total number of units in the new block are assigned

blkunt(numblk+1)=nunct

!----- UNTNUM array is shifted allowing for the new units.

mfinal=mdec-numunt-1
DO  i=decunt,decunt-mfinal,-1
  untnum(i)=untnum(i-nunct)
END DO

!----- The new units are inserted.

IF (edflg == 10) THEN
  
  DO  i=1,nunct
    untnum(numunt+i)=inunt+1
    UNIT=untnum(numunt+i)
    untyp(UNIT)=typtem(i)
    blkin(UNIT)=tin(i)
    blkout(UNIT)=tout(i)
    blkpar(UNIT)=tpar(i)
    DO  j=1,blkin(UNIT)
      indxin(UNIT,j)=xintem(i,j)
    END DO
    DO  j=1,blkout(UNIT)
      indxot(UNIT,j)=xottem(i,j)
    END DO
    DO  j=1,blkpar(UNIT)
      pvalue(UNIT,j)=valtem(i,j)
    END DO
    inunt=inunt+1
  END DO
  GO TO 120
  
ELSE
  edflg=1
  DO  i=1,nunct
    untnum(numunt+i)=inunt+1
    UNIT=untnum(numunt+i)
    CALL crunit(crtflg,edflg,rtn_crunit,STATUS)
    IF(rtn_crunit == 2) THEN
      rtn_insblk = 2
      RETURN
    END IF
109 IF (STATUS == 1) then
      rtn_insblk = 2
      RETURN
    END IF
    inunt=inunt+1
  END DO
END IF

120     WRITE(*,4) numblk+1,nunct,inunt+1-nunct
4       FORMAT(1X,'Block number',i4,' is being inserted. ',/,  &
    ' It has ',i4,' units beginning with unit number',i4,'.',/,  &
    ' All blocks in the simulation will be renumbered',//,  &
    '  ---> Note: The inputs and outputs for this unit should',/,  &
    ' be checked to be certain they are correct. This also',/,  &
    'holds true for the state, boundary and reported variables.')

edflg=4

!----- The subroutine that computes the new state vector is called.

300     CALL recalc(rtn_recalc)
IF(rtn_recalc == 1) THEN
  rtn_insblk = 2
  RETURN
END IF
CALL fsave(crtflg,edflg,rtn_fsave)
IF(rtn_fsave == 1) THEN
  rtn_insblk = 1
elseif(rtn_fsave == 2) then
  rtn_insblk = 2
endif

RETURN
END SUBROUTINE insblk
!=======================================================================

!      SUBROUTINE INCK2

!----- This subroutine returns superblock and unit information to the
!      subroutine that inserts blocks (INSBLK).

!=======================================================================

SUBROUTINE inck2(nsup,nn)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

INTEGER, INTENT(IN OUT)                   :: nsup
integer                                   :: lcnt,mcnt,maxsct

l=1
lcnt=1
mcnt=1
maxsct=simcnt(1)

DO  nn=1,decsup
  IF (nsup == nn) GO TO 200
  DO  m=lcnt,maxsct
    DO  n=1,blkunt(m)
      l=l+1
    END DO
    mcnt=mcnt+1
  END DO
  lcnt=mcnt
  maxsct=simcnt(nn+1)+(lcnt-1)
END DO

nn=0
200    RETURN
END SUBROUTINE inck2
!=======================================================================

!    SUBROUTINE RECALC

!      This subroutine recalculates the state vector position of the
!      variables, after an insert has occured. This assures that any
!      new variables insert as inputs or outputs are included.
!      This subroutine is used to call two other subroutines that
!      recalculate the state vector positions of the boundary and
!      reported variables.

!=======================================================================

SUBROUTINE recalc(rtn_recalc)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

REAL,dimension(maxstv)  :: nstate
INTEGER                 :: ipres,mpres
INTEGER                 :: iflow,mflow
INTEGER                 :: itemp,mtemp
INTEGER                 :: icntr,mcntr
INTEGER                 :: iothr,mothr
INTEGER                 :: ienrg,menrg
INTEGER                 :: ipowr,mpowr
INTEGER                 :: iahum,mahum

!------ Temporary array values are set to zero.

DO  i=1,maxstv
  nstate(i)=0
END DO

!------ Old values are temporarily stored in alternate variables
!       allowing for the computation of the new number of state
!       variables.

ipres=pres
iflow=flow
itemp=temp
icntr=cntr
iothr=othr
ienrg=enrg
ipowr=powr
iahum=ahum

!------ Maximum variable count assigned.

mpres=pres
mflow=flow
mtemp=temp
mcntr=cntr
mothr=othr
menrg=enrg
mpowr=powr
mahum=ahum
pres=0
flow=0
temp=0
cntr=0
othr=0
enrg=0
powr=0
ahum=0

!------ The new total variable count of each variable is computed.

DO  ii=1,decunt
  CALL typinf(ii,rtn_typinf)
  IF(rtn_typinf == 1) THEN
    rtn_recalc = 1
    RETURN
  END IF
  
  DO  jj=1,blkin(untnum(ii))
    IF (intype(untnum(ii),jj) == 1.AND.indxin(untnum(ii),jj)  &
         >= pres) pres=indxin(untnum(ii),jj)
    IF (intype(untnum(ii),jj) == 2.AND.indxin(untnum(ii),jj)  &
         >= flow) flow=indxin(untnum(ii),jj)
    IF (intype(untnum(ii),jj) == 3.AND.indxin(untnum(ii),jj)  &
         >= temp)  temp=indxin(untnum(ii),jj)
    IF (intype(untnum(ii),jj) == 4.AND.indxin(untnum(ii),jj)  &
         >= cntr) cntr=indxin(untnum(ii),jj)
    IF (intype(untnum(ii),jj) == 5.AND.indxin(untnum(ii),jj)  &
         >= othr) othr=indxin(untnum(ii),jj)
    IF (intype(untnum(ii),jj) == 6.AND.indxin(untnum(ii),jj)  &
         >= enrg) enrg=indxin(untnum(ii),jj)
    IF (intype(untnum(ii),jj) == 7.AND.indxin(untnum(ii),jj)  &
         >= powr) powr=indxin(untnum(ii),jj)
    IF (intype(untnum(ii),jj) == 8.AND.indxin(untnum(ii),jj)  &
         >= ahum) ahum=indxin(untnum(ii),jj)
  END DO
  DO  jj=1,blkout(untnum(ii))
    IF (ottype(untnum(ii),jj) == 1.AND.indxot(untnum(ii),jj)  &
         >= pres) pres=indxot(untnum(ii),jj)
    IF (ottype(untnum(ii),jj) == 2.AND.indxot(untnum(ii),jj)  &
         >= flow) flow=indxot(untnum(ii),jj)
    IF (ottype(untnum(ii),jj) == 3.AND.indxot(untnum(ii),jj)  &
         >= temp) temp=indxot(untnum(ii),jj)
    IF (ottype(untnum(ii),jj) == 4.AND.indxot(untnum(ii),jj)  &
         >= cntr) cntr=indxot(untnum(ii),jj)
    IF (ottype(untnum(ii),jj) == 5.AND.indxot(untnum(ii),jj)  &
         >= othr) othr=indxot(untnum(ii),jj)
    IF (ottype(untnum(ii),jj) == 6.AND.indxot(untnum(ii),jj)  &
         >= enrg) enrg=indxot(untnum(ii),jj)
    IF (ottype(untnum(ii),jj) == 7.AND.indxot(untnum(ii),jj)  &
         >= powr) powr=indxot(untnum(ii),jj)
    IF (ottype(untnum(ii),jj) == 8.AND.indxot(untnum(ii),jj)  &
         >= ahum) ahum=indxot(untnum(ii),jj)
  END DO
END DO

IF (mpres > pres) mpres=pres
IF (mflow > flow) mflow=flow
IF (mtemp > temp) mtemp=temp
IF (mcntr > cntr) mcntr=cntr
IF (mothr > othr) mothr=othr
IF (menrg > enrg) menrg=enrg
IF (mpowr > powr) mpowr=powr
IF (mahum > ahum) mahum=ahum

DO  i=1,mpres
  nstate(i)=state(i)
END DO
DO  i=1,mflow
  nstate(pres+i)=state(i+ipres)
END DO
DO  i=1,mtemp
  nstate(pres+flow+i)=state(i+ipres+iflow)
END DO
DO  i=1,mcntr
  nstate(pres+flow+temp+i)=state(i+ipres+iflow+itemp)
END DO
DO  i=1,mothr
  nstate(pres+flow+temp+cntr+i)=state(i+ipres+iflow+itemp+icntr)
END DO
DO  i=1,menrg
  nstate(pres+flow+temp+cntr+othr+i)=state(i+ipres+iflow+itemp+ icntr+iothr)
END DO
DO  i=1,mpowr
  nstate(pres+flow+temp+cntr+othr+enrg+i)=state(i+ipres+iflow+  &
      itemp+icntr+iothr+ienrg)
END DO
DO  i=1,mahum
  nstate(pres+flow+temp+cntr+othr+enrg+powr+i)=state(i+ipres+  &
      iflow+itemp+icntr+iothr+ienrg+ipowr)
END DO
DO  i=1,pres+flow+temp+cntr+othr+enrg+powr+ahum
  state(i)=nstate(i)
END DO

!----- Call to subroutine REBND, used to recalculate state vector
!      positions for the boundary variables.

CALL rebnd(ipres,iflow,itemp,icntr,iothr,ienrg,ipowr,iahum)

!----- Call to subroutine REREPT used to recalculate state vector
!      positions for the reported variables.

CALL rerept(ipres,iflow,itemp,icntr,iothr,ienrg,ipowr,iahum)
RETURN
END SUBROUTINE recalc
!=======================================================================

!     SUBROUTINE TYPINF

!----- This subroutine provides the information on the type of each
!      of the input and output indeces. This is used to recalculate
!      the total number of each variable.

!=======================================================================

SUBROUTINE typinf(ii,rtn_typinf)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

icheck = untyp(untnum(ii))
CALL typein(3)
IF (istat /= 0) GO TO 102

DO  i=1,nnout
  ottype(untnum(ii),i) = iodent(i)
END DO
DO  i=1,nnin
  intype(untnum(ii),i) = iident(i)
END DO
RETURN
102     CALL rite(' Error reading TYPES file. ')
CALL holdit
rtn_typinf = 1
RETURN
END SUBROUTINE typinf
!=======================================================================

!     SUBROUTINE REBND
!------- This subroutine calculates the new position in the state
!        vector of the boundry variables.

!=======================================================================

SUBROUTINE rebnd(ipres,iflow,itemp,icntr,iothr,ienrg,ipowr, iahum)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

INTEGER, INTENT(IN)                      :: ipres
INTEGER, INTENT(IN)                      :: iflow
INTEGER, INTENT(IN)                      :: itemp
INTEGER, INTENT(IN)                      :: icntr
INTEGER, INTENT(IN)                      :: iothr
INTEGER, INTENT(IN)                      :: ienrg
INTEGER, INTENT(IN)                      :: ipowr
INTEGER, INTENT(IN OUT)                  :: iahum
INTEGER,dimension(mdmbnd)                :: nbndry
integer                                  :: max

DO  k=1,bcnt
  
  IF (bndry(k) > pres.AND.bndry(k) <= ipres) THEN
    nbndry(k)=0
    CYCLE
  END IF
  
  IF (bndry(k) > ipres.OR.ipres <= 0) GO TO 210
  nbndry(k)=bndry(k)
  CYCLE
  
  210     IF ((bndry(k)-ipres) > flow.AND.(bndry(k)-ipres) <= iflow) THEN
    nbndry(k)=0
    CYCLE
  END IF
  
  IF (bndry(k) > (ipres+iflow).OR.iflow <= 0) GO TO 220
  nbndry(k)=(bndry(k)-ipres)+pres
  CYCLE
  
  220     IF ((bndry(k)-ipres-iflow) > temp.AND.(bndry(k)-ipres-iflow)  &
         <= itemp) THEN
    nbndry(k)=0
    CYCLE
  END IF
  
  IF (bndry(k) > (ipres+iflow+itemp).OR.itemp <= 0) GO TO 230
  nbndry(k)=(bndry(k)-ipres-iflow)+pres+flow
  CYCLE
  
  230     IF ((bndry(k)-ipres-iflow-itemp) > cntr.AND.(bndry(k)-ipres-  &
        iflow-itemp) <= icntr) THEN
    nbndry(k)=0
    CYCLE
  END IF
  
  IF (bndry(k) > (ipres+iflow+itemp+icntr).OR.icntr <= 0) GO TO 240
  nbndry(k)=(bndry(k)-ipres-iflow-itemp)+pres+flow+temp
  CYCLE
  
  240     IF ((bndry(k)-ipres-iflow-itemp-icntr) > othr.AND.(bndry(k)  &
        -ipres-iflow-itemp-icntr) <= iothr) THEN
    nbndry(k)=0
    CYCLE
  END IF
  
  IF (bndry(k) > (ipres+iflow+itemp+icntr+iothr).OR.iothr <= 0) GO TO 250
  nbndry(k)=(bndry(k)-ipres-iflow-itemp-icntr)+pres+flow+ temp+cntr
  CYCLE
  
  250     IF ((bndry(k)-ipres-iflow-itemp-icntr-iothr) > enrg.AND.  &
        (bndry(k)-ipres-iflow-itemp-icntr-iothr) <= ienrg) THEN
    nbndry(k)=0
    CYCLE
  END IF
  
  IF (bndry(k) > (ipres+iflow+itemp+icntr+iothr+ienrg).OR.ienrg  &
       <= 0) GO TO 260
  nbndry(k)=(bndry(k)-ipres-iflow-itemp-icntr-iothr)+pres+flow+ temp+cntr+othr
  CYCLE
  
  260     IF ((bndry(k)-ipres-iflow-itemp-icntr-iothr-ienrg) > powr.AND.  &
        (bndry(k)-ipres-iflow-itemp-icntr-iothr-ienrg) <= ipowr) THEN
    nbndry(k)=0
    CYCLE
  END IF
  
  IF (bndry(k) > (ipres+iflow+itemp+icntr+iothr+ienrg+ipowr).OR.  &
      ipowr <= 0) GO TO 270
  nbndry(k)=(bndry(k)-ipres-iflow-itemp-icntr-iothr-ienrg)+pres+  &
      flow+temp+cntr+othr+enrg
  CYCLE
  
  270     IF ((bndry(k)-ipres-iflow-itemp-icntr-iothr-ienrg-ipowr) >  &
        ahum.AND.(bndry(k)-ipres-iflow-itemp-icntr-iothr-ienrg-ipowr)  &
         <= iahum) THEN
    nbndry(k)=0
    CYCLE
  END IF
  
  IF (bndry(k) > (ipres+iflow+itemp+icntr+iothr+ienrg+ipowr+  &
      iahum).OR.iahum <= 0) CYCLE
  nbndry(k)=(bndry(k)-ipres-iflow-itemp-icntr-iothr-ienrg-ipowr)+  &
      pres+flow+temp+cntr+othr+enrg+powr
  
END DO

MAX=bcnt
297     DO  i=1,MAX
  IF (nbndry(i) == 0.AND.i <= bcnt) THEN
    bcnt=bcnt-1
    DO  j=i,MAX
      nbndry(j)=nbndry(j+1)
    END DO
    GO TO 297
  END IF
END DO

DO  i=1,bcnt
  bndry(i)=nbndry(i)
END DO

RETURN
END SUBROUTINE rebnd
!=======================================================================

!     SUBROUTINE REREPT
!------- This subroutine calculates the new position in the state
!        vector of the reported variables.
!=======================================================================

SUBROUTINE rerept(ipres,iflow,itemp,icntr,iothr,ienrg,ipowr, iahum)

use hvacsim_par
use hvaccomm
implicit none

INTEGER, INTENT(IN)                      :: ipres
INTEGER, INTENT(IN)                      :: iflow
INTEGER, INTENT(IN)                      :: itemp
INTEGER, INTENT(IN)                      :: icntr
INTEGER, INTENT(IN)                      :: iothr
INTEGER, INTENT(IN)                      :: ienrg
INTEGER, INTENT(IN)                      :: ipowr
INTEGER, INTENT(IN OUT)                  :: iahum
INTEGER,dimension(maxsbk,mdmrep)         :: nrepva,nntyp,nnindx
integer                                  :: i,j,k,m,max

DO  k=1,decsup
  DO  m=1,scnt(k)
    
    IF (repvar(k,m) > pres.AND.repvar(k,m) <= ipres) THEN
      nrepva(k,m)=0
      CYCLE
    END IF
    IF (repvar(k,m) > ipres.OR.ipres <= 0) GO TO 210
    nrepva(k,m)=repvar(k,m)
    nnindx(k,m)=nindx(k,m)
    nntyp(k,m)=ntyp(k,m)
    CYCLE
    
    210     IF ((repvar(k,m)-ipres) > flow.AND.(repvar(k,m)-ipres) <=  &
          iflow) THEN
      nrepva(k,m)=0
      CYCLE
    END IF
    IF (repvar(k,m) > (ipres+iflow).OR.iflow <= 0) GO TO 220
    nrepva(k,m)=(repvar(k,m)-ipres)+pres
    nnindx(k,m)=nindx(k,m)
    nntyp(k,m)=ntyp(k,m)
    CYCLE
    
    220     IF ((repvar(k,m)-ipres-iflow) > temp.AND.(repvar(k,m)-ipres-  &
          iflow) <= itemp) THEN
      nrepva(k,m)=0
      CYCLE
    END IF
    IF (repvar(k,m) > (ipres+iflow+itemp).OR.itemp <= 0) GO TO 230
    nrepva(k,m)=(repvar(k,m)-ipres-iflow)+pres+flow
    nnindx(k,m)=nindx(k,m)
    nntyp(k,m)=ntyp(k,m)
    CYCLE
    
    230     IF ((repvar(k,m)-ipres-iflow-itemp) > cntr.AND.(repvar(k,m)-  &
          ipres-iflow-itemp) <= icntr) THEN
      nrepva(k,m)=0
      CYCLE
    END IF
    IF (repvar(k,m) > (ipres+iflow+itemp+icntr).OR.icntr <= 0) GO TO 240
    nrepva(k,m)=(repvar(k,m)-ipres-iflow-itemp)+pres+flow+temp
    nnindx(k,m)=nindx(k,m)
    nntyp(k,m)=ntyp(k,m)
    CYCLE
    
    240     IF ((repvar(k,m)-ipres-iflow-itemp-icntr) > othr.AND.  &
          (repvar(k,m)-ipres-iflow-itemp-icntr) <= iothr) THEN
      nrepva(k,m)=0
      CYCLE
    END IF
    IF (repvar(k,m) > (ipres+iflow+itemp+icntr+iothr).OR.  &
        iothr <= 0) GO TO 250
    nrepva(k,m)=(repvar(k,m)-ipres-iflow-itemp-icntr)+pres+flow+ temp+cntr
    nnindx(k,m)=nindx(k,m)
    nntyp(k,m)=ntyp(k,m)
    CYCLE
    
    250     IF ((repvar(k,m)-ipres-iflow-itemp-icntr-iothr) > enrg.AND.  &
          (repvar(k,m)-ipres-iflow-itemp-icntr-iothr) <= ienrg) THEN
      nrepva(k,m)=0
      CYCLE
    END IF
    IF (repvar(k,m) > (ipres+iflow+itemp+icntr+iothr+ienrg).OR.  &
        ienrg <= 0) GO TO 260
    nrepva(k,m)=(repvar(k,m)-ipres-iflow-itemp-icntr-iothr)+pres  &
        +flow+temp+cntr+othr
    nnindx(k,m)=nindx(k,m)
    nntyp(k,m)=ntyp(k,m)
    CYCLE
    
    260     IF ((repvar(k,m)-ipres-iflow-itemp-icntr-iothr-ienrg) > powr  &
          .AND.(repvar(k,m)-ipres-iflow-itemp-icntr-iothr-ienrg)  <= ipowr) THEN
      nrepva(k,m)=0
      CYCLE
    END IF
    IF (repvar(k,m) > (ipres+iflow+itemp+icntr+iothr+ienrg+ipowr)  &
        .OR.ipowr <= 0) GO TO 270
    nrepva(k,m)=(repvar(k,m)-ipres-iflow-itemp-icntr-iothr-ienrg)  &
        +pres+flow+temp+cntr+othr+enrg
    nnindx(k,m)=nindx(k,m)
    nntyp(k,m)=ntyp(k,m)
    CYCLE
    
    270     IF ((repvar(k,m)-ipres-iflow-itemp-icntr-iothr-ienrg-ipowr) >  &
          ahum.AND.(repvar(k,m)-ipres-iflow-itemp-icntr-iothr-ienrg  &
          -ipowr) <= iahum) THEN
      nrepva(k,m)=0
      CYCLE
    END IF
    IF (repvar(k,m) > (ipres+iflow+itemp+icntr+iothr+ienrg+ipowr+  &
        iahum).OR.iahum <= 0) CYCLE
    nrepva(k,m)=(repvar(k,m)-ipres-iflow-itemp-icntr-iothr-ienrg-  &
        ipowr)+pres+flow+temp+cntr+othr+enrg+powr
    nnindx(k,m)=nindx(k,m)
    nntyp(k,m)=ntyp(k,m)
    
  END DO
END DO

DO  k=1,decsup
  MAX=scnt(k)
  297        DO  i=1,MAX
    IF (nrepva(k,i) == 0.AND.i <= scnt(k)) THEN
      scnt(k)=scnt(k)-1
      scnt2(k)=scnt2(k)-1
      scnt3(k)=scnt3(k)-1
      DO  j=i,MAX
        nrepva(k,j)=nrepva(k,j+1)
        nntyp(k,j)=nntyp(k,j+1)
        nnindx(k,j)=nnindx(k,j+1)
      END DO
      GO TO 297
    END IF
  END DO
END DO

DO  k=1,decsup
  DO  m=1,scnt(k)
    repvar(k,m)=nrepva(k,m)
    ntyp(k,m)=nntyp(k,m)
    nindx(k,m)=nnindx(k,m)
  END DO
END DO
RETURN
END SUBROUTINE rerept
!=======================================================================

!       SUBROUTINE DELSIM

!------ This subroutine allows the user to delete a UNit, BLock,
!       or SUperblock from a simulation.

!=======================================================================

SUBROUTINE delsim(crtflg,edflg,rtn_delsim)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

100     CALL rite(' Enter the line number of the item that is to ')
CALL rite(' be DELETED:')
CALL rite(' ')
CALL rite(' 1 UNit')
CALL rite(' 2 BLock')
!CC     CALL RITE(' 3 SUperblock (Not yet available)')
CALL rite(' ')
CALL rite(' or Carriage Return to return to previous menu    ')

101     messag=' '
CALL datain(item,anumbr,3.,0.,inswer,line,report,2,0,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  go to 101
elseif(rtn_reedt == 2) then
  rtn_delsim = 2
  return
END IF
!        CALL REEDT(ITEM,CRTFLG,EDFLG,*101,*450)
SELECT CASE ( inswer )
  CASE (    1)
    GO TO 110
  CASE (    2)
    GO TO 120
  CASE (    3)
    GO TO 130
END SELECT
IF (inswer == 0) GO TO 300
110     CALL delunt(crtflg,edflg,rtn_delunt)
IF(rtn_delunt == 1) THEN
  go to 100
elseif(rtn_delunt == 2) then
  rtn_delsim = 2
  return
END IF
120     CALL delblk(crtflg,edflg,rtn_delblk)
IF(rtn_delblk == 1) THEN
  go to 100
elseif(rtn_delblk == 2) then
  rtn_delsim = 2
  return
END IF

!------ DELSIM NOT IMPLEMENTED YET

130     GO TO 100
300     rtn_delsim = 1
RETURN
END SUBROUTINE delsim
!=======================================================================

!     SUBROUTINE DELUNT

!------ This subroutine is used to delete a unit from a simulation.

!=======================================================================

SUBROUTINE delunt(crtflg,edflg,rtn_delunt)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer               :: nunt,nblk,nposit,nsup

100     messag=' Enter the Unit number to be deleted'
CALL datain(item,anumbr,REAL(maxunt),1.,nunt,line,report, 2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_delunt = 2
  RETURN
END IF
102     WRITE(*,1)nunt,untyp(nunt)
1       FORMAT(' Delete Unit',i3,' Type',i3)
CALL okay(item,anumbr,REAL(maxunt),1.,intnum,same,abort)
IF (abort) then
  rtn_delunt = 2
  RETURN
endif
IF (.NOT.same) THEN
  nunt=intnum
  GO TO 102
END IF

!----- A call is made to find the superblock and block locations of
!      the unit to be deleted.

call delchk(nunt,nblk,nsup,nposit)
IF (nblk == 0) THEN
  CALL rite(' That unit is not in the simulation')
  GO TO 100
END IF
IF (blkunt(nblk) <= 1) THEN
  CALL rite(' There is only ONE unit in this block and it may not be deleted')
  GO TO 300
END IF
WRITE(*,2)nunt
2       FORMAT(' NOTE: All unit numbers greater than ',i3,' have been',  &
    ' decreased by one.')

!----- Total unit number in the specified block number is decremented.

blkunt(nblk)=blkunt(nblk)-1

!----- Arrays are shifted to fill in the unit that was deleted.

DO  i=nposit,decunt-1
  untnum(i)=untnum(i+1)
END DO

decunt=decunt-1
DO  i=1,200
  DO  l=1,decunt
    IF (untnum(l) == nunt+1) THEN
      untyp(untnum(l)-1)=untyp(untnum(l))
      blkin(untnum(l)-1)=blkin(untnum(l))
      blkout(untnum(l)-1)=blkout(untnum(l))
      blkpar(untnum(l)-1)=blkpar(untnum(l))
      DO  j=1,blkin(untnum(l)-1)
        indxin(untnum(l)-1,j)=indxin(untnum(l),j)
      END DO
      DO  j=1,blkout(untnum(l)-1)
        indxot(untnum(l)-1,j)=indxot(untnum(l),j)
      END DO
      DO  j=1,blkpar(untnum(l)-1)
        pvalue(untnum(l)-1,j)=pvalue(untnum(l),j)
      END DO
      untnum(l)=untnum(l)-1
    END IF
  END DO
  nunt=nunt+1
END DO

300     CALL recalc(rtn_recalc)
IF(rtn_recalc == 1) THEN
  rtn_delunt = 2
  RETURN
END IF
CALL fsave(crtflg,edflg,rtn_fsave)
if(rtn_fsave == 1) then
  rtn_delunt = 1
elseif(rtn_fsave == 2) then
  rtn_delunt = 2
endif

RETURN
END SUBROUTINE delunt
!=======================================================================

!     SUBROUTINE DELBLK

!------ This subroutine is used to delete a block from a simulation.

!=======================================================================

SUBROUTINE delblk(crtflg,edflg,rtn_delblk)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

INTEGER,dimension(maxunt)          :: tmpunt
integer                            :: nblk,mblk,nposit,nunt,mblkun,maxtmp,mcount

100     messag=' Enter the Block number to be deleted'
CALL datain(item,anumbr,REAL(maxblk),1.,nblk,line,report, 2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_delblk = 2
  RETURN
END IF
102     WRITE(*,1)nblk
1       FORMAT(' Delete Block',i3)
CALL okay(item,anumbr,REAL(maxblk),1.,intnum,same,abort)
IF (abort) then
  rtn_delblk = 2
  RETURN
endif
IF (.NOT.same) THEN
  nblk=intnum
  GO TO 102
END IF

!----- A call is made to find superblock and unit information for the
!      block that is being deleted.

call delck2(nblk,mblk,nposit,nunt)
!CALL delck2(simcnt,blkunt,untnum,decsup,nblk,mblk,nposit,nunt)

!----- NPOSIT defines the superblock location of the block.

IF (nposit == 0) THEN
  CALL rite(' That block is not in the simulation')
  GO TO 100
END IF
IF (simcnt(nposit) <= 1) THEN
  CALL rite(' There is only ONE block in this superblock and it may not be deleted !')
  GO TO 300
END IF

!----- The unit numbers that are deleted with the block are stored
!      in a temporary array.

DO  i=nunt,blkunt(nblk)+nunt-1
  tmpunt(i+1-nunt)=untnum(i)
END DO

!----- Arrays are shifted to fill in the units that were deleted when
!      the block was deleted.
!----- Array storing the unit numbers is shifted (UNTNUM).

DO  i=nunt,decunt-blkunt(nblk)
  untnum(i)=untnum(i+blkunt(nblk))
END DO
DO  i=decunt-blkunt(nblk)+1,decunt
  untnum(i)=0
END DO

!----- New maximum unit count is redefined

decunt=decunt-blkunt(nblk)
mblkun=blkunt(nblk)
150     maxtmp=0
mcount = 1
DO  i=1,mblkun
  IF (tmpunt(i) > maxtmp) THEN
    maxtmp=tmpunt(i)
    mcount=i
  END IF
END DO
tmpunt(mcount)=0
IF (maxtmp == 0) GO TO 151

DO  i=1,200
  DO  l=1,decunt
    IF (untnum(l) == maxtmp+1) THEN
      untyp(untnum(l)-1)=untyp(untnum(l))
      blkin(untnum(l)-1)=blkin(untnum(l))
      blkout(untnum(l)-1)=blkout(untnum(l))
      blkpar(untnum(l)-1)=blkpar(untnum(l))
      DO  j=1,blkin(untnum(l)-1)
        indxin(untnum(l)-1,j)=indxin(untnum(l),j)
      END DO
      DO  j=1,blkout(untnum(l)-1)
        indxot(untnum(l)-1,j)=indxot(untnum(l),j)
      END DO
      DO  j=1,blkpar(untnum(l)-1)
        pvalue(untnum(l)-1,j)=pvalue(untnum(l),j)
      END DO
      untnum(l)=untnum(l)-1
    END IF
  END DO
  maxtmp=maxtmp+1
END DO

GO TO 150

!----- The number of blocks in the simulation is decremented.

151     decblk=decblk-1

!----- The array containing the unit count for the blocks is shifted.

DO  i=nblk,decblk
  blkunt(i)=blkunt(i+1)
END DO

!----- The array containing the number of blocks per superblock
!      is adjusted to reflect the correct number of blocks

simcnt(nposit)=simcnt(nposit)-1
300     CALL recalc(rtn_recalc)
IF(rtn_recalc == 1) THEN
  rtn_delblk = 2
  RETURN
END IF
CALL fsave(crtflg,edflg,rtn_fsave)
if(rtn_fsave == 1) then
  rtn_delblk = 1
elseif(rtn_fsave == 2) then
  rtn_delblk = 2
endif
RETURN
END SUBROUTINE delblk
!=======================================================================

!      SUBROUTINE DELCHK

!----- This subroutine checks the unit number that is to be deleted
!      and finds the superblock and block where the unit is located.

!=======================================================================

SUBROUTINE delchk(nunt,m,nn,l)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

INTEGER, INTENT(IN OUT)                   :: nunt
integer                                   :: lcnt,mcnt,maxsct
l=1
lcnt=1
mcnt=1
maxsct=simcnt(1)

DO  nn=1,decsup
  DO  m=lcnt,maxsct
    DO  n=1,blkunt(m)
      IF (untnum(l) == nunt) GO TO 200
      l=l+1
    END DO
    mcnt=mcnt+1
  END DO
  lcnt=mcnt
  maxsct=simcnt(nn+1)+(lcnt-1)
END DO

m=0
200     RETURN
END SUBROUTINE delchk
!=======================================================================

!      SUBROUTINE DELCK2

!----- This subroutine checks the block number that is to be deleted
!      and finds the superblock and unit information required for the
!      deletion module.

!=======================================================================

SUBROUTINE delck2(nblk,m,nn,l)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

INTEGER                                   :: nblk,lcnt,mcnt,maxsct

l=1
lcnt=1
mcnt=1
maxsct=simcnt(1)

DO  nn=1,decsup
  DO  m=lcnt,maxsct
    IF (m == nblk) GO TO 200
    DO  n=1,blkunt(m)
      l=l+1
    END DO
    mcnt=mcnt+1
  END DO
  lcnt=mcnt
  maxsct=simcnt(nn+1)+(lcnt-1)
END DO

nn=0
200     RETURN
END SUBROUTINE delck2
!=======================================================================

!       SUBROUTINE REPSIM

!------ This subroutine allows the user to replace a UNit, BLock,
!       Superblock, in a simulation.

!=======================================================================

SUBROUTINE repsim(crtflg,edflg,rtn_repsim)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CALL repunt(crtflg,edflg,rtn_repunt)
if(rtn_repunt == 1) then
  rtn_repsim = 1
elseif(rtn_repunt == 2) then
  rtn_repsim = 2
endif
RETURN
END SUBROUTINE repsim
!=======================================================================

!     SUBROUTINE REPUNT

!------ This subroutine is used to replace a unit in a simulation.

!=======================================================================

SUBROUTINE repunt(crtflg,edflg,rtn_repunt)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer               :: nunt,nblk,nsup,nposit

100     messag=' Enter the Unit number to be replaced'
CALL datain(item,anumbr,REAL(maxunt),1.,nunt,line,report, 2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_repunt = 2
  RETURN
END IF
102     WRITE(*,1)nunt,untyp(nunt)
1       FORMAT(' Replace Unit',i3,' Type',i3)
CALL okay(item,anumbr,REAL(maxunt),1.,intnum,same,abort)
IF (abort) then
  rtn_repunt = 2
  RETURN
endif
IF (.NOT.same) THEN
  nunt=intnum
  GO TO 102
END IF

!----- A call is made to find the superblock and block locations of
!      the unit to be replaced.

CALL delchk(nunt,nblk,nsup,nposit)
!CALL delchk(simcnt,blkunt,untnum,decsup,nunt,nblk,nsup,nposit)
IF (nblk == 0) THEN
  CALL rite(' That unit is not in the simulation')
  GO TO 100
END IF
UNIT=untnum(nposit)

105     CALL rite(' Select one of the following: ')
CALL rite(' ')
CALL rite(' 1 Read the Unit from a file. ')
CALL rite(' 2 Create the new Unit at this time. ')

messag=' '
CALL datain(item,anumbr,2.,1.,inswer,line,report,2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 105
ELSE IF(rtn_reedt == 2) THEN
  rtn_repunt = 2
  RETURN
END IF

IF (inswer == 1) THEN
  
!---- EDFLG set to 9 for the call to the read unit module.
  
  edflg=9
  CALL rite(' ')
  CALL rite(' What is the unit that is to be read in ?')
  CALL readin(crtflg,edflg,rtn_readin)
  IF(rtn_readin == 2) THEN
    rtn_repunt = 2
    RETURN
  END IF
  106       edflg=4
  CALL rite('  ---> Note: The inputs and outputs for this unit')
  CALL rite(' should be checked to be certain they are')
  CALL rite(' correct. This also holds true for the variable ')
  CALL rite(' initial values, boundary and reported variables.')
  GO TO 300
  
ELSE IF (inswer == 2) THEN
  edflg=1
  CALL crunit(crtflg,edflg,rtn_crunit,STATUS)
  IF(rtn_crunit == 2) THEN
    rtn_repunt = 2
    RETURN
  END IF
  107       edflg=4
  IF (STATUS == 1) then
    rtn_repunt = 2
    RETURN
  endif
  CALL rite('  ---> Note: The inputs and outputs for this unit')
  CALL rite(' should be checked to be certain they are')
  CALL rite(' correct. This also holds true for the variable')
  CALL rite(' initial values, boundary and reported variables.')
  GO TO 300
END IF
GO TO 105

!----- The subroutine that computes the new state vector is called.

300     CALL recalc(rtn_recalc)
IF(rtn_recalc == 1) THEN
  rtn_repunt = 2
  RETURN
END IF
CALL fsave(crtflg,edflg,rtn_fsave)
if(rtn_fsave == 1) then
  rtn_repunt = 1
elseif(rtn_fsave == 2) then
  rtn_repunt = 2
endif
301     RETURN
END SUBROUTINE repunt
!=======================================================================

!       HELP MODULE

!------ This subroutine provides a help listing that may be
!       accessed by almost any portion of the program.

!=======================================================================

SUBROUTINE help(crtflg,edflg,abort)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

lincnt=1
CALL rite(' The currently available commands are as follows: ')
CALL scroll(lincnt)
CALL rite(' ')
CALL scroll(lincnt)
CALL rite(' CReate - this allows the user to create a UNit,  ')
CALL scroll(lincnt)
CALL rite(' BLock, or SImulation. The ability to create')
CALL scroll(lincnt)
CALL rite(' units and blocks separately will allow the user  ')
CALL scroll(lincnt)
CALL rite(' to insert them into a previously constructed ')
CALL scroll(lincnt)
CALL rite(' simulation. (See EDit.)')
CALL scroll(lincnt)
CALL rite(' ')
CALL scroll(lincnt)
CALL rite(' EDit - this provides the user with the ablity to ')
CALL scroll(lincnt)
CALL rite(' make changes to any portion of a UNit or any     ')
CALL scroll(lincnt)
CALL rite(' SImulation. Most importantly it will allow for   ')
CALL scroll(lincnt)
CALL rite(' insertion and deletion of BLocks and UNits into  ')
CALL scroll(lincnt)
CALL rite(' a simulation. The units and blocks may be read   ')
CALL scroll(lincnt)
CALL rite(' from a file or created at the time insertion.    ')
CALL scroll(lincnt)
CALL rite(' ')
CALL scroll(lincnt)
CALL rite(' VIew - the SImulation structure, reported and    ')
CALL scroll(lincnt)
CALL rite(' boundary variables, error tolerances and initial ')
CALL scroll(lincnt)
CALL rite(' variable values may be viewed through the use    ')
CALL scroll(lincnt)
CALL rite(' of this command. The structure of BLocks and     ')
CALL scroll(lincnt)
CALL rite(' UNits may also be viewed with this command.      ')
CALL scroll(lincnt)
CALL rite(' ')
CALL scroll(lincnt)
CALL rite(' ABort - this command will return the user back   ')
CALL scroll(lincnt)
CALL rite(' to the main menu. If editing, the item that was  ')
CALL scroll(lincnt)
CALL rite(' being edited will not be saved. ')
CALL scroll(lincnt)
CALL rite(' ')
CALL scroll(lincnt)
CALL rite(' Note: To implement any of the above commands     ')
CALL scroll(lincnt)
CALL rite(' simply type the CAPTILIZED letters.              ')
100     messag=' Hit the Carriage Return to continue'
CALL datain(item,anumbr,0.,0.,inswer,line,report,2,3,0,messag)
IF (inswer /= 0) GO TO 100
RETURN
END SUBROUTINE help
!=======================================================================

!     SUBROUTINE EXTBLK

!=======================================================================

SUBROUTINE extblk(crtflg,edflg,rtn_extblk)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer               :: nblk,lstunt,mblk,nposit,nunt

100     messag=' Enter the number of the block to be saved'
CALL datain(item,anumbr,REAL(maxblk),1.,nblk,line,report, 2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_extblk=2
  RETURN
END IF
102     WRITE(*,1)nblk
1       FORMAT(' Save Block',i3)
CALL okay(item,anumbr,REAL(maxblk),1.,intnum,same,abort)
IF (abort) THEN
  rtn_extblk=2
  RETURN
ENDIF
IF (.NOT.same) THEN
  nblk=intnum
  GO TO 102
END IF

!----- A call is made to find the first unit and number of units for
!      the block being saved.
CALL delck2(nblk,mblk,nposit,nunt)

!----- NPOSIT defines the superblock location of the block.

IF (nposit == 0) THEN
  CALL rite(' That block is not in the simulation')
  GO TO 100
END IF

!------first unit

lstunt = nunt

!------number of units

icnt = blkunt(nblk) + 1

edflg = 2
irdflg = 0
CALL opnfil(crtflg,edflg,abort)
IF (abort) then
  rtn_extblk=2
  RETURN
endif  
CALL savblk(lstunt,crtflg,edflg,rtn_savblk)
IF(rtn_savblk == 2) THEN
  rtn_extblk=2
  RETURN
END IF
edflg = 4

300  rtn_extblk=1
RETURN
END SUBROUTINE extblk

