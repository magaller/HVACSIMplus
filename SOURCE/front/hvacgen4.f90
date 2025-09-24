!=======================================================================
! HVACGEN SOURCE FILE #4 OF 6, VERSION 5.0
!=======================================================================

! hvacgen4.f90

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

!       VIEW MODULE

!------ This module transfers control to the proper subroutines
!       as indicated by the user's selection of viewing a
!       SImulation,BLock,UNit.
!------ CHOICE is the variable which indicates what the user
!       would like to viewing, i.e., a SImulation, BLock, Unit.

!=======================================================================

SUBROUTINE view(crtflg,edflg,abort)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

!----------------------------------------------------------------------

CHARACTER (LEN=2) :: choice

!messag=' '
!CALL datain(item,anumbr,0.,0.,intnum,line,report,1,0,1,messag)
!choice=item
!IF (choice == ' ') THEN
  100       WRITE(*,40401)
  40401     FORMAT(1X,'View a:',//, ' SImulation',/,' BLock',/, ' UNit')
  messag=' Enter Selection'
  CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
  CALL revew(item,crtflg,edflg,rtn_revew)
  IF(rtn_revew == 1) THEN
    GO TO 100
  ELSE IF(rtn_revew == 2) THEN
    GO TO 300
  END IF
  choice=item

IF (choice == 'SI'.OR.choice == 'si') THEN
  IF (crtflg /= 1.AND.crtflg /= 2.AND.crtflg /= 3.AND. crtflg /= 4) THEN
    CALL vewsim(crtflg,edflg,rtn_vewsim)
    IF(rtn_vewsim == 1) THEN
      GO TO 300
    ELSE IF(rtn_vewsim == 2) THEN
      edflg=0
      crtflg=0
      abort = .true.
      RETURN
    END IF
  END IF
ELSE IF (choice == 'BL'.OR.choice == 'bl') THEN
  IF (crtflg /= 1.AND.crtflg /= 2.AND.crtflg /= 3.AND.crtflg  /= 4) THEN
    CALL vewblk(crtflg,edflg,rtn_vewblk)
    IF(rtn_vewblk == 1) THEN
      GO TO 300
    ELSE IF(rtn_vewblk == 2) THEN
      edflg=0
      crtflg=0
      abort = .true.
      RETURN
    END IF
  END IF
ELSE IF (choice == 'UN'.OR.choice == 'un') THEN
  CALL vewunt(crtflg,edflg,rtn_vewunt)
  IF(rtn_vewunt == 1) THEN
    GO TO 300
  ELSE IF(rtn_vewunt == 2) THEN
    edflg=0
    crtflg=0
    abort = .true.
    RETURN
  END IF
END IF

300     IF (edflg == 5.OR.edflg == 6.OR.edflg == 7.OR.edflg == 8.OR.  &
      edflg == 10.OR.edflg == 11) THEN
  edflg=0
  crtflg=0
ELSE
  RETURN
END IF
!400     EDFLG=0
!        CRTFLG=0
RETURN
END SUBROUTINE view
!=======================================================================

!      SUBROUTINE VEWUNT

!----- This subroutine allows the user to view the inputs, outputs,
!      and parameters for a unit.

!----- NUNIT is the unit number assigned by the user for that is to
!      be viewed. This is set to a default value if viewing from the
!      CReate unit module.

!----- SCROLL is a subroutine that controls the screen scrolling.

!----- LNUMBR is used to count the number of lines written the console.

!----- CRTFLG indicates calls from the unit,block,superblock, and
!      simulation level; 1,2,3, and 4 respectively.

!----- EDFLG is similar but this indicates calls from the edit module.
!      This is also used for calls from the view module, 5,6,7, and 8,
!      indicating unit,block,superblock and simultion levels re-
!      spectively.
!=======================================================================

SUBROUTINE vewunt(crtflg,edflg,rtn_vewunt)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=15) :: iprmpt,oprmpt
integer            :: nnunit,lnumbr

!  If the unit number is negative this implies a minimum of processing
!  because the view all function is being invoked.

IF (edflg < 0) THEN
  nunit = -edflg
  GO TO 103
END IF

!------ Entering the unit number which the user would like to view.

!------ EDFLG is set to 8 by the VIEW module this allows for
!       viewing units in a simulation when call by the view simulation
!       structure subroutine.

IF (edflg == 8) GO TO 100
IF (crtflg == 0.AND.edflg == 0) GO TO 120
IF (blkflg == 1) GO TO 100
IF (supflg == 1) GO TO 100
IF (simflg == 1) GO TO 100

!------ Call from edit a simulation? If so jump to unit number request.

IF (edflg == 4.OR.edflg == 10.OR.edflg == 11) GO TO 100

nunit=1
untnum(nunit)=1
IF (edflg == 1) untnum(nunit)=UNIT
IF (crtflg /= 1.AND.crtflg /= 2.AND.crtflg /= 3.AND.  &
    crtflg /= 4.AND.edflg /= 1.AND.edflg /= 2.AND. edflg /= 3.AND.edflg /= 4)  &
    GO TO 120

GO TO 122
120       IF (crtflg == 0.AND.edflg == 0) edflg=5
CALL readin(crtflg,edflg,rtn_readin)
IF(rtn_readin == 2) THEN
  rtn_vewunt = 2
  RETURN
END IF
121     untnum(nunit)=UNIT

!------ If the call was not from create or edit then the information is
!       only viewed.

122     untyp(untnum(nunit))=icheck
GO TO 103
100     messag= ' Enter the UNit number to view (carriage return to skip)'

IF (edflg == 1.OR.edflg == 4) messag=  &
    ' Enter the UNit number to edit (carriage return to skip)'
CALL datain(item,anumbr,REAL(maxunt),0.,nnunit,line,report, 2,3,0,messag)
CALL revew(item,crtflg,edflg,rtn_revew)
IF(rtn_revew == 1) THEN
  rtn_vewunt = 1
  GO TO 100
ELSE IF(rtn_revew == 2) THEN
  rtn_vewunt = 2
  RETURN
END IF
IF (nnunit == 0) then
  rtn_vewunt = 1
  RETURN
endif

!------ Finding the unit in the simulation

DO  nunit=1,maxunt
  IF (untnum(nunit) == nnunit) GO TO 103
END DO

CALL rite(' That UNit was not found')
GO TO 100

!------ Information is read in from the types data file TYPAR.DAT

103     icheck = untyp(untnum(nunit))
CALL typein(4)
IF (istat == 2) GO TO 104
IF (istat /= 0) GO TO 112

!------ The information for the types is displayed.

WRITE(*,1) untnum(nunit),untyp(untnum(nunit))
IF (svmode) WRITE(luview,1) untnum(nunit),untyp(untnum(nunit))
1       FORMAT(/,' UNIT ',i3,5X,' TYPE ',i3)
lincnt=1
WRITE(*,2) dscrib
IF (svmode)  WRITE(luview,2) dscrib

2       FORMAT(1X,a72)
IF (edflg >= 0) CALL scroll(lincnt)
lnumbr=1
WRITE(*,3)lnumbr
IF (svmode)  WRITE(luview,3)lnumbr
3       FORMAT(/,i2,5X,'INPUTS:')
IF (edflg >= 0) CALL scroll(lincnt)

DO  j=1,blkin(untnum(nunit))
  CALL promptd(iident,j,iprmpt)
  WRITE(*,4) iprmpt,indxin(untnum(nunit),j),imssge(j)
  IF (svmode) WRITE(luview,4) iprmpt,indxin(untnum(nunit),j),imssge(j)
  4         FORMAT(7X,a15,1X,i3,' - ',a50)
  IF (edflg >= 0) CALL scroll(lincnt)
END DO

lnumbr=lnumbr+1
WRITE(*,5)lnumbr
IF (svmode) WRITE(luview,5)lnumbr
5       FORMAT(/,i2,5X,'OUTPUTS:')
IF (edflg >= 0) CALL scroll(lincnt)

DO  i=1,blkout(untnum(nunit))
  CALL promptd(iodent,i,oprmpt)
  WRITE(*,4)oprmpt,indxot(untnum(nunit),i),omssge(i)
  IF (svmode) WRITE(luview,4)oprmpt,indxot(untnum(nunit),i),omssge(i)
  IF (edflg >= 0) CALL scroll(lincnt)
END DO

lnumbr=lnumbr+1
WRITE(*,6)lnumbr
IF (svmode) WRITE(luview,6)lnumbr
6       FORMAT(/,i2,5X,'PARAMETERS:')
IF (edflg >= 0) CALL scroll(lincnt)

DO  k=1,blkpar(untnum(nunit))
  WRITE(*,7) pvalue(untnum(nunit),k),pmssge(k)
  IF (svmode) WRITE(luview,7) pvalue(untnum(nunit),k),pmssge(k)
  7         FORMAT(7X,g15.6,1X,a55)
  IF (edflg >= 0) CALL scroll(lincnt)
END DO

IF (edflg == 8.OR.edflg == 10.OR.edflg == 11) GO TO 100

!------ Continue? If call from the EDit module skip this and return
!       there.

IF (edflg == 1.OR.edflg == 4.OR.edflg == 8.OR.edflg == 10.OR.  &
    edflg == 11.OR.edflg < 0) GO TO 300
CALL holdit
300     RETURN
104     CALL rite(' That type does not exist.')
CALL holdit
GO TO 300
112     CALL rite(' READ ERROR in TYPES data file.')
CALL holdit
return
END SUBROUTINE vewunt
!=======================================================================

!       SUBROUTINE VEWBLK

!-----  This subroutine allows the user to view a block if it
!       was created as a block alone, i.e. with an extension
!       BLK

!=======================================================================

SUBROUTINE vewblk(crtflg,edflg,rtn_vewblk)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

!------ If EDFLG is set to 11 by the view block module then
!       the information is read into temporary storage for use
!       by this subroutine.

edflg=11
CALL readin(crtflg,edflg,rtn_readin)
IF(rtn_readin == 2) THEN
  rtn_vewblk = 2
  RETURN
END IF
100     lincnt=1
CALL rite(' ')
CALL rite(' BLOCK')
CALL rite(' ')

DO  i=1,nunct
  icheck=typtem(i)
  CALL typein(2)
  WRITE(*,1)numtem(i),typtem(i),dscrib
  IF (svmode) WRITE(luview,1)numtem(i),typtem(i),dscrib
  1         FORMAT(5X,'UNIT',i3,5X,'TYPE',i3,' - ',a50)
END DO

!------ The information read in from the read block subroutine is
!       assigned to the proper arrays for use by the view a unit
!       subroutine. The information was stored in the temporary arrays
!       found in the common block INSCOM.

DO  i=1,nunct
  untnum(i)=numtem(i)
  untyp(i)=typtem(i)
  blkin(untnum(i))=tin(i)
  blkout(untnum(i))=tout(i)
  blkpar(untnum(i))=tpar(i)
  DO  j=1,blkin(untnum(i))
    indxin(untnum(i),j)=xintem(i,j)
  END DO
  DO  j=1,blkout(untnum(i))
    indxot(untnum(i),j)=xottem(i,j)
  END DO
  DO  j=1,blkpar(untnum(i))
    pvalue(untnum(i),j)=valtem(i,j)
  END DO
END DO

CALL vewunt(crtflg,edflg,rtn_vewunt)
IF(rtn_vewunt == 1) THEN
  rtn_vewblk = 1
  RETURN
else if(rtn_vewunt == 2) then
  rtn_vewblk = 2
  return
END IF
RETURN
END SUBROUTINE vewblk
!=======================================================================

!     SUBROUTINE VEWSIM

!---- This subroutine provides the menu that contains the different
!     view options.

!=======================================================================

SUBROUTINE vewsim(crtflg,edflg,rtn_vewsim)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=2) :: answer

IF (edflg /= 4.AND.crtflg /= 4) THEN
  
!------ EDFLG set to 8 indicating a call from VIEWSIM
  
  edflg=8
  CALL readin(crtflg,edflg,rtn_readin)
  IF(rtn_readin == 2) THEN
    rtn_vewsim = 2
    RETURN
  END IF
END IF
120  WRITE(*,40201)
  40201    FORMAT('  What part of the simulation would you like to view:'  &
      ,//,' ALl the simulation information (for documentation)'  &
      ,//,' STructure (superblock,block, and unit Information)'  &
      ,//,' VAriable initial values'  &
      ,//,' ERror tolerances, variable scan and freeze options'  &
      ,//,' BOundary variables' ,//,' REported variables'  &
      ,//,' COntinue with the previous menu')
  messag=' '
  CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
  CALL revew(item,crtflg,edflg,rtn_revew)
  IF(rtn_revew == 1) THEN
    GO TO 120
  ELSE IF(rtn_revew == 2) THEN
    rtn_vewsim = 2
    RETURN
  END IF
  answer=item
IF (answer == 'AL'.OR.answer == 'al') THEN
  CALL vewall
ELSE IF (answer == 'ST'.OR.answer == 'st') THEN
  CALL struct(crtflg,edflg,rtn_struct)
  if(rtn_struct == 1) then
     go to 120
  elseif(rtn_struct == 2) then
     rtn_vewsim = 2
     return
  endif

ELSE IF (answer == 'VA'.OR.answer == 'va') THEN
  CALL varval(crtflg,edflg,rtn_varval)
  if(rtn_varval == 1) then
     go to 120
  elseif(rtn_varval == 2) then
     rtn_vewsim = 2
     return
  endif

ELSE IF (answer == 'ER'.OR.answer == 'er') THEN
  CALL errors(crtflg,edflg,rtn_errors)
  if(rtn_errors == 1) then
     go to 120
  elseif(rtn_errors == 2) then
     rtn_vewsim = 2
     return
  endif

ELSE IF (answer == 'BO'.OR.answer == 'bo') THEN
  CALL bound(crtflg,edflg,rtn_bound)
  if(rtn_bound == 1) then
     go to 120
  elseif(rtn_bound == 2) then
     rtn_vewsim = 2
     return
  endif

ELSE IF (answer == 'RE'.OR.answer == 're') THEN
  CALL rptvar(crtflg,edflg,rtn_rptvar)
  if(rtn_rptvar == 1) then
     go to 120
  elseif(rtn_rptvar == 2) then
     rtn_vewsim = 2
     return
  endif

ELSE IF (answer == 'CO'.OR.answer == 'co') THEN
  rtn_vewsim = 1
  RETURN
END IF

IF (intnum == 0) then
  rtn_vewsim = 1
  RETURN
endif
GO TO 120
END SUBROUTINE vewsim
!=======================================================================

!      SUBROUTINE STRUCT

!----- This subroutine allows viewing of the structure of a simulation

!=======================================================================

SUBROUTINE struct(crtflg,edflg,rtn_struct)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

INTEGER :: maxsct,lcnt,mcnt

l=1
lcnt=1
mcnt=1
maxsct=simcnt(1)
IF (crtflg == 4) decsup=simsup-1
WRITE(*,40301) titlsm
IF (svmode) WRITE(luview,40301) titlsm
40301   FORMAT(/,1X,a72)
lincnt=1

DO  nn=1,decsup
  WRITE(*,40302) nn
  IF (svmode) WRITE(luview,40302) nn
  40302     FORMAT(/,' SUPERBLOCK',i2)
  IF (edflg >= 0) CALL scroll(lincnt)
  DO  m=lcnt,maxsct
    WRITE(*,40303) m
    IF (svmode) WRITE(luview,40303) m
    40303       FORMAT(5X,' BLOCK',i2)
    IF (edflg >= 0) CALL scroll(lincnt)
    DO  n=1,blkunt(m)
      icheck=untyp(untnum(l))
      CALL typein(2)
      WRITE(*,4) untnum(l),untyp(untnum(l)),dscrib
      IF (svmode) WRITE(luview,4) untnum(l),untyp(untnum(l)),dscrib
      4             FORMAT(10X,'UNIT',i3,5X,'TYPE',i3,' - ',a45)
      IF (edflg >= 0) CALL scroll(lincnt)
      l=l+1
    END DO
    mcnt=mcnt+1
  END DO
  lcnt=mcnt
  maxsct=simcnt(nn+1)+(lcnt-1)
END DO

IF (edflg == 4.OR.edflg < 0) GO TO 300
IF (edflg == 8) CALL vewunt(crtflg,edflg,rtn_vewunt)
IF(rtn_vewunt == 2) THEN
  rtn_struct = 2
  RETURN
END IF
300       IF (edflg < 0)edflg = l-1
rtn_struct = 1
RETURN
END SUBROUTINE struct
!=======================================================================

!     SUBROUTINE VARVAL

!------ This subroutine allows the user to view the initial variable
!       values for all the inputs.
!=======================================================================

SUBROUTINE varval(crtflg,edflg,rtn_varval)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

WRITE(*,1)
IF (svmode) WRITE(luview,1)
1       FORMAT(' Initial Variable Values:',/)
lincnt=1

IF (pres <= 0) GO TO 150
DO  i=1,pres
  WRITE(*,2) i,state(i),iounit(1)
  IF (svmode) WRITE(luview,2) i,state(i),iounit(1)
  2         FORMAT(' PRESSURE       ',i3,' -> ',g15.6,a20)
  IF (edflg > 0) CALL scroll(lincnt)
END DO

150       IF (flow <= 0) GO TO 250
DO  j=1,flow
  WRITE(*,3) j,state(j+pres),iounit(2)
  IF (svmode) WRITE(luview,3) j,state(j+pres),iounit(2)
  3         FORMAT(' FLOW           ',i3,' -> ',g15.6,a20)
  IF (edflg > 0) CALL scroll(lincnt)
END DO

250       IF (temp <= 0) GO TO 350
DO  k=1,temp
  WRITE(*,4) k,state(k+pres+flow),iounit(3)
  IF (svmode) WRITE(luview,4) k,state(k+pres+flow),iounit(3)
  4         FORMAT(' TEMPERATURE    ',i3,' -> ',g15.6,a20)
  IF (edflg > 0) CALL scroll(lincnt)
END DO

350       IF (cntr <= 0) GO TO 451
DO  l=1,cntr
  WRITE(*,5) l,state(l+pres+flow+temp),iounit(4)
  IF (svmode) WRITE(luview,5) l,state(l+pres+flow+temp),iounit(4)
  5         FORMAT(' CONTROL        ',i3,' -> ',g15.6,a20)
  IF (edflg > 0) CALL scroll(lincnt)
END DO

451       IF (othr <= 0) GO TO 550
DO  m=1,othr
  WRITE(*,6) m,state(m+pres+flow+temp+cntr),iounit(5)
  IF (svmode) WRITE(luview,6) m,state(m+pres+flow+temp+cntr),iounit(5)
  6         FORMAT(' OTHER          ',i3,' -> ',g15.6,a20)
  IF (edflg > 0) CALL scroll(lincnt)
END DO

550       IF (enrg <= 0) GO TO 650
DO  n=1,enrg
  WRITE(*,7) n,state(n+pres+flow+temp+cntr+othr),iounit(6)
  IF (svmode ) WRITE(luview,7)  &
      n,state(n+pres+flow+temp+cntr+othr),iounit(6)
  7         FORMAT(' ENERGY         ',i3,' -> ',g15.6,a20)
  IF (edflg > 0) CALL scroll(lincnt)
END DO

650     IF (powr <= 0) GO TO 750
DO  ii=1,powr
  WRITE(*,8) ii,state(ii+pres+flow+temp+cntr+othr+enrg), iounit(7)
  IF (svmode) WRITE(luview,8) ii,state(ii+pres+flow+temp+cntr+othr+enrg),  &
      iounit(7)
  8         FORMAT(' POWER          ',i3,' -> ',g15.6,a20)
  IF (edflg > 0) CALL scroll(lincnt)
END DO

750     IF (ahum <= 0) GO TO 850
DO  jj=1,ahum
  WRITE(*,9) jj,state(jj+pres+flow+temp+cntr+othr+enrg+powr), iounit(8)
  IF (svmode) WRITE(luview,9)  &
      jj,state(jj+pres+flow+temp+cntr+othr+enrg+powr), iounit(8)
  9         FORMAT(' HUMIDITY       ',i3,' -> ',g15.6,a20)
  IF (edflg > 0) CALL scroll(lincnt)
END DO

850     IF (edflg == 4.OR.edflg < 0) GO TO 300
900     messag=' Hit the Carriage Return to continue'
CALL datain(item,anumbr,1.,0.,intnum,line,report,2,3,0,messag)
CALL revew(item,crtflg,edflg,rtn_revew)
IF(rtn_revew == 1) THEN
  GO TO 900
ELSE IF(rtn_revew == 2) THEN
  rtn_varval = 2
  RETURN
END IF

IF (intnum == 0) GO TO 300
GO TO 900
300     rtn_varval = 1
RETURN
END SUBROUTINE varval
!=======================================================================

!     SUBROUTINE BOUND

!------- This subroutine is used to view the boundary information for a
!        simulation.

!=======================================================================

SUBROUTINE bound(crtflg,edflg,rtn_bound)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

WRITE(*,1)
IF (svmode) WRITE(luview,1)
1       FORMAT(' The following are Boundary Variables in the ',  &
    'simulation:',/)
lincnt=1
IF (edflg >= 0) CALL scroll(lincnt)

DO  k=1,bcnt
  
  IF (bndry(k) > pres.OR.pres <= 0) GO TO 210
  WRITE(*,2) bndry(k)
  IF (svmode) WRITE(luview,2) bndry(k)
  2       FORMAT(' PRESSURE      ',i3)
  IF (edflg >= 0) CALL scroll(lincnt)
  CYCLE
  
  210     IF (bndry(k) > (pres+flow).OR.flow <= 0) GO TO 220
  WRITE(*,3) (bndry(k)-pres)
  IF (svmode) WRITE(luview,3) (bndry(k)-pres)
  3       FORMAT(' FLOW          ',i3)
  IF (edflg >= 0) CALL scroll(lincnt)
  CYCLE
  
  220     IF (bndry(k) > (pres+flow+temp).OR.temp <= 0) GO TO 230
  WRITE(*,4) (bndry(k)-pres-flow)
  IF (svmode) WRITE(luview,4) (bndry(k)-pres-flow)
  4       FORMAT(' TEMPERATURE   ',i3)
  IF (edflg >= 0) CALL scroll(lincnt)
  CYCLE
  
  230     IF (bndry(k) > (pres+flow+temp+cntr).OR.cntr <= 0) GO TO 240
  WRITE(*,5) (bndry(k)-pres-flow-temp)
  IF (svmode) WRITE(luview,5) (bndry(k)-pres-flow-temp)
  5       FORMAT(' CONTROL       ',i3)
  IF (edflg >= 0) CALL scroll(lincnt)
  CYCLE
  
  240     IF (bndry(k) > (pres+flow+temp+cntr+othr).OR.othr <= 0) GO TO 250
  WRITE(*,6) (bndry(k)-pres-flow-temp-cntr)
  IF (svmode) WRITE(luview,6) (bndry(k)-pres-flow-temp-cntr)
  6       FORMAT(' OTHER          ',i3)
  IF (edflg >= 0) CALL scroll(lincnt)
  CYCLE
  
  250     IF (bndry(k) > (pres+flow+temp+cntr+othr+enrg).OR.enrg <= 0)  &
      GO TO 260
  WRITE(*,7) (bndry(k)-pres-flow-temp-cntr-othr)
  IF (svmode) WRITE(luview,7) (bndry(k)-pres-flow-temp-cntr-othr)
  7       FORMAT(' ENERGY        ',i3)
  IF (edflg >= 0) CALL scroll(lincnt)
  CYCLE
  
  260     IF (bndry(k) > (pres+flow+temp+cntr+othr+enrg+powr).OR.powr  &
       <= 0) GO TO 270
  WRITE(*,8) (bndry(k)-pres-flow-temp-cntr-othr-enrg)
  IF (svmode) WRITE(luview,8) (bndry(k)-pres-flow-temp-cntr-othr-enrg)
  8       FORMAT(' POWER         ',i3)
  IF (edflg >= 0) CALL scroll(lincnt)
  CYCLE
  
  270     IF (bndry(k) > (pres+flow+temp+cntr+othr+enrg+powr+ahum)  &
      .OR.ahum <= 0) CYCLE
  WRITE(*,9) (bndry(k)-pres-flow-temp-cntr-othr-enrg-powr)
  IF (svmode) WRITE(luview,9)  &
      (bndry(k)-pres-flow-temp-cntr-othr-enrg-powr)
  9       FORMAT(' HUMIDITY      ',i3)
  IF (edflg >= 0) CALL scroll(lincnt)
END DO

IF (edflg < 0) GO TO 300
295     messag=' Hit the Carriage Return to continue'
CALL datain(item,anumbr,1.,0.,intnum,line,report,2,3,0,messag)
CALL revew(item,crtflg,edflg,rtn_revew)
IF(rtn_revew == 1) THEN
  GO TO 295
ELSE IF(rtn_revew == 2) THEN
  rtn_bound = 2
  RETURN
END IF
IF (intnum == 0) GO TO 300
GO TO 295
300     rtn_bound = 1
RETURN
END SUBROUTINE bound
!=======================================================================

!     SUBROUTINE RPTVAR

!------- This subroutine is used to view the reported variables in a
!        simulation.

!=======================================================================

SUBROUTINE rptvar(crtflg,edflg,rtn_rptvar)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

lincnt=1

DO  k=1,decsup
  
  WRITE(*,1)k,repinv(k)
  IF (svmode) WRITE(luview,1)k,repinv(k)
  1       FORMAT(/,' SUPERBLOCK',i2,5X,' REPORTING INTERVAL',g15.6)
  IF (edflg >= 0) CALL scroll(lincnt)
  
  DO  m=1,scnt(k)
    
    IF (repvar(k,m) > pres.OR.pres <= 0) GO TO 210
    WRITE(*,2) repvar(k,m)
    IF (svmode) WRITE(luview,2) repvar(k,m)
    2       FORMAT(' PRESSURE       ',i3)
    IF (edflg >= 0) CALL scroll(lincnt)
    CYCLE
    
    210     IF (repvar(k,m) > (pres+flow).OR.flow <= 0) GO TO 220
    WRITE(*,3) (repvar(k,m)-pres)
    IF (svmode) WRITE(luview,3) (repvar(k,m)-pres)
    3       FORMAT(' FLOW           ',i3)
    IF (edflg >= 0) CALL scroll(lincnt)
    CYCLE
    
    220     IF (repvar(k,m) > (pres+flow+temp).OR.temp <= 0) GO TO 230
    WRITE(*,4) (repvar(k,m)-pres-flow)
    IF (svmode) WRITE(luview,4) (repvar(k,m)-pres-flow)
    4       FORMAT(' TEMPERATURE    ',i3)
    IF (edflg >= 0) CALL scroll(lincnt)
    CYCLE
    
    230     IF (repvar(k,m) > (pres+flow+temp+cntr).OR.cntr <= 0) GO TO 240
    WRITE(*,5) (repvar(k,m)-pres-flow-temp)
    IF (svmode) WRITE(luview,5) (repvar(k,m)-pres-flow-temp)
    5       FORMAT(' CONTROL        ',i3)
    IF (edflg >= 0) CALL scroll(lincnt)
    CYCLE
    
    240     IF (repvar(k,m) > (pres+flow+temp+cntr+othr).OR.othr <= 0)  &
        GO TO 250
    WRITE(*,6) (repvar(k,m)-pres-flow-temp-cntr)
    IF (svmode) WRITE(luview,6) (repvar(k,m)-pres-flow-temp-cntr)
    6       FORMAT(' OTHER          ',i3)
    IF (edflg >= 0) CALL scroll(lincnt)
    CYCLE
    
    250     IF (repvar(k,m) > (pres+flow+temp+cntr+othr+enrg).OR.enrg <= 0)  &
        GO TO 260
    WRITE(*,7) (repvar(k,m)-pres-flow-temp-cntr-othr)
    IF (svmode) WRITE(luview,7) (repvar(k,m)-pres-flow-temp-cntr-othr)
    7       FORMAT(' ENERGY         ',i3)
    IF (edflg >= 0) CALL scroll(lincnt)
    CYCLE
    
    260     IF (repvar(k,m) > (pres+flow+temp+cntr+othr+enrg+powr).OR.powr  &
         <= 0) GO TO 270
    WRITE(*,8) (repvar(k,m)-pres-flow-temp-cntr-othr-enrg)
    IF (svmode) WRITE(luview,8) (repvar(k,m)-pres-flow-temp-cntr-othr-enrg)
    8       FORMAT(' POWER          ',i3)
    IF (edflg >= 0) CALL scroll(lincnt)
    CYCLE
    
    270     IF (repvar(k,m) > (pres+flow+temp+cntr+othr+enrg+powr+ahum)  &
        .OR.ahum <= 0) CYCLE
    WRITE(*,9) (repvar(k,m)-pres-flow-temp-cntr-othr-enrg-powr)
    IF (svmode) WRITE(luview,9)  &
        (repvar(k,m)-pres-flow-temp-cntr-othr-enrg-powr)
    9       FORMAT(' HUMIDITY       ',i3)
    IF (edflg >= 0) CALL scroll(lincnt)
    
  END DO
END DO

IF (edflg < 0) GO TO 300
295     messag=' Hit the Carriage Return to continue'
CALL datain(item,anumbr,1.,0.,intnum,line,report,2,3,0,messag)
CALL revew(item,crtflg,edflg,rtn_revew)
IF(rtn_revew == 1) THEN
  GO TO 295
ELSE IF(rtn_revew == 2) THEN
  rtn_rptvar = 2
  RETURN
END IF
IF (intnum == 0) GO TO 300
GO TO 295
300     rtn_rptvar = 1
RETURN
END SUBROUTINE rptvar
!=======================================================================

!     SUBROUTINE ERRORS

!------- This subroutine allows the user to view the error tolerances
!        and the freeze and scan options.

!=======================================================================

SUBROUTINE errors(crtflg,edflg,rtn_errors)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer              :: lnumbr

lincnt=1
lnumbr=1
WRITE(*,1)
IF (svmode) WRITE(luview,1)
1       FORMAT(' Simulation Error Tolerances:',/)
WRITE(*,2) lnumbr,rtolx,atolx
IF (svmode) WRITE(luview,2) lnumbr,rtolx,atolx
2       FORMAT(i2,5X,' RTOLX=',g15.6,5X,' ATOLX=',g15.6)
IF (edflg >= 0) CALL scroll(lincnt)
WRITE(*,3) xtol,ttime
IF (svmode) WRITE(luview,3) xtol,ttime
3       FORMAT(7X,' XTOL=',1X,g15.6,5X,' TTIME=',g15.6)
IF (edflg >= 0) CALL scroll(lincnt)

!------- The freeze and scan options are displayed for each superblock.

DO  k=1,decsup
  lnumbr=lnumbr+1
  WRITE(*,4)k
  IF (svmode) WRITE(luview,4)k
  4         FORMAT(/,' SUPERBLOCK',i2)
  IF (edflg >= 0) CALL scroll(lincnt)
  WRITE(*,5)lnumbr,freeze(k),scan(k)
  IF (svmode) WRITE(luview,5)lnumbr,freeze(k),scan(k)
  5         FORMAT(i2,5X,' FREEZE OPTION',i2,5X,' SCAN OPTION',i2)
  IF (edflg >= 0) CALL scroll(lincnt)
END DO

IF (edflg == 4.OR.edflg < 0) GO TO 300
200     messag=' Hit the Carriage Return to continue'
CALL datain(item,anumbr,1.,0.,intnum,line,report,2,3,0,messag)
CALL revew(item,crtflg,edflg,rtn_revew)
IF(rtn_revew == 1) THEN
  GO TO 200
ELSE IF(rtn_revew == 2) THEN
  rtn_errors = 2
  RETURN
END IF
IF (intnum == 0) GO TO 300
GO TO 200
300     rtn_errors = 1
RETURN
END SUBROUTINE errors
!=======================================================================

!       SUBROUTINE REVEW

!------ This subroutine is called by the VIEW module and handles
!       all calls to other selected modules. This subroutine will
!       only allow calls to the  HELP modules from the VIEW module.
!------ POINTR is described in REWORD.

!=======================================================================

SUBROUTINE revew(item,crtflg,edflg,rtn_revew)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=1) :: sure
CHARACTER (LEN=2) :: answer

answer=item
IF (answer == 'AB'.OR.answer == 'ab') GO TO 250
IF (answer == 'HE'.OR.answer == 'he') GO TO 101
RETURN

!------ A call is made to HELP with program control transferred there.

101     CALL help(crtflg,edflg,abort)
        rtn_revew = 1
        RETURN
250     messag=' Are you sure you wish to ABORT, Y/N'
CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
sure=item
IF (sure == 'Y'.OR.sure == 'y') then
  rtn_revew = 2
else
  rtn_revew = 1
endif

RETURN
END SUBROUTINE revew
!=======================================================================

!     SUBROUTINE VEWALL

!--- This routine is used to document a simulation by displaying all at
!  once a complete description of a simulation. There is no pause after
!  each screenful. The purpose of this is to allow the output to be
!  captured by a printer rather than being read as it passes by on the
!  screen.

!=======================================================================

SUBROUTINE vewall

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer              :: io,nunits
CHARACTER (LEN=1)    :: answer
character (len=40)   :: viewsave_file

PRINT *
PRINT *, ' Save the model setup? (y/n) '
READ (*,FMT='(A1)')  answer
IF (answer == 'y' .or. answer == 'Y' .or. answer == ' ' ) THEN
  svmode = .true.
ELSE
  svmode = .false.
END IF

if (svmode) then
  print *, 'use the default file name (viewsave.txt) (y/n)?'
  read (*,fmt='(a1)')  answer
  if (answer == 'y' .or. answer == 'Y'  .or. answer == ' ' ) then
    open (unit = luview, file = 'viewsave.txt',status='new', iostat=io)
    if (io /= 0 ) then
      close(luview)
      print *,'**** abort:  viewsave.txt already exists ****** '
      return
    end if
  else
    print *, 'Enter a full file name including an extension'
    read (*,fmt='(a40)') viewsave_file
    open (unit = luview, file = viewsave_file,status='new', iostat=io)
    if (io /= 0 ) then
      close(luview)
      print *, '**** abort:  ', viewsave_file,' already exists ****** '
      return
    end if
  endif
end if

edflg = -1
crtflg = 0
CALL struct(crtflg,edflg,rtn_struct)
IF(rtn_struct == 2) THEN
  CALL holdit
  RETURN
END IF
100   nunits = edflg

DO  n = 1,nunits
  edflg = -n
  CALL rite (' --------------------------------------------------')
  CALL vewunt(crtflg,edflg,rtn_vewunt)
  IF(rtn_vewunt == 2) THEN
    CALL holdit
    RETURN
  END IF
END DO

edflg = -1
CALL rite (' ----------------------------------------------------')
CALL varval(crtflg,edflg,rtn_varval)
IF(rtn_varval == 2) THEN
  CALL holdit
  RETURN
END IF
300   CALL rite (' ----------------------------------------------------')
CALL errors(crtflg,edflg,rtn_errors)
IF(rtn_errors == 2) THEN
  CALL holdit
  RETURN
END IF
400   CALL rite (' ----------------------------------------------------')
CALL bound(crtflg,edflg,rtn_bound)
IF(rtn_bound == 2) THEN
  CALL holdit
  RETURN
END IF
500   CALL rite (' ----------------------------------------------------')
CALL rite(' The following are the Reported Variables:')
CALL rite(' ')
CALL rptvar(crtflg,edflg,rtn_rptvar)
IF(rtn_rptvar == 2) THEN
  CALL holdit
  RETURN
END IF
900   CALL holdit
RETURN
END SUBROUTINE vewall

