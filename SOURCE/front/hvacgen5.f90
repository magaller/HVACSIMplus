!=======================================================================
! HVACGEN SOURCE FILE #5 OF 6, VERSION 5.0
!=======================================================================

! hvacgen5.f90

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

!       EDIT MODULE

!------ This module transfers control to the proper subroutines
!       as indicated by the user's selection of editting a SImulation,
!       BLock, UNit.
!------ CHOICE is the variable which indicates what the user
!       would like to editting, i.e., a SImulation, BLock, Unit.
!=======================================================================

SUBROUTINE edit(crtflg,edflg,abort)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=2) :: choice

abort = .false.
100     WRITE(*,50101)
50101   FORMAT(1X,'Edit a:',//,' SImulation',/,' UNit')
!50101   FORMAT(1X,'Edit a:',//, ' SImulation',/,&
!      & ' SUperblock',/, ' BLock',/,' UNit')
messag=' Enter Selection'
CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  edflg=0
  abort = .true.
  RETURN
END IF
choice=item

110     IF (choice == 'SI'.OR.choice == 'si') THEN
  IF (crtflg /= 1.AND.crtflg /= 2.AND.crtflg /= 3.AND.crtflg  &
       /= 4) CALL edsim(crtflg,edflg,abort)
!ELSE IF (choice == 'SU'.OR.choice == 'su') THEN
!  RETURN
!ELSE IF (choice == 'BL'.OR.choice == 'bl') THEN
!  RETURN
ELSE IF (choice == 'UN'.OR.choice == 'un') THEN
  CALL edunt(crtflg,edflg,abort)
END IF

edflg=0
RETURN
400     edflg=0
abort = .true.
RETURN
END SUBROUTINE edit
!=======================================================================

!      SUBROUTINE EDUNT

!----- This subroutine allows the user to edit a unit.

!----- EDFLG indicates a call from the EDUNT module, used by other
!      modules. This value set to 1 in EDUNT.

!----- INSWER used to determine what line the user would like to edit.
!=======================================================================

SUBROUTINE edunt(crtflg,edflg,abort)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=15),dimension(8) :: varcat
CHARACTER (LEN=15)              :: iprmpt,oprmpt

DATA varcat /' PRESSURE',' FLOW',' TEMPERATURE',' CONTROL',  &
    ' OTHER',' ENERGY',' POWER',' HUMIDITY'/

abort = .false.
change = .false.
rtn_vewunt = 0                      ! 4/5/07
IF (edflg == 4) GO TO 100
edflg=1
IF (crtflg == 1.OR.crtflg == 2.OR.crtflg == 3.OR.crtflg == 4) GO TO 100
CALL readin(crtflg,edflg,rtn_readin)
IF(rtn_readin == 2) THEN
  abort = .true.
  RETURN
END IF
100     CALL vewunt(crtflg,edflg,rtn_vewunt)
! write(*,*)' rtn_vewunt = ', rtn_vewunt   !  debug
IF(rtn_vewunt == 1) THEN
  GO TO 290
ELSE IF(rtn_vewunt == 2) THEN
  abort = .true.
  RETURN
END IF

!------ The line that is to be edited is entered.

150     messag=' Enter the line number to edit or hit CR to end editing'
CALL datain(item,anumbr,3.,0.,inswer,line,report,2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 150
ELSE IF(rtn_reedt == 2) THEN
  abort = .true.
  RETURN
END IF

!------ The input indeces are changed if the user desires

IF (inswer == 1) THEN
  WRITE(*,1)
  1         FORMAT(/,' INPUTS',/)
  
  DO  j=1,blkin(untnum(nunit))
    iprmpt = varcat(iident(j))
    223         WRITE(*,4) iprmpt,indxin(untnum(nunit),j),imssge(j)
    4           FORMAT(7X,a15,1X,i3,' - ',a50)
    CALL okay(item,anumbr,REAL(maxunt*minoiu),0.,intnum,same, abort)
    IF (abort) RETURN
    IF (.NOT.same) THEN
      indxin(untnum(nunit),j) = intnum
      change = .true.
      GO TO 223
    END IF
  END DO
  
!------ The output indeces are changed if the user desires.
  
ELSE IF (inswer == 2) THEN
  WRITE(*,3)
  3         FORMAT(/,' OUTPUTS',/)
  DO  i=1,blkout(untnum(nunit))
    oprmpt = varcat(iodent(i))
    233         WRITE(*,4) oprmpt,indxot(untnum(nunit),i),omssge(i)
    CALL okay(item,anumbr,REAL(maxunt*minoiu),0.,intnum,same, abort)
    IF (abort) RETURN
    IF (.NOT.same) THEN
      indxot(untnum(nunit),i) = intnum
      change = .true.
      GO TO 233
    END IF
  END DO
  
!------ The parameter values are changed if desired.
  
ELSE IF (inswer == 3) THEN
  WRITE(*,5)
  5         FORMAT(/,' PARAMETERS',/)
  
  DO  k=1,blkpar(untnum(nunit))
    243         WRITE(*,7) pvalue(untnum(nunit),k),pmssge(k)
    7           FORMAT(7X,g15.6,1X,a55)
    CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
    IF (abort) RETURN
    IF (.NOT.same) THEN
      pvalue(untnum(nunit),k) = anumbr
      change = .true.
      GO TO 243
    END IF
  END DO
  
ELSE IF (inswer == 0) THEN
  IF (crtflg == 1.OR.crtflg == 2.OR.crtflg == 3.OR.crtflg == 4) RETURN
  IF (edflg == 1) GO TO 290
  GO TO 100
END IF

GO TO 150
290    IF (edflg == 4) THEN
  CALL recalc(rtn_recalc)
  IF(rtn_recalc == 1) THEN
    abort = .true.
    RETURN
  END IF
END IF
IF (change) CALL fsave(crtflg,edflg,rtn_fsave)
if(rtn_fsave == 2) then
    abort = .true.
endif
RETURN
END SUBROUTINE edunt
!=======================================================================

!     SUBROUTINE EDSIM

!---- This subroutine provides the menu that contains the different
!     simulation editting options.

!=======================================================================

SUBROUTINE edsim(crtflg,edflg,abort)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=2) :: answer

abort = .false.
reask = .false.
edflg=4
CALL readin(crtflg,edflg,rtn_readin)
IF(rtn_readin == 2) THEN
  abort = .true.
  RETURN
END IF
120   messag=' '
CALL datain(item,anumbr,0.,0.,intnum,line,report,1,0,1,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)

IF(rtn_reedt == 1) THEN
   GO TO 120
ELSE IF(rtn_reedt == 2) THEN
  abort = .true.
  RETURN
END IF
  CALL rite(' What portion of the simulation would you like to edit:')
  CALL rite(' ')
  CALL rite(' TItle (change the title of the simulation)')
  CALL rite(' ')
  CALL rite(' STructure (superblock,block,and unit information)')
  CALL rite(' ')
  CALL rite(' VAriable initial values')
  CALL rite(' ')
  CALL rite(' BOundary variables')
  CALL rite(' ')
  CALL rite(' REported variables')
  CALL rite(' ')
  CALL rite(' ERror tolerances, variable freeze and scan option')
  CALL rite(' ')
  CALL rite(' DElete (block, or unit)')
  CALL rite(' ')
  CALL rite(' INsert or replace (block, or unit)')
  CALL rite(' ')
  CALL rite(' SAve in work file  (block)')
  CALL rite(' ')
  CALL rite(' COntinue with the previous menu')
  CALL rite(' ')
  messag=' '
  CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
  CALL reedt(item,crtflg,edflg,rtn_reedt)
  IF(rtn_reedt == 1) THEN
    GO TO 120
  ELSE IF(rtn_reedt == 2) THEN
    abort = .true.
    RETURN
  END IF
  answer=item

IF (answer == 'TI'.OR.answer == 'ti') THEN
  CALL edtitl(crtflg,edflg,rtn_edtitl)
  IF(rtn_edtitl == 1) THEN
    GO TO 120
  ELSE IF(rtn_edtitl == 2) THEN
    abort = .true.
    RETURN
  END IF
ELSE IF (answer == 'ST'.OR.answer == 'st') THEN
  CALL edstr(crtflg,edflg,rtn_edstr)
  IF(rtn_edstr == 1) THEN
    GO TO 120
  ELSE IF(rtn_edstr == 2) THEN
    abort = .true.
    RETURN
  END IF
ELSE IF (answer == 'VA'.OR.answer == 'va') THEN
  CALL edval(crtflg,edflg,rtn_edval)
  IF(rtn_edval == 1) THEN
    GO TO 120
  ELSE IF(rtn_edval == 2) THEN
    abort = .true.
    RETURN
  END IF
ELSE IF (answer == 'ER'.OR.answer == 'er') THEN
  CALL ederr(crtflg,edflg,rtn_ederr)
  IF(rtn_ederr == 1) THEN
    GO TO 120
  ELSE IF(rtn_ederr == 2) THEN
    abort = .true.
    RETURN
  END IF
ELSE IF (answer == 'BO'.OR.answer == 'bo') THEN
  CALL edbnd(crtflg,edflg,rtn_edbnd)
  IF(rtn_edbnd == 1) THEN
    GO TO 120
  ELSE IF(rtn_edbnd == 2) THEN
    abort = .true.
    RETURN
  END IF
ELSE IF (answer == 'RE'.OR.answer == 're') THEN
  CALL edrep(crtflg,edflg,rtn_edrep)
  IF(rtn_edrep == 1) THEN
    GO TO 120
  ELSE IF(rtn_edrep == 2) THEN
    abort = .true.
    RETURN
  END IF
ELSE IF (answer == 'DE'.OR.answer == 'de') THEN
  CALL delsim(crtflg,edflg,rtn_delsim)
  IF(rtn_delsim == 1) THEN
    GO TO 120
  ELSE IF(rtn_delsim == 2) THEN
    abort = .true.
    RETURN
  END IF
ELSE IF (answer == 'IN'.OR.answer == 'in') THEN
  CALL inssim(crtflg,edflg,rtn_inssim)
  IF(rtn_inssim == 1) THEN
    GO TO 120
  ELSE IF(rtn_inssim == 2) THEN
    abort = .true.
    RETURN
  END IF
ELSE IF (answer == 'SA'.OR.answer == 'sa') THEN
  CALL extblk(crtflg,edflg,rtn_extblk)
  IF(rtn_extblk == 1) THEN
    GO TO 120
  ELSE IF(rtn_extblk == 2) THEN
    abort = .true.
    RETURN
  END IF
END IF

IF (reask) GO TO 120
IF (abort) RETURN
IF (answer == 'CO'.OR.answer == 'co') GO TO 300
300   RETURN
END SUBROUTINE edsim
!=======================================================================

!     SUBROUTINE EDTITL

!----- This subroutine allows the user to edit the tile of a
!      simulation.
!=======================================================================

SUBROUTINE edtitl(crtflg,edflg,rtn_edtitl)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=1) :: answer

change = .false.
90     CALL rite(' The current title for the simulation is:')
CALL rite(' ')
WRITE(*,1)titlsm
1       FORMAT(1X,a75)
CALL rite(' ')
100     messag=' Is this title correct, Y/N ?'
iok = -13
CALL datain(item,anumbr,0.,0.,iok,line,report,2,0,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_edtitl = 2
  RETURN
END IF
answer=item

IF (iok == 0.AND.(answer == 'N'.OR.answer == 'n')) THEN
  200       messag=' Enter the new title for the simulation'
  CALL datain(item,anumbr,0.,0.,intnum,line,report,3,2,1,messag)
  CALL reedt(item,crtflg,edflg,rtn_reedt)
  IF(rtn_reedt == 1) THEN
    GO TO 200
  ELSE IF(rtn_reedt == 2) THEN
    rtn_edtitl = 2
    RETURN
  END IF
  titlsm=line
  change = .true.
  GO TO 90
ELSE
  IF (change) CALL fsave(crtflg,edflg,rtn_fsave)
    if(rtn_fsave == 1) then
      rtn_edtitl = 1
      return
    elseif(rtn_fsave == 2) then
      rtn_edtitl = 2
      return
    endif
END IF

RETURN
END SUBROUTINE edtitl
!=======================================================================

!      SUBROUTINE EDSTR

!----- This subroutine allows the user to edit the simulation
!      structure.
!=======================================================================

SUBROUTINE edstr(crtflg,edflg,rtn_edstr)

implicit none
include 'hvacgen.inc'

CALL struct(crtflg,edflg,rtn_struct)
IF(rtn_struct == 2) THEN
  rtn_edstr = 2
  RETURN
END IF
100     CALL edunt(crtflg,edflg,abort)
IF (abort) THEN
  rtn_edstr = 2
  RETURN
END IF
RETURN
END SUBROUTINE edstr
!=======================================================================

!     SUBROUTINE EDVAL

!----- This allows the user to edit the initial variable values in the
!      state vector.

!=======================================================================

SUBROUTINE edval(crtflg,edflg, rtn_edval)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

change = .false.
lincnt=1
100     WRITE(*,50001)
50001   FORMAT(' Enter the line number of the variable category to ',  &
    'edit:' ,/,' 1 PRESSURE           5 OTHER'  &
    ,/,' 2 FLOW               6 ENERGY' ,/,' 3 TEMPERATURE        7 POWER'  &
    ,/,' 4 CONTROL            8 HUMIDITY'  &
    ,/,' or hit the Carriage Return to return to previous menu.')
messag=' '
CALL datain(item,anumbr,8.,0.,inswer,line,report,2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_edval = 2
  RETURN
END IF

IF (inswer == 1) THEN
  IF (pres <= 0) GO TO 500
  DO  i=1,pres
    203         WRITE(*,2) i,state(i),iounit(1)
    2           FORMAT(' PRESSURE      ',i3,' -> ',g15.6,a20)
    CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
    IF (abort) THEN
      rtn_edval = 2
      RETURN
    END IF
    IF (.NOT.same) THEN
      state(i)=anumbr
      change = .true.
      GO TO 203
    END IF
  END DO
  
ELSE IF (inswer == 2) THEN
  IF (flow <= 0) GO TO 500
  DO  i=1,flow
    213        WRITE(*,3) i,state(i+pres),iounit(2)
    3           FORMAT(' FLOW          ',i3,' -> ',g15.6,a20)
    CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
    IF (abort) THEN
      rtn_edval = 2
      RETURN
    END IF
    IF (.NOT.same) THEN
      state(i+pres)=anumbr
      change = .true.
      GO TO 213
    END IF
  END DO
  
ELSE IF (inswer == 3) THEN
  IF (temp <= 0) GO TO 500
  DO  i=1,temp
    223         WRITE(*,4) i,state(i+pres+flow),iounit(3)
    4           FORMAT(' TEMPERATURE   ',i3,' -> ',g15.6,a20)
    CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
    IF (abort) THEN
      rtn_edval = 2
      RETURN
    END IF
    IF (.NOT.same) THEN
      state(i+pres+flow)=anumbr
      change = .true.
      GO TO 223
    END IF
  END DO
  
ELSE IF (inswer == 4) THEN
  IF (cntr <= 0) GO TO 500
  DO  i=1,cntr
    233         WRITE(*,5) i,state(i+pres+flow+temp),iounit(4)
    5           FORMAT(' CONTROL       ',i3,' -> ',g15.6,a20)
    CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
    IF (abort) THEN
      rtn_edval = 2
      RETURN
    END IF
    IF (.NOT.same) THEN
      state(i+pres+flow+temp)=anumbr
      change = .true.
      GO TO 233
    END IF
  END DO
  
ELSE IF (inswer == 5) THEN
  IF (othr <= 0) GO TO 500
  DO  i=1,othr
    243         WRITE(*,6) i,state(i+pres+flow+temp+cntr),iounit(5)
    6           FORMAT(' OTHER          ',i3,' -> ',g15.6,a20)
    CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
    IF (abort) THEN
      rtn_edval = 2
      RETURN
    END IF
    IF (.NOT.same) THEN
      state(i+pres+flow+temp+cntr)=anumbr
      change = .true.
      GO TO 243
    END IF
  END DO
  
ELSE IF (inswer == 6) THEN
  IF (enrg <= 0) GO TO 500
  DO  i=1,enrg
    253         WRITE(*,7) i,state(i+pres+flow+temp+cntr+othr),iounit(6)
    7           FORMAT(' ENERGY       ',i3,' -> ',g15.6,a20)
    CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
    IF (abort) THEN
      rtn_edval = 2
      RETURN
    END IF
    IF (.NOT.same) THEN
      state(i+pres+flow+temp+cntr+othr)=anumbr
      change = .true.
      GO TO 253
    END IF
  END DO
  
ELSE IF (inswer == 7) THEN
  IF (powr <= 0) GO TO 500
  DO  i=1,powr
    263         WRITE(*,8) i,state(i+pres+flow+temp+cntr+othr+enrg), iounit(7)
    8           FORMAT(' POWER         ',i3,' -> ',g15.6,a20)
    CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
    IF (abort) THEN
      rtn_edval = 2
      RETURN
    END IF
    IF (.NOT.same) THEN
      state(i+pres+flow+temp+cntr+othr+enrg)=anumbr
      change = .true.
      GO TO 263
    END IF
  END DO
  
ELSE IF (inswer == 8) THEN
  IF (ahum <= 0) GO TO 500
  DO  i=1,ahum
    273         WRITE(*,9) i,state(i+pres+flow+temp+cntr+othr+enrg+powr),  &
        iounit(8)
    9           FORMAT(' HUMIDITY      ',i3,' -> ',g15.6,a20)
    CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
    IF (abort) THEN
      rtn_edval = 2
      RETURN
    END IF
    IF (.NOT.same) THEN
      state(i+pres+flow+temp+cntr+othr+enrg+powr)=anumbr
      change = .true.
      GO TO 273
    END IF
  END DO
  
ELSE IF (inswer == 0) THEN
  IF (change) CALL fsave(crtflg,edflg,rtn_fsave)
  IF(rtn_fsave == 1) THEN
    rtn_edval = 1
    RETURN
  elseif(rtn_fsave == 2) then
    rtn_edval = 2
    return
  END IF

ELSE
  GO TO 100
END IF

GO TO 100
500     CALL rite(' There are no variables in that category.')
CALL rite(' ')
GO TO 100
END SUBROUTINE edval
!=======================================================================

!       SUBROUTINE EDBND

!------ This subroutine allows the user to change the boundary variables
!       that were set up during the create process.

!=======================================================================

SUBROUTINE edbnd(crtflg,edflg,rtn_edbnd)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

! Display boundary variables (using view routine) and get choice------

change = .false.
100     CALL bound(0,0,rtn_bound)
IF(rtn_bound == 2) THEN
  rtn_edbnd = 2
  RETURN
END IF
101     WRITE(*,50501)
50501   FORMAT(  &
    ' Enter the line number of the boundary variable edit option',  &
    /,' 1 INSERT a boundary variable', /,' 2 DELETE a boundary variable',  &
    /,' or hit the Carriage Return to return to the previous menu.')
messag=' '
CALL datain(item,anumbr,3.,0.,inswer,line,report,2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_edbnd = 2
  RETURN
END IF
IF (inswer == 1) THEN
  change = .true.
  CALL insert(crtflg,edflg,rtn_insert)
  IF(rtn_insert == 1) THEN
    GO TO 100                        
  ELSE IF(rtn_insert == 2) THEN
    rtn_edbnd = 2                    
    RETURN                           
  END IF                             
ELSE IF (inswer == 2) THEN
  change = .true.
  CALL DELETE(crtflg,edflg,rtn_delete)
  IF(rtn_delete == 1) THEN
    GO TO 100                          
  ELSE IF(rtn_delete == 2) THEN
    rtn_edbnd = 2                      
    RETURN                             
  END IF                               
END IF
IF (change) CALL fsave(crtflg,edflg,rtn_fsave)
IF(rtn_fsave == 1) THEN
  rtn_edbnd = 1
  RETURN
elseif(rtn_fsave == 2) then
  rtn_edbnd = 2
  return
END IF
RETURN
END SUBROUTINE edbnd
!=======================================================================

!     SUBROUTINE INSERT

!------ This program allows the user to insert a boundary variable.
!------ This subroutine calls a module OKAY that asks the user for
!       the new input.

!=======================================================================

SUBROUTINE insert(crtflg,edflg,rtn_insert)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer              :: ibcnt

100     CALL rite(' Enter the line number of the variable that you   ')
CALL rite(' wish to insert:                                  ')
CALL rite(' ')
CALL rite(' 1 PRESSURE           5 OTHER                     ')
CALL rite(' 2 FLOW               6 ENERGY                    ')
CALL rite(' 3 TEMPERATURE        7 POWER                     ')
CALL rite(' 4 CONTROL            8 HUMIDITY                  ')
CALL rite(' or hit the Carriage Return to return to previous ')
CALL rite(' menu.                                            ')

messag=' '
CALL datain(item,anumbr,8.,0.,inswer,line,report,2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_insert = 2
  RETURN
END IF
ibcnt=bcnt
SELECT CASE ( inswer )
  CASE (    1)
    GO TO 200
  CASE (    2)
    GO TO 210
  CASE (    3)
    GO TO 220
  CASE (    4)
    GO TO 230
  CASE (    5)
    GO TO 240
  CASE (    6)
    GO TO 250
  CASE (    7)
    GO TO 260
  CASE (    8)
    GO TO 270
END SELECT
IF (inswer == 0)RETURN
GO TO 100

!------- Pressure boundary variables

200     IF (pres <= 0) GO TO 500
201     messag=' Enter the index of a new PRESSURE boundary variable'
CALL datain(item,anumbr,REAL(pres),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 201
ELSE IF(rtn_reedt == 2) THEN
  rtn_insert = 2
  RETURN
END IF
202     WRITE(*,2) intnum
2       FORMAT(' PRESSURE       ',i3)
CALL okay(item,anumbr,REAL(pres),1.,intnum,same,abort)
IF (abort) THEN
  rtn_insert = 2
  RETURN
END IF
IF (.NOT.same) GO TO 202
IF (bcnt == 0) THEN
  bndry(1) = intnum
  GO TO 290
ELSE
  DO  k=1,ibcnt
    IF (bndry(k) <= pres) THEN
      IF (k == ibcnt) THEN
        bndry(k+1)=intnum
        GO TO 290
      END IF
    ELSE
      DO  j=ibcnt+1,k+1,-1
        bndry(j)=bndry(j-1)
      END DO
      bndry(k)=intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_insert = 1
RETURN

!------- Flow boundary variables

210     IF (flow <= 0) GO TO 500
211     messag=' Enter the index of a new FLOW boundary variable'
CALL datain(item,anumbr,REAL(flow),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 211
ELSE IF(rtn_reedt == 2) THEN
  rtn_insert = 2
  RETURN
END IF
212     WRITE(*,3) intnum
3       FORMAT(' FLOW           ',i3)
CALL okay(item,anumbr,REAL(flow),1.,intnum,same,abort)
IF (abort) THEN
  rtn_insert = 2
  RETURN
END IF
IF (.NOT.same) GO TO 212
IF (bcnt == 0) THEN
  bndry(1) = intnum+pres
  GO TO 290
ELSE
  DO  k=1,ibcnt
    IF (bndry(k) <= pres+flow) THEN
      IF (k == ibcnt) THEN
        bndry(k+1)=pres+intnum
        GO TO 290
      END IF
    ELSE
      DO  j=ibcnt+1,k+1,-1
        bndry(j)=bndry(j-1)
      END DO
      bndry(k)=pres+intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_insert = 1
RETURN

!------- Temperature boundary variables

220     IF (temp <= 0) GO TO 500
221     messag=' Enter the index of a new TEMPERATURE boundary variable'
CALL datain(item,anumbr,REAL(temp),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 221
ELSE IF(rtn_reedt == 2) THEN
  rtn_insert = 2
  RETURN
END IF
222     WRITE(*,4) intnum
4       FORMAT(' TEMPERATURE    ',i3)
CALL okay(item,anumbr,REAL(temp),1.,intnum,same,abort)
IF (abort) THEN
  rtn_insert = 2
  RETURN
END IF
IF (.NOT.same) GO TO 222
IF (bcnt == 0) THEN
  bndry(1) = intnum+pres+flow
  GO TO 290
ELSE
  DO  k=1,ibcnt
    IF (bndry(k) <= pres+flow+temp) THEN
      IF (k == ibcnt) THEN
        bndry(k+1)=pres+flow+intnum
        GO TO 290
      END IF
    ELSE
      DO  j=ibcnt+1,k+1,-1
        bndry(j)=bndry(j-1)
      END DO
      bndry(k)=pres+flow+intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_insert = 1
RETURN

!------- Control boundary variables

230     IF (cntr <= 0) GO TO 500
231     messag=' Enter the index of a new CONTROL boundary variable'
CALL datain(item,anumbr,REAL(cntr),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 231
ELSE IF(rtn_reedt == 2) THEN
  rtn_insert = 2
  RETURN
END IF
232     WRITE(*,5) intnum
5       FORMAT(' CONTROL        ',i3)
CALL okay(item,anumbr,REAL(cntr),1.,intnum,same,abort)
IF (abort) THEN
  rtn_insert = 2
  RETURN
END IF
IF (.NOT.same) GO TO 232
IF (bcnt == 0) THEN
  bndry(1) = intnum+pres+flow+temp
  GO TO 290
ELSE
  DO  k=1,ibcnt
    IF (bndry(k) <= pres+flow+temp+cntr) THEN
      IF (k == ibcnt) THEN
        bndry(k+1)=pres+flow+temp+intnum
        GO TO 290
      END IF
    ELSE
      DO  j=ibcnt+1,k+1,-1
        bndry(j)=bndry(j-1)
      END DO
      bndry(k)=pres+flow+temp+intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_insert = 1
RETURN

!------- Other boundary variables

240     IF (othr <= 0) GO TO 500
241     messag=' Enter the index of a new OTHER boundary variable'
CALL datain(item,anumbr,REAL(othr),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 241
ELSE IF(rtn_reedt == 2) THEN
  rtn_insert = 2
  RETURN
END IF
242     WRITE(*,6) intnum
6       FORMAT(' OTHER          ',i3)
CALL okay(item,anumbr,REAL(othr),1.,intnum,same,abort)
IF (abort) THEN
  rtn_insert = 2
  RETURN
END IF
IF (.NOT.same) GO TO 242
IF (bcnt == 0) THEN
  bndry(1) = intnum+pres+flow+temp+cntr
  GO TO 290
ELSE
  DO  k=1,ibcnt
    IF (bndry(k) <= pres+flow+temp+cntr+othr) THEN
      IF (k == ibcnt) THEN
        bndry(k+1)=pres+flow+temp+cntr+intnum
        GO TO 290
      END IF
    ELSE
      DO  j=ibcnt+1,k+1,-1
        bndry(j)=bndry(j-1)
      END DO
      bndry(k)=pres+flow+temp+cntr+intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_insert = 1
RETURN

!------- Energy boundary variables

250     IF (enrg <= 0) GO TO 500
251     messag=' Enter the index of a new ENERGY boundary variable'
CALL datain(item,anumbr,REAL(enrg),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 251
ELSE IF(rtn_reedt == 2) THEN
  rtn_insert = 2
  RETURN
END IF
252     WRITE(*,7) intnum
7       FORMAT(' ENERGY         ',i3)
CALL okay(item,anumbr,REAL(enrg),1.,intnum,same,abort)
IF (abort) THEN
  rtn_insert = 2
  RETURN
END IF
IF (.NOT.same) GO TO 252
IF (bcnt == 0) THEN
  bndry(1) = intnum+pres+flow+temp+cntr+othr
  GO TO 290
ELSE
  DO  k=1,ibcnt
    IF (bndry(k) <= pres+flow+temp+cntr+othr+enrg) THEN
      IF (k == ibcnt) THEN
        bndry(k+1)=pres+flow+temp+cntr+othr+intnum
        GO TO 290
      END IF
    ELSE
      DO  j=ibcnt+1,k+1,-1
        bndry(j)=bndry(j-1)
      END DO
      bndry(k)=pres+flow+temp+cntr+othr+intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_insert = 1
RETURN

!------- Power boundary variables

260     IF (powr <= 0) GO TO 500
261     messag=' Enter the index of a new POWER boundary variable'
CALL datain(item,anumbr,REAL(powr),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 261
ELSE IF(rtn_reedt == 2) THEN
  rtn_insert = 2
  RETURN
END IF
262     WRITE(*,8) intnum
8       FORMAT(' POWER          ',i3)
CALL okay(item,anumbr,REAL(powr),1.,intnum,same,abort)
IF (abort) THEN
  rtn_insert = 2
  RETURN
END IF
IF (.NOT.same) GO TO 262
IF (bcnt == 0) THEN
  bndry(1) = intnum+pres+flow+temp+cntr+othr+enrg
  GO TO 290
ELSE
  DO  k=1,ibcnt
    IF (bndry(k) <= pres+flow+temp+cntr+othr+enrg+powr) THEN
      IF (k == ibcnt) THEN
        bndry(k+1)=pres+flow+temp+cntr+othr+enrg+intnum
        GO TO 290
      END IF
    ELSE
      DO  j=ibcnt+1,k+1,-1
        bndry(j)=bndry(j-1)
      END DO
      bndry(k)=pres+flow+temp+cntr+othr+enrg+intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_insert = 1
RETURN

!-------  Humidity boundary variables

270     IF (ahum <= 0) GO TO 500
271     messag=' Enter the index of a new HUMIDITY boundary variable'
CALL datain(item,anumbr,REAL(ahum),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 271
ELSE IF(rtn_reedt == 2) THEN
  rtn_insert = 2
  RETURN
END IF
272     WRITE(*,9) intnum
9       FORMAT(' HUMIDITY       ',i3)
CALL okay(item,anumbr,REAL(ahum),1.,intnum,same,abort)
IF (abort) THEN
  rtn_insert = 2
  RETURN
END IF
IF (.NOT.same) GO TO 272
IF (bcnt == 0) THEN
  bndry(1) = intnum+pres+flow+temp+cntr+othr+enrg+powr
  GO TO 290
ELSE
  DO  k=1,ibcnt
    IF (bndry(k) <= pres+flow+temp+cntr+othr+enrg+powr+ahum) THEN
      IF (k == ibcnt) THEN
        bndry(k+1)=pres+flow+temp+cntr+othr+enrg+powr+intnum
        GO TO 290
      END IF
    ELSE
      DO  j=ibcnt+1,k+1,-1
        bndry(j)=bndry(j-1)
      END DO
      bndry(k)=pres+flow+temp+cntr+othr+enrg+powr+intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_insert = 1
RETURN
290     bcnt=bcnt+1
rtn_insert = 1
RETURN
500     CALL rite(' There are no variables in that category.')
CALL rite(' ')
GO TO 100
END SUBROUTINE insert
!=======================================================================

!     SUBROUTINE DELETE

!------ This program allows the user to delete a boundary variable.
!------ This subroutine calls a module OKAY that asks the user for
!       the new input.

!=======================================================================

SUBROUTINE DELETE(crtflg,edflg,rtn_delete)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'
integer              :: ibcnt

100     CALL rite(' Enter the line number of the category of         ')
CALL rite(' variable from which to delete:                   ')
CALL rite(' ')
CALL rite(' 1 PRESSURE        5 OTHER                        ')
CALL rite(' 2 FLOW            6 ENERGY                       ')
CALL rite(' 3 TEMPERATURE     7 POWER                        ')
CALL rite(' 4 CONTROL         8 HUMIDITY                     ')
CALL rite(' or hit the Carriage Return to return to previous ')
CALL rite(' menu.                                            ')

messag=' '
CALL datain(item,anumbr,8.,0.,inswer,line,report,2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_delete = 2
  RETURN
END IF
ibcnt=bcnt
SELECT CASE ( inswer )
  CASE (    1)
    GO TO 200
  CASE (    2)
    GO TO 210
  CASE (    3)
    GO TO 220
  CASE (    4)
    GO TO 230
  CASE (    5)
    GO TO 240
  CASE (    6)
    GO TO 250
  CASE (    7)
    GO TO 260
  CASE (    8)
    GO TO 270
END SELECT
IF (inswer == 0) RETURN
GO TO 100

!------- Pressure boundary variables

200     IF (pres <= 0) GO TO 500
201     messag= ' Enter the index of the PRESSURE boundary variable to delete'
CALL datain(item,anumbr,REAL(pres),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 201
ELSE IF(rtn_reedt == 2) THEN
  rtn_delete = 2
  RETURN
END IF
202     WRITE(*,2) intnum
2       FORMAT(' PRESSURE       ',i3)
CALL okay(item,anumbr,REAL(pres),0.,intnum,same,abort)
IF (abort) THEN
  rtn_delete = 2
  RETURN
END IF
IF (.NOT.same) GO TO 202

DO  k=1,ibcnt
  IF (bndry(k) == intnum) THEN
    DO  j=k+1,ibcnt
      bndry(j-1)=bndry(j)
    END DO
    GO TO 290
  END IF
END DO
rtn_delete = 1
RETURN

!------- Flow boundary variables

210     IF (flow <= 0) GO TO 500
211     messag= ' Enter the index of the FLOW boundary variable to delete'
CALL datain(item,anumbr,REAL(flow),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 211
ELSE IF(rtn_reedt == 2) THEN
  rtn_delete = 2
  RETURN
END IF
212     WRITE(*,3) intnum
3       FORMAT(' FLOW           ',i3)
CALL okay(item,anumbr,REAL(flow),0.,intnum,same,abort)
IF (abort) THEN
  rtn_delete = 2
  RETURN
END IF
IF (.NOT.same) GO TO 212

DO  k=1,ibcnt
  IF (bndry(k) == pres+intnum) THEN
    DO  j=k+1,ibcnt
      bndry(j-1)=bndry(j)
    END DO
    GO TO 290
  END IF
END DO
rtn_delete = 1
RETURN

!------- Temperature boundary variables

220     IF (temp <= 0) GO TO 500
221     messag=  &
    ' Enter the index of the TEMPERATURE boundary variable to delete'
CALL datain(item,anumbr,REAL(temp),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 221
ELSE IF(rtn_reedt == 2) THEN
  rtn_delete = 2
  RETURN
END IF
222     WRITE(*,4) intnum
4       FORMAT(' TEMPERATURE    ',i3)
CALL okay(item,anumbr,REAL(temp),0.,intnum,same,abort)
IF (abort) THEN
  rtn_delete = 2
  RETURN
END IF
IF (.NOT.same) GO TO 222

DO  k=1,ibcnt
  IF (bndry(k) == pres+flow+intnum) THEN
    DO  j=k+1,ibcnt
      bndry(j-1)=bndry(j)
    END DO
    GO TO 290
  END IF
END DO
rtn_delete = 1
RETURN

!------- Control boundary variables

230     IF (cntr <= 0) GO TO 500
231     messag= ' Enter the index of the CONTROL boundary variable to delete'
CALL datain(item,anumbr,REAL(cntr),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 231
ELSE IF(rtn_reedt == 2) THEN
  rtn_delete = 2
  RETURN
END IF
232     WRITE(*,5) intnum
5       FORMAT(' CONTROL        ',i3)
CALL okay(item,anumbr,REAL(cntr),0.,intnum,same,abort)
IF (abort) THEN
  rtn_delete = 2
  RETURN
END IF
IF (.NOT.same) GO TO 232

DO  k=1,ibcnt
  IF (bndry(k) == pres+flow+temp+intnum) THEN
    DO  j=k+1,ibcnt
      bndry(j-1)=bndry(j)
    END DO
    GO TO 290
  END IF
END DO
rtn_delete = 1
RETURN

!------- Other boundary variables

240     IF (othr <= 0) GO TO 500
241     messag= ' Enter the index of the OTHER boundary variable to delete'
CALL datain(item,anumbr,REAL(othr),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 241
ELSE IF(rtn_reedt == 2) THEN
  rtn_delete = 2
  RETURN
END IF
242     WRITE(*,6) intnum
6       FORMAT(' OTHER          ',i3)
CALL okay(item,anumbr,REAL(othr),0.,intnum,same,abort)
IF (abort) THEN
  rtn_delete = 2
  RETURN
END IF
IF (.NOT.same) GO TO 242

DO  k=1,ibcnt
  IF (bndry(k) == pres+flow+temp+cntr+intnum) THEN
    DO  j=k+1,ibcnt
      bndry(j-1)=bndry(j)
    END DO
    GO TO 290
  END IF
END DO
rtn_delete = 1
RETURN

!------- Energy boundary variables

250     IF (enrg <= 0) GO TO 500
251     messag= ' Enter the index of the ENERGY boundary variable to delete'
CALL datain(item,anumbr,REAL(enrg),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 251
ELSE IF(rtn_reedt == 2) THEN
  rtn_delete = 2
  RETURN
END IF
252     WRITE(*,7) intnum
7       FORMAT(' ENERGY         ',i3)
CALL okay(item,anumbr,REAL(enrg),0.,intnum,same,abort)
IF (abort) THEN
  rtn_delete = 2
  RETURN
END IF
IF (.NOT.same) GO TO 252

DO  k=1,ibcnt
  IF (bndry(k) == pres+flow+temp+cntr+othr+intnum) THEN
    DO  j=k+1,ibcnt
      bndry(j-1)=bndry(j)
    END DO
    GO TO 290
  END IF
END DO
rtn_delete = 1
RETURN

!------- Power boundary variables

260     IF (powr <= 0) GO TO 500
261     messag= ' Enter the index of the POWER boundary variable to delete'
CALL datain(item,anumbr,REAL(powr),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 261
ELSE IF(rtn_reedt == 2) THEN
  rtn_delete = 2
  RETURN
END IF
262     WRITE(*,8) intnum
8       FORMAT(' POWER          ',i3)
CALL okay(item,anumbr,REAL(powr),0.,intnum,same,abort)
IF (abort) THEN
  rtn_delete = 2
  RETURN
END IF
IF (.NOT.same) GO TO 262

DO  k=1,ibcnt
  IF (bndry(k) == pres+flow+temp+cntr+othr+enrg+intnum) THEN
    DO  j=k+1,ibcnt
      bndry(j-1)=bndry(j)
    END DO
    GO TO 290
  END IF
END DO
rtn_delete = 1
RETURN

!-------  Humidity boundary variables

270     IF (ahum <= 0) GO TO 500
271     messag= ' Enter the index of the HUMIDITY boundary variable to delete'
CALL datain(item,anumbr,REAL(ahum),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 271
ELSE IF(rtn_reedt == 2) THEN
  rtn_delete = 2
  RETURN
END IF
272     WRITE(*,9) intnum
9       FORMAT(' HUMIDITY       ',i3)
CALL okay(item,anumbr,REAL(ahum),0.,intnum,same,abort)
IF (abort) THEN
  rtn_delete = 2
  RETURN
END IF
IF (.NOT.same) GO TO 272

DO  k=1,ibcnt
  IF (bndry(k) == pres+flow+temp+cntr+othr+enrg+powr+intnum) THEN
    DO  j=k+1,ibcnt
      bndry(j-1)=bndry(j)
    END DO
    GO TO 290
  END IF
END DO
rtn_delete = 1
RETURN
290     bcnt=bcnt-1
rtn_delete = 1
RETURN
500     CALL rite(' There are no variables in that category.')
CALL rite(' ')
GO TO 100
END SUBROUTINE DELETE
!=======================================================================

!       SUBROUTINE EDREP

!------ This subroutine allows the user to change the reported variables
!       that were set up during the create process.

!=======================================================================

SUBROUTINE edrep(crtflg,edflg,rtn_edrep)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

change = .false.
100     CALL rptvar(0,0,rtn_rptvar)
IF(rtn_rptvar == 2) THEN
  rtn_edrep = 2
  RETURN
END IF
101     WRITE(*,50501)
50501   FORMAT(  &
    ' Enter the line number of the reported variable edit option',  &
    /,' 1 CHANGE reporting interval', /,' 2 INSERT a reported variable',  &
    /,' 3 DELETE a reported variable',  &
    /,' or hit the Carriage Return to return to the previous menu.')

messag=' '
CALL datain(item,anumbr,3.,0.,inswer,line,report,2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_edrep = 2
  RETURN
END IF
IF (inswer == 1) THEN
  change = .true.
  CALL rpchng(crtflg,edflg,rtn_rpchng)
  IF(rtn_rpchng == 1) THEN
    GO TO 100
  ELSE IF(rtn_rpchng == 2) THEN
    rtn_edrep = 2
    RETURN
  END IF
ELSE IF (inswer == 2) THEN
  change = .true.
  CALL rpinrt(crtflg,edflg,rtn_rpinrt)
  IF(rtn_rpinrt == 1) THEN
    GO TO 100
  ELSE IF(rtn_rpinrt == 2) THEN
    rtn_edrep = 2
    RETURN
  END IF
ELSE IF (inswer == 3) THEN
  change = .true.
  CALL rpdelt(crtflg,edflg,rtn_rpdelt)
  IF(rtn_rpdelt == 1) THEN
    GO TO 100
  ELSE IF(rtn_rpdelt == 2) THEN
    rtn_edrep = 2
    RETURN
  END IF
END IF
IF (change) CALL fsave(crtflg,edflg,rtn_fsave)
if(rtn_fsave == 1) then
  rtn_edrep = 1
  return
elseif(rtn_fsave == 2)then
  rtn_edrep = 2
  return
endif
RETURN
END SUBROUTINE edrep
!=======================================================================

!     SUBROUTINE RPCHNG

!------ This program allows the user to simply change the index of
!       the reported variable.
!------ This subroutine calls a module OKAY that asks the user for
!       the new input.

!=======================================================================

SUBROUTINE rpchng(crtflg,edflg,rtn_rpchng)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

DO  nsuper = 1, decsup
  500       WRITE(*,10) nsuper,repinv(nsuper)
  10        FORMAT(' Reporting interval for SUPERBLOCK',i2,' is',g15.6,  &
      ' sec.')
  CALL okay(item,anumbr,1.e35,0.,intnum,same,abort)
  IF (abort) THEN
    rtn_rpchng = 2
    RETURN
  END IF
  
  IF (.NOT.same) THEN
    repinv(nsuper)=anumbr
    GO TO 500
  END IF
END DO
rtn_rpchng = 1
RETURN
END SUBROUTINE rpchng
!=======================================================================

!     SUBROUTINE RPINRT

!------ This program allows the user to insert a reported variable.
!------ This subroutine calls a module OKAY that asks the user for
!       the new input.

!=======================================================================

SUBROUTINE rpinrt(crtflg,edflg,rtn_rpinrt)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

100     IF (decsup == 1) THEN
  nsuper = 1
ELSE
  CALL rite(' Enter the Superblock number that is to be edited ')
  messag=' or hit the Carriage Return to return to previous menu'
  CALL datain(item,anumbr,REAL(decsup),0.,nsuper,line,report, 2,3,0,messag)
  CALL reedt(item,crtflg,edflg,rtn_reedt)
  IF(rtn_reedt == 1) THEN
    GO TO 100
  ELSE IF(rtn_reedt == 2) THEN
    rtn_rpinrt = 2
    RETURN
  END IF
  IF (nsuper == 0) then
    rtn_rpinrt = 1
    RETURN
  endif
END IF

101     CALL rite(' Enter the line of the type of variable that you  ')
CALL rite(' wish to insert:                                  ')
CALL rite(' ')
CALL rite(' 1 PRESSURE        4 CONTROL       7 POWER        ')
CALL rite(' 2 FLOW            5 OTHER         8 HUMIDITY     ')
CALL rite(' 3 TEMPERATURE     6 ENERGY                       ')
CALL rite(' or hit the Carriage Return to return to the      ')
CALL rite(' previous menu                                    ')

messag=' '
CALL datain(item,anumbr,8.,0.,inswer,line,report,2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpinrt = 2
  RETURN
END IF
IF (inswer == 0) then
  rtn_rpinrt = 1
  RETURN
endif
iscnt=scnt(nsuper)
SELECT CASE ( inswer )
  CASE (    1)
    GO TO 200
  CASE (    2)
    GO TO 210
  CASE (    3)
    GO TO 220
  CASE (    4)
    GO TO 230
  CASE (    5)
    GO TO 240
  CASE (    6)
    GO TO 250
  CASE (    7)
    GO TO 260
  CASE (    8)
    GO TO 270
END SELECT
GO TO 101

!------- Pressure reported variables

200     IF (pres <= 0) GO TO 500
201     messag=' Enter the index of a new PRESSURE reported variable'
CALL datain(item,anumbr,REAL(pres),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 201
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpinrt = 2
  RETURN
END IF
202     WRITE(*,2) intnum
2       FORMAT(' PRESSURE       ',i3)
CALL okay(item,anumbr,REAL(pres),1.,intnum,same,abort)
IF (abort) THEN
  rtn_rpinrt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 202
IF (iscnt == 0) THEN
  repvar(nsuper,1) = intnum
  ntyp(nsuper,1) = 1
  nindx(nsuper,1) = intnum
  GO TO 290
ELSE
  DO  k=1,iscnt
    IF (repvar(nsuper,k) <= pres) THEN
      IF (k == iscnt) THEN
        repvar(nsuper,k+1)=intnum
        ntyp(nsuper,k+1)=1
        nindx(nsuper,k+1)=intnum
        GO TO 290
      END IF
    ELSE
      DO  j=iscnt+1,k+1,-1
        repvar(nsuper,j)=repvar(nsuper,j-1)
        ntyp(nsuper,j)=ntyp(nsuper,j-1)
        nindx(nsuper,j)=nindx(nsuper,j-1)
      END DO
      repvar(nsuper,k)=intnum
      ntyp(nsuper,k)=1
      nindx(nsuper,k)=intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_rpinrt = 1
RETURN

!------- Flow reported variables

210      IF (flow <= 0) GO TO 500
211      messag=' Enter the index of a new FLOW reported variable'
CALL datain(item,anumbr,REAL(flow),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 211
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpinrt = 2
  RETURN
END IF
212      WRITE(*,3) intnum
3        FORMAT(' FLOW           ',i3)
CALL okay(item,anumbr,REAL(flow),1.,intnum,same,abort)
IF (abort) THEN
  rtn_rpinrt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 212
IF (iscnt == 0) THEN
  repvar(nsuper,1) = intnum+pres
  ntyp(nsuper,1) = 2
  nindx(nsuper,1) = intnum
  GO TO 290
ELSE
  DO  k=1,iscnt
    IF (repvar(nsuper,k) <= pres+flow) THEN
      IF (k == iscnt) THEN
        repvar(nsuper,k+1)=pres+intnum
        ntyp(nsuper,k+1)=2
        nindx(nsuper,k+1)=intnum
        GO TO 290
      END IF
    ELSE
      DO  j=iscnt+1,k+1,-1
        repvar(nsuper,j)=repvar(nsuper,j-1)
        ntyp(nsuper,j)=ntyp(nsuper,j-1)
        nindx(nsuper,j)=nindx(nsuper,j-1)
      END DO
      repvar(nsuper,k)=pres+intnum
      ntyp(nsuper,k)=2
      nindx(nsuper,k)=intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_rpinrt = 1
RETURN

!------- Temperature reported variables

220      IF (temp <= 0) GO TO 500
221      messag= ' Enter the index of a new TEMPERATURE reported variable'
CALL datain(item,anumbr,REAL(temp),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 221
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpinrt = 2
  RETURN
END IF
222      WRITE(*,4) intnum
4        FORMAT(' TEMPERATURE    ',i3)
CALL okay(item,anumbr,REAL(temp),1.,intnum,same,abort)
IF (abort) THEN
  rtn_rpinrt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 222
IF (iscnt == 0) THEN
  repvar(nsuper,1) = intnum+pres+flow
  ntyp(nsuper,1) = 3
  nindx(nsuper,1) = intnum
  GO TO 290
ELSE
  DO  k=1,iscnt
    IF (repvar(nsuper,k) <= pres+flow+temp) THEN
      IF (k == iscnt) THEN
        repvar(nsuper,k+1)=pres+flow+intnum
        ntyp(nsuper,k+1)=3
        nindx(nsuper,k+1)=intnum
        GO TO 290
      END IF
    ELSE
      DO  j=iscnt+1,k+1,-1
        repvar(nsuper,j)=repvar(nsuper,j-1)
        ntyp(nsuper,j)=ntyp(nsuper,j-1)
        nindx(nsuper,j)=nindx(nsuper,j-1)
      END DO
      repvar(nsuper,k)=pres+flow+intnum
      ntyp(nsuper,k)=3
      nindx(nsuper,k)=intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_rpinrt = 1
RETURN

!------- Control reported variables

230      IF (cntr <= 0) GO TO 500
231      messag=' Enter the index of a new CONTROL reported variable'
CALL datain(item,anumbr,REAL(cntr),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 231
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpinrt = 2
  RETURN
END IF
232      WRITE(*,5) intnum
5        FORMAT(' CONTROL        ',i3)
CALL okay(item,anumbr,REAL(cntr),1.,intnum,same,abort)
IF (abort) THEN
  rtn_rpinrt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 232
IF (iscnt == 0) THEN
  repvar(nsuper,1) = intnum+pres+flow+temp
  ntyp(nsuper,1) = 4
  nindx(nsuper,1) = intnum
  GO TO 290
ELSE
  DO  k=1,iscnt
    IF (repvar(nsuper,k) <= pres+flow+temp+cntr) THEN
      IF (k == iscnt) THEN
        repvar(nsuper,k+1)=pres+flow+temp+intnum
        ntyp(nsuper,k+1)=4
        nindx(nsuper,k+1)=intnum
        GO TO 290
      END IF
    ELSE
      DO  j=iscnt+1,k+1,-1
        repvar(nsuper,j)=repvar(nsuper,j-1)
        ntyp(nsuper,j)=ntyp(nsuper,j-1)
        nindx(nsuper,j)=nindx(nsuper,j-1)
      END DO
      repvar(nsuper,k)=pres+flow+temp+intnum
      ntyp(nsuper,k)=4
      nindx(nsuper,k)=intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_rpinrt = 1
RETURN

!------- Other reported variables

240      IF (othr <= 0) GO TO 500
241      messag=' Enter the index of a new OTHER reported variable'
CALL datain(item,anumbr,REAL(othr),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 241
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpinrt = 2
  RETURN
END IF
242      WRITE(*,6) intnum
6        FORMAT(' OTHER          ',i3)
CALL okay(item,anumbr,REAL(othr),1.,intnum,same,abort)
IF (abort) THEN
  rtn_rpinrt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 242
IF (iscnt == 0) THEN
  repvar(nsuper,1) = intnum+pres+flow+temp+cntr
  ntyp(nsuper,1) = 5
  nindx(nsuper,1) = intnum
  GO TO 290
ELSE
  DO  k=1,iscnt
    IF (repvar(nsuper,k) <= pres+flow+temp+cntr+othr) THEN
      IF (k == iscnt) THEN
        repvar(nsuper,k+1)=pres+flow+temp+cntr+intnum
        ntyp(nsuper,k+1)=5
        nindx(nsuper,k+1)=intnum
        GO TO 290
      END IF
    ELSE
      DO  j=iscnt+1,k+1,-1
        repvar(nsuper,j)=repvar(nsuper,j-1)
        ntyp(nsuper,j)=ntyp(nsuper,j-1)
        nindx(nsuper,j)=nindx(nsuper,j-1)
      END DO
      repvar(nsuper,k)=pres+flow+temp+cntr+intnum
      ntyp(nsuper,k)=5
      nindx(nsuper,k)=intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_rpinrt = 1
RETURN

!------- Energy reported variables

250      IF (enrg <= 0) GO TO 500
251      messag=' Enter the index of a new ENERGY reported variable'
CALL datain(item,anumbr,REAL(enrg),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 251
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpinrt = 2
  RETURN
END IF
252      WRITE(*,7) intnum
7        FORMAT(' ENERGY         ',i3)
CALL okay(item,anumbr,REAL(enrg),1.,intnum,same,abort)
IF (abort) THEN
  rtn_rpinrt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 252
IF (iscnt == 0) THEN
  repvar(nsuper,1) = intnum+pres+flow+temp+cntr+othr
  ntyp(nsuper,1) = 6
  nindx(nsuper,1) = intnum
  GO TO 290
ELSE
  DO  k=1,iscnt
    IF (repvar(nsuper,k) <= pres+flow+temp+cntr+othr+enrg) THEN
      IF (k == iscnt) THEN
        repvar(nsuper,k+1)=pres+flow+temp+cntr+othr+intnum
        ntyp(nsuper,k+1)=6
        nindx(nsuper,k+1)=intnum
        GO TO 290
      END IF
    ELSE
      DO  j=iscnt+1,k+1,-1
        repvar(nsuper,j)=repvar(nsuper,j-1)
        ntyp(nsuper,j)=ntyp(nsuper,j-1)
        nindx(nsuper,j)=nindx(nsuper,j-1)
      END DO
      repvar(nsuper,k)=pres+flow+temp+cntr+othr+intnum
      ntyp(nsuper,k)=6
      nindx(nsuper,k)=intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_rpinrt = 1
RETURN

!------- Power reported variables

260      IF (powr <= 0) GO TO 500
261      messag=' Enter the index of a new POWER reported variable'
CALL datain(item,anumbr,REAL(powr),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 261
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpinrt = 2
  RETURN
END IF
262      WRITE(*,8) intnum
8        FORMAT(' POWER          ',i3)
CALL okay(item,anumbr,REAL(powr),1.,intnum,same,abort)
IF (abort) THEN
  rtn_rpinrt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 262
IF (iscnt == 0) THEN
  repvar(nsuper,1) = intnum+pres+flow+temp+cntr+othr+enrg
  ntyp(nsuper,1) = 7
  nindx(nsuper,1) = intnum
  GO TO 290
ELSE
  DO  k=1,iscnt
    IF (repvar(nsuper,k) <= pres+flow+temp+cntr+othr+enrg+ powr) THEN
      IF (k == iscnt) THEN
        repvar(nsuper,k+1)=pres+flow+temp+cntr+othr+enrg+intnum
        ntyp(nsuper,k+1)=7
        nindx(nsuper,k+1)=intnum
        GO TO 290
      END IF
    ELSE
      DO  j=iscnt+1,k+1,-1
        repvar(nsuper,j)=repvar(nsuper,j-1)
        ntyp(nsuper,j)=ntyp(nsuper,j-1)
        nindx(nsuper,j)=nindx(nsuper,j-1)
      END DO
      repvar(nsuper,k)=pres+flow+temp+cntr+othr+enrg+intnum
      ntyp(nsuper,k)=7
      nindx(nsuper,k)=intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_rpinrt = 1
RETURN

!-------  Humidity reported variables

270     IF (ahum <= 0) GO TO 500
271     messag=' Enter the index of a new HUMIDITY reported variable'
CALL datain(item,anumbr,REAL(ahum),1.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 271
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpinrt = 2
  RETURN
END IF
272     WRITE(*,9) intnum
9       FORMAT(' HUMIDITY       ',i3)
CALL okay(item,anumbr,REAL(ahum),1.,intnum,same,abort)
IF (abort) THEN
  rtn_rpinrt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 272
IF (iscnt == 0) THEN
  repvar(nsuper,1) = intnum+pres+flow+temp+cntr+othr+enrg+powr
  ntyp(nsuper,1) = 8
  nindx(nsuper,1) = intnum
  GO TO 290
ELSE
  DO  k=1,iscnt
    IF (repvar(nsuper,k) <= pres+flow+temp+cntr+othr+enrg+powr+ ahum) THEN
      IF (k == iscnt) THEN
        repvar(nsuper,k+1)=pres+flow+temp+cntr+othr+enrg+powr+ intnum
        ntyp(nsuper,k+1)=8
        nindx(nsuper,k+1)=intnum
        GO TO 290
      END IF
    ELSE
      DO  j=iscnt+1,k+1,-1
        repvar(nsuper,j)=repvar(nsuper,j-1)
        ntyp(nsuper,j)=ntyp(nsuper,j-1)
        nindx(nsuper,j)=nindx(nsuper,j-1)
      END DO
      repvar(nsuper,k)=pres+flow+temp+cntr+othr+enrg+powr+intnum
      ntyp(nsuper,k)=8
      nindx(nsuper,k)=intnum
      GO TO 290
    END IF
  END DO
END IF
rtn_rpinrt = 1
RETURN

290     scnt(nsuper)=scnt(nsuper)+1
scnt2(nsuper)=scnt(nsuper)
scnt3(nsuper)=scnt(nsuper)
rtn_rpinrt = 1
RETURN
500     CALL rite(' There are no variables in that category.')
CALL rite(' ')
GO TO 100
END SUBROUTINE rpinrt
!=======================================================================

!     SUBROUTINE RPDELT

!------ This program allows the user to delete a reported variable.
!------ This subroutine calls a module OKAY that asks the user for
!       the new input.

!=======================================================================

SUBROUTINE rpdelt(crtflg,edflg,rtn_rpdelt)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

100     IF (decsup == 1) THEN
  nsuper = 1
ELSE
  CALL rite(' Enter the Superblock number that is to be edited ')
  messag=' or hit the Carriage Return to return to previous menu'
  CALL datain(item,anumbr,REAL(decsup),0.,nsuper,line,report, 2,3,0,messag)
  CALL reedt(item,crtflg,edflg,rtn_reedt)
  IF(rtn_reedt == 1) THEN                       
    GO TO 100                                   
  ELSE IF(rtn_reedt == 2) THEN                  
    rtn_rpdelt = 2                              
    RETURN                                      
  END IF                                        
  IF (nsuper == 0) then
    rtn_rpdelt = 1
    RETURN
  endif
END IF

101     CALL rite(' Enter the line of the type of variable that you  ')
CALL rite(' wish to delete:                                  ')
CALL rite(' ')
CALL rite(' 1 PRESSURE        4 CONTROL   7 POWER            ')
CALL rite(' 2 FLOW            5 OTHER     8 HUMIDITY         ')
CALL rite(' 3 TEMPERATURE     6 ENERGY                       ')
CALL rite(' or hit the Carriage Return to return to the      ')
CALL rite(' previous menu                                    ')

messag=' '
CALL datain(item,anumbr,8.,0.,inswer,line,report,2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpdelt = 2
  RETURN
END IF
IF (inswer == 0) then
  rtn_rpdelt = 1
  RETURN
endif
iscnt=scnt(nsuper)
SELECT CASE ( inswer )
  CASE (    1)
    GO TO 200
  CASE (    2)
    GO TO 210
  CASE (    3)
    GO TO 220
  CASE (    4)
    GO TO 230
  CASE (    5)
    GO TO 240
  CASE (    6)
    GO TO 250
  CASE (    7)
    GO TO 260
  CASE (    8)
    GO TO 270
END SELECT
GO TO 101

!------- Pressure reported variables

200     IF (pres <= 0) GO TO 500
201     messag= ' Enter the index of the PRESSURE reported variable to delete'
CALL datain(item,anumbr,REAL(pres),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 201
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpdelt = 2
  RETURN
END IF
202     WRITE(*,2) intnum
2       FORMAT(' PRESSURE       ',i3)
CALL okay(item,anumbr,REAL(pres),0.,intnum,same,abort)
IF (abort) THEN
  rtn_rpdelt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 202

DO  k=1,iscnt
  IF (repvar(nsuper,k) == intnum) THEN
    DO  j=k+1,iscnt
      repvar(nsuper,j-1)=repvar(nsuper,j)
      ntyp(nsuper,j-1)=ntyp(nsuper,j)
      nindx(nsuper,j-1)=nindx(nsuper,j)
    END DO
    GO TO 290
  END IF
END DO
rtn_rpdelt = 1
RETURN

!------- Flow reported variables

210     IF (flow <= 0) GO TO 500
211     messag= ' Enter the index of the FLOW reported variable to delete'
CALL datain(item,anumbr,REAL(flow),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 211
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpdelt = 2
  RETURN
END IF
212     WRITE(*,3) intnum
3       FORMAT(' FLOW           ',i3)
CALL okay(item,anumbr,REAL(flow),0.,intnum,same,abort)
IF (abort) THEN
  rtn_rpdelt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 212

DO  k=1,iscnt
  IF (repvar(nsuper,k) == pres+intnum) THEN
    DO  j=k+1,iscnt
      repvar(nsuper,j-1)=repvar(nsuper,j)
      ntyp(nsuper,j-1)=ntyp(nsuper,j)
      nindx(nsuper,j-1)=nindx(nsuper,j)
    END DO
    GO TO 290
  END IF
END DO
rtn_rpdelt = 1
RETURN

!------- Temperature reported variables

220     IF (temp <= 0) GO TO 500
221     messag=  &
    ' Enter the index of the TEMPERATURE reported variable to delete'
CALL datain(item,anumbr,REAL(temp),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 221
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpdelt = 2
  RETURN
END IF
222     WRITE(*,4) intnum
4       FORMAT(' TEMPERATURE    ',i3)
CALL okay(item,anumbr,REAL(temp),0.,intnum,same,abort)
IF (abort) THEN
  rtn_rpdelt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 222

DO  k=1,iscnt
  IF (repvar(nsuper,k) == pres+flow+intnum) THEN
    DO  j=k+1,iscnt
      repvar(nsuper,j-1)=repvar(nsuper,j)
      ntyp(nsuper,j-1)=ntyp(nsuper,j)
      nindx(nsuper,j-1)=nindx(nsuper,j)
    END DO
    GO TO 290
  END IF
END DO
rtn_rpdelt = 1
RETURN

!------- Control reported variables

230     IF (cntr <= 0) GO TO 500
231     messag= ' Enter the index of the CONTROL reported variable to delete'
CALL datain(item,anumbr,REAL(cntr),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 231
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpdelt = 2
  RETURN
END IF
232     WRITE(*,5) intnum
5       FORMAT(' CONTROL        ',i3)
CALL okay(item,anumbr,REAL(cntr),0.,intnum,same,abort)
IF (abort) THEN
  rtn_rpdelt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 232

DO  k=1,iscnt
  IF (repvar(nsuper,k) == pres+flow+temp+intnum) THEN
    DO  j=k+1,iscnt
      repvar(nsuper,j-1)=repvar(nsuper,j)
      ntyp(nsuper,j-1)=ntyp(nsuper,j)
      nindx(nsuper,j-1)=nindx(nsuper,j)
    END DO
    GO TO 290
  END IF
END DO
rtn_rpdelt = 1
RETURN

!------- OTHR reported variables

240     IF (othr <= 0) GO TO 500
241     messag= ' Enter the index of the OTHER reported variable to delete'
CALL datain(item,anumbr,REAL(othr),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 241
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpdelt = 2
  RETURN
END IF
242     WRITE(*,6) intnum
6       FORMAT(' OTHER          ',i3)
CALL okay(item,anumbr,REAL(othr),0.,intnum,same,abort)
IF (abort) THEN
  rtn_rpdelt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 242

DO  k=1,iscnt
  IF (repvar(nsuper,k) == pres+flow+temp+cntr+intnum) THEN
    DO  j=k+1,iscnt
      repvar(nsuper,j-1)=repvar(nsuper,j)
      ntyp(nsuper,j-1)=ntyp(nsuper,j)
      nindx(nsuper,j-1)=nindx(nsuper,j)
    END DO
    GO TO 290
  END IF
END DO
rtn_rpdelt = 1
RETURN

!------- Energy reported variables

250     IF (enrg <= 0) GO TO 500
251     messag= ' Enter the index of the ENERGY reported variable to delete'
CALL datain(item,anumbr,REAL(enrg),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 251
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpdelt = 2
  RETURN
END IF
252     WRITE(*,7) intnum
7       FORMAT(' ENERGY         ',i3)
CALL okay(item,anumbr,REAL(enrg),0.,intnum,same,abort)
IF (abort) THEN
  rtn_rpdelt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 252

DO  k=1,iscnt
  IF (repvar(nsuper,k) == pres+flow+temp+cntr+othr+intnum) THEN
    DO  j=k+1,iscnt
      repvar(nsuper,j-1)=repvar(nsuper,j)
      ntyp(nsuper,j-1)=ntyp(nsuper,j)
      nindx(nsuper,j-1)=nindx(nsuper,j)
    END DO
    GO TO 290
  END IF
END DO
rtn_rpdelt = 1
RETURN

!------- Power reported variables

260     IF (powr <= 0) GO TO 500
261     messag= ' Enter the index of the POWER reported variable to delete'
CALL datain(item,anumbr,REAL(powr),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 261
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpdelt = 2
  RETURN
END IF
262     WRITE(*,8) intnum
8       FORMAT(' POWER          ',i3)
CALL okay(item,anumbr,REAL(powr),0.,intnum,same,abort)
IF (abort) THEN
  rtn_rpdelt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 262

DO  k=1,iscnt
  IF (repvar(nsuper,k) == pres+flow+temp+cntr+othr+enrg+ intnum) THEN
    DO  j=k+1,iscnt
      repvar(nsuper,j-1)=repvar(nsuper,j)
      ntyp(nsuper,j-1)=ntyp(nsuper,j)
      nindx(nsuper,j-1)=nindx(nsuper,j)
    END DO
    GO TO 290
  END IF
END DO
rtn_rpdelt = 1
RETURN

!-------  Humidity reported variables

270     IF (ahum <= 0) GO TO 500
271     messag= ' Enter the index of the HUMIDITY reported variable to delete'
CALL datain(item,anumbr,REAL(ahum),0.,intnum,line,report,2,3,0, messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 271
ELSE IF(rtn_reedt == 2) THEN
  rtn_rpdelt = 2
  RETURN
END IF
272     WRITE(*,9) intnum
9       FORMAT(' HUMIDITY       ',i3)
CALL okay(item,anumbr,REAL(ahum),0.,intnum,same,abort)
IF (abort) THEN
  rtn_rpdelt = 2
  RETURN
END IF
IF (.NOT.same) GO TO 272

DO  k=1,iscnt
  IF (repvar(nsuper,k) == pres+flow+temp+cntr+othr+enrg+powr+ intnum) THEN
    DO  j=k+1,iscnt
      repvar(nsuper,j-1)=repvar(nsuper,j)
      ntyp(nsuper,j-1)=ntyp(nsuper,j)
      nindx(nsuper,j-1)=nindx(nsuper,j)
    END DO
    GO TO 290
  END IF
END DO
rtn_rpdelt = 1
RETURN

290     scnt(nsuper)=scnt(nsuper)-1
scnt2(nsuper)=scnt(nsuper)
scnt3(nsuper)=scnt(nsuper)
rtn_rpdelt = 1
RETURN
450     rtn_rpdelt = 2
RETURN
500     CALL rite(' There are no variables in that category.')
CALL rite(' ')
GO TO 100
END SUBROUTINE rpdelt
!=======================================================================

!     SUBROUTINE EDERR

!------- This subroutine allows the user to edit the error tolerances
!        and the freeze and scan options.

!=======================================================================

SUBROUTINE ederr(crtflg,edflg,rtn_ederr)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

100     CALL rite(' Edit either:')
CALL rite(' ')
CALL rite(' 1 ERROR tolerances')
CALL rite(' 2 FREEZE and SCAN options')
CALL rite(' or hit the Carriage Return to return to the')
CALL rite(' previous menu')

messag=' '
CALL datain(item,anumbr,2.,0.,inswer,line,report,2,3,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF(rtn_reedt == 1) THEN
  GO TO 100
ELSE IF(rtn_reedt == 2) THEN
  rtn_ederr = 2
  RETURN
END IF
SELECT CASE ( inswer )
  CASE (    1)
    GO TO 101
  CASE (    2)
    GO TO 150
END SELECT
GO TO 300

101     WRITE(*,2) rtolx
2       FORMAT(' RTOLX=',g15.6)
CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
IF (abort) THEN
  rtn_ederr = 2
  RETURN
END IF
IF (.NOT.same) THEN
  rtolx=anumbr
  GO TO 101
END IF

102     WRITE(*,3) atolx
3       FORMAT(' ATOLX=',g15.6)
CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
IF (abort) THEN
  rtn_ederr = 2
  RETURN
END IF
IF (.NOT.same) THEN
  atolx=anumbr
  GO TO 102
END IF

103     WRITE(*,4) xtol
4       FORMAT(' XTOL=',g15.6)
CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
IF (abort) THEN
  rtn_ederr = 2
  RETURN
END IF
IF (.NOT.same) THEN
  xtol=anumbr
  GO TO 103
END IF

104     WRITE(*,5) ttime
5       FORMAT(' TTIME=',g15.6)
CALL okay(item,anumbr,1.e35,-1.e35,intnum,same,abort)
IF (abort) THEN
  rtn_ederr = 2
  RETURN
END IF
IF (.NOT.same) THEN
  ttime=anumbr
  GO TO 104
END IF

!------- The freeze and scan options are displayed for each superblock.

150     DO  k=1,decsup
  WRITE(*,6)k
  6         FORMAT(/,' SUPERBLOCK',i2,/)
  201       WRITE(*,7)freeze(k)
  7         FORMAT(' FREEZE OPTION',i2)
  CALL okay(item,anumbr,2.,0.,intnum,same,abort)
  IF (abort) THEN
    rtn_ederr = 2
    RETURN
  END IF
  IF (.NOT.same) THEN
    freeze(k)=intnum
    GO TO 201
  END IF
  202       WRITE(*,8)scan(k)
  8         FORMAT(' SCAN OPTION',i2)
  CALL okay(item,anumbr,1.,0.,intnum,same,abort)
  IF (abort) THEN
    rtn_ederr = 2
    RETURN
  END IF
  IF (.NOT.same) THEN
    scan(k)=intnum
    GO TO 202
  END IF
END DO

GO TO 100
300     CALL fsave(crtflg,edflg,rtn_fsave)
if(rtn_fsave == 1) then
  rtn_ederr = 1
elseif(rtn_fsave == 2) then
  rtn_ederr = 2
endif
RETURN
END SUBROUTINE ederr
!=======================================================================

!       SUBROUTINE REEDT

!------ This subroutine is called by the EDIT module and handles
!       all calls to other selected modules. This subroutine will
!       only allow calls to the  HELP module from the EDIT module.

!=======================================================================

SUBROUTINE reedt(item,crtflg,edflg,rtn_reedt)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=1) :: sure
CHARACTER (LEN=2) :: answer

answer=item
IF (answer == 'AB' .OR. answer == 'ab') GO TO 250
IF (answer == 'HE' .OR. answer == 'he') GO TO 101
RETURN

!------ A call is made to HELP with program control transferred there.

101     CALL help(crtflg,edflg,abort)
rtn_reedt = 1
200     RETURN
250     messag=' Are you sure you wish to ABORT, Y/N'
CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
sure=item
IF (sure /= 'Y' .OR. sure /= 'y') GO TO 200
rtn_reedt = 2
RETURN
END SUBROUTINE reedt

