!=======================================================================
! HVACGEN SOURCE FILE #1 OF 6; VERSION 5.0
!=======================================================================

! hvacgen1.f90

! This source code, when compiled with a FORTRAN 90/95 compiler, and linked
! with the compiled source code from the other five source code files,
! will produce the program HVACGEN which is the "front-end" program for
! the HVACSIM+ dynamic simulation package for building heating,
! ventilation, and air conditions systems and controls, plus other
! building systems. This program produces simulation work files which
! must be processed by the program SLIMCON to produce a model definition
! file which describes the structure and characteristics of the system
! to be simulated by the main simulation program, MODSIM.

! - Original concepts based on the program SIMCON, written by C.R. Hill,
!   which was a much simplified version of HVACGEN, but performed the
!   same function. Where SIMCON only allowed the creation of work files,
!   HVACGEN allows the editing and well as creation of work files.

! - HVACGEN was designed by:   W.B. May and D.B. Harris
!                              Building 226, Room B122
!                              National Institue of Standards and
!                                Technology (formerly NBS)
!                              Gaithersburg MD, 20899
!                              U.S.A.

! - HVACGEN was coded by:      D.B. Harris
! - The first complete version was 1.3 and was finished MARCH 1984.

! - changes/improvements

! - version 1.4, APRIL 1984: W.B. MAY, improved section for file name
!                entry after editing of work file (HVACGEN3.FTN)

! - version 1.5, SEPT. 1984: W.B. MAY, general cleanup, removal of dead
!                code, improved command parser, improved type
!                data file access, ability to change units for
!                input/output variables in data statement,
!                types file contents list, variable index
!                zeroing on insert unit and insert block.
! - version 1.6, OCT.  1984: W.B. MAY, consolidation of all references
!                to the types information file into a single
!                subroutine, TYPEIN. This subroutine has an
!                option to allow direct access file reads.
!                General cleanup of non-standard (77) fortran.
!                Improved boundary and reported variable
!                editing. Standardization of common blocks.
!                Installation of all common array dimensions
!                in parameter statements. Improvements to
!                editing of boundary conditions and reported
!                            variables.
! - version 1.7, APR.  1985: Removal of many alternate return statements
!                to improve clarity. Changed length of output
!                message in DATAIN. Revised OKAY. Made string
!                argument to RITE of variable length. Removed
!                character variables from common SIMINF and
!                VEWINF and placed in SIMIN2 and VEWIN2.
!                Display of types with TYPES was revised to
!                pause after each screen full. Write statements
!                were revised to keep output lines less than
!                72 characters; previously some long lines were
!                output. Insert Block and Unit functions were
!                revised to allow a positive or negative offset
!                to be added to the variable indices when read
!                from a work file. A VIEW SIMULATION <filename>
!                ALL function was added to allow a simulation
!                work file to be documented with one command.
!                A SAVE BLOCK function was added to the EDIT
!                Simulation section. This allows block work
!                files to be created from blocks in simulation
!                work files. The EDIT STRUCTURE function was
!                modified to allow sequential editing of units
!                without invokation of the EDIT STRUCTURE command
!                each time.
! - version 1.8, AUG. 16, 1985:
!                C. PARK and D. R. Clark
!                TYPEIN subroutine  was modified.

! - version 1.9, April 30, 1986:
!                C. PARK
!                (1) The parameter value of MRPTIS was changed to 30
!                    to be consistent with SLIMCON and MODSIM.
!                (2) The subroutine RPINRT was revised as
!                     NINDX(NSUPER,1)=INTNUM
!                     for all catagories.
!                (3) The numerical constant for ANUMAX which is an
!                    argument of both DATAIN and OKAY subroutines
!                    was changed into a real varialbe.
! - version 1.9A, Sept. 18, 1986:
!                C. PARK
!                FORMAT statements for INDEX numbers were changed
!                from I2 to I3.
! - version 1.9B, February 17, 1987:
!                C. PARK
!                BLKUNT(MUNTIB) was changed into BLKUNT(MAXBLK).
!                Set  MAXPAR=1200 and MAXDEQ=90.
!                Errors in entering the values of error tolerances,
!                integration period, and boundary variable indices
!                were corrected when a simulation file was initially
!                created.  Because of the corrections, the subroutine
!                CRSIM2 was modified and the parameter MDMBND=MAXBND*8
!                was used.

! - version 2.0, April 7, 1989:
!                C. Park, National Inst. of Standards and Technology,
!                         U.S.A.
!                P. Haves, Oxford University, United Kingdom

!                (1) All READs and WRITEs on unit LUC now use unit *.
!                (2) Changes in PARAMETER values.
!                (3) Component identification number can go up to 300.

! - version 2.1  January 30, 1990:
!                C. Park
!                The displays on a screen for VIEW ALL can be saved on
!                the file, VIEWSAVE.TXT as a text file.

! - version 2.2. October 11, 1990:
!                C. Park
!                The dimension of SIMCNT was changed into MAXSBK+1.
!                Added comments on the simulation setup file (*.sim).

! - version 2.3  February 12, 1992:
!                C. Park
!                The term "abolute humidity" is changed into "humidity"
!                The term "OTHR" is changed into "OTHR"

! - version 3.0  September 14, 1992:
!                Cheol Park
!                (1) The INCLUDE statement ( non-ANSI standard ) for
!                    PARAMETERs is used.
!                (2) The order of INPUTs and OUTPUTs in the TYPAR.DAT
!                    file is reversed and a # sign is inserted between
!                    the INPUTs and OUTPUTs for clarity.
!                (3) The direct access file DIRECT.DAT which is created
!                    from the sequential TYPAR.DAT is saved for
!                    debugging the error due to the format of  TYPAR.DAT.

! - version 3.1  March 4, 1993:
!                Minor changes were made based on the comments from
!                Zhu and Iwami, Nagoya University, Japan

! - version 3.2  March 31, 1994:
!                The file "HVACGEN.PAR" contains COMMON blocks and type
!                statements in addition to PARAMETER statements.

! - version 3.3  July 7, 1996:
!                Changes made by P Haves, Loughborough University, UK to
!                support the "simulation test-bed" produced by ASHRAE
!                825-RP
!                (1) Long file names supported (46 chars + '.' + 3 char
!                    extension = 50 chars)
!                (2) SUBROUTINE DATAIN: use a local copy of argument
!                    SWITCH since calling routines use constants for
!                    corresponding argument in CALL statement ("access
!                    violation" error on VAX-VMS)
!                (3) PARAMETER MAXTYP set in hvacsim.par
!                (4) SUBROUTINE TYPEIN: tries to open TYPES data file
!                    TYPARF in current directory or file TYPARFF
!                    specified with full path
!                (5) SUBROUTINEs BOUND, RPTVAR and EDVAL: one extra
!                    space introduced between the variable type
!                    (PRESSURE etc) and the index to match the spacing
!                    produced by SUBROUTINE VEWUNT. This facilitates
!                    using a text editor to search for occurancies of a
!                    particular simulation variable in different parts
!                    of a text file generated using the VIEW facility.

! - version 3.4  June 23, 1997:

!                The file "hvacgen.par" is replaced with "hvacsim.par",
!                and declarations are restored in each subroutines.

! - version 4.0  December 20, 1999  C.P.

!                An include file, hvaccomm.inc, is added to replace
!                most of declarations in subprograms for clarity.

! - version 5.0  March 9, 2007 Cheol Park

!                At first, code was converted using TO_F90 by Alan Miller
!                into Fortran 90 code, then alternate subroutine returns
!                were revised for Fortran 90 compatible.
!                Minor modifications also were made.

! ----------------------------------------------------------------------

!  DISCLAIMER

!  This program is furnished by the government and used by any recipient
!  with the express understanding that the United States Government
!  makes no warranty, expressed or implied, concerning the accuracy,
!  completeness, reliablilty, usability, or suitability for any
!  particular purpose of the information and data contained in
!  this program or furnished in connection therewith, and the United
!  States shall be under no liability whatsoever to any person by reason
!  of any use made thereof.  This program belongs to the government.
!  Therefore, the recipient further agrees not to assert any
!  proprietary rights therein or to represent this program to anyone
!  as other than a government program.

!=======================================================================

!        MAIN PROGRAM

!=======================================================================
PROGRAM hvacgen
!=======================================================================
!------- This portion of the "front end" routine is used to
!        transfer control to the VIEW, HELP, EDIT, and CREATE
!        modules. This module as all other modules
!        will receive user input through the use of the Command
!        Processor Module (COPMOD). This module will also allow the
!        user to exit the "front end" program at any time.

!------- LUF is the logical unit for the data files
!------- LUFB is a logical unit number that may be used for
!        a later option.
!------- LUTYPR is the logical unit of the types file TYPAR.DAT.

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=8)  :: vers
CHARACTER (LEN=20) :: update

 crtflg = 0
 edflg  = 0

vers = '5.0    '
update = 'March 9, 2007'

do

WRITE(*,1)
1   FORMAT(/////)
WRITE(*,2) vers,update
2   FORMAT (' HVACGEN - Simulation GENeration Program',  &
            //,' Version ',a8,1X,a20,//)
CALL rite(' Choose from the list below:')
WRITE(*,3)
3   FORMAT (//,' CReate (SImulation, BLock, UNit)',&
            //,' EDit (SImulation, UNit)', &
            //,' VIew (SImulation, BLock, UNit)',  &
            //,' HElp',                            &
            //,' ENd',/)

messag=' Selection ?'
CALL datain(item,anumbr,0.,0.,intnum,line,report,2,1,1,messag)
CALL remain(item,crtflg,edflg,abort)
enddo
END PROGRAM hvacgen

!=======================================================================

!        MODULE DATAIN

!        This module is called by any module or subroutine that
!        requires input from the console.
!        It is responsible for handling all calls to the necessary
!        subroutines which are used to input data at the keyboard
!        and handle all errors that may occur with such input.
!------- This subroutine examines the variable STCHEK, which is
!        passed by the first calling program, and it is compared
!        with the variable REPORT, which is passed by subroutine
!        COPMOD. If the two are not equal an error message is
!        printed and the original calling program is instructed
!        to ask the question again. SWITCH is set to 2 so that
!        a new line will be read in on the re-enter cycle.
!
!        This is specified by the user when calling the subroutine
!        DATAIN.
!
!------- ANUMAX is a maximum range that may be placed on, ANUMBR.
!------- ANUMIN is a minimum range on ANUMBR.
!------- INTNUM is the integer value of ANUMBR.
!------- STCHEK is a variable that is passed by the first calling
!        program that specifies what should have been entered at
!        the key board. This may take on any of the values for
!        REPORT found in subroutine COPMOD except 4, as this
!        value is not applicable as an input request. Also note
!        that if 0 is used no error message shall be printed and
!        there will be no way of re-entering the information entered
!------- RESET should be set to 1 by the calling program when entering
!        a word or phrase at the console, all other times it should be
!        set to 0. When the value is 1, RESET sets the value of ANUMBR
!        back to 0 avoiding any possible errors. If STCHEK is set to
!        0 this has no effect.
!------- MESSAG allows the user to display a message during input
!        from the console. If no message is desired this should be
!        to a ' ', space.
!------- All other variables are defined in subroutine COPMOD.

!=======================================================================

SUBROUTINE datain(item,anumbr,anumax,anumin,intnum,line,  &
    report,switch,stchek,reset,messag)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

REAL, INTENT(IN )                     :: anumax
REAL, INTENT(IN )                     :: anumin
INTEGER, INTENT(IN)                      :: switch
INTEGER, INTENT(IN )                  :: stchek
INTEGER, INTENT(IN )                  :: reset
LOGICAL :: reword
INTEGER :: localswi
character(len=1),dimension(72)           :: dash='-'

!------- Local copy of SWITCH - prevents potential runtime error since
!        corresponding arguments in CALL statements are constants

localswi=switch

!------- MESSAG is examined to determine if it should be displayed

100      WRITE(*,fmt='(72a1)') dash             ! 6/6/07
         WRITE(*,fmt='(a72)') messag
!100      WRITE(*,fmt='(a72)') messag

CALL copmod(item,anumbr,line,report,localswi)
IF (.NOT.reword(item)) THEN
  IF (reset == 1) anumbr=0.
  intnum=anumbr
  
!------- STCHEK and REPORT are compared to determine if the information
!        entered was the same as that requested.
!        If not appropriate error messages are displayed and the user
!        is requested to re-enter the same information.
  
  IF (stchek /= 0) THEN
    IF (stchek == report) THEN
      IF (anumbr > anumax.OR.anumbr < anumin) THEN
        CALL rite(' Number out of allowable range')
        localswi=2
        GO TO 100
      END IF
      RETURN
    ELSE IF (report == 1) THEN
      CALL rite(' A Word was entered - A Number is required')
    ELSE IF (report == 2.OR.report == 4) THEN
      CALL rite (' End of line encountered - Reenter information')
    ELSE IF (report == 3) THEN
      CALL rite(' A Number was entered - A Word is required')
    ELSE IF (report == 5) THEN
      CALL rite(' A Reserved word was entered, Re-enter the')
      CALL rite(' information requested.')
    ELSE
      GO TO 100
    END IF

    IF (stchek == 2) then
      localswi=3
    else
      localswi=2
    endif

    GO TO 100
  END IF
END IF
RETURN
END SUBROUTINE datain
!=======================================================================

!       SUBROUTINE COPMOD

!------ This subroutine is used to break down an inputted line
!       from the console into words and numbers. This is accomplished
!       by reading the line character by character up until a space
!       is detected. During the character read if at any time an
!       alpha-numeric symbol is detected the group of characters
!       is considered a word. All other cases result in a number
!       and a special check is performed on a group of characters
!       that contain an "E" to determine if the group is a word
!       or an exponential number.

!------ CHARC is the variable assigned to the characters in the line
!       that is read in from the console.
!------ ITEM is a character variable assigned to a word or number
!       if it is assigned to a number it is decoded before being
!       passed back to the calling program.
!------ LINE is the character variable assigned to the line read
!       in from the keyboard.
!------ I counting variable.
!------ J also a counting variable.
!------ TEST variable which is assigned to the variable CHARAC
!       used during a call to the subroutine CHECK.
!------ TEST1 and TEST2 are variables used to check the characters
!       before and after CHARAC, these are used during the
!       exponential check.
!------ STATUS is defined in subroutine CHECK.
!------ STAT1 and STAT2 return the conclusion of the numerical
!       and alpha-numeric symbol test from the subroutine CHECK. They
!       are used in conjunction with TEST1 and TEST2 respectively.
!------ REPORT a variable that is assigned a value of 1,2,3, or 4 where
!       1 is the end of a word,, 2 indicates the return of an
!       entire inputted line to the calling program, 3 is the end
!       of a number, 4 is the end of a line, and 5 indicates a reserve
!       word was entered and the information should be re-entered.
!       This value is passed back to the calling program.
!------ ANUMBR is the variable that is assigned to ITEM after
!       decoding and represents a number when passed back to
!       calling program.
!------ SWITCH is a variable used by the calling program which
!       is either 1,2 or 3; 1 instructs COPMOD to read from the
!       existing line on the console, 2 instructs it to read
!       a new line, and 3 results in the entire line being read in
!       and passed back to the calling program.
!------ DBCHEK is a variable which is used as a double check to
!       make certain that a number does not contain any other
!       alpha-numeric symbols with the exception of an "E". If
!       a number does than it is then treated as a word.
!------ ECHECK is used during the exponential checking procedure
!       to be certain that the number does not contain more than
!       one "E", if it does it is then considered a word.
!==================================================================

SUBROUTINE copmod(item,anumbr,line,report,switch)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

INTEGER, INTENT(IN OUT)      :: switch
CHARACTER (LEN=1)            :: charac,charc1,charc2
CHARACTER (LEN=80)           :: origin
INTEGER                      :: stat1,stat2, dbchek,echeck
LOGICAL                      :: NUMBER
integer                      :: nperod
DATA origin(1:40)  /'                                        '/,  &
     origin(41:80) /'                                        '/
SAVE

!---------Initialize flags and buffers
i = 1
anumbr = 0.
dbchek=0
echeck=0
nperod=0
NUMBER = .true.
item='                                                     '
line = origin
j=1

!----------Read in a new line if requested

IF (switch /= 1) THEN
  write(*,fmt='(a4)',advance='no') ' => '
  READ(*,fmt='(a80)') origin
  line = origin
  IF (i > 80) GO TO 106
  i=1
END IF

!-------- CHARAC assigned to a character in LINE

100     charac=line(i:i)
IF (.NOT.NUMBER) GO TO 104

!--------- Number checks

charc1=line(i+1:i+1)
IF (i /= 1) THEN
  charc2=line(i-1:i-1)
ELSE
  charc2 = charac
END IF
IF (i > 80) GO TO 106
IF (charc1 == '-'.AND.charac == '-'.OR.charc1 == '-'.AND.  &
      charac == '+'.OR.charc1 == '+'.AND.charac == '+'.OR.charc1 ==  &
      '+'.AND.charac == '-') THEN
  NUMBER = .false.
  GO TO 104
END IF
IF (charac == '-'.OR.charac == '+') THEN
  IF (charc1 > '9'.OR.charc1 < '0'.AND.charc1 /= '.'.OR.charc2  &
         >= '0'.AND.charc2 <= '9'.OR.charc2 == '.') THEN
    NUMBER = .false.
    GO TO 104
  END IF
END IF

IF (charac == '.') nperod=nperod+1
IF (nperod >= 2) THEN
  NUMBER = .false.
  GO TO 104
END IF

!------ Examination of CHARAC is performed to determine if it
!       is numerical or an alpha-numeric symbol.

CALL check(charac,STATUS)
IF (STATUS == 0.AND.charac /= 'E'.AND.charac /= 'e') dbchek=1
IF (STATUS == 0) THEN
  NUMBER = .false.
  
!----- Examination of CHARAC to determine if an "E" is part of
!      a word or if it is used in an exponent.
!----- ECHECK is used for the multiple "E" check
  
  IF (charac == 'E' .OR. charac == 'e') THEN
    echeck=echeck+1
    IF (echeck >= 2) THEN
      NUMBER = .false.
      GO TO 104
    END IF
    
    IF ((charac == 'E'.AND.charc1 == 'E'.OR.charac  &
           == 'E'.AND.charc2 == 'E') .OR.  &
          (charac == 'e'.AND.charc1 == 'e'.OR.charac  &
           == 'e'.AND.charc2 == 'e')) THEN
      dbchek=1
    END IF
    CALL check(charc1,stat1)
    CALL check(charc2,stat2)
    IF (stat1 /= 1.OR.stat2 /= 1) THEN
      NUMBER = .false.
    ELSE
      NUMBER = .true.
    END IF
  END IF
END IF

!----- Check for end of word/number

104   IF (i > 80) GO TO 106
IF (charac == ' '.OR.charac == ','.OR.j == 50) THEN
  IF ((dbchek == 1.AND.charac /= 'E') .OR.  &
      (dbchek == 1.AND.charac /= 'e'))  NUMBER = .false.
  IF (.NOT.NUMBER) THEN
    report=1
  ELSE IF (NUMBER) THEN
    
!------ ITEM is decoded and re-assigned to ANUMBR before being passed
!       back to the calling program.
    
    IF (j == 1.AND.item(j:j) == ' ') THEN
      anumbr = 0.
    ELSE
      IF (nperod == 0.AND.echeck == 0.AND.j /= 1.AND.j /= 50) item(j:j) = '.'
      READ(item,fmt='(f50.0)') anumbr
    END IF
    report=3
  END IF
  GO TO 107
END IF
item(j:j)=line(i:i)
i=i+1
j=j+1
GO TO 100

106   report=4

107   i=i+1
j=j+1
IF (switch == 3) report=2
RETURN
END SUBROUTINE copmod
!=======================================================================

!     SUBROUTINE CHECK

!-------- This subroutine determines whether "a" character is a
!         number or an alpha-numeric symbol. The decimal point
!         and plus and minus sign are considered numbers in this
!         subroutine.
!         Checking an exponent is done by the calling program.

!-------- TEST is the character that is tested
!-------- STATUS is a number that is passed back to the calling
!         program; a value of 1 signifies a number, a value of 0
!         signifies an alpha-numeric symbol.
!=======================================================================

SUBROUTINE check(test,STATUS)

implicit none
CHARACTER (LEN=1), INTENT(IN OUT)        :: test
INTEGER, INTENT(OUT)                     :: STATUS

IF (test == ' '.OR.test == ',') THEN
  STATUS = 1
ELSE IF (test >= '0'.AND.test <= '9') THEN
  STATUS = 1
ELSE IF (test == '+'.OR.test == '-'.OR.test == '.') THEN
  STATUS = 1
ELSE
  STATUS = 0
END IF

RETURN
END SUBROUTINE check
!=======================================================================

!        FUNCTION REWORD

!------- This subroutine is called after a call has been made
!        to the subroutine COPMOD and its purpose is to check
!        the latest input to see if it is a reserved word used
!        for calling other modules or performing other functions

!------- RECHEK is assigned to the word that is passed by the
!        calling program. It is used to determine if the passed word
!        is a reserved word or not.
!------- STAT is a variable that may have the value 0 or 1.
!        1 indicates that the word that was checked was a reserved
!        word and 0 indicates that it was not.
!------- EXAM is an array used to store the different reserve words
!        as allocated by the data statement.

!======================================================================

LOGICAL FUNCTION reword(item)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=2)                :: rechek
CHARACTER (LEN=6), dimension(14) :: exam =  &  !------- Reserve word data
      (/'CREATE','EDIT  ','VIEW  ','HELP  ', 'ABORT ','END   ','TYPES ', &
        'create','edit  ','view  ','help  ', 'abort ','end   ','types '/)

!------- Check to see if RECHEK is a reserve word

rechek=item
DO  i=1,14
  IF (rechek == exam(i)(1:2)) THEN
    reword = .true.
    RETURN
  END IF
END DO
reword = .false.

RETURN
END FUNCTION reword
!======================================================================

!       SUBROUTINE REMAIN

!------ This subroutine is called by the MAIN module and its
!       purpose is to transfer control to the requested module.
!       If a module is selected that is nonexistent then the program
!       will print an error message and transfer control back to
!       the MAIN module.
!------ POINTR is used as a stack point counter to keep track
!       of how many subroutines have been called and also
!       tells the program the position of the user away from the
!       MAIN calling module.
!=======================================================================

SUBROUTINE remain(item,crtflg,edflg,abort)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

CHARACTER (LEN=2) :: answer

answer=item
IF (answer == 'CR'.OR.answer == 'cr') THEN
  CALL create(crtflg,edflg,abort)
ELSE IF (answer == 'ED'.OR.answer == 'ed') THEN
  CALL edit(crtflg,edflg,abort)
ELSE IF (answer == 'VI'.OR.answer == 'vi') THEN
  CALL view(crtflg,edflg,abort)
ELSE IF (answer == 'HE'.OR.answer == 'he') THEN
  CALL help(crtflg,edflg,abort)
ELSE IF (answer == 'EN'.OR.answer == 'en') THEN
  close(lufb,status='delete')
  STOP
ELSE
  CALL rite(' Choose only from the given list.')
  CALL holdit
END IF
RETURN
END SUBROUTINE remain
!=======================================================================

!      SUBROUTINE HOLDIT

!      This routine causes a pause

!=======================================================================

SUBROUTINE holdit

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

messag=' Push the Carriage Return to continue'
CALL datain(item,anumbr,1.,0.,intnum,line,report,2,3,0,messag)
RETURN
END SUBROUTINE holdit
!=======================================================================

!     SUBROUTINE SCROLL

!---- This subroutine uses a variable LINCNT to determine when a screen
! full of information has been displayed. Display is paused after each
! screen.

!=======================================================================

SUBROUTINE scroll(lincnt)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

lincnt=lincnt+1
IF (lincnt <= nlscre) RETURN
100   messag=' Hit the CR to scroll to next screen'
CALL datain(item,anumbr,1.,0.,intnum,line,report,2,3,0,messag)
IF (intnum == 0) THEN
  lincnt=1
  RETURN
ELSE
  GO TO 100
END IF
RETURN
END SUBROUTINE scroll
!=======================================================================

!       SUBROUTINE OKAY

!------- This subroutine prompts the user for the correct value and
!        allows the user to make changes to it.

!=======================================================================

SUBROUTINE okay(item,anumbr,anumax,anumin,intnum,same,abort)

use hvacsim_par
use hvaccomm
implicit none
include 'hvacgen.inc'

REAL, INTENT(IN)           :: anumax
REAL, INTENT(IN)           :: anumin
CHARACTER (LEN=1)              :: answer

abort = .false.
200      messag=' OK?, Y/N '

!------- IOK is set to an arbitrary value so that the IF statements
!        will respond correctly,i.e., YES is equal to 'Y' or just
!        hitting the carriage return.

iok=-13
CALL datain(item,anumbr,1.,0.,iok,line,report,2,0,0,messag)
CALL reedt(item,crtflg,edflg,rtn_reedt)
IF (rtn_reedt == 2) GO TO 200
!IF (abort) RETURN
answer=item
IF (iok == 0.AND.(answer == 'N'.OR.answer == 'n')) THEN
  same = .false.
  201  messag=' Enter the new value'
  CALL datain(item,anumbr,anumax,anumin,intnum,line,report, 2,3,0,messag)
  CALL reedt(item,crtflg,edflg,rtn_reedt)
    IF (rtn_reedt == 2) GO TO 201
ELSE IF ((answer == 'Y'.OR.answer == 'y').OR.iok == 0) THEN
  same = .true.
  RETURN
ELSE
  GO TO 200
END IF
RETURN
END SUBROUTINE okay
!=======================================================================

!     SUBROUTINE RITE

!=======================================================================

SUBROUTINE rite(output)

implicit none
CHARACTER (LEN=*), INTENT(IN)            :: output
CHARACTER (LEN=80)                       :: linea

linea = output
WRITE(*,fmt='(a80)') linea

RETURN
END SUBROUTINE rite
!=======================================================================

!     SUBROUTINE PROMPTD

!------ This subroutine is used to provide the index labels for
!       the input and output information of each of the types.
!=======================================================================

SUBROUTINE promptd(ident,i,prmpt)

use hvacsim_par
implicit none

INTEGER, INTENT(IN OUT)                  :: ident(minoiu)
INTEGER, INTENT(IN)                      :: i
CHARACTER (LEN=15), INTENT(OUT)          :: prmpt
CHARACTER (LEN=15), dimension(8)         :: varcat = &
   (/' PRESSURE   ',' FLOW       ',' TEMPERATURE',' CONTROL    ',  &
     ' OTHER      ',' ENERGY     ',' POWER      ',' HUMIDITY   '/)

prmpt = varcat(ident(i))
RETURN
END SUBROUTINE promptd

