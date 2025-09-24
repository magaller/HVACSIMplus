                            =========
                             READ.ME
                            =========


  DISCLAIMER

  This program is furnished by the government and used by any
  recipient with the express understanding that the United
  States Government makes no warranty, expressed or implied,
  concerning the accuracy, completeness, reliability, usability,
  or suitability for any particular purpose of the information and data
  contained in this program or furnished in connection therewith,
  and the United States shall be under no liability whatsoever
  to any person by reason of any use made thereof.  This program
  belongs to the government.  Therefore, the recipient
  further agrees not to assert any proprietary rights therein or
  to represent this program to anyone as other than a government
  program.


System Requirements

* PC system computer or compatible computer
* 256 megabytes of memory 
* 100 megabytes of hard disk space 
* Fortran 90/95 compiler for development of source code
* PDF file reader for reading manuals 

Installing HVACSIM+    

1. Create a folder/directory (e.g., HVACSIM20) on a hard drive.
2. Copy all files on the CD into the folder.
3. Make a working subfolder (e.g., WORK) in the folder.
4. Copy the files in the folders, BIN and DATA, into the working folder for an initial run.
5. The folder DOC contains manuals in PDF format, and the folder SAMPLE contains the inputs for sample runs.

Manuals

* NBSIR 84-2996.pdf
Clark, D.R., HVACSIM+ Building Systems and Equipment Simulation Program Reference Manual, NBSIR 84-2996, NIST, Jan. 1985.
* NBSIR 85-3243.pdf
Clark, D.R., and W.B. May, HVACSIM+ Building Systems and Equipment Simulation Program Users Guide, NBSIR 85-3243, NIST, Sept.1985.
* NBSIR 86-3331.pdf
Park, C., D.R. Clark, and  G.E. Kelly, HVACSIM+ Building Systems and Equipment Simulation Program Building Loads Calculation, NBSIR 86-3331, NIST, Feb.1986.
* NISTIR 7514.pdf
      Park, C., HVACSIM+ User’s Guide Update, NISTIR 7514, NIST, July 2008.

To use HVACSIM+ version 20.0, read NISTIR 7514.pdf first, then refer to other manuals for details.

Note

(1) Every component model as TYPEn should include at least one parameter variable, PAR.

(2) Whether used or not, the boundary data file (default file name: HVACSIM.BND) must be present during simulations.

(3) As a default, the simulation summary file (default file name: HVACSIM.SUM) contains simulation outputs up to 3600 seconds.  If desired, change the value of WRTMAX in the file MODSIM_HEAD.F90.






