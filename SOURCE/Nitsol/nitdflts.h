
c  This include file contains defaults for various parameters that
c  control the nonlinear iterations.

c  Note that use of parameter statements here is not necessary, it
c  merely facilitates conversion to single precision:  merely change
c  each double precision declaration to real, and edit the definitions
c  of the constants in the first two parameter statements.

      real one,       two,       rtfiv
      parameter      ( one=1.0e0, two=2.0e0, rtfiv=2.23606797749978981 )
  
      real tenth,        half,        fournines
      parameter      ( tenth=0.10e0, half=0.50e0, fournines=one-1.0e-4 )
  
      real DFLT_CHOICE1_EXP
      parameter      ( DFLT_CHOICE1_EXP=(one+rtfiv)*half )
  
      real DFLT_CHOICE2_EXP
      parameter      ( DFLT_CHOICE2_EXP=two )
  
      real DFLT_CHOICE2_COEF
      parameter      ( DFLT_CHOICE2_COEF=one )
  
      real DFLT_ETA_CUTOFF
      parameter      ( DFLT_ETA_CUTOFF=tenth )
  
      real DFLT_ETA_MAX
      parameter      ( DFLT_ETA_MAX=fournines )
  
      real DFLT_THMIN
      parameter      ( DFLT_THMIN=tenth )
  
      real DFLT_THMAX
      parameter      ( DFLT_THMAX=half )
  
      real DFLT_ETA_FIXED
      parameter      ( DFLT_ETA_FIXED=tenth )
c	  parameter		 ( DFLT_ETA_FIXED=0.0001 )

      integer     DFLT_PRLVL
      parameter ( DFLT_PRLVL=0 )

      integer     STDOUT
      parameter ( STDOUT=6 )
