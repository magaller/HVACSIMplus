!***********************************************************************
!
!     SNSQ
!
!     Source :  CMLIB/SNLSE, National Bureau of Standards
!
! ----------------------------------------------------------------------
!
!      SUBROUTINE SNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,
!     &                 ML,MU,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,
!     &                 NJEV,R,LR,QTF,WA1,WA2,WA3,WA4)
!
!***BEGIN PROLOGUE  SNSQ
!***DATE WRITTEN   800301   (YYMMDD)
!***REVISION DATE  820801   (YYMMDD)
!***CATEGORY NO.  F2A
!***KEYWORDS  NONLINEAR SQUARE SYSTEM,POWELL HYBRID METHOD,ZERO
!***AUTHOR  HIEBERT, K. L., (SNLA)
!***PURPOSE  SNSQ finds to find a zero of a system of N nonlinear
!            functions in N variables by a modification of the Powell
!            hybrid method.  This code is the combination of the MINPACK
!            codes (Argonne) HYBRD and HYBRDJ.
!***DESCRIPTION
! 1. Purpose.
!
!       The purpose of SNSQ is to find a zero of a system of N non-
!       linear functions in N variables by a modification of the Powell
!       hybrid method.  The user must provide a subroutine which calcu-
!       lates the functions.  The user has the option of either to
!       provide a subroutine which calculates the Jacobian or to let the
!       code calculate it by a forward-difference approximation.
!       This code is the combination of the MINPACK codes (Argonne)
!       HYBRD and HYBRDJ.
!
!
! 2. Subroutine and Type Statements.
!
!       SUBROUTINE SNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,
!      *                 ML,MU,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,
!      *                 NJEV,R,LR,QTF,WA1,WA2,WA3,WA4)
!       INTEGER IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,NJEV,LR
!       REAL XTOL,EPSFCN,FACTOR
!       REAL X(N),FVEC(N),DIAG(N),FJAC(LDFJAC,N),R(LR),QTF(N),
!      *     WA1(N),WA2(N),WA3(N),WA4(N)
!       EXTERNAL FCN,JAC
!
!
! 3. Parameters.
!
!       Parameters designated as input parameters must be specified on
!       entry to SNSQ and are not changed on exit, while parameters
!       designated as output parameters need not be specified on entry
!       and are set to appropriate values on exit from SNSQ.
!
!       FCN is the name of the user-supplied subroutine which calculates
!         the functions.  FCN must be declared in an EXTERNAL statement
!         in the user calling program, and should be written as follows.
!
!         SUBROUTINE FCN(N,X,FVEC,IFLAG)
!         INTEGER N,IFLAG
!         REAL X(N),FVEC(N)
!         ----------
!         Calculate the functions at X and
!         return this vector in FVEC.
!         ----------
!         RETURN
!         END
!
!         The value of IFLAG should not be changed by FCN unless the
!         user wants to terminate execution of SNSQ.  In this case, set
!         IFLAG to a negative integer.
!
!       JAC is the name of the user-supplied subroutine which calculates
!         the Jacobian.  If IOPT=1, then JAC must be declared in an
!         EXTERNAL statement in the user calling program, and should be
!         written as follows.
!
!         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
!         INTEGER N,LDFJAC,IFLAG
!         REAL X(N),FVEC(N),FJAC(LDFJAC,N)
!         ----------
!         Calculate the Jacobian at X and return this
!         matrix in FJAC.  FVEC contains the function
!         values at X and should not be altered.
!         ----------
!         RETURN
!         END
!
!         The value of IFLAG should not be changed by JAC unless the
!         user wants to terminate execution of SNSQ.  In this case, set
!         IFLAG to a negative integer.
!
!         If IOPT=2, JACcan be ignored (treat it as a dummy argument).
!
!       IOPT is an input variable which specifies how the Jacobian will
!         be calculated.  If IOPT=1, then the user must supply the
!         Jacobian through the subroutine JAC.  If IOPT=2, then the
!         code will approximate the Jacobian by forward-differencing.
!
!       N is a positive integer input variable set to the number of
!         functions and variables.
!
!       X is an array of length N.  On input, X must contain an initial
!         estimate of the solution vector.  On output, X contains the
!         final estimate of the solution vector.
!
!       FVE! is an output array of length N which contains the functions
!         evaluated at the output X.
!
!       FJAC is an output N by N array which contains the orthogonal
!         matrix Q produced by the QR factorization of the final approx-
!         imate Jacobian.
!
!       LDFJAC is a positive integer input variable not less than N
!         which specifies the leading dimension of the array FJAC.
!
!       XTOL is a non-negative input variable.  Termination occurs when
!         the relative error between two consecutive iterates is at most
!         XTOL.  Therefore, XTOL measures the relative error desired in
!         the approximate solution.  Section 4 contains more details
!         about XTOL.
!
!       MAXFEV is a positive integer input variable.  Termination occurs
!         when the number of calls to FCN is at least MAXFEV by the end
!         of an iteration.
!
!       ML is a non-negative integer input variable which specifies the
!         number of subdiagonals within the band of the Jacobian matrix.
!         If the Jacobian is not banded or IOPT=1, set ML to at
!         least N - 1.
!
!       MU is a non-negative integer input variable which specifies the
!         number of superdiagonals within the band of the Jacobian
!         matrix.  If the Jacobian is not banded or IOPT=1, set MU to at
!         least N - 1.
!
!       EPSFCN is an input variable used in determining a suitable step
!         for the forward-difference approximation.  This approximation
!         assumes that the relative errors in the functions are of the
!         order of EPSFCN.  If EPSFCN is less than the machine preci-
!         sion, it is assumed that the relative errors in the functions
!         are of the order of the machine precision.  If IOPT=1, then
!         EPSFCN can be ignored (treat it as a dummy argument).
!
!       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
!         internally set.  If MODE = 2, DIAG must contain positive
!         entries that serve as implicit (multiplicative) scale factors
!         for the variables.
!
!       MODE is an integer input variable.  If MODE = 1, the variables
!         will be scaled internally.  If MODE = 2, the scaling is speci-
!         fied by the input DIAG.  Other values of MODE are equivalent
!         to MODE = 1.
!
!       FACTOR is a positive input variable used in determining the ini-
!         tial step bound.  This bound is set to the product of FACTOR
!         and the Euclidean norm of DIAG*X if nonzero, or else to FACTOR
!         itself.  In most cases FACTOR should lie in the interval
!         (.1,100.).  100. is a generally recommended value.
!
!       NPRINT is an integer input variable that enables controlled
!         printing of iterates if it is positive.  In this case, FCN is
!         called with IFLAG = 0 at the beginning of the first iteration
!         and every NPRINT iteration thereafter and immediately prior
!         to return, with X and FVE! available for printing. Appropriate
!         print statements must be added to FCN(see example).  If NPRINT
!         is not positive, no special calls of FCN with IFLAG = 0 are
!         made.
!
!       INFO is an integer output variable.  If the user has terminated
!         execution, INFO is set to the (negative) value of IFLAG.  See
!         description of FCN and JAC. Otherwise, INFO is set as follows.
!
!         INFO = 0  improper input parameters.
!
!         INFO = 1  relative error between two consecutive iterates is
!                   at most XTOL.
!
!         INFO = 2  number of calls to FCN has reached or exceeded
!                   MAXFEV.
!
!         INFO = 3  XTOL is too small.  No further improvement in the
!                   approximate solution X is possible.
!
!         INFO = 4  iteration is not making good progress, as measured
!                   by the improvement from the last five Jacobian eval-
!                   uations.
!
!         INFO = 5  iteration is not making good progress, as measured
!                   by the improvement from the last ten iterations.
!
!         Sections 4 and 5 contain more details about INFO.
!
!       NFEV is an integer output variable set to the number of calls to
!         FCN.
!
!       NJEV is an integer output variable set to the number of calls to
!         JAC. (If IOPT=2, then NJEV is set to zero.)
!
!       R is an output array of length LR which contains the upper
!         triangular matrix produced by the QR factorization of the
!         final approximate Jacobian, stored rowwise.
!
!       LR is a positive integer input variable not less than
!         (N*(N+1))/2.
!
!       QTF is an output array of length N which contains the vector
!         (Q TRANSPOSE)*FVEC.
!
!       WA1, WA2, WA3, and WA4 are work arrays of length N.
!
!
! 4. Successful Completion.
!
!       The accuracy of SNSQ is controlled by the convergence parameter
!       XTOL.  This parameter is used in a test which makes a comparison
!       between the approximation X and a solution XSOL.  SNSQ termi-
!       nates when the test is satisfied.  If the convergence parameter
!       is less than the machine precision (as defined by the function
!       R1MACH(4)), then SNSQ only attempts to satisfy the test
!       defined by the machine precision.  Further progress is not
!       usually possible.
!
!       The test assumes that the functions are reasonably well behaved,
!       and, if the Jacobian is supplied by the user, that the functions
!       and the Jacobian are coded consistently.  If these conditions
!       are not satisfied, then SNSQ may incorrectly indicate conver-
!       gence.  The coding of the Jacobian can be checked by the
!       subroutine CHKDER. If the Jacobian is coded correctly or IOPT=2,
!       then the validity of the answer can be checked, for example, by
!       rerunning SNSQ with a tighter tolerance.
!
!       Convergence Test.  If ENORM(Z) denotes the Euclidean norm of a
!         vector Z and D is the diagonal matrix whose entries are
!         defined by the array DIAG, then this test attempts to guaran-
!         tee that
!
!               ENORM(D*(X-XSOL)) .LE. XTOL*ENORM(D*XSOL).
!
!         If this condition is satisfied with XTOL = 10**(-K), then the
!         larger components of D*X have K significant decimal digits and
!         INFO is set to 1.  There is a danger that the smaller compo-
!         nents of D*X may have large relative errors, but the fast rate
!         of convergence of SNSQ usually avoids this possibility.
!         Unless high precision solutions are required, the recommended
!         value for XTOL is the square root of the machine precision.
!
!
! 5. Unsuccessful Completion.
!
!       Unsuccessful termination of SNSQ can be due to improper input
!       parameters, arithmetic interrupts, an excessive number of func-
!       tion evaluations, or lack of good progress.
!
!       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1,
!         or IOPT .GT. 2, or N .LE. 0, or LDFJAC .LT. N, or
!         XTOL .LT. 0.E0, or MAXFEV .LE. 0, or ML .LT. 0, or MU .LT. 0,
!         or FACTOR .LE. 0.E0, or LR .LT. (N*(N+1))/2.
!
!       Arithmetic Interrupts.  If these interrupts occur in the FCN
!         subroutine during an early stage of the computation, they may
!         be caused by an unacceptable choice of X by SNSQ.  In this
!         case, it may be possible to remedy the situation by rerunning
!         SNSQ with a smaller value of FACTOR.
!
!       Excessive Number of Function Evaluations.  A reasonable value
!         for MAXFEV is 100*(N+1) for IOPT=1 and 200*(N+1) for IOPT=2.
!         If the number of calls to FCN reaches MAXFEV, then this
!         indicates that the routine is converging very slowly as
!         measured by the progress of FVEC, and INFO is set to 2.  This
!         situation should be unusual because, as indicated below, lack
!         of good progress is usually diagnosed earlier by SNSQ,
!         causing termination with INFO = 4 or INFO = 5.
!
!       Lack of Good Progress.  SNSQ searches for a zero of the system
!         by minimizing the sum of the squares of the functions.  In so
!         doing, it can become trapped in a region where the minimum
!         does not correspond to a zero of the system and, in this situ-
!         ation, the iteration eventually fails to make good progress.
!         In particular, this will happen if the system does not have a
!         zero.  If the system has a zero, rerunning SNSQ from a dif-
!         ferent starting point may be helpful.
!
!
! 6. Characteristics of the Algorithm.
!
!       SNSQ is a modification of the Powell hybrid method.  Two of its
!       main characteristics involve the choice of the correction as a
!       convex combination of the Newton and scaled gradient directions,
!       and the updating of the Jacobian by the rank-1 method of Broy-
!       den.  The choice of the correction guarantees (under reasonable
!       conditions) global convergence for starting points far from the
!       solution and a fast rate of convergence.  The Jacobian is
!       calculated at the starting point by either the user-supplied
!       subroutine or a forward-difference approximation, but it is not
!       recalculated until the rank-1 method fails to produce satis-
!       factory progress.
!
!       Timing.  The time required by SNSQ to solve a given problem
!         depends on N, the behavior of the functions, the accuracy
!         requested, and the starting point.  The number of arithmetic
!         operations needed by SNSQ is about 11.5*(N**2) to process
!         each evaluation of the functions (call to FCN) and 1.3*(N**3)
!         to process each evaluation of the Jacobian (call to JAC,
!         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
!         the timing of SNSQ will be strongly influenced by the time
!         spent in FCN and JAC.
!
!       Storage.  SNSQ requires (3*N**2 + 17*N)/2 single precision
!         storage locations, in addition to the storage required by the
!         program.  There are no internally declared storage arrays.
!
!
! 7. Example.
!
!       The problem is to determine the values of X(1), X(2), ..., X(9),
!       which solve the system of tridiagonal equations
!
!       (3-2*X(1))*X(1)           -2*X(2)                   = -1
!               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
!                                   -X(8) + (3-2*X(9))*X(9) = -1
! !     **********
!
!       PROGRAM TEST(INPUT,OUTPUT,TAPE6=OUTPUT)
! !
! !     Driver for SNSQ example.
! !
!       INTEGER J,IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR,
!      *        NWRITE
!       REAL XTOL,EPSFCN,FACTOR,FNORM
!       REAL X(9),FVEC(9),DIAG(9),FJAC(9,9),R(45),QTF(9),
!      *     WA1(9),WA2(9),WA3(9),WA4(9)
!       REAL ENORM,R1MACH
!       EXTERNAL FCN
!       DATA NWRITE /6/
! !
!       IOPT = 2
!       N = 9
! !
! !     The following starting values provide a rough solution.
! !
!       DO 10 J = 1, 9
!          X(J) = -1.E0
!    10    CONTINUE
! !
!       LDFJAC = 9
!       LR = 45
! !
! !     Set XTOL to the square root of the machine precision.
! !     Unless high precision solutions are required,
! !     this is the recommended setting.
! !
!       XTOL = SQRT(R1MACH(4))
! !
!       MAXFEV = 2000
!       ML = 1
!       MU = 1
!       EPSFCN = 0.E0
!       MODE = 2
!       DO 20 J = 1, 9
!          DIAG(J) = 1.E0
!    20    CONTINUE
!       FACTOR = 1.E2
!       NPRINT = 0
! !
!       CALL SNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,ML,MU,
!      *           EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
!      *           R,LR,QTF,WA1,WA2,WA3,WA4)
!       FNORM = ENORM(N,FVEC)
!       WRITE (NWRITE,1000) FNORM,NFEV,INFO,(X(J),J=1,N)
!       STOP
!  1000 FORMAT (5X,31H FINAL L2 NORM OF THE RESIDUALS,E15.7 //
!      *        5X,31H NUMBER OF FUNCTION EVALUATIONS,I10 //
!      *        5X,15H EXIT PARAMETER,16X,I10 //
!      *        5X,27H FINAL APPROXIMATE SOLUTION // (5X,3E15.7))
!       END
!       SUBROUTINE FCN(N,X,FVEC,IFLAG)
!       INTEGER N,IFLAG
!       REAL X(N),FVEC(N)
!       INTEGER K
!       REAL ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
!       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
! !
!       IF (IFLAG .NE. 0) GO TO 5
! !
! !     Insert print statements here when NPRINT is positive.
! !
!       RETURN
!     5 CONTINUE
!       DO 10 K = 1, N
!          TEMP = (THREE - TWO*X(K))*X(K)
!          TEMP1 = ZERO
!          IF (K .NE. 1) TEMP1 = X(K-1)
!          TEMP2 = ZERO
!          IF (K .NE. N) TEMP2 = X(K+1)
!          FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
!    10    CONTINUE
!       RETURN
!       END
!
!       Results obtained with different compilers or machines
!       may be slightly different.
!
!       FINAL L2 NORM OF THE RESIDUALS  0.1192636E-07
!
!       NUMBER OF FUNCTION EVALUATIONS        14
!
!       EXIT PARAMETER                         1
!
!       FINAL APPROXIMATE SOLUTION
!
!       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
!       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
!       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
!***REFERENCES  POWELL, M. J. D.
!                 A HYBRID METHOD FOR NONLINEAR EQUATIONS.
!                 NUMERICAL METHODS FOR NONLINEAR ALGEBRAIC EQUATIONS,
!                 P. RABINOWITZ, EDITOR.  GORDON AND BREACH, 1970.
!***ROUTINES CALLED  DOGLEG,ENORM,FDJAC1,QFORM,QRFAC,R1MACH,R1MPYQ,
!                    R1UPDT,XERROR
!***END PROLOGUE  SNSQ
! **********************************************************************
!
      subroutine snsq(fcn,jac,iopt,n,x,fvec,fjac,ldfjac,xtol,maxfev,&
                      ml,mu,epsfcn,diag,mode,factor,nprint,info,nfev,&
                      njev,r,lr,qtf,wa1,wa2,wa3,wa4,isblk)
!
! ----------------------------------------------------------------------
!
!     SNSQ : Find a zero of a system of N nonlinear functions in
!            N variables by a modification of the Powell hybrid methos.
!
!            The original program was written by K.L. Hiebert
!
!            A minor modification was made by C. R. Hill and C. Park
!            August 27, 1984
!
!            Updated on : Oct. 25, 1984 C.P.
!   
!            Modified by P. Haves, Oxford Univ., U.K. to include
!            simulation times in warning/error messages
!            April 12, 1989
!                    
!            Updated to Fortran 90.  Dec. 27, 2006 Cheol Park, NIST
!
! **********************************************************************
!
      use modsim_head
      implicit none

      external fcn,jac
      logical                  :: jeval,sing
      integer                  :: isblk,njev,nfev,info,nprint,mu,&
                                  ml,maxfev,iopt
      integer                  :: mode,n,lr,ldfjac,iflag,j,iset,iter,ncsuc,&
                                  ncfail,nslow1,nslow2,i,l,jm1
      integer,dimension(1)     :: iwa
      real                     :: factor,epsfcn,epsmch,r1mach,xnorm,&
                                  fnorm,enorm,delta
      real                     :: sum,temp,pnorm,fnorm1,actred,prered,&
                                  ratio
      real                     :: xtol
      real                     :: one=1.0,p1=0.1,p5=0.5,p001=1.0e-3,&
                                  p0001=1.0e-4,zero=0.0
      real,dimension(n)        :: x,fvec,diag,qtf,wa1,wa2,wa3,wa4
      real,dimension(lr)       :: r
      real,dimension(ldfjac,n) :: fjac

      epsmch = r1mach(4)
      xnorm=0.
!
      info = 0
      iflag = 0
      nfev = 0
      njev = 0
!
!     Check the input parameters for errors.
!
      if (iopt < 1 .or. iopt > 2 .or.&
          n <= 0 .or. xtol < zero .or. maxfev <= 0&
          .or. ml < 0 .or. mu < 0 .or. factor <= zero&
          .or. ldfjac < n .or. lr < (n*(n + 1))/2) go to 300
      if (mode == 2) then
        do j = 1, n
          if (diag(j) <= zero) go to 300
        end do
      end if
!
!     Evaluate the function at the starting ponit
!     and calculate its norm.
!
      iflag = 1
      iset=1
      call fcn(n,x,fvec,iflag,isblk,iset)
      nfev = 1
      if (iflag < 0) go to 300
      fnorm = enorm(n,fvec)
!
!     Initialize iteration counter and monitors.
!
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
!
!     Beginning of the outer loop.
!
   30 continue
         jeval = .true.
!
!        Calculate the jacobian matrix.
!
         if (iopt == 2)  then
!
!           Code approximates the jacobian.
!
            iflag = 2
            iset=-1
            call fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,&
                        wa2,isblk,iset)
            nfev = nfev + min(ml+mu+1,n)
         else
!
!           User supplies jacobian.
!
!*           call jac(n,x,fvec,fjac,ldfjac,iflag)
!*           njev = njev+1
         end if
!
         if (iflag < 0) go to 300
!
!        Compute the qr factorization of the jacobian.
!
         call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
!
!        On the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
         if (iter == 1) then
           if (mode /= 2) then
             do j = 1, n
                diag(j) = wa2(j)                       
                if (wa2(j) == zero) then
                  diag(j) = one
                end if
             end do                                    
           end if
!                                                          
!          On the first iteration, calculate the norm of the scaled x
!          and initialize the step bound delta.            
!                                                          
           do j = 1, n                                     
              wa3(j) = diag(j)*x(j)                        
           end do                                          
           xnorm = enorm(n,wa3)                            
           delta = factor*xnorm                            
           if (delta == zero) then
             delta = factor
           end if
         end if
!
!        Form (q transpose)*fvec and store in qtf.
!
         do i = 1, n
            qtf(i) = fvec(i)
         end do
         do j = 1, n
            if (fjac(j,j) == zero) cycle
            sum = zero
            do i = j, n
               sum = sum + fjac(i,j)*qtf(i)
            end do
            temp = -sum/fjac(j,j)
            do i = j, n
               qtf(i) = qtf(i) + fjac(i,j)*temp
            end do
         end do
!
!        Copy the triangular factor of the qr factorization into r.
!
         sing = .false.
         do j = 1, n
            l = j
            jm1 = j - 1
            if (jm1 >= 1) then
              do i = 1, jm1
                 r(l) = fjac(i,j)      
                 l = l + n - i         
              end do                   
            end if
            r(l) = wa1(j)
            if (wa1(j) == zero) then
              sing = .true.
            end if
         end do
!
!        Accumulate the orthogonal factor in fjac.
!
         call qform(n,n,fjac,ldfjac,wa1)
!
!        Rescale if necessary.
!
         if (mode /= 2) then
           do j = 1, n
              diag(j) = amax1(diag(j),wa2(j))        
           end do                                    
         end if
!
!        Beginning of the inner loop.
!
  180    continue
!
!           If requested, call fcn to enable printing of iterates.
!
            if (nprint > 0) then
              iflag = 0
              iset=0                                                    
              if (mod(iter-1,nprint) == 0) then
                call fcn(n,x,fvec,iflag,isblk,iset)
              end if
              if (iflag < 0) go to 300                                  
            end if
!
!           Determine the direction p.
!
            call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
!
!           Store the direction p and x + p. calculate the norm of p.
!
            do j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
            end do
            pnorm = enorm(n,wa3)
!
!           On the first iteration, adjust the initial step bound.
!
            if (iter == 1) then
              delta = amin1(delta,pnorm)
            end if
!
!           Evaluate the function at x + p and calculate its norm.
!
            iflag = 1
            iset=0
            call fcn(n,wa2,wa4,iflag,isblk,iset)
            nfev = nfev + 1
            if (iflag < 0) go to 300
            fnorm1 = enorm(n,wa4)
!
!           Compute the scaled actual reduction.
!
            actred = -one
            if (fnorm1 < fnorm) then
              actred = one - (fnorm1/fnorm)**2
            end if
!
!           Compute the scaled predicted reduction.
!
            l = 1
            do i = 1, n
               sum = zero
               do j = i, n
                  sum = sum + r(l)*wa1(j)
                  l = l + 1
               end do
               wa3(i) = qtf(i) + sum
            end do
            temp = enorm(n,wa3)
            prered = zero
            if (temp < fnorm) then
              prered = one - (temp/fnorm)**2
            end if
!
!           Compute the ratio of the actual to the predicted
!           reduction.
!
            ratio = zero
            if (prered > zero) then
              ratio = actred/prered
            end if
!
!           Update the step bound.
!
            if (ratio >= p1) then
              ncfail = 0
              ncsuc = ncsuc + 1
              if (ratio >= p5 .or. ncsuc > 1) then
                 delta = amax1(delta,pnorm/p5)
              end if
              if (abs(ratio-one) <= p1) then
                delta = pnorm/p5
              end if
            else
              ncsuc = 0
              ncfail = ncfail + 1
              delta = p5*delta
            end if
!
!           Test for successful iteration.
!
            if (ratio >= p0001) then
!
!             Successful iteration. update x, fvec, and their norms.
!                                                                           
              do j = 1, n                                                   
                 x(j) = wa2(j)                                              
                 wa2(j) = diag(j)*x(j)                                      
                 fvec(j) = wa4(j)                                           
              end do                                                        
              xnorm = enorm(n,wa2)                                          
              fnorm = fnorm1                                                
              iter = iter + 1                                               
            end if
!
!           Determine the progress of the iteration.
!
            nslow1 = nslow1 + 1
            if (actred >= p001) nslow1 = 0
            if (jeval) nslow2 = nslow2 + 1
            if (actred >= p1) nslow2 = 0
!
!           Test for convergence.
!
            if (delta <= xtol*xnorm .or. fnorm == zero) info = 1
            if (info /= 0) go to 300
!
!           Tests for termination and stringent tolerances.
!
            if (nfev >= maxfev) info = 2
            if (p1*amax1(p1*delta,pnorm) <= epsmch*xnorm) info = 3
            if (nslow2 == 5) info = 4
            if (nslow1 == 10) info = 5
            if (info /= 0) go to 300
!
!           Criterion for recalculating jacobian.
!
            if (ncfail /= 2) then
!
!             Calculate the rank one modification to the jacobian
!             and update qtf if necessary.                                  
!                                                                           
              do j = 1, n                                                   
                 sum = zero                                                 
                 do i = 1, n                                                
                    sum = sum + fjac(i,j)*wa4(i)                            
                 end do                                                     
                 wa2(j) = (sum - wa3(j))/pnorm                              
                 wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)                  
                 if (ratio >= p0001) then
                   qtf(j) = sum
                 end if
              end do                                                        
!                                                                           
!             Compute the qr factorization of the updated jacobian.         
!                                                                           
              call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)                        
              call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)                          
              call r1mpyq(1,n,qtf,1,wa2,wa3)                                
!                                                                           
!             End of the inner loop.                                        
!                                                                           
              jeval = .false.                                               
              go to 180                                                     
            end if
!
!        End of the outer loop.
!
         go to 30
  300 continue
!
!     Termination, either normal or user imposed.
!
      if (iflag < 0) info = iflag
      iflag = 0
      iset=0
      if (nprint > 0) call fcn(n,x,fvec,iflag,isblk,iset)
      if (info < 0) print 1000, time, isblk       ! 3/1/96
      if (info == 0) print 2000, time, isblk
      if (info == 2) print 3000, time, isblk
      if (info == 3) print 4000, time, isblk
      if (info > 4) print 5000, time, isblk
!
1000  format(' time = ',f8.1,' snsq - execuation terminated because',&
                               ' user set iflag negative.')
2000  format(' time = ',f8.1,' snsq - invalid input parameter.')
3000  format(' time = ',f8.1,' snsq - too many function evaluation.')
4000  format(' time = ',f8.1,' snsq - xtol too small.',&
                               '  no further improvement possible.')
5000  format(' time = ',f8.1,' snsq - iteration not making good',&
                               ' progress.','isblk=',i3)    ! 3/1/96
!
      return
      end subroutine snsq

! **********************************************************************
!
      subroutine snsq1(fcn,jac,iopt,n,x,fvec,fjac,ldfjac,xtol,maxfev,&
                       ml,mu,epsfcn,diag,mode,factor,nprint,info,nfev,&
                       njev,r,lr,qtf,wa1,wa2,wa3,wa4,iblk,isblk)
!
! ----------------------------------------------------------------------
!
!     SNSQ1 : This program is essentially same as SNSQ.
!           Called by subroutine BLOCK. Calls FDJAC2, while SNSQ calls
!           FDJAC1. Note that the argument of FCN has IBLK and ISBLK.
!
!           The original program was written by K.L. Hiebert
!
!           A minor modification was made by C. R. Hill and C. Park
!           August 13, 1984
!
!           Updated for writing diagonastic information by D.C.,
!           March 8, 1985
!
!           Modified by P. Haves, Oxford Univ., U.K. to include
!           simulation times in warning/error messages
!           April 12, 1989
!                    
!            Updated to Fortran 90.  Dec. 27, 2006 Cheol Park, NIST
!
! **********************************************************************

      use modsim_head
      implicit none

      external fcn,jac
      logical                  :: jeval,sing
      integer                  :: iblk,isblk,njev,nfev,info,nprint,mu,&
                                  ml,maxfev,iopt
      integer                  :: mode,n,lr,ldfjac,iflag,j,iset,iter,ncsuc,&
                                  ncfail,nslow1,nslow2,i,l,jm1
      integer,dimension(1)     :: iwa
      real                     :: factor,epsfcn,epsmch,r1mach,xnorm,&
                                  fnorm,enorm,delta,xtnorm
      real                     :: sum,temp,pnorm,fnorm1,actred,prered,&
                                  ratio
      real                     :: xtol
      real                     :: one=1.0,p1=0.1,p5=0.5,p001=1.0e-3,&
                                  p0001=1.0e-4,zero=0.0
      real,dimension(n)        :: x,fvec,diag,qtf,wa1,wa2,wa3,wa4
      real,dimension(lr)       :: r
      real,dimension(ldfjac,n) :: fjac

      epsmch = r1mach(4)
      xnorm=0.
!
      info = 0
      iflag = 0
      nfev = 0
      njev = 0
!
!     Check the input parameters for errors.
!
      if (iopt < 1 .or. iopt > 2 .or.&
          n <= 0 .or. xtol < zero .or. maxfev <= 0&
          .or. ml < 0 .or. mu < 0 .or. factor <= zero&
          .or. ldfjac < n .or. lr < (n*(n + 1))/2) go to 300
      if (mode == 2) then
        do j = 1, n
           if (diag(j) <= zero) go to 300
        end do
      end if   
!
!     Evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      call fcn(n,x,fvec,iflag,iblk,isblk)
      nfev = 1
      if (iflag < 0) go to 300
      fnorm = enorm(n,fvec)
!
!     Initialize iteration counter and monitors.
!
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
!
!     Beginning of the outer loop.
!
   30 continue
         jeval = .true.
!
!        Calculate the jacobian matrix.
!
         if (iopt == 2) then
!
!        Code approximates the jacobian.
!
   31       iflag = 2
            call fdjac2(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,&
                        wa2,iblk,isblk)
            nfev = nfev + min(ml+mu+1,n)
!**debug
            if(nprint==2) then
               write(ifile3,666)                                      
  666 format('----------- jacobian ---------------')
               do i=1,n
                 write(ifile3,667) (fjac(i,j),j=1,n)
  667 format(1x,1p5g15.6)
               end do
            endif                                                     
!**end debug
         else
!
!          User supplies jacobian.
!
!*           call jac(n,x,fvec,fjac,ldfjac,iflag)
!*           njev = njev+1

         end if
   32    if (iflag < 0) go to 300
!
!        Compute the qr factorization of the jacobian.
!
         call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
!
!        On the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
         if (iter == 1) then
           if (mode /= 2) then
             do j = 1, n
                diag(j) = wa2(j)             
                if (wa2(j) == zero) then     
                  diag(j) = one              
                end if                       
             end do                          
           end if
!                                                   
!          On the first iteration, calculate the norm of the scaled x
!          and initialize the step bound delta.     
!                                                   
           do j = 1, n                              
              wa3(j) = diag(j)*x(j)                 
           end do                                   
           xnorm = enorm(n,wa3)                     
           delta = factor*xnorm                     
           if (delta == zero) then                  
             delta = factor                         
           end if                                   
         end if
!
!        Form (q transpose)*fvec and store in qtf.
!
         do i = 1, n
            qtf(i) = fvec(i)
         end do
         do j = 1, n
            if (fjac(j,j) == zero) cycle
            sum = zero
            do i = j, n
               sum = sum + fjac(i,j)*qtf(i)
            end do
            temp = -sum/fjac(j,j)
            do i = j, n
               qtf(i) = qtf(i) + fjac(i,j)*temp
            end do
         end do
!
!        Copy the triangular factor of the qr factorization into r.
!
         sing = .false.
         do j = 1, n
            l = j
            jm1 = j - 1
            if (jm1 >= 1) then
              do i = 1, jm1
                 r(l) = fjac(i,j)         
                 l = l + n - i            
              end do                      
            end if
            r(l) = wa1(j)
            if (wa1(j) == zero) then
              sing = .true.
            end if  
         end do
!
!        Accumulate the orthogonal factor in fjac.
!
         call qform(n,n,fjac,ldfjac,wa1)
!
!        Rescale if necessary.
!
         if (mode /= 2) then
           do j = 1, n
              diag(j) = amax1(diag(j),wa2(j))       
           end do                                   
         end if
!
!        Beginning of the inner loop.
!
  180    continue
!
!           If requested, call fcn to enable printing of iterates.
!
            if (nprint > 0) then
              iflag = 0
              if (mod(iter-1,nprint) == 0) then                    
                 call fcn(n,x,fvec,iflag,iblk,isblk)               
              end if                                               
              if (iflag < 0) go to 300                             
            end if
!
!           Determine the direction p.
!
            call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
!
!           Store the direction p and x + p. calculate the norm of p.
!
            do j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
            end do
            pnorm = enorm(n,wa3)
!
!           On the first iteration, adjust the initial step bound.
!
            if (iter == 1) then
              delta = amin1(delta,pnorm)
            end if
!
!           Evaluate the function at x + p and calculate its norm.
!
            iflag = 1
            iset=0
            call fcn(n,wa2,wa4,iflag,iblk,isblk)
            nfev = nfev + 1
            if (iflag < 0) go to 300
            fnorm1 = enorm(n,wa4)
!
!           Compute the scaled actual reduction.
!
            actred = -one
            if (fnorm1 < fnorm) then
              actred = one - (fnorm1/fnorm)**2
            end if
!
!           Compute the scaled predicted reduction.
!
            l = 1
            do i = 1, n
               sum = zero
               do j = i, n
                  sum = sum + r(l)*wa1(j)
                  l = l + 1
               end do
               wa3(i) = qtf(i) + sum
            end do
            temp = enorm(n,wa3)
            prered = zero
            if (temp < fnorm) then
              prered = one - (temp/fnorm)**2
            end if
!
!           Compute the ratio of the actual to the predicted
!           reduction.
!
            ratio = zero
            if (prered > zero) then
              ratio = actred/prered
            end if
!
!           Update the step bound.
!
            if (ratio >= p1) then
               ncfail = 0
               ncsuc = ncsuc + 1
               if (ratio >= p5 .or. ncsuc > 1) then
                 delta = amax1(delta,pnorm/p5)
               end if
               if (abs(ratio-one) <= p1) then
                 delta = pnorm/p5
               end if
            else
               ncsuc = 0
               ncfail = ncfail + 1
               delta = p5*delta
            end if
!
!           Test for successful iteration.
!
            if (ratio >= p0001) then
!
!             Successful iteration. update x, fvec, and their norms.
!                                                                              
              do j = 1, n                                                      
                 x(j) = wa2(j)                                                 
                 wa2(j) = diag(j)*x(j)                                         
                 fvec(j) = wa4(j)                                              
              end do                                                           
              xnorm = enorm(n,wa2)                                             
              fnorm = fnorm1                                                   
              iter = iter + 1                                                  
            end if
!
!           Determine the progress of the iteration.
!
            nslow1 = nslow1 + 1
            if (actred >= p001) nslow1 = 0
            if (jeval) nslow2 = nslow2 + 1
            if (actred >= p1) nslow2 = 0
!
!           Test for convergence.
!
            if (delta <= xtol*xnorm .or. fnorm == zero) then
              info = 1
            end if
!**debug
            if(nprint==2) then
              xtnorm=xtol*xnorm                                            
              write(ifile3,668) delta,xtnorm                               
  668 format(' delta, xtnorm at convergence test:',1p2g15.6)
            endif                                                          
!**end debug
            if (info /= 0) go to 300
!
!           Tests for termination and stringent tolerances.
!
            if (nfev >= maxfev) info = 2
            if (p1*amax1(p1*delta,pnorm) <= epsmch*xnorm) info = 3
            if (nslow2 == 5) info = 4
            if (nslow1 == 10) info = 5
            if (info /= 0) go to 300
!
!           Criterion for recalculating jacobian.
!
            if (ncfail /= 2) then

!             Calculate the rank one modification to the jacobian
!             and update qtf if necessary.                                
!                                                                         
              do j = 1, n                                                 
                 sum = zero                                               
                 do i = 1, n                                              
                    sum = sum + fjac(i,j)*wa4(i)                          
                 end do                                                   
                 wa2(j) = (sum - wa3(j))/pnorm                            
                 wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)                
                 if (ratio >= p0001) then
                   qtf(j) = sum
                 end if
              end do                                                      
!                                                                         
!             Compute the qr factorization of the updated jacobian.       
!                                                                         
              call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)                      
              call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)                        
              call r1mpyq(1,n,qtf,1,wa2,wa3)                              
!                                                                         
!             End of the inner loop.                                      
!                                                                         
              jeval = .false.                                             
              go to 180                                                   
            end if
!
!        End of the outer loop.
!
         go to 30
  300 continue
!
!     Termination, either normal or user imposed.
!
      if (iflag < 0) info = iflag
      iflag = 0
      if (nprint > 0) call fcn(n,x,fvec,iflag,iblk,isblk)
      if (info < 0) print 1000, time, iblk       ! 3/1/96
      if (info == 0) print 2000, time, iblk
      if (info == 2) print 3000, time, iblk
      if (info == 3) print 4000, time, iblk
      if (info > 4) print 5000, time, iblk
!
1000  format(' time = ',f8.1,' snsq - execuation terminated because',&
                               ' user set iflag negative.')
2000  format(' time = ',f8.1,' snsq - invalid input parameter.')
3000  format(' time = ',f8.1,' snsq - too many function evaluation.')
4000  format(' time = ',f8.1,' snsq - xtol too small.',&
                               '  no further improvement possible.')
5000  format(' time = ',f8.1,' snsq - iteration not making good',&
                               ' progress.','iblk=',i3)   ! 3/1/96
!
      return
      end subroutine snsq1

! **********************************************************************
!
      subroutine snsq2(fcn,jac,iopt,n,x,fvec,fjac,ldfjac,xtol,maxfev,&
                       ml,mu,epsfcn,diag,mode,factor,nprint,info,nfev,&
                       njev,r,lr,qtf,wa1,wa2,wa3,wa4,iblk,isblk)
!
! ----------------------------------------------------------------------
!
!     SNSQ2 : This program is abbreviated version of SNSQ.
!             Only one iteration is performed and no convergence tests
!             are made. Called by BLOCK when the Jacobian matrix for
!             a superblock is being calculated.
!             SNSQ2 calls FDJAC2.
!             Note that the argument of FCN has IBLK and ISBLK.
!
!            The original program was written by K.L. Hiebert
!
!           A minor modification was made by C. R. Hill and C. Park
!           August 13, 1984
!
!           Updated on : Oct. 25, 1984 C.P.
!
!           Updated to Fortran 90.  Dec. 27, 2006 Cheol Park, NIST
!
! **********************************************************************
!
      use modsim_head
      implicit none

      external fcn,jac
      logical                  :: jeval,sing
      integer                  :: iblk,isblk,njev,nfev,info,nprint,mu,&
                                  ml,maxfev,iopt
      integer                  :: mode,n,lr,ldfjac,iflag,j,iset,iter,ncsuc,&
                                  ncfail,nslow1,nslow2,i,l,jm1
      integer,dimension(1)     :: iwa
      real                     :: factor,epsfcn,epsmch,r1mach,xnorm,&
                                  fnorm,enorm,delta
      real                     :: sum,temp,pnorm,fnorm1,actred,prered,&
                                  ratio
      real                     :: xtol
      real                     :: one=1.0,p1=0.1,p5=0.5,p001=1.0e-3,&
                                  p0001=1.0e-4,zero=0.0
      real,dimension(n)        :: x,fvec,diag,qtf,wa1,wa2,wa3,wa4
      real,dimension(lr)       :: r
      real,dimension(ldfjac,n) :: fjac

!     Evaluate the function at the starting ponit
!     and calculate its norm.
!
      iflag = 1
      call fcn(n,x,fvec,iflag,iblk,isblk)
!
!     Calculate the jacobian matrix.
!        code approximates the jacobian.
!
      iflag = 2
      call fdjac2(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,&
                  wa2,iblk,isblk)
!
!     Compute the qr factorization of the jacobian.
!
      call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
!
!     Calculate scaling matrix and step bound.
!
      do j = 1, n
        diag(j) = wa2(j)
        if (wa2(j) == zero) then
          diag(j) = one
        end if
        wa3(j) = diag(j)*x(j)
      end do
      xnorm = enorm(n,wa3)
      delta = factor*xnorm
      if (delta == zero) then
        delta = factor
      end if
!
!     Form (q transpose)*fvec and store in qtf.
!
      do i = 1, n
        qtf(i) = fvec(i)
      end do

      do j = 1, n
        if (fjac(j,j) == zero) cycle
          sum = zero
          do i = j, n
            sum = sum + fjac(i,j)*qtf(i)
          end do
          temp = -sum/fjac(j,j)
          do i = j, n
            qtf(i) = qtf(i) + fjac(i,j)*temp
          end do
      end do
!
!     Copy the triangular factor of the qr factorization into r.
!
      do j = 1, n
        l = j
        jm1 = j - 1
        if (jm1 >= 1) then
          do i = 1, jm1
            r(l) = fjac(i,j)      
            l = l + n - i         
          end do                  
        end if
        r(l) = wa1(j)
        if (wa1(j) == zero) then
          sing = .true.
        end if
      end do
!
!     Accumulate the orthogonal factor in fjac.
!
      call qform(n,n,fjac,ldfjac,wa1)
!
!     Determine the direction p.
!
      call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
!
!     Store the direction p and x + p. calculate the norm of p.
!
      do j = 1, n
        x(j)=x(j)-wa1(j)
      end do
!
      return
      end subroutine snsq2

