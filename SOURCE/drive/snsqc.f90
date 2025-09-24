!***********************************************************************
!
      subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,&
                        wa1,wa2,isblk,iset)
!
! ----------------------------------------------------------------------
!
!     FDJAC1: This is called by SNSQ. ISBLK and ISET were added.
!             Minor but important modifications by Dan Clark, Dec. 1984
!
!***BEGIN PROLOGUE  FDJAC1
!***REFER TO  SNSQ,SNSQE
!              Subroutine FDJAC1
!
!     This subroutine computes a forward-difference approximation
!     to the N by N Jacobian matrix associated with a specified
!     problem of N functions in N VARIABLES. If the Jacobian has
!     a banded form, then function evaluations are saved by only
!     approximating the nonzero terms.
!
!     The subroutine statement is
!
!       SUBROUTINE FDJAC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
!                         WA1,WA2)
!
!     where
!
!       FCN is the name of the user-supplied subroutine which
!         calculates the functions. FCN must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
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
!         The value of IFLAG should not be changed by FCN unless
!         the user wants to terminate execution of FDJAC1.
!         In this case set IFLAG to a negative integer.
!
!       N Is a positive integer input variable set to the number
!         of functions and variables.
!
!       X is an input array of length N.
!
!       FVEC is an input array of length N which must contain the
!         functions evaluated at X.
!
!       FJAC is an output N by N array which contains the
!         approximation to the Jacobian matrix evaluated at X.
!
!       LDFJAC is a positive integer input variable not less than N
!         which specifies the leading dimension of the array FJAC.
!
!       IFLAG is an integer variable which can be used to terminate
!         the execution of FDJAC1. See description of FCN.
!
!       ML is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         Jacobian matrix. If the Jacobian is not banded, set
!         ML to at least N - 1.
!
!       EPSFCN is an input variable used in determining a suitable
!         step length for the forward-difference approximation. This
!         approximation assumes that the relative errors in the
!         functions are of the order of EPSFCN. If EPSFCN is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       MU is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         Jacobian matrix. If the Jacobian is not banded, set
!         MU to at least N - 1.
!
!       WA1 and WA2 are work arrays of length N. If ML + MU + 1 is at
!         least N, then the Jacobian is considered dense, and WA2 is
!         not referenced.
!
!     Subprograms called
!
!       MINPACK-supplied ... R1MACH
!
!       FORTRAN-supplied ... ABS,AMAX1,SQRT
!
!     MINPACK. Version of June 1979.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
!     Updated to Fortran 90.  Dec. 27, 2006 Cheol Park, NIST
!
!     **********
!***routines called  r1mach
!***end prologue  fdjac1
!***********************************************************************

      use modsim_head
      implicit none

      integer                       :: n,ldfjac,iflag,ml,mu
      integer                       :: iset,isblk
      real                          :: epsfcn
      real,dimension(n)             :: x,fvec,wa1,wa2
      real,dimension(ldfjac,n)      :: fjac
      integer                       :: i,j,k,msum
      real                          :: eps,epsmch,h,temp
      real                          :: r1mach
      real                          :: zero=0.0

!
!     First executable statement FDJAC1
!
      epsmch = r1mach(4)
!
      eps = sqrt(amax1(epsfcn,epsmch))
      msum = ml + mu + 1
      if (msum >= n) then
!
!        Computation of dense approximate Jacobian.
!
         do j = 1, n
            temp = x(j)
            h = eps*abs(temp)*isign(issolv(isblk,j))
            if (h == zero) h = eps*isign(issolv(isblk,j))
            x(j) = temp + h
            call fcn(n,x,wa1,iflag,isblk,iset)
            if (iflag < 0) exit
            x(j) = temp
            do i = 1, n
               fjac(i,j) = (wa1(i) - fvec(i))/h
            end do
         end do
      else 
!
!        Computation of banded approximate Jacobian.
!
         do k = 1, msum
            do j = k, n, msum
               wa2(j) = x(j)
               h = eps*abs(wa2(j))
               if (h == zero) h = eps
               x(j) = wa2(j) + h
            end do
            call fcn(n,x,wa1,iflag,isblk,iset)
            if (iflag < 0) exit
            do j = k, n, msum
               x(j) = wa2(j)
               h = eps*abs(wa2(j))
               if (h == zero) h = eps
               do i = 1, n
                  fjac(i,j) = zero
                  if (i >= j - mu .and. i <= j + ml)&
                     fjac(i,j) = (wa1(i) - fvec(i))/h
               end do
            end do
         end do
      end if

      return
!
!     Last statement of FDJAC1
!
      end subroutine fdjac1

!***********************************************************************
!
      subroutine fdjac2(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,&
                        wa1,wa2,iblk,isblk)
!
! ----------------------------------------------------------------------
!
!     FDJAC2: This program is essentially same as FDJAC1 except for
!             change in argument. The change is needed to be called
!             by SNSQ2. See the description in FDJAC1.
!
! **********************************************************************

      use modsim_head
      implicit none

      integer                       :: n,ldfjac,iflag,ml,mu
      integer                       :: iset,isblk,iblk
      real                          :: epsfcn
      real,dimension(n)             :: x,fvec,wa1,wa2
      real,dimension(ldfjac,n)      :: fjac
      integer                       :: i,j,k,msum
      real                          :: eps,epsmch,h,temp
      real                          :: r1mach
      real                          :: zero=0.0

!
!     First executable statement FDJAC2
!
      epsmch = r1mach(4)
!
      eps = sqrt(amax1(epsfcn,epsmch))
      msum = ml + mu + 1
      if (msum >= n) then
!
!        Computation of dense approximate Jacobian.
!
         do j = 1, n
            temp = x(j)
            h = eps*abs(temp)*isign(isolve(iblk,j))
            if (h == zero) h = eps*isign(isolve(iblk,j))
            x(j) = temp + h
            call fcn(n,x,wa1,iflag,iblk,isblk)
            if (iflag < 0) exit
            x(j) = temp
            do i = 1, n
               fjac(i,j) = (wa1(i) - fvec(i))/h
            end do
         end do
      else  
!
!        Computation of banded approximate Jacobian.
!
         do k = 1, msum
            do j = k, n, msum
               wa2(j) = x(j)
               h = eps*abs(wa2(j))
               if (h == zero) h = eps
               x(j) = wa2(j) + h
            end do
            call fcn(n,x,wa1,iflag,iblk,isblk)
            if (iflag < 0) exit
            do j = k, n, msum
               x(j) = wa2(j)
               h = eps*abs(wa2(j))
               if (h == zero) h = eps
               do i = 1, n
                  fjac(i,j) = zero
                  if (i >= j - mu .and. i <= j + ml)&
                     fjac(i,j) = (wa1(i) - fvec(i))/h
               end do
            end do
         end do
      end if

      return
!
!     Last statement of FDJAC2
!
      end subroutine fdjac2

!***********************************************************************
!
      real function r1mach(i)
!
! ----------------------------------------------------------------------
!
!**BEGIN PROLOGUE R1MACH
!**DATE WRITTEN 790101  (YYMMDD)
!**REVISION DATE 801001 (YYMMDD)
!**CATEGORY NO. Q
!**KEYWORDS MACHINE CONSTANTS
!**AUTHOR  FOX, P. A., (BELL LABS)
!          HALL, A. D., (BELL LABS)
!          SCHRYER, N. L., (BELL LABS)
!**PURPOSE  Returns single precision machine dependent constants
!**DESCRIPTION
!    R1MACH can be used to obtain machine-dependent parameters
!    for the local machine environment.  It is a function
!    subroutine with one (input) argument, and can be called
!    as follows, for example
!
!         A = R1MACH(I)
!
!    where I=1,...,5.  The (output) value of A above is
!    determined by the (input) value of I.  The results for
!    various values of I are discussed below.
!
!  Single-Precision Machine Constants
!  R1MACH(1) = B**(EMIN-1), the smallest positive magnitude..
!  R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!  R1MACH(3) = B**(-T), the smallest relative spacing.
!  R1MACH(4) = B**(1-T), the largest relative spacing.
!  R1MACH(5) = LOG10(B)
!
!  where  B = Base
!         T = number of base-B digits in mantissa
!      EMIN = smallest (most negative) exponent
!      EMAX = largest exponent
!
!  For example, the single-precision machine constants of an IBM PC are:
!
!          B = 2
!          T = 23
!       EMIN = 127
!       EMAX = 128
!
!  To alter this function for a particular environment,
!  the desired set of data statements should be activated by
!  removing the C from column 1.
!
!  Where possible, octal or hexadecimal constants have been used
!  to specify the constants exactly which has in some cases
!  required the use of equivalent integer arrays.
!
!**REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L., *FRAMEWORK FOR
!                A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE-
!                MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978,
!                PP. 177-188.
!**ROUTINES CALLED XERROR
!**END PROLOGUE R1MACH
!
!     Updated to Fortran 90.  Dec. 27, 2006 Cheol Park, NIST
!
!***********************************************************************

      integer,dimension(2)     :: small,large,right,diver,log10
      real,dimension(5)        :: rmach

!      INTEGER SMALL(2)
!      INTEGER LARGE(2)
!      INTEGER RIGHT(2)
!      INTEGER DIVER(2)
!      INTEGER LOG10(2)
!      real rmach(5)

      equivalence (rmach(1),small(1))
      equivalence (rmach(2),large(1))
      equivalence (rmach(3),right(1))
      equivalence (rmach(4),diver(1))
      equivalence (rmach(5),log10(1))
!
!     Machine constants for the IBM PC or alike
!
!*      data rmach(1) / z'007fffff'/
!*      data rmach(2) / z'7f7fffff'/
!*      data rmach(3) / z'34000000'/
!*      data rmach(4) / z'34800000'/
!*      data rmach(5) / 0.3/

      data rmach(1) / 1.2e-38/
      data rmach(2) / 3.4e+38/
      data rmach(3) / 0.119209e-06/
      data rmach(4) / 0.238419e-06/
      data rmach(5) / 0.30103/
!
!     Machine constants for the SUN-3
!
!      data small(1) /    8388608 /
!      data large(1) / 2139095039 /
!      data right(1) /  864026624 /
!      data diver(1) /  872415232 /
!      data log10(1) / 1050288283 /
!
!     Machine constants used for the examples given in 
!     NBSIR 85-3243 and NBSIR 86-3331 are as follows:
!
!*    data rmach(1) / 5.3976e-79/
!*    data rmach(2) / 7.2370e+75/
!*    data rmach(3) / 5.9605e-8/
!*    data rmach(4) / 9.5367e-7/
!
!     ***first executable statement r1mach
!
!
      r1mach = rmach(i)
!
      return
      end function r1mach

