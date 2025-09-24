!***********************************************************************
!
      subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
!
!----------------------------------------------------------------------
!
!**begin prologue  dogleg
!**refer to  snsq,snsqe
!    **********
!     Subroutine DOGLEG
!
!     Given an M by N matrix A, an N by N nonsingular DIAGONAL
!     matrix D, an M-vector B, and a positive number DELTA, the
!     problem is to determine the convex combination X of the
!     Gauss-Newton and scaled gradient directions that minimizes
!     (A*X - B) in the least squares sense, subject to the
!     restriction that the Euclidean norm of D*X be at most DELTA.
!
!     This subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     QR factorization of A. That is, if A = Q*R, where Q has
!     orthogonal columns and R is an upper triangular matrix,
!     then DOGLEG expects the full upper triangle of R and
!     the first N components of (Q TRANSPOSE)*B.
!
!     The subroutine statement is
!
!       SUBROUTINE DOGLEG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
!
!     where
!
!       N is a positive integer input variable set to the order of R.
!
!       R is an input array of length LR which must contain the upper
!         triangular matrix R stored by rows.
!
!       LR is a positive integer input variable not less than
!         (N*(N+1))/2.
!
!       DIAG is an input array of length N which must contain the
!         diagonal elements of the matrix D.
!
!       QTB is an input array of length N which must contain the first
!         N elements of the vector (Q TRANSPOSE)*B.
!
!       DELTA is a positive input variable which specifies an upper
!         bound on the Euclidean norm of D*X.
!
!       X is an output array of length N which contains the desired
!         convex combination of the Gauss-Newton direction and the
!         scaled gradient direction.
!
!       WA1 and WA2 are work arrays of length N.
!
!     Subprograms called
!
!       MINPACK-supplied ... R1MACH,ENORM
!
!       FORTRAN-supplied ... ABS,AMAX1,AMIN1,SQRT
!
!     MINPACK. Version of July 1978.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
!     **********
!     Updated to Fortran 90.  Dec. 27, 2006 Cheol Park, NIST
!
!***routines called  enorm,r1mach
!***end prologue  dogleg
!
!***********************************************************************

      implicit none
      integer                :: n,lr
      real                   :: delta
      real,dimension(lr)     :: r
      real,dimension(n)      :: diag,qtb,x,wa1,wa2
      integer                :: i,j,jj,jp1,k,l
      real                   :: alpha,bnorm,epsmch,gnorm,qnorm,&
                                sgnorm,sum,temp
      real                   :: r1mach,enorm
      real                   :: one=1.0,zero=0.0

!
!     First executable statement DOGLEG
!
      epsmch = r1mach(4)
!
!     First, calculate the Gauss-Newton direction.
!
      jj = (n*(n + 1))/2 + 1
      do k = 1, n
         j = n - k + 1
         jp1 = j + 1
         jj = jj - k
         l = jj + 1
         sum = zero
         if (n >= jp1) then
           do i = jp1, n
              sum = sum + r(l)*x(i)      
              l = l + 1                  
           end do                        
         end if
         temp = r(jj)
         if (temp == zero) then
           l = j
           do i = 1, j                                  
              temp = amax1(temp,abs(r(l)))              
              l = l + n - i                             
           end do                                       
           temp = epsmch*temp                           
           if (temp == zero) then
             temp = epsmch
           end if
         end if
         x(j) = (qtb(j) - sum)/temp
      end do
!
!     Test whether the Gauss-Newton direction is acceptable.
!
      do j = 1, n
         wa1(j) = zero
         wa2(j) = diag(j)*x(j)
      end do
      qnorm = enorm(n,wa2)
      if (qnorm <= delta) return
!     The Gauss-Newton direction is not acceptable.
!     Next, calculate the scaled gradient direction.
!
      l = 1
      do j = 1, n
         temp = qtb(j)
         do i = j, n
            wa1(i) = wa1(i) + r(l)*temp
            l = l + 1
         end do
         wa1(j) = wa1(j)/diag(j)
      end do
!
!     Calculate the norm of the scaled gradient direction,
!     Normalize, and rescale the gradient.
!
      gnorm = enorm(n,wa1)
      sgnorm = zero
      alpha = delta/qnorm
      if (gnorm /= zero) then
        do j = 1, n
           wa1(j) = (wa1(j)/gnorm)/diag(j)                                  
        end do                                                              
!                                                                           
!       Calculate the point along the scaled gradient                       
!       at which the quadratic is minimized.                                
!                                                                           
        l = 1                                                               
        do j = 1, n                                                         
           sum = zero                                                       
           do i = j, n                                                      
              sum = sum + r(l)*wa1(i)                                       
              l = l + 1                                                     
           end do                                                           
           wa2(j) = sum                                                     
        end do                                                              
        temp = enorm(n,wa2)                                                 
        sgnorm = (gnorm/temp)/temp                                          
!                                                                           
!       Test whether the scaled gradient direction is acceptable.           
!                                                                           
        alpha = zero                                                        
        if (sgnorm < delta) then
!                                                                           
!         The scaled gradient direction is not acceptable.                  
!         Finally, calculate the point along the dogleg                     
!         at which the quadratic is minimized.                              
!                                                                           
          bnorm = enorm(n,qtb)                                              
          temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)                 
          temp = temp - (delta/qnorm)*(sgnorm/delta)**2&
                 + sqrt((temp-(delta/qnorm))**2&
                 +(one-(delta/qnorm)**2)*(one-(sgnorm/delta)**2))
          alpha = ((delta/qnorm)*(one - (sgnorm/delta)**2))/temp            
        end if
      end if
!
!     Form appropriate convex combination of the Gauss-Newton
!     direction and the scaled gradient direction.
!
      temp = (one - alpha)*amin1(sgnorm,delta)
      do j = 1, n
         x(j) = temp*wa1(j) + alpha*x(j)
      end do

      return
!
!     Last statement of subroutine DOGLEG.
!
      end subroutine dogleg

!***********************************************************************
!
      real function enorm(n,x)
!
! ----------------------------------------------------------------------
!
!***begin prologue  enorm
!***refer to  snls1,snls1e,snsq,snsqe
!     **********
!     Function ENORM
!
!     Given an N-vector X, this function calculates the
!     Euclidean norm of X.
!
!     The Euclidean norm is computed by accumulating the sum of
!     squares in three different sums. The sums of squares for the
!     small and large components are scaled so that no overflows
!     occur. Non-destructive underflows are permitted. Underflows
!     and overflows do not occur in the computation of the unscaled
!     sum of squares for the intermediate components.
!     The definitions of small, intermediate and large components
!     depend on two constants, RDWARF and RGIANT. The main
!     restrictions on these constants are that RDWARF**2 not
!     underflow and RGIANT**2 not overflow. The constants
!     given here are suitable for every known computer.
!
!     The function statement is
!
!       REAL FUNCTION ENORM(N,X)
!
!     where
!
!       N is a positive integer input variable.
!
!       X is an input array of length N.
!
!     Subprograms called
!
!       FORTRAN-supplied ... ABS,SQRT
!
!     MINPACK. Version of October 1979.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
!
!     Minor modification by D.R. Clark, Oct. 1984, to prevent underflow
!
!     Updated to Fortran 90.  Dec. 27, 2006 Cheol Park, NIST
!
!     **********
!***routines called  (none)
!***end prologue  enorm
!***********************************************************************
!
      implicit none
      integer               :: n
      real,dimension(n)     :: x
      integer               :: i
      real                  :: agiant,floatn,s1,s2,s3,xabs,x1max,x3max
      real                  :: one=1.0,zero=0.0
      real                  :: rdwarf=3.834e-20,rgiant=1.304e19

!
!     First executable statement ENORM
!
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant/floatn

      big_loop: do i = 1, n
         xabs = abs(x(i))

         if (xabs <= rdwarf .or. xabs >= agiant) then
            if (xabs > rdwarf) then
!
!              Sum for large components.
!
               if (xabs > x1max) then
                  s1 = one + s1*(x1max/xabs)**2
                  x1max = xabs
               else
                  s1 = s1 + (xabs/x1max)**2
               end if
            else
!
!              Sum for small components.
!
               if (xabs > x3max) then
                  s3 = one + s3*(x3max/xabs)**2        
                  x3max = xabs                         
               else
                  if (xabs /= zero) then               
                    s3 = s3 + (xabs/x3max)**2          
                  end if
               end if
            end if
         else
!
!           Sum for intermediate components.
!
            s2 = s2 + xabs**2
         end if
      end do big_loop

!
!     Calculation of norm.
!
      if (s1 /= zero) then
         enorm = x1max*sqrt(s1+(s2/x1max)/x1max)
         return
      end if

      if (s2 /= zero) then
!***modification : if s2>=x3max, small components are negligible for
!*** all reasonable values of n. comment out four lines:
!***            if (s2 >= x3max)
!***     *         enorm = sqrt(s2*(one+(x3max/s2)*(x3max*s3)))
!***            if (s2 < x3max)
!***     *         enorm = sqrt(x3max*((s2/x3max)+(x3max*s3)))
!***
!*** insert modified lines:
         if (s2 >= x3max) then
             enorm = sqrt(s2)
         else
             enorm = sqrt(x3max*((s2/x3max)+(x3max*s3)))
         endif
!*** end of modification
      else
            enorm = x3max*sqrt(s3)
      end if

      return
!
!     Last statement of FUNCTION ENORM.
!
      end function enorm

!***********************************************************************
!
      subroutine qform(m,n,q,ldq,wa)
!
! ----------------------------------------------------------------------
!
!***begin prologue  qform
!***refer to  snsq,snsqe
!     **********
!     Subroutine QFORM
!
!     This subroutine proceeds from the computed QR factorization of
!     an M by N matrix A to accumulate the M by M orthogonal matrix
!     Q from its factored form.
!
!     The subroutine statement is
!
!       SUBROUTINE QFORM(M,N,Q,LDQ,WA)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of A and the order of Q.
!
!       N is a positive integer input variable set to the number
!         of columns of A.
!
!       Q is an M by M array. On input the full lower trapezoid in
!         the first min(M,N) columns of Q contains the factored form.
!         On output Q has been accumulated into a square matrix.
!
!       LDQ is a positive integer input variable not less than M
!         which specifies the leading dimension of the array Q.
!
!       WA is a work array of length M.
!
!     Subprograms called
!
!       FORTRAN-supplied ... MIN0
!
!     MINPACK. Version of January 1979.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
!     Updated to Fortran 90.  Dec. 27, 2006 Cheol Park, NIST
!
!     **********
!***routines called  (none)
!***end prologue  qform
!***********************************************************************
!
      implicit none
      integer               :: m,n,ldq
      real,dimension(ldq,m) :: q
      real,dimension(m)     :: wa
      integer               :: i,j,jm1,k,l,minmn,np1
      real                  :: one=1.0,sum,temp,zero=0.0

!
!     First executable statement QFORM
!
      minmn = min(m,n)
      if (minmn >= 2) then
        do j = 2, minmn
           jm1 = j - 1                
           do i = 1, jm1              
              q(i,j) = zero           
           end do                     
        end do                        
      end if
!
!     Initialize remaining columns to those of the identity matrix.
!
      np1 = n + 1
      if (m >= np1) then
        do j = np1, m
           do i = 1, m                
              q(i,j) = zero           
           end do                     
           q(j,j) = one               
        end do                        
      end if
!
!     Accumulate Q from its factored form.
!
      do l = 1, minmn
         k = minmn - l + 1
         do i = k, m
            wa(i) = q(i,k)
            q(i,k) = zero
         end do
         q(k,k) = one
         if (wa(k) /= zero) then
           do j = k, m
              sum = zero                                  
              do i = k, m                                 
                 sum = sum + q(i,j)*wa(i)                 
              end do                                      
              temp = sum/wa(k)                            
              do i = k, m                                 
                 q(i,j) = q(i,j) - temp*wa(i)             
              end do                                      
           end do                                         
         end if
      end do
      return
!
!     Last statement of QFORM.
!
      end subroutine qform

!***********************************************************************
!
      subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,sigma,acnorm,wa)
!
! ----------------------------------------------------------------------
!
!***begin prologue  qrfac
!***refer to  snls1,snls1e,snsq,snsqe
!***routines called  enorm,r1mach
!***description
!     Subroutine QRFAC
!
!     This subroutine uses Householder transformations with column
!     pivoting (optional) to compute a QR factorization of the
!     M by N matrix A. That is, QRFAC determines an orthogonal
!     matrix Q, a permutation matrix P, and an upper trapezoidal
!     matrix R with diagonal elements of nonincreasing magnitude,
!     such that A*P = Q*R. The Householder transformation for
!     column K, K = 1,2,...,MIN(M,N), is of the form
!
!                           T
!           I - (1/U(K))*U*U
!
!     where U has zeros in the first K-1 positions. The form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding LINPACK subroutine.
!
!     The subroutine statement is
!
!       SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of A.
!
!       N is a positive integer input variable set to the number
!         of columns of A.
!
!       A is an M by N array. On input A contains the matrix for
!         which the QR factorization is to be computed. On output
!         the strict upper trapezoidal part of A contains the strict
!         upper trapezoidal part of R, and the lower trapezoidal
!         part of A contains a factored form of Q (the non-trivial
!         elements of the U vectors described above).
!
!       LDA is a positive integer input variable not less than M
!         which specifies the leading dimension of the array A.
!
!       PIVOT is a logical input variable. If pivot is set .TRUE.,
!         then column pivoting is enforced. If pivot is set .FALSE.,
!         then no column pivoting is done.
!
!       IPVT is an integer output array of length LIPVT. IPVT
!         defines the permutation matrix P such that A*P = Q*R.
!         Column J of P is column IPVT(J) of the identity matrix.
!         If pivot is .FALSE., IPVT is not referenced.
!
!       LIPVT is a positive integer input variable. If PIVOT is
!             .FALSE., then LIPVT may be as small as 1. If PIVOT is
!             .TRUE., then LIPVT must be at least N.
!
!       SIGMA is an output array of length N which contains the
!         diagonal elements of R.
!
!       ACNORM is an output array of length N which contains the
!         norms of the corresponding columns of the input matrix A.
!         If this information is not needed, then ACNORM can coincide
!         with SIGMA.
!
!       WA is a work array of length N. If pivot is .FALSE., then WA
!         can coincide with SIGMA.
!
!     Subprograms called
!
!       MINPACK-Supplied ... R1MACH,ENORM
!       FORTRAN-Supplied ... AMAX1,SQRT,MIN0
!
!     MINPACK. Version of December 1978.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
!     Updated to Fortran 90.  Dec. 27, 2006 Cheol Park, NIST
!
!***end prologue  qrfac
!***********************************************************************
!
      implicit none
      integer                 :: m,n,lda,lipvt
      integer,dimension(lipvt):: ipvt
      logical                 :: pivot
      real,dimension(lda,n)   :: a
      real,dimension(n)       :: sigma,acnorm,wa
      integer                 :: i,j,jp1,k,kmax,minmn
      real                    :: ajnorm,epsmch,sum,temp
      real                    :: r1mach,enorm
      real                    :: one=1.0,p05=5.0e-2,zero=0.0

!
!     First executable statement QRFAC
!
      epsmch = r1mach(4)
!
!     Compute the initial column norms and initialize several arrays.
!
      do j = 1, n
         acnorm(j) = enorm(m,a(1,j))
         sigma(j) = acnorm(j)
         wa(j) = sigma(j)
         if (pivot) then
           ipvt(j) = j
         end if
      end do
!
!     Reduce A to R with Householder transformations.
!
      minmn = min(m,n)
      big_loop: do j = 1, minmn
         if (pivot) then
!
!          Bring the column of largest norm into the pivot position.
!                                                                       
           kmax = j                                                     
           do k = j, n                                                  
              if (sigma(k) > sigma(kmax)) kmax = k                      
           end do                                                       
           if (kmax /= j) then
             do i = 1, m                                                
                temp = a(i,j)                                           
                a(i,j) = a(i,kmax)                                      
                a(i,kmax) = temp                                        
             end do                                                     
             sigma(kmax) = sigma(j)                                     
             wa(kmax) = wa(j)                                           
             k = ipvt(j)                                                
             ipvt(j) = ipvt(kmax)                                       
             ipvt(kmax) = k                                             
           end if
         end if
!
!        Compute the Householder transformation to reduce the
!        J-th column of A to A multiple of the J-th unit vector.
!
         ajnorm = enorm(m-j+1,a(j,j))
         if (ajnorm /= zero) then
           if (a(j,j) < zero) then
             ajnorm = -ajnorm
           end if
           do i = j, m                                                        
              a(i,j) = a(i,j)/ajnorm                                          
           end do                                                             
           a(j,j) = a(j,j) + one                                              
!                                                                             
!          Apply the transformation to the remaining columns                  
!          and update teh norms.                                              
!                                                                             
           jp1 = j + 1                                                        
           if (n >= jp1) then
             do k = jp1, n                                                    
                sum = zero                                                    
                do i = j, m                                                   
                   sum = sum + a(i,j)*a(i,k)                                  
                end do                                                        
                temp = sum/a(j,j)                                             
                do i = j, m                                                   
                   a(i,k) = a(i,k) - temp*a(i,j)                              
                end do                                                        
                if (pivot .and. sigma(k) /= zero) then
                  temp = a(j,k)/sigma(k)
                  sigma(k) = sigma(k)*sqrt(amax1(zero,one-temp**2))           
                  if (p05*(sigma(k)/wa(k))**2 <= epsmch) then
                    sigma(k) = enorm(m-j,a(jp1,k))                            
                    wa(k) = sigma(k)                                          
                  end if                                                      
                end if
             end do
           end if
         end if
         sigma(j) = -ajnorm
      end do big_loop
      return
!
!     Last statement of QRFAC.
!
      end subroutine qrfac

!***********************************************************************
!
      subroutine r1mpyq(m,n,a,lda,v,w)
!
! ----------------------------------------------------------------------
!
!***begin prologue  r1mpyq
!***refer to  snsq,snsqe
!     **********
!     Subroutine R1MPYQ
!
!     Given an M by N matrix A, this subroutine computes A*Q where
!     Q is the product of 2*(N - 1) transformations
!
!           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
!
!     and GV(I), GW(I) are Givens rotations in the (I,N) plane which
!     eliminate elements in the I-th and N-th planes, respectively.
!     Q itself is not given, rather the information to recover the
!     GV, GW rotations is supplied.
!
!     The subroutine statement is
!
!       SUBROUTINE R1MPYQ(M,N,A,LDA,V,W)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of A.
!
!       N is a positive integer input variable set to the number
!         of columns of A.
!
!       A is an M by N ARRAY. On input A must contain the matrix
!         to be postmultiplied by the orthogonal matrix Q
!         described above. On output A*Q has replaced A.
!
!       LDA is a positive integer input variable not less than M
!         which specifies the leading dimension of the array A.
!
!       V is an input array of length N. V(I) must contain the
!         information necessary to recover the Givens rotation GV(I)
!         described above.
!
!       W is an input array of length N. W(I) must contain the
!         information necessary to recover the Givens rotation GW(I)
!         described above.
!
!     Subroutines called
!
!       FORTRAN-Supplied ... abs,sqrt
!
!     MINPACK. Version of December 1978.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
!     Updated to Fortran 90.  Dec. 27, 2006 Cheol Park, NIST
!
!     **********
!***routines called  (none)
!***end prologue  r1mpyq
!***********************************************************************

      implicit none
      integer                :: m,n,lda
      real,dimension(lda,n)  :: a
      real,dimension(n)      :: v,w
      integer                :: i,j,nmj,nm1
      real                   :: cos,sin,temp
      real                   :: one=1.0

!
!     First executable statement R1MPYQ
!
      nm1 = n - 1
      if (nm1 < 1) return
      do nmj = 1, nm1
         j = n - nmj
         if (abs(v(j)) > one) then
           cos = one/v(j)
           sin = sqrt(one-cos**2)
         else
           sin = v(j)
           cos = sqrt(one-sin**2)
         endif
         do i = 1, m
            temp = cos*a(i,j) - sin*a(i,n)
            a(i,n) = sin*a(i,j) + cos*a(i,n)
            a(i,j) = temp
         end do
      end do
!
!     Apply the second set of Givens rotations to A
!
      do j = 1, nm1
         if (abs(w(j)) > one) then
           cos = one/w(j)
           sin = sqrt(one-cos**2)
         else
           sin = w(j)
           cos = sqrt(one-sin**2)
         endif
         do i = 1, m
            temp = cos*a(i,j) + sin*a(i,n)
            a(i,n) = -sin*a(i,j) + cos*a(i,n)
            a(i,j) = temp
         end do
      end do

      return
!
!     Last statement of R1MPYQ.
!
      end subroutine r1mpyq

!***********************************************************************
!
      subroutine r1updt(m,n,s,ls,u,v,w,sing)
!
! ----------------------------------------------------------------------
!
!***begin prologue  r1updt
!***refer to  snsq,snsqe
!     **********
!     Subroutine R1UPDT
!
!     Given an M by N lower trapezoidal matrix S, an M-vector U,
!     and an N-vector V, the problem is to determine an
!     orthogonal matrix Q such that
!
!                   T
!           (S + U*V )*Q
!
!     is again lower trapezoidal.
!
!     This subroutine determines Q as the product of 2*(N - 1)
!     transformations
!
!           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
!
!     where GV(I), GW(I) are Givens rotations in the (I,N) plane
!     which eliminate elements in the I-th and N-th planes,
!     respectively. Q Itself is not accumulated, rather the
!     information to recover the GV, GW rotations is returned.
!
!     The subroutine statement is
!
!       SUBROUTINE R1UPDT(M,N,S,LS,U,V,W,SING)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of S.
!
!       N is a positive integer input variable set to the number
!         of columns of S. N must not exceed M.
!
!       S is an array of length LS. On input S must contain the lower
!         trapezoidal matrix S stored by columns. On output S contains
!         the lower trapezoidal matrix produced as described above.
!
!       LS is a positive integer input variable not less than
!         (N*(2*M-N+1))/2.
!
!       U is an input array of length M which must contain the
!         vector U.
!
!       V is an array of length N. On input V must contain the vector
!         V. On output V(I) contains the information necessary to
!         recover the Givens rotation GV(I) described above.
!
!       W is an output array of length M. W(I) contains information
!         necessary to recover the Givens rotation GW(I) described
!         above.
!
!       SING is a logical output variable. SING is set .TRUE. if any
!         of the diagonal elements of the output S are zero. Otherwise
!         SING is set .FALSE.
!
!     Subprograms called
!
!       MINPACK-supplied ... R1MACH
!       FORTRAN-supplied ... ABS,SQRT
!
!     MINPACK. Version of December 1978.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More,
!     John L. Nazareth
!     Updated to Fortran 90.  Dec. 27, 2006 Cheol Park, NIST
!
!     **********
!***routines called  r1mach
!***end prologue  r1updt
!***********************************************************************

      implicit none
      integer                 :: m,n,ls
      logical                 :: sing
      real,dimension(ls)      :: s
      real,dimension(m)       :: u,w
      real,dimension(n)       :: v
      integer                 :: i,j,jj,l,nmj,nm1
      real                    :: cos,cotan,giant,sin,tan,tau,temp
      real                    :: r1mach
      real                    :: one=1.0,p5=5.0e-1,p25=2.5e-1,zero=0.0

!
!     First executable statement R1UPDT
!
      giant = r1mach(2)
!
!     Initialize the diagonal element pointer.
!
      jj = (n*(2*m - n + 1))/2 - (m - n)
!
!     Move the nontrivial part of the last column of S into W.
!
      l = jj
      do i = n, m
         w(i) = s(l)
         l = l + 1
      end do
!
!     Rotate the vector V into a multiple of the N-th unit vector
!     in such a way that a spike is introduced into W.
!
      nm1 = n - 1
      if (nm1 >= 1) then
        do nmj = 1, nm1
           j = n - nmj                                                    
           jj = jj - (m - j + 1)                                          
           w(j) = zero                                                    
           if (v(j) /= zero) then
!                                                                         
!            Determine a Givens rotation which eliminates the             
!            J-th element of V.                                           
!                                                                         
             if (abs(v(n)) < abs(v(j))) then
                cotan = v(n)/v(j)                                         
                sin = p5/sqrt(p25+p25*cotan**2)                           
                cos = sin*cotan                                           
                tau = one                                                 
                if (abs(cos)*giant > one) tau = one/cos                   
             else
                tan = v(j)/v(n)                                           
                cos = p5/sqrt(p25+p25*tan**2)                             
                sin = cos*tan                                             
                tau = sin                                                 
             end if
!                                                                         
!            Apply the transformation to V and store the information      
!            necessary to recover the Givens rotation.                    
!                                                                         
             v(n) = sin*v(j) + cos*v(n)                                   
             v(j) = tau                                                   
!                                                                         
!            Apply the transformation to S and extend the spike in W.     
!                                                                         
             l = jj                                                       
             do i = j, m                                                  
                temp = cos*s(l) - sin*w(i)                                
                w(i) = sin*s(l) + cos*w(i)                                
                s(l) = temp                                               
                l = l + 1                                                 
             end do                                                       
           end if
        end do                                                            
      end if
!
!     Add the spike from the rank 1 update to W
!
      do i = 1, m
         w(i) = w(i) + v(n)*u(i)
      end do
!
!     Eliminate the spike.
!
      sing = .false.
      if (nm1 >= 1) then
        do j = 1, nm1
           if (w(j) /= zero) then
!                                                                             
!            Determine a Givens rotation which eliminates the                 
!            J-th element of the spike.                                       
!                                                                             
             if (abs(s(jj)) < abs(w(j))) then
                cotan = s(jj)/w(j)                                            
                sin = p5/sqrt(p25+p25*cotan**2)                               
                cos = sin*cotan                                               
                tau = one                                                     
                if (abs(cos)*giant > one) tau = one/cos                       
             else
                tan = w(j)/s(jj)                                              
                cos = p5/sqrt(p25+p25*tan**2)                                 
                sin = cos*tan                                                 
                tau = sin                                                     
             end if
!                                                                             
!            Apply the transformation to S and reduce the spike in W.         
!                                                                             
             l = jj                                                           
             do i = j, m                                                      
                temp = cos*s(l) + sin*w(i)                                    
                w(i) = -sin*s(l) + cos*w(i)                                   
                s(l) = temp                                                   
                l = l + 1                                                     
             end do                                                           
!                                                                             
!            Store the information necessary to recover the                   
!            Givens rotation.                                                 
!                                                                             
             w(j) = tau                                                       
           end if
!                                                                             
!          Test for zero diagonal elements in the output S.                   
!                                                                             
           if (s(jj) == zero) sing = .true.                                   
           jj = jj + (m - j + 1)                                              
        end do
      end if
!
!     Move W back into the last column of the output S.
!
      l = jj
      do i = n, m
         s(l) = w(i)
         l = l + 1
      end do
      if (s(jj) == zero) then
        sing = .true.
      end if
      return
!
!     Last statement of R1UPDT.
!
      end subroutine r1updt

