	subroutine jacv(n, xcur, fcur, ijob, v, z, rpar, ipar, itrmjv,
     &				nfe,iblk,isblk)
c ------------------------------------------------------------------------
c
c user-supplied subroutine for evaluating J*v or P(inverse)*v, where J is 
c the Jacobian of f and P is a right preconditioning operator. If neither 
c analytic J*v evaluations nor right preconditioning is used, this can be 
c a dummy subroutine; if right preconditioning is used but not analytic
c J*v evaluations, this need only evaluate P(inverse)*v.
c
c ------------------------------------------------------------------------
c Explanation: 
c
c  n       = dimension of the problem.
c
c  xcur    = vector of length n, initial guess on input and final 
c            approximate solution on output. 
c
c  fcur    = vector of length n, value of f at xcur. 
c
c
c
c  rpar    = real parameter/work array passed to the f and jacv routines. 
c
c  ipar    = integer parameter/work array passed to the f and jacv routines. 
c
c
c  ijob    = ijob is an integer flag indicating which product is desired
c              0 => z = J*v
c              1 => z = P(inverse)*v 
c
c  nfe     = number of function evaluations.
c
c  v       = vector to be multiplied. 
c
c  z       = desired product. 
c
c  itrmjv  = termination flag; values suppose have the following meanings: 
c              0 => normal termination; desired product evaluated. 
c              2 => failure to produce P(inverse)*v. 
c			Note: In the current implementation, itrmjv is not active
c
c ------------------------------------------------------------------------
	use modsim_head
	use nist_solver
	implicit none
	external fnc
	integer i,j
	integer n,ijob,ipar(*),itrmjv
	real xcur(n),fcur(n),v(n),z(n),rpar(n)
	integer nfe
	real m(n,n),minv(n,n)

      integer iflag
      real epsfcn
      real x(n),fvec(n),fjac(n,n),wa1(n),wa2(n)
	integer iblk,isblk

	real timep(10),minvp(500,500,10)
	integer iblkp,isblkp,itrmjvp
	logical compute_pc
	integer pc_status(10)
	common /jacvsaved/ iblkp,isblkp,timep,minvp,itrmjvp,
     &		compute_pc,pc_status

	integer iterm_rec(10),itrmks_rec(10)
	common /info_saved/ iterm_rec,itrmks_rec

	itrmjv=0
!-------------------------------------------------------------------------
! initialize
!-------------------------------------------------------------------------
	x(1:n)=xcur(1:n)
	fvec(1:n)=fcur(1:n)
	epsfcn=1e-5
!-------------------------------------------------------------------------
	compute_pc=.false.

		if (itrmks_rec(iblk)>=3 .or. 
     &		iterm_rec(iblk)>=5 .or. iterm_rec(iblk)==1) then
			compute_pc=.true. 
			pc_status(iblk)=1	
			iterm_rec(iblk)=0
			itrmks_rec(iblk)=0
		end if
			
!-------------------------------------------------------------------------
! determine whether to update preconditioner
!-------------------------------------------------------------------------
	if (.not. compute_pc) then
		itrmjv=itrmjvp	!if first call fail, will activate z=v
		if (itrmjv>0) goto 100	
		if (pc_status(iblk)==1) then		
			minv(1:n,1:n)=minvp(1:n,1:n,iblk)
		end if
		goto 100
	else
!-------------------------------------------------------------------------
! first call at every time step for ever block
!-------------------------------------------------------------------------
	itrmjv=0
	if (ijob==1) then !right preconditioning
!-------------------------------------------------------------------------	
! use inverse of jacobian directly
!-------------------------------------------------------------------------
		iflag=2
          call fdjac2(fnc,n,x,fvec,fjac,n,iflag,n,n,epsfcn,wa1,
     &                wa2,iblk,isblk)
		nfe=nfe+n
		m(1:n,1:n)=fjac(1:n,1:n)
		call inverse(m,minv,n)
	endif

	minvp(1:n,1:n,iblk)=minv(1:n,1:n)

30	continue

	itrmjvp=itrmjv	! record weather first call success

	endif


100   continue

	if (pc_status(iblk) .ne. 1) then
		z=v
		goto 110
	endif

	if (itrmjv==0) then
		z=matmul(minv,v)
	else
		z=v
	endif

110   continue

	timep(iblk)=time

	return
	end

!*********************************************************************

	subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
	implicit none 
	integer n
	real a(n,n), c(n,n)
	real L(n,n), U(n,n), b(n), d(n), x(n)
	real coeff
	integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
	L=0.0
	U=0.0
	b=0.0

! step 1: forward elimination
	do k=1, n-1
		do i=k+1,n
			coeff=a(i,k)/a(k,k)
			L(i,k) = coeff
			do j=k+1,n
				a(i,j) = a(i,j)-coeff*a(k,j)
			end do
		end do
	end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
	do i=1,n
		L(i,i) = 1.0
	end do
! U matrix is the upper triangular part of A
	do j=1,n
		do i=1,j
			U(i,j) = a(i,j)
		end do
	end do

! Step 3: compute columns of the inverse matrix C
	do k=1,n
		b(k)=1.0
		d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
		do i=2,n
			d(i)=b(i)
			do j=1,i-1
				d(i) = d(i) - L(i,j)*d(j)
			end do
		end do
! Step 3b: Solve Ux=d using the back substitution
		x(n)=d(n)/U(n,n)
		do i = n-1,1,-1
			x(i) = d(i)
			do j=n,i+1,-1
				x(i)=x(i)-U(i,j)*x(j)
			end do
			x(i) = x(i)/u(i,i)
		end do
! Step 3c: fill the solutions x(n) into column k of C
		do i=1,n
			c(i,k) = x(i)
		end do
		b(k)=0.0
	end do
	end subroutine inverse