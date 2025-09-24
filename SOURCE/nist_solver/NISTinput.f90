!**************************************************************************
	
	module nist_solver

!   Residual tolerance
	real    :: user_ftol
!   Slover (0:snsq,1:nitsol)
	integer :: solver
!   NITSOL inputs
	integer :: nitinput(10)

	end module nist_solver

!**************************************************************************
	subroutine NISTinput
!--------------------------------------------------------------------------
!   Purpose: Read inputs for nitsol
!--------------------------------------------------------------------------
	use nist_solver
	logical res
    
    ! MAG 211109 add code to check if file exists, if not use nist solver
    inquire(file='solver.par', exist=res)
    if(res == true) then
!	Solver selection and tolerance
		open (unit=102,file='solver.par',status='old',action='read')
		read (102,*) solver
		read (102,*) user_ftol
		close (102)
    else
        solver = 0	!0:SNSQ, 1:NITSOL
        user_ftol = 0;!Residual
    end if
    
!	NITSOL inputs
    inquire(file='nitinput.par', exist=res)
    if(res == true) then
		open (unit=103,file='nitinput.par',status='old',action='read')	
		read (103,*) nitinput(1:10)
		close (103)	
    else 
        nitinput = (/0,   0,   0,   0,   1,   0,   0,   0,   0,   0/)
    end if


    end subroutine NISTinput

! NOTES:
!    c Optional every-user input:
!c
!
!c    input(1) = nnimax = maximum number of nonlinear iterations (default 200).
!c 
!
!c    input(2) = ijacv = flag for determining the method of J*v evaluation:
!
!c                 0 => finite-difference evaluation (default) 
!
!c                 1 => analytic evaluation
!c
!c    input(3) = ikrysl = flag for determining the Krylov solver: 
!
!c                 0 => GMRES (default)
!
!c                 1 => BiCGSTAB
!
!c                 2 => TFQMR
!
!c
!
!c               For brief descriptions of the solvers plus references, 
!
!c               see the subroutines nitgm, nitstb, and nittfq. 
!
!c
!
!c    input(4) = kdmax = maximum Krylov subspace dimension when GMRES is used 
!
!c               (default 20). 
!
!c
!
!c    input(5) = irpre = flag for right preconditioning: 
!
!c                 0 => no right preconditioning
!
!c                 1 => right preconditioning
!c
!c Optional experienced user input:
!
!c
!
!c    input(6) = iksmax = maximum allowable number of iterations per call 
!
!c               to the Krylov solver routine (default 1000). 
!
!c
!
!c    input(7) = iresup = residual update flag when GMRES is used; on 
!
!c               restarts, the residual is updated as follows: 
!
!c                 0 => linear combination (default) 
!
!c                 1 => direct evaluation
!
!c               The first is cheap (one n-vector saxpy) but may lose 
!
!c               accuracy with extreme residual reduction; the second 
!
!c               retains accuracy better but costs one J*v product per 
!
!c               restart. 
!
!c
!
!c    input(8) = ifdord = order of the finite-difference formula (sometimes) 
!
!c               used when input(2) = ijacv = 0. When input(2) = ijacv = 0, 
!
!c               this must be 0, 1, 2, or 4 on input; otherwise, it is 
!
!c               irrelevant. With input(2) = ijacv = 0, the precise 
!
!c               meaning is as follows: 
!
!c
!
!c               If GMRES is used, then ifdord matters only if input(7) = 
!
!c               iresup = 1, in which case it determines the order of 
!
!c               the finite-difference formula used in evaluating the 
!
!c               initial residual at each GMRES restart (default 2); if 
!
!c               ifdord = 0 on input, then it is set to 2 below. NOTE: This 
!
!c               only affects initial residuals at restarts; first-order 
!
!c               differences are always used within each GMRES cycle. Using 
!
!c               higher-order differences at restarts only should give 
!
!c               the same accuracy as if higher-order differences were 
!
!c               used throughout; see K. Turner and H. F. Walker, "Efficient 
!
!c               high accuracy solutions with GMRES(m)," SIAM J. Sci. 
!
!c               Stat. Comput., 13 (1992), pp. 815--825. 
!
!c               
!
!c               If BiCGSTAB or TFQMR is used, then ifdord determines the 
!
!c               order of the finite-difference formula used at each 
!
!c               iteration (default 1); if ifdord = 0 on input, then it 
!
!c               is set to 1 below. 
!
!c
!
!c    input(9) = ibtmax = maximum allowable number of backtracks (step 
!
!c               reductions) per call to nitbt (default 10). 
!
!c
!
!c               USAGE NOTE: Backtracking can be turned off by setting 
!
!c		ibtmax = -1. Other negative values of ibtmax are not 
!
!c               valid. 
!
!c
!
!c    input(10) = ieta = flag determining the forcing term eta as follows: 
!
!c                 0 => abs( ||fcur|| - ||fprev+Jprev*sprev|| )/||fprev||
!
!c                      (default) 
!
!c                 1 => (||fcur||/||fprev||)**2 
!
!c                 2 => gamma*(||fcur||/||fprev||)**alpha 
!
!c                      for user-supplied gamma in (0,1] and alpha in (1,2] 
!
!c                 3 => fixed (constant) eta in (0,1), either 0.1 (default) 
!
!c		       or specified by the user (see USAGE NOTE below) 
!
!c               Here, fcur = current f, fprev = previous f, etc. The Krylov 
!
!c               iterations are terminated when an iterate s satisfies 
!
!c               an inexact Newton condition ||F + J*s|| .le. eta*||F||.
!
!c
!
!c               USAGE NOTE: If input(10) = ieta = 2, then alpha and gamma 
!
!c               must be set in common block nitparam.h as described below. 
!
!c		If input(10) = ieta = 3, then the desired constant eta may 
!
!c		be similarly set in nitparam.h if a value other than the 
!
!c		default of 0.1 is desired. 
!**************************************************************************