subroutine intel_direct_solver(n,nnz,row,col,mat,x,rhs)
implicit none
include 'mkl_rci.fi'

integer, intent(in) :: n, nnz
integer, intent(in) :: row(n+1), col(nnz)
real *8, intent(in) :: mat(nnz), rhs(n)
real *8, intent(inout) :: x(n)

integer pt(64)

      INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
      INTEGER iparm(64)
      INTEGER i, idum
      REAL*8 waltime1, waltime2, ddum
      real *8 x2(n)
      
      !DATA nrhs /1/, maxfct /1/, mnum /1/

      nrhs=1
      maxfct=1
      mnum=1

      iparm(:)=0
      iparm(1) = 1 ! no solver default
      iparm(2) = 2 ! fill-in reordering from METIS
      iparm(3) = 1 ! numbers of processors
      iparm(4) = 0 ! no iterative-direct algorithm
      iparm(5) = 0 ! no user fill-in reducing permutation
      iparm(6) = 0 ! =0 solution on the first n compoments of x
      iparm(7) = 0 ! not in use
      iparm(8) = 9 ! numbers of iterative refinement steps
      iparm(9) = 0 ! not in use
      iparm(10) = 13 ! perturbe the pivot elements with 1E-13
      iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      iparm(12) = 0 ! not in use
      iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
      iparm(14) = 0 ! Output: number of perturbed pivots
      iparm(15) = 0 ! not in use
      iparm(16) = 0 ! not in use
      iparm(17) = 0 ! not in use
      iparm(18) = -1 ! Output: number of nonzeros in the factor LU
      iparm(19) = -1 ! Output: Mflops for LU factorization
      iparm(20) = 0 ! Output: Numbers of CG Iterations
      error = 0 ! initialize error flag
      msglvl = 0 ! print statistical information
      mtype = 11 ! real unsymmetric

      pt(:)=0
      

      phase = 11 ! only reordering and symbolic factorization
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, mat, row, col, &
      idum, nrhs, iparm, msglvl, ddum, ddum, error)
      !WRITE(*,*) 'Reordering completed ... '
      IF (error .NE. 0) THEN
         WRITE(*,*) 'E R R O R: in reordering, error =', error
         STOP
      END IF
      !WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
      !WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
      
      phase = 22 ! only factorization
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, mat, row, col, &
      idum, nrhs, iparm, msglvl, ddum, ddum, error)
      !WRITE(*,*) 'Factorization completed ... '
      IF (error .NE. 0) THEN
         WRITE(*,*) 'E R R O R: in factorization, error =', error
         STOP 
      ENDIF
      
      iparm(8) = 2 ! max numbers of iterative refinement steps
      phase = 33 ! only factorization
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, mat, row, col, &
      idum, nrhs, iparm, msglvl, rhs, x, error)
      !WRITE(*,*) 'Solve completed ... '
      !WRITE(*,*) 'The solution of the system is '

      phase = -1 ! release internal memory
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
     idum, nrhs, iparm, msglvl, ddum, ddum, error)
     
      call mkl_dcsrgemv('n', n, mat, row, col, x, x2)
      !print *, 'norm2(residue)=', norm2(x2-rhs)/norm2(rhs)
     
     return
     
     end subroutine intel_direct_solver

