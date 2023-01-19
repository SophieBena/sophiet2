/** 
 Algorithm NNLS (Non-negative least-squares)
 Given an m by n matrix A, and an m-vector B, computes an n-vector X, that solves the least squares problem A * X = B, subject to X>=0

 @param a
 On entry, a[ 0... N ][ 0 ... M ] contains the M by N matrix A.
 On exit, a[][] contains the product matrix Q*A, where Q is an m by n orthogonal matrix generated implicitly by this function.

 @param m
 Matrix dimension

 @param n
 Matrix dimension

 @param b
 On entry, b[] must contain the m-vector B.
 On exit, b[] contains Q*B

 @param x
 On exit, x[] will contain the solution vector

 @param rnorm
 On exit, rnorm contains the Euclidean norm of the residual vector.

 @param w
 An n-array of working space, wp[].
 On exit, wp[] will contain the dual solution vector. wp[i]=0.0 for all i in set p and wp[i]<=0.0 for all i in set z.

 @param zz
 An m-array of working space, zz[].

 @param index
 An n-array of working space, index[].

 @return
 True if succesful. False if iteration count exceeded 3*N or initial conditions are invalid.
 */
int nnlsWithConstraintMatrix(double **a, int m, int n, double *b, double *x, double *rnorm,
    double *w, double *zz, int *index, int debug);
