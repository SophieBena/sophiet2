#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#undef MAX
#undef MIN

#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define MIN(a,b)  ((a) < (b) ? (a) : (b))

#define SQUARE(x) ((x)*(x))

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/**
 Householder decomposition of a vector and application of the resulting transformation to column or row vectors of a matrix.

 Algorithm H12 (Lawson and Hanson, page 57). See documentation below.

 @param mode
   Calculation mode, i.e. algorithm H1 or H2

 @param p
   index of the pivot element of the vector stored in v

 @param l
   range of indices to zero elements of vector stored in v

 @param m
   TBD

 @param v
   matrix containing the vector v stored as v[][vCol]

 @param vCol
   index number of column of v that contains vector

 @param h
   extra storage for the transformed pth element of v (u[p])

 @param c
  matrix containing the set of column vectors to apply the transformation, i.e. c[][cColStart] -> c[][cColEnd]

 @param cColStart
   starting column

 @param cColEnd
   ending column

 Algorithm H12 (page 57)
 Householder decomposition of a vector and application of the resulting
 transformation to column or row vectors of a matrix.

 If mode = HOUSEHOLDER_CALC_MODE_H1

 Householder decomposition of a vector v into:

 -        -
 |  v[0]  |
 |    .   |
 | v[p-1] |
 |  y[p]  |
 Q v  =  | v[p+1] |  =  y (stored in v)
 |    .   |
 | v[l-1] |
 |    0   |
 |    .   |
 |    0   |
 -        -

 with y[p] = -s*(v[p]^2 + sum(l<=i<=m) v[i]^2)^1/2. On return, the pth
 element of v is used to store s, and the storage location h is used to
 contain the quantity v[p] - s. The vector v may be stored in a column
 of the input matrix (in which case one uses the householderCol* routines
 and indicates the designated column, vCol), or as a row of the input matrix
 (in which case one uses the householderRow* routines and indicates the row,
 vRow). The elements of v are in either case stored in a pointer to pointer
 to double matrix.

 If mode = HOUSEHOLDER_CALC_MODE_H1 or HOUSEHOLDER_CALC_MODE_H2

 The transformation is applied to the set of vectors stored in the matrix c.
 These may be column vectors beginning with cColStart and ending with
 cColEnd, for which the householder*Col routines are suitable. Alternatively,
 the transformation may be applied to row vectors beginning with cRowStart
 and ending with cRowEnd, in which case the householder*Row routines are
 appropriate. Note that if end < start, c is ignored (in effect the
 identity transformation is performed). The information used to tranform
 c is provided in the elements of v and the scaler h.

 For both modes the user must input the pivot element p, and the indices
 of the elements of the vector v that are to be zeroed (H1 mode) or
 transformed in c (H1 mode, if end >= start, or H2 mode).

 Note that all vectors and matrices have zero-based indices and that
 on input 0 <= p < l and that l <= m. The elemnets of v indexed l through
 m are zero in H1 or transformed in H1/2.
 */

#define HOUSEHOLDER_CALC_MODE_H1     1
#define HOUSEHOLDER_CALC_MODE_H2     2

static void householderColCol(int mode, int p, int l, int m, double **v, int vCol, double *h,
	double **c, int cColStart, int cColEnd) {
	double b, s;
	int i, j;

	if (0 > p || p >= l || l > m) return;
	switch (mode) {

		case HOUSEHOLDER_CALC_MODE_H1:
			for (i=l, s=SQUARE(v[p][vCol]); i<=m; i++) s += SQUARE(v[i][vCol]);
			s = sqrt(s);
			if (v[p][vCol] > 0.0) s *= -1.0;
			*h = v[p][vCol] - s; v[p][vCol] = s;

		case HOUSEHOLDER_CALC_MODE_H2:
			b = v[p][vCol]*(*h);
			if (b == 0.0) return;

			for (j=cColStart; j<=cColEnd; j++) {
				for (i=l, s=c[p][j]*(*h); i<=m; i++) s += c[i][j]*v[i][vCol];
				s /= b;
				c[p][j] += s*(*h);
				for (i=l; i<=m; i++) c[i][j] += s*v[i][vCol];
			}

	}
}

static void householderColRow(int mode, int p, int l, int m, double **v, int vCol, double *h,
    double **c, int cRowStart, int cRowEnd) {
	double b, s;
	int i, j;

	if (0 > p || p >= l || l > m) return;
	switch (mode) {

		case HOUSEHOLDER_CALC_MODE_H1:
			for (i=l, s=SQUARE(v[p][vCol]); i<=m; i++) s += SQUARE(v[i][vCol]);
			s = sqrt(s);
			if (v[p][vCol] > 0.0) s *= -1.0;
			*h = v[p][vCol] - s; v[p][vCol] = s;

		case HOUSEHOLDER_CALC_MODE_H2:
			b = v[p][vCol]*(*h);
			if (b == 0.0) return;

			for (j=cRowStart; j<=cRowEnd; j++) {
				for (i=l, s=c[j][p]*(*h); i<=m; i++) s += c[j][i]*v[i][vCol];
				s /= b;
				c[j][p] += s*(*h);
				for (i=l; i<=m; i++) c[j][i] += s*v[i][vCol];
			}

	}
}

/**
 SUBROUTINE G1 (A,B,CTERM,STERM,SIG)

 COMPUTE ORTHOGONAL ROTATION MATRIX..

 The original version of this code was developed by
 Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratoryc  1973 JUN 12, and published in the book
 "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.

 Revised FEB 1995 to accompany reprinting of the book by SIAM.

 COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))
 (-S,C)         (-S,C)(B)   (   0          )
 COMPUTE SIG = SQRT(A**2+B**2)
 SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT
 SIG MAY BE IN THE SAME LOCATION AS A OR B .
 */
static void G1(double a, double b, double *cterm, double *sterm, double *sig) {
	if (fabs(a) > fabs(b)) {
		double xr = b/a;
		double yr = sqrt(1.0 + xr*xr);
		*cterm = (a > 0.0) ? 1.0/yr : -1.0/yr;
		*sterm = (*cterm)*xr;
		*sig   = fabs(a)*yr;
	} else if (b != 0.0) {
		double xr = a/b;
		double yr = sqrt(1.0 + xr*xr);
		*sterm = (b > 0.0) ? 1.0/yr : -1.0/yr;
		*cterm = (*sterm)*xr;
		*sig   = fabs(b)*yr;
	} else {
		*sig   = 0.0;
		*cterm = 0.0;
		*sterm = 1.0;
	}
}

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
    double *w, double *zz, int *index, int debug) {

  	if (m <= 0 || n <= 0) return FALSE;

	for (int i=0; i<n; i++) { x[i] = 0.0; index[i] = i; }
	int itmax = 3*n;
	int iter  = 0;
	int iz2   = n - 1;
	int iz1   = 0;
	int nsetp = 0;
	int npp1  = 0;

	int mainLoop = TRUE;
	while (mainLoop) {
		if (debug) printf("npp1 %d, nsetp %d, iter %d, iz2 %d, iz1 %d\n", npp1, nsetp, iter, iz2, iz1);

		// quit if all coefficients are already in the solution or if m cols of a have been triangularized.
		if (iz1 > iz2 || nsetp >= m) {
			double sm = 0.0;
		    if (npp1 < m) for (int i=npp1; i<m; i++) sm += b[i]*b[i];
		    else          for (int i=0; i<n; i++) w[i] = 0.0;
		    *rnorm = sqrt(sm);
		    return TRUE;
		}

		// compute components of the dual (negative gradient) vector w().
		for (int iz=iz1; iz<=iz2; iz++) {
			double sm = 0.0;
			for (int l=npp1; l<m; l++) sm += a[l][index[iz]]*b[l];
			w[index[iz]] = sm;
		}

		// find largest positive w(j).
		int wLoop = TRUE;
		int izmax = iz1;
		double up;
		while (wLoop) {
			double wmax = 0.0;
			for (int iz=iz1; iz<=iz2; iz++) {
				if (w[index[iz]] > wmax) {
					wmax  = w[index[iz]];
					izmax = iz;
				}
			}

			// if wmax <= 0.0 then terminate. This condition indicates satisfaction of the kuhn-tucker conditions.
			if (wmax <= 0.0) {
				double sm = 0.0;
		        if (npp1 < m) for (int i=npp1; i<m; i++) sm += b[i]*b[i];
		        else          for (int i=0; i<n; i++) w[i] = 0.0;
		        *rnorm = sqrt(sm);
		        return TRUE;
			}

			// the sign of w(index[izmax]) is ok for index[izmax] to be moved to set p.
			// begin the transformation and check new diagonal element to avoid near linear dependence.
			// original FORTRAN code: call h12 (1,npp1,npp1+1,m,a(1,index[izmax]),1,up,dummy,1,1,0)
			double asave = a[npp1][index[izmax]];
			householderColCol(HOUSEHOLDER_CALC_MODE_H1, npp1, npp1+1, m-1, a,
				index[izmax], &up, NULL, 0, -1);

			double unorm = 0.0;
			if (nsetp != 0) for (int l=0; l<nsetp; l++) unorm += a[l][index[izmax]]*a[l][index[izmax]];
			unorm = sqrt(unorm);

			if (fabs(a[npp1][index[izmax]]*0.01) > DBL_EPSILON*fabs(unorm)) {
				// col index[izmax] is sufficiently independent.  copy b into zz, update zz and solve for ztest ( = proposed new value for x(index[izmax]) ).
				// original FORTRAN code: call h12 (2,npp1,npp1+1,m,a(1,index[izmax]),1,up,zz,1,1,1)
				for (int l=0; l<m; l++) zz[l] = b[l];
				householderColRow(HOUSEHOLDER_CALC_MODE_H2, npp1, npp1+1, m-1, a,
				    index[izmax], &up, &zz, 0, 0);

	            if (zz[npp1]/a[npp1][index[izmax]] > 0.0) break; // out of wLoop
			}

			// reject index[izmax] as a candidate to be moved from set z to set p.
			// restore a(npp1,index[izmax]), set w(index[izmax])=0., and loop back to test dual coeffs again.
			a[npp1][index[izmax]] = asave;
			w[index[izmax]]       = 0.0;
		} // end while on wloop

		// The index index(izmax) has been selected to be moved from set z to set p.
		// Update b, update indices, apply householder transformations to cols in new set z,  zero subdiagonal elts in col index(izmax),  set w(index(izmax))=0.
		for (int l=0; l<m; l++) b[l] = zz[l];
		int indexTemp = index[izmax];
		index[izmax] = index[iz1];
		index[iz1]   = indexTemp;
		nsetp        = npp1+1;
	    iz1++;
		npp1++;

		// call h12 (2,nsetp,npp1,m,a(1,indexTemp),1,up,a(1,index[jz]),1,mda,1)
		if (iz1 <= iz2) for (int jz=iz1; jz<=iz2; jz++)
			householderColCol(HOUSEHOLDER_CALC_MODE_H2, nsetp-1, npp1, m-1, a,
				indexTemp, &up, a, index[jz], index[jz]);

		if (nsetp != m) for (int l=npp1; l<m; l++) a[l][indexTemp] = 0.0;
		w[indexTemp] = 0.0;

		//solve the triangular system. store the solution temporarily in zz().
		for (int l=0; l<nsetp; l++) {
			if (l != 0) for (int ii=0; ii<=(nsetp-(l+1)); ii++)
			    zz[ii] -= a[ii][index[nsetp-(l+1)+1]]*zz[nsetp-(l+1)+1];
			zz[nsetp-(l+1)] /= a[nsetp-(l+1)][index[nsetp-(l+1)]];
		}

		int secondaryLoop = TRUE;
		while (secondaryLoop) {
			iter++;
			if (iter > itmax) {
				double sm = 0.0;
		        if (npp1 < m) for (int i=npp1; i<m; i++) sm += b[i]*b[i];
		        else          for (int i=0; i<n; i++) w[i] = 0.0;
		        *rnorm = sqrt(sm);
		        return FALSE;
			}

			// see if all new constrained coeffs are feasible. if not compute alpha.
			double alpha = 2.0;
			int alphaIndex = 0;
			for (int ip=0; ip<nsetp; ip++) {
				if (zz[ip] <= 0.0) {
					double t = -x[index[ip]]/(zz[ip] - x[index[ip]]);
					if (alpha > t) {
						alpha      = t;
						alphaIndex = ip-1;
					}
				}
			}

			// if all new constrained coeffs are feasible then alpha will still = 2. if so exit from secondary loop to main loop.
			if (alpha == 2.0) break; // out of the while loop

			// otherwise use alpha which will be between 0. and 1. to interpolate between the old x and the new zz.
			for (int ip=0; ip<nsetp; ip++) x[index[ip]] += alpha*(zz[ip]-x[index[ip]]);

			// modify a and b and the index arrays to move coefficient i from set p to set z.
			int loop = TRUE;
			int superAlphaIndex = index[alphaIndex+1];
			while (loop) {
				x[superAlphaIndex] = 0.0;

				if (alphaIndex != (nsetp-1)) {
					alphaIndex++;
					for (int j=(alphaIndex+1); j<nsetp; j++) {
						index[j-1] = index[j];
						// call g1 (a(j-1,index[j]),a(j,index[j]),cc,ss,a(j-1,index[j]))
						double cc, ss;
						G1(a[j-1][index[j]], a[j][index[j]], &cc, &ss, &a[j-1][index[j]]);
						a[j][index[j]] = 0.0;
						for (int l=0; l<n; l++) {
							if (l != index[j]) {
								// apply procedure g2 (cc,ss,a(j-1,l),a(j,l))
								double temp = a[j-1][l];
								a[j-1][l] =  cc*temp + ss*a[j][l];
								a[j][l]   = -ss*temp + cc*a[j][l];
							}
						}

						// apply procedure g2 (cc,ss,b(j-1),b(j))
						double temp = b[j-1];
						b[j-1] =  cc*temp + ss*b[j];
						b[j]   = -ss*temp + cc*b[j];
					}
				}

				npp1 = nsetp-1;
				nsetp--;
				iz1--;
				index[iz1] = superAlphaIndex;

				// see if the remaining coeffs in set p are feasible.  they should be because of the way alpha was determined.
				// if any are infeasible it is due to round-off error.  any that are nonpositive will be set to zero and moved from set p to set z.
				loop = FALSE;
				for (int jj=0; jj<nsetp; jj++) {
					superAlphaIndex = index[jj];
					if (x[superAlphaIndex] <= 0.0) { loop = TRUE; break; }
				}
			}

			// copy b( ) into zz( ).  then solve again and loop back.
			for (int i=0; i<m; i++) zz[i] = b[i];
			for (int l=0; l<nsetp; l++) {
				if (l != 0) for (int ii=0; ii<=(nsetp-(l+1)); ii++)
				    zz[ii] -= a[ii][index[nsetp-(l+1)+1]]*zz[nsetp-(l+1)+1];
				zz[nsetp-(l+1)] /= a[nsetp-(l+1)][index[nsetp-(l+1)]];
			}

		} // end secondary while loop

		for (int ip=0; ip<nsetp; ip++) x[index[ip]] = zz[ip];
		//all new coeffs are positive.  loop back to beginning.
	} // end primary while loop

	double sm = 0.0;
	if (npp1 < m) for (int i=npp1; i<m; i++) sm += b[i]*b[i];
	else          for (int i=0; i<n; i++) w[i] = 0.0;
	*rnorm = sqrt(sm);
	return FALSE;
}
