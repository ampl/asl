/* dl7upd.f -- translated by f2c (version 20160102).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Subroutine */ int dl7upd_(doublereal *beta, doublereal *gamma, doublereal *
	l, doublereal *lambda, doublereal *lplus, integer *n, doublereal *w, 
	doublereal *z__)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal a, b;
    static integer i__, j, k;
    static doublereal s, bj, gj;
    static integer ij, jj;
    static doublereal lj, wj, nu, zj;
    static integer jp1, nm1, np1;
    static doublereal eta, lij, ljj, theta;


/*  ***  COMPUTE LPLUS = SECANT UPDATE OF L  *** */

/*  ***  PARAMETER DECLARATIONS  *** */

/*     DIMENSION L(N*(N+1)/2), LPLUS(N*(N+1)/2) */

/* --------------------------  PARAMETER USAGE  -------------------------- */

/*   BETA = SCRATCH VECTOR. */
/*  GAMMA = SCRATCH VECTOR. */
/*      L (INPUT) LOWER TRIANGULAR MATRIX, STORED ROWWISE. */
/* LAMBDA = SCRATCH VECTOR. */
/*  LPLUS (OUTPUT) LOWER TRIANGULAR MATRIX, STORED ROWWISE, WHICH MAY */
/*             OCCUPY THE SAME STORAGE AS  L. */
/*      N (INPUT) LENGTH OF VECTOR PARAMETERS AND ORDER OF MATRICES. */
/*      W (INPUT, DESTROYED ON OUTPUT) RIGHT SINGULAR VECTOR OF RANK 1 */
/*             CORRECTION TO  L. */
/*      Z (INPUT, DESTROYED ON OUTPUT) LEFT SINGULAR VECTOR OF RANK 1 */
/*             CORRECTION TO  L. */

/* -------------------------------  NOTES  ------------------------------- */

/*  ***  APPLICATION AND USAGE RESTRICTIONS  *** */

/*        THIS ROUTINE UPDATES THE CHOLESKY FACTOR  L  OF A SYMMETRIC */
/*     POSITIVE DEFINITE MATRIX TO WHICH A SECANT UPDATE IS BEING */
/*     APPLIED -- IT COMPUTES A CHOLESKY FACTOR  LPLUS  OF */
/*     L * (I + Z*W**T) * (I + W*Z**T) * L**T.  IT IS ASSUMED THAT  W */
/*     AND  Z  HAVE BEEN CHOSEN SO THAT THE UPDATED MATRIX IS STRICTLY */
/*     POSITIVE DEFINITE. */

/*  ***  ALGORITHM NOTES  *** */

/*        THIS CODE USES RECURRENCE 3 OF REF. 1 (WITH D(J) = 1 FOR ALL J) */
/*     TO COMPUTE  LPLUS  OF THE FORM  L * (I + Z*W**T) * Q,  WHERE  Q */
/*     IS AN ORTHOGONAL MATRIX THAT MAKES THE RESULT LOWER TRIANGULAR. */
/*        LPLUS MAY HAVE SOME NEGATIVE DIAGONAL ELEMENTS. */

/*  ***  REFERENCES  *** */

/* 1.  GOLDFARB, D. (1976), FACTORIZED VARIABLE METRIC METHODS FOR UNCON- */
/*             STRAINED OPTIMIZATION, MATH. COMPUT. 30, PP. 796-811. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY (FALL 1979). */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH SUPPORTED */
/*     BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS MCS-7600324 AND */
/*     MCS-7906671. */

/* ------------------------  EXTERNAL QUANTITIES  ------------------------ */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/* --------------------------  LOCAL VARIABLES  -------------------------- */


/*  ***  DATA INITIALIZATIONS  *** */

/* /6 */
/*     DATA ONE/1.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --l;
    --lplus;
    --z__;
    --w;
    --lambda;
    --gamma;
    --beta;

    /* Function Body */
    nu = 1.;
    eta = 0.;
    if (*n <= 1) {
	goto L30;
    }
    nm1 = *n - 1;

/*  ***  TEMPORARILY STORE S(J) = SUM OVER K = J+1 TO N OF W(K)**2 IN */
/*  ***  LAMBDA(J). */

    s = 0.;
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n - i__;
/* Computing 2nd power */
	d__1 = w[j + 1];
	s += d__1 * d__1;
	lambda[j] = s;
/* L10: */
    }

/*  ***  COMPUTE LAMBDA, GAMMA, AND BETA BY GOLDFARB*S RECURRENCE 3. */

    i__1 = nm1;
    for (j = 1; j <= i__1; ++j) {
	wj = w[j];
	a = nu * z__[j] - eta * wj;
	theta = a * wj + 1.;
	s = a * lambda[j];
/* Computing 2nd power */
	d__1 = theta;
	lj = sqrt(d__1 * d__1 + a * s);
	if (theta > 0.) {
	    lj = -lj;
	}
	lambda[j] = lj;
	b = theta * wj + s;
	gamma[j] = b * nu / lj;
	beta[j] = (a - b * eta) / lj;
	nu = -nu / lj;
/* Computing 2nd power */
	d__1 = a;
	eta = -(eta + d__1 * d__1 / (theta - lj)) / lj;
/* L20: */
    }
L30:
    lambda[*n] = (nu * z__[*n] - eta * w[*n]) * w[*n] + 1.;

/*  ***  UPDATE L, GRADUALLY OVERWRITING  W  AND  Z  WITH  L*W  AND  L*Z. */

    np1 = *n + 1;
    jj = *n * (*n + 1) / 2;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	j = np1 - k;
	lj = lambda[j];
	ljj = l[jj];
	lplus[jj] = lj * ljj;
	wj = w[j];
	w[j] = ljj * wj;
	zj = z__[j];
	z__[j] = ljj * zj;
	if (k == 1) {
	    goto L50;
	}
	bj = beta[j];
	gj = gamma[j];
	ij = jj + j;
	jp1 = j + 1;
	i__2 = *n;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	    lij = l[ij];
	    lplus[ij] = lj * lij + bj * w[i__] + gj * z__[i__];
	    w[i__] += lij * wj;
	    z__[i__] += lij * zj;
	    ij += i__;
/* L40: */
	}
L50:
	jj -= j;
/* L60: */
    }

/* L999: */
    return 0;
/*  ***  LAST CARD OF DL7UPD FOLLOWS  *** */
} /* dl7upd_ */

