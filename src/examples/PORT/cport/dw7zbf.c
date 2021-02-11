/* dw7zbf.f -- translated by f2c (version 20160102).
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

/* Subroutine */ int dw7zbf_(doublereal *l, integer *n, doublereal *s, 
	doublereal *w, doublereal *y, doublereal *z__)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal cs, cy, ys, shs, theta, epsrt;
    extern /* Subroutine */ int dl7ivm_(integer *, doublereal *, doublereal *,
	     doublereal *);
    extern doublereal dd7tpr_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int dl7tvm_(integer *, doublereal *, doublereal *,
	     doublereal *);


/*  ***  COMPUTE  Y  AND  Z  FOR  DL7UPD  CORRESPONDING TO BFGS UPDATE. */

/*     DIMENSION L(N*(N+1)/2) */

/* --------------------------  PARAMETER USAGE  -------------------------- */

/* L (I/O) CHOLESKY FACTOR OF HESSIAN, A LOWER TRIANG. MATRIX STORED */
/*             COMPACTLY BY ROWS. */
/* N (INPUT) ORDER OF  L  AND LENGTH OF  S,  W,  Y,  Z. */
/* S (INPUT) THE STEP JUST TAKEN. */
/* W (OUTPUT) RIGHT SINGULAR VECTOR OF RANK 1 CORRECTION TO L. */
/* Y (INPUT) CHANGE IN GRADIENTS CORRESPONDING TO S. */
/* Z (OUTPUT) LEFT SINGULAR VECTOR OF RANK 1 CORRECTION TO L. */

/* -------------------------------  NOTES  ------------------------------- */

/*  ***  ALGORITHM NOTES  *** */

/*        WHEN  S  IS COMPUTED IN CERTAIN WAYS, E.G. BY  GQTSTP  OR */
/*     DBLDOG,  IT IS POSSIBLE TO SAVE N**2/2 OPERATIONS SINCE  (L**T)*S */
/*     OR  L*(L**T)*S IS THEN KNOWN. */
/*        IF THE BFGS UPDATE TO L*(L**T) WOULD REDUCE ITS DETERMINANT TO */
/*     LESS THAN EPS TIMES ITS OLD VALUE, THEN THIS ROUTINE IN EFFECT */
/*     REPLACES  Y  BY  THETA*Y + (1 - THETA)*L*(L**T)*S,  WHERE  THETA */
/*     (BETWEEN 0 AND 1) IS CHOSEN TO MAKE THE REDUCTION FACTOR = EPS. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY (FALL 1979). */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH SUPPORTED */
/*     BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS MCS-7600324 AND */
/*     MCS-7906671. */

/* ------------------------  EXTERNAL QUANTITIES  ------------------------ */

/*  ***  FUNCTIONS AND SUBROUTINES CALLED  *** */

/* DD7TPR RETURNS INNER PRODUCT OF TWO VECTORS. */
/* DL7IVM MULTIPLIES L**-1 TIMES A VECTOR. */
/* DL7TVM MULTIPLIES L**T TIMES A VECTOR. */

/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/* --------------------------  LOCAL VARIABLES  -------------------------- */


/*  ***  DATA INITIALIZATIONS  *** */

/* /6 */
/*     DATA EPS/0.1D+0/, ONE/1.D+0/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --l;
    --z__;
    --y;
    --w;
    --s;

    /* Function Body */
    dl7tvm_(n, &w[1], &l[1], &s[1]);
    shs = dd7tpr_(n, &w[1], &w[1]);
    ys = dd7tpr_(n, &y[1], &s[1]);
    if (ys >= shs * .1) {
	goto L10;
    }
    theta = shs * .90000000000000002 / (shs - ys);
    epsrt = sqrt(.1);
    cy = theta / (shs * epsrt);
    cs = ((theta - 1.) / epsrt + 1.) / shs;
    goto L20;
L10:
    cy = 1. / (sqrt(ys) * sqrt(shs));
    cs = 1. / shs;
L20:
    dl7ivm_(n, &z__[1], &l[1], &y[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	z__[i__] = cy * z__[i__] - cs * w[i__];
    }

/* L999: */
    return 0;
/*  ***  LAST CARD OF DW7ZBF FOLLOWS  *** */
} /* dw7zbf_ */

