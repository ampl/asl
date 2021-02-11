/* dn2gb.f -- translated by f2c (version 20160102).
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

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int dn2gb_(integer *n, integer *p, doublereal *x, doublereal 
	*b, S_fp calcr, S_fp calcj, integer *iv, integer *liv, integer *lv, 
	doublereal *v, integer *uiparm, doublereal *urparm, U_fp ufparm)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer d1, n1, n2, r1, nf, dr1, rd1, iv1;
    extern /* Subroutine */ int drn2gb_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), divset_(integer *, integer *, 
	    integer *, integer *, doublereal *);


/*  ***  VERSION OF NL2SOL THAT HANDLES SIMPLE BOUNDS ON X  *** */

/*  ***  PARAMETERS  *** */

/* /6 */
/*     INTEGER IV(LIV), UIPARM(1) */
/*     DOUBLE PRECISION X(P), B(2,P), V(LV), URPARM(1) */
/* /7 */
/* / */

/*  ***  DISCUSSION  *** */

/*        NOTE... NL2SOL (MENTIONED BELOW) IS A CODE FOR SOLVING */
/*     NONLINEAR LEAST-SQUARES PROBLEMS.  IT IS DESCRIBED IN */
/*     ACM TRANS. MATH. SOFTWARE, VOL. 7 (1981), PP. 369-383 */
/*     (AN ADAPTIVE NONLINEAR LEAST-SQUARES ALGORITHM, BY J.E. DENNIS, */
/*     D.M. GAY, AND R.E. WELSCH). */

/*        LIV GIVES THE LENGTH OF IV.  IT MUST BE AT LEAST 82 + 4*P. */
/*     IF NOT, THEN  DN2GB RETURNS WITH IV(1) = 15.  WHEN  DN2GB */
/*     RETURNS, THE MINIMUM ACCEPTABLE VALUE OF LIV IS STORED IN */
/*     IV(LASTIV) = IV(44), (PROVIDED THAT LIV .GE. 44). */

/*        LV GIVES THE LENGTH OF V.  THE MINIMUM VALUE FOR LV IS */
/*     LV0 = 105 + P*(N + 2*P + 21) + 2*N.  IF LV IS SMALLER THAN THIS, */
/*     THEN  DN2GB RETURNS WITH IV(1) = 16.  WHEN  DN2GB RETURNS, THE */
/*     MINIMUM ACCEPTABLE VALUE OF LV IS STORED IN IV(LASTV) = IV(45) */
/*     (PROVIDED LIV .GE. 45). */

/*        RETURN CODES AND CONVERGENCE TOLERANCES ARE THE SAME AS FOR */
/*     NL2SOL, WITH SOME SMALL EXTENSIONS... IV(1) = 15 MEANS LIV WAS */
/*     TOO SMALL.   IV(1) = 16 MEANS LV WAS TOO SMALL. */

/*        THERE ARE TWO NEW V INPUT COMPONENTS...  V(LMAXS) = V(36) AND */
/*     V(SCTOL) = V(37) SERVE WHERE V(LMAX0) AND V(RFCTOL) FORMERLY DID */
/*     IN THE SINGULAR CONVERGENCE TEST -- SEE THE NL2SOL DOCUMENTATION. */

/*  ***  BOUNDS  *** */

/*     THE BOUNDS  B(1,I) .LE. X(I) .LE. B(2,I), I = 1(1)P, ARE ENFORCED. */

/*  ***  DEFAULT VALUES  *** */

/*        DEFAULT VALUES ARE PROVIDED BY SUBROUTINE DIVSET, RATHER THAN */
/*     DFAULT.  THE CALLING SEQUENCE IS... */
/*             CALL DIVSET(1, IV, LIV, LV, V) */
/*     THE FIRST PARAMETER IS AN INTEGER 1.  IF LIV AND LV ARE LARGE */
/*     ENOUGH FOR DIVSET, THEN DIVSET SETS IV(1) TO 12.  OTHERWISE IT */
/*     SETS IV(1) TO 15 OR 16.  CALLING  DN2GB WITH IV(1) = 0 CAUSES ALL */
/*     DEFAULT VALUES TO BE USED FOR THE INPUT COMPONENTS OF IV AND V. */
/*        IF YOU FIRST CALL DIVSET, THEN SET IV(1) TO 13 AND CALL  DN2GB, */
/*     THEN STORAGE ALLOCATION ONLY WILL BE PERFORMED.  IN PARTICULAR, */
/*     IV(D) = IV(27), IV(J) = IV(70), AND IV(R) = IV(61) WILL BE SET */
/*     TO THE FIRST SUBSCRIPT IN V OF THE SCALE VECTOR, THE JACOBIAN */
/*     MATRIX, AND THE RESIDUAL VECTOR RESPECTIVELY, PROVIDED LIV AND LV */
/*     ARE LARGE ENOUGH.  IF SO, THEN  DN2GB RETURNS WITH IV(1) = 14. */
/*     WHEN CALLED WITH IV(1) = 14,  DN2GB ASSUMES THAT STORAGE HAS */
/*     BEEN ALLOCATED, AND IT BEGINS THE MINIMIZATION ALGORITHM. */

/*  ***  SCALE VECTOR  *** */

/*        ONE DIFFERENCE WITH NL2SOL IS THAT THE SCALE VECTOR D IS */
/*     STORED IN V, STARTING AT A DIFFERENT SUBSCRIPT.  THE STARTING */
/*     SUBSCRIPT VALUE IS STILL STORED IN IV(D) = IV(27).  THE */
/*     DISCUSSION OF DEFAULT VALUES ABOVE TELLS HOW TO HAVE IV(D) SET */
/*     BEFORE THE ALGORITHM IS STARTED. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY. */

/*  ***  EXTERNAL SUBROUTINES  *** */

/* DIVSET.... PROVIDES DEFAULT IV AND V INPUT COMPONENTS. */
/* DRN2GB... CARRIES OUT OPTIMIZATION ITERATIONS. */

/*  ***  LOCAL VARIABLES  *** */


/*  ***  IV COMPONENTS  *** */

/* /6 */
/*     DATA D/27/, J/70/, NEXTV/47/, NFCALL/6/, NFGCAL/7/, R/61/, */
/*    1     REGD0/82/, TOOBIG/2/, VNEED/4/ */
/* /7 */
/* / */
/* ---------------------------------  BODY  ------------------------------ */

    /* Parameter adjustments */
    b -= 3;
    --x;
    --iv;
    --v;
    --uiparm;
    --urparm;

    /* Function Body */
    if (iv[1] == 0) {
	divset_(&c__1, &iv[1], liv, lv, &v[1]);
    }
    iv1 = iv[1];
    if (iv1 == 14) {
	goto L10;
    }
    if (iv1 > 2 && iv1 < 12) {
	goto L10;
    }
    if (iv1 == 12) {
	iv[1] = 13;
    }
    if (iv[1] == 13) {
	iv[4] = iv[4] + *p + *n * (*p + 2);
    }
    drn2gb_(&b[3], &x[1], &v[1], &iv[1], liv, lv, n, n, &n1, &n2, p, &v[1], &
	    v[1], &v[1], &x[1]);
    if (iv[1] != 14) {
	goto L999;
    }

/*  ***  STORAGE ALLOCATION  *** */

    iv[27] = iv[47];
    iv[61] = iv[27] + *p;
    iv[82] = iv[61] + *n;
    iv[70] = iv[82] + *n;
    iv[47] = iv[70] + *n * *p;
    if (iv1 == 13) {
	goto L999;
    }

L10:
    d1 = iv[27];
    dr1 = iv[70];
    r1 = iv[61];
    rd1 = iv[82];

L20:
    drn2gb_(&b[3], &v[d1], &v[dr1], &iv[1], liv, lv, n, n, &n1, &n2, p, &v[r1]
	    , &v[rd1], &v[1], &x[1]);
    if ((i__1 = iv[1] - 2) < 0) {
	goto L30;
    } else if (i__1 == 0) {
	goto L50;
    } else {
	goto L999;
    }

/*  ***  NEW FUNCTION VALUE (R VALUE) NEEDED  *** */

L30:
    nf = iv[6];
    (*calcr)(n, p, &x[1], &nf, &v[r1], &uiparm[1], &urparm[1], (U_fp)ufparm);
    if (nf > 0) {
	goto L40;
    }
    iv[2] = 1;
    goto L20;
L40:
    if (iv[1] > 0) {
	goto L20;
    }

/*  ***  COMPUTE DR = GRADIENT OF R COMPONENTS  *** */

L50:
    (*calcj)(n, p, &x[1], &iv[7], &v[dr1], &uiparm[1], &urparm[1], (U_fp)
	    ufparm);
    if (iv[7] == 0) {
	iv[2] = 1;
    }
    goto L20;

L999:
    return 0;

/*  ***  LAST CARD OF  DN2GB FOLLOWS  *** */
} /* dn2gb_ */

