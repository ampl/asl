/* dmngb.f -- translated by f2c (version 20160102).
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

static integer c__2 = 2;

/* Subroutine */ int dmngb_(integer *n, doublereal *d__, doublereal *x, 
	doublereal *b, S_fp calcf, S_fp calcg, integer *iv, integer *liv, 
	integer *lv, doublereal *v, integer *uiparm, doublereal *urparm, U_fp 
	ufparm)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal f;
    static integer g1, nf, iv1;
    extern /* Subroutine */ int drmngb_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *), divset_(integer *, 
	    integer *, integer *, integer *, doublereal *);


/*  ***  MINIMIZE GENERAL SIMPLY BOUNDED OBJECTIVE FUNCTION USING  *** */
/*  ***  ANALYTIC GRADIENT AND HESSIAN APPROX. FROM SECANT UPDATE  *** */

/* /6S */
/*     INTEGER IV(LIV), UIPARM(1) */
/*     DOUBLE PRECISION D(N), X(N), B(2,N), V(LV), URPARM(1) */
/* /7S */
/* / */
/*     DIMENSION IV(59 + N), V(71 + N*(N+21)/2), UIPARM(*), URPARM(*) */

/*  ***  DISCUSSION  *** */

/*        THIS ROUTINE IS LIKE  DMNG, EXCEPT FOR THE EXTRA PARAMETER B, */
/*     AN ARRAY OF LOWER AND UPPER BOUNDS ON X...  DMNGB ENFORCES THE */
/*     CONSTRAINTS THAT  B(1,I) .LE. X(I) .LE. B(2,I), I = 1(1)N. */
/*     (INSTEAD OF CALLING DRMNG,  DMNGB CALLS DRMNGB.) */
/* . */

/* ----------------------------  DECLARATIONS  --------------------------- */


/* DIVSET.... SUPPLIES DEFAULT IV AND V INPUT COMPONENTS. */
/* DRMNGB... REVERSE-COMMUNICATION ROUTINE THAT CARRIES OUT  DMNG ALGO- */
/*             RITHM. */


/*  ***  SUBSCRIPTS FOR IV   *** */


/* /6 */
/*     DATA NEXTV/47/, NFCALL/6/, NFGCAL/7/, G/28/, TOOBIG/2/, VNEED/4/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    b -= 3;
    --x;
    --d__;
    --iv;
    --v;
    --uiparm;
    --urparm;

    /* Function Body */
    if (iv[1] == 0) {
	divset_(&c__2, &iv[1], liv, lv, &v[1]);
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
	iv[4] += *n;
    }
    drmngb_(&b[3], &d__[1], &f, &v[1], &iv[1], liv, lv, n, &v[1], &x[1]);
    if (iv[1] != 14) {
	goto L999;
    }

/*  ***  STORAGE ALLOCATION */

    iv[28] = iv[47];
    iv[47] = iv[28] + *n;
    if (iv1 == 13) {
	goto L999;
    }

L10:
    g1 = iv[28];

L20:
    drmngb_(&b[3], &d__[1], &f, &v[g1], &iv[1], liv, lv, n, &v[1], &x[1]);
    if ((i__1 = iv[1] - 2) < 0) {
	goto L30;
    } else if (i__1 == 0) {
	goto L40;
    } else {
	goto L999;
    }

L30:
    nf = iv[6];
    (*calcf)(n, &x[1], &nf, &f, &uiparm[1], &urparm[1], (U_fp)ufparm);
    if (nf <= 0) {
	iv[2] = 1;
    }
    goto L20;

L40:
    nf = iv[7];
    (*calcg)(n, &x[1], &nf, &v[g1], &uiparm[1], &urparm[1], (U_fp)ufparm);
    if (nf <= 0) {
	iv[2] = 1;
    }
    goto L20;

L999:
    return 0;
/*  ***  LAST CARD OF  DMNGB FOLLOWS  *** */
} /* dmngb_ */

