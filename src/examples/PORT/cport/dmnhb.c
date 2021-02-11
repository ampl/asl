/* dmnhb.f -- translated by f2c (version 20160102).
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

/* Subroutine */ int dmnhb_(integer *n, doublereal *d__, doublereal *x, 
	doublereal *b, S_fp calcf, S_fp calcgh, integer *iv, integer *liv, 
	integer *lv, doublereal *v, integer *uiparm, doublereal *urparm, U_fp 
	ufparm)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal f;
    static integer g1, h1, lh, nf, iv1;
    extern /* Subroutine */ int drmnhb_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *), 
	    divset_(integer *, integer *, integer *, integer *, doublereal *);


/*  ***  MINIMIZE GENERAL SIMPLY BOUNDED OBJECTIVE FUNCTION USING   *** */
/*  ***  (ANALYTIC) GRADIENT AND HESSIAN PROVIDED BY THE CALLER.    *** */

/* /6S */
/*     INTEGER IV(LIV), UIPARM(1) */
/*     DOUBLE PRECISION B(2,N), D(N), X(N), V(LV), URPARM(1) */
/* /7S */
/* / */
/*     DIMENSION IV(59 + 3*N), V(78 + N*(N+15)) */

/* ------------------------------  DISCUSSION  --------------------------- */

/*        THIS ROUTINE IS LIKE  DMNGB, EXCEPT THAT THE SUBROUTINE PARA- */
/*     METER CALCG OF  DMNGB (WHICH COMPUTES THE GRADIENT OF THE OBJEC- */
/*     TIVE FUNCTION) IS REPLACED BY THE SUBROUTINE PARAMETER CALCGH, */
/*     WHICH COMPUTES BOTH THE GRADIENT AND (LOWER TRIANGLE OF THE) */
/*     HESSIAN OF THE OBJECTIVE FUNCTION.  THE CALLING SEQUENCE IS... */
/*             CALL CALCGH(N, X, NF, G, H, UIPARM, URPARM, UFPARM) */
/*     PARAMETERS N, X, NF, G, UIPARM, URPARM, AND UFPARM ARE THE SAME */
/*     AS FOR  DMNGB, WHILE H IS AN ARRAY OF LENGTH N*(N+1)/2 IN WHICH */
/*     CALCGH MUST STORE THE LOWER TRIANGLE OF THE HESSIAN AT X.  START- */
/*     ING AT H(1), CALCGH MUST STORE THE HESSIAN ENTRIES IN THE ORDER */
/*     (1,1), (2,1), (2,2), (3,1), (3,2), (3,3), ... */
/*        THE VALUE PRINTED (BY DITSUM) IN THE COLUMN LABELLED STPPAR */
/*     IS THE LEVENBERG-MARQUARDT USED IN COMPUTING THE CURRENT STEP. */
/*     ZERO MEANS A FULL NEWTON STEP.  IF THE SPECIAL CASE DESCRIBED IN */
/*     REF. 1 IS DETECTED, THEN STPPAR IS NEGATED.  THE VALUE PRINTED */
/*     IN THE COLUMN LABELLED NPRELDF IS ZERO IF THE CURRENT HESSIAN */
/*     IS NOT POSITIVE DEFINITE. */
/*        IT SOMETIMES PROVES WORTHWHILE TO LET D BE DETERMINED FROM THE */
/*     DIAGONAL OF THE HESSIAN MATRIX BY SETTING IV(DTYPE) = 1 AND */
/*     V(DINIT) = 0.  THE FOLLOWING IV AND V COMPONENTS ARE RELEVANT... */

/* IV(DTOL)..... IV(59) GIVES THE STARTING SUBSCRIPT IN V OF THE DTOL */
/*             ARRAY USED WHEN D IS UPDATED.  (IV(DTOL) CAN BE */
/*             INITIALIZED BY CALLING  DMNHB WITH IV(1) = 13.) */
/* IV(DTYPE).... IV(16) TELLS HOW THE SCALE VECTOR D SHOULD BE CHOSEN. */
/*             IV(DTYPE) .LE. 0 MEANS THAT D SHOULD NOT BE UPDATED, AND */
/*             IV(DTYPE) .GE. 1 MEANS THAT D SHOULD BE UPDATED AS */
/*             DESCRIBED BELOW WITH V(DFAC).  DEFAULT = 0. */
/* V(DFAC)..... V(41) AND THE DTOL AND D0 ARRAYS (SEE V(DTINIT) AND */
/*             V(D0INIT)) ARE USED IN UPDATING THE SCALE VECTOR D WHEN */
/*             IV(DTYPE) .GT. 0.  (D IS INITIALIZED ACCORDING TO */
/*             V(DINIT), DESCRIBED IN  DMNG.)  LET */
/*                  D1(I) = MAX(SQRT(ABS(H(I,I))), V(DFAC)*D(I)), */
/*             WHERE H(I,I) IS THE I-TH DIAGONAL ELEMENT OF THE CURRENT */
/*             HESSIAN.  IF IV(DTYPE) = 1, THEN D(I) IS SET TO D1(I) */
/*             UNLESS D1(I) .LT. DTOL(I), IN WHICH CASE D(I) IS SET TO */
/*                  MAX(D0(I), DTOL(I)). */
/*             IF IV(DTYPE) .GE. 2, THEN D IS UPDATED DURING THE FIRST */
/*             ITERATION AS FOR IV(DTYPE) = 1 (AFTER ANY INITIALIZATION */
/*             DUE TO V(DINIT)) AND IS LEFT UNCHANGED THEREAFTER. */
/*             DEFAULT = 0.6. */
/* V(DTINIT)... V(39), IF POSITIVE, IS THE VALUE TO WHICH ALL COMPONENTS */
/*             OF THE DTOL ARRAY (SEE V(DFAC)) ARE INITIALIZED.  IF */
/*             V(DTINIT) = 0, THEN IT IS ASSUMED THAT THE CALLER HAS */
/*             STORED DTOL IN V STARTING AT V(IV(DTOL)). */
/*             DEFAULT = 10**-6. */
/* V(D0INIT)... V(40), IF POSITIVE, IS THE VALUE TO WHICH ALL COMPONENTS */
/*             OF THE D0 VECTOR (SEE V(DFAC)) ARE INITIALIZED.  IF */
/*             V(DFAC) = 0, THEN IT IS ASSUMED THAT THE CALLER HAS */
/*             STORED D0 IN V STARTING AT V(IV(DTOL)+N).  DEFAULT = 1.0. */

/*  ***  REFERENCE  *** */

/* 1. GAY, D.M. (1981), COMPUTING OPTIMAL LOCALLY CONSTRAINED STEPS, */
/*         SIAM J. SCI. STATIST. COMPUT. 2, PP. 186-197. */
/* . */
/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY (WINTER, SPRING 1983). */

/* ----------------------------  DECLARATIONS  --------------------------- */


/* DIVSET.... PROVIDES DEFAULT INPUT VALUES FOR IV AND V. */
/* DRMNHB... REVERSE-COMMUNICATION ROUTINE THAT DOES  DMNHB ALGORITHM. */


/*  ***  SUBSCRIPTS FOR IV   *** */


/* /6 */
/*     DATA NEXTV/47/, NFCALL/6/, NFGCAL/7/, G/28/, H/56/, TOOBIG/2/, */
/*    1     VNEED/4/ */
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
    lh = *n * (*n + 1) / 2;
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
	iv[4] += *n * (*n + 3) / 2;
    }
    drmnhb_(&b[3], &d__[1], &f, &v[1], &v[1], &iv[1], &lh, liv, lv, n, &v[1], 
	    &x[1]);
    if (iv[1] != 14) {
	goto L999;
    }

/*  ***  STORAGE ALLOCATION */

    iv[28] = iv[47];
    iv[56] = iv[28] + *n;
    iv[47] = iv[56] + *n * (*n + 1) / 2;
    if (iv1 == 13) {
	goto L999;
    }

L10:
    g1 = iv[28];
    h1 = iv[56];

L20:
    drmnhb_(&b[3], &d__[1], &f, &v[g1], &v[h1], &iv[1], &lh, liv, lv, n, &v[1]
	    , &x[1]);
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
    (*calcgh)(n, &x[1], &nf, &v[g1], &v[h1], &uiparm[1], &urparm[1], (U_fp)
	    ufparm);
    if (nf <= 0) {
	iv[2] = 1;
    }
    goto L20;

L999:
    return 0;
/*  ***  LAST CARD OF  DMNHB FOLLOWS  *** */
} /* dmnhb_ */

