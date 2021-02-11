/* drmnhb.f -- translated by f2c (version 20160102).
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
static integer c_n1 = -1;
static doublereal c_b38 = 1.;
static doublereal c_b43 = -1.;

/* Subroutine */ int drmnhb_(doublereal *b, doublereal *d__, doublereal *fx, 
	doublereal *g, doublereal *h__, integer *iv, integer *lh, integer *
	liv, integer *lv, integer *n, doublereal *v, doublereal *x)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal t;
    static integer w1;
    static doublereal gi;
    static integer x01, x11;
    static doublereal xi;
    static integer dg1, td1, tg1, ipi, ipn, nn1o2, temp0, temp1, ipiv2, step0,
	     step1, dummy;
    extern logical stopx_(integer *);
    extern /* Subroutine */ int dd7dup_(doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, doublereal *), dg7qsb_(
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal dd7tpr_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int da7sst_(integer *, integer *, integer *, 
	    doublereal *), dv7scp_(integer *, doublereal *, doublereal *);
    extern doublereal dv2nrm_(integer *, doublereal *);
    extern /* Subroutine */ int ds7ipr_(integer *, integer *, doublereal *), 
	    dv7ipr_(integer *, integer *, doublereal *), ds7lvm_(integer *, 
	    doublereal *, doublereal *, doublereal *), dv2axy_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dv7cpy_(
	    integer *, doublereal *, doublereal *), dv7vmp_(integer *, 
	    doublereal *, doublereal *, doublereal *, integer *), i7pnvr_(
	    integer *, integer *, integer *), dparck_(integer *, doublereal *,
	     integer *, integer *, integer *, integer *, doublereal *);
    extern doublereal drldst_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    extern /* Subroutine */ int divset_(integer *, integer *, integer *, 
	    integer *, doublereal *), ditsum_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *);
    static integer lstgst, rstrst;


/*  ***  CARRY OUT  DMNHB (SIMPLY BOUNDED MINIMIZATION) ITERATIONS, */
/*  ***  USING HESSIAN MATRIX PROVIDED BY THE CALLER. */

/*  ***  PARAMETER DECLARATIONS  *** */


/* --------------------------  PARAMETER USAGE  -------------------------- */

/* D.... SCALE VECTOR. */
/* FX... FUNCTION VALUE. */
/* G.... GRADIENT VECTOR. */
/* H.... LOWER TRIANGLE OF THE HESSIAN, STORED ROWWISE. */
/* IV... INTEGER VALUE ARRAY. */
/* LH... LENGTH OF H = P*(P+1)/2. */
/* LIV.. LENGTH OF IV (AT LEAST 59 + 3*N). */
/* LV... LENGTH OF V (AT LEAST 78 + N*(N+27)/2). */
/* N.... NUMBER OF VARIABLES (COMPONENTS IN X AND G). */
/* V.... FLOATING-POINT VALUE ARRAY. */
/* X.... PARAMETER VECTOR. */

/*  ***  DISCUSSION  *** */

/*        PARAMETERS IV, N, V, AND X ARE THE SAME AS THE CORRESPONDING */
/*     ONES TO  DMNHB (WHICH SEE), EXCEPT THAT V CAN BE SHORTER (SINCE */
/*     THE PART OF V THAT  DMNHB USES FOR STORING G AND H IS NOT NEEDED). */
/*     MOREOVER, COMPARED WITH  DMNHB, IV(1) MAY HAVE THE TWO ADDITIONAL */
/*     OUTPUT VALUES 1 AND 2, WHICH ARE EXPLAINED BELOW, AS IS THE USE */
/*     OF IV(TOOBIG) AND IV(NFGCAL).  THE VALUE IV(G), WHICH IS AN */
/*     OUTPUT VALUE FROM  DMNHB, IS NOT REFERENCED BY DRMNHB OR THE */
/*     SUBROUTINES IT CALLS. */

/* IV(1) = 1 MEANS THE CALLER SHOULD SET FX TO F(X), THE FUNCTION VALUE */
/*             AT X, AND CALL DRMNHB AGAIN, HAVING CHANGED NONE OF THE */
/*             OTHER PARAMETERS.  AN EXCEPTION OCCURS IF F(X) CANNOT BE */
/*             COMPUTED (E.G. IF OVERFLOW WOULD OCCUR), WHICH MAY HAPPEN */
/*             BECAUSE OF AN OVERSIZED STEP.  IN THIS CASE THE CALLER */
/*             SHOULD SET IV(TOOBIG) = IV(2) TO 1, WHICH WILL CAUSE */
/*             DRMNHB TO IGNORE FX AND TRY A SMALLER STEP.  THE PARA- */
/*             METER NF THAT  DMNH PASSES TO CALCF (FOR POSSIBLE USE BY */
/*             CALCGH) IS A COPY OF IV(NFCALL) = IV(6). */
/* IV(1) = 2 MEANS THE CALLER SHOULD SET G TO G(X), THE GRADIENT OF F AT */
/*             X, AND H TO THE LOWER TRIANGLE OF H(X), THE HESSIAN OF F */
/*             AT X, AND CALL DRMNHB AGAIN, HAVING CHANGED NONE OF THE */
/*             OTHER PARAMETERS EXCEPT PERHAPS THE SCALE VECTOR D. */
/*                  THE PARAMETER NF THAT  DMNHB PASSES TO CALCG IS */
/*             IV(NFGCAL) = IV(7).  IF G(X) AND H(X) CANNOT BE EVALUATED, */
/*             THEN THE CALLER MAY SET IV(NFGCAL) TO 0, IN WHICH CASE */
/*             DRMNHB WILL RETURN WITH IV(1) = 65. */
/*                  NOTE -- DRMNHB OVERWRITES H WITH THE LOWER TRIANGLE */
/*             OF  DIAG(D)**-1 * H(X) * DIAG(D)**-1. */
/* . */
/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY (WINTER, SPRING 1983). */

/*        (SEE  DMNG AND  DMNH FOR REFERENCES.) */

/* +++++++++++++++++++++++++++  DECLARATIONS  ++++++++++++++++++++++++++++ */

/*  ***  LOCAL VARIABLES  *** */


/*     ***  CONSTANTS  *** */


/*  ***  NO INTRINSIC FUNCTIONS  *** */

/*  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  *** */


/* DA7SST.... ASSESSES CANDIDATE STEP. */
/* DIVSET.... PROVIDES DEFAULT IV AND V INPUT VALUES. */
/* DD7TPR... RETURNS INNER PRODUCT OF TWO VECTORS. */
/* DD7DUP.... UPDATES SCALE VECTOR D. */
/* DG7QSB... COMPUTES APPROXIMATE OPTIMAL BOUNDED STEP. */
/* I7PNVR... INVERTS PERMUTATION ARRAY. */
/* DITSUM.... PRINTS ITERATION SUMMARY AND INFO ON INITIAL AND FINAL X. */
/* DPARCK.... CHECKS VALIDITY OF INPUT IV AND V VALUES. */
/* DRLDST... COMPUTES V(RELDX) = RELATIVE STEP SIZE. */
/* DS7IPR... APPLIES PERMUTATION TO LOWER TRIANG. OF SYM. MATRIX. */
/* DS7LVM... MULTIPLIES SYMMETRIC MATRIX TIMES VECTOR, GIVEN THE LOWER */
/*             TRIANGLE OF THE MATRIX. */
/* STOPX.... RETURNS .TRUE. IF THE BREAK KEY HAS BEEN PRESSED. */
/* DV2NRM... RETURNS THE 2-NORM OF A VECTOR. */
/* DV2AXY.... COMPUTES SCALAR TIMES ONE VECTOR PLUS ANOTHER. */
/* DV7CPY.... COPIES ONE VECTOR TO ANOTHER. */
/* DV7IPR... APPLIES PERMUTATION TO VECTOR. */
/* DV7SCP... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR. */
/* DV7VMP... MULTIPLIES (OR DIVIDES) TWO VECTORS COMPONENTWISE. */

/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/*  ***  IV SUBSCRIPT VALUES  *** */

/*  ***  (NOTE THAT NC AND N0 ARE STORED IN IV(G0) AND IV(STLSTG) RESP.) */

/* /6 */
/*     DATA CNVCOD/55/, DG/37/, DTOL/59/, DTYPE/16/, IRC/29/, IVNEED/3/, */
/*    1     KAGQT/33/, LMAT/42/, MODE/35/, MODEL/5/, MXFCAL/17/, */
/*    2     MXITER/18/, N0/41/, NC/48/, NEXTIV/46/, NEXTV/47/, NFCALL/6/, */
/*    3     NFGCAL/7/, NGCALL/30/, NITER/31/, PERM/58/, RADINC/8/, */
/*    4     RESTOR/9/, STEP/40/, STGLIM/11/, TOOBIG/2/, VNEED/4/, W/34/, */
/*    5     XIRC/13/, X0/43/ */
/* /7 */
/* / */

/*  ***  V SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA DGNORM/1/, DINIT/38/, DSTNRM/2/, DTINIT/39/, D0INIT/40/, */
/*    1     F/10/, F0/13/, FDIF/11/, GTSTEP/4/, INCFAC/23/, LMAX0/35/, */
/*    2     LMAXS/36/, PHMXFC/21/, PREDUC/7/, RADFAC/16/, RADIUS/8/, */
/*    3     RAD0/9/, RELDX/17/, STPPAR/5/, TUNER4/29/, TUNER5/30/ */
/* /7 */
/* / */

/* /6 */
/*     DATA NEGONE/-1.D+0/, ONE/1.D+0/, ONEP2/1.2D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --h__;
    --iv;
    --v;
    --x;
    --g;
    --d__;
    b -= 3;

    /* Function Body */
    i__ = iv[1];
    if (i__ == 1) {
	goto L50;
    }
    if (i__ == 2) {
	goto L60;
    }

/*  ***  CHECK VALIDITY OF IV AND V INPUT VALUES  *** */

    if (iv[1] == 0) {
	divset_(&c__2, &iv[1], liv, lv, &v[1]);
    }
    if (iv[1] < 12) {
	goto L10;
    }
    if (iv[1] > 13) {
	goto L10;
    }
    iv[4] = iv[4] + *n * (*n + 27) / 2 + 7;
    iv[3] += *n * 3;
L10:
    dparck_(&c__2, &d__[1], &iv[1], liv, lv, n, &v[1]);
    i__ = iv[1] - 2;
    if (i__ > 12) {
	goto L999;
    }
    nn1o2 = *n * (*n + 1) / 2;
    if (*lh >= nn1o2) {
	switch (i__) {
	    case 1:  goto L250;
	    case 2:  goto L250;
	    case 3:  goto L250;
	    case 4:  goto L250;
	    case 5:  goto L250;
	    case 6:  goto L250;
	    case 7:  goto L190;
	    case 8:  goto L150;
	    case 9:  goto L190;
	    case 10:  goto L20;
	    case 11:  goto L20;
	    case 12:  goto L30;
	}
    }
    iv[1] = 81;
    goto L440;

/*  ***  STORAGE ALLOCATION  *** */

L20:
    iv[59] = iv[42] + nn1o2;
    iv[43] = iv[59] + (*n << 1);
    iv[40] = iv[43] + (*n << 1);
    iv[37] = iv[40] + *n * 3;
    iv[34] = iv[37] + (*n << 1);
    iv[47] = iv[34] + (*n << 2) + 7;
    iv[46] = iv[58] + *n * 3;
    if (iv[1] != 13) {
	goto L30;
    }
    iv[1] = 14;
    goto L999;

/*  ***  INITIALIZATION  *** */

L30:
    iv[31] = 0;
    iv[6] = 1;
    iv[30] = 1;
    iv[7] = 1;
    iv[35] = -1;
    iv[5] = 1;
    iv[11] = 1;
    iv[2] = 0;
    iv[55] = 0;
    iv[8] = 0;
    iv[48] = *n;
    v[9] = 0.;
    v[5] = 0.;
    if (v[38] >= 0.) {
	dv7scp_(n, &d__[1], &v[38]);
    }
    k = iv[59];
    if (v[39] > 0.) {
	dv7scp_(n, &v[k], &v[39]);
    }
    k += *n;
    if (v[40] > 0.) {
	dv7scp_(n, &v[k], &v[40]);
    }

/*  ***  CHECK CONSISTENCY OF B AND INITIALIZE IP ARRAY  *** */

    ipi = iv[58];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iv[ipi] = i__;
	++ipi;
	if (b[(i__ << 1) + 1] > b[(i__ << 1) + 2]) {
	    goto L420;
	}
/* L40: */
    }

/*  ***  GET INITIAL FUNCTION VALUE  *** */

    iv[1] = 1;
    goto L450;

L50:
    v[10] = *fx;
    if (iv[35] >= 0) {
	goto L250;
    }
    v[13] = *fx;
    iv[1] = 2;
    if (iv[2] == 0) {
	goto L999;
    }
    iv[1] = 63;
    goto L440;

/*  ***  MAKE SURE GRADIENT COULD BE COMPUTED  *** */

L60:
    if (iv[2] == 0) {
	goto L70;
    }
    iv[1] = 65;
    goto L440;

/*  ***  UPDATE THE SCALE VECTOR D  *** */

L70:
    dg1 = iv[37];
    if (iv[16] <= 0) {
	goto L90;
    }
    k = dg1;
    j = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j += i__;
	v[k] = h__[j];
	++k;
/* L80: */
    }
    dd7dup_(&d__[1], &v[dg1], &iv[1], liv, lv, n, &v[1]);

/*  ***  COMPUTE SCALED GRADIENT AND ITS NORM  *** */

L90:
    dg1 = iv[37];
    dv7vmp_(n, &v[dg1], &g[1], &d__[1], &c_n1);

/*  ***  COMPUTE SCALED HESSIAN  *** */

    k = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = 1. / d__[i__];
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    h__[k] = t * h__[k] / d__[j];
	    ++k;
/* L100: */
	}
/* L110: */
    }

/*  ***  CHOOSE INITIAL PERMUTATION  *** */

    ipi = iv[58];
    ipn = ipi + *n;
    ipiv2 = ipn - 1;
/*     *** INVERT OLD PERMUTATION ARRAY *** */
    i7pnvr_(n, &iv[ipn], &iv[ipi]);
    k = iv[48];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (b[(i__ << 1) + 1] >= b[(i__ << 1) + 2]) {
	    goto L120;
	}
	xi = x[i__];
	gi = g[i__];
	if (xi <= b[(i__ << 1) + 1] && gi > 0.) {
	    goto L120;
	}
	if (xi >= b[(i__ << 1) + 2] && gi < 0.) {
	    goto L120;
	}
	iv[ipi] = i__;
	++ipi;
	j = ipiv2 + i__;
/*           *** DISALLOW CONVERGENCE IF X(I) HAS JUST BEEN FREED *** */
	if (iv[j] > k) {
	    iv[55] = 0;
	}
	goto L130;
L120:
	--ipn;
	iv[ipn] = i__;
L130:
	;
    }
    iv[48] = ipn - iv[58];

/*  ***  PERMUTE SCALED GRADIENT AND HESSIAN ACCORDINGLY  *** */

    ipi = iv[58];
    ds7ipr_(n, &iv[ipi], &h__[1]);
    dv7ipr_(n, &iv[ipi], &v[dg1]);
    v[1] = 0.;
    if (iv[48] > 0) {
	v[1] = dv2nrm_(&iv[48], &v[dg1]);
    }

    if (iv[55] != 0) {
	goto L430;
    }
    if (iv[35] == 0) {
	goto L380;
    }

/*  ***  ALLOW FIRST STEP TO HAVE SCALED 2-NORM AT MOST V(LMAX0)  *** */

    v[8] = v[35] / (v[21] + 1.);

    iv[35] = 0;


/* -----------------------------  MAIN LOOP  ----------------------------- */


/*  ***  PRINT ITERATION SUMMARY, CHECK ITERATION LIMIT  *** */

L140:
    ditsum_(&d__[1], &g[1], &iv[1], liv, lv, n, &v[1], &x[1]);
L150:
    k = iv[31];
    if (k < iv[18]) {
	goto L160;
    }
    iv[1] = 10;
    goto L440;

L160:
    iv[31] = k + 1;

/*  ***  INITIALIZE FOR START OF NEXT ITERATION  *** */

    x01 = iv[43];
    v[13] = v[10];
    iv[29] = 4;
    iv[33] = -1;

/*     ***  COPY X TO X0  *** */

    dv7cpy_(n, &v[x01], &x[1]);

/*  ***  UPDATE RADIUS  *** */

    if (k == 0) {
	goto L180;
    }
    step1 = iv[40];
    k = step1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[k] = d__[i__] * v[k];
	++k;
/* L170: */
    }
    t = v[16] * dv2nrm_(n, &v[step1]);
    if (v[16] < 1. || t > v[8]) {
	v[8] = t;
    }

/*  ***  CHECK STOPX AND FUNCTION EVALUATION LIMIT  *** */

L180:
    if (! stopx_(&dummy)) {
	goto L200;
    }
    iv[1] = 11;
    goto L210;

/*     ***  COME HERE WHEN RESTARTING AFTER FUNC. EVAL. LIMIT OR STOPX. */

L190:
    if (v[10] >= v[13]) {
	goto L200;
    }
    v[16] = 1.;
    k = iv[31];
    goto L160;

L200:
    if (iv[6] < iv[17]) {
	goto L220;
    }
    iv[1] = 9;
L210:
    if (v[10] >= v[13]) {
	goto L440;
    }

/*        ***  IN CASE OF STOPX OR FUNCTION EVALUATION LIMIT WITH */
/*        ***  IMPROVED V(F), EVALUATE THE GRADIENT AT X. */

    iv[55] = iv[1];
    goto L370;

/* . . . . . . . . . . . . .  COMPUTE CANDIDATE STEP  . . . . . . . . . . */

L220:
    step1 = iv[40];
    l = iv[42];
    w1 = iv[34];
    ipi = iv[58];
    ipn = ipi + *n;
    ipiv2 = ipn + *n;
    tg1 = iv[37];
    td1 = tg1 + *n;
    x01 = iv[43];
    x11 = x01 + *n;
    dg7qsb_(&b[3], &d__[1], &h__[1], &g[1], &iv[ipi], &iv[ipn], &iv[ipiv2], &
	    iv[33], &v[l], lv, n, &iv[41], &iv[48], &v[step1], &v[td1], &v[
	    tg1], &v[1], &v[w1], &v[x11], &v[x01]);
    if (iv[29] != 6) {
	goto L230;
    }
    if (iv[9] != 2) {
	goto L250;
    }
    rstrst = 2;
    goto L260;

/*  ***  CHECK WHETHER EVALUATING F(X0 + STEP) LOOKS WORTHWHILE  *** */

L230:
    iv[2] = 0;
    if (v[2] <= 0.) {
	goto L250;
    }
    if (iv[29] != 5) {
	goto L240;
    }
    if (v[16] <= 1.) {
	goto L240;
    }
    if (v[7] > v[11] * 1.2) {
	goto L240;
    }
    if (iv[9] != 2) {
	goto L250;
    }
    rstrst = 0;
    goto L260;

/*  ***  COMPUTE F(X0 + STEP)  *** */

L240:
    dv2axy_(n, &x[1], &c_b38, &v[step1], &v[x01]);
    ++iv[6];
    iv[1] = 1;
    goto L450;

/* . . . . . . . . . . . . .  ASSESS CANDIDATE STEP  . . . . . . . . . . . */

L250:
    rstrst = 3;
L260:
    x01 = iv[43];
    v[17] = drldst_(n, &d__[1], &x[1], &v[x01]);
    da7sst_(&iv[1], liv, lv, &v[1]);
    step1 = iv[40];
    lstgst = step1 + (*n << 1);
    i__ = iv[9] + 1;
    switch (i__) {
	case 1:  goto L300;
	case 2:  goto L270;
	case 3:  goto L280;
	case 4:  goto L290;
    }
L270:
    dv7cpy_(n, &x[1], &v[x01]);
    goto L300;
L280:
    dv7cpy_(n, &v[lstgst], &x[1]);
    goto L300;
L290:
    dv7cpy_(n, &x[1], &v[lstgst]);
    dv2axy_(n, &v[step1], &c_b43, &v[x01], &x[1]);
    v[17] = drldst_(n, &d__[1], &x[1], &v[x01]);
    iv[9] = rstrst;

L300:
    k = iv[29];
    switch (k) {
	case 1:  goto L310;
	case 2:  goto L340;
	case 3:  goto L340;
	case 4:  goto L340;
	case 5:  goto L310;
	case 6:  goto L320;
	case 7:  goto L330;
	case 8:  goto L330;
	case 9:  goto L330;
	case 10:  goto L330;
	case 11:  goto L330;
	case 12:  goto L330;
	case 13:  goto L410;
	case 14:  goto L380;
    }

/*     ***  RECOMPUTE STEP WITH NEW RADIUS  *** */

L310:
    v[8] = v[16] * v[2];
    goto L180;

/*  ***  COMPUTE STEP OF LENGTH V(LMAXS) FOR SINGULAR CONVERGENCE TEST. */

L320:
    v[8] = v[36];
    goto L220;

/*  ***  CONVERGENCE OR FALSE CONVERGENCE  *** */

L330:
    iv[55] = k - 4;
    if (v[10] >= v[13]) {
	goto L430;
    }
    if (iv[13] == 14) {
	goto L430;
    }
    iv[13] = 14;

/* . . . . . . . . . . . .  PROCESS ACCEPTABLE STEP  . . . . . . . . . . . */

L340:
    if (iv[29] != 3) {
	goto L370;
    }
    temp1 = lstgst;

/*     ***  PREPARE FOR GRADIENT TESTS  *** */
/*     ***  SET  TEMP1 = HESSIAN * STEP + G(X0) */
/*     ***             = DIAG(D) * (H * STEP + G(X0)) */

    k = temp1;
    step0 = step1 - 1;
    ipi = iv[58];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iv[ipi];
	++ipi;
	step1 = step0 + j;
	v[k] = d__[j] * v[step1];
	++k;
/* L350: */
    }
/*        USE X0 VECTOR AS TEMPORARY. */
    ds7lvm_(n, &v[x01], &h__[1], &v[temp1]);
    temp0 = temp1 - 1;
    ipi = iv[58];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iv[ipi];
	++ipi;
	temp1 = temp0 + j;
	v[temp1] = d__[j] * v[x01] + g[j];
	++x01;
/* L360: */
    }

/*  ***  COMPUTE GRADIENT AND HESSIAN  *** */

L370:
    ++iv[30];
    iv[2] = 0;
    iv[1] = 2;
    goto L450;

L380:
    iv[1] = 2;
    if (iv[29] != 3) {
	goto L140;
    }

/*  ***  SET V(RADFAC) BY GRADIENT TESTS  *** */

    step1 = iv[40];
/*     *** TEMP1 = STLSTG *** */
    temp1 = step1 + (*n << 1);

/*     ***  SET  TEMP1 = DIAG(D)**-1 * (HESSIAN*STEP + (G(X0)-G(X)))  *** */

    k = temp1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v[k] = (v[k] - g[i__]) / d__[i__];
	++k;
/* L390: */
    }

/*     ***  DO GRADIENT TESTS  *** */

    if (dv2nrm_(n, &v[temp1]) <= v[1] * v[29]) {
	goto L400;
    }
    if (dd7tpr_(n, &g[1], &v[step1]) >= v[4] * v[30]) {
	goto L140;
    }
L400:
    v[16] = v[23];
    goto L140;

/* . . . . . . . . . . . . . .  MISC. DETAILS  . . . . . . . . . . . . . . */

/*  ***  BAD PARAMETERS TO ASSESS  *** */

L410:
    iv[1] = 64;
    goto L440;

/*  ***  INCONSISTENT B  *** */

L420:
    iv[1] = 82;
    goto L440;

/*  ***  PRINT SUMMARY OF FINAL ITERATION AND OTHER REQUESTED ITEMS  *** */

L430:
    iv[1] = iv[55];
    iv[55] = 0;
L440:
    ditsum_(&d__[1], &g[1], &iv[1], liv, lv, n, &v[1], &x[1]);
    goto L999;

/*  ***  PROJECT X INTO FEASIBLE REGION (PRIOR TO COMPUTING F OR G)  *** */

L450:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (x[i__] < b[(i__ << 1) + 1]) {
	    x[i__] = b[(i__ << 1) + 1];
	}
	if (x[i__] > b[(i__ << 1) + 2]) {
	    x[i__] = b[(i__ << 1) + 2];
	}
/* L460: */
    }

L999:
    return 0;

/*  ***  LAST CARD OF DRMNHB FOLLOWS  *** */
} /* drmnhb_ */

