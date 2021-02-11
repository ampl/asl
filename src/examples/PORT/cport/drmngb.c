/* drmngb.f -- translated by f2c (version 20160102).
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
static doublereal c_b16 = 0.;
static logical c_false = FALSE_;
static integer c_n1 = -1;
static doublereal c_b40 = 1.;
static doublereal c_b45 = -1.;

/* Subroutine */ int drmngb_(doublereal *b, doublereal *d__, doublereal *fx, 
	doublereal *g, integer *iv, integer *liv, integer *lv, integer *n, 
	doublereal *v, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal t;
    static integer z__, i1, n1, w1, g01;
    static doublereal gi;
    static integer x01;
    static doublereal xi;
    static integer dg1, td1, tg1, np1, ipi, ipn, temp0, temp1, step1, dummy;
    extern /* Subroutine */ int dd7dgb_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern logical stopx_(integer *);
    extern /* Subroutine */ int dl7upd_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *);
    static integer dstep1;
    extern /* Subroutine */ int dw7zbf_(doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *);
    extern doublereal dd7tpr_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int da7sst_(integer *, integer *, integer *, 
	    doublereal *), i7shft_(integer *, integer *, integer *), dl7vml_(
	    integer *, doublereal *, doublereal *, doublereal *);
    extern doublereal dv2nrm_(integer *, doublereal *);
    extern /* Subroutine */ int dq7rsh_(integer *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *), dv7scp_(integer *, 
	    doublereal *, doublereal *), dv7ipr_(integer *, integer *, 
	    doublereal *), dv7cpy_(integer *, doublereal *, doublereal *), 
	    dl7tvm_(integer *, doublereal *, doublereal *, doublereal *), 
	    dv2axy_(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dv7vmp_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    static integer nwtst1;
    extern /* Subroutine */ int dparck_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *);
    extern doublereal drldst_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    extern /* Subroutine */ int divset_(integer *, integer *, integer *, 
	    integer *, doublereal *), ditsum_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *);
    static integer lstgst, rstrst;


/*  ***  CARRY OUT  DMNGB (SIMPLY BOUNDED MINIMIZATION) ITERATIONS, */
/*  ***  USING DOUBLE-DOGLEG/BFGS STEPS. */

/*  ***  PARAMETER DECLARATIONS  *** */


/* --------------------------  PARAMETER USAGE  -------------------------- */

/* B.... VECTOR OF LOWER AND UPPER BOUNDS ON X. */
/* D.... SCALE VECTOR. */
/* FX... FUNCTION VALUE. */
/* G.... GRADIENT VECTOR. */
/* IV... INTEGER VALUE ARRAY. */
/* LIV.. LENGTH OF IV (AT LEAST 59) + N. */
/* LV... LENGTH OF V (AT LEAST 71 + N*(N+19)/2). */
/* N.... NUMBER OF VARIABLES (COMPONENTS IN X AND G). */
/* V.... FLOATING-POINT VALUE ARRAY. */
/* X.... VECTOR OF PARAMETERS TO BE OPTIMIZED. */

/*  ***  DISCUSSION  *** */

/*        PARAMETERS IV, N, V, AND X ARE THE SAME AS THE CORRESPONDING */
/*     ONES TO  DMNGB (WHICH SEE), EXCEPT THAT V CAN BE SHORTER (SINCE */
/*     THE PART OF V THAT  DMNGB USES FOR STORING G IS NOT NEEDED). */
/*     MOREOVER, COMPARED WITH  DMNGB, IV(1) MAY HAVE THE TWO ADDITIONAL */
/*     OUTPUT VALUES 1 AND 2, WHICH ARE EXPLAINED BELOW, AS IS THE USE */
/*     OF IV(TOOBIG) AND IV(NFGCAL).  THE VALUE IV(G), WHICH IS AN */
/*     OUTPUT VALUE FROM  DMNGB (AND SMSNOB), IS NOT REFERENCED BY */
/*     DRMNGB OR THE SUBROUTINES IT CALLS. */
/*        FX AND G NEED NOT HAVE BEEN INITIALIZED WHEN DRMNGB IS CALLED */
/*     WITH IV(1) = 12, 13, OR 14. */

/* IV(1) = 1 MEANS THE CALLER SHOULD SET FX TO F(X), THE FUNCTION VALUE */
/*             AT X, AND CALL DRMNGB AGAIN, HAVING CHANGED NONE OF THE */
/*             OTHER PARAMETERS.  AN EXCEPTION OCCURS IF F(X) CANNOT BE */
/*             (E.G. IF OVERFLOW WOULD OCCUR), WHICH MAY HAPPEN BECAUSE */
/*             OF AN OVERSIZED STEP.  IN THIS CASE THE CALLER SHOULD SET */
/*             IV(TOOBIG) = IV(2) TO 1, WHICH WILL CAUSE DRMNGB TO IG- */
/*             NORE FX AND TRY A SMALLER STEP.  THE PARAMETER NF THAT */
/*              DMNGB PASSES TO CALCF (FOR POSSIBLE USE BY CALCG) IS A */
/*             COPY OF IV(NFCALL) = IV(6). */
/* IV(1) = 2 MEANS THE CALLER SHOULD SET G TO G(X), THE GRADIENT VECTOR */
/*             OF F AT X, AND CALL DRMNGB AGAIN, HAVING CHANGED NONE OF */
/*             THE OTHER PARAMETERS EXCEPT POSSIBLY THE SCALE VECTOR D */
/*             WHEN IV(DTYPE) = 0.  THE PARAMETER NF THAT  DMNGB PASSES */
/*             TO CALCG IS IV(NFGCAL) = IV(7).  IF G(X) CANNOT BE */
/*             EVALUATED, THEN THE CALLER MAY SET IV(NFGCAL) TO 0, IN */
/*             WHICH CASE DRMNGB WILL RETURN WITH IV(1) = 65. */
/* . */
/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY (DECEMBER 1979).  REVISED SEPT. 1982. */
/*     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH SUPPORTED */
/*     IN PART BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS */
/*     MCS-7600324 AND MCS-7906671. */

/*        (SEE  DMNG FOR REFERENCES.) */

/* +++++++++++++++++++++++++++  DECLARATIONS  ++++++++++++++++++++++++++++ */

/*  ***  LOCAL VARIABLES  *** */


/*     ***  CONSTANTS  *** */


/*  ***  NO INTRINSIC FUNCTIONS  *** */

/*  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  *** */


/* DA7SST.... ASSESSES CANDIDATE STEP. */
/* DD7DGB... COMPUTES SIMPLY BOUNDED DOUBLE-DOGLEG (CANDIDATE) STEP. */
/* DIVSET.... SUPPLIES DEFAULT IV AND V INPUT COMPONENTS. */
/* DD7TPR... RETURNS INNER PRODUCT OF TWO VECTORS. */
/* I7SHFT... CYCLICALLLY SHIFTS AN ARRAY OF INTEGERS. */
/* DITSUM.... PRINTS ITERATION SUMMARY AND INFO ON INITIAL AND FINAL X. */
/* DL7TVM... MULTIPLIES TRANSPOSE OF LOWER TRIANGLE TIMES VECTOR. */
/* LUPDT.... UPDATES CHOLESKY FACTOR OF HESSIAN APPROXIMATION. */
/* DL7VML.... MULTIPLIES LOWER TRIANGLE TIMES VECTOR. */
/* DPARCK.... CHECKS VALIDITY OF INPUT IV AND V VALUES. */
/* DQ7RSH... CYCLICALLY SHIFTS CHOLESKY FACTOR. */
/* DRLDST... COMPUTES V(RELDX) = RELATIVE STEP SIZE. */
/* STOPX.... RETURNS .TRUE. IF THE BREAK KEY HAS BEEN PRESSED. */
/* DV2NRM... RETURNS THE 2-NORM OF A VECTOR. */
/* DV2AXY.... COMPUTES SCALAR TIMES ONE VECTOR PLUS ANOTHER. */
/* DV7CPY.... COPIES ONE VECTOR TO ANOTHER. */
/* DV7IPR... CYCLICALLY SHIFTS A FLOATING-POINT ARRAY. */
/* DV7SCP... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR. */
/* DV7VMP... MULTIPLIES VECTOR BY VECTOR RAISED TO POWER (COMPONENTWISE). */
/* DW7ZBF... COMPUTES W AND Z FOR DL7UPD CORRESPONDING TO BFGS UPDATE. */

/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/*  ***  IV SUBSCRIPT VALUES  *** */

/*  ***  (NOTE THAT NC IS STORED IN IV(G0)) *** */

/* /6 */
/*     DATA CNVCOD/55/, DG/37/, INITH/25/, IRC/29/, IVNEED/3/, KAGQT/33/, */
/*    1     MODE/35/, MODEL/5/, MXFCAL/17/, MXITER/18/, NC/48/, */
/*    2     NEXTIV/46/, NEXTV/47/, NFCALL/6/, NFGCAL/7/, NGCALL/30/, */
/*    3     NITER/31/, NWTSTP/34/, PERM/58/, RADINC/8/, RESTOR/9/, */
/*    4     STEP/40/, STGLIM/11/, STLSTG/41/, TOOBIG/2/, XIRC/13/, X0/43/ */
/* /7 */
/* / */

/*  ***  V SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA DGNORM/1/, DINIT/38/, DSTNRM/2/, F/10/, F0/13/, FDIF/11/, */
/*    1     GTSTEP/4/, INCFAC/23/, LMAT/42/, LMAX0/35/, LMAXS/36/, */
/*    2     PREDUC/7/, RADFAC/16/, RADIUS/8/, RAD0/9/, RELDX/17/, */
/*    3     TUNER4/29/, TUNER5/30/, VNEED/4/ */
/* /7 */
/* / */

/* /6 */
/*     DATA NEGONE/-1.D+0/, ONE/1.D+0/, ONEP2/1.2D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --iv;
    --v;
    --x;
    --g;
    --d__;
    b -= 3;

    /* Function Body */
    i__ = iv[1];
    if (i__ == 1) {
	goto L70;
    }
    if (i__ == 2) {
	goto L80;
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
    iv[4] += *n * (*n + 19) / 2;
    iv[3] += *n;
L10:
    dparck_(&c__2, &d__[1], &iv[1], liv, lv, n, &v[1]);
    i__ = iv[1] - 2;
    if (i__ > 12) {
	goto L999;
    }
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

/*  ***  STORAGE ALLOCATION  *** */

L20:
    l = iv[42];
    iv[43] = l + *n * (*n + 1) / 2;
    iv[40] = iv[43] + (*n << 1);
    iv[41] = iv[40] + (*n << 1);
    iv[34] = iv[41] + *n;
    iv[37] = iv[34] + (*n << 1);
    iv[47] = iv[37] + (*n << 1);
    iv[46] = iv[58] + *n;
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

/*  ***  CHECK CONSISTENCY OF B AND INITIALIZE IP ARRAY  *** */

    ipi = iv[58];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iv[ipi] = i__;
	++ipi;
	if (b[(i__ << 1) + 1] > b[(i__ << 1) + 2]) {
	    goto L410;
	}
/* L40: */
    }

    if (v[38] >= 0.) {
	dv7scp_(n, &d__[1], &v[38]);
    }
    if (iv[25] != 1) {
	goto L60;
    }

/*     ***  SET THE INITIAL HESSIAN APPROXIMATION TO DIAG(D)**-2  *** */

    l = iv[42];
    i__1 = *n * (*n + 1) / 2;
    dv7scp_(&i__1, &v[l], &c_b16);
    k = l - 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k += i__;
	t = d__[i__];
	if (t <= 0.) {
	    t = 1.;
	}
	v[k] = t;
/* L50: */
    }

/*  ***  GET INITIAL FUNCTION VALUE  *** */

L60:
    iv[1] = 1;
    goto L440;

L70:
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
    goto L430;

/*  ***  MAKE SURE GRADIENT COULD BE COMPUTED  *** */

L80:
    if (iv[2] == 0) {
	goto L90;
    }
    iv[1] = 65;
    goto L430;

/*  ***  CHOOSE INITIAL PERMUTATION  *** */

L90:
    ipi = iv[58];
    ipn = ipi + *n;
    n1 = *n;
    np1 = *n + 1;
    l = iv[42];
    w1 = iv[34] + *n;
    k = *n - iv[48];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--ipn;
	j = iv[ipn];
	if (b[(j << 1) + 1] >= b[(j << 1) + 2]) {
	    goto L100;
	}
	xi = x[j];
	gi = g[j];
	if (xi <= b[(j << 1) + 1] && gi > 0.) {
	    goto L100;
	}
	if (xi >= b[(j << 1) + 2] && gi < 0.) {
	    goto L100;
	}
/*           *** DISALLOW CONVERGENCE IF X(J) HAS JUST BEEN FREED *** */
	if (i__ <= k) {
	    iv[55] = 0;
	}
	goto L120;
L100:
	i1 = np1 - i__;
	if (i1 >= n1) {
	    goto L110;
	}
	i7shft_(&n1, &i1, &iv[ipi]);
	dq7rsh_(&i1, &n1, &c_false, &g[1], &v[l], &v[w1]);
L110:
	--n1;
L120:
	;
    }

    iv[48] = n1;
    v[1] = 0.;
    if (n1 <= 0) {
	goto L130;
    }
    dg1 = iv[37];
    dv7vmp_(n, &v[dg1], &g[1], &d__[1], &c_n1);
    dv7ipr_(n, &iv[ipi], &v[dg1]);
    v[1] = dv2nrm_(&n1, &v[dg1]);
L130:
    if (iv[55] != 0) {
	goto L420;
    }
    if (iv[35] == 0) {
	goto L370;
    }

/*  ***  ALLOW FIRST STEP TO HAVE SCALED 2-NORM AT MOST V(LMAX0)  *** */

    v[8] = v[35];

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
    goto L430;

/*  ***  UPDATE RADIUS  *** */

L160:
    iv[31] = k + 1;
    if (k == 0) {
	goto L170;
    }
    t = v[16] * v[2];
    if (v[16] < 1. || t > v[8]) {
	v[8] = t;
    }

/*  ***  INITIALIZE FOR START OF NEXT ITERATION  *** */

L170:
    x01 = iv[43];
    v[13] = v[10];
    iv[29] = 4;
    iv[33] = -1;

/*     ***  COPY X TO X0  *** */

    dv7cpy_(n, &v[x01], &x[1]);

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
	goto L430;
    }

/*        ***  IN CASE OF STOPX OR FUNCTION EVALUATION LIMIT WITH */
/*        ***  IMPROVED V(F), EVALUATE THE GRADIENT AT X. */

    iv[55] = iv[1];
    goto L360;

/* . . . . . . . . . . . . .  COMPUTE CANDIDATE STEP  . . . . . . . . . . */

L220:
    step1 = iv[40];
    dg1 = iv[37];
    nwtst1 = iv[34];
    w1 = nwtst1 + *n;
    dstep1 = step1 + *n;
    ipi = iv[58];
    l = iv[42];
    tg1 = dg1 + *n;
    x01 = iv[43];
    td1 = x01 + *n;
    dd7dgb_(&b[3], &d__[1], &v[dg1], &v[dstep1], &g[1], &iv[ipi], &iv[33], &v[
	    l], lv, n, &iv[48], &v[nwtst1], &v[step1], &v[td1], &v[tg1], &v[1]
	    , &v[w1], &v[x01]);
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
    dv2axy_(n, &x[1], &c_b40, &v[step1], &v[x01]);
    ++iv[6];
    iv[1] = 1;
    goto L440;

/* . . . . . . . . . . . . .  ASSESS CANDIDATE STEP  . . . . . . . . . . . */

L250:
    rstrst = 3;
L260:
    x01 = iv[43];
    v[17] = drldst_(n, &d__[1], &x[1], &v[x01]);
    da7sst_(&iv[1], liv, lv, &v[1]);
    step1 = iv[40];
    lstgst = iv[41];
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
    dv2axy_(n, &v[step1], &c_b45, &v[x01], &x[1]);
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
	case 13:  goto L400;
	case 14:  goto L370;
    }

/*     ***  RECOMPUTE STEP WITH CHANGED RADIUS  *** */

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
	goto L420;
    }
    if (iv[13] == 14) {
	goto L420;
    }
    iv[13] = 14;

/* . . . . . . . . . . . .  PROCESS ACCEPTABLE STEP  . . . . . . . . . . . */

L340:
    x01 = iv[43];
    step1 = iv[40];
    dv2axy_(n, &v[step1], &c_b45, &v[x01], &x[1]);
    if (iv[29] != 3) {
	goto L360;
    }

/*     ***  SET  TEMP1 = HESSIAN * STEP  FOR USE IN GRADIENT TESTS  *** */

/*     ***  USE X0 AS TEMPORARY... */

    ipi = iv[58];
    dv7cpy_(n, &v[x01], &v[step1]);
    dv7ipr_(n, &iv[ipi], &v[x01]);
    l = iv[42];
    dl7tvm_(n, &v[x01], &v[l], &v[x01]);
    dl7vml_(n, &v[x01], &v[l], &v[x01]);

/*        *** UNPERMUTE X0 INTO TEMP1 *** */

    temp1 = iv[41];
    temp0 = temp1 - 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iv[ipi];
	++ipi;
	k = temp0 + j;
	v[k] = v[x01];
	++x01;
/* L350: */
    }

/*  ***  SAVE OLD GRADIENT, COMPUTE NEW ONE  *** */

L360:
    g01 = iv[34] + *n;
    dv7cpy_(n, &v[g01], &g[1]);
    ++iv[30];
    iv[2] = 0;
    iv[1] = 2;
    goto L999;

/*  ***  INITIALIZATIONS -- G0 = G - G0, ETC.  *** */

L370:
    g01 = iv[34] + *n;
    dv2axy_(n, &v[g01], &c_b45, &v[g01], &g[1]);
    step1 = iv[40];
    temp1 = iv[41];
    if (iv[29] != 3) {
	goto L390;
    }

/*  ***  SET V(RADFAC) BY GRADIENT TESTS  *** */

/*     ***  SET  TEMP1 = DIAG(D)**-1 * (HESSIAN*STEP + (G(X0)-G(X)))  *** */

    dv2axy_(n, &v[temp1], &c_b45, &v[g01], &v[temp1]);
    dv7vmp_(n, &v[temp1], &v[temp1], &d__[1], &c_n1);

/*        ***  DO GRADIENT TESTS  *** */

    if (dv2nrm_(n, &v[temp1]) <= v[1] * v[29]) {
	goto L380;
    }
    if (dd7tpr_(n, &g[1], &v[step1]) >= v[4] * v[30]) {
	goto L390;
    }
L380:
    v[16] = v[23];

/*  ***  UPDATE H, LOOP  *** */

L390:
    w1 = iv[34];
    z__ = iv[43];
    l = iv[42];
    ipi = iv[58];
    dv7ipr_(n, &iv[ipi], &v[step1]);
    dv7ipr_(n, &iv[ipi], &v[g01]);
    dw7zbf_(&v[l], n, &v[step1], &v[w1], &v[g01], &v[z__]);

/*     ** USE THE N-VECTORS STARTING AT V(STEP1) AND V(G01) FOR SCRATCH.. */
    dl7upd_(&v[temp1], &v[step1], &v[l], &v[g01], &v[l], n, &v[w1], &v[z__]);
    iv[1] = 2;
    goto L140;

/* . . . . . . . . . . . . . .  MISC. DETAILS  . . . . . . . . . . . . . . */

/*  ***  BAD PARAMETERS TO ASSESS  *** */

L400:
    iv[1] = 64;
    goto L430;

/*  ***  INCONSISTENT B  *** */

L410:
    iv[1] = 82;
    goto L430;

/*  ***  PRINT SUMMARY OF FINAL ITERATION AND OTHER REQUESTED ITEMS  *** */

L420:
    iv[1] = iv[55];
    iv[55] = 0;
L430:
    ditsum_(&d__[1], &g[1], &iv[1], liv, lv, n, &v[1], &x[1]);
    goto L999;

/*  ***  PROJECT X INTO FEASIBLE REGION (PRIOR TO COMPUTING F OR G)  *** */

L440:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (x[i__] < b[(i__ << 1) + 1]) {
	    x[i__] = b[(i__ << 1) + 1];
	}
	if (x[i__] > b[(i__ << 1) + 2]) {
	    x[i__] = b[(i__ << 1) + 2];
	}
/* L450: */
    }

L999:
    return 0;

/*  ***  LAST CARD OF DRMNGB FOLLOWS  *** */
} /* drmngb_ */

