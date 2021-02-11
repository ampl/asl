/* drn2gb.f -- translated by f2c (version 20160102).
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
static doublereal c_b14 = 0.;
static integer c__0 = 0;
static logical c_true = TRUE_;

/* Subroutine */ int drn2gb_(doublereal *b, doublereal *d__, doublereal *dr, 
	integer *iv, integer *liv, integer *lv, integer *n, integer *nd, 
	integer *n1, integer *n2, integer *p, doublereal *r__, doublereal *rd,
	 doublereal *v, doublereal *x)
{
    /* System generated locals */
    integer dr_dim1, dr_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, l;
    static doublereal t;
    static integer g1, y1, gi, lh, nn, yi, rd1, iv1, qtr1, rmat1, jtol1;
    extern /* Subroutine */ int dg7itb_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *), dq7rad_(integer *, 
	    integer *, integer *, doublereal *, logical *, doublereal *, 
	    doublereal *, doublereal *), dd7upd_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *), dq7apl_(integer *, integer *,
	     integer *, doublereal *, doublereal *, integer *);
    extern doublereal dd7tpr_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int dl7vml_(integer *, doublereal *, doublereal *,
	     doublereal *), dv7scp_(integer *, doublereal *, doublereal *);
    extern doublereal dv2nrm_(integer *, doublereal *);
    extern /* Subroutine */ int dv7cpy_(integer *, doublereal *, doublereal *)
	    , dr7tvm_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer ivmode;
    extern /* Subroutine */ int divset_(integer *, integer *, integer *, 
	    integer *, doublereal *), ditsum_(doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *);


/*  ***  REVISED ITERATION DRIVER FOR NL2SOL WITH SIMPLE BOUNDS  *** */


/* --------------------------  PARAMETER USAGE  -------------------------- */

/* B........ BOUNDS ON X. */
/* D........ SCALE VECTOR. */
/* DR....... DERIVATIVES OF R AT X. */
/* IV....... INTEGER VALUES ARRAY. */
/* LIV...... LENGTH OF IV... LIV MUST BE AT LEAST 4*P + 82. */
/* LV....... LENGTH OF V...  LV  MUST BE AT LEAST 105 + P*(2*P+20). */
/* N........ TOTAL NUMBER OF RESIDUALS. */
/* ND....... MAX. NO. OF RESIDUALS PASSED ON ONE CALL. */
/* N1....... LOWEST  ROW INDEX FOR RESIDUALS SUPPLIED THIS TIME. */
/* N2....... HIGHEST ROW INDEX FOR RESIDUALS SUPPLIED THIS TIME. */
/* P........ NUMBER OF PARAMETERS (COMPONENTS OF X) BEING ESTIMATED. */
/* R........ RESIDUALS. */
/* V........ FLOATING-POINT VALUES ARRAY. */
/* X........ PARAMETER VECTOR BEING ESTIMATED (INPUT = INITIAL GUESS, */
/*             OUTPUT = BEST VALUE FOUND). */

/*  ***  DISCUSSION  *** */

/*     THIS ROUTINE CARRIES OUT ITERATIONS FOR SOLVING NONLINEAR */
/*  LEAST SQUARES PROBLEMS.  IT IS SIMILAR TO  DRN2G, EXCEPT THAT */
/*  THIS ROUTINE ENFORCES THE BOUNDS  B(1,I) .LE. X(I) .LE. B(2,I), */
/*  I = 1(1)P. */

/*  ***  GENERAL  *** */

/*     CODED BY DAVID M. GAY. */

/* +++++++++++++++++++++++++++++  DECLARATIONS  ++++++++++++++++++++++++++ */

/*  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  *** */


/* DIVSET.... PROVIDES DEFAULT IV AND V INPUT COMPONENTS. */
/* DD7TPR... COMPUTES INNER PRODUCT OF TWO VECTORS. */
/* DD7UPD...  UPDATES SCALE VECTOR D. */
/* DG7ITB... PERFORMS BASIC MINIMIZATION ALGORITHM. */
/* DITSUM.... PRINTS ITERATION SUMMARY, INFO ABOUT INITIAL AND FINAL X. */
/* DL7VML.... COMPUTES L * V, V = VECTOR, L = LOWER TRIANGULAR MATRIX. */
/* DQ7APL... APPLIES QR TRANSFORMATIONS STORED BY DQ7RAD. */
/* DQ7RAD.... ADDS A NEW BLOCK OF ROWS TO QR DECOMPOSITION. */
/* DR7TVM... MULT. VECTOR BY TRANS. OF UPPER TRIANG. MATRIX FROM QR FACT. */
/* DV7CPY.... COPIES ONE VECTOR TO ANOTHER. */
/* DV7SCP... SETS ALL ELEMENTS OF A VECTOR TO A SCALAR. */
/* DV2NRM... RETURNS THE 2-NORM OF A VECTOR. */


/*  ***  LOCAL VARIABLES  *** */



/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/*  ***  IV SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA DTYPE/16/, G/28/, JCN/66/, JTOL/59/, MODE/35/, NEXTV/47/, */
/*    1     NF0/68/, NF00/81/, NF1/69/, NFCALL/6/, NFCOV/52/, NFGCAL/7/, */
/*    2     QTR/77/, RDREQ/57/, RESTOR/9/, REGD/67/, RMAT/78/, TOOBIG/2/, */
/*    3     VNEED/4/ */
/* /7 */
/* / */

/*  ***  V SUBSCRIPT VALUES  *** */

/* /6 */
/*     DATA DINIT/38/, DTINIT/39/, D0INIT/40/, F/10/, RLIMIT/46/ */
/* /7 */
/* / */
/* /6 */
/*     DATA HALF/0.5D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    /* Parameter adjustments */
    --iv;
    --v;
    --rd;
    --r__;
    --x;
    dr_dim1 = *nd;
    dr_offset = 1 + dr_dim1;
    dr -= dr_offset;
    --d__;
    b -= 3;

    /* Function Body */
    lh = *p * (*p + 1) / 2;
    if (iv[1] == 0) {
	divset_(&c__1, &iv[1], liv, lv, &v[1]);
    }
    iv1 = iv[1];
    if (iv1 > 2) {
	goto L10;
    }
    nn = *n2 - *n1 + 1;
    iv[9] = 0;
    i__ = iv1 + 4;
    if (iv[2] == 0) {
	switch (i__) {
	    case 1:  goto L150;
	    case 2:  goto L130;
	    case 3:  goto L150;
	    case 4:  goto L120;
	    case 5:  goto L120;
	    case 6:  goto L150;
	}
    }
    if (i__ != 5) {
	iv[1] = 2;
    }
    goto L40;

/*  ***  FRESH START OR RESTART -- CHECK INPUT INTEGERS  *** */

L10:
    if (*nd <= 0) {
	goto L220;
    }
    if (*p <= 0) {
	goto L220;
    }
    if (*n <= 0) {
	goto L220;
    }
    if (iv1 == 14) {
	goto L30;
    }
    if (iv1 > 16) {
	goto L270;
    }
    if (iv1 < 12) {
	goto L40;
    }
    if (iv1 == 12) {
	iv[1] = 13;
    }
    if (iv[1] != 13) {
	goto L20;
    }
    iv[4] += *p * (*p + 15) / 2;
L20:
    dg7itb_(&b[3], &d__[1], &x[1], &iv[1], liv, lv, p, p, &v[1], &x[1], &x[1])
	    ;
    if (iv[1] != 14) {
	goto L999;
    }

/*  ***  STORAGE ALLOCATION  *** */

    iv[28] = iv[47];
    iv[66] = iv[28] + (*p << 1);
    iv[78] = iv[66] + *p;
    iv[77] = iv[78] + lh;
    iv[59] = iv[77] + (*p << 1);
    iv[47] = iv[59] + (*p << 1);
/*  ***  TURN OFF COVARIANCE COMPUTATION  *** */
    iv[57] = 0;
    if (iv1 == 13) {
	goto L999;
    }

L30:
    jtol1 = iv[59];
    if (v[38] >= 0.) {
	dv7scp_(p, &d__[1], &v[38]);
    }
    if (v[39] > 0.) {
	dv7scp_(p, &v[jtol1], &v[39]);
    }
    i__ = jtol1 + *p;
    if (v[40] > 0.) {
	dv7scp_(p, &v[i__], &v[40]);
    }
    iv[68] = 0;
    iv[69] = 0;
    if (*nd >= *n) {
	goto L40;
    }

/*  ***  SPECIAL CASE HANDLING OF FIRST FUNCTION AND GRADIENT EVALUATION */
/*  ***  -- ASK FOR BOTH RESIDUAL AND JACOBIAN AT ONCE */

    g1 = iv[28];
    y1 = g1 + *p;
    dg7itb_(&b[3], &d__[1], &v[g1], &iv[1], liv, lv, p, p, &v[1], &x[1], &v[
	    y1]);
    if (iv[1] != 1) {
	goto L260;
    }
    v[10] = 0.;
    dv7scp_(p, &v[g1], &c_b14);
    iv[1] = -1;
    qtr1 = iv[77];
    dv7scp_(p, &v[qtr1], &c_b14);
    iv[67] = 0;
    rmat1 = iv[78];
    goto L100;

L40:
    g1 = iv[28];
    y1 = g1 + *p;
    dg7itb_(&b[3], &d__[1], &v[g1], &iv[1], liv, lv, p, p, &v[1], &x[1], &v[
	    y1]);
    if ((i__1 = iv[1] - 2) < 0) {
	goto L50;
    } else if (i__1 == 0) {
	goto L60;
    } else {
	goto L260;
    }

L50:
    v[10] = 0.;
    if (iv[69] == 0) {
	goto L240;
    }
    if (iv[9] != 2) {
	goto L240;
    }
    iv[68] = iv[69];
    dv7cpy_(n, &rd[1], &r__[1]);
    iv[67] = 0;
    goto L240;

L60:
    dv7scp_(p, &v[g1], &c_b14);
    if (iv[35] > 0) {
	goto L230;
    }
    rmat1 = iv[78];
    qtr1 = iv[77];
    rd1 = qtr1 + *p;
    dv7scp_(p, &v[qtr1], &c_b14);
    iv[67] = 0;
    if (*nd < *n) {
	goto L90;
    }
    if (*n1 != 1) {
	goto L90;
    }
    if (iv[35] < 0) {
	goto L100;
    }
    if (iv[69] == iv[7]) {
	goto L70;
    }
    if (iv[68] != iv[7]) {
	goto L90;
    }
    dv7cpy_(n, &r__[1], &rd[1]);
    goto L80;
L70:
    dv7cpy_(n, &rd[1], &r__[1]);
L80:
    dq7apl_(nd, n, p, &dr[dr_offset], &rd[1], &c__0);
    i__1 = min(*n,*p);
    dr7tvm_(nd, &i__1, &v[y1], &v[rd1], &dr[dr_offset], &rd[1]);
    iv[67] = 0;
    goto L110;

L90:
    iv[1] = -2;
    if (iv[35] < 0) {
	iv[1] = -3;
    }
L100:
    dv7scp_(p, &v[y1], &c_b14);
L110:
    dv7scp_(&lh, &v[rmat1], &c_b14);
    goto L240;

/*  ***  COMPUTE F(X)  *** */

L120:
    t = dv2nrm_(&nn, &r__[1]);
    if (t > v[46]) {
	goto L210;
    }
/* Computing 2nd power */
    d__1 = t;
    v[10] += d__1 * d__1 * .5;
    if (*n2 < *n) {
	goto L250;
    }
    if (*n1 == 1) {
	iv[69] = iv[6];
    }
    goto L40;

/*  ***  COMPUTE Y  *** */

L130:
    y1 = iv[28] + *p;
    yi = y1;
    i__1 = *p;
    for (l = 1; l <= i__1; ++l) {
	v[yi] += dd7tpr_(&nn, &dr[l * dr_dim1 + 1], &r__[1]);
	++yi;
/* L140: */
    }
    if (*n2 < *n) {
	goto L250;
    }
    iv[1] = 2;
    if (*n1 > 1) {
	iv[1] = -3;
    }
    goto L240;

/*  ***  COMPUTE GRADIENT INFORMATION  *** */

L150:
    g1 = iv[28];
    ivmode = iv[35];
    if (ivmode < 0) {
	goto L170;
    }
    if (ivmode == 0) {
	goto L180;
    }
    iv[1] = 2;

/*  ***  COMPUTE GRADIENT ONLY (FOR USE IN COVARIANCE COMPUTATION)  *** */

    gi = g1;
    i__1 = *p;
    for (l = 1; l <= i__1; ++l) {
	v[gi] += dd7tpr_(&nn, &r__[1], &dr[l * dr_dim1 + 1]);
	++gi;
/* L160: */
    }
    goto L200;

/*  *** COMPUTE INITIAL FUNCTION VALUE WHEN ND .LT. N *** */

L170:
    if (*n <= *nd) {
	goto L180;
    }
    t = dv2nrm_(&nn, &r__[1]);
    if (t > v[46]) {
	goto L210;
    }
/* Computing 2nd power */
    d__1 = t;
    v[10] += d__1 * d__1 * .5;

/*  ***  UPDATE D IF DESIRED  *** */

L180:
    if (iv[16] > 0) {
	dd7upd_(&d__[1], &dr[dr_offset], &iv[1], liv, lv, n, nd, &nn, n2, p, &
		v[1]);
    }

/*  ***  COMPUTE RMAT AND QTR  *** */

    qtr1 = iv[77];
    rmat1 = iv[78];
    dq7rad_(&nn, nd, p, &v[qtr1], &c_true, &v[rmat1], &dr[dr_offset], &r__[1])
	    ;
    iv[69] = 0;
    if (*n1 > 1) {
	goto L200;
    }
    if (*n2 < *n) {
	goto L250;
    }

/*  ***  SAVE DIAGONAL OF R FOR COMPUTING Y LATER  *** */

    rd1 = qtr1 + *p;
    l = rmat1 - 1;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l += i__;
	v[rd1] = v[l];
	++rd1;
/* L190: */
    }

L200:
    if (*n2 < *n) {
	goto L250;
    }
    if (ivmode > 0) {
	goto L40;
    }
    iv[81] = iv[7];

/*  ***  COMPUTE G FROM RMAT AND QTR  *** */

    dl7vml_(p, &v[g1], &v[rmat1], &v[qtr1]);
    iv[1] = 2;
    if (ivmode == 0) {
	goto L40;
    }
    if (*n <= *nd) {
	goto L40;
    }

/*  ***  FINISH SPECIAL CASE HANDLING OF FIRST FUNCTION AND GRADIENT */

    y1 = g1 + *p;
    iv[1] = 1;
    dg7itb_(&b[3], &d__[1], &v[g1], &iv[1], liv, lv, p, p, &v[1], &x[1], &v[
	    y1]);
    if (iv[1] != 2) {
	goto L260;
    }
    goto L40;

/*  ***  MISC. DETAILS  *** */

/*     ***  X IS OUT OF RANGE (OVERSIZE STEP)  *** */

L210:
    iv[2] = 1;
    goto L40;

/*     ***  BAD N, ND, OR P  *** */

L220:
    iv[1] = 66;
    goto L270;

/*  ***  RECORD EXTRA EVALUATIONS FOR FINITE-DIFFERENCE HESSIAN  *** */

L230:
    ++iv[52];
    ++iv[6];
    iv[7] = iv[6];
    iv[1] = -1;

/*  ***  RETURN FOR MORE FUNCTION OR GRADIENT INFORMATION  *** */

L240:
    *n2 = 0;
L250:
    *n1 = *n2 + 1;
    *n2 += *nd;
    if (*n2 > *n) {
	*n2 = *n;
    }
    goto L999;

/*  ***  PRINT SUMMARY OF FINAL ITERATION AND OTHER REQUESTED ITEMS  *** */

L260:
    g1 = iv[28];
L270:
    ditsum_(&d__[1], &v[g1], &iv[1], liv, lv, p, &v[1], &x[1]);

L999:
    return 0;
/*  ***  LAST CARD OF DRN2GB FOLLOWS  *** */
} /* drn2gb_ */

