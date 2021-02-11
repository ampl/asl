/* df7dhb.f -- translated by f2c (version 20160102).
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

static doublereal c_b3 = 0.;

/* Subroutine */ int df7dhb_(doublereal *b, doublereal *d__, doublereal *g, 
	integer *irt, integer *iv, integer *liv, integer *lv, integer *p, 
	doublereal *v, doublereal *x)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, k, l, m;
    static doublereal t, xm;
    static integer mm1;
    static doublereal xm1, del;
    static integer hmi, hes, hpi, hpm;
    static doublereal del0;
    static integer stp0, kind, mm1o2, pp1o2, stpi, stpm, newm1, gsave1;
    extern /* Subroutine */ int dv7scp_(integer *, doublereal *, doublereal *)
	    , dv7cpy_(integer *, doublereal *, doublereal *);
    static logical offsid;


/*  ***  COMPUTE FINITE-DIFFERENCE HESSIAN, STORE IT IN V STARTING */
/*  ***  AT V(IV(FDH)) = V(-IV(H)).  HONOR SIMPLE BOUNDS IN B. */

/*  ***  IF IV(COVREQ) .GE. 0 THEN DF7DHB USES GRADIENT DIFFERENCES, */
/*  ***  OTHERWISE FUNCTION DIFFERENCES.  STORAGE IN V IS AS IN DG7LIT. */

/* IRT VALUES... */
/*     1 = COMPUTE FUNCTION VALUE, I.E., V(F). */
/*     2 = COMPUTE G. */
/*     3 = DONE. */


/*  ***  PARAMETER DECLARATIONS  *** */


/*  ***  LOCAL VARIABLES  *** */


/*  ***  EXTERNAL SUBROUTINES  *** */


/* DV7CPY.... COPY ONE VECTOR TO ANOTHER. */
/* DV7SCP... COPY SCALAR TO ALL COMPONENTS OF A VECTOR. */

/*  ***  SUBSCRIPTS FOR IV AND V  *** */


/* /6 */
/*     DATA HALF/0.5D+0/, HLIM/0.1D+0/, ONE/1.D+0/, TWO/2.D+0/, */
/*    1     ZERO/0.D+0/ */
/* /7 */
/* / */

/* /6 */
/*     DATA COVREQ/15/, DELTA/52/, DELTA0/44/, DLTFDC/42/, F/10/, */
/*    1     FDH/74/, FX/53/, H/56/, KAGQT/33/, MODE/35/, NFGCAL/7/, */
/*    2     SAVEI/63/, SWITCH/12/, TOOBIG/2/, W/65/, XMSAVE/51/ */
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
    *irt = 4;
    kind = iv[15];
    m = iv[35];
    if (m > 0) {
	goto L10;
    }
    hes = abs(iv[56]);
    iv[56] = -hes;
    iv[74] = 0;
    iv[33] = -1;
    v[53] = v[10];
/*        *** SUPPLY ZEROS IN CASE B(1,I) = B(2,I) FOR SOME I *** */
    i__1 = *p * (*p + 1) / 2;
    dv7scp_(&i__1, &v[hes], &c_b3);
L10:
    if (m > *p) {
	goto L999;
    }
    if (kind < 0) {
	goto L120;
    }

/*  ***  COMPUTE FINITE-DIFFERENCE HESSIAN USING BOTH FUNCTION AND */
/*  ***  GRADIENT VALUES. */

    gsave1 = iv[65] + *p;
    if (m > 0) {
	goto L20;
    }
/*        ***  FIRST CALL ON DF7DHB.  SET GSAVE = G, TAKE FIRST STEP  *** */
    dv7cpy_(p, &v[gsave1], &g[1]);
    iv[12] = iv[7];
    goto L80;

L20:
    del = v[52];
    x[m] = v[51];
    if (iv[2] == 0) {
	goto L30;
    }

/*     ***  HANDLE OVERSIZE V(DELTA)  *** */

/* Computing MAX */
    d__2 = 1. / d__[m], d__3 = (d__1 = x[m], abs(d__1));
    del0 = v[44] * max(d__2,d__3);
    del *= .5;
    if ((d__1 = del / del0, abs(d__1)) <= .1) {
	goto L140;
    }

L30:
    hes = -iv[56];

/*  ***  SET  G = (G - GSAVE)/DEL  *** */

    del = 1. / del;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = del * (g[i__] - v[gsave1]);
	++gsave1;
/* L40: */
    }

/*  ***  ADD G AS NEW COL. TO FINITE-DIFF. HESSIAN MATRIX  *** */

    k = hes + m * (m - 1) / 2;
    l = k + m - 2;
    if (m == 1) {
	goto L60;
    }

/*  ***  SET  H(I,M) = 0.5 * (H(I,M) + G(I))  FOR I = 1 TO M-1  *** */

    mm1 = m - 1;
    i__1 = mm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (b[(i__ << 1) + 1] < b[(i__ << 1) + 2]) {
	    v[k] = (v[k] + g[i__]) * .5;
	}
	++k;
/* L50: */
    }

/*  ***  ADD  H(I,M) = G(I)  FOR I = M TO P  *** */

L60:
    ++l;
    i__1 = *p;
    for (i__ = m; i__ <= i__1; ++i__) {
	if (b[(i__ << 1) + 1] < b[(i__ << 1) + 2]) {
	    v[l] = g[i__];
	}
	l += i__;
/* L70: */
    }

L80:
    ++m;
    iv[35] = m;
    if (m > *p) {
	goto L340;
    }
    if (b[(m << 1) + 1] >= b[(m << 1) + 2]) {
	goto L80;
    }

/*  ***  CHOOSE NEXT FINITE-DIFFERENCE STEP, RETURN TO GET G THERE  *** */

/* Computing MAX */
    d__2 = 1. / d__[m], d__3 = (d__1 = x[m], abs(d__1));
    del = v[44] * max(d__2,d__3);
    xm = x[m];
    if (xm < 0.) {
	goto L90;
    }
    xm1 = xm + del;
    if (xm1 <= b[(m << 1) + 2]) {
	goto L110;
    }
    xm1 = xm - del;
    if (xm1 >= b[(m << 1) + 1]) {
	goto L100;
    }
    goto L280;
L90:
    xm1 = xm - del;
    if (xm1 >= b[(m << 1) + 1]) {
	goto L100;
    }
    xm1 = xm + del;
    if (xm1 <= b[(m << 1) + 2]) {
	goto L110;
    }
    goto L280;

L100:
    del = -del;
L110:
    v[51] = xm;
    x[m] = xm1;
    v[52] = del;
    *irt = 2;
    goto L999;

/*  ***  COMPUTE FINITE-DIFFERENCE HESSIAN USING FUNCTION VALUES ONLY. */

L120:
    stp0 = iv[65] + *p - 1;
    mm1 = m - 1;
    mm1o2 = m * mm1 / 2;
    hes = -iv[56];
    if (m > 0) {
	goto L130;
    }
/*        ***  FIRST CALL ON DF7DHB.  *** */
    iv[63] = 0;
    goto L240;

L130:
    if (iv[2] == 0) {
	goto L150;
    }
/*        ***  PUNT IN THE EVENT OF AN OVERSIZE STEP  *** */
L140:
    iv[74] = -2;
    goto L350;
L150:
    i__ = iv[63];
    if (i__ > 0) {
	goto L190;
    }

/*  ***  SAVE F(X + STP(M)*E(M)) IN H(P,M)  *** */

    pp1o2 = *p * (*p - 1) / 2;
    hpm = hes + pp1o2 + mm1;
    v[hpm] = v[10];

/*  ***  START COMPUTING ROW M OF THE FINITE-DIFFERENCE HESSIAN H.  *** */

    newm1 = 1;
    goto L260;
L160:
    hmi = hes + mm1o2;
    if (mm1 == 0) {
	goto L180;
    }
    hpi = hes + pp1o2;
    i__1 = mm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = 0.;
	if (b[(i__ << 1) + 1] < b[(i__ << 1) + 2]) {
	    t = v[53] - (v[10] + v[hpi]);
	}
	v[hmi] = t;
	++hmi;
	++hpi;
/* L170: */
    }
L180:
    v[hmi] = v[10] - v[53] * 2.;
    if (offsid) {
	v[hmi] = v[53] - v[10] * 2.;
    }

/*  ***  COMPUTE FUNCTION VALUES NEEDED TO COMPLETE ROW M OF H.  *** */

    i__ = 0;
    goto L200;

L190:
    x[i__] = v[52];

/*  ***  FINISH COMPUTING H(M,I)  *** */

    stpi = stp0 + i__;
    hmi = hes + mm1o2 + i__ - 1;
    stpm = stp0 + m;
    v[hmi] = (v[hmi] + v[10]) / (v[stpi] * v[stpm]);
L200:
    ++i__;
    if (i__ > m) {
	goto L230;
    }
    if (b[(i__ << 1) + 1] < b[(i__ << 1) + 2]) {
	goto L210;
    }
    goto L200;

L210:
    iv[63] = i__;
    stpi = stp0 + i__;
    v[52] = x[i__];
    x[i__] += v[stpi];
    *irt = 1;
    if (i__ < m) {
	goto L999;
    }
    newm1 = 2;
    goto L260;
L220:
    x[m] = v[51] - del;
    if (offsid) {
	x[m] = v[51] + del * 2.;
    }
    goto L999;

L230:
    iv[63] = 0;
    x[m] = v[51];

L240:
    ++m;
    iv[35] = m;
    if (m > *p) {
	goto L330;
    }
    if (b[(m << 1) + 1] < b[(m << 1) + 2]) {
	goto L250;
    }
    goto L240;

/*  ***  PREPARE TO COMPUTE ROW M OF THE FINITE-DIFFERENCE HESSIAN H. */
/*  ***  COMPUTE M-TH STEP SIZE STP(M), THEN RETURN TO OBTAIN */
/*  ***  F(X + STP(M)*E(M)), WHERE E(M) = M-TH STD. UNIT VECTOR. */

L250:
    v[51] = x[m];
    newm1 = 3;
L260:
    xm = v[51];
/* Computing MAX */
    d__1 = 1. / d__[m], d__2 = abs(xm);
    del = v[42] * max(d__1,d__2);
    xm1 = xm + del;
    offsid = FALSE_;
    if (xm1 <= b[(m << 1) + 2]) {
	goto L270;
    }
    offsid = TRUE_;
    xm1 = xm - del;
    if (xm - del * 2. >= b[(m << 1) + 1]) {
	goto L300;
    }
    goto L280;
L270:
    if (xm - del >= b[(m << 1) + 1]) {
	goto L290;
    }
    offsid = TRUE_;
    if (xm + del * 2. <= b[(m << 1) + 2]) {
	goto L310;
    }

L280:
    iv[74] = -2;
    goto L350;

L290:
    if (xm >= 0.) {
	goto L310;
    }
    xm1 = xm - del;
L300:
    del = -del;
L310:
    switch (newm1) {
	case 1:  goto L160;
	case 2:  goto L220;
	case 3:  goto L320;
    }
L320:
    x[m] = xm1;
    stpm = stp0 + m;
    v[stpm] = del;
    *irt = 1;
    goto L999;

/*  ***  HANDLE SPECIAL CASE OF B(1,P) = B(2,P) -- CLEAR SCRATCH VALUES */
/*  ***  FROM LAST ROW OF FDH... */

L330:
    if (b[(*p << 1) + 1] < b[(*p << 1) + 2]) {
	goto L340;
    }
    i__ = hes + *p * (*p - 1) / 2;
    dv7scp_(p, &v[i__], &c_b3);

/*  ***  RESTORE V(F), ETC.  *** */

L340:
    iv[74] = hes;
L350:
    v[10] = v[53];
    *irt = 3;
    if (kind < 0) {
	goto L999;
    }
    iv[7] = iv[12];
    gsave1 = iv[65] + *p;
    dv7cpy_(p, &g[1], &v[gsave1]);
    goto L999;

L999:
    return 0;
/*  ***  LAST LINE OF DF7DHB FOLLOWS  *** */
} /* df7dhb_ */

