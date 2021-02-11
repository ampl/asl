/* dd7dgb.f -- translated by f2c (version 20160102).
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

static integer c__3 = 3;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c_n1 = -1;
static logical c_false = FALSE_;
static doublereal c_b23 = 1.;

/* Subroutine */ int dd7dgb_(doublereal *b, doublereal *d__, doublereal *dig, 
	doublereal *dst, doublereal *g, integer *ipiv, integer *ka, 
	doublereal *l, integer *lv, integer *p, integer *pc, doublereal *
	nwtst, doublereal *step, doublereal *td, doublereal *tg, doublereal *
	v, doublereal *w, doublereal *x0)
{
    /* Initialized data */

    static doublereal meps2 = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;
    static integer p1;
    static doublereal t1, t2, ti, xi, x0i, rad;
    static integer p1m1;
    static doublereal nred, pred, gnorm;
    extern /* Subroutine */ int dd7dog_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal dr7mdc_(integer *);
    extern /* Subroutine */ int dv7shf_(integer *, integer *, doublereal *), 
	    dl7ivm_(integer *, doublereal *, doublereal *, doublereal *);
    static doublereal gnorm0;
    extern doublereal dd7tpr_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int i7shft_(integer *, integer *, integer *), 
	    dl7vml_(integer *, doublereal *, doublereal *, doublereal *), 
	    dv7scp_(integer *, doublereal *, doublereal *);
    extern doublereal dv2nrm_(integer *, doublereal *);
    extern /* Subroutine */ int dl7itv_(integer *, doublereal *, doublereal *,
	     doublereal *), dq7rsh_(integer *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *), dv7ipr_(integer *, 
	    integer *, doublereal *), dv7cpy_(integer *, doublereal *, 
	    doublereal *), dl7tvm_(integer *, doublereal *, doublereal *, 
	    doublereal *), dv2axy_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), dv7vmp_(integer *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal ghinvg, dnwtst;


/*  ***  COMPUTE DOUBLE-DOGLEG STEP, SUBJECT TO SIMPLE BOUNDS ON X  *** */


/*     DIMENSION L(P*(P+1)/2) */


/*  ***  LOCAL VARIABLES  *** */


/*  ***  V SUBSCRIPTS  *** */


/* /6 */
/*     DATA DGNORM/1/, DST0/3/, DSTNRM/2/, GRDFAC/45/, GTHG/44/, */
/*    1     GTSTEP/4/, NREDUC/6/, NWTFAC/46/, PREDUC/7/, RADIUS/8/, */
/*    2     STPPAR/5/ */
/* /7 */
/* / */
/* /6 */
/*     DATA HALF/0.5D+0/, ONE/1.D+0/, TWO/2.D+0/, ZERO/0.D+0/ */
/* /7 */
/* / */
    /* Parameter adjustments */
    --l;
    --v;
    --x0;
    --w;
    --tg;
    --td;
    --step;
    --nwtst;
    --ipiv;
    --g;
    --dst;
    --dig;
    --d__;
    b -= 3;

    /* Function Body */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    if (meps2 <= 0.) {
	meps2 = 2. * dr7mdc_(&c__3);
    }
    gnorm0 = v[1];
    v[2] = 0.;
    if (*ka < 0) {
	goto L10;
    }
    dnwtst = v[3];
    nred = v[6];
L10:
    pred = 0.;
    v[5] = 0.;
    rad = v[8];
    if (*pc > 0) {
	goto L20;
    }
    dnwtst = 0.;
    dv7scp_(p, &step[1], &c_b5);
    goto L140;

L20:
    p1 = *pc;
    dv7cpy_(p, &td[1], &d__[1]);
    dv7ipr_(p, &ipiv[1], &td[1]);
    dv7scp_(pc, &dst[1], &c_b5);
    dv7cpy_(p, &tg[1], &g[1]);
    dv7ipr_(p, &ipiv[1], &tg[1]);

L30:
    dl7ivm_(&p1, &nwtst[1], &l[1], &tg[1]);
    ghinvg = dd7tpr_(&p1, &nwtst[1], &nwtst[1]);
    v[6] = ghinvg * .5;
    dl7itv_(&p1, &nwtst[1], &l[1], &nwtst[1]);
    dv7vmp_(&p1, &step[1], &nwtst[1], &td[1], &c__1);
    v[3] = dv2nrm_(pc, &step[1]);
    if (*ka >= 0) {
	goto L40;
    }
    *ka = 0;
    dnwtst = v[3];
    nred = v[6];
L40:
    v[8] = rad - v[2];
    if (v[8] <= 0.) {
	goto L100;
    }
    dv7vmp_(&p1, &dig[1], &tg[1], &td[1], &c_n1);
    gnorm = dv2nrm_(&p1, &dig[1]);
    if (gnorm <= 0.) {
	goto L100;
    }
    v[1] = gnorm;
    dv7vmp_(&p1, &dig[1], &dig[1], &td[1], &c_n1);
    dl7tvm_(&p1, &w[1], &l[1], &dig[1]);
    v[44] = dv2nrm_(&p1, &w[1]);
    ++(*ka);
    dd7dog_(&dig[1], lv, &p1, &nwtst[1], &step[1], &v[1]);

/*     ***  FIND T SUCH THAT X - T*STEP IS STILL FEASIBLE. */

    t = 1.;
    k = 0;
    i__1 = p1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = ipiv[i__];
	x0i = x0[j] + dst[i__] / td[i__];
	xi = x0i + step[i__];
	if (xi < b[(j << 1) + 1]) {
	    goto L50;
	}
	if (xi <= b[(j << 1) + 2]) {
	    goto L70;
	}
	ti = (b[(j << 1) + 2] - x0i) / step[i__];
	j = i__;
	goto L60;
L50:
	ti = (b[(j << 1) + 1] - x0i) / step[i__];
	j = -i__;
L60:
	if (t <= ti) {
	    goto L70;
	}
	k = j;
	t = ti;
L70:
	;
    }

/*  ***  UPDATE DST, TG, AND PRED  *** */

    dv7vmp_(&p1, &step[1], &step[1], &td[1], &c__1);
    dv2axy_(&p1, &dst[1], &t, &step[1], &dst[1]);
    v[2] = dv2nrm_(pc, &dst[1]);
    t1 = t * v[45];
    t2 = t * v[46];
/* Computing 2nd power */
    d__1 = v[44] * t1;
    pred = pred - t1 * gnorm * ((t2 + 1.) * gnorm) - t2 * (t2 * .5 + 1.) * 
	    ghinvg - d__1 * d__1 * .5;
    if (k == 0) {
	goto L100;
    }
    dl7vml_(&p1, &w[1], &l[1], &w[1]);
    t2 = 1. - t2;
    i__1 = p1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L80: */
	tg[i__] = t2 * tg[i__] - t1 * w[i__];
    }

/*     ***  PERMUTE L, ETC. IF NECESSARY  *** */

    p1m1 = p1 - 1;
    j = abs(k);
    if (j == p1) {
	goto L90;
    }
    dq7rsh_(&j, &p1, &c_false, &tg[1], &l[1], &w[1]);
    i7shft_(&p1, &j, &ipiv[1]);
    dv7shf_(&p1, &j, &tg[1]);
    dv7shf_(&p1, &j, &td[1]);
    dv7shf_(&p1, &j, &dst[1]);
L90:
    if (k < 0) {
	ipiv[p1] = -ipiv[p1];
    }
    p1 = p1m1;
    if (p1 > 0) {
	goto L30;
    }

/*     ***  UNSCALE STEP, UPDATE X AND DIHDI  *** */

L100:
    dv7scp_(p, &step[1], &c_b5);
    i__1 = *pc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = (i__2 = ipiv[i__], abs(i__2));
	step[j] = dst[i__] / td[i__];
/* L110: */
    }

/*  ***  FUDGE STEP TO ENSURE THAT IT FORCES APPROPRIATE COMPONENTS */
/*  ***  TO THEIR BOUNDS  *** */

    if (p1 >= *pc) {
	goto L140;
    }
    dv2axy_(p, &td[1], &c_b23, &step[1], &x0[1]);
    k = p1 + 1;
    i__1 = *pc;
    for (i__ = k; i__ <= i__1; ++i__) {
	j = ipiv[i__];
	t = meps2;
	if (j > 0) {
	    goto L120;
	}
	t = -t;
	j = -j;
	ipiv[i__] = j;
L120:
/* Computing MAX */
	d__3 = (d__1 = td[j], abs(d__1)), d__4 = (d__2 = x0[j], abs(d__2));
	t *= max(d__3,d__4);
	step[j] += t;
/* L130: */
    }

L140:
    v[1] = gnorm0;
    v[6] = nred;
    v[7] = pred;
    v[8] = rad;
    v[3] = dnwtst;
    v[4] = dd7tpr_(p, &step[1], &g[1]);

/* L999: */
    return 0;
/*  ***  LAST LINE OF DD7DGB FOLLOWS  *** */
} /* dd7dgb_ */

