/* ds7bqn.f -- translated by f2c (version 20160102).
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
static logical c_false = FALSE_;

/* Subroutine */ int ds7bqn_(doublereal *b, doublereal *d__, doublereal *dst, 
	integer *ipiv, integer *ipiv1, integer *ipiv2, integer *kb, 
	doublereal *l, integer *lv, integer *ns, integer *p, integer *p1, 
	doublereal *step, doublereal *td, doublereal *tg, doublereal *v, 
	doublereal *w, doublereal *x, doublereal *x0)
{
    /* Initialized data */

    static doublereal fudge = 1.0001;
    static doublereal half = .5;
    static doublereal meps2 = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t;
    static integer p0;
    static doublereal t1, dx, ti, xi;
    static integer p1m1;
    static doublereal gts, dst0, dst1, alpha;
    extern doublereal dr7mdc_(integer *);
    extern /* Subroutine */ int dv7shf_(integer *, integer *, doublereal *), 
	    dl7ivm_(integer *, doublereal *, doublereal *, doublereal *);
    extern doublereal dd7tpr_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int i7shft_(integer *, integer *, integer *), 
	    dv7scp_(integer *, doublereal *, doublereal *);
    extern doublereal dv2nrm_(integer *, doublereal *);
    extern /* Subroutine */ int dl7itv_(integer *, doublereal *, doublereal *,
	     doublereal *), dq7rsh_(integer *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *), dv7ipr_(integer *, 
	    integer *, doublereal *), dv7cpy_(integer *, doublereal *, 
	    doublereal *), dv2axy_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal dstmin, dstmax;


/*  ***  COMPUTE BOUNDED MODIFIED NEWTON STEP  *** */

/*     DIMENSION L(P*(P+1)/2) */


/*  ***  LOCAL VARIABLES  *** */


/*  ***  V SUBSCRIPTS  *** */


/* /6 */
/*     DATA DSTNRM/2/, GTSTEP/4/, PHMNFC/20/, PHMXFC/21/, PREDUC/7/, */
/*    1     RADIUS/8/, STPPAR/5/ */
/* /7 */
/* / */

    /* Parameter adjustments */
    --l;
    --v;
    --x0;
    --x;
    --w;
    --tg;
    --td;
    --step;
    --ipiv2;
    --ipiv1;
    --ipiv;
    --dst;
    --d__;
    b -= 3;

    /* Function Body */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    dstmax = fudge * (one + v[21]) * v[8];
    dstmin = (one + v[20]) * v[8];
    dst1 = zero;
    if (meps2 <= zero) {
	meps2 = two * dr7mdc_(&c__3);
    }
    p0 = *p1;
    *ns = 0;
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ipiv1[i__] = i__;
	ipiv2[i__] = i__;
/* L10: */
    }
    i__1 = *p1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	w[i__] = -step[i__] * td[i__];
    }
    alpha = abs(v[5]);
    v[7] = zero;
    gts = -v[4];
    if (*kb < 0) {
	dv7scp_(p, &dst[1], &zero);
    }
    *kb = 1;

/*     ***  -W = D TIMES RESTRICTED NEWTON STEP FROM X + DST/D. */

/*     ***  FIND T SUCH THAT X - T*W IS STILL FEASIBLE. */

L30:
    t = one;
    k = 0;
    i__1 = *p1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = ipiv[i__];
	dx = w[i__] / d__[j];
	xi = x[j] - dx;
	if (xi < b[(j << 1) + 1]) {
	    goto L40;
	}
	if (xi <= b[(j << 1) + 2]) {
	    goto L60;
	}
	ti = (x[j] - b[(j << 1) + 2]) / dx;
	k = i__;
	goto L50;
L40:
	ti = (x[j] - b[(j << 1) + 1]) / dx;
	k = -i__;
L50:
	if (t <= ti) {
	    goto L60;
	}
	t = ti;
L60:
	;
    }

    if (*p > *p1) {
	i__1 = *p - *p1;
	dv7cpy_(&i__1, &step[*p1 + 1], &dst[*p1 + 1]);
    }
    d__1 = -t;
    dv2axy_(p1, &step[1], &d__1, &w[1], &dst[1]);
    dst0 = dst1;
    dst1 = dv2nrm_(p, &step[1]);

/*  ***  CHECK FOR OVERSIZE STEP  *** */

    if (dst1 <= dstmax) {
	goto L80;
    }
    if (*p1 >= p0) {
	goto L70;
    }
    if (dst0 < dstmin) {
	*kb = 0;
    }
    goto L110;

L70:
    k = 0;

/*  ***  UPDATE DST, TG, AND V(PREDUC)  *** */

L80:
    v[2] = dst1;
    dv7cpy_(p1, &dst[1], &step[1]);
    t1 = one - t;
    i__1 = *p1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L90: */
	tg[i__] = t1 * tg[i__];
    }
    if (alpha > zero) {
	d__1 = t * alpha;
	dv2axy_(p1, &tg[1], &d__1, &w[1], &tg[1]);
    }
    v[7] += t * ((one - half * t) * gts + half * alpha * t * dd7tpr_(p1, &w[1]
	    , &w[1]));
    if (k == 0) {
	goto L110;
    }

/*     ***  PERMUTE L, ETC. IF NECESSARY  *** */

    p1m1 = *p1 - 1;
    j = abs(k);
    if (j == *p1) {
	goto L100;
    }
    ++(*ns);
    ipiv2[*p1] = j;
    dq7rsh_(&j, p1, &c_false, &tg[1], &l[1], &w[1]);
    i7shft_(p1, &j, &ipiv[1]);
    i7shft_(p1, &j, &ipiv1[1]);
    dv7shf_(p1, &j, &tg[1]);
    dv7shf_(p1, &j, &dst[1]);
L100:
    if (k < 0) {
	ipiv[*p1] = -ipiv[*p1];
    }
    *p1 = p1m1;
    if (*p1 <= 0) {
	goto L110;
    }
    dl7ivm_(p1, &w[1], &l[1], &tg[1]);
    gts = dd7tpr_(p1, &w[1], &w[1]);
    dl7itv_(p1, &w[1], &l[1], &w[1]);
    goto L30;

/*     ***  UNSCALE STEP  *** */

L110:
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = (i__2 = ipiv[i__], abs(i__2));
	step[j] = dst[i__] / d__[j];
/* L120: */
    }

/*  ***  FUDGE STEP TO ENSURE THAT IT FORCES APPROPRIATE COMPONENTS */
/*  ***  TO THEIR BOUNDS  *** */

    if (*p1 >= p0) {
	goto L150;
    }
    k = *p1 + 1;
    i__1 = p0;
    for (i__ = k; i__ <= i__1; ++i__) {
	j = ipiv[i__];
	t = meps2;
	if (j > 0) {
	    goto L130;
	}
	t = -t;
	j = -j;
	ipiv[i__] = j;
L130:
/* Computing MAX */
	d__3 = (d__1 = x[j], abs(d__1)), d__4 = (d__2 = x0[j], abs(d__2));
	t *= max(d__3,d__4);
	step[j] += t;
/* L140: */
    }

L150:
    dv2axy_(p, &x[1], &one, &step[1], &x0[1]);
    if (*ns > 0) {
	dv7ipr_(&p0, &ipiv1[1], &td[1]);
    }
/* L999: */
    return 0;
/*  ***  LAST LINE OF DS7BQN FOLLOWS  *** */
} /* ds7bqn_ */

