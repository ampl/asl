/* dq7rsh.f -- translated by f2c (version 20160102).
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

/* Subroutine */ int dq7rsh_(integer *k, integer *p, logical *havqtr, 
	doublereal *qtr, doublereal *r__, doublereal *w)
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal a, b;
    static integer i__, j;
    static doublereal t, x, y, z__;
    static integer i1, j1, k1;
    static doublereal wj;
    static integer jm1, km1, jp1, pm1;
    extern /* Subroutine */ int dh2rfa_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *);
    extern doublereal dh2rfg_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern /* Subroutine */ int dv7cpy_(integer *, doublereal *, doublereal *)
	    ;


/*  ***  PERMUTE COLUMN K OF R TO COLUMN P, MODIFY QTR ACCORDINGLY  *** */

/*     DIMSNSION R(P*(P+1)/2) */


/*  ***  LOCAL VARIABLES  *** */


    /* Parameter adjustments */
    --w;
    --qtr;
    --r__;

    /* Function Body */

/* +++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++ */

    if (*k >= *p) {
	goto L999;
    }
    km1 = *k - 1;
    k1 = *k * km1 / 2;
    dv7cpy_(k, &w[1], &r__[k1 + 1]);
    wj = w[*k];
    pm1 = *p - 1;
    j1 = k1 + km1;
    i__1 = pm1;
    for (j = *k; j <= i__1; ++j) {
	jm1 = j - 1;
	jp1 = j + 1;
	if (jm1 > 0) {
	    dv7cpy_(&jm1, &r__[k1 + 1], &r__[j1 + 2]);
	}
	j1 += jp1;
	k1 += j;
	a = r__[j1];
	b = r__[j1 + 1];
	if (b != zero) {
	    goto L10;
	}
	r__[k1] = a;
	x = zero;
	z__ = zero;
	goto L40;
L10:
	r__[k1] = dh2rfg_(&a, &b, &x, &y, &z__);
	if (j == pm1) {
	    goto L30;
	}
	i1 = j1;
	i__2 = pm1;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	    i1 += i__;
	    dh2rfa_(&c__1, &r__[i1], &r__[i1 + 1], &x, &y, &z__);
/* L20: */
	}
L30:
	if (*havqtr) {
	    dh2rfa_(&c__1, &qtr[j], &qtr[jp1], &x, &y, &z__);
	}
L40:
	t = x * wj;
	w[j] = wj + t;
	wj = t * z__;
/* L50: */
    }
    w[*p] = wj;
    dv7cpy_(p, &r__[k1 + 1], &w[1]);
L999:
    return 0;
} /* dq7rsh_ */

