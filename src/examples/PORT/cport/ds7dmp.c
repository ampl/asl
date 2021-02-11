/* ds7dmp.f -- translated by f2c (version 20160102).
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

/* Subroutine */ int ds7dmp_(integer *n, doublereal *x, doublereal *y, 
	doublereal *z__, integer *k)
{
    /* Initialized data */

    static doublereal one = 1.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l;
    static doublereal t;


/* ***  SET X = DIAG(Z)**K * Y * DIAG(Z)**K */
/* ***  FOR X, Y = COMPACTLY STORED LOWER TRIANG. MATRICES */
/* ***  K = 1 OR -1. */

/* /6S */
/*     DOUBLE PRECISION X(1), Y(1), Z(N) */
/* /7S */
/* / */
    /* Parameter adjustments */
    --z__;
    --x;
    --y;

    /* Function Body */

    l = 1;
    if (*k >= 0) {
	goto L30;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = one / z__[i__];
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    x[l] = t * y[l] / z__[j];
	    ++l;
/* L10: */
	}
/* L20: */
    }
    goto L999;

L30:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t = z__[i__];
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    x[l] = t * y[l] * z__[j];
	    ++l;
/* L40: */
	}
/* L50: */
    }
L999:
    return 0;
/*  ***  LAST CARD OF DS7DMP FOLLOWS  *** */
} /* ds7dmp_ */

