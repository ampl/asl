/* ds7ipr.f -- translated by f2c (version 20160102).
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

/* Subroutine */ int ds7ipr_(integer *p, integer *ip, doublereal *h__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal t;
    static integer j1, k1, kk, jm, km, kmj;


/*  APPLY THE PERMUTATION DEFINED BY IP TO THE ROWS AND COLUMNS OF THE */
/*  P X P SYMMETRIC MATRIX WHOSE LOWER TRIANGLE IS STORED COMPACTLY IN H. */
/*  THUS H.OUTPUT(I,J) = H.INPUT(IP(I), IP(J)). */



/* ***  BODY  *** */

    /* Parameter adjustments */
    --ip;
    --h__;

    /* Function Body */
    i__1 = *p;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = ip[i__];
	if (j == i__) {
	    goto L90;
	}
	ip[i__] = abs(j);
	if (j < 0) {
	    goto L90;
	}
	k = i__;
L10:
	j1 = j;
	k1 = k;
	if (j <= k) {
	    goto L20;
	}
	j1 = k;
	k1 = j;
L20:
	kmj = k1 - j1;
	l = j1 - 1;
	jm = j1 * l / 2;
	km = k1 * (k1 - 1) / 2;
	if (l <= 0) {
	    goto L40;
	}
	i__2 = l;
	for (m = 1; m <= i__2; ++m) {
	    ++jm;
	    t = h__[jm];
	    ++km;
	    h__[jm] = h__[km];
	    h__[km] = t;
/* L30: */
	}
L40:
	++km;
	kk = km + kmj;
	++jm;
	t = h__[jm];
	h__[jm] = h__[kk];
	h__[kk] = t;
	j1 = l;
	l = kmj - 1;
	if (l <= 0) {
	    goto L60;
	}
	i__2 = l;
	for (m = 1; m <= i__2; ++m) {
	    jm = jm + j1 + m;
	    t = h__[jm];
	    ++km;
	    h__[jm] = h__[km];
	    h__[km] = t;
/* L50: */
	}
L60:
	if (k1 >= *p) {
	    goto L80;
	}
	l = *p - k1;
	--k1;
	km = kk;
	i__2 = l;
	for (m = 1; m <= i__2; ++m) {
	    km = km + k1 + m;
	    jm = km - kmj;
	    t = h__[jm];
	    h__[jm] = h__[km];
	    h__[km] = t;
/* L70: */
	}
L80:
	k = j;
	j = ip[k];
	ip[k] = -j;
	if (j > i__) {
	    goto L10;
	}
L90:
	;
    }
/* L999: */
    return 0;
/*  ***  LAST LINE OF DS7IPR FOLLOWS  *** */
} /* ds7ipr_ */

