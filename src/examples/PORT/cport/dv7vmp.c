/* dv7vmp.f -- translated by f2c (version 20160102).
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

/* Subroutine */ int dv7vmp_(integer *n, doublereal *x, doublereal *y, 
	doublereal *z__, integer *k)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/* ***  SET X(I) = Y(I) * Z(I)**K, 1 .LE. I .LE. N (FOR K = 1 OR -1)  *** */


    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    if (*k >= 0) {
	goto L20;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	x[i__] = y[i__] / z__[i__];
    }
    goto L999;

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	x[i__] = y[i__] * z__[i__];
    }
L999:
    return 0;
/*  ***  LAST CARD OF DV7VMP FOLLOWS  *** */
} /* dv7vmp_ */

