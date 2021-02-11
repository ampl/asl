/* dh2rfg.f -- translated by f2c (version 20160102).
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

doublereal dh2rfg_(doublereal *a, doublereal *b, doublereal *x, doublereal *y,
	 doublereal *z__)
{
    /* Initialized data */

    static doublereal zero = 0.;

    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__, t, a1, b1;


/*  ***  DETERMINE X, Y, Z SO  I + (1,Z)**T * (X,Y)  IS A 2X2 */
/*  ***  HOUSEHOLDER REFLECTION SENDING (A,B)**T INTO (C,0)**T, */
/*  ***  WHERE  C = -SIGN(A)*SQRT(A**2 + B**2)  IS THE VALUE DH2RFG */
/*  ***  RETURNS. */


/* /+ */
/* / */

/*  ***  BODY  *** */

    if (*b != zero) {
	goto L10;
    }
    *x = zero;
    *y = zero;
    *z__ = zero;
    ret_val = *a;
    goto L999;
L10:
    t = abs(*a) + abs(*b);
    a1 = *a / t;
    b1 = *b / t;
/* Computing 2nd power */
    d__1 = a1;
/* Computing 2nd power */
    d__2 = b1;
    c__ = sqrt(d__1 * d__1 + d__2 * d__2);
    if (a1 > zero) {
	c__ = -c__;
    }
    a1 -= c__;
    *z__ = b1 / a1;
    *x = a1 / c__;
    *y = b1 / c__;
    ret_val = t * c__;
L999:
    return ret_val;
/*  ***  LAST LINE OF DH2RFG FOLLOWS  *** */
} /* dh2rfg_ */

