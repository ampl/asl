/* i7pnvr.f -- translated by f2c (version 20160102).
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

/* Subroutine */ int i7pnvr_(integer *n, integer *x, integer *y)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;


/*  ***  SET PERMUTATION VECTOR X TO INVERSE OF Y  *** */


    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = y[i__];
	x[j] = i__;
/* L10: */
    }

/* L999: */
    return 0;
/*  ***  LAST LINE OF I7PNVR FOLLOWS  *** */
} /* i7pnvr_ */

