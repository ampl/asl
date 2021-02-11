/* dr7tvm.f -- translated by f2c (version 20160102).
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

/* Subroutine */ int dr7tvm_(integer *n, integer *p, doublereal *y, 
	doublereal *d__, doublereal *u, doublereal *x)
{
    /* System generated locals */
    integer u_dim1, u_offset, i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal t;
    static integer ii, pl, pp1;
    extern doublereal dd7tpr_(integer *, doublereal *, doublereal *);


/*  ***  SET Y TO R*X, WHERE R IS THE UPPER TRIANGULAR MATRIX WHOSE */
/*  ***  DIAGONAL IS IN D AND WHOSE STRICT UPPER TRIANGLE IS IN U. */

/*  ***  X AND Y MAY SHARE STORAGE. */



/*  ***  LOCAL VARIABLES  *** */


/*  ***  BODY  *** */

    /* Parameter adjustments */
    --x;
    u_dim1 = *n;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    --d__;
    --y;

    /* Function Body */
    pl = min(*n,*p);
    pp1 = pl + 1;
    i__1 = pl;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = pp1 - ii;
	t = x[i__] * d__[i__];
	if (i__ > 1) {
	    i__2 = i__ - 1;
	    t += dd7tpr_(&i__2, &u[i__ * u_dim1 + 1], &x[1]);
	}
	y[i__] = t;
/* L10: */
    }
/* L999: */
    return 0;
/*  ***  LAST LINE OF DR7TVM FOLLOWS  *** */
} /* dr7tvm_ */

