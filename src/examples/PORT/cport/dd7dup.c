/* dd7dup.f -- translated by f2c (version 20160102).
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

/* Subroutine */ int dd7dup_(doublereal *d__, doublereal *hdiag, integer *iv, 
	integer *liv, integer *lv, integer *n, doublereal *v)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal t;
    static integer d0i;
    static doublereal vdfac;
    static integer dtoli;


/*  ***  UPDATE SCALE VECTOR D FOR  DMNH  *** */

/*  ***  PARAMETER DECLARATIONS  *** */


/*  ***  LOCAL VARIABLES  *** */


/*  ***  INTRINSIC FUNCTIONS  *** */
/* /+ */
/* / */
/*  ***  SUBSCRIPTS FOR IV AND V  *** */

/* /6 */
/*     DATA DFAC/41/, DTOL/59/, DTYPE/16/, NITER/31/ */
/* /7 */
/* / */

/* -------------------------------  BODY  -------------------------------- */

    /* Parameter adjustments */
    --iv;
    --v;
    --hdiag;
    --d__;

    /* Function Body */
    i__ = iv[16];
    if (i__ == 1) {
	goto L10;
    }
    if (iv[31] > 0) {
	goto L999;
    }

L10:
    dtoli = iv[59];
    d0i = dtoli + *n;
    vdfac = v[41];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = sqrt((d__1 = hdiag[i__], abs(d__1))), d__3 = vdfac * d__[i__];
	t = max(d__2,d__3);
	if (t < v[dtoli]) {
/* Computing MAX */
	    d__1 = v[dtoli], d__2 = v[d0i];
	    t = max(d__1,d__2);
	}
	d__[i__] = t;
	++dtoli;
	++d0i;
/* L20: */
    }

L999:
    return 0;
/*  ***  LAST CARD OF DD7DUP FOLLOWS  *** */
} /* dd7dup_ */

