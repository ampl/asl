/****************************************************************
Copyright (C) 2014 AMPL Optimization, Inc.; written by David M. Gay.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notice appear in all copies and that
both that the copyright notice and this permission notice and warranty
disclaimer appear in supporting documentation.

The author and AMPL Optimization, Inc. disclaim all warranties with
regard to this software, including all implied warranties of
merchantability and fitness.  In no event shall the author be liable
for any special, indirect or consequential damages or any damages
whatsoever resulting from loss of use, data or profits, whether in an
action of contract, negligence or other tortious action, arising out
of or in connection with the use or performance of this software.
****************************************************************/

/*
 This file is source for bspline.dll, which is discussed in
 "The AMPL Modeling Language -- an Aid to Formulating and Solving
 Optimization Problems'', pp. 95-116 of "Numerical Analysis and
 Optimization", edited by Mehiddin Al-Baali, Lucio Grandinetti, and
 Anton Purnama, Springer Proceedings in Matrhematics and Statistics
 vol. 134, 2015.  Draft: http://ampl.com/REFS/muscat14.pdf .

 For rules to make *.dll files on various systems, see the sample
 makefiles in http://ampl.com/netlib/ampl/solvers/funclink .
*/

#include <string.h>
#include "funcadd.h"

 static real*
squawk(arglist *al, const char *fmt, ...)
{
	AmplExports *ae;
	char buf[1024], *s;
	int n;
	va_list ap;

	ae = al->AE;
	va_start(ap, fmt);
	n = vsnprintf(buf, sizeof(buf), fmt, ap);
	va_end(ap);
	al->Errmsg = s = (char*)TempMem(al->TMI, n + 1);
	if (n < sizeof(buf))
		memcpy(s, buf, n+1);
	else {
		va_start(ap, fmt);
		vsnprintf(s, n+1, fmt, ap);
		va_end(ap);
		}
	return 0;
	}

 static real *
argcheck(arglist *al, const char *who, int *nkp, real **wp)
	/* check arguments, return breakpoint b with b[0] <= x < b[1] */
{
	int i, j, nargs, nd, nd1, ni, ni2, nk;
	real *b, *b0, *be, *ra, *w, x;

	/* nd = spline degree (order of spline >= 0) */
	/* ni >= 1 = number of intervals where x may lie */
	/* nw = ni + nd		# number of basis-function weights */
	/* nk = ni + 2*nd + 1	# number of knots */
	/* nargs = 2 + nw + nk = 3 + 2*ni + 3*nd */

	/* ni = (nargs - 3*nd - 3)/2 */
	/* ni >= 1 ==> nargs >= 5 + 3*nd */

	if ((nargs = al->n) < 5)
		return squawk(al, "%s: too few arguments (%d)\n", who, nargs);
	ra = al->ra;
	nd = ra[0];
	ni2 = nargs - 3*nd - 3;
	if (nd < 0 || ni2 <= 0 || ni2 & 1)
		return squawk(al, "%s(nd, ...): bad (nd = %d, ni2 = %d)\n",
			who, nd, ni2);
	ni = ni2 >> 1;
	nk = ni + 2*nd + 1;
	x = ra[1];
	b0 = ra + 2;	/* start of breakpoints */
	w = b0 + nk;
	be = w - 1;	/* last breakpoint */
	for(b = b0; b < be; ++b) {
		if (b[0] >= b[1])
			return squawk(al, "%s: Breakpoints not strictly increasing: "
					"b[%d] = %.g, b[%d] = %.g.", who,
					(int)(b-b0), b[0], (int)(b-b0) + 1, b[1]);
		}
	for(b = b0; b < be && x > *b; ++b);
	if (b > b0 && x < *b)
		--b;
	if ((i = b - b0) < nd) {
		if (i == nd - 1) {
			++b; /* be slightly tolerant */
			++i;
			}
		else
			return squawk(al, "Too few breakpoints (%d) "
					"before x = %.g", i, x);
		}
	nd1 = nd + 1;
	if ((j = w - b) <= nd1) {
		if (j == nd1) {
			--b; /* be slightly tolerant */
			--i;
			}
		else
			return squawk(al, "Too few breakpoints (%d) "
					"after x = %.g", j, x);
		}
	*nkp = nk;
	*wp = w + i - nd;
	return b;
	}

 static real
bspline0(arglist *al)
{
	/* function value only -- no derivatives */

	AmplExports *ae;
	int i, j, nd, nk;
	real alpha, *b, *f, f0, fj, *ra, *w, x;
	size_t L;

	if (!(b = argcheck(al, "spline", &nk, &w)))
		return 0.;
	if (al->derivs) {
		al->Errmsg = "'bspline0'(...) is not available.";
		return 0.;
		}
	ae = al->AE;
	ra = al->ra;
	nd = ra[0];
	x = ra[1];
	L = (nd + 1)*sizeof(real);
	f = (real*)TempMem(al->TMI, L);
	*f = 1.;
	for(i = 1; i <= nd; ++i) {
		alpha = (b[1]-x)/(b[1] - b[1-i]);
		f0 = *f;
		*f = f0*alpha;
		for(j = 2; j <= i; ++j) {
			fj = (1. - alpha)*f0;
			alpha = (b[j] - x)/(b[j] - b[j-i]);
			f[j-1] = fj + alpha*(f0 = f[j-1]);
			}
		f[i] = (1. - alpha)*f0;
		}
	/* Now f[i], i = 0 ... nd, are B-spline basis-function values, summing to 1. */
	/* See de Boor, p. 131, (5). */
	x = f[0]*w[0];
	for(i = 1; i <= nd; ++i)
		x += f[i]*w[i];
	return x;
	}

 static real
bspline(arglist *al)
{
	/* variant of bspline0() that computes first derivatives */

	AmplExports *ae;
	int i, j, j1, nd, nk;
	real a, af0, alpha, *b, *b0, *d, *d0, dx, *f, f0, fj, oma, *ra, rv, t, *w, x;
	size_t L;

	if (!(b = argcheck(al, "spline", &nk, &w)))
		return 0.;
	if (al->hes) {
		al->Errmsg = "\"bspline\"(...) is not available.";
		return 0;
		}
	ae = al->AE;
	ra = al->ra;
	nd = ra[0];
	x = ra[1];
	L = (((nd+1)*(nd+2))>>1)*sizeof(real);
	f = (real*)TempMem(al->TMI, L);
	*f = 1.;
	for(i = 1; i <= nd; ++i) {
		alpha = (b[1]-x)/(b[1] - b[1-i]);
		f0 = *f;
		f[i] = f0*alpha;
		for(j = 1; j < i;) {
			j1 = j++;
			fj = (1. - alpha)*f0;
			alpha = (b[j] - x)/(b[j] - b[j-i]);
			f[i+j1] = fj + alpha*(f0 = f[j1]);
			}
		f += i;
		f[i] = (1. - alpha)*f0;
		}
	/* Now f[i], i = 0 ... nd, are B-spline basis-function values, summing to 1. */
	/* See de Boor, p. 131, (5). */
	rv = f[0]*w[0];
	for(i = 1; i <= nd; ++i)
		rv += f[i]*w[i];
	if (!(d0 = al->derivs))
		goto ret;

	/* First derivative computation, using backward AD "by hand" */
	/* -- not the recommended procedure in general, but it makes */
	/* this code self-contained.  We visit operations in reverse */
	/* order, overwriting intermediate values by their adjoints  */
	/* (i.e., partials of the final result with respect to the   */
	/* intermediate values). */

	memset(d0, 0, al->n*sizeof(real));
	b0 = ra + 2;
	d = d0 + 2 + (b-b0);
	memcpy(d + nk - nd, f, (nd+1)*sizeof(real));	/* w partials */
	memcpy(f, w, (nd+1)*sizeof(real));		/* f adjoints */
	dx = 0.;
	for(i = nd; i >= 1; --i) {
		alpha = (b[i] - x)/(t = b[i] - b[0]);
		oma = 1. - alpha;
		a = f[i];
		f -= i;
		f0 = f[i-1];
		dx += af0 = a*f0/t;
		d[0] -= af0*alpha;
		d[i] -= af0*oma;
		f[i-1] = a*oma;
		for(j = i; j > 1; j = j1) {
			j1 = j - 1;
			a = f[i+j1];
			f[j1] += a*alpha;
			af0 = a*f0/t;
			dx -= af0;
			d[j-i] += af0*alpha;
			d[j] += af0*oma;
			alpha = (b[j1] - x)/(t = b[j1] - b[j1-i]);
			oma = 1. - alpha;
			f0 = f[j1-1];
			f[j1-1] = a*oma;
			af0 = a*f0/t;
			dx += af0;
			d[j1-i] -= af0*alpha;
			d[j1] -= af0*oma;
			}
		a = f[i];
		f[0] += a*alpha;
		af0 = a*f0/t;
		dx -= af0;
		d[1-i] += af0*alpha;
		d[1] += af0*oma;
		}
	d0[1] = dx;
 ret:
	return rv;
	}

 static real
splinew(arglist *al)
{
	/* Computation based on http://en.wikipedia.org/wiki/De_Boor%27s_algorithm */
	AmplExports *ae;
	int i, j, j1, nd, nd1, ndi, nk;
	real alpha, *b, *f, *ra, *w, x;
	size_t L;

	if (!(b = argcheck(al, "splinew", &nk, &w)))
		return 0.;
	if (al->derivs) {
		al->Errmsg = "splinew'(...) is not available.";
		return 0.;
		}
	ae = al->AE;
	ra = al->ra;
	nd = ra[0];
	x = ra[1];
	L = (nd + 1)*sizeof(real);
	f = (real*)TempMem(al->TMI, L);
	memcpy(f, w, (nd+1)*sizeof(real));
	nd1 = nd + 1;
	for(i = 1; i <= nd; ++i) {
		ndi = nd1 - i;
		for(j1 = 0, j = i-nd; j <= 0; ++j, ++j1) {
			alpha = (x - b[j])/(b[j + ndi] - b[j]);
			f[j1] = (1. - alpha)*f[j1] + alpha*f[j1+1];
			}
		}
	return f[0];
	}

 void
funcadd(AmplExports *ae)
{
	addfunc("bspline", (rfunc)bspline,   0, -5, 0);
	addfunc("bspline0", (rfunc)bspline0, 0, -5, 0);
	addfunc("splinew", (rfunc)splinew,   0, -5, 0);
	}
