/*******************************************************************
Copyright (C) 2017, 2018, 2019, 2020 AMPL Optimization, Inc.; written by David M. Gay.

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
*******************************************************************/

#define PSHVREAD
#include "asl.h"
#include "psinfo.h"
#include "opno2.h"
#include "errchk.h"

#ifdef __cplusplus
 extern "C" {
#endif

 extern void conpgrd_ew_ASL(EvalWorkspace*, int nc, real *X, real *G, fint *nerror);
 extern void conpgrd_nomap_ew_ASL(EvalWorkspace*, int nc, real *X, real *G, fint *nerror);
 extern real conpival_ew_ASL(EvalWorkspace*, int nc, real *X, fint *ne);
 extern real conpival_nomap_ew_ASL(EvalWorkspace*, int nc, real *X, fint *ne);
 extern void conpval_ew_ASL(EvalWorkspace*, real *X, real *F, fint *nerror);
 extern real eval2_ASL(int *o, EvalWorkspace*);
 extern void jacpval_ew_ASL(EvalWorkspace*, real *X, real *JAC, fint *nerror);
 extern int  lconpval_ew_ASL(EvalWorkspace*, int nc, real *X, fint *ne);
 extern void objpgrd_ew_ASL(EvalWorkspace*, int nobj, real *X, real *G, fint *nerror);
 extern real objpval_ew_ASL(EvalWorkspace*, int nobj, real *X, fint *nerror);

#ifdef __cplusplus
	}
#endif

#ifdef No_dtoa

 static real
Round(real x, int prec)
{
	real scale;
	int flip;

	if (!x)
		return x;
	flip = 0;
	if (x < 0.) {
		x = -x;
		flip = 1;
		}
	if (!prec)
		x = floor(x + 0.5);
	else if (prec > 0) {
		scale = mypow(10., (real)prec);
		x = floor(x*scale + 0.5) / scale;
		}
	else {
		scale = mypow(10., -(real)prec);
		x = scale*floor(x/scale + 0.5);
		}
	return flip ? -x : x;
	}
#else

 static real
Round(real x, int prec)
{
	char *b, *s, sbuf[400], *se;
	int decpt, L, sign;
	char buf[96];

	if (!x)
		return x;
	s = dtoa_r(x, 3, prec, &decpt, &sign, &se, sbuf, sizeof(sbuf));
	if (decpt == 9999) {
 zreturn:
		return x;
		}
	L = se - s;
	if (L <= 0) {
		x = 0.;
		goto zreturn;
		}
	if (L > 80)
		se = s + 80;
	b = buf;
	if (sign)
	*b++ = '-';
	*b++ = '.';
	while(s < se)
		*b++ = *s++;
	*b = 0;
	if (decpt)
		snprintf(b, buf + sizeof(buf) - b, "e%d", decpt);
	return strtod(buf, (char **)0);
	}
#endif

 typedef struct
jb_st {
	jmp_buf jb;
	} jb_st;

 static int
rcompj(const void *a, const void *b, void *v)
{
	jb_st *J;
	real t = *(real *)a - *(real *)b;

	if (!t) {
		J = (jb_st*)v;
		longjmp(J->jb, 1);
		}
	return t < 0 ? -1 : 1;
	}

#define introuble(who, a, jv) introuble_ASL(ew, who, a, jv)
#define introuble2(who, a, b, jv) introuble2_ASL(ew, who, a, b, jv)
#define zero_div(L,s) zero_div_ASL(ew, L, s)

 real
eval2_ASL(int *o, EvalWorkspace *ew)
{
	Condptrs *cp;
	Eresult *r;
	GOps *g;
	Minmaxptrs *mmp;
	TMInfo T, *T1, *T1prev;
	U rv, tv;
	arglist *al;
	char buf[32], *s;
	const char **sa;
	derpblock **db, **db1;
	func_info *fi;
	int i, j, k, n, nr, ns, *o1, **pop, prec, sign, wd, wdf, z;
	jb_st J;
	plterm *p, **pp;
	real L, R, *bs, rbuf[128], *rp, t, t0, t1, *w;
	tfinfo **ptfi, *tfi;
	void **v;
	static real Le10;

	rv.d = 0.;
	if (!o)
		goto done;
	w = ew->w;
	wd = ew->wantderiv;

 top:
	switch(*o) {
	  case OPRET:
		rv.d = w[o[1]];
		goto done;
	  case OPPLUS0:
		w[o[1]] = w[o[2]] + w[o[3]];
		o += 4;
		goto top;
	  case OPPLUS01:
	  case OPPLUS10:
	  case OPPLUS2:
		w[o[2]] = w[o[3]] + w[o[4]];
		o += 5;
		goto top;
	  case OPMINUS0:
		w[o[1]] = w[o[2]] - w[o[3]];
		o += 4;
		goto top;
	  case OPMINUS01:
	  case OPMINUS10:
	  case OPMINUS2:
		w[o[2]] = w[o[3]] - w[o[4]];
		o += 5;
		goto top;
	  case OPMULT0:
		w[o[1]] = w[o[2]] * w[o[3]];
		o += 4;
		goto top;
	  case OPMULT01:
	  case OPMULT10:
	  case OPMULT2:
		w[o[2]] = w[o[3]] * w[o[4]];
		o += 5;
		goto top;
	  case OPDIV0:
		L = w[o[2]];
		if (!(R = w[o[3]])
#ifdef WANT_INFNAN
		 && !L
#endif
			)
			zero_div(L, "/");
		w[o[1]] = L / R;
		o += 4;
		goto top;
	  case OPDIV10:
		L = w[o[3]];
		if (!(R = w[o[4]])
#ifdef WANT_INFNAN
		 && !L
#endif
			)
			zero_div(L, "/");
		r = (Eresult*)(w + o[2]);
		r->O = L / R;
		if (wd)
			r->dL = 1. / R;
		o += 5;
		goto top;
	  case OPDIV01:
		L = w[o[3]];
		if (!(R = w[o[4]])
#ifdef WANT_INFNAN
		 && !L
#endif
			)
			zero_div(L, "/");
		r = (Eresult*)(w + o[2]);
		r->O = L / R;
		if (wd) {
			rv.d = 1. / R;
			r->dL = -r->O * rv.d; /* dR */
			r->dL2 = -2. * rv.d * r->dL; /* dR2 */
			}
		o += 5;
		goto top;
	  case OPDIV2:
		L = w[o[3]];
		if (!(R = w[o[4]])
#ifdef WANT_INFNAN
		 && !L
#endif
			)
			zero_div(L, "/");
		r = (Eresult*)(w + o[2]);
		r->O = L / R;
		if (wd) {
			r->dL = t = 1. / R;
			r->dLR = -t*t;
			r->dR = -r->O * t;
			r->dL2 = -2. * t * r->dR; /* dR2 */
			}
		o += 5;
		goto top;
	  case n_OPREM0:
		L = w[o[2]];
		R = w[o[3]];
		rv.d = fmod(L,R);
		if (errchk(rv))
			introuble2("fmod",L,R,1);
		w[o[1]] = rv.d;
		o += 4;
		goto top;
	  case nOPREM10:
		L = w[o[3]];
		R = w[o[4]];
		rv.d = fmod(L,R);
		if (errchk(rv))
			introuble2("fmod",L,R,1);
		r = (Eresult*)(w + o[2]);
		r->O = rv.d;
		o += 5;
		goto top;
	  case nOPREM01:
	  case nOPREM2:
		L = w[o[3]];
		R = w[o[4]];
		rv.d = fmod(L,R);
		if (errchk(rv))
			introuble2("fmod",L,R,1);
		r = (Eresult*)(w + o[2]);
		r->O = rv.d;
		if (wd) {
			t = -L / R;
			r->dL = t >= 0. ? floor(t) : ceil(t);	/* dR */
			}
		o += 5;
		goto top;
	  case n_OPPOW0:
		rv.d = w[i = o[1]] = mypow(L = w[o[2]], R = w[o[3]]);
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				return w[i] = Infinity;
				}
#endif
			introuble2("pow",L,R,1);
			}
		o += 4;
		goto top;
	  case OP1POW_g:
		g = (GOps*)(w + o[1]);
		rv.d = mypow(L = w[o[2]], R = w[o[3]]);
		o += 4;
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				g->dL = negInfinity;
				g->O = g->dL2 = Infinity;
				goto top;
				}
#endif
			introuble2("pow",L,R,1);
			}
		g->O = rv.d;
		if (wd) {
			if (L) {
				g->dL = R * (rv.d/L);
				g->dL2 = (R - 1.) * g->dL / L;
				}
			else if (R > 1.) {
				g->dL = 0.;
				if (R >= 2.) {
					g->dL2 = R > 2. ? 0. : 2.;
					goto top;
					}
#ifdef WANT_INFNAN
				e->dL2 = Infinity;
#else
				introuble2("pow\"",L,R,3);
#endif
				}
			else if (R == 1.) {
				g->dL = 1.;
				g->dL2 = 0.;
				}
			else if (R == 0.)
				g->dL = g->dL2 = 0.;
			else { /* 0 < R < 1 */
#ifdef WANT_INFNAN
				g->dL = Infinity;
				g->dL2 = negInfinity;
#else
				goto badpowder;
#endif
				}
			}
		goto top;
	  case nOPPOW10:
		r = (Eresult*)(w + o[2]);
		rv.d = mypow(L = w[o[3]], R = w[o[4]]);
		o += 5;
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				r->dL = negInfinity;
				r->O = r->dL2 = Infinity;
				goto top;
				}
#endif
			introuble2("pow",L,R,1);
			}
		r->O = rv.d;
		if (wd) {
			if (L) {
				r->dL = R * (rv.d/L);
				r->dL2 = (R - 1.) * r->dL / L;
				}
			else if (R > 1.) {
				r->dL = 0.;
				if (R >= 2.) {
					r->dL2 = R > 2. ? 0. : 2.;
					goto top;
					}
#ifdef WANT_INFNAN
				e->dL2 = Infinity;
#else
				introuble2("pow\"",L,R,3);
#endif
				}
			else if (R == 1.) {
				r->dL = 1.;
				r->dL2 = 0.;
				}
			else if (R == 0.)
				r->dL = r->dL2 = 0.;
			else { /* 0 < R < 1 */
#ifdef WANT_INFNAN
				r->dL = Infinity;
				r->dL2 = negInfinity;
#else
				goto badpowder;
#endif
				}
			}
		goto top;
	  case OPCPOW_g:
		g = (GOps*)(w + o[1]);
		g->O = rv.d = mypow(L = w[o[2]], R = w[o[3]]);
		o += 4;
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				g->dL = negInfinity;
				g->O = g->dL2 = Infinity;
				goto top;
				}
#endif
 badpow:
			introuble2("pow",L,R,1);
			break; /* not reached */
			}
		if (wd) {
			if (L > 0.) {
				g->dL = (t = log(L)) * rv.d;
				g->dL2 = t * g->dL;
				}
			else if (L != 0.)
				goto badpowder;
			else {
				g->dL = g->dL2 = 0.;
#ifndef WANT_INFNAM
				if (R < 2.)
					introuble2("pow\"",L,R,3);
#endif
				}
			}
		goto top;
	  case nOPPOW01:
		r = (Eresult*)(w + o[2]);
		r->O = rv.d = mypow(L = w[o[3]], R = w[o[4]]);
		o += 5;
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				r->dL = negInfinity;
				r->O = r->dL2 = Infinity;
				goto top;
				}
#endif
			goto badpow;
			}
		if (wd) {
			if (L > 0.) {
				r->dL = (t = log(L)) * rv.d;
				r->dL2 = t * r->dL;
				}
			else if (L != 0.)
				goto badpowder;
			else {
				r->dL = r->dL2 = 0.;
#ifndef WANT_INFNAM
				if (R < 2.)
					introuble2("pow\"",L,R,3);
#endif
				}
			}
		goto top;
	  case n_OPPOW2:
		r = (Eresult*)(w + o[2]);
		r->O = rv.d = mypow(L = w[o[3]], R = w[o[4]]);
		o += 5;
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				if (wd) {
					r->dL = r->dR = negInfinity;
					r->dL2 = r->dLR = r->dR2 = Infinity;
					}
				r->O = Infinity;
				goto top;
				}
#endif
			introuble2("pow",L,R,1);
			}
		if (wd) {
			if (L > 0.) {
				t1 = rv.d / L;
				t = log(L);
				r->dL = R*t1;
				r->dR = t * rv.d;
				r->dL2 = (R - 1.) * (r->dL / L);
				r->dLR = t1 * (1. + R*t);
				r->dR2 = t * r->dR;
				goto top;
				}
			r->dL = r->dR = r->dL2 = r->dLR = r->dR2 = 0.;
			if (L != 0.) {
 badpowder:
				introuble2("pow'",L,R,2);
				break; /* not reached */
				}
			if (R > 1.) {
				if (R == 2.)
					r->dL2 = 2;
				else if (R < 2.) {
					r->dL2 = Infinity;
#ifndef WANT_INFNAN
					introuble2("pow\"",L,R,3);
#endif
					}
				goto top;
				}
			if (R == 1.)
				r->dL = 1.;
			else { /* 0 < R < 1 */
#ifdef WANT_INFNAN
				r->dL = Infinity;
				r->dL2 = r->dLR = negInfinity;
#else
				goto badpowder;
#endif
				}
			}
		goto top;
	  case nOPPOW1i:
		r = (Eresult*)(w + o[2]);
		r->O = rv.d = mypow(L = w[o[3]], R = w[o[4]]);
		o += 5;
		/* R = constant integer */
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				if (wd)
					r->dL = r->dR = r->dL2 = r->dR2 = r->dLR = Infinity;
				r->O = Infinity;
				goto top;
				}
#endif
			introuble2("pow",L,R,1);
			}
		if (wd) {
			if (L != 0.) {
				r->dL = R * (rv.d/L);
				r->dL2 = (R - 1.) * r->dL / L;
				}
			else {
				if (R > 1.)
					r->dL = r->dL2 = 0.;
				else /* R < 0 */
#ifdef WANT_INFNAN
					r->dL = r->dL2 = Infinity;
#else
					goto badpowder;
#endif
				}
			}
		goto top;
	  case n_OPLESS0:
		rv.d = w[o[2]] - w[o[3]];
		if (rv.d < 0.)
			rv.d = 0.;
		w[o[1]] = rv.d;
		o += 4;
		goto top;
#ifdef X64_bit_pointers
	  case OPLESS10align:
		db = (derpblock**)(o + 6);
		goto more_OPLESS10;
#endif
	  case nOPLESS10:
		db = (derpblock**)(o + 5);
alignarg(more_OPLESS10:)
		r = (Eresult*)(w + o[2]);
		db1 = (derpblock**)&r->dL2;
		rv.d = w[o[3]] - w[o[4]];
		if (rv.d < 0.) {
			rv.d = r->dL = 0.;
			*db1 = db[1];
			}
		else {
			r->dL = 1.;
			*db1 = db[0];
			}
		r->O = rv.d;
		o = (int*)&db[2];
		goto top;
#ifdef X64_bit_pointers
	  case OPLESS01align:
		db = (derpblock**)(o + 6);
		goto more_OPLESS01;
#endif
	  case nOPLESS01:
		db = (derpblock**)(o + 5);
alignarg(more_OPLESS01:)
		r = (Eresult*)(w + o[2]);
		db1 = (derpblock**)&r->dL2;
		rv.d = w[o[3]] - w[o[4]];
		if (rv.d < 0.) {
			rv.d = r->dL = 0.;
			*db1 = db[1];
			}
		else {
			r->dL = -1.;
			*db1 = db[0];
			}
		r->O = rv.d;
		o = (int*)&db[2];
		goto top;
#ifdef X64_bit_pointers
	  case OPLESS2align:
		db = (derpblock**)(o + 6);
		goto more_OPLESS2;
#endif
	  case nOPLESS2:
		db = (derpblock**)(o + 5);
alignarg(more_OPLESS2:)
		r = (Eresult*)(w + o[2]);
		db1 = (derpblock**)&r->dR;
		rv.d = w[o[3]] - w[o[4]];
		if (rv.d < 0.) {
			rv.d = r->dL = r->dL2 = 0.;
			*db1 = db[1];
			}
		else {
			r->dL = 1.;
			r->dL2 = -1.;
			*db1 = db[0];
			}
		r->O = rv.d;
		o = (int*)&db[2];
		goto top;
	  case OPMINLIST0:
		i = 3;
		rv.d = w[o[3]];
		k = i + o[2];
		while(++i < k) {
			t = w[o[i]];
			if (rv.d > t)
				rv.d = t;
			}
		w[o[1]] = rv.d;
		o += k;
		goto top;
#ifdef X64_bit_pointers
	  case OPMINLISTalign:
		o1 = o + 5;
		goto have_o1;
#endif
	  case OPMINLIST1:
		o1 = o + 4;
alignarg(have_o1:)
		n = o[3];
		rv.d = w[*o1];
		i = j = 0;
		while(++i < n) {
			t = w[o1[i]];
			if (rv.d > t) {
				rv.d = t;
				j = i;
				}
			}
 finish_xlist:
		mmp = (Minmaxptrs*)&o1[n];
		w[i = o[2]] = rv.d;
		o = (int*)&mmp[n] + 2;
		*(Minmaxptrs**)&w[i+4] = mmp += j;
		*(derpblock**)&w[i+5] = mmp->db;
		goto top;
	  case OPMAXLIST0:
		i = 3;
		rv.d = w[o[3]];
		k = i + o[2];
		while(++i < k) {
			t = w[o[i]];
			if (rv.d < t)
				rv.d = t;
			}
		w[o[1]] = rv.d;
		o += k;
		goto top;
#ifdef X64_bit_pointers
	  case OPMAXLISTalign:
		o1 = o + 5;
		goto have_o1x;
#endif
	  case OPMAXLIST1:
		o1 = o + 4;
alignarg(have_o1x:)
		n = o[3];
		rv.d = w[*o1];
		i = j = 0;
		while(++i < n) {
			t = w[o1[i]];
			if (rv.d < t) {
				rv.d = t;
				j = i;
				}
			}
		goto finish_xlist;
	  case n_FLOOR:
		w[o[1]] = floor(w[o[2]]);
		o += 3;
		goto top;
	  case n_CEIL:
		w[o[1]] = ceil(w[o[2]]);
		o += 3;
		goto top;
	  case n_OPABS0:
		if ((rv.d = w[o[2]]) < 0.)
			rv.d = -rv.d;
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPABS_g:
		t = 1.;
		if ((rv.d = w[o[2]]) < 0.) {
			rv.d = -rv.d;
			t = -1.;
			}
		g = (GOps*)(w + o[1]);
		g->O = rv.d;
		g->dL = t;
		o += 3;
		goto top;
	  case n_OPABS1:
		t = 1.;
		if ((rv.d = w[o[3]]) < 0.) {
			rv.d = -rv.d;
			t = -1.;
			}
		r = (Eresult*)(w + o[2]);
		r->O = rv.d;
		r->dL = t;
		o += 4;
		goto top;
	  case OPUMINUS0:
		w[o[1]] = -w[o[2]];
		o += 3;
		goto top;
	  case OPUMINUS1:
		w[o[2]] = -w[o[3]];
		o += 4;
		goto top;
	  Opalign(OP_ORalign)
	  case n_OPOR:
		pop = (int**)&o[3];
		if (w[o[2]] != 0.) {
			w[o[1]] = 1.;
			o = *pop;
			goto top;
			}
		o = (int*)&pop[1];
		goto top;
	  Opalign(OP_ANDalign)
	  case n_OPAND:
		pop = (int**)&o[3];
		if (w[o[2]] == 0.) {
			w[o[1]] = 0.;
			o = *pop;
			goto top;
			}
		o = (int*)&pop[1];
		goto top;
	  case n_OPLT:
	  case n_OPNOTATMOST:
		w[o[1]] = w[o[2]] < w[o[3]] ? 1. : 0.;
		o += 4;
		goto top;
	  case n_OPLE:
	  case n_OPATLEAST:
		w[o[1]] = w[o[2]] <= w[o[3]] ? 1. : 0.;
		o += 4;
		goto top;
	  case OPEQ:
	  case n_OPEXACTLY:
		w[o[1]] = w[o[2]] == w[o[3]] ? 1. : 0.;
		o += 4;
		goto top;
	  case n_OPGE:
	  case n_OPATMOST:
		w[o[1]] = w[o[2]] >= w[o[3]] ? 1. : 0.;
		o += 4;
		goto top;
	  case n_OPGT:
	  case n_OPNOTATLEAST:
		w[o[1]] = w[o[2]] > w[o[3]] ? 1. : 0.;
		o += 4;
		goto top;
	  case n_OPNE:
	  case n_OPNOTEXACTLY:
		w[o[1]] = w[o[2]] != w[o[3]] ? 1. : 0.;
		o += 4;
		goto top;
	  case n_OPNOT:
		w[o[1]] = w[o[2]] == 0. ? 1. : 0.;
		o += 3;
		goto top;
#ifdef X64_bit_pointers
	  case OP_IMPELSE_align:
	  case OPIF0align:
		++o;
#endif
	  case n_OPIMPELSE:
	  case nOPIF0:	/* and OPIFSYM */
		pop = (int**)&o[2];
		if (w[o[1]] == 0.)
			++pop;
		o = *pop;
		goto top;
#ifdef X64_bit_pointers
	  case OPIF1align:
	  case OPIF11align:
	  case OPIF12align:
	  case OPIF13align:
		cp = (Condptrs*)&o[6];
		goto more_if;
#endif
	  case nOPIF1:
	  case nOPIF11:
	  case nOPIF12:
	  case nOPIF13:
		cp = (Condptrs*)&o[5];
alignarg(more_if:)
		v = (void**)&w[o[4]];
/* //		v[2] = *(int**)&cp[2];	*/
		if (w[o[3]] == 0.)
			++cp;
		v[0] = cp->db;
		v[1] = cp;
		o = cp->e;
		goto top;
	  case OPtanh0:
		if ((t = w[o[2]]) >= 175.)
			rv.d = 1.;
		else if (t <= -175.)
			rv.d = -1.;
		else {
			rv.d = tanh(t);
			if (errchk(rv))
				introuble("tanh",t,1);
			}
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPtanh_g:
		g = (GOps*)(w + o[1]);
		if ((t = w[o[2]]) >= 175.) {
			g->O = 1.;
			if (wd)
				g->dL2 = g->dL = 0.;
			}
		else if (t <= -175.) {
			g->O = -1.;
			if (wd)
				g->dL2 = g->dL = 0.;
			}
		else {
			rv.d = tanh(t);
			if (errchk(rv))
				introuble("tanh",t,1);
			if (wd) {
				tv.d = cosh(t);
				if (errchk(tv))
					introuble("tanh'", t, 2);
				else {
					t1 = 1. / tv.d;
					g->dL2 = -(rv.d + rv.d) * (g->dL = t1*t1);
					}
				}
			g->O = rv.d;
			}
		o += 3;
		goto top;
	  case OPtanh1:
		r = (Eresult*)(w + o[2]);
		if ((t = w[o[3]]) >= 175.) {
			r->O = 1.;
			if (wd)
				r->dL = r->dL2 = 0.;
			}
		else if (t <= -175.) {
			r->O = -1.;
			if (wd)
				r->dL = r->dL2 = 0.;
			}
		else {
			rv.d = tanh(t);
			if (errchk(rv))
				introuble("tanh",t,1);
			if (wd) {
				tv.d = cosh(t);
				if (errchk(tv))
					introuble("tanh'", t, 2);
				else {
					t1 = 1. / tv.d;
					r->dL2 = -(rv.d + rv.d) * (r->dL = t1*t1);
					}
				}
			r->O = rv.d;
			}
		o += 4;
		goto top;
	  case OP_tan0:
		rv.d = tan(w[o[2]]);
		if (errchk(rv))
			introuble("tan",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPtan_g:
		rv.d = tan(t = w[o[2]]);
		if (errchk(rv))
			introuble("tan",t,1);
		g = (GOps*)(w + o[1]);
		if (wd) {
			tv.d = cos(t);
			if (errchk(tv) || !tv.d)
				introuble("tan'", t, 2);
			else {
				t1 = 1. / tv.d;
				g->dL2 = (rv.d + rv.d) * (g->dL = t1*t1);
				}
			}
		g->O = rv.d;
		o += 3;
		goto top;
	  case OP_tan1:
		rv.d = tan(t = w[o[3]]);
		if (errchk(rv))
			introuble("tan",t,1);
		r = (Eresult*)(w + o[2]);
		if (wd) {
			tv.d = cos(t);
			if (errchk(tv) || !tv.d)
				introuble("tan'", t, 2);
			else {
				t1 = 1. / tv.d;
				r->dL2 = (rv.d + rv.d) * (r->dL = t1*t1);
				}
			}
		r->O = rv.d;
		o += 4;
		goto top;
	  case OP_sqrt0:
		t = w[o[2]];
		if (t < 0.) {
 badsqrt:
			introuble("sqrt",t,1);
			}
		rv.d = sqrt(t);
		if (errchk(rv))
			goto badsqrt;
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPsqrt_g:
		t = w[o[2]];
		if (t < 0.)
			goto badsqrt;
		g = (GOps*)(w + o[1]);
		rv.d = sqrt(t);
		if (errchk(rv))
			goto badsqrt;
		if (wd) {
			if (rv.d <= 0.)
				introuble("sqrt'",t,2);
			else {
				g->dL = 0.5 / rv.d;
				g->dL2 = -0.5 * g->dL / t;
				}
			}
		g->O = rv.d;
		o += 3;
		goto top;
	  case OP_sqrt1:
		t = w[o[3]];
		if (t < 0.)
			goto badsqrt;
		r = (Eresult*)(w + o[2]);
		rv.d = sqrt(t);
		if (errchk(rv))
			goto badsqrt;
		if (wd) {
			if (rv.d <= 0.)
				introuble("sqrt'",t,2);
			else {
				r->dL = 0.5 / rv.d;
				r->dL2 = -0.5 * r->dL / t;
				}
			}
		r->O = rv.d;
		o += 4;
		goto top;
	  case OP_sinh0:
		rv.d = sinh(w[o[2]]);
		if (errchk(rv))
			introuble("sinh",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPsinh_g:
		rv.d = sinh(t = w[o[2]]);
		if (errchk(rv))
			introuble("sinh",t,1);
		g = (GOps*)(w + o[1]);
		if (wd) {
			tv.d = cosh(t);
			if (errchk(tv))
				introuble("sinh'", t, 2);
			g->dL = tv.d;
			g->dL2 = rv.d;
			}
		g->O = rv.d;
		o += 3;
		goto top;
	  case OP_sinh1:
		rv.d = sinh(t = w[o[3]]);
		if (errchk(rv))
			introuble("sinh",t,1);
		r = (Eresult*)(w + o[2]);
		if (wd) {
			tv.d = cosh(t);
			if (errchk(tv))
				introuble("sinh'", t, 2);
			r->dL = tv.d;
			r->dL2 = rv.d;
			}
		r->O = rv.d;
		o += 4;
		goto top;
	  case OP_sin0:
		rv.d = sin(w[o[2]]);
		if (errchk(rv))
			introuble("sin",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPsin_g:
		rv.d = sin(t = w[o[2]]);
		if (errchk(rv))
			introuble("sin",t,1);
		g = (GOps*)(w + o[1]);
		if (wd) {
			tv.d = cos(t);
			if (errchk(tv))
				introuble("sin'", t, 2);
			g->dL = tv.d;
			g->dL2 = -rv.d;
			}
		g->O = rv.d;
		o += 3;
		goto top;
	  case OP_sin1:
		rv.d = sin(t = w[o[3]]);
		if (errchk(rv))
			introuble("sin",t,1);
		r = (Eresult*)(w + o[2]);
		if (wd) {
			tv.d = cos(t);
			if (errchk(tv))
				introuble("sin'", t, 2);
			r->dL = tv.d;
			r->dL2 = -rv.d;
			}
		r->O = rv.d;
		o += 4;
		goto top;
	  case OP_log100:
		rv.d = log10(w[o[2]]);
		if (errchk(rv))
			introuble("log10",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPlog10_g:
		rv.d = log10(t = w[o[2]]);
		if (errchk(rv))
			introuble("log10",t,1);
		g = (GOps*)(w + o[1]);
		if (wd) {
			if (!Le10)
				Le10 = 1. / log(10.);
			g->dL = Le10 / t;
			g->dL2 = -g->dL / t;
			}
		g->O = rv.d;
		o += 3;
		goto top;
	  case OP_log101:
		rv.d = log10(t = w[o[3]]);
		if (errchk(rv))
			introuble("log10",t,1);
		r = (Eresult*)(w + o[2]);
		if (wd) {
			if (!Le10)
				Le10 = 1. / log(10.);
			r->dL = Le10 / t;
			r->dL2 = -r->dL / t;
			}
		r->O = rv.d;
		o += 4;
		goto top;
	  case OP_log0:
		rv.d = log(w[o[2]]);
		if (errchk(rv))
			introuble("log",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPlog_g:
		rv.d = log(t = w[o[2]]);
		if (errchk(rv))
			introuble("log",t,1);
		g = (GOps*)(w + o[1]);
		if (wd) {
			t = g->dL = 1. / t;
			g->dL2 = -t*t;
			}
		g->O = rv.d;
		o += 3;
		goto top;
	  case OP_log1:
		rv.d = log(t = w[o[3]]);
		if (errchk(rv))
			introuble("log",t,1);
		r = (Eresult*)(w + o[2]);
		if (wd) {
			r->dL = t1 = 1. / t;
			r->dL2 = -t1*t1;
			}
		r->O = rv.d;
		o += 4;
		goto top;
	  case OP_exp0:
		rv.d = exp(t = w[o[2]]);
		if (errchk(rv)) {
			if (t >= 0.)
				introuble("exp",t,1);
			else {
				errno_set(0);
				rv.d = 0.;
				}
			}
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OP_exp1:
		rv.d = exp(t = w[o[3]]);
		if (errchk(rv)) {
			if (t >= 0.)
				introuble("exp",t,1);
			else {
				errno_set(0);
				rv.d = 0.;
				}
			}
		w[o[2]] = rv.d;
		o += 4;
		goto top;
	  case OPexp_g:
		rv.d = exp(t = w[o[2]]);
		if (errchk(rv)) {
			if (t >= 0.)
				introuble("exp",t,1);
			else {
				errno_set(0);
				rv.d = 0.;
				}
			}
		g = (GOps*)(w + o[1]);
		g->O = g->dL = g->dL2 = rv.d;
		o += 3;
		goto top;
	  case OP_cosh0:
		rv.d = cosh(w[o[2]]);
		if (errchk(rv))
			introuble("cosh",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPcosh_g:
		rv.d = cosh(t = w[o[2]]);
		if (errchk(rv))
			introuble("cosh",t,1);
		g = (GOps*)(w + o[1]);
		if (wd) {
			tv.d = sinh(t);
			if (errchk(tv))
				introuble("cosh'", t, 2);
			g->dL = tv.d;
			g->dL2 = rv.d;
			}
		g->O = rv.d;
		o += 3;
		goto top;
	  case OP_cosh1:
		rv.d = cosh(t = w[o[3]]);
		if (errchk(rv))
			introuble("cosh",t,1);
		r = (Eresult*)(w + o[2]);
		if (wd) {
			tv.d = sinh(t);
			if (errchk(tv))
				introuble("cosh'", t, 2);
			r->dL = tv.d;
			r->dL2 = rv.d;
			}
		r->O = rv.d;
		o += 4;
		goto top;
	  case OP_cos0:
		rv.d = cos(w[o[2]]);
		if (errchk(rv))
			introuble("cos",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPcos_g:
		rv.d = cos(t = w[o[2]]);
		if (errchk(rv))
			introuble("cos",t,1);
		g = (GOps*)(w + o[1]);
		if (wd) {
			tv.d = sin(t);
			if (errchk(tv))
				introuble("cos'", t, 2);
			g->dL = -tv.d;
			g->dL2 = -rv.d;
			}
		g->O = rv.d;
		o += 3;
		goto top;
	  case OP_cos1:
		rv.d = cos(t = w[o[3]]);
		if (errchk(rv))
			introuble("cos",t,1);
		r = (Eresult*)(w + o[2]);
		if (wd) {
			tv.d = sin(t);
			if (errchk(tv))
				introuble("cos'", t, 2);
			r->dL = -tv.d;
			r->dL2 = -rv.d;
			}
		r->O = rv.d;
		o += 4;
		goto top;
	  case OP_atanh0:
		t = w[o[2]];
		r = (Eresult*)(w + o[1]);
		o += 3;
		if (t <= -1. || t >= 1.) {
 bad_atanh:
			errno_set(EDOM);
			rv.d = 0.;
			introuble("atanh",t,1);
			}
		rv.d = 0.5*log((1. + t) / (1. - t));
		if (errchk(rv))
			goto bad_atanh;
		r->O = rv.d;
		goto top;
	  case OPatanh_g:
		t = w[o[2]];
		g = (GOps*)(w + o[1]);
		o += 3;
		if (t <= -1. || t >= 1.) {
			r = (Eresult*)g;
			goto bad_atanh;
			}
		rv.d = 0.5*log((1. + t) / (1. - t));
		if (errchk(rv)) {
			r = (Eresult*)g;
			goto bad_atanh;
			}
		g->O = rv.d;
		if (wd) {
			g->dL = t1 = 1. / (1. - t*t);
			g->dL2 = (t+t)*t1*t1;
			}
		goto top;
	  case OP_atanh1:
		t = w[o[3]];
		r = (Eresult*)(w + o[2]);
		o += 4;
		if (t <= -1. || t >= 1.)
			goto bad_atanh;
		rv.d = 0.5*log((1. + t) / (1. - t));
		if (errchk(rv))
			goto bad_atanh;
		r->O = rv.d;
		if (wd) {
			r->dL = t1 = 1. / (1. - t*t);
			r->dL2 = (t+t)*t1*t1;
			}
		goto top;
	  case OP_atan20:
		rv.d = atan2(L = w[o[2]], R = w[o[3]]);
		if (errchk(rv))
			introuble2("atan2",L,R,1);
		w[o[1]] = rv.d;
		o += 4;
		goto top;
	  case OP_atan210:
		r = (Eresult*)(w + o[2]);
		r->O = rv.d = atan2(L = w[o[3]], R = w[o[4]]);
		if (errchk(rv))
			introuble2("atan2",L,R,1);
		if (wd) {
			t = 1. / (L*L + R*R);
			r->dL =  t * R;
			t *= t;
			t1 = L*R;
			r->dL2 = -(t * (t1+t1));
			}
		o += 5;
		goto top;
	  case OPatan210_g:
		g = (GOps*)(w + o[1]);
		g->O = rv.d = atan2(L = w[o[2]], R = w[o[3]]);
		if (errchk(rv))
			introuble2("atan2",L,R,1);
		if (wd) {
			t = 1. / (L*L + R*R);
			g->dL =  t * R;
			t *= t;
			t1 = L*R;
			g->dL2 = -(t * (t1+t1));
			}
		o += 4;
		goto top;
	  case OP_atan201:
		r = (Eresult*)(w + o[2]);
		r->O = rv.d = atan2(L = w[o[3]], R = w[o[4]]);
		if (errchk(rv))
			introuble2("atan2",L,R,1);
		if (wd) {
			t = 1. / (L*L + R*R);
			r->dL = -t * L;
			t *= t;
			t1 = L*R;
			r->dL2 = t * (t1+t1);
			}
		o += 5;
		goto top;
	  case OPatan201_g:
		g = (GOps*)(w + o[1]);
		g->O = rv.d = atan2(L = w[o[2]], R = w[o[3]]);
		if (errchk(rv))
			introuble2("atan2",L,R,1);
		if (wd) {
			t = 1. / (L*L + R*R);
			g->dL = -t * L;
			t *= t;
			t1 = L*R;
			g->dL2 = t * (t1+t1);
			}
		o += 4;
		goto top;
	  case OP_atan22:
		r = (Eresult*)(w + o[2]);
		r->O = rv.d = atan2(L = w[o[3]], R = w[o[4]]);
		if (errchk(rv))
			introuble2("atan2",L,R,1);
		if (wd) {
			t = 1. / (L*L + R*R);
			r->dL =  t * R;
			r->dR = -t * L;
			t *= t;
			t1 = L*R;
			r->dL2 = -(r->dR2 = t * (t1+t1));
			r->dLR = t * (L*L - R*R);
			}
		o += 5;
		goto top;
	  case OP_atan0:
		rv.d = atan(w[o[2]]);
		if (errchk(rv))
			introuble("atan",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPatan_g:
		g = (GOps*)(w + o[1]);
		g->O = rv.d = atan(t = w[o[2]]);
		if (errchk(rv))
			introuble("atan",t,1);
		if (wd) {
			g->dL = t1 = 1. / (1. + t*t);
			g->dL2 = -(t+t)*t1*t1;
			}
		o += 3;
		goto top;
	  case OP_atan1:
		r = (Eresult*)(w + o[2]);
		r->O = rv.d = atan(t = w[o[3]]);
		if (errchk(rv))
			introuble("atan",t,1);
		if (wd) {
			r->dL = t1 = 1. / (1. + t*t);
			r->dL2 = -(t+t)*t1*t1;
			}
		o += 4;
		goto top;
	  case OP_asinh0:
		t = t0 = w[o[2]];
		if ((sign = t < 0.))
			t = -t;
		rv.d = log(t + (t1 = sqrt(t*t + 1.)));
		if (errchk(rv))
			introuble("asinh",t0,1);
		if (sign)
			rv.d = -rv.d;
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPasinh_g:
		t = t0 = w[o[2]];
		rv.d = log(t + (t1 = sqrt(t0 = t*t + 1.)));
		if ((sign = t < 0.)) {
			t = -t;
			t0 = -t0;
			}
		if (errchk(rv))
			introuble("asinh",t0,1);
		if (sign)
			rv.d = -rv.d;
		g = (GOps*)(w + o[1]);
		g->O = rv.d;
		if (wd)
			g->dL2 = -(t/t0) * (g->dL = 1. / t1);
		o += 3;
		goto top;
	  case OP_asinh1:
		t = t0 = w[o[3]];
		rv.d = log(t + (t1 = sqrt(t0 = t*t + 1.)));
		if ((sign = t < 0.)) {
			t = -t;
			t0 = -t0;
			}
		if (errchk(rv))
			introuble("asinh",t0,1);
		if (sign)
			rv.d = -rv.d;
		r = (Eresult*)(w + o[2]);
		r->O = rv.d;
		if (wd)
			r->dL2 = -(t/t0) * (r->dL = 1. / t1);
		o += 4;
		goto top;
	  case OP_asin0:
		rv.d = asin(t = w[o[2]]);
		if (errchk(rv))
			introuble("asin",t,1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPasin_g:
		rv.d = asin(t = w[o[2]]);
		if (errchk(rv))
			introuble("asin",t,1);
		g = (GOps*)(w + o[1]);
		g->O = rv.d;
		if (wd) {
			if ((t1 = 1. - t*t) <= 0.)
				introuble("asin'",t,2);
			g->dL2 = t * (g->dL = 1. / sqrt(t1)) / t1;
			}
		o += 3;
		goto top;
	  case OP_asin1:
		rv.d = asin(t = w[o[3]]);
		if (errchk(rv))
			introuble("asin",t,1);
		r = (Eresult*)(w + o[2]);
		r->O = rv.d;
		if (wd) {
			if ((t1 = 1. - t*t) <= 0.)
				introuble("asin'",t,2);
			r->dL2 = t * (r->dL = 1. / sqrt(t1)) / t1;
			}
		o += 4;
		goto top;
	  case OP_acosh0:
		if ((t = w[o[2]]) < 1.) {
 bad_acosh:
			errno_set(EDOM);
			rv.d = t1 = 0.;
			introuble("acosh",t,1);
			}
		rv.d = log(t + sqrt(t*t - 1.));
		if (errchk(rv))
			goto bad_acosh;
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPacosh_g:
		if ((t = w[o[2]]) < 1.)
			goto bad_acosh;
		rv.d = log(t + (t1 = sqrt(t0 = t*t - 1.)));
		if (errchk(rv))
			goto bad_acosh;
		g = (GOps*)(w + o[1]);
		g->O = rv.d;
		if (wd) {
			if (t1 <= 0.)
				introuble("acosh'",t,2);
			g->dL2 = -t * (g->dL = 1. / t1) / t0;
			}
		o += 3;
		goto top;
	  case OP_acosh1:
		if ((t = w[o[3]]) < 1.)
			goto bad_acosh;
		rv.d = log(t + (t1 = sqrt(t0 = t*t - 1.)));
		if (errchk(rv))
			goto bad_acosh;
		r = (Eresult*)(w + o[2]);
		r->O = rv.d;
		if (wd) {
			if (t1 <= 0.)
				introuble("acosh'",t,2);
			r->dL2 = -t * (r->dL = 1. / t1) / t0;
			}
		o += 4;
		goto top;
	  case OP_acos0:
		rv.d = acos(t = w[o[2]]);
		if (errchk(rv))
			introuble("acos",t,1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case OPacos_g:
		rv.d = acos(t = w[o[2]]);
		if (errchk(rv))
			introuble("acos",t,1);
		g = (GOps*)(w + o[1]);
		g->O = rv.d;
		if (wd) {
			if ((t1 = 1. - t*t) <= 0.)
				introuble("acos'",t,2);
			g->dL2 = t * (g->dL = -1. / sqrt(t1)) / t1;
			}
		o += 3;
		goto top;
	  case OP_acos1:
		rv.d = acos(t = w[o[3]]);
		if (errchk(rv))
			introuble("acos",t,1);
		r = (Eresult*)(w + o[2]);
		r->O = rv.d;
		if (wd) {
			if ((t1 = 1. - t*t) <= 0.)
				introuble("acos'",t,2);
			r->dL2 = t * (r->dL = -1. / sqrt(t1)) / t1;
			}
		o += 4;
		goto top;
	  case OPSUMLIST1:
		++o;
	  case OPSUMLIST0:
		rp = w + o[1];
		o1 = o + 3;
		o = o1 + o[2];
		t = 0.;
		while(o1 < o)
			t += w[*o1++];
		*rp = t;
		goto top;
	  case OPintdiv:
		L = w[o[2]];
		if (!(R = w[o[3]]))
			zero_div(L, " div ");
		w[o[1]] = (L /= R) >= 0. ? floor(L) : ceil(L);
		o += 4;
		goto top;
	  case OP_precision:
		R = w[o[3]];
		g_fmtp(buf, w[o[2]], (int)R);
		w[o[1]] = strtod(buf, (char**)0);
		o += 4;
		goto top;
	  case OP_round:
		w[o[1]] = Round(w[o[2]], w[o[3]]);
		o += 4;
		goto top;
	  case OP_trunc:
		L = w[o[2]];
		if (!(R = w[o[3]]))
			rv.d = L >= 0. ? floor(L) : ceil(L);
		else {
			rv.d = Round(L, prec = (int)R);
			if (rv.d != L) {
				R = 0.5*mypow(10., (real)-prec);
				rv.d = Round(L > 0 ? L - R : L + R, prec);
				}
			}
		w[o[1]] = rv.d;
		o += 4;
		goto top;
	  case n_OPCOUNT:
		n = o[2];
		j = 3;
		rv.d = 0.;
		for(k = j + n; j < k; ++j)
			if (w[o[j]] != 0.)
				++rv.d;
		w[o[1]] = rv.d;
		o += k;
		goto top;
	  case n_OPNUMBEROF:
	  case n_OPNUMBEROFs:	/* TEMPORARY */
		n = o[2];
		t = w[o[3]];
		j = 4;
		rv.d = 0.;
		for(k = j + n; j < k; ++j)
			if (w[o[j]] == t)
				++rv.d;
		w[o[1]] = rv.d;
		o += k;
		goto top;
#ifdef X64_bit_pointers
	  case OP_PLTERM0align:
		pp = (plterm**)(o+4);
		goto more_plterm0;
#endif
	  case n_OPPLTERM0:
		pp = (plterm**)(o+3);
alignarg(more_plterm0:)
		i = o[1];
		R = w[o[2]];
		p = *pp;
		o = (int*)&pp[1];
		n = p->n;
		z = p->z;
		bs = p->bs;
		bs += z;
		z >>= 1;
		if (R >= 0) {
			n -= z;
			if (n <= 1 || R <= bs[1]) {
				w[i] = R**bs;
				goto top;
				}
			for(t = bs[0]*bs[1]; --n > 1 && R > bs[3]; bs += 2)
				t += (bs[3]-bs[1])*bs[2];
			w[i] = t + (R-bs[1])*(w[i+1] = bs[2]);
			goto top;
			}
		if (z <= 0) {
			w[i] = R*(w[i+1] = bs[0]);
			goto top;
			}
		for(t = bs[0]*bs[-1]; --z > 0 && R < bs[-3]; bs -= 2)
			t += bs[-2]*(bs[-3] - bs[-1]);
		w[i] = t + (R - bs[-1])*bs[-2];
		goto top;
#ifdef X64_bit_pointers
	  case OP_PLTERM1align:
		pp = (plterm**)(o+5);
		goto more_plterm;
#endif
	  case n_OPPLTERM1:
		pp = (plterm**)(o+4);
alignarg(more_plterm:)
		r = (Eresult*)(w + o[2]);
		R = w[o[3]];
		p = *pp;
		o = (int*)&pp[1];
		n = p->n;
		z = p->z;
		bs = p->bs + z;
		z >>= 1;
		if (R >= 0) {
			n -= z;
			if (n <= 1 || R <= bs[1]) {
				r->O = R*(r->dL = *bs);
				goto top;
				}
			for(t = bs[0]*bs[1]; --n > 1 && R > bs[3]; bs += 2)
				t += (bs[3]-bs[1])*bs[2];
			r->O = t + (R-bs[1])*(r->dL = bs[2]);
			goto top;
			}
		if (z <= 0) {
			r->O = R*(r->dL = bs[0]);
			goto top;
			}
		for(t = bs[0]*bs[-1]; --z > 0 && R < bs[-3]; bs -= 2)
			t += bs[-2]*(bs[-3] - bs[-1]);
		r->O = t + (R - bs[-1])*(r->dL = bs[-2]);
		goto top;
	  case n_OPANDLIST:
		k = 3 + o[2];
		for(i = 3; i < k; ++i) {
			if (!w[o[i]]) {
				w[o[1]] = 0.;
				o += k;
				goto top;
				}
			}
		w[o[1]] = 1.;
		o += k;
		goto top;
	  case n_OPORLIST:
		k = 3 + o[2];
		for(i = 3; i < k; ++i) {
			if (w[o[i]]) {
				w[o[1]] = 1.;
				o += k;
				goto top;
				}
			}
		w[o[1]] = 0.;
		o += k;
		goto top;
	  case n_OP_IFF:
		i = w[o[2]] != 0.;
		j = w[o[3]] != 0.;
		w[o[1]] = i == j ? 1. : 0.;
		o += 4;
		goto top;
	  case n_OPALLDIFF:
		n = o[2];
		bs = rbuf;
		if (n > sizeof(rbuf)/sizeof(real))
			bs = (real*)Malloc(n*sizeof(real));
		for(i = 0; i < n; ++i)
			bs[i] = w[o[i+3]];
		t = 1.;
		if (setjmp(J.jb)) {
			t = 0.;
			goto alldiff_done;
			}
		qsortv(bs, n, sizeof(real), rcompj, &J);
 alldiff_done:	if (bs != rbuf)
			free(bs);
		w[o[1]] = t;
		o += n + 3;
		goto top;
	  case n_OPSOMESAME:
		n = o[2];
		bs = rbuf;
		if (n > sizeof(rbuf)/sizeof(real))
			bs = (real*)Malloc(n*sizeof(real));
		for(i = 0; i < n; ++i)
			bs[i] = w[o[i+3]];
		t = 0.;
		if (setjmp(J.jb)) {
			t = 1.;
			goto alldiff_done;
			}
		qsortv(bs, n, sizeof(real), rcompj, &J);
		goto alldiff_done;
	  case OP_2POW0:
		L = w[o[2]];
		w[o[1]] = L*L;
		o += 3;
		goto top;
	  case OP2POW_g:
		g = (GOps*)(w + o[1]);
		L = w[o[2]];
		r = (Eresult*)(w + o[1]);
		g->O = L*L;
		g->dL = L + L;
		o += 3;
		goto top;
	  case OP_2POW1:
		L = w[o[3]];
		r = (Eresult*)(w + o[2]);
		r->O = L*L;
		r->dL = L + L;
		o += 4;
		goto top;
#ifdef X64_bit_pointers
	  case OPCPOW0align:
		bs = (real*)&o[4];
		goto more_CPOW0;
	  case OPCPOW1align:
		bs = (real*)&o[5];
		goto more_CPOW1;
#endif
	  case nOPCPOW0:
		bs = (real*)&o[3];
alignarg(more_CPOW0:)
		o1 = (int*)&bs[1];
		wdf = 0;
		goto more_CPOW;
	  case nOPCPOW1:
		bs = (real*)&o[4];
alignarg(more_CPOW1:)
		o1 = (int*)&bs[1];
		wdf = wd;
		++o;
 more_CPOW:
		L = w[o[2]];
		R = *bs;
		r = (Eresult*)(w + o[1]);
		o = o1;
		rv.d = mypow(L, R);
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				r->O = r->dL2 = Infinity;
				r->dL2 = negInfinity;
				goto top;
				}
#endif
			introuble2("pow",L,R,1);
			}
		r->O = rv.d;
		if (wd) {
			if (L > 0.) {
				t1 = rv.d / L;
				r->dL = R*t1;
				r->dL2 = (R - 1.) * (r->dL / L);
				goto top;
				}
			r->dL = r->dL2 = 0.;
			if (L != 0.)
				goto badpowder;
			if (R > 1.) {
				if (R == 2.)
					r->dL2 = 2;
				else if (R < 2.) {
					r->dL2 = Infinity;
#ifndef WANT_INFNAN
					introuble2("pow\"",L,R,3);
#endif
					}
				goto top;
				}
			if (R == 1.)
				r->dL = 1.;
			else { /* 0 < R < 1 */
#ifdef WANT_INFNAN
				r->dL = Infinity;
				r->dL2 = negInfinity;
#else
				goto badpowder;
#endif
				}
			}
		r->O = rv.d;
		goto top;

#ifdef X64_bit_pointers
	  case OP_FUNCALL0align:
		ptfi = (tfinfo**)(o+3);
		wdf = 0;
		goto more_func;

	  case OP_FUNCALL1align:
		ptfi = (tfinfo**)(o+4);
		goto more_func1;
#endif
	  case OP_FUNCALL1:
		ptfi = (tfinfo**)(o+3);
alignarg(more_func1:)
		++o;
		wdf = wd;
		goto more_func;
	  case OP_FUNCALL0:
		ptfi = (tfinfo**)(o+2);
		wdf = 0;
 more_func:
		k = o[1];
		tfi = *ptfi;
		al = ew->al + tfi->ali;
		o = (int*)(ptfi+1);
		rp = al->ra;
		if ((nr = al->nr) > 0) {
			for(i = 0; i < nr; ++i)
				*rp++ = w[*o++];
			}
		if ((ns = al->nsin) > 0) {
			sa = al->sa;
			for(i = 0; i < ns; ++i)
				*sa++ = *(char**)&w[*o++];
			}
		rp = w + k + 4;
		if (wdf && nr) {
			memset(al->derivs = rp, 0, sizeof(real)*nr*(nr+3)/2);
			al->hes = rp + nr;
			al->dig = tfi->wd;
			}
		else {
			al->derivs = al->hes = 0;
			al->dig = 0;
			}
		al->Errmsg = 0;
		T.u.prev = 0;
		al->TMI = &T;
		fi = tfi->fi;
		w[k] = (*fi->funcp)(al);
		errno_set(0);
		if ((s = al->Errmsg))
			fintrouble_ASL(ew, fi, s, &T);
		for(T1 = T.u.prev; T1; T1 = T1prev) {
			T1prev = T1->u.prev;
			free(T1);
			}
		goto top;

	  case n_OPHOL:
		*(char**)&w[o[1]] = (char*)(o+3);
		o += o[2];
		goto top;

	  /*case nOPVARVAL: accessed directly */
	  case OPCOPY0:
		w[o[1]] = w[o[2]];
		o += 3;
		goto top;

	  case OPCOPY1:
	  case OPCOPY1a:
		w[o[2]] = w[o[3]];
		o += 4;
		goto top;

	  case OP_COPYSYM:
		*(char**)&w[o[1]] = *(char**)&w[o[2]];
		o += 3;
		goto top;

	  case OP_GOTO:
	  case OPGOTO2:
	  case OP_NEXTBLK:
		o = *(int**)(o+1);	/* for chaining blocks of code */
		goto top;

	  case OPGOTOF:
		o = (int*)&((int**)(o+1))[1];
		goto top;

	  case OPGOTOF2:
	  case OPGOTOF2n:
		pop = (int**)&o[1];
		o = (int*)&pop[2];
		goto top;

	  case OPGOTOMM:
		o += 2;
		goto top;

	  case OPRETB:
		++o;
		goto top;

#ifdef X64_bit_pointers
	  case OP_GOTOalign:
	  case OPGOTO2align:
	  case OP_NEXTBLKalign:
		o = *(int**)(o+2);
		goto top;

	  case OPGOTOFalign:
		o = (int*)&((int**)(o+2))[1];
		goto top;

	  case OPGOTOF2align:
	  case OPGOTOF2nalign:
		pop = (int**)&o[2];
		o = (int*)&pop[2];
		goto top;
	  case OPGOTOBalign:
		++o;
#endif
	  case OPGOTOB:
		pop = (int**)(o+1);
		o = (int*)&pop[1];
		goto top;

	  case OPVARREF:
		o += 3;
		goto top;

	  default:
		fprintf(Stderr, "\nUnexpected opno %d in eval2_ASL()\n", *o);
		fflush(Stderr);
		exit(1);
	  }
 done:
	return rv.d;
	}

#ifdef DERPDEBUG
/*DEBUG*/ /*extern*/ int derpzork, derpzork1, derpzork2, derpzork3;
#endif

 static void
derpropa(derpblock *db, uint a0, real *s, real *w, real f)
{
	derp *d, *de;
	real t;
	size_t n;
	uint a;

#ifdef DERPDEBUG
	printf("\nStarting derpropa call %d\n", ++derpzork2);
	if (derpzork2 == derpzork3)
		printf("derpzork2 = %d\n", derpzork2);
#endif

	d = db->d0;
	s[d->b] = f;
	for(;;) {
		for(de = db->de; d < de; ++d) {
#ifdef DERPDEBUG
/*DEBUG*/		if (++derpzork == derpzork1)
/*DEBUG*/			printf("derpzork = %d\n", derpzork);
#endif
			t = s[d->b] * w[d->c];
			if ((a = d->a) >= a0) {
#ifdef DERPDEBUG
/*DEBUG*/			printf("%d\ts[%d] := s[%d]*w[%d] = %g * %g = %g\n",
/*DEBUG*/				derpzork, d->a, d->b, d->c, s[d->b], w[d->c], t);
#endif
				s[a] = t;
				}
			else {
				s[a] += t;
#ifdef DERPDEBUG
/*DEBUG*/			printf("%d\ts[%d] += s[%d]*w[%d] = %g * %g ==> %g\n",
/*DEBUG*/				derpzork, d->a, d->b, d->c, s[d->b], w[d->c], s[a]);
#endif
				}
			}
		if (!(n = db->nxt))
			break;
		if (!(db = db->next))
			db = *(derpblock**)&w[n];
		d = db->d0;
		}
	}

 static real
copeval(EvalWorkspace *ew, ps_func *f)
{
	psb_elem *b, *be;
	real t;

	t = 0.;
	for(b = f->pi.b, be = f->pi.be; b < be; b++)
		t += eval2_ASL(b->o.e, ew);
	return t;
	}

 static real
cogeval(EvalWorkspace *ew, ps_func *f)
{
	Varval *V;
	psb_elem *b, *be;
	psg_elem *g, *ge;
	lincoef *lc, *lce;
	linpart *L;
	real t, t1, *w;


	V = (Varval*)(w = ew->w);
	t = 0.;
	for(g = f->g, ge = f->ge; g < ge; g++) {
		t1 = g->g0;
		if ((L = g->L)) {
			lc = L->lc;
			for(lce = lc + L->n; lc < lce; ++lc)
				t1 += lc->coef * V[lc->varno].O;
			}
		if ((b = g->pi.b)) {
			for(be = g->pi.be; b < be; b++)
				t1 += eval2_ASL(b->o.e, ew);
			}
		w[g->gm] = t1;
		t += g->scale * eval2_ASL(g->o, ew);
		}
	return t;
	}

 static void
psderprop(EvalWorkspace *ew, Psbinfo *pi)
{
	ASL_pfgh *asl;
	derp *d, *de;
	derpblock *db, **pdb;
	int *ce, *cee, *ov, *ove;
	linarg *la, **lap;
	real *oc, *s, *sdv, t, *w;
	size_t n;
	uint a, a0;

	asl = (ASL_pfgh*)ew->asl;
	a0 = asl->i.maxvar;
	w = ew->w;
	s = ew->derps;
	sdv = s + asl->i.defvar0;
	if ((ce = pi->ce)) {
		cee = ce + *ce;
		do sdv[*++ce] = 0.;
		   while(ce < cee);
		}
	pdb = pi->pdb;
	db = *pdb;
	d = db->d0;
	s[d->b] = 1.;
#ifdef DERPDEBUG
	printf("\nStarting psderprop call %d\n", ++derpzork2);
	if (derpzork2 == derpzork3)
		printf("derpzork2 = %d\n", derpzork2);
#endif
	for(;;) {
		for(;;) {
			d = db->d0;
			for(de = db->de; d < de; ++d) {
#ifdef DERPDEBUG
/*DEBUG*/		if (++derpzork == derpzork1)
/*DEBUG*/			printf("derpzork = %d\n", derpzork);
#endif
				t = s[d->b] * w[d->c];
				if ((a = d->a) >= a0) {
#ifdef DERPDEBUG
/*DEBUG*/			printf("%d\ts[%d] := s[%d]*w[%d] = %g * %g = %g\n",
/*DEBUG*/				derpzork, d->a, d->b, d->c, s[d->b], w[d->c], t);
#endif
					s[a] = t;
					}
				else {
					s[a] += t;
#ifdef DERPDEBUG
/*DEBUG*/			printf("%d\ts[%d] += s[%d]*w[%d] = %g * %g ==> %g\n",
/*DEBUG*/				derpzork, d->a, d->b, d->c, s[d->b], w[d->c], s[a]);
#endif
					}
				}
			if (!(n = db->nxt))
				break;
			if (!(db = db->next))
				db = *(derpblock**)&w[n];
			}
		if (!(db = *++pdb))
			break;
		}
	if ((lap = pi->lap)) {
		while((la = *lap++))
			if ((t = s[la->u.v])) {
				oc = la->oc;
				ov = la->ov;
				for(ove = ov + la->nnz; ov < ove; ++ov, ++oc) {
#ifdef DERPDEBUG
/*DEBUG*/				if (++derpzork == derpzork1)
/*DEBUG*/					printf("derpzork = %d\n", derpzork);
/*DEBUG*/				printf("%d\ts[%d] = %g += %g*%g ==> %g\n",
/*DEBUG*/					derpzork, *ov, s[*ov], t, *oc, s[*ov] + t**oc);
#endif
					s[*ov] += t * *oc;
					}
				}
		}
	}

 static void
psgcomp(EvalWorkspace *ew, ps_func *f)
{
	GOps *u, *u1;
	int i, j, n, nov, *o, *ov;
	lincoef *lc, *lce;
	linpart *L;
	psg_elem *g, *ge;
	real *oc, *s, t, t1, t2, *us, *w;

	w = ew->w;
	s = ew->derps;
	ew->npsgcomp++;
	us = ew->unopscr;
	for(g = f->g, ge = f->ge; g < ge; g++) {
		nov = g->nov;
		ov = g->ov;
		u = (GOps*)(w + g->gm);
		for(j = 0; j < nov; ++j)
			s[ov[j]] = 0;
		if ((L = g->L)) {
			for(lc = L->lc, lce = lc + L->n; lc < lce; ++lc)
				s[lc->varno] = lc->coef;
			}
		if (g->pi.pdb)
			psderprop(ew, &g->pi);
		oc = (real*)(u + 1);
		for(j = 0; j < nov; ++j)
			oc[j] = s[ov[j]];

		/* compute u->dL and u->dL2 */

		o = g->g.i;
		u1 = (GOps*)(w + *o);
		t = u1->dL;
		if ((n = g->nu) == 1)
			t2 = u1->dL2;
		else {
			us[0] = 1.;
			i = 1;
			do {
				us[i] = t;
				u1 = (GOps*)&w[o[i]];
				t *= u1->dL;
				}
				while(++i < n);
			t1 = u1->dL;
			t2 = u1->dL2 * us[--i];
			for(;;) {
				u1 = (GOps*)&w[o[--i]];
				t2 += us[i] * u1->dL2 * (t1*t1);
				if (i <= 0)
					break;
				t1 *= u1->dL;
				}
			}
		u->dL = t *= g->scale;
		u->dL2 = t2 * g->scale;
		}
	}

 static void
addgr(EvalWorkspace *ew, ps_func *p)
{
	GOps *u;
	int i, n, *ov;
	psg_elem *g, *ge;
	real *oc, *s, t, *w;

	s = ew->derps;
	w = ew->w;
	g = p->g;
	ge = p->ge;
	do {
		u = (GOps*)(w + g->gm);
		if ((t = u->dL)) {
			oc = (real*)(u + 1);
			ov = g->ov;
			n = g->nov;
			for(i = 0; i < n; ++i)
				s[ov[i]] += t*oc[i];
			}
		} while(++g < ge);
	}

 void
dv_comp_ASL(EvalWorkspace *ew, int n0, int n)
{
	ASL_pfgh *asl;
	Varval *V, *Vd;
	cexp *c, *c1;
	ograd *og;
	int *dvsp0, i, j, k, n1, *ndvsp;
	lincoef *lc, *lce;
	linpart *L;
	real t;

	asl = (ASL_pfgh*)ew->asl;
	dvsp0 = asl->P.dvsp0;
	ndvsp = asl->P.ndvsp;
	c = cexps;
	V = (Varval*)ew->w;
	Vd = (Varval*)ew->dv;
	for(i = n0, n1 = n0 + n; i < n1; ++i) {
		if ((k = ndvsp[i])) {
			k += j = dvsp0[i];
			do {
				c1 = c + j;
				ew->cv_index = j + 1;
				V[c1->varno].O = eval2_ASL(c1->o.e, ew);
				} while(++j < k);
			}
		c1 = c + i;
		ew->cv_index = i + 1;
		t = eval2_ASL(c1->o.e, ew);
		if ((L = c1->lp)) {
			for(lc = L->lc, lce = lc + L->n; lc < lce; ++lc)
				t += lc->coef * V[lc->varno].O;
			}
		else if (!c1->db && (og = asl->P.dv[i].ll)) {
			if (og->varno < 0) {
				t += og->coef;
				og = og->next;
				}
			while(og) {
				t += og->coef*V[og->varno].O;
				og = og->next;
				}
			}
		Vd[i].O = t;
		}
	ew->cv_index = 0;
	}

 static void
funnelset(EvalWorkspace *ew, cexp **dvf)
{
	ASL_pfgh *asl;
	cexp *ce;
	derp *d, *de;
	derpblock *db;
	int *vr, *vre;
	real *s, *w;
	size_t n;
	uint a;

	asl = (ASL_pfgh*)ew->asl;
	a = asl->i.maxvar;
	s = ew->derps;
	w = ew->w;
	for(; (ce = *dvf); ++dvf) {
		if ((vr = ce->vref)) {
			n = vr[0];
			vr += 3;
			vre = vr + n;
			do s[*vr] = 0.;
			   while(++vr < vre);
			}
		derpropa(ce->dbf, a, s, w, 1.);
		db = ce->db;
		d = db->d0;
		for(;;) {
			for(de = db->de; d < de; ++d)
				w[d->c] = s[d->a];
			if (!(n = db->nxt))
				break;
			if (!(db = db->next))
				db = *(derpblock**)&w[n];
			d = db->d0;
			}
		}
	}

 void
conpval_ew_ASL(EvalWorkspace *ew, real *X, real *F, fint *nerror)
{
	ASL *a;
	ASL_pfgh *asl;
	Jmp_buf err_jmp0;
	cgrad *gr, **gr0;
	int *c1, *cm, i, j, j1, je, kv, nc1, *o, *vmi;
	ps_func *p, *p0;
	real *cscale, f, *vscale;
	size_t *ncxval, nx;

	asl = (ASL_pfgh*)(a = ew->asl);
	ASL_CHECK(a, ASL_read_pfgh, "conpval");
	if (nerror && *nerror >= 0) {
		ew->err_jmpw = &err_jmp0;
		i = setjmp(err_jmp0.jb);
		if ((*nerror = i))
			return;
		}
	ew->wantderiv = want_derivs;
	errno = 0;	/* in case f77 set errno opening files */
	++ew->stats.conval;
	cm = asl->i.cmap;
	j = n_conjac[0];
	if (!(ew->x0kind & ASL_x_known)) {
		ew->co_index = cm ? cm[j] : j;
		xp_check_ASL(ew, X);
		}
	if (ew->x0kind & ASL_need_concom) {
		dv_comp_ASL(ew, comb, comc);
		ew->x0kind &= ~ASL_need_concom;
		}
	je = n_conjac[1];
	if ((c1 = asl->i.c_cexp1st_) && (nc1 = c1[je] - c1[j]))
		dv_comp_ASL(ew, c1[j], nc1);
	if (!(gr0 = asl->i.Cgrad0))
		asl->i.Cgrad0 = gr0 = asl->i.Cgrad_;
	p0 = asl->P.cps;
	cscale = asl->i.cscale;
	kv = 0;
	vmi = 0;
	if ((vscale = asl->i.vscale))
		kv = 2;
	if (asl->i.vmap) {
		vmi = get_vminv_ASL(a);
		++kv;
		}
	nx = ew->nxval;
	ncxval = ew->ncxval;
	for(; j < je; j++) {
		i = j;
		if (cm)
			i = cm[i];
		ew->co_index = i;
		p = p0 + i;
		if (p->pi.b) {
			f = copeval(ew, p);
			if (p->g)
				f += cogeval(ew, p);
			}
		else if (p->g)
			f = cogeval(ew, p);
		else if ((o = con_de[i].o.e))
			f = ew->w[o[1]];
		else
			f = 0.;
		ncxval[i] = nx;
		if (!F)
			continue;
		gr = gr0[i];
		switch(kv) {
		  case 3:
			for(; gr; gr = gr->next) {
				j1 = vmi[gr->varno];
				f += X[j1] * vscale[j1] * gr->coef;
				}
			break;
		  case 2:
			for(; gr; gr = gr->next) {
				j1 = gr->varno;
				f += X[j1] * vscale[j1] * gr->coef;
				}
			break;
		  case 1:
			for(; gr; gr = gr->next)
				f += X[vmi[gr->varno]] * gr->coef;
			break;
		  case 0:
			for(; gr; gr = gr->next)
				f += X[gr->varno] * gr->coef;
		  }
		if (cscale)
			f *= cscale[j];
		*F++ = f;
		}
	ew->x0kind |= ASL_have_conval;
	ew->err_jmpw = 0;
	}

 void
jacpval_ew_ASL(EvalWorkspace *ew, real *X, real *G, fint *nerror)
{
	ASL *a;
	ASL_pfgh *asl;
	Jmp_buf err_jmp0;
	cgrad *gr, **gr0, *gr1;
	fint ne0;
	int *cm, i, i1, j, k, *vmi;
	ps_func *p, *p0;
	real *cscale, *s, t, *vscale;
	size_t nx;
	static char who[] = "jacpval";

	asl = (ASL_pfgh*)(a = ew->asl);
	ASL_CHECK(a, ASL_read_pfgh, who);
	if (!want_derivs)
		No_derivs_ASL(who);
	ne0 = -1;
	s = ew->derps;
	if (nerror && (ne0 = *nerror) >= 0) {
		ew->err_jmpw = &err_jmp0;
		i = setjmp(err_jmp0.jb);
		if ((*nerror = i))
			return;
		}
	errno = 0;	/* in case f77 set errno opening files */
	++ew->stats.jacval;
	cm = asl->i.cmap;
	j = n_conjac[0];
	ew->co_index = cm ? cm[j] : j;
	if ((!(ew->x0kind & ASL_x_known) && xp_check_ASL(ew,X))
	 || !(ew->x0kind & ASL_have_conval)) {
		if (!(ew->x0kind & ASL_x_known)) {
			ew->x0kind |= ASL_x_known;
			conpval_ew_ASL(ew,X,0,nerror);
			ew->x0kind &= ~ASL_x_known;
			}
		else
			conpval_ew_ASL(ew,X,0,nerror);
		if (ne0 >= 0 && *nerror)
			return;
		}
	nx = ew->nxval;
	gr0 = asl->i.Cgrad0;
	cscale = asl->i.cscale;
	vscale = asl->i.vscale;
	vmi = 0;
	if (asl->i.vmap)
		vmi = get_vminv_ASL(a);
	k = n_conjac[1];
	if (ew->Derrs)
		deriv_errchk_ASL(ew, j, k-j, 2);
	if (ew->x0kind & ASL_need_comba) {
		funnelset(ew, asl->I.dvfb);
		ew->x0kind &= ~ASL_need_comba;
		}
	if (ew->x0kind & ASL_need_comca) {
		funnelset(ew, asl->I.dvfc);
		ew->x0kind &= ~ASL_need_comca;
		}
	p0 = asl->P.cps;
	for(; j < k; ++j) {
		i = j;
		if (cm)
			i = cm[i];
		p = p0 + i;
		ew->ncxval[i] = nx;
		if (p->g)
			psgcomp(ew, p);
		for(gr = gr1 = gr0[i]; gr; gr = gr->next)
			s[gr->varno] = gr->coef;
		if (p->pi.pdb)
			psderprop(ew, &p->pi);
		if (p->g)
			addgr(ew, p);
		if (vscale) {
			gr = gr1;
			if (vmi)
				for(; gr; gr = gr->next) {
					i1 = gr->varno;
					s[i1] *= vscale[vmi[i1]];
					}
			else
				for(; gr; gr = gr->next) {
					i1 = gr->varno;
					s[i1] *= vscale[i1];
					}
			}
		gr = gr1;
		if (cscale)
			for(t = cscale[j]; gr; gr = gr->next)
				G[gr->goff] = t*s[gr->varno];
		else
			for(; gr; gr = gr->next)
				G[gr->goff] = s[gr->varno];
		}
	ew->err_jmpw = 0;
	}

 int
jacpdim_ASL(ASL *asl, const char *stub, fint *M, fint *N, fint *NO, fint *NZ,
	fint *MXROW, fint *MXCOL, ftnlen stub_len)
{
	FILE *nl;

	nl = jac_dim_ASL(asl, stub, M, N, NO, NZ, MXROW, MXCOL, stub_len);
	if (!nl)
		return ASL_readerr_nofile;
	X0 = (real *)M1alloc(n_var*sizeof(real));
	return pfgh_read_ASL(asl, nl, ASL_findgroups | ASL_return_read_err);
	}

/******** objpval, objpgrd ********/

 static void
NNOBJ_chk(ASL *asl, int i, const char *who)
{
	ASL_CHECK(asl, ASL_read_pfgh, who);
	if (i < 0 || i >= n_obj) {
		fprintf(Stderr,
			"%s: got NOBJ = %d; expected 0 <= NOBJ < %d\n",
			who, i, n_obj);
		exit(1);
		}
	}

 real
objpval_ew_ASL(EvalWorkspace *ew, int i, real *X, fint *nerror)
{
	ASL *a;
	ASL_pfgh *asl;
	Jmp_buf err_jmp0;
	int *c1, ij, j1, kv, nc1, *o, *vmi;
	ograd *gr;
	ps_func *p;
	real f, *vscale;

	asl = (ASL_pfgh*)(a = ew->asl);
	NNOBJ_chk(a, i, "objpval");
	if (nerror && *nerror >= 0) {
		ew->err_jmpw = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			return 0.;
		}
	ew->wantderiv = want_derivs;
	errno = 0;	/* in case f77 set errno opening files */
	++ew->stats.objval;
	ew->co_index = -(i + 1);
	if (!(ew->x0kind & ASL_x_known))
		xp_check_ASL(ew, X);
	if (ew->x0kind & ASL_need_objcom) {
		dv_comp_ASL(ew, combc, como);
		ew->x0kind &= ~ASL_need_objcom;
		}
	if ((c1 = asl->i.o_cexp1st_) && (nc1 = c1[i+1] - c1[i]))
		dv_comp_ASL(ew, c1[i], nc1);
	p = asl->P.ops + i;
	if (p->pi.b) {
		f = copeval(ew, p);
		if (p->g)
			f += cogeval(ew, p);
		}
	else if (p->g)
		f = cogeval(ew, p);
	else if ((o = obj_de[i].o.e))
		f = ew->w[o[1]];
	else
		f = 0.;
	ew->noxval[i] = ew->nxval;
	gr = Ograd[i];
	kv = 0;
	vmi = 0;
	if ((vscale = asl->i.vscale))
		kv = 2;
	if (asl->i.vmap) {
		vmi = get_vminv_ASL(a);
		++kv;
		}
	switch(kv) {
	  case 3:
		for(; gr; gr = gr->next) {
			j1 = vmi[gr->varno];
			f += X[j1] * vscale[j1] * gr->coef;
			}
		break;
	  case 2:
		for(; gr; gr = gr->next) {
			j1 = gr->varno;
			f += X[j1] * vscale[j1] * gr->coef;
			}
		break;
	  case 1:
		for(; gr; gr = gr->next)
			f += X[vmi[gr->varno]] * gr->coef;
		break;
	  case 0:
		for(; gr; gr = gr->next)
			f += X[gr->varno] * gr->coef;
	  }
	ew->err_jmpw = 0;
	return f;
	}

 void
objpgrd_ew_ASL(EvalWorkspace *ew, int i, real *X, real *G, fint *nerror)
{
	ASL *a;
	ASL_pfgh *asl;
	Jmp_buf err_jmp0;
	fint ne0;
	int ij, j, *vmi, *z;
	ograd *gr, *gr0;
	ps_func *p;
	real *s, *vscale;
	static char who[] = "objpgrd";

	asl = (ASL_pfgh*)(a = ew->asl);
	NNOBJ_chk(a, i, who);
	if (!want_derivs)
		No_derivs_ASL(who);
	p = asl->P.ops + i;
	ne0 = -1;
	if (nerror && (ne0 = *nerror) >= 0) {
		ew->err_jmpw = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			return;
		}
	errno = 0;	/* in case f77 set errno opening files */
	++ew->stats.objgrd;
	if (!(ew->x0kind & ASL_x_known)) {
		ew->co_index = -(i + 1);
		xp_check_ASL(ew, X);
		}
	if (ew->noxval[i] != ew->nxval) {
		if (!(ew->x0kind & ASL_x_known)) {
			ew->x0kind |= ASL_x_known;
			objpval_ew_ASL(ew, i, X, nerror);
			ew->x0kind &= ~ASL_x_known;
			}
		else
			objpval_ew_ASL(ew, i, X, nerror);
		if (ne0 >= 0 && *nerror)
			return;
		}
	if (ew->Derrs)
		deriv_errchk_ASL(ew, -(i+1), 1, 2);
	if (ew->x0kind & ASL_need_comba) {
		funnelset(ew, asl->I.dvfb);
		ew->x0kind &= ~ASL_need_comba;
		}
	if (ew->x0kind & ASL_need_comoa) {
		funnelset(ew, asl->I.dvfo);
		ew->x0kind &= ~ASL_need_comoa;
		}
	if (p->g)
		psgcomp(ew, p);
	gr0 = Ograd[i];
	s = ew->derps;
	for(gr = gr0; gr; gr = gr->next)
		s[gr->varno] = gr->coef;
	if (p->pi.pdb)
		psderprop(ew, &p->pi);
	if (p->g)
		addgr(ew, p);
	ew->noxgval[i] = ew->nxval;
	if (!G)
		return;
	if (zerograds) {	/* sparse gradients */
		z = zerograds[i];
		while((j = *z++) >= 0)
			G[j] = 0;
		}
	gr = gr0;
	vmi = 0;
	if (asl->i.vmap)
		vmi = get_vminv_ASL(a);
	if ((vscale = asl->i.vscale)) {
		if (vmi)
			for(; gr; gr = gr->next) {
				j = vmi[i = gr->varno];
				G[j] = s[i] * vscale[j];
				}
		else
			for(; gr; gr = gr->next) {
				i = gr->varno;
				G[i] = s[i] * vscale[i];
				}
		}
	else if (vmi)
		for(; gr; gr = gr->next) {
			i = gr->varno;
			G[vmi[i]] = s[i];
			}
	else
		for(; gr; gr = gr->next) {
			i = gr->varno;
			G[i] = s[i];
			}
	ew->err_jmpw = 0;
	}

 static void
INchk(ASL *asl, const char *who, int i, int ix)
{
	ASL_CHECK(asl, ASL_read_pfgh, who);
	if (i < 0 || i >= ix) {
		fprintf(Stderr, "%s: got I = %d; expected 0 <= I < %d\n",
			who, i, ix);
		exit(1);
		}
	}

 static real
cpval(EvalWorkspace *ew, int i, real *X, fint *nerror)
{
	ASL_pfgh *asl;
	Jmp_buf err_jmp0;
	int L, *c1, nc, nc1, *o;
	ps_func *p;
	real f;

	asl = (ASL_pfgh*)ew->asl;
	if (nerror && *nerror >= 0) {
		ew->err_jmpw = &err_jmp0;
		L = setjmp(err_jmp0.jb);
		if ((*nerror = L))
			return 0.;
		}
	ew->wantderiv = want_derivs;
	errno = 0;	/* in case f77 set errno opening files */
	ew->co_index = i;
	if (!(ew->x0kind & ASL_x_known))
		xp_check_ASL(ew, X);
	if (ew->x0kind & ASL_need_concom) {
		dv_comp_ASL(ew, comb, comc);
		ew->x0kind &= ~ASL_need_concom;
		}
	if ((c1 = asl->i.c_cexp1st_) && (nc1 = c1[i+1] - c1[i]))
		dv_comp_ASL(ew, c1[i], nc1);
	if (i >= (nc = asl->i.n_con0)) {
		o = lcon_de[i-nc].o.e;
		f = eval2_ASL(o, ew);
		goto done;
		}
	p = asl->P.cps + i;
	if (p->pi.b) {
		f = copeval(ew, p);
		if (p->g)
			f += cogeval(ew, p);
		}
	else if (p->g)
		f = cogeval(ew, p);
	else if ((o = con_de[i].o.e))
		f = ew->w[o[1]];
	else
		f = 0.;
 done:
	ew->ncxval[i] = ew->nxval;
	ew->err_jmpw = 0;
	return f;
	}

 static real
Conivalp(EvalWorkspace *ew, int i, real *X, fint *nerror)
{
	ASL_pfgh *asl;
	cgrad *gr, **gr0;
	int j1, kv, *vmi;
	real f, *vscale;

	asl = (ASL_pfgh*)ew->asl;
	++ew->stats.conival;
	if (i < asl->i.n_con0)
		f = cpval(ew, i, X, nerror);
	else
		f = 0.;
	kv = 0;
	vmi = 0;
	if ((vscale = asl->i.vscale))
		kv = 2;
	if (asl->i.vmap) {
		vmi = get_vminv_ASL((ASL*)asl);
		++kv;
		}
	if (!(gr0 = asl->i.Cgrad0))
		asl->i.Cgrad0 = gr0 = asl->i.Cgrad_;
	gr = gr0[i];
	switch(kv) {
	  case 3:
		for(; gr; gr = gr->next) {
			j1 = vmi[gr->varno];
			f += X[j1] * vscale[j1] * gr->coef;
			}
		break;
	  case 2:
		for(; gr; gr = gr->next) {
			j1 = gr->varno;
			f += X[j1] * vscale[j1] * gr->coef;
			}
		break;
	  case 1:
		for(; gr; gr = gr->next)
			f += X[vmi[gr->varno]] * gr->coef;
		break;
	  case 0:
		for(; gr; gr = gr->next)
			f += X[gr->varno] * gr->coef;
	  }
	return f;
	}

 real
conpival_nomap_ew_ASL(EvalWorkspace *ew, int i, real *X, fint *nerror)
{
	ASL *a = ew->asl;
	INchk(a, "conpival_nomap", i, a->i.n_con0);
	return  Conivalp(ew, i, X, nerror);
	}

 real
conpival_ew_ASL(EvalWorkspace *ew, int i, real *X, fint *nerror)
{
	ASL *a;
	int *cm;

	a = ew->asl;
	INchk(a, "conpival", i, a->i.n_con_);
	if ((cm = a->i.cmap))
		i = cm[i];
	return  Conivalp(ew, i, X, nerror);
	}

 int
lconpval_ew_ASL(EvalWorkspace *ew, int i, real *X, fint *nerror)
{
	ASL *a;
	real f;

	a = ew->asl;
	INchk(a, "lconpival", i, a->i.n_lcon_);
	++ew->stats.lconval;
	f = cpval(ew, i + a->i.n_con0, X, nerror);
	return f != 0.;
	}

 static void
Congrdp(EvalWorkspace *ew, int i, real *X, real *G, fint *nerror)
{
	ASL_pfgh *asl;
	Jmp_buf err_jmp0;
	cgrad *gr, *gr0;
	fint ne0;
	int i0, j, *vmi;
	ps_func *p;
	real *s, *vscale;

	ne0 = -1;
	if (nerror && (ne0 = *nerror) >= 0) {
		ew->err_jmpw = &err_jmp0;
		i0 = setjmp(err_jmp0.jb);
		if ((*nerror = i0))
			return;
		}
	errno = 0;	/* in case f77 set errno opening files */
	++ew->stats.conigrd;
	asl = (ASL_pfgh*)ew->asl;
	if (!(ew->x0kind & ASL_x_known)) {
		ew->co_index = -(i + 1);
		xp_check_ASL(ew, X);
		}
	if (ew->ncxval[i] != ew->nxval
	 && (!(ew->x0kind & ASL_have_conval)
	     || i < n_conjac[0] || i >= n_conjac[1])) {
		if (!(ew->x0kind & ASL_x_known)) {
			ew->x0kind |= ASL_x_known;
			conpival_ew_ASL(ew,i,X,nerror);
			ew->x0kind &= ~ASL_x_known;
			}
		else
			conpival_ew_ASL(ew,i,X,nerror);
		if (ne0 >= 0 && *nerror)
			return;
		}
	if (ew->Derrs)
		deriv_errchk_ASL(ew, i, 1, 2);
	if (ew->x0kind & ASL_need_comba) {
		funnelset(ew, asl->I.dvfb);
		ew->x0kind &= ~ASL_need_comba;
		}
	if (ew->x0kind & ASL_need_comca) {
		funnelset(ew, asl->I.dvfc);
		ew->x0kind &= ~ASL_need_comca;
		}
	s = ew->derps;
	p = asl->P.cps + i;
	gr0 = asl->i.Cgrad0[i];
	if (p->g)
		psgcomp(ew, p);
	for(gr = gr0; gr; gr = gr->next)
		s[gr->varno] = gr->coef;
	if (p->pi.pdb)
		psderprop(ew, &p->pi);
	if (p->g)
		addgr(ew, p);
	ew->ncxgval[i] = ew->nxval;
	if (!G)
		return;
	vmi = 0;
	if (asl->i.vmap)
		vmi = get_vminv_ASL((ASL*)asl);
	if ((vscale = asl->i.vscale)) {
		gr = gr0;
		if (vmi)
			for(; gr; gr = gr->next) {
				i0 = gr->varno;
				s[i0] *= vscale[vmi[i0]];
				}
		else
			for(; gr; gr = gr->next) {
				i0 = gr->varno;
				s[i0] *= vscale[i0];
				}
		}
	gr = gr0;
	i0 = 0;
	switch(asl->i.congrd_mode) {
	  case 1:
		for(; gr; gr = gr->next)
			G[i0++] = s[gr->varno];
		break;
	  case 2:
		for(; gr; gr = gr->next)
			G[gr->goff] = s[gr->varno];
		break;
	  default:
		if (vmi) {
			for(; gr; gr = gr->next) {
				i = vmi[j = gr->varno];
				while(i0 < i)
					G[i0++] = 0;
				G[i] = s[j];
				i0 = i + 1;
				}
			}
		else
			for(; gr; gr = gr->next) {
				i = gr->varno;
				while(i0 < i)
					G[i0++] = 0;
				G[i] = s[i];
				i0 = i + 1;
				}
		i = n_var;
		while(i0 < i)
			G[i0++] = 0;
	  }
	ew->err_jmpw = 0;
	}

 void
conpgrd_nomap_ew_ASL(EvalWorkspace *ew, int i, real *X, real *G, fint *nerror)
{
	ASL *a;
	static char who[] = "conpgrd_nomap";

	a = ew->asl;
	INchk(a, who, i, a->i.n_con0);
	if (!ew->want_derivs)
		No_derivs_ASL(who);
	Congrdp(ew, i, X, G, nerror);
	}

 void
conpgrd_ew_ASL(EvalWorkspace *ew, int i, real *X, real *G, fint *nerror)
{
	ASL *a;
	int *cm;
	static char who[] = "conpgrd";

	a = ew->asl;
	INchk(a, who, i, a->i.n_con_);
	if (!ew->want_derivs)
		No_derivs_ASL(who);
	if ((cm = a->i.cmap))
		i = cm[i];
	Congrdp(ew, i, X, G, nerror);
	}

 static void
xpsgchk(EvalWorkspace *ew, ps_func *f0, size_t *xv, size_t *xvg, int n, size_t nx,
	real (*ev)(EvalWorkspace *ew, int i, real *X, fint *nerror),
	void (*gv)(EvalWorkspace *ew, int i, real *X, real *G, fint *nerror),
	real *y, int oxk)
{
	int i, i1, i2, needeval;
	ps_func *f;

	i1 = i2 = -1;
	for(i = 0; i < n; i++)
		if (y[i] != 0.) {
			if (i1 < 0)
				i1 = i;
			i2 = i;
			if ((needeval = xv[i] != nx))
				(*ev)(ew, i, ew->Lastx, 0);
			f = f0 + i;
			if (f->g && xvg[i] != nx)
				(*gv)(ew, i, ew->Lastx, 0, 0);
			}
	if (i1 >= 0 && ew->Derrs) {
		ew->x0kind = oxk;
		for(i = i1; i <= i2; i = i1) {
			i1 = i + 1;
			if (y[i]) {
				while(i1 <= i2 && y[i1])
					++i1;
				deriv_errchk_ASL(ew, i, i1-i, 2);
				}
			}
		ew->x0kind |= ASL_x_known;
		}
	}

 void
xpsg_check_ASL(EvalWorkspace *ew, int nobj, real *ow, real *y)
{
	ASL_pfgh *asl;
	int i, nc, no, nz, oxk, tno;
	ps_func *f;
	real t, *x;
	size_t nx, *xv, *xvg;

	asl = (ASL_pfgh*)ew->asl;
	nc = nlc;
	no = nlo;
	if (ew->x0kind & ASL_first_x) {
		if (!(x = X0))
			memset(x = ew->Lastx, 0, n_var*sizeof(real));
		if (y) {
			for(i = 0; i < nc; ++i) {
				if (y[i]) {
					ew->co_index = i;
					goto chk;
					}
				}
			}
		if (ow) {
			for(i = 0; i < no; ++i) {
				if (ow[i]) {
					ew->co_index = -(i+1);
					goto chk;
					}
				}
			}
		if (nobj >= 0 && nobj < no)
			ew->co_index = -(nobj + 1);
		else if (no)
			ew->co_index = -1;
		else
			ew->co_index = 0;
 chk:
		xp_check_ASL(ew, x);
		}
	tno = -1;
	if (!no) {
		if (!nc)
			return;
		ow = 0;
		}
	else if (ow) {
		for(i = 0; i < no; ++i) {
			if ((t = ow[i])) {
				if (t != 1. || tno >= 0) {
					tno = -2;
					break;
					}
				tno = i;
				}
			}
		}
	else if (nobj >= 0 && nobj < nlo)
		tno = nobj;
	if (!(x = ew->oyow))
		ew->oyow = x = ew->oyow0;
	else {
		if (ew->onxval != ew->nxval)
			goto work;
		if (tno == -2) {
			if (memcmp(ow, x, no*sizeof(real)))
				goto work;
			}
		else if (ew->onobj != tno)
			goto work;
		if (nc) {
			if (y) {
				if (memcmp(y, x+no, nc*sizeof(real)))
					goto work;
				}
			else if (ew->nynz)
				goto work;
			}
		return;
		}
 work:
	ew->ihdcur = 0;
	ew->onxval = nx = ew->nxval;
	ew->onobj = tno;
	if (no) {
		if (ow)
			memcpy(x, ow, no*sizeof(real));
		else
			memset(x, 0, no*sizeof(real));
		x += no;
		}
	nz = 0;
	if (nc) {
		if (y) {
			for(i = 0; i < nc; ++i)
				if ((x[i] = y[i]))
					++nz;
			}
		else
			memset(x, 0, nc*sizeof(real));
		}
	else
		y = 0;
	ew->nynz = nz;
	oxk = ew->x0kind;
	ew->x0kind = oxk | ASL_x_known;
	if (y)
		xpsgchk(ew, asl->P.cps, ew->ncxval, ew->ncxgval, nc,
			nx, conpival_ew_ASL, conpgrd_ew_ASL, y, oxk);
	f = asl->P.ops;
	xv = ew->noxval;
	xvg = ew->noxgval;
	if (nobj >= 0 && nobj < n_obj) {
		if (nobj >= no)
			goto done;
		if (ow) {
			ow += nobj;
			if (*ow == 0.)
				goto done;
			}
		if (xv[nobj] != nx)
			objpval_ew_ASL(ew, nobj, ew->Lastx, 0);
		f += nobj;
		if (!f->g)
			goto done;
		if (xvg[nobj] != nx)
			objpgrd_ew_ASL(ew, nobj, ew->Lastx, 0, 0);
		}
	else if (ow && no)
		xpsgchk(ew, f, xv, xvg, no, nx, objpval_ew_ASL, objpgrd_ew_ASL, ow, oxk);
 done:
	ew->x0kind = oxk;
	return;
	}
