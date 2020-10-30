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

#include "nlp.h"
#include "opno.hd"
#include "errchk.h"

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

 static real
eval1(int *o, EvalWorkspace *ew)
{
	TMInfo T, *T1, *T1prev;
	U rv, tv;
	arglist *al;
	char buf[32], *s;
	const char **sa;
	derpblock **db, **db1;
	func_info *fi;
	int i, j, k, n, *o1, **pop, prec, sign, wd, wdf, z;
	jb_st J;
	plterm *p;
	real L, R, *bs, rbuf[128], *rp, t, t0, t1, *w;
	static real Le10;

	w = ew->w;
	wd = ew->wantderiv;
	rv.d = 0.;

 top:
	switch(*o) {
	  case nOPRET:
		rv.d = w[o[1]];
		goto done;
	  case nOPPLUS:
		w[o[1]] = w[o[2]] + w[o[3]];
		o += 4;
		goto top;
	  case nOPMINUS:
		w[o[1]] = w[o[2]] - w[o[3]];
		o += 4;
		goto top;
	  case nOPMULT:
		w[o[1]] = w[o[2]] * w[o[3]];
		o += 4;
		goto top;
	  case nOPDIV0:
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
	  case nOPDIV1:
		L = w[o[2]];
		if (!(R = w[o[3]])
#ifdef WANT_INFNAN
		 && !L
#endif
			)
			zero_div(L, "/");
		w[i = o[1]] = L / R;
		w[i+1] = 1. / R;
		o += 4;
		goto top;
	  case nOPDIV2:
		L = w[o[2]];
		if (!(R = w[o[3]])
#ifdef WANT_INFNAN
		 && !L
#endif
			)
			zero_div(L, "/");
		w[i = o[1]] = rv.d = L / R;
		w[i+1] = -rv.d / R;
		o += 4;
		goto top;
	  case nOPDIV3:
		L = w[o[2]];
		if (!(R = w[o[3]])
#ifdef WANT_INFNAN
		 && !L
#endif
			)
			zero_div(L, "/");
		w[i = o[1]] = rv.d = L / R;
		w[i+2] = -rv.d * (w[i+1] = 1. / R);
		o += 4;
		goto top;
	  case nOPREM0:
		L = w[o[2]];
		R = w[o[3]];
		rv.d = fmod(L,R);
		if (errchk(rv))
			introuble2("fmod",L,R,1);
		w[o[1]] = rv.d;
		o += 4;
		goto top;
	  case nOPREM1:
		L = w[o[2]];
		R = w[o[3]];
		rv.d = fmod(L,R);
		if (errchk(rv))
			introuble2("fmod",L,R,1);
		w[i = o[1]] = rv.d;
		t = -L / R;
		w[i+1] = t >= 0. ? floor(t) : ceil(t);
		o += 4;
		goto top;
	  case nOPPOW0:
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
	  case nOPPOW1:
		rv.d = w[i = o[1]] = mypow(L = w[o[2]], R = w[o[3]]);
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				if (wd)
					w[i] = w[i+1] = w[i+2] = Infinity;
				return w[i] = Infinity;
				}
#endif
			introuble2("pow",L,R,1);
			}
		if (wd) {
			if (L > 0.)
				w[i+1] = R * (rv.d/L);
			else if (L != 0.) {
	 bad:
				introuble2("pow'",L,R,2);
				}
			else {
				if (R > 1.)
					w[i+1] = 0.;
				else if (R == 1.) {
					w[i+1] = 1.;
					}
				else
#ifdef WANT_INFNAN
					w[i+1] = Infinity;
#else
					goto bad;
#endif
				}
			}
		o += 4;
		goto top;
	  case nOPPOW2:
		rv.d = w[i = o[1]] = mypow(L = w[o[2]], R = w[o[3]]);
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				if (wd)
					w[i] = w[i+1] = Infinity;
				return w[i] = Infinity;
				}
#endif
			goto bad;
			}
		if (wd) {
			if (L > 0.)
				w[i+1] = log(L) * rv.d;
			else if (L != 0.)
				goto bad;
			else {
				if (R > 1.)
					w[i+1] = 0.;
				else if (R == 1.)
					w[i+1] = 0.;
				else
#ifdef WANT_INFNAN
					w[i+1] = negInfinity;
#else
					goto bad;
#endif
				}
			}
		o += 4;
		goto top;
	  case nOPPOW3:
		rv.d = w[i = o[1]] = mypow(L = w[o[2]], R = w[o[3]]);
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				if (wd)
					w[i] = w[i+1] = w[i+2] = Infinity;
				return w[i] = Infinity;
				}
#endif
			introuble2("pow",L,R,1);
			}
		if (wd) {
			if (L > 0.) {
				w[i+1] = R * (rv.d/L);
				w[i+2] = log(L) * rv.d;
				}
			else if (L != 0.)
				goto bad;
			else {
				if (R > 1.)
					w[i+1] = w[i+2] = 0.;
				else if (R == 1.) {
					w[i+1] = 1.;
					w[i+2] = 0.;
					}
				else
#ifdef WANT_INFNAN
					{
					w[i+1] = Infinity;
					w[i+2] = negInfinity;
					}
#else
					goto bad;
#endif
				}
			}
		o += 4;
		goto top;
	  case nOPPOW4:
		rv.d = w[i = o[1]] = mypow(L = w[o[2]], R = w[o[3]]);
		/* R = constant integer */
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				if (wd)
					w[i] = w[i+1] = w[i+2] = Infinity;
				return w[i] = Infinity;
				}
#endif
			introuble2("pow",L,R,1);
			}
		if (wd) {
			if (L != 0.)
				w[i+1] = R * (rv.d/L);
			else {
				if (R > 1.)
					w[i+1] = 0.;
				else if (R == 1.) {
					w[i+1] = 1.;
					}
				else
#ifdef WANT_INFNAN
					w[i+1] = Infinity;
#else
					goto bad;
#endif
				}
			}
		o += 4;
		goto top;
	  case nOPLESS0:
		rv.d = w[o[2]] - w[o[3]];
		if (rv.d < 0.)
			rv.d = 0.;
		w[o[1]] = rv.d;
		o += 4;
		goto top;
	  Opalign(OPLESSalign)
	  case nOPLESS1:
		i = o[1];
		db = (derpblock**)(o + 4);
		db1 = (derpblock**)&w[i+1];
		if ((rv.d = w[o[2]] - w[o[3]]) >= 0.)
			*db1 = db[0];
		else {
			rv.d = 0.;
			*db1 = db[1];
			}
		w[i] = rv.d;
		o = (int*)&db[2];
		goto top;
	  case nOPMINLIST0:
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
	  Opalign(OPMINLIST1align)
	  case nOPMINLIST1:
		n = o[2];
		o1 = o + 3;
		rv.d = w[*o1];
		i = j = 0;
		while(++i < n) {
			t = w[o1[i]];
			if (rv.d > t) {
				rv.d = t;
				j = i;
				}
			}
		w[i = o[1]] = rv.d;
		o1 += n;
		*(derpblock**)&w[i+1] = ((derpblock**)o1)[j];
		o = (int*)&((derpblock**)o1)[n];
		goto top;
	  case nOPMAXLIST0:
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
	  Opalign(OPMAXLIST1align)
	  case nOPMAXLIST1:
		n = o[2];
		o1 = o + 3;
		rv.d = w[*o1];
		i = j = 0;
		while(++i < n) {
			t = w[o1[i]];
			if (rv.d < t) {
				rv.d = t;
				j = i;
				}
			}
		w[i = o[1]] = rv.d;
		o1 += n;
		*(derpblock**)&w[i+1] = ((derpblock**)o1)[j];
		o = (int*)&((derpblock**)o1)[n];
		goto top;
	  case nFLOOR:
		w[o[1]] = floor(w[o[2]]);
		o += 3;
		goto top;
	  case nCEIL:
		w[o[1]] = ceil(w[o[2]]);
		o += 3;
		goto top;
	  case nOPABS0:
		if ((rv.d = w[o[2]]) < 0.)
			rv.d = -rv.d;
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case nOPABS1:
		t = 1.;
		if ((rv.d = w[o[2]]) < 0.) {
			rv.d = -rv.d;
			t = -1.;
			}
		w[i = o[1]] = rv.d;
		w[i+1] = t;
		o += 3;
		goto top;
	  case nOPUMINUS:
		w[o[1]] = -w[o[2]];
		o += 3;
		goto top;
	  Opalign(OPORalign)
	  case nOPOR:
		pop = (int**)&o[3];
		if (w[o[2]] != 0.) {
			w[o[1]] = 1.;
			o = *pop;
			goto top;
			}
		o = (int*)&pop[1];
		goto top;
	  Opalign(OPANDalign)
	  case nOPAND:
		pop = (int**)&o[3];
		if (w[o[2]] == 0.) {
			w[o[1]] = 0.;
			o = *pop;
			goto top;
			}
		o = (int*)&pop[1];
		goto top;
	  case nOPLT:
	  case nOPNOTATMOST:
		w[o[1]] = w[o[2]] < w[o[3]] ? 1. : 0.;
		o += 4;
		goto top;
	  case nOPLE:
	  case nOPATLEAST:
		w[o[1]] = w[o[2]] <= w[o[3]] ? 1. : 0.;
		o += 4;
		goto top;
	  case nOPEQ:
	  case nOPEXACTLY:
		w[o[1]] = w[o[2]] == w[o[3]] ? 1. : 0.;
		o += 4;
		goto top;
	  case nOPGE:
	  case nOPATMOST:
		w[o[1]] = w[o[2]] >= w[o[3]] ? 1. : 0.;
		o += 4;
		goto top;
	  case nOPGT:
	  case nOPNOTATLEAST:
		w[o[1]] = w[o[2]] > w[o[3]] ? 1. : 0.;
		o += 4;
		goto top;
	  case nOPNE:
	  case nOPNOTEXACTLY:
		w[o[1]] = w[o[2]] != w[o[3]] ? 1. : 0.;
		o += 4;
		goto top;
	  case nOPNOT:
		w[o[1]] = w[o[2]] == 0. ? 1. : 0.;
		o += 3;
		goto top;
#ifdef X64_bit_pointers
	  case OPIMPELSE_align:
	  case OPIFnl0align:
		++o;
		/* no break */
#endif
	  case nOPIMPELSE:
	  case nOPIFnl0:	/* and OPIFSYM */
		pop = (int**)&o[5];
		if (w[o[2]] == 0.)
			++pop;
		o = *pop;
		goto top;
	  Opalign(OPIFnl1align)
	  case nOPIFnl1:
		pop = (int**)&o[6];
		if (w[o[2]] == 0.)
			++pop;
		*(derpblock**)&w[o[5]] = *(derpblock**)&pop[5];
		o = *pop;
		goto top;
	  case nOP_tanh0:
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
	  case nOP_tanh1:
		i = o[1];
		if ((t = w[o[2]]) >= 175.) {
			rv.d = 1.;
			w[i+1] = 0.;
			}
		else if (t <= -175.) {
			rv.d = -1.;
			w[i+1] = 0.;
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
					w[i+1] = t1*t1;
					}
				}
			}
		w[i] = rv.d;
		o += 3;
		goto top;
	  case nOP_tan0:
		rv.d = tan(w[o[2]]);
		if (errchk(rv))
			introuble("tan",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case nOP_tan1:
		rv.d = tan(t = w[o[2]]);
		if (errchk(rv))
			introuble("tan",t,1);
		i = o[1];
		if (wd) {
			tv.d = cos(t);
			if (errchk(tv) || !tv.d)
				introuble("tan'", t, 2);
			else {
				t1 = 1. / tv.d;
				w[i+1] = t1*t1;
				}
			}
		w[i] = rv.d;
		o += 3;
		goto top;
	  case nOP_sqrt0:
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
	  case nOP_sqrt1:
		t = w[o[2]];
		if (t < 0.)
			goto badsqrt;
		i = o[1];
		rv.d = sqrt(t);
		if (errchk(rv))
			goto badsqrt;
		if (wd) {
			if (rv.d <= 0.)
				introuble("sqrt'",t,2);
			else
				w[i+1] = 0.5 / rv.d;
			}
		w[i] = rv.d;
		o += 3;
		goto top;
	  case nOP_sinh0:
		rv.d = sinh(w[o[2]]);
		if (errchk(rv))
			introuble("sinh",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case nOP_sinh1:
		rv.d = sinh(t = w[o[2]]);
		if (errchk(rv))
			introuble("sinh",t,1);
		i = o[1];
		if (wd) {
			tv.d = cosh(t);
			if (errchk(tv))
				introuble("sinh'", t, 2);
			w[i+1] = tv.d;
			}
		w[i] = rv.d;
		o += 3;
		goto top;
	  case nOP_sin0:
		rv.d = sin(w[o[2]]);
		if (errchk(rv))
			introuble("sin",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case nOP_sin1:
		rv.d = sin(t = w[o[2]]);
		if (errchk(rv))
			introuble("sin",t,1);
		i = o[1];
		if (wd) {
			tv.d = cos(t);
			if (errchk(tv))
				introuble("sin'", t, 2);
			w[i+1] = tv.d;
			}
		w[i] = rv.d;
		o += 3;
		goto top;
	  case nOP_log100:
		rv.d = log10(w[o[2]]);
		if (errchk(rv))
			introuble("log10",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case nOP_log101:
		rv.d = log10(t = w[o[2]]);
		if (errchk(rv))
			introuble("log10",t,1);
		i = o[1];
		if (wd) {
			if (!Le10)
				Le10 = 1. / log(10.);
			w[i+1] = Le10 / t;
			}
		w[i] = rv.d;
		o += 3;
		goto top;
	  case nOP_log0:
		rv.d = log(w[o[2]]);
		if (errchk(rv))
			introuble("log",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case nOP_log1:
		rv.d = log(t = w[o[2]]);
		if (errchk(rv))
			introuble("log",t,1);
		i = o[1];
		if (wd)
			w[i+1] = 1. / t;
		w[i] = rv.d;
		o += 3;
		goto top;
	  case nOP_exp:
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
	  case nOP_cosh0:
		rv.d = cosh(w[o[2]]);
		if (errchk(rv))
			introuble("cosh",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case nOP_cosh1:
		rv.d = cosh(t = w[o[2]]);
		if (errchk(rv))
			introuble("cosh",t,1);
		i = o[1];
		if (wd) {
			tv.d = sinh(t);
			if (errchk(tv))
				introuble("cosh'", t, 2);
			w[i+1] = tv.d;
			}
		w[i] = rv.d;
		o += 3;
		goto top;
	  case nOP_cos0:
		rv.d = cos(w[o[2]]);
		if (errchk(rv))
			introuble("cos",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case nOP_cos1:
		rv.d = cos(t = w[o[2]]);
		if (errchk(rv))
			introuble("cos",t,1);
		i = o[1];
		if (wd) {
			tv.d = sin(t);
			if (errchk(tv))
				introuble("cos'", t, 2);
			w[i+1] = -tv.d;
			}
		w[i] = rv.d;
		o += 3;
		goto top;
	  case nOP_atanh0:
		t = w[o[2]];
		if (t <= -1. || t >= 1.) {
 bad_atan:
			errno_set(EDOM);
			rv.d = 0.;
			introuble("atanh",t,1);
			}
		rv.d = 0.5*log((1. + t) / (1. - t));
		if (errchk(rv))
			goto bad_atan;
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case nOP_atanh1:
		t = w[o[2]];
		if (t <= -1. || t >= 1.)
			goto bad_atan;
		rv.d = 0.5*log((1. + t) / (1. - t));
		if (errchk(rv))
			goto bad_atan;
		w[i = o[1]] = rv.d;
		if (wd)
			w[i+1] = 1. / (1. - t*t);
		o += 3;
		goto top;
	  case nOP_atan20:
		rv.d = atan2(L = w[o[2]], R = w[o[3]]);
		if (errchk(rv))
			introuble2("atan2",L,R,1);
		w[o[1]] = rv.d;
		o += 4;
		goto top;
	  case nOP_atan21:
		rv.d = atan2(L = w[o[2]], R = w[o[3]]);
		if (errchk(rv))
			introuble2("atan2",L,R,1);
		w[i = o[1]] = rv.d;
		if (wd) {
			if ((t = L) < 0.)
				t = -t;
			if ((t1 = R) < 0.)
				t1 = -t1;
			if (t1 >= t) {
				t = L / R;
				t1 = 1. / (1. + t*t);
				w[i+1] = t1 /= R;
				}
			else {
				t = R / L;
				t1 = -1. / (1. + t*t);
				w[i+1] = -t*t1/L;
				}
			}
		o += 4;
		goto top;
	  case nOP_atan22:
		rv.d = atan2(L = w[o[2]], R = w[o[3]]);
		if (errchk(rv))
			introuble2("atan2",L,R,1);
		w[i = o[1]] = rv.d;
		if (wd) {
			if ((t = L) < 0.)
				t = -t;
			if ((t1 = R) < 0.)
				t1 = -t1;
			if (t1 >= t) {
				t = L / R;
				t1 = 1. / (1. + t*t);
				w[i+1] = -t*t1/R;
				}
			else {
				t = R / L;
				t1 = -1. / (1. + t*t);
				w[i+1] = t1 / L;
				}
			}
		o += 4;
		goto top;
	  case nOP_atan23:
		rv.d = atan2(L = w[o[2]], R = w[o[3]]);
		if (errchk(rv))
			introuble2("atan2",L,R,1);
		w[i = o[1]] = rv.d;
		if (wd) {
			if ((t = L) < 0.)
				t = -t;
			if ((t1 = R) < 0.)
				t1 = -t1;
			if (t1 >= t) {
				t = L / R;
				t1 = 1. / (1. + t*t);
				w[i+1] = t1 /= R;
				w[i+2] = -t*t1;
				}
			else {
				t = R / L;
				t1 = -1. / (1. + t*t);
				w[i+2] = t1 /= L;
				w[i+1] = -t*t1;
				}
			}
		o += 4;
		goto top;
	  case nOP_atan0:
		rv.d = atan(w[o[2]]);
		if (errchk(rv))
			introuble("atan",w[o[2]],1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case nOP_atan1:
		rv.d = atan(t = w[o[2]]);
		if (errchk(rv))
			introuble("atan",t,1);
		i = o[1];
		if (wd)
			w[i+1] = 1. / (1. + t*t);
		w[i] = rv.d;
		o += 3;
		goto top;
	  case nOP_asinh0:
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
	  case nOP_asinh1:
		t = t0 = w[o[2]];
		if ((sign = t < 0.))
			t = -t;
		rv.d = log(t + (t1 = sqrt(t*t + 1.)));
		if (errchk(rv))
			introuble("asinh",t0,1);
		if (sign)
			rv.d = -rv.d;
		w[i = o[1]] = rv.d;
		if (wd)
			w[i+1] = 1. / t1;
		o += 3;
		goto top;
	  case nOP_asin0:
		rv.d = asin(t = w[o[2]]);
		if (errchk(rv))
			introuble("asin",t,1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case nOP_asin1:
		rv.d = asin(t = w[o[2]]);
		if (errchk(rv))
			introuble("asin",t,1);
		w[i = o[1]] = rv.d;
		if (wd) {
			if ((t1 = 1. - t*t) <= 0.)
				introuble("asin'",t,2);
			w[i+1] = 1. / sqrt(t1);
			}
		o += 3;
		goto top;
	  case nOP_acosh0:
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
	  case nOP_acosh1:
		if ((t = w[o[2]]) < 1.)
			goto bad_acosh;
		rv.d = log(t + (t1 = sqrt(t*t - 1.)));
		if (errchk(rv))
			goto bad_acosh;
		w[i = o[1]] = rv.d;
		if (wd) {
			if (t1 <= 0.)
				introuble("acosh'",t,2);
			w[i+1] = 1. / t1;
		}
		o += 3;
		goto top;
	  case nOP_acos0:
		rv.d = acos(t = w[o[2]]);
		if (errchk(rv))
			introuble("acos",t,1);
		w[o[1]] = rv.d;
		o += 3;
		goto top;
	  case nOP_acos1:
		rv.d = acos(t = w[o[2]]);
		if (errchk(rv))
			introuble("acos",t,1);
		w[i = o[1]] = rv.d;
		if (wd) {
			if ((t1 = 1. - t*t) <= 0.)
				introuble("acos'",t,2);
			w[i+1] = -1. / sqrt(t1);
			}
		o += 3;
		goto top;
	  case nOPSUMLIST:
		n = o[2];
		j = 3;
		rv.d = 0.;
		for(k = j + n; j < k; ++j)
			rv.d += w[o[j]];
		w[o[1]] = rv.d;
		o += k;
		goto top;
	  case nOPintDIV:
		L = w[o[2]];
		if (!(R = w[o[3]]))
			zero_div(L, " div ");
		w[o[1]] = (L /= R) >= 0. ? floor(L) : ceil(L);
		o += 4;
		goto top;
	  case nOPprecision:
		R = w[o[3]];
		g_fmtp(buf, w[o[2]], (int)R);
		w[o[1]] = strtod(buf, (char**)0);
		o += 4;
		goto top;
	  case nOPround:
		w[o[1]] = Round(w[o[2]], w[o[3]]);
		o += 4;
		goto top;
	  case nOPtrunc:
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
	  case nOPCOUNT:
		n = o[2];
		j = 3;
		rv.d = 0.;
		for(k = j + n; j < k; ++j)
			if (w[o[j]] != 0.)
				++rv.d;
		w[o[1]] = rv.d;
		o += k;
		goto top;
	  case nOPNUMBEROF:
	  case nOPNUMBEROFs:	/* TEMPORARY */
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
	  Opalign(OPPLTERM0align)
	  case nOPPLTERM0:
		i = o[1];
		R = w[o[2]];
		o += 3;
		p = (plterm*)o;
		n = p->n;
		z = p->z;
		bs = p->bs;
		o = (int*)&bs[2*n-1];
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
		for(t = bs[0]*bs[-1]; --z > 0 && R < bs[-3]; bs -= 2) {
			t += bs[-2]*(bs[-3] - bs[-1]);
			}
		w[i] = t + (R - bs[-1])*bs[-2];
		goto top;
	  Opalign(OPPLTERM1align)
	  case nOPPLTERM1:
		i = o[1];
		R = w[o[2]];
		o += 3;
		p = (plterm*)o;
		n = p->n;
		z = p->z;
		bs = p->bs;
		o = (int*)&bs[2*n-1];
		bs += z;
		z >>= 1;
		if (R >= 0) {
			n -= z;
			if (n <= 1 || R <= bs[1]) {
				w[i] = R*(w[i+1] = *bs);
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
		for(t = bs[0]*bs[-1]; --z > 0 && R < bs[-3]; bs -= 2) {
			t += bs[-2]*(bs[-3] - bs[-1]);
			}
		w[i] = t + (R - bs[-1])*(w[i+1] = bs[-2]);
		goto top;
	  case nOPANDLIST:
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
	  case nOPORLIST:
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
	  case nOP_IFF:
		i = w[o[2]] != 0.;
		j = w[o[3]] != 0.;
		w[o[1]] = i == j ? 1. : 0.;
		o += 4;
		goto top;
	  case nOPALLDIFF:
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
	  case nOPSOMESAME:
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
	  case OP2POW0:
		L = w[o[2]];
		w[o[1]] = L*L;
		o += 3;
		goto top;
	  case OP2POW1:
		L = w[o[2]];
		w[i = o[1]] = L*L;
		w[i+1] = L + L;
		o += 3;
		goto top;
#ifdef X64_bit_pointers
	  case OPCPOWalign:
		R = w[o[1]];
		rp = (real*)&o[2];
		goto more_CPOW;
#endif
	  case nOPCPOW:
		rp = (real*)&o[1];
		R = w[o[5]];
 alignarg(more_CPOW:)
		L = rp[0];
		rv.d = mypow(L, R);
		t = 0.;
		if (errchk(rv)) {
#ifdef WANT_INFNAN
			if (!L && R < 0.) {
				errno_set(0);
				rv.d = Infinity;
				goto cpow_store;
				}
#endif
			introuble2("pow",L,R,1);
			}
		if (wd) {
			t = rp[1];
			if (L > 0.)
				t *= rv.d;
			else if (L < 0.) {
	#if defined(WANT_INFNAN) && defined(QNaN0)
				U u;
				u.ui[0] = QNaN0;
				u.ui[1] = QNaN1;
				t = u.r;
	#else
				introuble2("pow'",L,R,2);
	#endif
				}
			}
#ifdef WANT_INFNAN
 cpow_store:
#endif
		w[i = o[6]] = rv.d;
		w[i+1] = t;
		o += 7;
		goto top;

	  case OPFUNCALL0:
		wdf = 0;
		goto more_func;

	  case OPFUNCALL1:
		wdf = wd;
 more_func:
		k = o[1];
		al = ew->al + o[2];
		fi = (func_info*)al->f;
		o += 3;
		if ((n = al->nr)) {
			rp = al->ra;
			for(i = 0; i < n; ++i)
				*rp++ = w[*o++];
			al->derivs = 0;
			if (wdf) {
				al->derivs = rp = w + k + 1;
				for(i = 0; i < n; ++i)
					rp[i] = 0.;
				}
			}
		if ((n = al->n - n)) {
			sa = al->sa;
			for(i = 0; i < n; ++i)
				*sa++ = *(char**)&w[*o++];
			}
		al->Errmsg = 0;
		T.u.prev = 0;
		al->TMI = &T;
		w[k] = (*fi->funcp)(al);
		errno_set(0);
		if ((s = al->Errmsg))
			fintrouble_ASL(ew, fi, s, &T);
		for(T1 = T.u.prev; T1; T1 = T1prev) {
			T1prev = T1->u.prev;
			free(T1);
			}
		goto top;

	  case nOPHOL:
		*(char**)&w[o[1]] = (char*)(o+3);
		o += o[2];
		goto top;

	  /*case nOPVARVAL: accessed directly */
	  case OPCOPY:
		w[o[1]] = w[o[2]];
		o += 3;
		goto top;

	  case OPCOPYSYM:
		*(char**)&w[o[1]] = *(char**)&w[o[2]];
		o += 3;
		goto top;

	  case OPGOTO:
	  case OPNEXTBLK:
		o = *(int**)(o+1);	/* for chaining blocks of code */
		goto top;

#ifdef X64_bit_pointers
	  case OPGOTOalign:
	  case OPNEXTBLKalign:
		o = *(int**)(o+2);
		goto top;
#endif

	  default:
		fprintf(Stderr, "\nUnexpected opno %d in eval1_ASL()\n", *o);
		fflush(Stderr);
		exit(1);
	  }
 done:
	return rv.d;
	}

 static void
comeval(EvalWorkspace *ew, int i, int ie)
{
	ASL_fg *asl;
	cexp *ce, *cx0;
	int *o;
	lincoef *lc, *lce;
	linpart *lp;
	real *dv, t, *w;

	asl = (ASL_fg*)ew->asl;
	cx0 = cexps;
	w = ew->w;
	dv = (real*)ew->dv;
	do {
		ce = cx0 + i;
		o = ce->o.e;
		if (*o == nOPRET)
			t = w[o[1]];
		else {
			ew->cv_index = i + 1;
			t = eval1(o, ew);
			}
		if ((lp = ce->lp)) {
			lc = lp->lc;
			lce = lc + lp->n;
			do t += lc->coef * w[lc->varno];
			   while(++lc < lce);
			}
		dv[i] = t;
		} while(++i < ie);
	ew->cv_index = 0;
	}

 static void
com1eval(EvalWorkspace *ew, int i, int n)
{
	ASL_fg *asl;
	cexp1 *ce, *cee;
	int *o;
	lincoef *lc, *lce;
	linpart *lp;
	real *dv, t, *w;

	asl = (ASL_fg*)ew->asl;
	ce = cexps1 + i;
	cee = ce + n;
	w = ew->w;
	dv = (real*)ew->dv1 + i;
	i += ncom0;
	do {
		ew->cv_index = ++i;
		o = ce->o.e;
		if (*o == nOPRET)
			t = w[o[1]];
		else
			t = eval1(o, ew);
		if ((lp = ce->lp)) {
			lc = lp->lc;
			lce = lc + lp->n;
			do t += lc->coef * w[lc->varno];
			   while(++lc < lce);
			}
		*dv++ = t;
		} while(++ce < cee);
	ew->cv_index= 0;
	}

#ifdef DERPDEBUG
/*DEBUG*/ int derpzork, derpzork0, derpzork1;
#endif

 static void
derpropa(derpblock *db, uint a0, real *s, real *w, real f)
{
	derp *d, *de;
	real t;
	size_t n;
	uint a;

	d = db->d0;
	s[d->b] = f;
#ifdef DERPDEBUG
	printf("derpropa call %d:\n", ++derpzork0);
#endif
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
		if (!(db = db->next)) {
			if (!(db = *(derpblock**)&w[n]))
				break;
			}
		d = db->d0;
		}
	}

 static void
dv_funnel(ASL_fg *asl, EvalWorkspace *ew, cexp **dvf)
{
	cexp *ce;
	derp *d, *de;
	derpblock *db, *db1;
	int n, *vr, *vre;
	real *s, *w;
	uint a;

	w = ew->w;
	s = ew->derps;
	a = asl->i.maxvar;
	while((ce = *dvf++)) {
		db = ce->db;
		if ((vr = ce->vref)) {
			n = vr[0];
			vr += 3;
			vre = vr + n;
			do s[*vr] = 0.;
			   while(++vr < vre);
			}
		derpropa(db, a, s, w, 1.);
		db = ce->dbf;
		d = db->d0;
		de = db->de;
		for(;;) {
			while(d < de) {
				w[d->c] = s[d->a];
				++d;
				}
			if (!db->nxt)
				break;
			if ((db1 = db->next))
				db = db1;
			else
				db = *(derpblock**)&w[db->nxt];
			d = db->d0;
			de = db->de;
			}
		}
	}

 static int
x0_check1(EvalWorkspace *ew, real *X)
{
	ASL_fg *asl;
	int i, nv, *vm;
	real *vscale, *w;

	asl = (ASL_fg*)ew->asl;
	if (x0len == 0) {
		ew->x0kind = 0;
		return 0;
		}
	if (ew->x0kind == ASL_first_x || memcmp(ew->Lastx, X, x0len)) {
		++ew->stats.newx;
		if (ew->Derrs)
			deriv_errclear_ASL(ew);
		memcpy(ew->Lastx, X, x0len);
		++ew->nxval;
		nv = n_var;
		vscale = asl->i.vscale;
		w = ew->w;
		i = 0;
		if ((vm = asl->i.vmap)) {
			if (vscale)
				for(; i < nv; ++i)
					w[vm[i]] = vscale[i] * X[i];
			else
				for(; i < nv; ++i)
					w[vm[i]] = X[i];
			}
		else if (vscale)
				for(; i < nv; ++i)
					w[i] = vscale[i] * X[i];
		/* else Lastx == w */
		ew->x0kind = asl->i.x0kindinit;
		if (comb)
			comeval(ew, 0, comb);
		return 1;
		}
	++ew->stats.oldx;
	return 0;
	}

 int
x1known_ew_ASL(EvalWorkspace *ew, real *X, fint *nerror)
{
	ASL_fg *asl;
	Jmp_buf err_jmp0;
	int ij, rc;

	asl = (ASL_fg*)ew->asl;
	rc = 1;
	if (ew->x0kind & ASL_xknown_ignore)
		goto ret;
	if (nerror && *nerror >= 0) {
		ew->err_jmpw = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	ew->co_index = nlo ? -1 : 0;
	rc = x0_check1(ew, X);
	ew->x0kind |= ASL_x_known;
 done:
	ew->err_jmpw = 0;
 ret:
	return rc;
	}

 real
obj1val_ew_ASL(EvalWorkspace *ew, int nobj, real *X, fint *nerror)
{
	ASL_fg *asl;
	Jmp_buf err_jmp0;
	cde *d;
	int *c1, ij, j1, kv, nc1, *o, *vmi;
	ograd *gr;
	real f, *vscale;

	asl = (ASL_fg*)ew->asl;
	if (nerror && *nerror >= 0) {
		ew->err_jmpw = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij)) {
			f = 0.;
			goto done;
			}
		}
	errno = 0;	/* in case f77 set errno opening files */
	++ew->stats.objval;
	ew->co_index = -(nobj + 1);
	if (!(ew->x0kind & ASL_x_known))
		x0_check1(ew,X);
	if (ew->x0kind & ASL_need_objcom) {
		if (combc < ncom0)
			comeval(ew, combc, ncom0);
		ew->x0kind &= ~ASL_need_objcom;
		}
	if ((c1 = asl->i.o_cexp1st_) && (nc1 = c1[nobj+1] - c1[nobj]))
		com1eval(ew, c1[nobj], nc1);
	d = obj_de + nobj;
	o = d->o.e;
	if (*o == nOPRET)
		f = ew->w[o[1]];
	else
		f = eval1(o, ew);
	ew->noxval[nobj] = ew->nxval;
	kv = 0;
	vmi = 0;
	if ((vscale = asl->i.vscale))
		kv = 2;
	if (asl->i.vmap) {
		vmi = get_vminv_ASL((ASL*)asl);
		++kv;
		}
	gr = Ograd[nobj];
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
 done:
	ew->err_jmpw = 0;
	return f;
	}

 real
con1ival_ew_ASL(EvalWorkspace *ew, int ncon, real *X, fint *nerror)
{
	ASL_fg *asl;
	Jmp_buf err_jmp0;
	cde *d;
	cgrad *gr, **gr0;
	int *c1, ij, j1, kv, nc1, *o, *vmi;
	real *cscale, f, *vscale;

	asl = (ASL_fg*)ew->asl;
	if (nerror && *nerror >= 0) {
		ew->err_jmpw = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij)) {
			f = 0.;
			goto done;
			}
		}
	errno = 0;	/* in case f77 set errno opening files */
	++ew->stats.conival;
	ew->co_index = ncon;
	if (!(ew->x0kind & ASL_x_known))
		x0_check1(ew,X);
	if (ew->x0kind & ASL_need_concom) {
		if (comb < combc)
			comeval(ew, comb, combc);
		ew->x0kind &= ~ASL_need_concom;
		}
	if ((c1 = asl->i.c_cexp1st_) && (nc1 = c1[ncon+1] - c1[ncon]))
		com1eval(ew, c1[ncon], nc1);
	d = con_de + ncon;
	o = d->o.e;
	if (*o == nOPRET)
		f = ew->w[o[1]];
	else
		f = eval1(o, ew);
	ew->ncxval[ncon] = ew->nxval;
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
	gr = gr0[ncon];
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
	if ((cscale = asl->i.cscale))
		f *= cscale[ncon];
 done:
	ew->err_jmpw = 0;
	return f;
	}

 void
con1val_ew_ASL(EvalWorkspace *ew, real *X, real *F, fint *nerror)
{
	ASL_fg *asl;
	Jmp_buf err_jmp0;
	cde *d, *d0;
	cgrad *gr, **gr0;
	int *c1, *cm, i, j, j1, k, kv, nc1, *o, *vmi;
	real *cscale, f, *vscale, *w;

	asl = (ASL_fg *)ew->asl;
	if (nerror && *nerror >= 0) {
		ew->err_jmpw = &err_jmp0;
		i = setjmp(err_jmp0.jb);
		if ((*nerror = i))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	++ew->stats.conval;
	j = n_conjac[0];
	cm = asl->i.cmap;
	if (!(ew->x0kind & ASL_x_known)) {
		ew->co_index = cm ? cm[j] : j;
		x0_check1(ew,X);
		}
	if (ew->x0kind & ASL_need_concom) {
		if (comb < combc)
			comeval(ew, comb, combc);
		ew->x0kind &= ~ASL_need_concom;
		}
	d0 = con_de;
	k = n_conjac[1];
	cscale = asl->i.cscale;
	c1 = asl->i.c_cexp1st_;
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
	w = ew->w;
	for(; j < k; ++j) {
		i = j;
		if (cm)
			i = cm[j];
		ew->co_index = i;
		d = d0 + i;
		if (c1 && (nc1 = c1[i+1] - c1[i]))
			com1eval(ew, c1[i], nc1);
		o = d->o.e;
		if (*o == nOPRET)
			f = w[o[1]];
		else
			f = eval1(o, ew);
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
 done:
	ew->err_jmpw = 0;
	}

 int
lcon1val_ew_ASL(EvalWorkspace *ew, int i, real *X, fint *nerror)
{
	ASL_fg *asl;
	Jmp_buf err_jmp0;
	cde *d;
	int *c1, ij, nc1, *o;
	real f;

	asl = (ASL_fg*)ew->asl;
	if (nerror && *nerror >= 0) {
		ew->err_jmpw = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij)) {
			f = 0.;
			goto done;
			}
		}
	errno = 0;	/* in case f77 set errno opening files */
	++ew->stats.lconval;
	ew->co_index = i + n_con;
	if (!(ew->x0kind & ASL_x_known))
		x0_check1(ew,X);
	if (ew->x0kind & ASL_need_concom) {
		if (comb < combc)
			comeval(ew, comb, combc);
		ew->x0kind &= ~ASL_need_concom;
		}
	if ((c1 = l_cexp1st) && (nc1 = c1[i+1] - c1[i]))
		com1eval(ew, c1[i], nc1);
	d = lcon_de + i;
	o = d->o.e;
	if (*o == nOPRET)
		f = ew->w[o[1]];
	else
		f = eval1(o, ew);
	ew->nlxval[i] = ew->nxval;
 done:
	ew->err_jmpw = 0;
	return f != 0.;
	}

 void
obj1grd_ew_ASL(EvalWorkspace *ew, int nobj, real *X, real *G, fint *nerror)
{
	ASL_fg *asl;
	Jmp_buf err_jmp0;
	cde *c;
	derpblock *db;
	fint ne0;
	int *dv, *dvr, i, ij, j, *vmi, *z;
	lincoef *lc, *lce;
	linpart *lp, **plp;
	ograd *gr, *gr1;
	real *s, t, *vscale, *w;
	static char who[] = "obj1grd";

	asl = (ASL_fg*)ew->asl;
	if (!want_derivs)
		No_derivs_ASL(who);
	ne0 = -1;
	if (nerror && (ne0 = *nerror) >= 0) {
		ew->err_jmpw = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	++ew->stats.objgrd;
	if (!(ew->x0kind & ASL_x_known)) {
		ew->co_index = -(nobj + 1);
		x0_check1(ew,X);
		}
	if (!ew->noxval || ew->noxval[nobj] != ew->nxval) {
		if (!(ew->x0kind & ASL_x_known)) {
			ew->x0kind |= ASL_x_known;
			obj1val_ew_ASL(ew, nobj, X, nerror);
			ew->x0kind &= ~ASL_x_known;
			}
		else
			obj1val_ew_ASL(ew, nobj, X, nerror);
		if (ne0 >= 0 && *nerror)
			goto done;
		}
	if (ew->Derrs)
		deriv_errchk_ASL(ew, -(nobj+1), 1, 2);
	if (ew->x0kind & ASL_need_comba) {
		dv_funnel(asl, ew, asl->I.dvfb);
		ew->x0kind &= ~ASL_need_comba;
		}
	if (ew->x0kind & ASL_need_comoa) {
		dv_funnel(asl, ew, asl->I.dvfo);
		ew->x0kind &= ~ASL_need_comoa;
		}
	s = ew->derps;
	w = ew->w;
	c = &obj_de[nobj];
	if (c->afn)
		memset(s + c->af1, 0, c->afn*sizeof(real));
	if ((dvr = c->dvref)) {
		dv = dvr + *dvr;
		do s[*dv] = 0.; while(--dv > dvr);
		}
	for(gr = gr1 = Ograd[nobj]; gr; gr = gr->next)
		s[gr->varno] = gr->coef;
	if ((db = c->db))
		derpropa(db, asl->i.maxvar, s, w, 1.);
	if ((plp = c->c1lp)) {
		while((lp = *plp++)) {
			if ((t = s[lp->a])) {
				lc = lp->lc;
				lce = lc + lp->n;
				do s[lc->varno] += t*lc->coef;
					while(++lc < lce);
				}
			}
		}
	if (zerograds) {	/* sparse gradients */
		z = zerograds[nobj];
		while((i = *z++) >= 0)
			G[i] = 0;
		}
	gr = gr1;
	vmi = 0;
	if (asl->i.vmap)
		vmi = get_vminv_ASL((ASL*)asl);
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
 done:
	ew->err_jmpw = 0;
	}

 void
con1grd_ew_ASL(EvalWorkspace *ew, int nc, real *X, real *G, fint *nerror)
{
	ASL_fg *asl;
	Jmp_buf err_jmp0;
	cde *c;
	cgrad *gr, *gr1;
	derpblock *db;
	int *dv, *dvr, i, i0, ij, j, *vmi;
	lincoef *lc, *lce;
	linpart *lp, **plp;
	real *s, t, *vscale, *w;
	static char who[] = "con1grd";

	asl = (ASL_fg*)ew->asl;
	if (!want_derivs)
		No_derivs_ASL(who);
	if (nerror && *nerror >= 0) {
		ew->err_jmpw = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			return;
		}
	errno = 0;	/* in case f77 set errno opening files */
	++ew->stats.conigrd;
	if (!(ew->x0kind & ASL_x_known)) {
		ew->co_index = nc;
		x0_check1(ew,X);
		}
	if ((!ew->ncxval || ew->ncxval[nc] != ew->nxval)
	 && (!(ew->x0kind & ASL_have_conval)
	     || nc < n_conjac[0] || nc >= n_conjac[1])) {
		if (!(ew->x0kind & ASL_x_known)) {
			ew->x0kind |= ASL_x_known;
			con1ival(ew,nc,X,nerror);
			ew->x0kind &= ~ASL_x_known;
			}
		else
			con1ival(ew,nc,X,nerror);
		if (nerror && *nerror)
			return;
		}
	if (ew->Derrs)
		deriv_errchk_ASL(ew, nc, 1, 2);
	if (ew->x0kind & ASL_need_comba) {
		dv_funnel(asl, ew, asl->I.dvfb);
		ew->x0kind &= ~ASL_need_comba;
		}
	if (ew->x0kind & ASL_need_comca) {
		dv_funnel(asl, ew, asl->I.dvfc);
		ew->x0kind &= ~ASL_need_comca;
		}
	s = ew->derps;
	w = ew->w;
	gr1 = asl->i.Cgrad0[nc];
	c = &con_de[nc];
	if (c->afn)
		memset(s + c->af1, 0, c->afn*sizeof(real));
	if ((dvr = c->dvref)) {
		dv = dvr + *dvr;
		do s[*dv] = 0.; while(--dv > dvr);
		}
	for(gr = gr1; gr; gr = gr->next)
		s[gr->varno] = gr->coef;
	if ((db = c->db))
		derpropa(db, asl->i.maxvar, s, w, 1.);
	if ((plp = c->c1lp)) {
		while((lp = *plp++)) {
			if ((t = s[lp->a])) {
				lc = lp->lc;
				lce = lc + lp->n;
				do s[lc->varno] += t*lc->coef;
					while(++lc < lce);
				}
			}
		}
	vmi = 0;
	if (asl->i.vmap)
		vmi = get_vminv_ASL((ASL*)asl);
	if ((vscale = asl->i.vscale)) {
		if (vmi)
			for(gr = gr1; gr; gr = gr->next) {
				i0 = gr->varno;
				s[i0] *= vscale[vmi[i0]];
				}
		else
			for(gr = gr1; gr; gr = gr->next) {
				i0 = gr->varno;
				s[i0] *= vscale[i0];
				}
		}
	gr = gr1;
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
jac1val_ew_ASL(EvalWorkspace *ew, real *X, real *G, fint *nerror)
{
	ASL_fg *asl;
	Jmp_buf err_jmp0;
	cde *c, *c0;
	cgrad *gr, **gr0, *gr1;
	derpblock *db;
	fint ne0;
	int *cm, *dv, *dvr, i, j, j1, k, maxvar, *vmi;
	lincoef *lc, *lce;
	linpart *lp, **plp;
	real *cscale, *s, t, *vscale, *w;
	static char who[] = "jac1val";

	asl = (ASL_fg*)ew->asl;
	maxvar = asl->i.maxvar;
	if (!want_derivs)
		No_derivs_ASL(who);
	ne0 = -1;
	if (nerror && (ne0 = *nerror) >= 0) {
		ew->err_jmpw = &err_jmp0;
		j = setjmp(err_jmp0.jb);
		if ((*nerror = j))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	++ew->stats.jacval;
	cm = asl->i.cmap;
	j = n_conjac[0];
	ew->co_index = cm ? cm[j] : j;
	if ((!(ew->x0kind & ASL_x_known) && x0_check1(ew,X))
	 || !(ew->x0kind & ASL_have_conval)) {
		if (!(ew->x0kind & ASL_x_known)) {
			ew->x0kind |= ASL_x_known;
			con1val_ew_ASL(ew, X, 0, nerror);
			ew->x0kind &= ~ASL_x_known;
			}
		else
			con1val_ew_ASL(ew, X, 0, nerror);
		if (ne0 >= 0 && *nerror)
			goto done;
		}
	k = n_conjac[1];
	if (ew->Derrs)
		deriv_errchk_ASL(ew, j, k-j, 2);
	if (asl->i.zap_J)
		memset(G, 0, asl->i.zap_J);
	s = ew->derps;
	w = ew->w;
	c0 = con_de;
	cscale = asl->i.cscale;
	vscale = asl->i.vscale;
	vmi = 0;
	if (asl->i.vmap)
		vmi = get_vminv_ASL((ASL*)asl);
	if (ew->x0kind & ASL_need_comba) {
		dv_funnel(asl, ew, asl->I.dvfb);
		ew->x0kind &= ~ASL_need_comba;
		}
	if (ew->x0kind & ASL_need_comca) {
		dv_funnel(asl, ew, asl->I.dvfc);
		ew->x0kind &= ~ASL_need_comca;
		}
	gr0 = asl->i.Cgrad0;
	for(; j < k; ++j) {
		j1 = j;
		if (cm)
			j1 = cm[j];
		c = &c0[j1];
		if (c->afn)
			memset(s + c->af1, 0, c->afn*sizeof(real));
		if ((dvr = c->dvref)) {
			dv = dvr + *dvr;
			do s[*dv] = 0.; while(--dv > dvr);
			}
		for(gr = gr1 = gr0[j1]; gr; gr = gr->next)
			s[gr->varno] = gr->coef;
		if ((db = c->db))
			derpropa(db, maxvar, s, w, 1.);
		if ((plp = c->c1lp)) {
			while((lp = *plp++)) {
				if ((t = s[lp->a])) {
					lc = lp->lc;
					lce = lc + lp->n;
					do s[lc->varno] += t*lc->coef;
						while(++lc < lce);
					}
				}
			}
		if (vscale) {
			gr = gr1;
			if (vmi)
				for(; gr; gr = gr->next) {
					i = gr->varno;
					s[i] *= vscale[vmi[i]];
					}
			else
				for(; gr; gr = gr->next) {
					i = gr->varno;
					s[i] *= vscale[i];
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
 done:
	ew->err_jmpw = 0;
	}
