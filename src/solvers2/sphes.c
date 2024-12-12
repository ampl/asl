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

#include "asl_pfgh.h"
#include "opno2.h"
#include "obj_adj.h"
#ifdef X64_bit_pointers
#define alignarg(x) x
#else
#define alignarg(x) /*nothing*/
#endif
#ifdef NO_PTHREADS
#define pthread_mutex_lock(t) /*nothing*/
#define pthread_mutex_unlock(t) /*nothing*/
#define PTHREADS(x) /*nothing*/
#else
#define PTHREADS(x) x
#endif
#undef exit
extern void wk_init_ASL(real*, int*, real);

 static void
hv_fwd(int *o, real *w, int *oend)	/* for determining sparsity */
{
	Condptrs *cp;
	Eresult *L, *R, *r;
	Minmaxptrs *mmp, *mmpe;
	int *d, *de, k, *o0, *o1, *opg1, **pop;
	plterm **pp;
	real dO, *rp;
	tfinfo **ptfi, *tfi;

	for(;;) {
		o0 = o;
		switch(*o) {

		case OPRET:
			return;

		case OP_GOTO:
		case OP_NEXTBLK:
		case OPGOTOF:
		case OPGOTOF2n:
			o = *(int**)(o+1);
			continue;

		case OPGOTO2:
		case OPGOTOF2:
			pop = (int**)(o+1);
			o = pop[1];
			continue;

/*		case Hv_timesR:	*/
/*		case Hv_binaryR:*/
		case OPDIV01:
		case OPMULT01:
		case nOPPOW01:
		case nOPREM01:
		case OP_atan201:
			r = (Eresult*)(w + o[2]);
			R = (Eresult*)(w + o[4]);
			r->dO = R->dO;
			o += 5;
			break;

#ifdef X64_bit_pointers
		case OPGOTO2align:
		case OPGOTOF2align:
			pop = (int**)(o+2);
			o = pop[1];
			continue;

		case OPGOTOF2nalign:
			o = *(int**)(o+2);
			continue;

		case OPLESS01align:
			o1 = o + 6 + 2*sizeof(derpblock*)/sizeof(int);
			goto more_OPLESS01;
#endif
		case nOPLESS01:
			o1 = o + 5 + 2*sizeof(derpblock*)/sizeof(int);
alignarg(more_OPLESS01:)
			r = (Eresult*)(w + o[2]);
			R = (Eresult*)(w + o[4]);
			r->dO = r->dL * R->dO;
			o = o1;
			break;

/*		case Hv_timesLR:*/
		case OPMULT2:
/*		case Hv_binaryLR:*/
		case n_OPPOW2:
		case OP_atan22:
		case OPDIV2:
		case nOPREM2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			r->dO = L->dO + R->dO;
			o += 5;
			break;

#ifdef X64_bit_pointers
		case OPLESS2align:
			o1 = o + 6 + 2*sizeof(derpblock*)/sizeof(int);
			goto more_OPLESS2;
#endif
		case nOPLESS2:
			o1 = o + 5 + 2*sizeof(derpblock*)/sizeof(int);
alignarg(more_OPLESS2:)
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			r->dO = L->dO + R->dO;
			o = o1;
			break;

/*		case Hv_timesL:	*/
/*		case Hv_unary:	*/
		case OP_2POW1:
		case n_OPABS1:
		case OP_acos1:
		case OP_acosh1:
		case OP_asin1:
		case OP_asinh1:
		case OP_atan1:
		case OP_atanh1:
		case OP_cos1:
		case OP_cosh1:
		case OP_exp1:
		case OP_log101:
		case OP_log1:
		case OP_sin1:
		case OP_sinh1:
		case OP_sqrt1:
		case OP_tan1:
		case OPtanh1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = L->dO;
			o += 4;
			break;

		case OPDIV10:
		case OPMULT10:
		case nOPREM10:
		case nOPPOW10:
		case OP_atan210:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = L->dO;
			o += 5;
			break;

#ifdef X64_bit_pointers
		case OP_GOTOalign:
		case OP_NEXTBLKalign:
		case OPGOTOFalign:
			o = *(int**)(o+2);
			continue;

		case OPLESS10align:
			o1 = o + 6 + 2*sizeof(derpblock*)/sizeof(int);
			goto more_OPLESS10;
#endif
		case nOPLESS10:
			o1 = o + 5 + 2*sizeof(derpblock*)/sizeof(int);
alignarg(more_OPLESS10:)
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = r->dL * L->dO;
			o = o1;
			break;

/*		case Hv_vararg:	*/
#ifdef X64_bit_pointers
		case OPMINLISTalign:
		case OPMAXLISTalign:
			opg1 = o + 5;
			goto more_maxminlist;
#endif
		case OPMINLIST1:
		case OPMAXLIST1:
			opg1 = o + 4;
alignarg(more_maxminlist:)
			r = (Eresult*)(w + o[2]);
			r->dO = 0.;
			mmp = (Minmaxptrs*)&opg1[k = o[3]];
			for(mmpe = mmp + k; mmp < mmpe; ++mmp) {
				if ((k = mmp->d) >= 0) {
					L = (Eresult*)(w + k);
					if ((r->dO = L->dO))
						break;
					}
				}
			o = (int*)mmpe + 2;
			break;

/*		case Hv_if:	*/
#ifdef X64_bit_pointers
		case OPIF11align:
		case OPIF13align:
			cp = (Condptrs*)&o[6];
			goto more_if13;
		case OPIF1align:
			cp = (Condptrs*)&o[6];
			goto more_if;
#endif
		case nOPIF11:
		case nOPIF13:
			cp = (Condptrs*)&o[5];
alignarg(more_if13:)
			r = (Eresult*)(w + o[2]);
			r->dO = 0.;
			if (cp->bder >= 0) {
				o = cp->f;
				if (*o == OPRET) {
					o = cp->b;
					if (*o == OPCOPY1 || *o == OPCOPY1a) {
						r = (Eresult*)(w + o[2]);
						L = (Eresult*)(w + o[3]);
						r->dO = L->dO;
						}
					}
				else
					hv_fwd(o, w, cp->b);
				}
			o = cp[1].f;
			continue;
		case nOPIF1:
			cp = (Condptrs*)&o[5];
alignarg(more_if:)
			r = (Eresult*)(w + o[2]);
			r->dO = 0.;
			if (cp->bder >= 0)
				hv_fwd(cp->f, w, cp->b);
			o = cp[1].f;
			continue;

/*		case Hv_plterm:*/
#ifdef X64_bit_pointers
		case OP_PLTERM1align:
			pp = (plterm**)(o+5);
			goto more_plterm;
#endif
		case n_OPPLTERM1:
			pp = (plterm**)(o+4);
alignarg(more_plterm:)
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = L->dO;
			o = (int*)&pp[1];
			break;

/*		case Hv_sumlist:	*/
		case OPSUMLIST1:
			r = (Eresult*)(w + o[2]);
			o1 = o + 4;
			o = o1 + o[3];
			dO = 0.;
			while(o1 < o) {
				if ((k = *o1++) >= 0) {
					L = (Eresult*)(w + k);
					dO += L->dO;
					}
				}
			r->dO = dO;
			break;

/*		case Hv_func:	*/
#ifdef X64_bit_pointers
		case OP_FUNCALL1align:
			ptfi = (tfinfo**)(o+4);
			goto more_func;
#endif
		case OP_FUNCALL1:
			ptfi = (tfinfo**)(o+3);
alignarg(more_func:)
			tfi = *ptfi;
			r = (Eresult*)(rp = w + o[2]);
			rp += 4;
			dO = 0.;
			for(d = tfi->doff, de = d + 2*tfi->nd; d < de; d += 2) {
				L = (Eresult*)(w + *d);
				dO += L->dO;
				}
			r->dO = dO;
			o = (int*)&ptfi[1] + tfi->n;
			break;

/*		case Hv_negate:	*/
		case OPUMINUS1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = L->dO;
			o += 4;
			break;

/*		case Hv_plusR:	*/
		case OPPLUS01:
			r = (Eresult*)(w + o[2]);
			R = (Eresult*)(w + o[4]);
			r->dO = R->dO;
			o += 5;
			break;

/*		case Hv_plusL:	*/
		case OPMINUS10:
		case OPPLUS10:
		case nOPPOW1i:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = L->dO;
			o += 5;
			break;

/*		case Hv_plusLR:	*/
		case OPPLUS2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			r->dO = L->dO + R->dO;
			o += 5;
			break;

/*		case Hv_minusR:	*/
		case OPMINUS01:
			r = (Eresult*)(w + o[2]);
			R = (Eresult*)(w + o[4]);
			r->dO = R->dO;
			o += 5;
			break;

/*		case Hv_minusLR: */
		case OPMINUS2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			r->dO = L->dO + R->dO;
			o += 5;
			break;

		case OPCOPY1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = L->dO;
			o += 4;
			break;

		case OPCOPY1a:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO += L->dO;
			o += 4;
			break;

#ifdef X64_bit_pointers
		case OPCPOW1align:
			rp = (real*)&o[5];
			goto more_OPCPOW1;
#endif

		case nOPCPOW1:
			rp = (real*)&o[4];
alignarg(more_OPCPOW1:)
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = L->dO;
			o = (int*)&rp[1];
			break;

		case OPVARREF:
			r = (Eresult*)&w[o[2]];
			o += 3;
			break;

#ifdef X64_bit_pointers
		case OPGOTOBalign:
			++o;
#endif
		case OPGOTOB:
			pop = (int**)(o+1);
			o = (int*)&pop[1];
			goto endchk;

		default:/*DEBUG*/
			fprintf(Stderr, "Bad *o = %d in hv_fwd\n", *o);
			exit(1);
			r = 0; /* not reached */
		  }
		r->aO = r->adO = 0.;
 endchk:
		if (o0 == oend)
			return;
		}
	}

 static void
func_back(int *o, real *w, tfinfo **ptfi)
{
	Eresult *r;
	int *d1, *d2, *de, *doff;
	real aO, adO, t;
	tfinfo *tfi;

	tfi = *ptfi;
	r = (Eresult*)(w + o[2]);
	doff = tfi->doff;
	aO = r->aO;
	adO = r->adO;
	de = doff + 2*tfi->nd;
	for(d1 = doff; d1 < de; d1 += 2) {
		r = (Eresult*)(w + *d1);
		r->adO += adO;
		r->aO += aO;
		t = adO*r->dO;
		for(d2 = doff; d2 < de; d2 += 2) {
			r = (Eresult*)(w + *d2);
			r->aO += t;
			}
		}
	}

static void
hv_back1(int *o, real *w, int *oend)
{
	Condptrs *cp;
	Eresult *r, *L, *R;
	Minmaxptrs *mmp, *mmp1;
	int k, *o1, *oe, *opg1;
	real adO, t1, t2;

 top:
	for(;;) {
		switch(*o) {
		case OPRET:
		case OPRETB:
			return;
#ifdef X64_bit_pointers
		case OPGOTOBalign:
			++o;
#endif
		case OPGOTOB:
			o = *(int**)(o + 1);
			goto top;
/*		case Hv_binaryR: */
		case OPDIV01:
		case nOPLESS01:
		case nOPPOW01:
		case nOPREM01:
		case OP_atan201:
			r = (Eresult*)(w + o[2]);
			R = (Eresult*)(w + o[4]);
			R->adO += r->adO;
			R->aO += r->aO  +  r->adO * R->dO;
			break;

/*		case Hv_binaryLR: */
		case n_OPPOW2:
		case OP_atan22:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			adO = r->adO;
			L->adO += adO;
			R->adO += adO;;
			t1 = adO * L->dO;
			t2 = adO * R->dO;
			L->aO  += r->aO + t1 + t2;
			R->aO  += r->aO + t1 + t2;
			break;

		case OPDIV2:
			/* dL2 is used for dR2, */
			/* since the true dL2 vanishes */
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			adO = r->adO;
			L->adO += adO;
			R->adO += adO;
			t1 = adO * L->dO;
			t2 = adO * R->dO;
			L->aO  += r->aO + t2;
			R->aO  += r->aO + t1 + t2;
			break;

		case nOPREM2:
		case nOPLESS2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			adO = r->adO;
			L->aO  += r->aO;
			R->aO  += r->aO;
			L->adO += adO;
			R->adO += adO;
			break;

/*		case Hv_unary: */
		case n_OPABS1:
		case nOPREM10:
		case nOPLESS10:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->adO += r->adO;
			L->aO  += r->aO;
			break;

/*		case Hv_unary: */
	alignarg(case OPCPOW1align:)
		case nOPCPOW1:
		case nOPPOW10:
		case nOPPOW1i:
		case OP_acos1:
		case OP_acosh1:
		case OP_asin1:
		case OP_asinh1:
		case OP_atan1:
		case OP_atan210:
		case OP_atanh1:
		case OP_cos1:
		case OP_cosh1:
		case OP_exp1:
		case OP_log101:
		case OP_log1:
		case OP_sin1:
		case OP_sinh1:
		case OP_sqrt1:
		case OP_tan1:
		case OPtanh1:
		case OP_2POW1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->adO += r->adO;
			L->aO  += r->aO  +  r->adO * L->dO;
			break;

/*		case Hv_vararg:	*/
#ifdef X64_bit_pointers
		case OPMINLISTalign:
		case OPMAXLISTalign:
			opg1 = o + 5;
			goto more_maxminlist;
#endif
		case OPMINLIST1:
		case OPMAXLIST1:
			opg1 = o + 4;
alignarg(more_maxminlist:)
			r = (Eresult*)(w + o[2]);
			if (r->adO || r->aO) {
				mmp = (Minmaxptrs*)&opg1[o[3]];
				mmp1 = mmp + o[3];
				opg1 = (int*)mmp;
				while(mmp1 > mmp) {
					--opg1;
					--mmp1;
					if (mmp1->d >= 0) {
						L = (Eresult*)(w + *opg1);
						if (L->dO) {
							L->aO += r->aO;
							L->adO += r->adO;
							if (mmp1->b)
								hv_back1(mmp1->b, w, mmp1->f);
							}
						}
					}
				}
			break;

/*		case Hv_if:	*/
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
			if (cp[1].bder >= 0)
				hv_back1(cp[1].b, w, cp[1].f);
			if (cp[0].bder >= 0)
				hv_back1(cp[0].b, w, cp[0].f);
			break;

#ifdef X64_bit_pointers
		case OP_PLTERM1align:
#endif
		case n_OPPLTERM1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->aO += r->aO;
			break;

		case OPSUMLIST1:
			r = (Eresult*)(w + o[2]);
			t1 = r->aO;
			t2 = r->adO;
			o1 = o + 4;
			oe = o1 + o[3];
			while(o1 < oe) {
				if ((k = *o1++) >= 0) {
					L = (Eresult*)(w + k);
					L->aO += t1;
					L->adO += t2;
					}
				}
			break;

/*		case Hv_func: */
		case OP_FUNCALL1:
			func_back(o, w, (tfinfo**)(o+3));
			break;

#ifdef X64_bit_pointers
		case OP_FUNCALL1align:
			func_back(o, w, (tfinfo**)(o+4));
			break;
#endif

/*		case Hv_plusL: */
		case OPPLUS10:
		case OPMINUS10:
/*		case Hv_negate: */
		case OPUMINUS1:
			L = (Eresult*)(w + o[3]);
 neg_end:
			r = (Eresult*)(w + o[2]);
			L->aO += r->aO;
			L->adO += r->adO;
			break;

/*		case Hv_plusR: */
		case OPPLUS01:
			L = (Eresult*)(w + o[4]);
			goto neg_end;

/*		case Hv_plusLR: */
		case OPPLUS2:
/*		case Hv_minusLR: */
		case OPMINUS2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			L->aO += t1 = r->aO;
			L->adO += t2 = r->adO;
			R->aO += t1;
			R->adO += t2;
			break;

/*		case Hv_minusR: */
		case OPMINUS01:
			L = (Eresult*)(w + o[4]);
			goto neg_end;

/*		case Hv_timesR: */
		case OPMULT01:
			r = (Eresult*)(w + o[2]);
			R = (Eresult*)(w + o[4]);
			R->adO += r->adO;
			R->aO += r->aO;
			break;

/*		case Hv_timesL: */
		case OPDIV10:
		case OPMULT10:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->adO += r->adO;
			L->aO  += r->aO;
			break;

/*		case Hv_timesLR: */
		case OPMULT2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			adO = r->adO;
			L->aO  += r->aO  +  adO * R->dO;
			R->aO  += r->aO  +  adO * L->dO;
			L->adO += adO;
			R->adO += adO;
			break;

		  case OPCOPY1:
		  case OPCOPY1a:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->adO += r->adO;
			L->aO += r->aO;
			break;

		  case OPVARREF:
			break;

		  default:/*DEBUG*/
			fprintf(Stderr, "Bad *o = %d in hv_back\n", *o);
			exit(1);
		  }
		if (o == oend)
			return;
		o -= o[1];
		}
	}

 static void
hv_back(int *o, real *w, real aO, real adO0, int *oend)
{
	Eresult *r;

	for(;;) {
		switch(*o) {
			case OPRET:
				return;
			case OPGOTOB:
				o = *(int**)(o + 1);
				continue;

#ifdef X64_bit_pointers
			case OPGOTOBalign:
				o = *(int**)(o + 2);
				continue;
#endif
			}
		break;
		}
	r = (Eresult*)(w + o[2]);
	r->aO = aO;
	r->adO = adO0;
	hv_back1(o, w, oend);
	}

 static void
hv_fwd0(EvalWorkspace *ew, cexp *c, Varval *v)	/* for determining sparsity */
{
	Varval *V, *v1;
	int i, *o, *ob;
	linarg *la;
	lincoef *lc, *lce;
	linpart *L;
	real *w, x;

	w = ew->w;
	V = (Varval*)w;
	v->aO = v->adO = 0;
	if ((o = c->o.f)) {
		hv_fwd(o, w, ob = c->o.b);
		v1 = (Varval*)(w + ob[2]);
		x = v1->dO;
		}
	else if ((o = c->o.e) && *o == OPRET && (i = o[1]) >= 0)
		x = w[i+1];
	else
		x = 0.;
	if ((la = c->la))
		x += V[la->u.v].dO;
	else if ((L = c->lp)) {
		lc = L->lc;
		lce = lc + L->n;
		for(lce = lc + L->n; lc < lce; ++lc)
			x += V[lc->varno].dO;
		}
	v->dO = x;
	}

 static void
pshv_prod1(EvalWorkspace *ew, range *r, int nobj, int ow, int y) /* for determining sparsity */
{
	ASL_pfgh *asl;
	Varval *V, *Vc, *v;
	cexp *c;
	int *cei, *cei0, *ceie, i, *o, *o1;
	linarg *la, **lap, **lape;
	lincoef *lc, *lce;
	linpart *L;
	psb_elem *b;
	real *s, *w;

	asl = (ASL_pfgh*)ew->asl;
	V = (Varval*)(w = ew->w);
	Vc = V + asl->i.defvar0;
	s = w + asl->P.dOscratch;
	lap = r->lap;
	lape = lap + r->n;
	while(lap < lape) {
		la = *lap++;
		v = V + la->u.v;
		v->dO = *s++;
		v->adO = v->aO = 0.;
		}
	if ((cei = cei0 = r->cei)) {
		i = *cei0++;
		ceie = (cei = cei0) + i;
		do {
			i = *cei++;
			hv_fwd0(ew, cexps + i, Vc + i);
			}
			while(cei < ceie);
		}
	for(b = r->refs; b; b = b->next) {
		if ((i = b->conno) < 0) {
			i = -2 - i;
			if (!ow && i != nobj)
				continue;
			}
		else if (!y)
			continue;
		if ((o = b->o.f)) {
			hv_fwd(o, w, o1 = b->o.b);
			hv_back(o1, w, 0., 1., o);
			}
		}
	while(cei > cei0) {
		i = *--cei;
		c = cexps + i;
		v = Vc + i;
		if (v->aO && (L = c->lp)) {
			if ((la = c->la))
				V[la->u.v].aO = 1;
			else
				for(lc = L->lc, lce = lc + L->n; lc < lce; ++lc)
					V[lc->varno].aO++;
			}
		if ((o = c->o.b))
			hv_back(o, w, 1., v->adO, c->o.f);
		}
	}

#ifdef __cplusplus
extern "C" {
static int compar(const void*, const void*);
}
#endif

 static int
compar(const void *a, const void *b)
{ return *(int*)a - *(int*)b; }

#undef nzc
#undef asl
#undef del_mblk
#define del_mblk(c) Del_mblk_ASL(a,(Char*)(c))

 static void
new_Hesoprod(EvalWorkspace *ew, int nov, int *ov, real *oc, int nR, int *Rov, real *Roc, real coef)
{
	ASL_pfgh *asl;
	Hesoprod *h, **hp;
	real *w;

	/*DEBUG*/ if (ew->hop_free >= ew->hop_free_end) {
		fprintf(Stderr, "\n**** hop_free botch in new_Hesoprod!\n");
		exit(1);
		}
	h = ew->hop_free++;
	w = ew->w;
	h->ov = ov;
	h->ove = ov + nov;
	h->oc = oc;
	h->rv = Rov;
	h->rve = Rov + nR;
	h->roc = Roc;
	h->coef = coef;
	asl = (ASL_pfgh*)ew->asl;
	hp = (Hesoprod**)(w + asl->P.otodo) + *Rov;
	h->next = *hp;
	*hp = h;
	}

 static uHeswork*
new_uhw(EvalWorkspace *ew, range *r)
{
	uHeswork *rv = ew->uhw_free;
	ew->uhw_free = (uHeswork*)&rv->ogp[r->n];
	if (ew->uhw_free > ew->uhw_free_end) {
		/*DEBUG*/ fprintf(Stderr, "\n**** botch in new_uhw!\n");
		exit(1);
		}
	return rv;
	}

 static fint
bothadj(EvalWorkspace *ew, SputInfo *spi)
{
	/* Adjust to compute both triangles of Hessian */
	ASL_pfgh *asl;
	fint *hr, *hre, *hrn, *hrn0;
	int kz;
	size_t *hcs, *ucs;
	size_t i, i0, i1, j, k, k0, L, n, n1, nz;
	ssize_t nod, *ulc, *uli,  *z, *z0, *z1;

	asl = (ASL_pfgh*)ew->asl;
	n = n_var;
	if ((nod = spi->nod) >= 0) {
		if (!nod)
			return 0;
		goto done;
		}
	n1 = n + 1;
	hcs = spi->hcolstartsZ;
	nod = nz = hcs[n] - hcs[0];
	hr = spi->hrownos - 1;
	i = i0 = Fortran;
	for(j = i + n; i < j; i++, hcs++) {
		hr += k = hcs[1] - hcs[0];
		if (k && *hr == i)
			--nod;
		}
	/* nod = number of off-diagonal elements in upper triangle */
	if (!(spi->nod = nod))
		return 0;	/* diagonal Hessian */
	nz += nod;
	L = nz*sizeof(fint) + (2*(nod+n1))*sizeof(size_t);
	if (sizeof(size_t) != sizeof(fint))
		L += n1*sizeof(fint);
	spi->khinfob = kz = htcl(L);
	spi->ulinc0 = uli = (ssize_t*)new_mblk(kz);
	spi->ulcopy0 = ulc = uli + n1;
	spi->hcs[1] = hcs = (size_t*)(ulc + 2*nod);
	hrn0 = (fint*)(hcs + n1);
	if (sizeof(size_t) != sizeof(fint))
		hrn0 += n1;
	spi->hrn[1] = hrn0;
	z = z0 = (ssize_t*)new_mblk(kz = htcl(n*sizeof(ssize_t)));
	z1 = z - Fortran;
	ucs = spi->hcs[0];
	hre = spi->hrn[0];
	for(i = i0; i < j; i++, ucs++) {
		hr = hre;
		hre += *z++ = ucs[1] - ucs[0];
		while(hr < hre)
			if ((k = *hr++) != i)
				z1[k]++;
		}
	ucs = spi->hcs[0];
	hre = spi->hrn[0];
	*uli++ = 0;
	for(i = k = i0; i < j; i++, ucs++) {
		hr = hre;
		hre += L = ucs[1] - ucs[0];
		*hcs++ = k;
		k0 = k - i0;
		hrn = hrn0 + k0;
		*uli++ = z1[i] - L;
		k += z1[i];
		z1[i] = k0 + L;
		while(hr < hre)
			if ((i1 = *hrn++ = *hr++) != i) {
				*ulc++ = k0++;
				hrn0[*ulc++ = z1[i1]++] = i;
				}
		}
	*hcs = k;
	spi->ulcend = ulc;
	Del_mblk_ASL((ASL*)asl, z0);
	spi->ulinc = spi->ulinc0;
	spi->ulcopy = spi->ulcopy0;
 done:
	spi->hrownos = spi->hrn[1];
	spi->hcolstartsZ = spi->hcs[1];
	return nod;
	}

 static void
upper_to_lower(EvalWorkspace *ew, SputInfo *spi, size_t nz)
{	/* convert upper to lower triangular */

	ASL_pfgh *asl;
	fint *hrownos, *rn;
	int k, k1, *u0, *utoL;
	size_t L, *cs, *hcolstarts;
	ssize_t f, i, j, j1, j2, n, n1, *rs, *z;

	asl = (ASL_pfgh*)ew->asl;
	f = Fortran;
	n = n_var;
	n1 = n + 1;
	hrownos = spi->hrownos;
	hcolstarts = spi->hcolstartsZ;
	L = nz*sizeof(fint) + n1*sizeof(size_t);
	if (sizeof(size_t) != sizeof(fint))
		L += n1*sizeof(fint);
	spi->khinfob = k = htcl(L);
	spi->hcolstartsZ = cs = (size_t*)new_mblk(k);
	spi->ulinc0 = (ssize_t*)cs;
	rn = (fint*)(cs + n1);
	if (sizeof(size_t) != sizeof(fint))
		rn += n1;
	spi->hrownos = rn;
	k = htcl((n+nz)*sizeof(ssize_t));
	rs = (ssize_t*)new_mblk(k);
	z = rs + n;
	memset(rs, 0, n*sizeof(size_t));
	for(i = 0; i < nz; i++)
		rs[hrownos[i]-f]++;
	for(i = j = 0; i < n; i++) {
		cs[i] = j + f;
		j1 = rs[i];
		rs[i] = j;
		j += j1;
		}
	cs[n] = nz + f;
	j1 = hcolstarts[1] - f;
	for(i = j = 0; i < nz; i++) {
		while(i >= j1)
			j1 = hcolstarts[++j + 1] - f;
		rn[z[i] = rs[hrownos[i]-f]++] = j + f;
		}
	for(i = j = 0; i < nz; i++) {
		if ((j1 = z[i]) <= i) {
			if (j1 < 0)
				z[i] = -(j1 + 1);
			continue;
			}
		j += 3;
		while((j2 = z[j1]) != i) {
			z[j1] = -(j2 + 1);
			j++;
			j1 = j2;
			}
		}
	if (j) {
		j += 2;
		k1 = htcl(j*sizeof(int));
		spi->uptolow = utoL = (int*)new_mblk(k1);
		*utoL++ = k1;
		for(i = 0; i < nz; i++) {
			if ((j = z[i]) <= i)
				continue;
			u0 = utoL++;
			*utoL++ = i;
			*utoL++ = j;
			while((j2 = z[j]) != i) {
				z[j] = -(j2 + 1);
				j = *utoL++ = j2;
				}
			*u0 = (utoL - u0) - 1;
			}
		*utoL = 0;
		}
	Del_mblk_ASL((ASL*)asl, rs);
	}

#ifdef MULTIPLE_THREADS
 typedef struct RangeInfo RangeInfo;
 struct
RangeInfo {
	range *r;
	real *hsave;	/* for results of pshv_prod */
	uHeswork *uhw;
	int i1;		/* 0 <= i1 < r->n indicates r->lap[i1] */
	};

 typedef struct {
	ASL_pfgh *asl;
	RangeInfo *ri, *ri0, *ria, *rie;
	PTHREADS(pthread_mutex_t mutex;)
	EvalWorkspace *ew;	/* original workspace */
	SputInfo *spi;
	uHeswork *uhw, **utodo;
	int *rperm;
	int i, iow, ir, iy, nobj, nr, nzc, start;
	int *tdone;
	linarg **lap, **lape;
	real *ow, *y;
	} Thpars;

 typedef struct {
	int tno;		/* thread number */
	int i;			/* currently working on x[i] */
	int i1;			/* currently working on ri->r->lap[i1] */
	Thpars *tp;
	EvalWorkspace *ew;	/* thread workspace and ASL pointer */
	RangeInfo *ri;		/* current */
	uHeswork *uhw;		/* current (n >= nv) */
	} Thpars1;

 static int
ttodo(Thpars *tp, Thpars1 *tp1)
{
	RangeInfo *ri;
	linarg **lap;
	range *r;
	uHeswork *uhw;

	if (tp->start) {
		tp->start = 0;
		tp->ri = tp->ri0;
		goto starting;
		}
	if ((lap = tp->lap) && lap < tp->lape) {
		ri = tp1->ri = tp->ri;
		tp1->uhw = 0;
		tp1->i = tp->i;
		tp1->i1 = lap - ri->r->lap;
		tp->lap = lap + 1;
		goto ret;
		}
	tp->ri++;
 starting:
	if (tp->ri >= tp->rie) {
		tp1->ri = 0;
		return 0;
		}
	tp1->ri = ri = tp->ri;
	tp1->i = tp->i;
	r = ri->r;
	if ((uhw = ri->uhw) /*r->n >= r->nv*/) {
		tp1->uhw = uhw;
		tp->lap = 0;
		tp1->i1 = ri->i1++;
		}
	else {
		tp1->uhw = 0;
		tp->lap = r->lap + 1;
		tp->lape = r->lap + r->n;
		tp1->i1 = 0;
		}
 ret:
	return 1;
	}

 static void*
tstart1(void *arg)
{
	ASL_pfgh *asl;
	EvalWorkspace *ew;
	Ogptrs *q, *qe;
	RangeInfo *ri;
	Thpars *tp;
	Thpars1 *tp1;
	Varval *V, *v;
	int i, i1, j, n, nobj, ow, *ui, *uie, y;
	linarg *la, **lap, **lape;
	range *r;
	real *hs, *s, *si, *w;
	uHeswork *uhw;

	tp1 = arg;
	tp = tp1->tp;
	ow = tp->iow;
	y = tp->iy;
	nobj = tp->nobj;
	ew = tp1->ew;
	asl = (ASL_pfgh*)ew->asl;
	w = ew->w;
	V = (Varval*)w;
	s = w + asl->P.dOscratch;
 top:
	i = tp1->i;
	ri = tp1->ri;
	r = ri->r;
	lap = r->lap;
	lape = lap + (n = r->n);
	ui = r->ui;
	uie = ui + r->nv;
	while(ui < uie) {
		v = V + *ui++;
		v->aO = v->dO = v->adO = 0.;
		}
	hs = ri->hsave;
	if ((uhw = tp1->uhw)) {
		q = uhw->ogp;
		qe = q + r->n;
		si = s;
		do {
			if (q->ov < q->ove && *q->ov == i)
				*si = 1;
			si++;
			} while(++q < qe);
		pshv_prod1(ew, r, nobj, ow, y);
		q = uhw->ogp;
		si = s;
		do {
			if (q->ov < q->ove && *q->ov == i) {
				*si = 0;
				++q->ov;
				++q->oc;
				}
			si++;
			} while(++q < qe);
		}
	else {
		i1 = tp1->i1;
		si = s + i1;
		*si = 1;
		pshv_prod1(ew, r, nobj, ow, y);
		*si = 0;
		lap += i1;
		hs += i1*n - i1*(i1-1)/2;
		}
	do {
		la = *lap++;
		*hs++ =  V[la->u.v].aO;
		}
		while(lap < lape);
	pthread_mutex_lock(&tp->mutex);
#ifdef ALLOW_OPENMP
#pragma omp critical
 {
#endif
	j = ttodo(tp,tp1);
	pthread_mutex_unlock(&tp->mutex);
#ifdef ALLOW_OPENMP
 }
#endif
	if (j)
		goto top;
	return 0;
	}

 static void*
tstart2(void *arg)
{
	ASL_pfgh *asl;
	EvalWorkspace *ew;
	Ogptrs *q, *qe;
	RangeInfo *ri;
	Thpars *tp;
	Thpars1 *tp1;
	Varval *V, *v;
	int i, i1, j, n, nobj, *ui, *uie;
	linarg *la, **lap, **lape;
	range *r;
	real *hs, *ow, *s, *si, *w, *y;
	uHeswork *uhw;

	tp1 = arg;
	tp = tp1->tp;
	ow = tp->ow;
	y = tp->y;
	nobj = tp->nobj;
	ew = tp1->ew;
	asl = (ASL_pfgh*)ew->asl;
	w = ew->w;
	V = (Varval*)w;
	s = w + asl->P.dOscratch;
 top:
	i = tp1->i;
	ri = tp1->ri;
	r = ri->r;
	n = r->n;
	si = s;
	lap = r->lap;
	lape = lap + n;
	ui = r->ui;
	uie = ui + r->nv;
	while(ui < uie) {
		v = V + *ui++;
		v->aO = v->dO = v->adO = 0.;
		}
	hs = ri->hsave;
	if ((uhw = tp1->uhw)) {
		q = uhw->ogp;
		qe = q + n;
		do {
			if (q->ov < q->ove && *q->ov == i)
				*si = *q->oc;
			si++;
			} while(++q < qe);
		pshv_prod_ASL(ew, r, nobj, ow, y);
		q = uhw->ogp;
		si = s;
		do {
			if (q->ov < q->ove && *q->ov == i) {
				*si = 0;
				++q->ov;
				++q->oc;
				}
			si++;
			} while(++q < qe);
		}
	else {
		if ((i1 = tp1->i1)) {
			hs += i1*n - i1*(i1-1)/2;
			lap += i1;
			}
		si = s + i1;
		*si = 1;
		pshv_prod_ASL(ew, r, nobj, ow, y);
		*si = 0;
		}
	do {
		la = *lap++;
		*hs++ =  V[la->u.v].aO;
		}
		while(lap < lape);
	pthread_mutex_lock(&tp->mutex);
#ifdef ALLOW_OPENMP
#pragma omp critical
 {
#endif
	j = ttodo(tp,tp1);
	pthread_mutex_unlock(&tp->mutex);
#ifdef ALLOW_OPENMP
 }
#endif
	if (j)
		goto top;
	return 0;
	}

 EvalWorkspace*
ewalloc3_ASL(ASL_pfgh *asl, EvalWorkspace *ew, int keep)
{
	/* Similar to ewalloc2_ASL, but omitting H0 and returning storage from */
	/* Malloc when keep == 0 (so we can can free the storage) rather than */
	/* M1alloc (keep == 1). */

	EvalWorkspace *rv;
	Varval *v, *ve;
	arglist *al;
	char *s;
	const char **sa;
	func_info *fi;
	int i, nc, nf, nlogc, nlv, no, nv, nv0, nvx;
	real *d, *lastx, *ra, *w;
	size_t extra, L, La, Le, Lh, Lo, Lu, Luw, Lx, Lxc, *nxc;
	tfinfo **ptfi, *tfi;
#ifdef NANDEBUG
	typedef union {real r; unsigned int x[2]; } U;
	U u0, *u, *ue;
#endif

	switch (asl->i.ASLtype) {
	  /* case ASL_read_fgh: */
	  /* case ASL_read_pfg: */
	  case ASL_read_pfgh:
		break;
	  default:
		fprintf(Stderr, "\nUnexpected ASLtype = %d in ewalloc3(()\n", asl->i.ASLtype);
		fflush(Stderr);
		exit(1);
		}
	nc = nclcon;
	nlogc = n_lcon;
	no = n_obj;
	nf = asl->i.nfinv;
	La = nf * sizeof(arglist);
	Lo = (nc + no)*sizeof(real);
	Lu = asl->I.gscrx * sizeof(real);
	Lx = x0len;
	Le = (sizeof(EvalWorkspace) + sizeof(real) - 1) & ~(sizeof(real)-1);
	Lxc = ((2*(nc + no)+ nlogc)*sizeof(size_t) + sizeof(real) - 1) & ~(sizeof(real)-1);
	Lh = ew ? 0 : (asl->I.nhop * sizeof(Hesoprod) + sizeof(real) - 1) & ~(sizeof(real)-1);
/*	Luw = (asl->I.uhlen + sizeof(real) - 1) & ~(sizeof(real)-1); */
	Luw = 0;
	if ((nlv = nlvc) < nlvo)
		nlv = nlvo;
	L = La + Le + Lh + Lo + Lu + Luw + Lx + asl->i.derplen
		+ asl->i.wlen + asl->i.numlen + Lxc
		+ asl->i.ra_max*sizeof(real) + asl->i.sa_max*sizeof(char*);
	if (L & (sizeof(real) - 1)) {
		extra = sizeof(real) - (L & (sizeof(real) - 1));
		L += extra;
		}
	rv = (EvalWorkspace*) (keep ? M1alloc(L) : Malloc(L));
	ACQUIRE_MBLK_LOCK(&asl->i, MemLock);
	++asl->i.n_ew0;
	*asl->i.pewthread = rv;
	asl->i.pewthread = &rv->ewthread;
	FREE_MBLK_LOCK(&asl->i, MemLock);
	s = (char*)rv + Le;
	memset(rv, 0, Le + Lxc);
	rv->asl = (ASL*)asl;
	rv->nlv = nlv;
	rv->nxval = 1;
	rv->x0kind = ASL_first_x;
	nxc = (size_t*)s;
	rv->ncxval = nxc;
	rv->ncxgval = nxc += nc;
	rv->nlxval = nxc += nc;
	rv->noxval = nxc += nlogc;
	rv->noxgval = nxc += no;
	s += Lxc;
	if (ew)
		rv->hop_free = rv->hop_free0 = 0;
	else {
		rv->hop_free0 = rv->hop_free = (Hesoprod*)s;
		rv->hop_free_end = rv->hop_free0 + asl->I.nhop;
		s += Lh;
		}
	rv->uhw_free = rv->uhw_free0 = rv->uhw_free_end = 0;
/*	rv->uhw_free = rv->uhw_free0 =(uHeswork*)s;	*/
/*	rv->uhw_free_end = (uHeswork*)((char*)rv->uhw_free0 + asl->I.uhlen);	*/
/*	s += Luw;					*/
	rv->H0 = 0;
	lastx = (real*)s;
	s += Lx;
	rv->unopscr = (real*)s;
	s += Lu;
	rv->oyow0 = (real*)s;
	s += Lo;
	memcpy(s, asl->i.numvals, asl->i.numlen);
	rv->w = w = (real*)(s += asl->i.numlen);
#ifdef NANDEBUG
	u0.x[0] = 0x12345678;
	u0.x[1] = 0xfff7abcd;	/* signalling NaN */
	u = (U*)w;
	ue = (U*)(w + asl->P.rtodo);
	while(u < ue)
		*u++ = u0;
	__fpu_control = _FPU_IEEE - _FPU_EXTENDED + _FPU_DOUBLE - _FPU_MASK_IM;
	_FPU_SETCW(__fpu_control);
	signal(SIGFPE, fpecatch_ASL);
#endif
	nv = asl->i.n_var_;

	/* The following loop only matters under obscure conditions when setting */
	/* up Hessian compuations involving defined variables having linear parts */
	/* and used linearly in a constraint or objective. */
	/* With default AMPL settings, this situation does not arise. */
	v = (Varval*)w;
	for(ve = v + nv; v < ve; ++v)
		v->dO = v->aO = 0.;

	memset(w + asl->P.rtodo, 0, asl->P.zaplen);
	rv->Lastx = lastx;
	nv0 = asl->i.n_var0;
	nvx = nv0 + asl->i.nsufext[ASL_Sufkind_var];
	rv->dv = (Varval*)w + nvx;
	rv->dv1 = (Varval*)rv->dv + combc + como;
	s += asl->i.wlen;
	rv->derps = d = (real*)s;
	if ((i = asl->i.maxvar - asl->i.defvar0) > 0)
		memset(d + asl->i.defvar0, 0, i*sizeof(real));
	s += asl->i.derplen;
/*	rv->derpzap = d + nvx + combc + como + comc1 + como1; */
	rv->wantderiv = want_derivs;
	rv->ndhmax = asl->P.ndhmax;
	if (asl->i.X0_) {
		memcpy(w, asl->i.X0_, nv0*sizeof(real));
		if (nv > nv0)
			memset(w + nv0, 0, (nv-nv0)*sizeof(real));
		}
	else
		memset(w, 0, nv*sizeof(real));
	if (nf) {
		ra = (real*)s;
		sa = (const char**)(ra + asl->i.ra_max);
		rv->al = al = (arglist*)(sa + asl->i.sa_max);
		ptfi = (tfinfo**)asl->i.invd;
		memset(al, 0, nf*sizeof(arglist));
		for(i = 0; i < nf; ++i, ++al) {
			tfi = ptfi[i];
			al->n = al->nin = tfi->n;
			if ((al->nr = tfi->nr))
				al->ra = ra;
			al->at = tfi->at;
			al->dig = tfi->wd;
			if ((al->nsin = tfi->n - tfi->nr))
				al->sa = sa;
			fi = tfi->fi;
			al->f = (function*)fi;
			al->funcinfo = fi->funcinfo;
			al->AE = asl->i.ae;
			}
		}
	if (asl->P.wkinit0)
		wk_init_ASL(w, asl->P.wkinit0, 0.);
	if (asl->P.wkinit2)
		wk_init_ASL(w, asl->P.wkinit2, 2.);
	if (asl->P.wkinitm1)
		wk_init_ASL(w, asl->P.wkinitm1, -1.);
	if (ew)
		memcpy(rv->w, ew->w, asl->i.wlen);
	return rv;
	}
#endif

 fint
sphes_setup_ew_ASL(EvalWorkspace *ew, SputInfo **pspi, int nobj, int ow, int y, int uptri)
{
	ASL *a;
	ASL_pfgh *asl;
	Hesoprod *hop, *hop1, **otodo, **otodoi, **otodoj;
	Objrep **por;
	Ogptrs *q, *qe;
	SputInfo *spi, *spi1;
#ifdef MULTIPLE_THREADS
	RangeInfo *ri, *ri0, *rie;
	Thpars tp;
	Thpars1 *tp1, *tpi;
	int m, maxit = 1, nr;
	real *hs, *hs0;
	size_t nrng;
#ifndef ALLOW_OPENMP
	pthread_t *T;
	int it;
	void *rc;
#endif
#endif
	Varval *V, *v;
	fint *hr, *hre, *hrownos, rv;
	int dassume, i, j, khinfo, kz;
	int n, n1, nhinfo, nlc0, no, noe, nov, nov1, nqslim, nzc, nzc0;
	int rfilter, robjno, *ov, *ov1, *ove, *ui, ushow, *zc, *zci;
	linarg *la, **lap, **lap1, **lape;
	ps_func *p, *pe;
	psb_elem *b;
	psg_elem *g, *ge;
	range *r, *r0, **rnext, **rp, **rtodo;
	real *s, *si, t, *w;
	size_t *hcolstarts, iz, jz, n1spi, *tf;
	uHeswork *uhw, *uhwi, **utodo, **utodoi, **utodoj;

	asl = (ASL_pfgh*)(a = ew->asl);
	ASL_CHECK(a, ASL_read_pfgh, "sphes_setup");
	w = ew->w;
	V = (Varval*)w;
	j = asl->i.maxvar;
	dassume = asl->P.sph_opts & 1;
	for(i = asl->i.defvar0 + asl->P.ncom; i < j; ++i) {
		v = V + i;
		v->dO = v->aO = v->adO = 0.;
		}
	n1 = n_var + 1;
	if (!pspi)
		pspi = &ew->Sputinfo;
	nlc0 = asl->i.nlc0;
	if (nobj >= 0 && nobj < n_obj) {
		robjno = -2 - nobj;
		rfilter = n_obj > 1 || (!y && nlc0 > 0);
		ow = 0;
		no = nobj;
		noe = no + 1;
		}
	else {
		robjno = nobj = -1;
		rfilter = (!ow && n_obj > 0) || (!y && nlc0 > 0);
		no = noe = 0;
		if (ow) {
			noe = n_obj;
			ow = 1;
			if ((por = asl->i.Or)) {
				for(i = 0; i < noe; ++i) {
					if (por[i]) {
						y = 1;
						break;
						}
					}
				}
			}
		}
	if (y)
		y = 1;
	if ((n = nlvo) < nlvc)
		n = nlvc;
	ushow = asl->P.sph_opts & 2;
	if ((spi = *pspi)) {
		if (spi->ow == ow && spi->y == y && spi->nobj == nobj
		 && spi->uptri == uptri)
			goto done;
		if (spi->ulinc0)
			del_mblk(spi->ulinc0);
		if ((ui = spi->uptolow))
			del_mblk(ui);
		del_mblk(spi);
		*pspi = 0;
		}
	ew->hes_setup_called = 3;
	otodo = otodoi = (Hesoprod**)(w + asl->P.otodo);
	rtodo = (range**)(w + asl->P.rtodo);
	memset(rtodo, 0,  asl->P.nran*sizeof(range*));
	utodo = utodoi = (uHeswork**)(w + asl->P.utodo);
	s = w + asl->P.dOscratch;
	nqslim = n >> 3;
	kz = htcl(2*sizeof(int)*n + asl->P.nran*sizeof(range*));
	rnext = (range**)new_mblk_ASL(a, kz);
	zc = (int*)(rnext + asl->P.nran);
	zci = zc + n;
	memset(zc, 0, n*sizeof(int));
	n1spi = sizeof(SputInfo) + n1*sizeof(size_t);
	if (sizeof(size_t) != sizeof(fint))
		n1spi += n1*sizeof(fint);
	khinfo = htcl((n1 + 30)*sizeof(fint) + n1spi);
	spi = (SputInfo*)new_mblk_ASL(a, khinfo);
	hcolstarts = (size_t*)(spi+1);
	hr = hrownos = (fint*)((char*)spi + n1spi);
	nhinfo = ((sizeof(Char*)<<khinfo) - n1spi) / sizeof(fint);
	hre = hr + nhinfo;
	r0 = (range*)&asl->P.rlist;
	ew->hop_free = ew->hop_free0;
	ew->uhw_free = ew->uhw_free0;
	for(r = asl->P.rlist.next; r != r0; r = r->rlist.next) {
		if (r->n <= 0)
			continue;
		if (rfilter) {
			for(b = r->refs; b; b = b->next) {
				if (b->conno >= 0) {
					if (y)
						goto keep;
					}
				else if (b->conno == robjno)
					goto keep;
				}
			continue;
			}
 keep:
		rp = rtodo + r->lasttermno;
		rnext[r->irange] = *rp;
		*rp = r;
		}
	if (asl->P.nobjgroups) {
		for(i = no; i < noe; i++) {
			p = asl->P.ops + i;
			g = p->g;
			for(ge = p->ge; g < ge; g++) {
				if ((nov = g->nov)) {
					ov = g->ov;
					new_Hesoprod(ew, nov, ov, 0, nov, ov, 0, 1.);
					}
				}
			}
		}
	if (asl->P.ncongroups && y) {
		p = asl->P.cps;
		for(pe = p + nlc0; p < pe; p++)
			for(g = p->g, ge = p->ge; g < ge; g++) {
				if ((nov = g->nov)) {
					ov = g->ov;
					new_Hesoprod(ew, nov, ov, 0, nov, ov, 0, 1.);
					}
				}
		}
#ifdef MULTIPLE_THREADS
	nrng = asl->P.nran;
	if ((m = asl->P.hessetupth) <= 0 || (asl->P.sph_opts & 5) || nrng <= 0) {
		m = 0;
		hs0 = 0;
		}
	else {
		if (m > nrng)
			m = nrng;
		tp.asl = asl;
#ifdef ALLOW_OPENMP
		/*omp_set_num_threads(m);*/
#else
		pthread_mutex_init(&tp.mutex, 0);
#endif
		tp.nobj = nobj;
		tp.iow = ow;
		tp.iy = y;
		tp.utodo = utodo;
		hs0 = hs = (real*)new_mblk(htcl(m*(PTHREADS(sizeof(pthread_t*) +)
						   sizeof(Thpars1))
					+ asl->P.thlen*sizeof(real)
					+ nrng*(sizeof(RangeInfo))));
		tp1 = (Thpars1*)(hs + asl->P.thlen);
		tp.ew = ew;
		tp.spi = spi;
#ifdef ALLOW_OPENMP
		tp.ri0 = ri0 = (RangeInfo*)(tp1 + m);
#else
		T = (pthread_t*)(tp1 + m);
		T[0] = pthread_self();
		tp.ri0 = ri0 = (RangeInfo*)(T + m);
#endif
		tp.rie = rie = ri0 + nrng;
		tp1->tno = 0;
		tp1->tp = &tp;
		tp1->ew = ew;
		for(i = 0; i < m; ++i) {
			tpi = tp1 + i;
			tpi->tno = i;
			tpi->tp = &tp;
			if (i)
				tpi->ew = ewalloc3_ASL(asl,0,0);
			}
		}
#endif
	for(i = 0; i < n; i++) {
		nzc = 0;
		rp = rtodo;
		uhwi = *utodoi;
		*utodoi++ = 0;
#ifdef MULTIPLE_THREADS /*{*/

		if (m) {
			ri = ri0;
			hs = hs0;
			while((r = *rp)) {
				rp = rnext + r->irange;
				nr = r->n;
				if (nr >= r->nv) {
					uhw = new_uhw(ew, r);
					uhw->next = uhwi;
					uhwi = uhw;
					uhw->r = r;
					uhw->uie = (uhw->ui = r->ui) + r->nv;
					q = uhw->ogp;
					lap = r->lap;
					lape = lap + nr;
					while(lap < lape) {
						la = *lap++;
						q->ove = (q->ov = la->ov) + la->nnz;
						q->oc = la->oc;
						++q;
						}
					}
				else {
					ri->r = r;
					ri->uhw = 0;
					ri->hsave = hs;
					hs += nr*(nr+1)/2;
					++ri;
					}
				}
			tp.ria = ri;
			while((uhw = uhwi)) {
				uhwi = uhw->next;
				r = ri->r = uhw->r;
				ri->uhw = uhw;
				ri->hsave = hs;
				hs += r->n;
				++ri;
				}
			if (tp.ri0 >= (tp.rie = ri))
				goto uhw_done;
			tp.i = i;
			tp.start = 1;
#ifdef ALLOW_OPENMP
#pragma omp parallel num_threads(m)
 {
			int it, jt;
			Thpars1 *tpi;

			it = omp_get_thread_num();
			tpi = tp1 + it;
#pragma omp critical
 {
			if ((jt = ttodo(&tp,tpi)) && maxit <= it)
				maxit = it + 1;
 }
			if (jt)
				tstart1(tpi);
 }
#else
			pthread_mutex_lock(&tp.mutex);
			for(it = 1; it < m; ++it) {
				tpi = tp1 + it;
				tpi->ri = 0;
				if (!ttodo(&tp,tpi))
					break;
				if ((j = pthread_create(&T[it], 0, &tstart1, tpi)))
					fprintf(Stderr, "Return %d from pthread_create for tstart1.\n", j);
				}
			maxit = it;
			j = ttodo(&tp,tp1);
			pthread_mutex_unlock(&tp.mutex);
			if (j) {
				tp1->tp = &tp;
				tstart1(tp1);
				}
			while(it > 1) {
				if ((j = pthread_join(T[--it], &rc)))
					fprintf(Stderr,
						"Return %d from pthread_join for thread %d\n",
						j, it);
				}
#endif
			for(ri = ri0; ri < tp.ria; ++ri) {
				r = ri->r;
				lap = r->lap;
				lape = lap + r->n;
				hs = ri->hsave;
				while(lap < lape) {
					lap1 = lap++;
					la = *lap1++;
					ov = la->ov;
					nov = la->nnz;
					if ((t = *hs++))
						new_Hesoprod(ew,nov,ov,0,nov,ov,0,t);
					while(lap1 < lape) {
					    la = *lap1++;
					    if ((t = *hs++)) {
						ov1 = la->ov;
						nov1 = la->nnz;
						new_Hesoprod(ew,nov,ov,0,nov1,ov1,0,t);
						new_Hesoprod(ew,nov1,ov1,0,nov,ov,0,t);
						}
					    }
					}
				}
			for(; ri < tp.rie; ++ri) {
				uhw = ri->uhw;
				r = uhw->r;
				lap = r->lap;
				lape = lap + r->n;
				hs = ri->hsave;
				kz = 0;
				nzc0 = nzc;
				do {
					if (*hs++) {
						++kz;
						la = *lap;
						ov = la->ov;
						for(ove = ov + la->nnz; ov < ove; ++ov)
							if ((j = *ov) <= i
							 && !zc[j]++)
								zci[nzc++] = j;
						}
					}
					while(++lap < lape);
				if (ushow) {
					if (r->n > r->nv)
						printf("obrange %d\t%d", r->n, r->nv);
					else
						printf("sqrange %d", r->n);
					printf("\t%d\t%d\n", kz, nzc-nzc0);
					}
				if ((ui = ++uhw->ui) < uhw->uie) {
					utodoj = utodo + *ui;
					uhw->next = *utodoj;
					*utodoj = uhw;
					}
				}
			}
		else
#endif /*}*/
		while((r = *rp)) {
			rp = rnext + r->irange;
			lap = r->lap;
			lape = lap + r->n;
			if (r->n >= r->nv) {
				uhw = new_uhw(ew, r);
				uhw->next = uhwi;
				uhwi = uhw;
				uhw->r = r;
				uhw->ui = ui = r->ui;
				uhw->uie = ui + r->nv;
				q = uhw->ogp;
				while(lap < lape) {
					la = *lap++;
					q->ove = (q->ov = la->ov) + la->nnz;
					++q;
					}
				}
			else {
				si = s;
				while(lap < lape) {
					*si = 1;
					pshv_prod1(ew, r, nobj, ow, y);
					*si++ = 0;
					lap1 = lap++;
					la = *lap1++;
					ov = la->ov;
					nov = la->nnz;
					v = V + la->u.v;
					if ((t = v->aO))
						new_Hesoprod(ew,nov,ov,0,nov,ov,0,t);
					while(lap1 < lape) {
					    la = *lap1++;
					    v = V + la->u.v;
					    if ((t = v->aO)) {
						ov1 = la->ov;
						nov1 = la->nnz;
						new_Hesoprod(ew,nov,ov,0,nov1,ov1,0,t);
						new_Hesoprod(ew,nov1,ov1,0,nov,ov,0,t);
						}
					    }
					}
				}
			}
#ifdef MULTIPLE_THREADS
		if (m)
			goto uhw_done;
#endif
		while((uhw = uhwi)) {
			uhwi = uhwi->next;
			si = s;
			r = uhw->r;
			q = uhw->ogp;
			qe = q + r->n;
			do {
				if (q->ov < q->ove && *q->ov == i)
					*si = 1.;
				si++;
				} while(++q < qe);
			lap = r->lap;
			lape = lap + r->n;
			if (dassume) {
				if (ushow) {
					if (r->n > r->nv)
						printf("obrange %d %d\n", r->n, r->nv);
					else
						printf("sqrange %d\n", r->n);
					}
				do {
					la = *lap++;
					ov = la->ov;
					for(ove = ov + la->nnz; ov < ove; ++ov)
						if ((j = *ov) <= i
						 && !zc[j]++)
							zci[nzc++] = j;
					}
					while(lap < lape);
				}
			else {
				pshv_prod1(ew, r, nobj, ow, y);
				kz = 0;
				nzc0 = nzc;
				do {
					la = *lap++;
					if (V[la->u.v].aO) {
						++kz;
						ov = la->ov;
						for(ove = ov + la->nnz; ov < ove; ++ov)
							if ((j = *ov) <= i
							 && !zc[j]++)
								zci[nzc++] = j;
						}
					}
					while(lap < lape);
				if (ushow) {
					if (r->n > r->nv)
						printf("obrange %d\t%d", r->n, r->nv);
					else
						printf("sqrange %d", r->n);
					printf("\t%d\t%d\n", kz, nzc-nzc0);
					}
				}

			q = uhw->ogp;
			si = s;
			do {
				if (q->ov < q->ove && *q->ov == i) {
					*si = 0;
					++q->ov;
					}
				si++;
				} while(++q < qe);
			if ((ui = ++uhw->ui) < uhw->uie) {
				utodoj = utodo + *ui;
				uhw->next = *utodoj;
				*utodoj = uhw;
				}
			}
#ifdef MULTIPLE_THREADS
 uhw_done:
#endif
		*rtodo++ = 0;	/* reset */
		hop1 = *otodoi;
		*otodoi++ = 0;
		while((hop = hop1)) {
			hop1 = hop->next;
			ov = hop->ov;
			ove = hop->ove;
			while((j = *ov) <= i) {
				if (!zc[j]++)
					zci[nzc++] = j;
				if (++ov >= ove)
					break;
				}
			ov1 = hop->rv;
			ove = hop->rve;
			if (++ov1 < ove) {
				hop->rv = ov1;
				otodoj = otodo + *ov1;
				hop->next = *otodoj;
				*otodoj = hop;
				}
			}
		hcolstarts[i] = hr - hrownos;
		if (nzc > hre - hr) {
			spi1 = (SputInfo*)new_mblk_ASL(a, ++khinfo);
			tf = (size_t*)(spi1+1);
			memcpy(tf, hcolstarts, (char*)hr - (char*)hcolstarts);
			del_mblk(spi);
			spi = spi1;
			hcolstarts = tf;
			hrownos = (fint*)((char*)spi1 + n1spi);
			hr = hrownos + hcolstarts[i];
			nhinfo = ((sizeof(Char*)<<khinfo) - n1spi) / sizeof(fint);
			hre = hrownos + nhinfo;
			}
		if (nzc > nqslim) {
			for(j = 0; j < n; j++)
				if (zc[j])
					zc[*hr++ = j] = 0;
			}
		else {
			if (nzc > 1)
				qsort(zci, nzc, sizeof(int), compar);
			for(j = 0; j < nzc; j++)
				zc[*hr++ = zci[j]] = 0;
			}
		}
	jz = hcolstarts[n] = hr - hrownos;
	for(i = n; ++i < n1; )
		hcolstarts[i] = jz;
	if ((j = Fortran)) {
		for(i = 0; i < n1; i++)
			hcolstarts[i] += j;
		iz = hcolstarts[n] - j;
		while(iz)
			hrownos[--iz] += j;
		}
	spi->hcs[0] = hcolstarts;
	spi->hrn[0] = hrownos;
	spi->nod = -1;
	spi->ulcend = 0;
	spi->khinfo = khinfo;
	spi->nobj = nobj;
	spi->ow = ow;
	spi->y = y;
	spi->uptri = uptri;
	*pspi = spi;
	spi->ulinc0 = spi->ulinc = 0;
	spi->ulcopy = 0;
	spi->uptolow = 0;
	del_mblk(rnext);
#ifdef MULTIPLE_THREADS
	if (hs0)
		del_mblk(hs0);
#endif
 done:
#ifdef MULTIPLE_THREADS
	asl->P.thusedsetup = maxit;
	if (asl->P.hessetupth > 0 && (asl->P.sph_opts & 16)) {
		printf("threads for sph_setup = %d\n", maxit);
		fflush(stdout);
		}
	asl->P.thused = 0;
#endif
	spi->hrownos = spi->hrn[0];
	spi->hcolstartsZ = hcolstarts = spi->hcs[0];
	rv = hcolstarts[n] - hcolstarts[0];
	switch(uptri) {
	  case 0:
		rv += bothadj(ew, spi);
		break;
	  case 2:
		upper_to_lower(ew, spi, rv);
	  }
	hcolstarts = spi->hcolstartsZ;
	if (sizeof(size_t) == sizeof(fint))
		spi->hcolstarts = (fint*)hcolstarts;
	else {
		spi->hcolstarts = hr = (fint*)(hcolstarts + n1);
		for(i = 0; i < n1; ++i)
			hr[i] = hcolstarts[i];
		}
	return rv;
	}

 void
sphes_ew_ASL(EvalWorkspace *ew, SputInfo **pspi, real *H, int nobj, real *ow, real *y)
{
	/* sparse upper triangle of Hessian */

	ASL *a;
	ASL_pfgh *asl;
	EvalWorkspace *ew0;
	Ogptrs *q, *qe;
	Hesoprod *hop, *hop1, **otodo, **otodoi, **otodoj;
	SputInfo*spi, *spi0;
#ifdef MULTIPLE_THREADS
	RangeInfo *ri, *ri0, *rie;
	Thpars tp;
	Thpars1 *tp1, *tpi;
	int m, maxit = 1, nr;
	real *hs, *hs0;
	size_t nrng;
#ifndef ALLOW_OPENMP
	pthread_t *T;
	int it;
	void *rc;
#endif
#endif
	Varval *V, *v;
	fint *hr;
	int i, j, k, n, no, noe, nov, nov1, *ov, *ov1, *ove, *ui, *uie, uptri;
	linarg *la, **lap, **lap1, **lape;
	ps_func *p, *pe;
	psg_elem *g, *ge;
	range *r, *r0, **rnext, **rp, **rtodo;
	real *Hi, *H0, *H00, *cscale, *oc, *oc1, *owi, *s, *si;
	real t, t1, *vsc, *vsc0, *vsc1, *w, *y1;
	size_t *hcs;
	ssize_t *ulc, *uli;
	uHeswork *uhw, *uhwi, **utodo, **utodoi, **utodoj;

	asl = (ASL_pfgh*)(a = ew->asl);
	ASL_CHECK(a, ASL_read_pfgh, "sputhes");
	++ew->stats.hesmat;
	if (ew->Derrs)
		deriv2_errchk_ASL(ew, 3);
	w = ew->w;
	V = (Varval*)w;
	xpsg_check_ASL(ew, nobj, ow, y);
	if (!pspi)
		pspi = &ew->Sputinfo;
	i = j = 0;
	if (y)
		j = 1;
	if (nobj >= 0 && nobj < n_obj) {
		no = nobj;
		noe = no + 1;
		owi = ow ? ow + no : &edag_one_ASL;
		ow = 0;
		}
	else {
		nobj = -1;
		no = noe = 0;
		if ((owi = ow)) {
			noe = n_obj;
			i = 1;
			}
		}
	if (ew->hes_setup_called != 3 || !(spi = *pspi)) {
		uptri = 0;
		if ((ew0 = asl->i.Ew0) && (spi0 = ew0->Sputinfo))
			uptri = spi0->uptri;
		sphes_setup_ew_ASL(ew, pspi, nobj, ow != 0, y != 0, uptri);
		}
	else if ((spi->nobj != nobj && nobj >= 0) || spi->ow < i || spi->y < j) {
		fprintf(Stderr,
		 "\nsphes() call inconsistent with previous sphsetup()\n");
		exit(1);
		}
	spi = *pspi;
	otodo = otodoi = (Hesoprod**)(w + asl->P.otodo);
	rtodo = (range**)(w + asl->P.rtodo);
	memset(rtodo, 0,  asl->P.nran*sizeof(range*));
	utodo = utodoi = (uHeswork**)(w + asl->P.utodo);
	s = w + asl->P.dOscratch;
	n = ew->nlv;
	Hi = H0 = ew->H0;
	memset(Hi, 0, n*sizeof(real) + asl->P.nran*sizeof(range*));
	rnext = (range**)(Hi + n);
	H0 -= Fortran;
	r0 = (range*)&asl->P.rlist;
	ew->hop_free = ew->hop_free0;
	ew->uhw_free = ew->uhw_free0;
	for(r = asl->P.rlist.next; r != r0; r = r->rlist.next) {
		if ((j = r->n) <= 0)
			continue;
		i = r->lasttermno;
		rp = rtodo + i;
		rnext[r->irange] = *rp;
		*rp = r;
		}
	if (asl->P.nobjgroups)
	    for(; no < noe; no++)
		if ((t = *owi++)) {
		    p = asl->P.ops + no;
		    g = p->g;
		    for(ge = p->ge; g < ge; g++) {
			if ((t1 = t*w[g->gm + 2])) {
				ov = g->ov;
				oc = w + g->gm + 3;
				nov = g->nov;
				for(i = 0; i < nov; ++i) {
					if (oc[i]) {
						nov -= i;
						ov += i;
						oc += i;
						new_Hesoprod(ew, nov, ov, oc, nov, ov, oc, t1);
						break;
						}
					}
				}
			}
		}
	if (asl->P.ncongroups && y) {
		cscale = asl->i.lscale;
		p = asl->P.cps;
		y1 = y;
		for(pe = p + asl->i.nlc0; p < pe; p++, y1++) {
			if ((t = cscale ? *cscale++ * *y1 : *y1)) {
				for(g = p->g, ge = p->ge; g < ge; g++) {
				    if ((t1 = t*w[g->gm + 2])) {
					ov = g->ov;
					oc = w + g->gm + 3;
					nov = g->nov;
					for(i = 0; i < nov; ++i)
					    if (oc[i]) {
						nov -= i;
						ov += i;
						oc += i;
						new_Hesoprod(ew, nov, ov, oc, nov, ov, oc, t1);
						break;
						}
					}
				    }
				}
			}
		}
	hcs = spi->hcs[0];
	hr  = spi->hrn[0];
	uli = spi->ulinc;
	H00 = H;
	vsc = asl->i.vscale;
	vsc0 = vsc - Fortran;
	vsc1 = vsc;
	if (ew->x0kind & ASL_need_funnel)
		funnelhes_ew_ASL(ew);
#ifdef MULTIPLE_THREADS
	if ((m = asl->P.hesevalth) < 1 || (asl->P.sph_opts & 8) || n <= 0) {
		m = 0;
		hs0 = 0;
		}
	else {
		if (m > n)
			m = n;
		tp.asl = asl;
#ifdef ALLOW_OPENMP
		/*omp_set_num_threads(m);*/
#else
		pthread_mutex_init(&tp.mutex, 0);
#endif
		tp.nobj = nobj;
		tp.ow = ow;
		tp.y = y;
		tp.utodo = utodo;
		nrng = asl->P.nran;
		if ((hs = asl->P.hesparwk) && m > asl->P.hesmaxth) {
			free(hs);
			*asl->P.hestofree = 0;
			hs0 = 0;
			goto alloc_anew;
			}
		if (!(hs0 = hs = asl->P.hesparwk)) {
 alloc_anew:
			asl->P.hesparwk =
			hs = (real*)M1alloc(m*(PTHREADS(sizeof(pthread_t*) +)
						sizeof(Thpars1))
					+ asl->P.thlen*sizeof(real)
					+ nrng*(sizeof(RangeInfo)));
			asl->P.hestofree = (real**)(asl->i.Mbnext - 1); /* kludge */
			asl->P.hesmaxth = m;
			}
		tp1 = (Thpars1*)(hs + asl->P.thlen);
		tp.ew = ew;
		tp.spi = spi;
#ifdef ALLOW_OPENMP
		tp.ri0 = ri0 = (RangeInfo*)(tp1 + m);
#else
		T = (pthread_t*)(tp1 + m);
		T[0] = pthread_self();
		tp.ri0 = ri0 = (RangeInfo*)(T + m);
#endif
		tp.rie = rie = ri0 + nrng;
		tp1->tno = 0;
		tp1->tp = &tp;
		tp1->ew = ew;
		if (!hs0) {
			hs0 = hs;
			for(i = 0; i < m; ++i) {
				tpi = tp1 + i;
				tpi->tno = i;
				tpi->tp = &tp;
				if (i)
					tpi->ew = ewalloc3_ASL(asl,ew,1);
				}
			}
		}
#endif
	for(i = 0; i < n; i++) {
		rp = rtodo;
		uhwi = *utodoi;
		*utodoi++ = 0;
#ifdef MULTIPLE_THREADS /*{*/
		if (m) {
			ri = ri0;
			hs = hs0;
			while((r = *rp)) {
				rp = rnext + r->irange;
				nr = r->n;
				if (nr >= r->nv) {
					uhw = new_uhw(ew, r);
					uhw->next = uhwi;
					uhwi = uhw;
					uhw->r = r;
					uhw->uie = (uhw->ui = r->ui) + r->nv;
					q = uhw->ogp;
					lap = r->lap;
					lape = lap + r->n;
					while(lap < lape) {
						la = *lap++;
						q->ove = (q->ov = la->ov) + la->nnz;
						q->oc = la->oc;
						++q;
						}
					}
				else {
					ri->r = r;
					ri->uhw = 0;
					ri->hsave = hs;
					hs += nr*(nr+1)/2;
					ri->i1 = 0;
					++ri;
					}
				}
			tp.ria = ri;
			while((uhw = uhwi)) {
				uhwi = uhw->next;
				r = ri->r = uhw->r;
				ri->uhw = uhw;
				ri->hsave = hs;
				hs += r->n;
				ri->i1 = 0;
				++ri;
				}
			if (tp.ri0 >= (tp.rie = ri))
				goto uhw_done;
			tp.i = i;
			tp.start = 1;
#ifdef ALLOW_OPENMP
#pragma omp parallel num_threads(m)
 {
			int it, jt;
			Thpars1 *tpi;

			it = omp_get_thread_num();
			tpi = tp1 + it;
#pragma omp critical
 {
			if ((jt = ttodo(&tp,tpi)) && maxit <= it)
				maxit = it + 1;
 }
			if (jt)
				tstart2(tpi);
 }
#else
			pthread_mutex_lock(&tp.mutex);
			for(it = 1; it < m; ++it) {
				tpi = tp1 + it;
				tpi->ri = 0;
				if (!ttodo(&tp,tpi))
					break;
				if ((j = pthread_create(&T[it], 0, &tstart2, tpi)))
					fprintf(Stderr, "Return %d from pthread_create.\n", j);
				}
			maxit = it;
			j = ttodo(&tp,tp1);
			pthread_mutex_unlock(&tp.mutex);
			if (j) {
				tp1->tp = &tp;
				tstart2(tp1);
				}
			while(it > 1) {
				if ((j = pthread_join(T[--it], &rc)))
					fprintf(Stderr,
						"Return %d from pthread_join for thread %d\n",
						j, it);
				}
#endif
			for(ri = ri0; ri < tp.ria; ++ri) {
				r = ri->r;
				lap = r->lap;
				lape = lap + r->n;
				hs = ri->hsave;
				while(lap < lape) {
					lap1 = lap++;
					la = *lap1++;
					ov = la->ov;
					oc = la->oc;
					nov = la->nnz;
					if ((t = *hs++))
						new_Hesoprod(ew,nov,ov,oc,nov,ov,oc,t);
					while(lap1 < lape) {
					    la = *lap1++;
					    if ((t = *hs++)) {
						ov1 = la->ov;
						oc1 = la->oc;
						nov1 = la->nnz;
						new_Hesoprod(ew,nov,ov,oc,nov1,ov1,oc1,t);
						new_Hesoprod(ew,nov1,ov1,oc1,nov,ov,oc,t);
						}
					    }
					}
				}
			for(; ri < tp.rie; ++ri) {
				uhw = ri->uhw;
				r = uhw->r;
				lap = r->lap;
				lape = lap + r->n;
				hs = ri->hsave;
				do {
					if ((t = *hs++)) {
						la = *lap;
						ov = la->ov;
						oc = la->oc;
						for(ov1 = ov + la->nnz; ov < ov1; ++ov, ++oc)
							if ((j = *ov) <= i)
								Hi[j] += t**oc;
						}
					}
					while(++lap < lape);
				if ((ui = ++uhw->ui) < uhw->uie) {
					utodoj = utodo + *ui;
					uhw->next = *utodoj;
					*utodoj = uhw;
					}
				}
			}
		else
#endif /*}*/
			{
			while((r = *rp)) {
				rp = rnext + r->irange;
				lap = r->lap;
				lape = lap + r->n;
				ui = r->ui;
				uie = ui + r->nv;
				while(ui < uie) {
					v = V + *ui++;
					v->aO = v->dO = v->adO = 0.;
					}
				if (r->n >= r->nv) {
					uhw = new_uhw(ew, r);
					uhw->next = uhwi;
					uhwi = uhw;
					uhw->r = r;
					uhw->ui = r->ui;
					uhw->uie = uie;
					q = uhw->ogp;
					while(lap < lape) {
						la = *lap++;
						q->ove = (q->ov = la->ov) + la->nnz;
						q->oc = la->oc;
						++q;
						}
					}
				else {
					si = s;
					while(lap < lape) {
						la = *lap++;
						lap1 = lap;
						ov = la->ov;
						oc = la->oc;
						nov = la->nnz;
						*si = 1;
						for(j = 0; j < nov; ++ j)
							V[ov[j]].dO = oc[j];
						pshv_prod_ASL(ew, r, nobj, ow, y);
						for(j = 0; j < nov; ++ j)
							V[ov[j]].dO = 0.;
						*si++ = 0;
						v = V + la->u.v;
						if ((t = v->aO))
							new_Hesoprod(ew, nov, ov, oc, nov, ov, oc, t);
						while(lap1 < lape) {
						    la = *lap1++;
						    v = V + la->u.v;
						    if ((t = v->aO)) {
							ov1 = la->ov;
							oc1 = la->oc;
							nov1 = la->nnz;
							new_Hesoprod(ew,nov,ov,oc,nov1,ov1,oc1,t);
							new_Hesoprod(ew,nov1,ov1,oc1,nov,ov,oc,t);
							}
						    }
						}
					}
				}
			while((uhw = uhwi)) {
				uhwi = uhwi->next;
				r = uhw->r;
				q = uhw->ogp;
				qe = q + r->n;
				si = s;
				do {
					if (q->ov < q->ove && *q->ov == i)
						*si = *q->oc;
					si++;
					} while(++q < qe);
				pshv_prod_ASL(ew, r, nobj, ow, y);
				lap = r->lap;
				lape = lap + r->n;
				do {
					la = *lap++;
					if ((t = V[la->u.v].aO)) {
						ov = la->ov;
						oc = la->oc;
						for(ov1 = ov + la->nnz; ov < ov1; ++ov, ++oc)
							if ((j = *ov) <= i)
								Hi[j] += t**oc;
						}
					}
					while(lap < lape);
				q = uhw->ogp;
				si = s;
				do {
					if (q->ov < q->ove && *q->ov == i) {
						*si = 0;
						++q->ov;
						++q->oc;
						}
					si++;
					} while(++q < qe);
				if ((ui = ++uhw->ui) < uhw->uie) {
					utodoj = utodo + *ui;
					uhw->next = *utodoj;
					*utodoj = uhw;
					}
				}
			}
#ifdef MULTIPLE_THREADS
 uhw_done:
#endif
		*rtodo++ = 0;	/* reset */
		hop1 = *otodoi;
		*otodoi++ = 0;
		while((hop = hop1)) {
			hop1 = hop->next;
			ov = hop->ov;
			ove = hop->ove;
			oc = hop->oc;
			oc1 = hop->roc;
			t = hop->coef * *oc1;
			while((j = *ov) <= i) {
				Hi[j] += t**oc;
				if (++ov >= ove)
					break;
				++oc;
				}
			if ((ov1 = hop->rv + 1) < hop->rve) {
				hop->rv = ov1;
				hop->roc = oc1 + 1;
				otodoj = otodo + *ov1;
				hop->next = *otodoj;
				*otodoj = hop;
				}
			}
		k = (int)(hcs[1] - hcs[0]);
		hcs++;
		if (uli)
			H += *uli++;
		if (vsc) {
			t = *vsc1++;
			while(--k >= 0) {
				j = (int)*hr++;
				*H++ = t * vsc0[j] * H0[j];
				H0[j] = 0;
				}
			}
		else {
			while(--k >= 0) {
				*H++ = H0[j = (int)*hr++];
				H0[j] = 0;
				}
			}
		}
#ifdef MULTIPLE_THREADS
	if (asl->P.hesevalth > 0 && (asl->P.sph_opts & 32) && asl->P.thused != maxit) {
		printf("threads for sphes = %d\n", maxit);
		fflush(stdout);
		}
	asl->P.thused = maxit;
#endif
	H = H00;
	if ((ulc = spi->ulcopy))
		for(uli = spi->ulcend; ulc < uli; ulc += 2)
			H[ulc[1]] = H[ulc[0]];
	else if ((ui = spi->uptolow))
		while((k = *++ui)) {
			t = H[j = *++ui];
			while(--k) {
				t1 = H[i = *++ui];
				H[i] = t;
				t = t1;
				}
			H[j] = t;
			}
	}

/* Variant of sphes that has a final nerror argument, working
   similarly to the final nerror argument to objval_(), etc. */

 void
sphese_ew_ASL(EvalWorkspace *ew, SputInfo **spi, real *H, int nobj, real *ow, real *y, fint *nerror)
{
	Jmp_buf **Jp, *Jsave, b;

	Jp = !nerror || *nerror >= 0 ? &ew->err_jmpw : &ew->err_jmpw1;
	Jsave = *Jp;
	*Jp = &b;
	*nerror = 0;
	if (setjmp(b.jb))
		*nerror = 1;
	else
		sphes_ew_ASL(ew, spi, H, nobj, ow, y);
	*Jp = Jsave;
	}
