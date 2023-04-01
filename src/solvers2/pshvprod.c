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

#include "jacpdim.h"
#include "opno2.h"

#ifdef __cplusplus
extern "C" {
#endif
#ifdef X64_bit_pointers
#define alignarg(x) x
#else
#define alignarg(x)
#endif

 static void
hv_fwd(int *o, real *w)
{
	Condptrs *cp;
	Eresult *L, *R, *r;
	Minmaxptrs *mmp;
	int *d, *de, k, *o1, **pop;
	plterm **pp;
	real dO, *rp;
	tfinfo **ptfi, *tfi;
	void **v;

 top:
	for(;;) {
		switch(*o) {

		case OPRET:
			return;

		case OP_GOTO:
		case OP_NEXTBLK:
		case OPGOTO2:
		case OPGOTOF:
		case OPGOTOF2:
		case OPGOTOF2n:
			if (!(o = *(int**)(o+1)))
				return;
			continue;

		case OPGOTOMM:
			r = (Eresult*)(rp = w + o[1]);
			mmp = *(Minmaxptrs**)&rp[4];
			if ((k = mmp->d) >= 0) {
				L = (Eresult*)(w + k);
				r->dO = L->dO;
				}
			else
				r->dO = 0;
			o += 2;
			break;

/*		case Hv_timesR:	*/
		case OPMULT01:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			r->dO = R->dO * L->O;
			o += 5;
			break;

/*		case Hv_binaryR:*/
		case OPDIV01:
		case nOPPOW01:
		case nOPREM01:
		case OP_atan201:
			r = (Eresult*)(w + o[2]);
			R = (Eresult*)(w + o[4]);
			r->dO = R->dO * r->dL;
			o += 5;
			break;

#ifdef X64_bit_pointers
		case OP_GOTOalign:
		case OP_NEXTBLKalign:
		case OPGOTO2align:
		case OPGOTOFalign:
		case OPGOTOF2align:
		case OPGOTOF2nalign:
			if (!(o = *(int**)(o+2)))
				return;
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
			r->dO = R->dO * r->dL;
			o = o1;
			break;

/*		case Hv_timesLR:*/
		case OPMULT2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			r->dO = L->dO*R->O + R->dO*L->O;
			o += 5;
			break;

/*		case Hv_binaryLR:*/
		case n_OPPOW2:
		case OP_atan22:
		case OPDIV2:
		case nOPREM2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			r->dO = L->dO*r->dL + R->dO*r->dR;
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
			r = (Eresult*)(w + o[2]);	/* r->dR holds a derpblock* */
			L = (Eresult*)(w + o[3]);	/* nOPLESS2 uses r->dL2 instead of r->dR */
			R = (Eresult*)(w + o[4]);
			r->dO = L->dO*r->dL + R->dO*r->dL2;
			o = o1;
			break;

/*		case Hv_timesL:	*/
		case OPMULT10:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			r->dO = L->dO * R->O;
			o += 5;
			break;

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
		case OP_log101:
		case OP_log1:
		case OP_sin1:
		case OP_sinh1:
		case OP_sqrt1:
		case OP_tan1:
		case OPtanh1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = L->dO * r->dL;
			o += 4;
			break;

		case OP_exp1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = L->dO * r->O;
			o += 4;
			break;

		case OPDIV10:
		case nOPPOW1i:
		case nOPPOW10:
		case OP_atan210:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = L->dO * r->dL;
			o += 5;
			break;

#ifdef X64_bit_pointers
		case OPLESS10align:
			o1 = o + 6 + 2*sizeof(derpblock*)/sizeof(int);
			goto more_OPLESS10;
#endif
		case nOPLESS10:
			o1 = o + 5 + 2*sizeof(derpblock*)/sizeof(int);
alignarg(more_OPLESS10:)
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = L->dO * r->dL;
			o = o1;
			break;

/*		case Hv_vararg:	*/
		case OPMINLIST1:
		case OPMAXLIST1:
#ifdef X64_bit_pointers
		case OPMINLISTalign:
		case OPMAXLISTalign:
#endif
			mmp = *(Minmaxptrs**)&w[o[2]+4];
			if (!(o = mmp->f))
				return;
			continue;

/*		case Hv_if:	*/
		case nOPIF1:
		case nOPIF11:
		case nOPIF12:
		case nOPIF13:
#ifdef X64_bit_pointers
		case OPIF1align:
		case OPIF11align:
		case OPIF12align:
		case OPIF13align:
#endif
			r = (Eresult*)(w + o[2]);
			r->dO = r->aO = r->adO = 0.;
			v = (void**)&w[o[4]];
			cp = (Condptrs*)v[1];
			o = cp->f;
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
			r->dO = r->dL * L->dO;
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
				dO += rp[d[1]] * L->dO;
				}
			r->dO = dO;
			o = (int*)&ptfi[1] + tfi->n;
			break;

/*		case Hv_negate:	*/
		case OPUMINUS1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = -L->dO;
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
		case OPPLUS10:
		case OPMINUS10:
		case nOPREM10:
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
			r->dO = -R->dO;
			o += 5;
			break;

/*		case Hv_minusLR: */
		case OPMINUS2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			r->dO = L->dO - R->dO;
			o += 5;
			break;

		case OPCOPY1:
		case OPCOPY1a:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = L->dO;
			o += 4;
			break;

#ifdef X64_bit_pointers
		case OPCPOW1align:
			rp = (real*)&o[5];
			goto moreOPCPOW1;
#endif
		case nOPCPOW1:
			rp = (real*)&o[4];
alignarg(moreOPCPOW1:)
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			r->dO = r->dL * L->dO;
			o = (int*)&rp[1];
			break;

		case OPVARREF:
			o += 3;
			goto top;

#ifdef X64_bit_pointers
		case OPGOTOBalign:
			++o;
#endif
		case OPGOTOB:
			pop = (int**)(o+1);
			o = (int*)&pop[1];
			goto top;

		default:/*DEBUG*/
			fprintf(Stderr, "bad *o = %d in hv_fwd\n", *o);
			exit(1);
			r = 0; /* not reached */
		  }
		r->aO = r->adO = 0.;
		}
	}

 static void
func_back(int *o, real *w, tfinfo **ptfi)
{
	Eresult *r;
	int *d1, *d2, *de, *doff, *fh;
	real aO, adO, *g, *h, t;
	tfinfo *tfi;

	tfi = *ptfi;
	r = (Eresult*)(w + o[2]);
	g = (real*)r + 4;
	h = g + tfi->nr;
	doff = tfi->doff;
	fh = tfi->fh;
	aO = r->aO;
	adO = r->adO;
	de = doff + 2*tfi->nd;
	for(d1 = doff; d1 < de; d1 += 2) {
		r = (Eresult*)(w + *d1);
		r->adO += (t = g[d1[1]]) * adO;
		r->aO += t * aO;
		t = adO*r->dO;
		for(d2 = doff; d2 < de; d2 += 2) {
			r = (Eresult*)(w + *d2);
			r->aO += t * h[*fh++];
			}
		}
	}

 static void
funnel_back(EvalWorkspace *ew, cexp *c, Varval *v, real t)
{
	Varval *V;
	hes_fun *hf;
	int *vp, *vp1, *vpe;
	ograd *og;
	real *g, *h, *W;
	real aO, adO;

	W = ew->w;
	V = (Varval*)W;
	aO = v->aO = t;
	adO = v->adO;
	hf = c->hfun;
	if ((og = hf->og)) {
		do {
			v = V + og->varno;
			v->adO += (t = og->coef) * adO;
			v->aO += t*aO;
			}
			while((og = og->next));
		return;
		}
	g = W + hf->grdhes;
	h = g + hf->nd;
	vp = hf->vp;
	vpe = vp + hf->n;
	do {
		v = V + *vp++;
		v->adO += (t = *g++) * adO;
		v->aO += t*aO;
		t = adO * v->dO;
		vp1 = hf->vp;
		do V[*vp1++].aO += t * *h++;
		   while(vp1 < vpe);
		}
		while(vp < vpe);
	}

 static void
hv_back(int *o, real *w, real aO0, real adO0)
{
	Condptrs *cp;
	Eresult *r, *L, *R;
	Minmaxptrs *mmp;
	int k, *o1, *oe;
	real adO, t1, t2;
	void **v;

	if (!aO0 && !adO0)
		return;
	r = (Eresult*)(w + o[2]);
	r->aO = aO0;
	r->adO = adO0;
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
		case nOPPOW01:
		case OP_atan201:
			r = (Eresult*)(w + o[2]);
			R = (Eresult*)(w + o[4]);
			R->adO += r->adO * r->dL;
			R->aO += r->aO * r->dL  +  r->adO * R->dO * r->dL2;
			break;

		case nOPLESS01:
			r = (Eresult*)(w + o[2]);
			R = (Eresult*)(w + o[4]);
			R->adO += r->adO * r->dL;
			R->aO += r->aO * r->dL;
			break;

/*		case Hv_binaryLR: */
		case n_OPPOW2:
		case OP_atan22:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			adO = r->adO;
			L->adO += adO * r->dL;
			R->adO += adO * r->dR;
			t1 = adO * L->dO;
			t2 = adO * R->dO;
			L->aO  += r->aO*r->dL + t1*r->dL2 + t2*r->dLR;
			R->aO  += r->aO*r->dR + t1*r->dLR + t2*r->dR2;
			break;

		case OPDIV2:
			/* dL2 is used for dR2, */
			/* since the true dL2 vanishes */
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			adO = r->adO;
			L->adO += adO * r->dL;
			R->adO += adO * r->dR;
			t1 = adO * L->dO;
			t2 = adO * R->dO;
			L->aO  += r->aO*r->dL + t2*r->dLR;
			R->aO  += r->aO*r->dR + t1*r->dLR + t2*r->dL2;
			break;

		case nOPREM2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			adO = r->adO;
			L->aO  += r->aO*r->dL;
			R->aO  += r->aO*r->dR;
			L->adO += adO * r->dL;
			R->adO += adO * r->dR;
			break;

		case nOPLESS2:
			r = (Eresult*)(w + o[2]);	/* r->dR holds a derpblock* */
			L = (Eresult*)(w + o[3]);	/* nOPLESS2 uses r->dL2 instead of r->dR */
			R = (Eresult*)(w + o[4]);
			adO = r->adO;
			L->aO  += r->aO*r->dL;
			R->aO  += r->aO*r->dL2;
			L->adO += adO * r->dL;
			R->adO += adO * r->dL2;
			break;

/*		case Hv_unary: */
		case n_OPABS1:
		case nOPREM01:
		case nOPREM10:
		case nOPLESS10:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->adO += r->adO * r->dL;
			L->aO  += r->aO * r->dL;
			break;

/*		case Hv_unary: */
		case nOPCPOW1:
		alignarg(case OPCPOW1align:)
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
		case OP_log101:
		case OP_log1:
		case OP_sin1:
		case OP_sinh1:
		case OP_sqrt1:
		case OP_tan1:
		case OPtanh1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->adO += r->adO * r->dL;
			L->aO  += r->aO * r->dL  +  r->adO * L->dO * r->dL2;
			break;

		case OP_exp1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->adO += r->adO * r->O;
			L->aO  += r->aO * r->O  +  r->adO * L->dO * r->O;
			break;

		case OP_2POW1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->adO += r->adO * r->dL;
			L->aO  += r->aO * r->dL  +  2. * r->adO * L->dO;
			break;

/*		case Hv_vararg:	*/
#ifdef X64_bit_pointers
		case OPMINLISTalign:
		case OPMAXLISTalign:
#endif
		case OPMINLIST1:
		case OPMAXLIST1:
			r = (Eresult*)(w + (k = o[2]));
			mmp = *(Minmaxptrs**)(w + k + 4);
			if ((k = mmp->d) >= 0) {
				L = (Eresult*)(w + k);
				L->aO += r->aO;
				L->adO += r->adO;
				}
			o = mmp->b;
			continue;

/*		case Hv_if:	*/
#ifdef X64_bit_pointers
		case OPIF1align:
		case OPIF11align:
		case OPIF12align:
		case OPIF13align:
#endif
		case nOPIF1:
		case nOPIF11:
		case nOPIF12:
		case nOPIF13:
			v = (void**)&w[o[4]];
			cp = (Condptrs*)v[1];
			o = cp->b;
			continue;

		case n_OPPLTERM1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->aO += r->dL * r->aO;
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

/*		case Hv_negate: */
		case OPUMINUS1:
			L = (Eresult*)(w + o[3]);
 neg_end:
			r = (Eresult*)(w + o[2]);
			L->aO -= r->aO;
			L->adO -= r->adO;
			break;

/*		case Hv_plusR: */
		case OPPLUS01:
			L = (Eresult*)(w + o[4]);
			goto plus_end;

/*		case Hv_plusL: */
		case OPPLUS10:
		case OPMINUS10:
			L = (Eresult*)(w + o[3]);
 plus_end:
			r = (Eresult*)(w + o[2]);
			L->aO += r->aO;
			L->adO += r->adO;
			break;

/*		case Hv_plusLR: */
		case OPPLUS2:
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

/*		case Hv_minusLR: */
		case OPMINUS2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			L->aO += t1 = r->aO;
			L->adO += t2 = r->adO;
			R->aO -= t1;
			R->adO -= t2;
			break;

/*		case Hv_timesR: */
		case OPMULT01:
			r = (Eresult*)(w + o[2]);
			R = (Eresult*)(w + o[4]);
			R->adO += r->adO * (t1 = w[o[3]]);
			R->aO += r->aO * t1;
			break;

/*		case Hv_timesL: */
		case OPDIV10:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->adO += r->adO * r->dL;
			L->aO  += r->aO * r->dL;
			break;

		case OPMULT10:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			L->adO += r->adO * R->O;
			L->aO  += r->aO * R->O;
			break;

/*		case Hv_timesLR: */
		case OPMULT2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			adO = r->adO;
			L->aO  += r->aO*R->O  +  adO * R->dO;
			R->aO  += r->aO*L->O  +  adO * L->dO;
			L->adO += adO * R->O;
			R->adO += adO * L->O;
			break;

		  case OPCOPY1:
		  case OPCOPY1a:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->adO += r->adO;
			L->aO += r->aO;
			break;

		  case OPVARREF:
		  case OP_PLTERM0align:
		  case OP_PLTERM1align:
			break;

		  default:/*DEBUG*/
			fprintf(Stderr, "bad *o = %d in hv_back\n", *o);
			exit(1);
		  }
		o -= o[1];
		}
	}

 static void
hv_fwd0(EvalWorkspace *ew, cexp *c, Varval *v)
{
	ASL_pfgh *asl;
	Varval *V, *v1;
	hes_fun *hf;
	int i, *o, *ob, *vp, *vpe;
	linarg *la;
	lincoef *lc, *lce;
	linpart *L;
	ograd *og;
	real *g, *w, x;

	w = ew->w;
	V = (Varval*)w;
	v->aO = v->adO = 0;
	if ((hf = c->hfun)) {
		x = 0;
		if ((og = hf->og))
			do x += og->coef * V[og->varno].dO;
			while((og = og->next));
		else {
			g = w + hf->grdhes;
			vp = hf->vp;
			vpe = vp + hf->n;
			do x += *g++ * V[*vp++].dO;
			   while(vp < vpe);
			}
		}
	else if ((o = c->o.f)) {
		hv_fwd(o, w);
		ob = c->o.b;
		v1 = (Varval*)(w + ob[2]);
		x = v1->dO;
		}
	else if ((o = c->o.e) && *o == OPRET && (i = o[1]) >= 0) {
		v1 = (Varval*)(w + i);
		x = v1->dO;
		}
	else
		x = 0.;
	if ((la = c->la)) {
		asl = (ASL_pfgh*)ew->asl;
		i = c - asl->I.cexps2_;
		x += asl->P.dv[i].scale * V[la->u.v].dO;
		}
	else if ((L = c->lp)) {
		lc = L->lc;
		lce = lc + L->n;
		for(lce = lc + L->n; lc < lce; ++lc)
			x += lc->coef * V[lc->varno].dO;
		}
	v->dO = x;
	}

 static void
hfg_fwd(Ops *O, real *w)
{
	Condptrs *cp;
	Eresult *r;
	Minmaxptrs *mmp;
	int *o, *o0, *oend, **pop;
	plterm **pp;
	real *rp;
	tfinfo **ptfi, *tfi;
	void **v;

	o = O->f;
	oend = O->b;
	for(;;) {
		o0 = o;
		switch(*o) {
/*		case Hv_timesL:	*/
		case OPMULT10:
		case OPDIV10:
		case nOPPOW1i:
/*		case Hv_timesR:	*/
		case OPMULT01:
/*		case Hv_binaryR:*/
		case OPDIV01:
		case nOPPOW01:
		case nOPLESS01:
		case OP_atan201:
/*		case Hv_binaryLR:*/
		case n_OPPOW2:
		case OP_atan22:
		case OPDIV2:
		case nOPREM2:
		case nOPLESS2:
/*		case Hv_timesLR:*/
		case OPMULT2:
/*		case Hv_plusR:	*/
		case OPPLUS01:
/*		case Hv_plusL:	*/
		case OPPLUS10:
		case OPMINUS10:
/*		case Hv_plusLR:	*/
		case OPPLUS2:
/*		case Hv_minusR:	*/
		case OPMINUS01:
/*		case Hv_minusLR: */
		case OPMINUS2:
			r = (Eresult*)&w[o[2]];
			o += 5;
			break;
		case nOPCPOW1:
			r = (Eresult*)&w[o[2]];
			rp = (real*)&o[4];
			o = (int*)&rp[1];
			break;

/*		case Hv_unary:	*/
		case n_OPABS1:
		case nOPLESS10:
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
		case OP_2POW1:
/*		case Hv_negate:	*/
		case OPUMINUS1:
		case OPCOPY1:
			r = (Eresult*)&w[o[2]];
			o += 4;
			break;

/*		case Hv_vararg:	*/
		case OPMINLIST1:
		case OPMAXLIST1:
#ifdef X64_bit_pointers
		case OPMINLISTalign:
		case OPMAXLISTalign:
#endif
			r = (Eresult*)&w[o[2]];
			mmp = *(Minmaxptrs**)&w[o[2]+4];
			o = mmp->f;
			break;

/*		case Hv_if:	*/
		case nOPIF1:
		case nOPIF11:
		case nOPIF12:
		case nOPIF13:
#ifdef X64_bit_pointers
		case OPIF1align:
		case OPIF11align:
		case OPIF12align:
		case OPIF13align:
#endif
			v = (void**)&w[o[4]];
			cp = (Condptrs*)v[1];
			o = cp->f;
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
			o = (int*)&pp[1];
			break;

/*		case Hv_sumlist:	*/
		case OPSUMLIST1:
			r = (Eresult*)(w + o[2]);
			o = o + o[3] + 4;
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
			r = (Eresult*)(w + o[2]);
			o = (int*)&ptfi[1] + tfi->n;
			break;

		case OP_GOTO:
		case OPGOTO2:
		case OP_NEXTBLK:
			o = *(int**)(o+1);
			continue;

		case OPGOTOF:
			o = *(int**)(o+1);
			continue;

		case OPGOTOF2:
		case OPGOTOF2n:
			pop = (int**)&o[1];
			o = (int*)&pop[2];
			continue;

		case OPGOTOMM:
			o += 2;
			continue;

#ifdef X64_bit_pointers
		case OP_GOTOalign:
		case OPGOTO2align:
		case OPGOTOFalign:
			o = *(int**)(o+2);
			continue;

		case OP_NEXTBLKalign:
			o = *(int**)(o+2);
			continue;

		case OPGOTOBalign:
			++o;

		case OPCPOW1align:
			r = (Eresult*)(w + o[2]);
			rp = (real*)&o[5];
			o = (int*)&rp[1];
			break;
#endif
		case OPGOTOB:
			pop = (int**)(o+1);
			o = (int*)&pop[1];
			continue;

		default:/*DEBUG*/
			fprintf(Stderr, "bad *o = %d in hfg_fwd\n", *o);
			exit(1);
			r = 0; /* not reached */
		}
		r->aO = 0.;
		if (o0 == oend)
			break;
		}
	}

 static void
hfg_back(Ops *O, real *w)
{
	Condptrs *cp;
	Eresult *r, *L, *R;
	Minmaxptrs *mmp;
	int k, *o, *o1, *oe, *oend;
	real t1;
	void **v;

	o = O->b;
	oend = O->f;
	r = (Eresult*)(w + o[2]);
	r->aO = 1.;
 top:
	for(;;) {
		switch(*o) {
		case OPRET:
			return;
		case OPGOTOB:
			o = *(int**)(o + 1);
			goto top;
#ifdef X64_bit_pointers
		case OPGOTOBalign:
			o = *(int**)(o + 2);
			goto top;
#endif
/*		case Hv_timesR: */
/*		case Hv_binaryR: */
		case OPDIV01:
		case nOPPOW01:
		case nOPLESS01:
		case OP_atan201:
			r = (Eresult*)(w + o[2]);
			R = (Eresult*)(w + o[4]);
			R->aO += r->aO * r->dL;
			break;

		case OPMULT01:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			R->aO += r->aO * L->O;
			break;

/*		case Hv_binaryLR: */
		case n_OPPOW2:
		case OP_atan22:
/*		case Hv_timesLR: */
		case OPDIV2:
		case nOPREM2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			L->aO  += r->aO*r->dL;
			R->aO  += r->aO*r->dR;
			break;

		case nOPLESS2:
			r = (Eresult*)(w + o[2]);	/* r->dR holds a derpblock* */
			L = (Eresult*)(w + o[3]);	/* nOPLESS2 uses r->dL2 instead of r->dR */
			R = (Eresult*)(w + o[4]);
			L->aO  += r->aO*r->dL;
			R->aO  += r->aO*r->dL2;
			break;

/*		case Hv_unary: */
/*		case Hv_timesL: */
		case n_OPABS1:
		case nOPREM10:
		case nOPLESS10:
		case nOPCPOW1:
		alignarg(case OPCPOW1align:)
		case nOPPOW10:
		case OP_acos1:
		case OP_acosh1:
		case OP_asin1:
		case OP_asinh1:
		case OP_atan1:
		case OP_atan210:
		case OP_atanh1:
		case OP_cos1:
		case OP_cosh1:
		case OP_log101:
		case OP_log1:
		case OP_sin1:
		case OP_sinh1:
		case OP_sqrt1:
		case OP_tan1:
		case OPtanh1:
		case OP_2POW1:
		case OPDIV10:
		case nOPPOW1i:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->aO  += r->aO * r->dL;
			break;

		case OP_exp1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->aO  += r->aO * r->O;
			break;

/*		case Hv_vararg:	*/
#ifdef X64_bit_pointers
		case OPMINLISTalign:
		case OPMAXLISTalign:
#endif
		case OPMINLIST1:
		case OPMAXLIST1:
			mmp = *(Minmaxptrs**)(w + o[2] + 4);
			if ((k = mmp->d) >= 0) {
				r = (Eresult*)(w + o[2]);
				L = (Eresult*)(w + k);
				L->aO += r->aO;
				}
			o = mmp->b;
			continue;

/*		case Hv_if:	*/
#ifdef X64_bit_pointers
		case OPIF1align:
		case OPIF11align:
		case OPIF12align:
		case OPIF13align:
#endif
		case nOPIF1:
		case nOPIF11:
		case nOPIF12:
		case nOPIF13:
			v = (void**)&w[o[4]];
			cp = (Condptrs*)v[1];
			if (cp->bder >= 0) {
				r = (Eresult*)(w + o[2]);
				L = (Eresult*)(w + cp->bder);
				L->aO += r->aO;
				}
			o = cp->b;
			continue;

		case n_OPPLTERM1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->aO += r->dL * r->aO;
			break;

		case OPSUMLIST1:
			r = (Eresult*)(w + o[2]);
			t1 = r->aO;
			o1 = o + 4;
			oe = o1 + o[3];
			while(o1 < oe) {
				if ((k = *o1++) >= 0) {
					L = (Eresult*)(w + k);
					L->aO += t1;
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

/*		case Hv_negate: */
		case OPUMINUS1:
			L = (Eresult*)(w + o[3]);
 neg_end:
			r = (Eresult*)(w + o[2]);
			L->aO -= r->aO;
			break;

/*		case Hv_plusR: */
		case OPPLUS01:
			L = (Eresult*)(w + o[4]);
			goto plus_end;

/*		case Hv_plusL: */
		case OPPLUS10:
		case OPMINUS10:
			L = (Eresult*)(w + o[3]);
 plus_end:
			r = (Eresult*)(w + o[2]);
			L->aO += r->aO;
			break;

/*		case Hv_plusLR: */
		case OPPLUS2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			L->aO += t1 = r->aO;
			R->aO += t1;
			break;

/*		case Hv_minusR: */
		case OPMINUS01:
			L = (Eresult*)(w + o[4]);
			goto neg_end;

/*		case Hv_minusLR: */
		case OPMINUS2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			L->aO += t1 = r->aO;
			R->aO -= t1;
			break;

		case OPMULT10:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			L->aO  += r->aO * R->O;
			break;

/*		case Hv_timesLR: */
		case OPMULT2:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			R = (Eresult*)(w + o[4]);
			L->aO  += r->aO*R->O;
			R->aO  += r->aO*L->O;
			break;

		  case OPCOPY1:
			r = (Eresult*)(w + o[2]);
			L = (Eresult*)(w + o[3]);
			L->aO += r->aO;
			break;

		  default:/*DEBUG*/
			fprintf(Stderr, "bad *o = %d in hvg_back\n", *o);
			exit(1);
		  }
		if (o == oend)
			break;
		o -= o[1];
		}
	}

 static void
funnelhes(EvalWorkspace *ew)
{
	ASL_pfgh *asl;
	Eresult *u;
	Varval *V, *v;
	cexp *c;
	hes_fun *hf;
	int *e, nd, *vp, *vp1, *vpe;
	real *g, *h, *w;

	asl = (ASL_pfgh*)ew->asl;
	ew->x0kind &= ~ASL_need_funnel;
	w = ew->w;
	V = (Varval*)w;
	for(hf = asl->I.hesthread; hf; hf = hf->hfthread) {
		if (hf->og)
			continue;
		g = w  + hf->grdhes;
		nd = hf->nd;
		h = g + nd;
		c = hf->c;
		vp = hf->vp;
		vpe = vp + nd;

		do V[*vp++].aO = 0.;
		   while(vp < vpe);

		hfg_fwd(&c->o, w);
		hfg_back(&c->o, w);

		vp = hf->vp;
		do {
			v = &V[*vp++];
			*g++ = v->aO;
			v->dO = v->aO = v->adO = 0;
			}
			while(vp < vpe);

		vp = hf->vp;
		vpe = vp + hf->n;
		do {
			v = &V[*vp++];
			v->dO = 1;
			if ((e = c->o.f)) {
				hv_fwd(e, w);
				hv_back(c->o.b, w, 0., 1.);
				}
			else if (*(e = c->o.e) != OPRET) {
				u = (Eresult*)&w[e[2]];
				u->aO = 0;
				u->adO = 1;
				}
			v->dO = 0;
			vp1 = hf->vp;
			do {
				v = &V[*vp1++];
				*h++ = v->aO;
				v->aO = v->adO = 0;
				}
				while(vp1 < vpe);
			}
			while(vp < vpe);
		}
	}

 static void
hvp0comp(EvalWorkspace *ew, real *hv, real *p, int nobj, real *ow, real *y)
	/* p = direction */
	/* y = Lagrange multipliers */
	/* hv = result */
{
	ASL_pfgh *asl;
	Varval *v, *vp, *x, *x0, *x1, *xe;
	cexp *c, *c1;
	int *dvsp0, *e, gm, i, i0, i1, j, n, nc, ndv, *ndvsp, no, noe, nov, *ov, *ove;
	linarg *la;
	lincoef *lc, *lce;
	linpart *L;
	ps_func *f, *f0;
	psb_elem *b, *be;
	psg_elem *g, *ge;
	real *cscale, g1, g2, *oc, *p1, t, t2, *W;

#ifdef IGNORE_BOGUS_WARNINGS
	c = 0;
	dvsp0 = ndvsp = 0;
	x1 = 0;
#endif

	asl = (ASL_pfgh*)ew->asl;
	W = ew->w;
	if (ew->x0kind & ASL_need_funnel)
		funnelhes(ew);
	if (nobj >= 0 && nobj < n_obj) {
		ow = ow ? ow + nobj : &edag_one_ASL;
		no = nobj;
		noe = no + 1;
		}
	else {
		no = noe = 0;
		if (ow)
			noe = n_obj;
		}
	x0 = (Varval*)W;
	vp = x0 + asl->i.defvar0;
	n = c_vars >= o_vars ? c_vars : o_vars;
	for(la = asl->P.lalist; la; la = la->lnext) {
		ov = la->ov;
		ove = ov + la->nnz;
		oc = la->oc;
		t = p[*ov]**oc++;
		while(++ov < ove)
			t += p[*ov]**oc++;
		x = x0 + la->u.v;
		x->dO = t;
		x->aO = x->adO = 0;
		}
	p1 = p;
	x = x0;
	for(xe = x + n; x < xe; x++) {
		x->dO = *p1++;
		x->aO = x->adO = 0.;
		}
	if ((ndv = asl->P.ncom)) {
		dvsp0 = asl->P.dvsp0;
		ndvsp = asl->P.ndvsp;
		x1 = (Varval*)ew->dv;
		c = cexps;
		for(i = 0; i < ndv; ++i) {
			if ((i1 = ndvsp[i])) {
				i1 += i0 = dvsp0[i];
				do hv_fwd0(ew, c + i0, &vp[i0]);
				   while(++i0 < i1);
				}
			hv_fwd0(ew, c + i, x1 + i);
			}
		}
	if (!y || (nc = n_con) <= 0)
		goto no_y;
	cscale = asl->i.lscale;
	f0 = asl->P.cps;
	for(i = 0; i < nc; ++i) {
	    if ((t2 = y[i])) {
		if (cscale)
			t2 *= cscale[i];
		f = f0 + i;
		for(b = f->pi.b, be = f->pi.be; b < be; b++) {
			if ((e = b->o.f)) {
				hv_fwd(e, W);
				hv_back(b->o.b, W, 0., t2);
				}
			else if (*(e = b->o.e) != OPRET) {
				x = (Varval*)&W[e[2]];
				x->aO = 0;
				x->adO = t2;
				}
			}
		for(g = f->g, ge = f->ge; g < ge; g++) {
			gm = g->gm;
			g1 = W[gm+1];
			for(b = g->pi.b, be = g->pi.be; b < be; b++) {
				if ((e = b->o.f)) {
					hv_fwd(e, W);
					hv_back(b->o.b, W, 0., t2*g1);
					}
				else if (*(e = b->o.e) != OPRET) {
					x = (Varval*)&W[e[2]];
					x->aO = 0;
					x->adO = t2*g1;
					}
				}
			if ((g2 = W[gm+2])) {
				oc = W + gm + 3;
				ov = g->ov;
				nov = g->nov;
				t = 0.;
				for(j = 0; j < nov; ++j)
					t += oc[j]*p[ov[j]];
				t *= t2*g2;
				for(j = 0; j < nov; ++j)
					x0[ov[j]].aO += t*oc[j];
				}
			}
		}
	    }
 no_y:
	for(; no < noe; no++) {
	    if ((t2 = *ow++)) {
		f = asl->P.ops + no;
		for(b = f->pi.b, be = f->pi.be; b < be; b++) {
			if ((e = b->o.f)) {
				hv_fwd(e, W);
				hv_back(b->o.b, W, 0., t2);
				}
			else if (*(e = b->o.e) != OPRET) {
				x = (Varval*)&W[e[2]];
				x->aO = 0;
				x->adO = t2;
				}
			}
		for(g = f->g, ge = f->ge; g < ge; g++) {
			gm = g->gm;
			g1 = W[gm+1];
			for(b = g->pi.b, be = g->pi.be; b < be; b++) {
				if ((e = b->o.f)) {
					hv_fwd(e, W);
					hv_back(b->o.b, W, 0., t2*g1);
					}
				else if (*(e = b->o.e) != OPRET) {
					x = (Varval*)&W[e[2]];
					x->aO = 0;
					x->adO = t2*g1;
					}
				}
			if ((g2 = W[gm+2])) {
				t = 0.;
				oc = W + gm + 3;
				ov = g->ov;
				nov = g->nov;
				for(j = 0; j < nov; ++j)
					t += oc[j]*p[ov[j]];
				t *= t2*g2;
				for(j = 0; j < nov; ++j)
					x0[ov[j]].aO += t*oc[j];
				}
			}
		}
	    }
	for(i = ndv; --i >= 0; ) {
		if ((i1 = ndvsp[i])) {
			i1 += i0 = dvsp0[i];
			do {
				v = &vp[--i1];
				c1 = c + i1;
				if ((t = v->aO) && (L = c1->lp))
					for(lc = L->lc, lce = lc + L->n; lc < lce; ++lc)
						x0[lc->varno].aO += t * lc->coef;
				if (c1->hfun)
					funnel_back(ew, c1, v, t);
				else if ((e = c1->o.b))
					hv_back(e, W, t, v->adO);
				else if (*(e = c1->o.e) != OPRET) {
					x = (Varval*)&W[e[2]];
					x->aO = t;
					x->adO = v->adO;
					}
				} while(i1 > i0);
			}
		x = x1 + i;
		c1 = c + i;
		if ((t = x->aO) && (L = c1->lp))
			for(lc = L->lc, lce = lc + L->n; lc < lce; ++lc)
				x0[lc->varno].aO += t * lc->coef;
		if (c1->hfun)
			funnel_back(ew, c1, x, t);
		else if ((e = c1->o.b))
			hv_back(e, W, t, x->adO);
		else if ((e = c1->o.e) && *e != OPRET) {
			v = (Varval*)&W[e[2]];
			v->aO = t;
			v->adO = x->adO;
			}
		}
	for(la = asl->P.lalist; la; la = la->lnext) {
		if ((t = x0[la->u.v].aO)) {
			ov = la->ov;
			ove = ov + la->nnz;
			oc = la->oc;
			do x0[*ov++].aO += t**oc++;
				while(ov < ove);
			}
		}
	x = x0;
	while(x < xe)
		*hv++ = (x++)->aO;
	}

 static real *	/* Compute vector x0 = mat(h)*y0,	*/
		/* where h = upper triang of mat(h).	*/
dtmul(int n, real *x0, real *h, real *y0)
{
	int i;
	real *hi, t, *x, *y, *y1, yi;

	y1 = y0;
	--h;
	for(i = 0; i < n; i++) {
		hi = ++h + i;
		yi = *y1++;
		t = yi**hi;
		x = x0;
		y = y0;
		while(h < hi) {
			t += *y++**h;
			*x++ += yi**h++;
			}
		*x = t;
		}
	return x0;
	}

 extern void hvpinit_nc_ASL(EvalWorkspace*, int, int, real*, real*);

 void
hvpcomp_ew_ASL(EvalWorkspace *ew, real *hv, real *p, int nobj, real *ow, real *y)
	/* p = direction */
	/* y = Lagrange multipliers */
	/* hv = result */
{
	ASL *a;
	ASL_pfgh *asl;
	Ihinfo *ihi;
	Varval *V, *v;
	int j, kp, kw, n, no, noe, nov, ns, nv, *ov, *ove, *ui, *uie;
	linarg *la, **lap, **lape;
	ps_func *ps, *pe;
	psg_elem *g, *ge;
	range *r;
	real *W, *cscale, *oc, *owi, t, t1, t2, *p0, *s, *w, *wi, *x;

	a = ew->asl;
	W = ew->w;
	V = (Varval*)W;
	ASL_CHECK(a, ASL_read_pfgh, "hvpcomp");
	if (ew->Derrs)
		deriv2_errchk_ASL(ew, 3);
	++ew->stats.hesvec;
	asl = (ASL_pfgh*)a;
	xpsg_check_ASL(ew, nobj, ow, y);
	nv = n_var;
	kp = htcl(nv*sizeof(real));
	p0 = 0;
	if ((s = asl->i.vscale)) {
		p0 = (real*)new_mblk(kp);
		for(n = 0; n < nv; n++)
			p0[n] = s[n] * p[n];
		p = p0;
		}
	if (!ew->ihdcur) {
		if (ew->ndhmax <= 0) {
			hvp0comp(ew,hv,p,nobj,ow,y);
			goto done;
			}
		if (!(n = ew->nhvprod))
			n = asl->P.ihdmin;
		if (n >= asl->P.ihdmin)
			hvpinit_nc_ASL(ew, ihd_limit, nobj, ow, y);
		}
	ew->nhvprod++;
	memset(hv, 0, nv*sizeof(real));
	for(la = asl->P.lalist; la; la = la->lnext) {
		ov = la->ov;
		ove = ov + la->nnz;
		oc = la->oc;
		t = p[*ov]**oc++;
		while(++ov < ove)
			t += p[*ov]**oc++;
		v = &V[la->u.v];
		v->dO = t;
		v->aO = v->adO = 0;
		}
	kw = kp + 1;
	w = (real*)new_mblk(kw);
	x = w + n_var;
	s = &W[asl->P.dOscratch];
	ns = 0;
	for(ihi = asl->P.ihi1; ihi && (r = ihi->r); ihi = ihi->next) {
		if (r->hest)
		    for(; r; r = r->rlist.prev) {
			n = r->n;
			nv = r->nv;
			wi = w;
			if (n < nv) {
				lap = r->lap;
				lape = lap + n;
				do {
					la = *lap++;
					ov = la->ov;
					ove = ov + la->nnz;
					oc = la->oc;
					t = p[*ov]**oc++;
					while(++ov < ove)
						t += p[*ov]**oc++;
					*wi++ = t;
					}
					while(lap < lape);
				wi = dtmul(n, x, W + r->hest, w);
				lap = r->lap;
				do if ((t = *wi++)) {
					la = *lap;
					ov = la->ov;
					ove = ov + la->nnz;
					oc = la->oc;
					do hv[*ov++] += t**oc++;
						while(ov < ove);
					}
					while(++lap < lape);
				}
			else {
				ui = r->ui;
				uie = ui + nv;
				do *wi++ = p[*ui++];
					while(ui < uie);
				wi = dtmul(nv, x, W + r->hest, w);
				ui = r->ui;
				do hv[*ui++] += *wi++;
					while(ui < uie);
				}
			}
		else
		    for(; r; r = r->rlist.prev) {
			n = r->n;
			if (ns < n)
				ns = n;
			wi = s;
			lap = r->lap;
			lape = lap + n;
			do {
				la = *lap++;
				ov = la->ov;
				ove = ov + la->nnz;
				oc = la->oc;
				t = p[*ov]**oc++;
				while(++ov < ove)
					t += p[*ov]**oc++;
				*wi++ = t;
				}
				while(lap < lape);
			pshv_prod_ASL(ew, r, nobj, ow, y);
			lap = r->lap;
			do {
				la = *lap++;
				if ((t = V[la->u.v].aO)) {
					ov = la->ov;
					ove = ov + la->nnz;
					oc = la->oc;
					do hv[*ov++] += t**oc++;
						while(ov < ove);
					}
				}
				while(lap < lape);
			}
		}
	del_mblk(w);
	wi = s + ns;
	while(wi > s)
		*--wi = 0.;
	if (asl->P.nobjgroups) {
	    if (nobj >= 0 && nobj < n_obj) {
		owi = ow ? ow + nobj : &edag_one_ASL;
		no = nobj;
		noe = no + 1;
		}
	    else {
		nobj = -1;
		no = noe = 0;
		if ((owi = ow))
			noe = n_obj;
		}
	    for(; no < noe; no++)
		if ((t = *owi++)) {
		    ps = asl->P.ops + no;
		    g = ps->g;
		    for(ge = ps->ge; g < ge; g++)
			if ((t2 = W[g->gm+2]) && (nov = g->nov) > 0) {
				oc = W + g->gm + 3;
				ov = g->ov;
				t1 = p[*ov] * oc[0];
				for(j = 1; j < nov; ++j)
					t1 += p[ov[j]] * oc[j];
				t2 *= t*t1;
				for(j = 0; j < nov; ++j)
					hv[ov[j]]  += t2 * oc[j];
				}
		}
	    }
	if (asl->P.ncongroups && y) {
		cscale = a->i.lscale;
		ps = asl->P.cps;
		for(pe = ps + n_con; ps < pe; ps++, y++)
		    if ((t = cscale ? *cscale++ * *y : *y))
			for(g = ps->g, ge = ps->ge; g < ge; g++)
			    if ((t2 = W[g->gm+2]) && (nov = g->nov) > 0) {
				ov = g->ov;
				oc = W + g->gm + 3;
				t1 = p[*ov] * oc[0];
				for(j = 1; j < nov; ++j)
					t1 += p[ov[j]] * oc[j];
				t2 *= t*t1;
				for(j = 0; j < nov; ++j)
					hv[ov[j]]  += t2 * oc[j];
				}
		}
	for(la = asl->P.lalist; la; la = la->lnext)
		if ((t = V[la->u.v].aO)) {
			ov = la->ov;
			ove = ov + la->nnz;
			oc = la->oc;
			do V[*ov++].aO += t**oc++;
				while(ov < ove);
			}
 done:
	if (p0) {
		del_mblk(p0);
		s = asl->i.vscale;
		w = hv + n_var;
		while(hv < w)
			*hv++ *= *s++;
		}
	}

/* Variant of hvpcomp_ew_ASL that has a final nerror argument, working
   similarly to the final nerror argument to objval_(), etc. */

 void
hvpcompe_ew_ASL(EvalWorkspace *ew, real *hv, real *p, int nobj, real *ow, real *y, fint *nerror)
{
	Jmp_buf **Jp, *Jsave, b;

	Jp = !nerror || *nerror >= 0 ? &ew->err_jmpw : &ew->err_jmpw1;
	Jsave = *Jp;
	*Jp = &b;
	*nerror = 0;
	if (setjmp(b.jb))
		*nerror = 1;
	else
		hvpcomp_ew_ASL(ew, hv, p, nobj, ow, y);
	*Jp = Jsave;
	}

 void
pshv_prod_ASL(EvalWorkspace *ew, range *r, int nobj, real *ow, real *y)
{
	ASL_pfgh *asl;
	Varval *V, *v, *vp, *x;
	cexp *c;
	int *cei, *cei0, *ceie, *e, i, *o;
	linarg *la, **lap, **lape;
	lincoef *lc, *lce;
	linpart *L;
	ps_func *p;
	psb_elem *b;
	psg_elem *g;
	real *cscale, *s, owi, t, *w;

	asl = (ASL_pfgh*)ew->asl;
	w = ew->w;
	V = (Varval*)w;
	vp = V + asl->i.defvar0;
	cscale = asl->i.lscale;
	owi = 1.;
	if (nobj >= 0 && nobj < n_obj) {
		if (ow) {
			if ((owi = ow[nobj]) == 0.)
				nobj = -1;
			ow = 0;
			}
		}
	if (ew->x0kind & ASL_need_funnel)
		funnelhes(ew);
	s = &w[asl->P.dOscratch];
	lap = r->lap;
	lape = lap + r->n;
	while(lap < lape) {
		v = &V[(*lap++)->u.v];
		v->dO = *s++;
		v->adO = v->aO = 0.;
		}
	if ((cei = cei0 = r->cei)) {
		i = *cei0++;
		ceie = (cei = cei0) + i;
		do {
			i = *cei++;
			hv_fwd0(ew, cexps + i, &vp[i]);
			}
			while(cei < ceie);
		}
	for(b = r->refs; b; b = b->next) {
		if ((i = b->conno) < 0) {
			i = -2 - i;
			if (i == nobj)
				t = owi;
			else if (ow) {
				if (!(t = ow[i]))
					continue;
				}
			else
				continue;
			p = asl->P.ops;
			}
		else {
			if (!y || !(t = y[i]))
				continue;
			if (cscale)
				t *= cscale[i];
			p = asl->P.cps;
			}
		if (b->groupno) {
			p += i;
			g = p->g + (b->groupno - 1);
			if (asl->P.pshv_g1)
				t *= w[g->gm+1];
			}
		if ((e = b->o.f)) {
			hv_fwd(e, w);
			hv_back(b->o.b, w, 0., t);
			}
		else if (*(e = b->o.e) != OPRET) {
			v = (Varval*)&w[e[2]];
			v->adO += t;
			}
		else {
			v = (Varval*)&w[e[1]];
			v->adO += t;
			}
		}
	while(cei > cei0) {
		i = *--cei;
		c = cexps + i;
		v = &vp[i];
		if ((t = v->aO) && (L = c->lp)) {
		    if ((la = c->la))
			V[la->u.v].aO += t * asl->P.dv[i].scale;
		    else {
			for(lc = L->lc, lce = lc + L->n; lc < lce; ++lc)
				V[lc->varno].aO += t * lc->coef;
			}
		    }
		if (c->hfun)
			funnel_back(ew, c, v, t);
		else if ((e = c->o.b))
			hv_back(e, w, t, v->adO);
		else if ((o = c->o.e) && *o != OPRET) {
			x = (Varval*)&w[e[2]];
			x->aO += t;
			x->adO += v->adO;
			}
		}
	}

 void
hvpcompd_ew_ASL(EvalWorkspace *ew, real *hv, real *p, int co)
	/* p = direction */
	/* hv = result */
	/* co >= 0: behave like hvpcomp_ASL with nobj = -1, ow = 0, y[i] = i == co ? 1. : 0. */
	/* co < 0: behave like hvpcomp_ASL with nobj = -1 - co, ow = 0, y = 0 */
{
	ASL *a;
	ASL_pfgh *asl;
	Varval *V, *v, *vp, *x, *x1;
	cexp *c, *c1;
	cgrad *cg, *cg0;
	int *dvsp0, *e, gm, i0, i1, j, kp, n, ndv, no, nov, *ndvsp, nx, *ov, *ove, oxk;
	linarg *la;
	lincoef *lc, *lce;
	linpart *L;
	ograd *og, *og0;
	ps_func *f;
	psb_elem *b, *be;
	psg_elem *g, *ge;
	real g1, g2, *oc, *s, t, t2, *p0, *w;
	varno_t i;

#ifdef IGNORE_BOGUS_WARNINGS
	c = 0;
	dvsp0 = ndvsp = 0;
	x1 = 0;
#endif
	asl = (ASL_pfgh*)(a = ew->asl);
	ASL_CHECK(a, ASL_read_pfgh, "hvpcompi");
	++ew->stats.hesvec;
	w = ew->w;
	V = (Varval*)w;
	vp = (Varval*)ew->dv;
	if (ew->x0kind & ASL_first_x) {
		if (!(s = X0))
			memset(s = ew->Lastx, 0, n_var*sizeof(real));
		xp_check_ASL(ew, s);
		ew->x0kind &= ~ASL_first_x;
		}
	nx = ew->nxval;
	oxk = ew->x0kind;
	ew->x0kind |= ASL_x_known;

	p0 = 0;
	cg0 = 0;
	og0 = 0;
	n = c_vars >= o_vars ? c_vars : o_vars;
	t2 = 1.;
	memset(hv, 0, n_var*sizeof(real));
	for(la = asl->P.lalist; la; la = la->lnext) {
		ov = la->ov;
		ove = ov + la->nnz;
		oc = la->oc;
		t = p[*ov]**oc++;
		while(++ov < ove)
			t += p[*ov]**oc++;
		x = &V[la->u.v];
		x->dO = t;
		x->aO = x->adO = 0;
		}
	if (co >= 0) {
		if (co >= nlc)
			return;
		f = asl->P.cps + co;
		if (ew->ncxval[co] != nx)
			conpival_ew_ASL(ew, co, ew->Lastx, 0);
		if (f->g && ew->ncxval[co] != nx)
			conpgrd_ew_ASL(ew, co, ew->Lastx, 0, 0);
		if ((s = asl->i.lscale))
			t2 = s[co];
		cg = cg0 = Cgrad[co];
		if ((s = asl->i.vscale)) {
			kp = htcl(n*sizeof(real));
			p0 = (real*)new_mblk(kp);
			for(; cg; cg = cg->next) {
				i = cg->varno;
				x = V + i;
				x->dO = p0[i] = p[i]*s[i];
				x->aO = x->adO = 0.;
				}
			p = p0;
			}
		else {
			for(; cg; cg = cg->next) {
				i = cg->varno;
				x = V + i;
				x->dO = p[i];
				x->aO = x->adO = 0.;
				}
			}
		}
	else {
		no = -1 - co;
		if (no >= nlo)
			return;
		f = asl->P.ops + no;
		if (ew->ncxval[no] != nx)
			objpval_ew_ASL(ew, no, ew->Lastx, 0);
		if (f->g && ew->noxval[no] != nx)
			objpgrd_ew_ASL(ew, no, ew->Lastx, 0, 0);
		og = og0 = Ograd[no];
		if ((s = asl->i.vscale)) {
			kp = htcl(n*sizeof(real));
			p0 = (real*)new_mblk(kp);
			for(; og; og = og->next) {
				i = og->varno;
				x = V + i;
				x->dO = p0[i] = p[i]*s[i];
				x->aO = x->adO = 0.;
				}
			p = p0;
			}
		else {
			for(; og; og = og->next) {
				i = og->varno;
				x = V + i;
				x->dO = p[i];
				x->aO = x->adO = 0.;
				}
			}
		}
	if (ew->Derrs) {
		ew->x0kind = oxk;
		deriv_errchk_ASL(ew, co, 1, 3);
		ew->x0kind |= ASL_x_known;
		}
	if ((ndv = asl->P.ncom)) {
		dvsp0 = asl->P.dvsp0;
		ndvsp = asl->P.ndvsp;
		x1 = (Varval*)ew->dv;
		c = cexps;
		for(j = 0; j < ndv; ++j) {
			if ((i1 = ndvsp[j])) {
				i1 += (i0 = dvsp0[j]);
				do hv_fwd0(ew, c + i0, &vp[i0]);
				   while(++i0 < i1);
				}
			hv_fwd0(ew, c + j, vp + j);
			}
		}
	for(b = f->pi.b, be = f->pi.be; b < be; b++) {
		if ((e = b->o.f)) {
			hv_fwd(e, w);
			hv_back(b->o.b, w, 0., t2);
			}
		else if (*(e = b->o.e) != OPRET) {
			x = (Varval*)&w[e[2]];
			x->aO = 0;
			x->adO = t2;
			}
		}
	for(g = f->g, ge = f->ge; g < ge; g++) {
		gm = g->gm;
		g1 = w[gm + 1];
		for(b = g->pi.b, be = g->pi.be; b < be; b++) {
			if ((e = b->o.f)) {
				hv_fwd(e, w);
				hv_back(b->o.b, w, 0., t2*g1);
				}
			else if (*(e = b->o.e) != OPRET) {
				x = (Varval*)&w[e[2]];
				x->aO = 0;
				x->adO = t2*g1;
				}
			}
		if ((g2 = w[gm+2]) && (nov = g->nov) > 0) {
			ov = g->ov;
			oc = w + gm + 3;
			t = oc[0] * p[*ov];
			for(j = 1; j < nov; ++j)
				t  += oc[j] * p[ov[j]];
			t *= t2*g2;
			for(j = 0; j < nov; ++j)
				V[ov[j]].aO += t*oc[j];
			}
		}
	for(j = ndv; --j >= 0; ) {
		if ((i1 = ndvsp[j])) {
			i1 += i0 = dvsp0[j];
			do {
				v = &vp[--i1];
				c1 = c + i1;
				if ((t = v->aO) && (L = c1->lp))
					for(lc = L->lc, lce = lc + L->n; lc < lce; ++lc)
						V[lc->varno].aO += t * lc->coef;
				if (c1->hfun)
					funnel_back(ew, c1, v, t);
				else if ((e = c1->o.b))
					hv_back(e, w, t, v->adO);
				else if (*(e = c1->o.e) != OPRET) {
					x = (Varval*)&w[e[2]];
					x->aO = t;
					x->adO = v->adO;
					}
				} while(i1 > i0);
			}
		c1 = c + j;
		x = x1 + j;
		if ((t = x->aO) && (L = c1->lp))
			for(lc = L->lc, lce = lc + L->n; lc < lce; ++lc)
				V[lc->varno].aO += t * lc->coef;
		x = x1 + j;
		if (c1->hfun)
			funnel_back(ew, c1, x, t);
		else if ((e = c1->o.b))
			hv_back(e, w, t, x->adO);
		else if (*(e = c1->o.e) != OPRET) {
			x = (Varval*)&w[e[2]];
			x->aO = t;
			x->adO = x->adO;
			}
		}
	for(la = asl->P.lalist; la; la = la->lnext)
		if ((t = V[la->u.v].aO)) {
			ov = la->ov;
			ove = ov + la->nnz;
			oc = la->oc;
			do V[*ov++].aO += t**oc++;
				while(ov < ove);
			}
	if ((cg = cg0)) {
		if (s) {
			while(cg) {
				i = cg->varno;
				hv[i] = s[i]*V[i].aO;
				cg = cg->next;
				}
			}
		else {
			while(cg) {
				i = cg->varno;
				hv[i] = V[i].aO;
				cg = cg->next;
				}
			}
		}
	else {
		og = og0;
		if (s) {
			while(og) {
				i = og->varno;
				hv[i] = s[i]*V[i].aO;
				og = og->next;
				}
			}
		else {
			while(og) {
				i = og->varno;
				hv[i] = V[i].aO;
				og = og->next;
				}
			}
		}
	if (p0)
		del_mblk(p0);
	}

/* Variant of hvpcompd_ew_ASL that has a final nerror argument, working
   similarly to the final nerror argument to objval_(), etc. */

 void
hvpcompde_ew_ASL(EvalWorkspace *ew, real *hv, real *p, int co, fint *nerror)
{
	Jmp_buf **Jp, *Jsave, b;

	Jp = !nerror || *nerror >= 0 ? &ew->err_jmpw : &ew->err_jmpw1;
	Jsave = *Jp;
	*Jp = &b;
	*nerror = 0;
	if (setjmp(b.jb))
		*nerror = 1;
	else
		hvpcompd_ew_ASL(ew, hv, p, co);
	*Jp = Jsave;
	}

 varno_t
hvpcomps_ew_ASL(EvalWorkspace *ew, real *hv, real *p, int co, varno_t nz, varno_t *z)
	/* p = direction */
	/* hv = result */
	/* co >= 0: behave like hvpcomp_ASL with nobj = -1, ow = 0, y[i] = i == co ? 1. : 0. */
	/* co < 0: behave like hvpcomp_ASL with nobj = -1 - co, ow = 0, y = 0 */
	/* Indices of up to nz nonzeros of the Hessian-vector product are stored in z. */
	/* The number of such nonzeros is returned (even if > nz). */
{
	ASL *a;
	ASL_pfgh *asl;
	Varval *V, *v, *vp, *x, *x1;
	cexp *c, *c1;
	cgrad *cg, *cg0;
	int *dvsp0, *e, gm, i0, i1, j, kp, n, ndv, *ndvsp, no, nov, nx, *ov, *ove, oxk;
	linarg *la;
	lincoef *lc, *lce;
	linpart *L;
	ograd *og, *og0;
	ps_func *f;
	psb_elem *b, *be;
	psg_elem *g, *ge;
	real g1, g2, *hve, *oc, *p0, *s, t, t2, *vscale, *w;
	varno_t F, i, rv, *ze;

#ifdef IGNORE_BOGUS_WARNINGS
	c = 0;
	dvsp0 = ndvsp = 0;
#endif
	asl = (ASL_pfgh*)(a = ew->asl);
	ASL_CHECK(a, ASL_read_pfgh, "hvpcompi");
	++ew->stats.hesvec;
	w = ew->w;
	V = (Varval*)w;
	vp = (Varval*)ew->dv;
	if (ew->x0kind & ASL_first_x) {
		if (!(s = X0))
			memset(s = ew->Lastx, 0, n_var*sizeof(real));
		xp_check_ASL(ew, s);
		ew->x0kind &= ~ASL_first_x;
		}
	nx = ew->nxval;
	oxk = ew->x0kind;
	ew->x0kind |= ASL_x_known;

	p0 = 0;
	cg0 = 0;
	og0 = 0;
	n = c_vars >= o_vars ? c_vars : o_vars;
	t2 = 1.;
	memset(hv, 0, n_var*sizeof(real));
	for(la = asl->P.lalist; la; la = la->lnext) {
		ov = la->ov;
		ove = ov + la->nnz;
		oc = la->oc;
		t = p[*ov]**oc++;
		while(++ov < ove)
			t += p[*ov]**oc++;
		x = &V[la->u.v];
		x->dO = t;
		x->aO = x->adO = 0;
		}
	no = -1 - co;
	if (co >= 0) {
		if (co >= nlc)
			return 0;
		f = asl->P.cps + co;
		if (ew->ncxval[co] != nx)
			conpival_ew_ASL(ew, co, ew->Lastx, 0);
		if (f->g && ew->ncxval[co] != nx)
			conpgrd_ew_ASL(ew, co, ew->Lastx, 0, 0);
		if ((s = asl->i.lscale))
			t2 = s[co];
		cg = cg0 = Cgrad[co];
		if ((vscale = asl->i.vscale)) {
			kp = htcl(n*sizeof(real));
			p0 = (real*)new_mblk(kp);
			for(; cg; cg = cg->next) {
				i = cg->varno;
				x = V + i;
				x->dO = p0[i] = p[i]*vscale[i];
				x->aO = x->adO = 0.;
				}
			p = p0;
			}
		else {
			for(; cg; cg = cg->next) {
				i = cg->varno;
				x = V + i;
				x->dO = p[i];
				x->aO = x->adO = 0.;
				}
			}
		}
	else {
		if (no >= nlo)
			return 0;
		f = asl->P.ops + no;
		if (ew->ncxval[no] != nx)
			objpval_ew_ASL(ew, no, ew->Lastx, 0);
		if (f->g && ew->noxval[no] != nx)
			objpgrd_ew_ASL(ew, no, ew->Lastx, 0, 0);
		og = og0 = Ograd[no];
		if ((vscale = asl->i.vscale)) {
			kp = htcl(n*sizeof(real));
			p0 = (real*)new_mblk(kp);
			for(; og; og = og->next) {
				i = og->varno;
				x = V + i;
				x->dO = p0[i] = p[i]*vscale[i];
				x->aO = x->adO = 0.;
				}
			p = p0;
			}
		else {
			for(; og; og = og->next) {
				i = og->varno;
				x = V + i;
				x->dO = p[i];
				x->aO = x->adO = 0.;
				}
			}
		}
	if (ew->Derrs) {
		ew->x0kind = oxk;
		deriv_errchk_ASL(ew, co, 1, 3);
		ew->x0kind |= ASL_x_known;
		}
	if ((ndv = asl->P.ncom)) {
		dvsp0 = asl->P.dvsp0;
		ndvsp = asl->P.dvsp0;
		c = cexps;
		for(j = 0; j < ndv; ++j) {
			if ((i1 = ndvsp[j])) {
				i1 += i0 = dvsp0[j];
				do hv_fwd0(ew, c + i0, &vp[i0]);
				   while(++i0 < i1);
				}
			hv_fwd0(ew, c + j, vp + j);
			}
		}
	for(b = f->pi.b, be = f->pi.be; b < be; b++) {
		if ((e = b->o.f)) {
			hv_fwd(e, w);
			hv_back(b->o.b, w, 0., t2);
			}
		else if (*(e = b->o.e) != OPRET) {
			x = (Varval*)&w[e[2]];
			x->aO = 0.;
			x->adO = t2;
			}
		}
	for(g = f->g, ge = f->ge; g < ge; g++) {
		gm = g->gm;
		g1 = w[gm + 1];
		for(b = g->pi.b, be = g->pi.be; b < be; b++) {
			if ((e = b->o.f)) {
				hv_fwd(e, w);
				hv_back(b->o.b, w, 0., t2*g1);
				}
			else if (*(e = b->o.e) != OPRET) {
				x = (Varval*)&w[e[2]];
				x->aO = 0.;
				x->adO = t2*g1;
				}
			}
		if ((g2 = w[gm+2]) && (nov = g->nov) > 0) {
			ov = g->ov;
			oc = w + gm + 3;
			t = oc[0] * p[*ov];
			for(j = 1; j < nov; ++j)
				t += oc[j] * p[ov[j]];
			t *= t2*g2;
			for(j = 0; j < nov; ++j)
				V[ov[j]].aO += t*oc[j];
			}
		}
	for(j = ndv; --j >= 0; ) {
		if ((i1 = ndvsp[j])) {
			i1 += i0 = dvsp0[j];
			do {
				v = &vp[--i1];
				c1 = c + i1;
				if ((t = v->aO) && (L = c1->lp))
					for(lc = L->lc, lce = lc + L->n; lc < lce; ++lc)
						V[lc->varno].aO += t * lc->coef;
				if (c1->hfun)
					funnel_back(ew, c1, v, t);
				else if ((e = c1->o.b))
					hv_back(e, w, t, v->adO);
				else if (*(e = c1->o.e) != OPRET) {
					x = (Varval*)&w[e[2]];
					x->aO = t;
					x->adO = v->adO;
					}
				} while(i1 > i0);
			}
		c1 = c + j;
		x1 = vp + j;
		if ((t = x1->aO) && (L = c1->lp))
			for(lc = L->lc, lce = lc + L->n; lc < lce; ++lc)
				V[lc->varno].aO += t * lc->coef;
		if (c1->hfun)
			funnel_back(ew, c1, x1, t);
		else if ((e = c1->o.b))
			hv_back(e, w, t, x1->adO);
		else if (*(e = c1->o.e) != OPRET) {
			x = (Varval*)&w[e[2]];
			x->aO = t;
			x->adO = x1->adO;
			}
		}
	for(la = asl->P.lalist; la; la = la->lnext) {
		if ((t = V[la->u.v].aO)) {
			ov = la->ov;
			ove = ov + la->nnz;
			oc = la->oc;
			do V[*ov++].aO += t**oc++;
				while(ov < ove);
			}
		}
	rv = 0;
	if ((ze = z))
		ze += nz;
	if ((hve = hv))
		hve += nz;
	F = Fortran;
	if ((cg = cg0)) {
		if (!hv) {
			while(cg) {
				++rv;
				if (z < ze)
					*z++ = cg->varno;
				cg = cg->next;
				}
			}
		else if (vscale) {
			while(cg) {
				++rv;
				i = cg->varno;
				if (z < ze)
					*z++ = F + i;
				if (hv < hve)
					*hv++ = vscale[i]*V[i].aO;
				cg = cg->next;
				}
			}
		else {
			while(cg) {
				++rv;
				i = cg->varno;
				if (z < ze)
					*z++ = F + i;
				if (hv < hve)
					*hv++ = V[i].aO;
				cg = cg->next;
				}
			}
		}
	else {
		og = Ograd[no];
		if (!hv) {
			while(og) {
				++rv;
				if (z < ze)
					*z++ = og->varno;
				og = og->next;
				}
			}
		else if (vscale) {
			while(og) {
				++rv;
				i = og->varno;
				if (z < ze)
					*z++ = F + i;
				if (hv < hve)
					*hv++ = vscale[i]*V[i].aO;
				og = og->next;
				}
			}
		else {
			while(og) {
				++rv;
				i = og->varno;
				if (z < ze)
					*z++ = F + i;
				if (hv < hve)
					*hv++ = V[i].aO;
				og = og->next;
				}
			}
		}
	if (p0)
		del_mblk(p0);
	return rv;
	}

/* Variant of hvpcomps_ew_ASL that has a final nerror argument, working
   similarly to the final nerror argument to objval_(), etc. */

 varno_t
hvpcompse_ew_ASL(EvalWorkspace *ew, real *hv, real *p, int co, varno_t nz, varno_t *z, fint *nerror)
{
	Jmp_buf **Jp, *Jsave, b;
	varno_t rv;

	Jp = !nerror || *nerror >= 0 ? &ew->err_jmpw : &ew->err_jmpw1;
	Jsave = *Jp;
	*Jp = &b;
	*nerror = 0;
	rv = 0;
	if (setjmp(b.jb))
		*nerror = 1;
	else
		rv = hvpcomps_ew_ASL(ew, hv, p, co, nz, z);
	*Jp = Jsave;
	return rv;
	}

#ifdef __cplusplus
}
#endif
