/*******************************************************************
Copyright (C) 2016, 2018, 2019 AMPL Optimization, Inc.; written by David M. Gay.

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

#include "asl.h"
#define SKIP_NL2_DEFINES
#include "nlp.h"
#include "nlp2.h"
#include "asl_pfg.h"
#include "asl_pfgh.h"
#undef cde
#include "obj_adj.h"
#include "opno.hd"
#include "opno2.h"

#ifndef LCADJ_GULP
#define LCADJ_GULP 1023
#endif
 enum { Gulp = LCADJ_GULP };

 enum { Lintype = 1, Contype = 2, BTtype = 3 };

 typedef struct
head {
	int type;
	int n;
	} head;

 typedef struct
Linform {
	head h;
	ograd *og;
	} Linform;

 typedef struct
Lconstr {
	head h; /* h.n gives constraint sense: 0 for <= , 1 for >= ,h 2 for == */
	real rhs;
	ograd *og;
	struct Lconstr *next;	/* for constraint && constraint && ... */
	} Lconstr;

 typedef struct
Btest {
	head h;	/* h.type = BTtype, h.n = varno */
	int bval;	/* 0 or 1 */
	Lconstr *lc[2];	/* "true" (==>) constraints and "false" (else) constraints */
	} Btest;

 typedef union Tchunk { union Tchunk *prev; real r; } Tchunk;

 typedef struct
LCADJ_Info {
	ASL *asl;
	Btest *bt;
	Lconstr **plc;
	Linform *lffree;
	head **w, **w0;
	ograd *freeog, **s;
	int *iego, *z, *zc;
	real *numvals, *r, *tfree, *tfree0;
	Tchunk *tchunks;	/* temp. memory (unlikly to be needed; freed after use) */
	Tchunk *tchunks0;	/* for defined var results */
	int nlv[6];
	int mxv;	/* asl->i.maxvar */
	int nnum;	/* asl->numlen / sizeof(real) = length of numvals array */
	int nv;		/* asl->i.n_var_ */
	int ntfree;	/* elements remaining in tfree */
	int n1og;	/* numbers of reals in one ograd */
	void *v;
	int *errinfo;
	Add_Indicator add_indic;
	head *(*Leval)(struct LCADJ_Info*, int*);
	} LCADJ_Info;

 static void
new_tchunk(LCADJ_Info *lci, int tneed)
{
	ASL *asl;
	Tchunk *tc;
	int k;

	asl = lci->asl;
	if (tneed < Gulp)
		tneed = Gulp;
	k = htcl(tneed*sizeof(real));
	tneed = ((sizeof(void*)<<k) / sizeof(Tchunk)) - 1;
	tc = (Tchunk*)new_mblk(k);
	tc->prev = lci->tchunks;
	lci->tchunks = tc;
	lci->tfree = &tc[1].r;
	lci->ntfree = tneed;
	}

 static real*
tmem(LCADJ_Info *lci, size_t L)
{
	int n = (L + sizeof(real) - 1) / sizeof(real);
	real *r;

	if (n > lci->ntfree)
		new_tchunk(lci, n);
	r = lci->tfree;
	lci->tfree = r + n;
	lci->ntfree -= n;
	return r;
	}

 static ograd *
new_og(LCADJ_Info *lci, int varno, real coef)
{
	int n1;
	ograd *og;

	if ((og = lci->freeog))
		lci->freeog = og->next;
	else {
		if (lci->ntfree < (n1 = lci->n1og))
			new_tchunk(lci, n1);
		og = (ograd*)lci->tfree;
		lci->tfree += n1;
		lci->ntfree -= n1;
		}
	og->next = 0;
	og->varno = varno;
	og->coef = coef;
	return og;
	}

 static ograd *
ogcopy(LCADJ_Info *lci, ograd *og)
{
	ograd *ogc, *ogp, *og1;

	if (!og)
		return 0;
	ogp = og1 = new_og(lci, og->varno, og->coef);
	while((og = og->next)) {
		ogp->next = ogc = new_og(lci, og->varno, og->coef);
		ogp = ogc;
		}
	return og1;
	}

 ograd *
getog(LCADJ_Info *lci, int i)
{
	Linform *lf;
	ograd *og;

	if (i < 0) {
		if ((i += lci->nnum) < 0)
			return 0;
		return new_og(lci, -1, lci->numvals[i]);
		}
	if (i < lci->nv)
		return new_og(lci, i, 1.);
	lf = (Linform*)lci->w[i];
	if (!lf || lf->h.type != Lintype)
		return 0;
	og = lf->og;
	if (i < lci->mxv)
		og = ogcopy(lci, og);
	else {
		lf->og = 0;
		*(Linform**)lf = lci->lffree;
		lci->lffree = lf;
		}
	return og;
	}

 ograd *
getog2(LCADJ_Info *lci, int i)
{
	Linform *lf;
	int i4;
	ograd *og;

	if (i < 0) {
		if ((i += lci->nnum) < 0)
			return 0;
		return new_og(lci, -1, lci->numvals[i]);
		}
	i4 = i >> 2;
	if (i4 < lci->nv)
		return new_og(lci, i4, 1.);
	lf = (Linform*)lci->w[i];
	if (!lf || lf->h.type != Lintype)
		return 0;
	og = lf->og;
	if (i < lci->mxv)
		og = ogcopy(lci, og);
	else {
		lf->og = 0;
		*(Linform**)lf = lci->lffree;
		lci->lffree = lf;
		}
	return og;
	}

 static Linform*
getlf(LCADJ_Info *lci, ograd *og)
{
	Linform *lf;

	if ((lf = lci->lffree))
		lci->lffree = *(Linform**)lf;
	else
		lf = (Linform*)tmem(lci, sizeof(Linform));
	lf->h.type = Lintype;
	lf->og = og;
	return lf;
	}

 static ograd *
og_free1(LCADJ_Info *lci, ograd *og)
{
	ograd *rv = og->next;
	og->next = lci->freeog;
	lci->freeog = og;
	return rv;
	}

 static Btest *
B_test(LCADJ_Info *lci, int t, int i, int j)
{
	ASL *asl;
	Btest *rv;
	int b;
	real *L, *U, rhs;

	if (i < 0) {
		b = j;
		j = i;
		i = b;
		if (t != nOPNE)
			t = 2*nOPEQ - t;
		}
	else if (j >= 0)
		return 0;
	if ((j += lci->nnum) < 0)
		return 0;
	if (i < lci->nlv[0] || (i >= lci->nlv[1] && i < lci->nlv[2])
	 || (i >= lci->nlv[3] && i < lci->nlv[4]) || i >= lci->nlv[5])
		return 0;
	asl = lci->asl;
	L = LUv;
	if ((U = Uvx)) {
		L += i;
		U += i;
		}
	else {
		L += 2*i;
		U = L + 1;
		}
	if (L[0] <= -1. || U[0] >= 2.)
		return 0;
	rhs = lci->numvals[j];
	switch(t) {
	  case nOPLE:
		if (rhs >= 1. || rhs < 0.)
			return 0;
		b = 0;
		break;
	  case nOPEQ:
		if (rhs == 0.)
			b = 0;
		else if (rhs == 1.)
			b = 1;
		else
			return 0;
		break;
	  case nOPGE:
		if (rhs <= 0. || rhs > 1.)
			return 0;
		b = 1;
		break;
	  case nOPNE:
		if (rhs == 0.)
			b = 1;
		else if (rhs == 1.)
			b = 0;
		else
			return 0;
		break;
	  default:
		return 0;
	  }
	rv = (Btest*)tmem(lci, sizeof(Btest));
	rv->h.type = BTtype;
	rv->h.n = i;
	rv->bval = b;
	rv->lc[0] = rv->lc[1] = 0;
	return rv;
	}

 static Btest *
B_test2(LCADJ_Info *lci, int t, int i, int j)
{
	switch(t) {
	  case n_OPLE:
		t = nOPLE;
		break;
	  case OPEQ:
		t = nOPEQ;
		break;
	  case n_OPGE:
		t = nOPGE;
		break;
	  case n_OPNE:
		t = nOPNE;
		break;
	  default:
		return 0;
	  }
	if (i > 0)
		i >>= 2;
	if (j > 0)
		j >>= 2;
	return B_test(lci, t, i, j);
	}

 static ograd *
ogsum(LCADJ_Info *lci, ograd *og1, ograd *og2)
{
	ograd **p, *rv;
	ssize_t d;

	rv = 0;	/* not necessary */
	p = &rv;
	for(;;) {
		if ((d = og1->varno - og2->varno) < 0) {
			*p = og1;
			if (!(og1 = *(p = &og1->next))) {
				*p = og2;
				break;
				}
			}
		else if (d > 0) {
			*p = og2;
			if (!(og2 = *(p = &og2->next))) {
				*p = og1;
				break;
				}
			}
		else {
			if ((og1->coef += og2->coef)) {
				*p = og1;
				og1 = *(p = &og1->next);
				}
			else
				og1 = og_free1(lci, og1);
			if (!(og2 = og_free1(lci, og2))) {
				*p = og1;
				break;
				}
			if (!og1) {
				*p = og2;
				break;
				}
			}
		}
	return rv;
	}

 static int
intcomp(const void *a, const void *b, void *v)
{
	return *(int*)a - *(int*)b;
	}

 static head *
Leval1(LCADJ_Info *lci, int *o)
{
	ASL *asl;
	Btest *bt;
	Lconstr *lc;
	Linform *lf, *lfd;
	head *rv, **w;
	int i, j, k, m, neg, nnum, nv, nz, **pop, *z, *zc;
	ograd *og, *og1, *og2, **pog;
	real *numvals, *r, t;

	asl = lci->asl;
	numvals = asl->i.numvals;
	nnum = lci->nnum;
	w = lci->w;
	nv = n_var;

	rv = 0;
 top:
	switch(*o) {

	  case nOPRET:
		if ((bt = lci->bt) && bt->lc[0])
			rv = &bt->h;
		else {
			if (!(og = getog(lci, o[1])))
				return 0;
			lf = getlf(lci, og);
			rv = &lf->h;
			}
		goto done;

	  case nOPPLUS:
		neg = 0;
 plusminus:
		if (!(og1 = getog(lci, o[2]))) {
 bad:
			return 0;
			}
		if (!(og2 = getog(lci, o[3])))
			goto bad;
		if (neg) {
			for(og = og2; og; og = og->next)
				og->coef = -og->coef;
			}
		if (!(og = ogsum(lci, og1, og2)))
			og = new_og(lci, -1, 0.);
		lf = getlf(lci, og);
		w[o[1]] = &lf->h;
		o += 4;
		goto top;

	  case nOPMINUS:
		neg = 1;
		goto plusminus;

	  case nOPMULT:
		j = o[3];
		if ((i = o[2]) < 0) {
			if ((i += nnum) < 0)
				goto bad;
			t = numvals[i];
			}
		else if (j < 0) {
			if ((j += nnum) < 0)
				goto bad;
			t = numvals[j];
			j = i;
			}
		else
			goto bad;
		if (!(og = getog(lci, j)))
			goto bad;
 finishmult:
		lf = getlf(lci, og);
		while(og) {
			og->coef *= t;
			og = og->next;
			}
		w[o[1]] = &lf->h;
		o += 4;
		goto top;

	  case nOPDIV0:
	  case nOPDIV1:
		if ((j = o[3]) >= 0 || (j += nnum) < 0)
			goto bad;
		if (!(t = numvals[j]))
			goto bad;
		if (!(og = getog(lci, o[2])))
			goto bad;
		t =  1. / t;
		goto finishmult;

	  case nOPUMINUS:
		if (!(og = getog(lci, o[2])))
			goto bad;
		lf = getlf(lci, og);
		do og->coef = -og->coef; while((og = og->next));
		w[o[1]] = &lf->h;
		o += 3;
		goto top;

	  case nOPLE:
		j = 0;
		goto finish_lincon;

	  case nOPGE:
		j = 1;
		goto finish_lincon;

	  case nOPNE:
		if (lci->bt)
			goto bad;
		if (!(bt = B_test(lci, o[0], o[2], o[3])))
			goto bad;
 have_bt:
		lci->bt = bt;
		w[o[1]] = &bt->h;
		o += 4;
		goto top;

	  case nOPEQ:
		j = 2;
 finish_lincon:
		if (!lci->bt && (bt = B_test(lci, o[0], o[2], o[3])))
			goto have_bt;
		if (!(og1 = getog(lci, o[2]))
		 || !(og2 = getog(lci, o[3])))
			goto bad;
		if (!lci->plc)
			goto bad;
		t = 0.;
		if (og1->varno < 0) {
			t = -og1->coef;
			og1 = og_free1(lci, og1);
			}
		if (og2->varno < 0) {
			t += og2->coef;
			og2 = og_free1(lci, og2);
			}
		for(og = og2; og; og = og->next)
			og->coef = -og->coef;
		if (!og1)
			og = og2;
		else if (!og2)
			og = og1;
		else
			og = ogsum(lci, og1, og2);
		lc = (Lconstr*)tmem(lci, sizeof(Lconstr));
		lc->h.type = Contype;
		lc->h.n = j;
		lc->rhs = t;
		lc->og = og;
		*lci->plc = lc;
		*(lci->plc = &lc->next) = 0;
		w[o[1]] = &lc->h;
		o += 4;
		goto top;

	  case nOPSUMLIST:
		k = o[2] + 3;
		t = 0.;
		if (!(zc = lci->zc)) {
			lci->zc = zc = (int*)new_mblk(htcl(nv*sizeof(int)));
			memset(zc, 0, nv*sizeof(int));
			}
		z = lci->z;
		r = lci->r;
		nz = 0;
		for(j = 3; j < k; ++j) {
			if ((i = o[j]) < 0) {
				if ((i += nnum) < 0)
					goto bad;
				t += numvals[i];
				}
			else if (i < nv) {
				if (!zc[i]) {
					z[nz] = i;
					zc[i] = ++nz;
					r[i] = 1.;
					}
				else if ((r[i] += 1.) == 0.) {
					z[--nz] = z[zc[i]-1];
					zc[i] = 0;
					}
				}
			else {
				lf = (Linform*)w[i];
				if (lf->h.type != Lintype) {
					while(nz > 0)
						zc[z[--nz]] = 0;
					goto bad;
					}
				lfd = lf;
				if (i < lci->mxv)
					lfd = 0;
				for(pog = &lf->og; (og = *pog); pog = &og->next) {
					if ((m = og->varno) < 0)
						t += og->coef;
					else if (!zc[m]) {
						z[nz] = m;
						zc[m] = ++nz;
						r[m] = og->coef;
						}
					else if ((r[m] += og->coef) == 0.) {
						z[--nz] = z[zc[m]-1];
						zc[m] = 0;
						}
					}
				if (lfd) {
					w[i] = 0;
					*pog = lci->freeog;
					lci->freeog = lf->og;
					lf->og = 0;
					*(Linform**)lf = lci->lffree;
					lci->lffree = lf;
					}
				}
			}
		if (nz > 1)
			qsortv(z, nz, sizeof(int), intcomp, NULL);
		og = 0;
		while(nz > 0) {
			zc[i = z[--nz]] = 0;
			(og1 = new_og(lci, i, r[i]))->next = og;
			og = og1;
			}
		if (t) {
			(og1 = new_og(lci, -1, t))->next = og;
			og = og1;
			}
		lf = getlf(lci, og);
		w[o[1]] = &lf->h;
		o += k;
		goto top;

	  case OPCOPY:
		w[o[1]] = w[o[2]];
		o += 3;
		goto top;

	  Opalign(OPORalign)
	  case nOPOR:
		if (lci->plc)
			goto bad;
		bt = (Btest*)w[o[2]];
		if (bt != lci->bt || !bt)
			goto bad;
		lci->plc = &bt->lc[0];
		w[o[1]] = 0;
		pop = (int**)&o[3];
		o = (int*)&pop[1];
		goto top;

	  Opalign(OPANDalign)
	  case nOPAND:
		if (!lci->plc)
			goto bad;
		pop = (int**)&o[3];
		o = (int*)&pop[1];
		goto top;

	  Opalign(OPIMPELSE_align)
	  case nOPIMPELSE:
		if (lci->plc || lci->iego)
			goto bad;
		bt = (Btest*)w[o[2]];
		if (bt != lci->bt || !bt)
			goto bad;
		lci->plc = &bt->lc[1];
		w[o[1]] = 0;
		pop = (int**)&o[5];
		o = pop[0];
		lci->iego = pop[1];
		goto top;

#ifdef X64_bit_pointers
	  case OPGOTOalign:
	  case OPNEXTBLKalign:
		++o;
#endif
	  case OPGOTO:
	  case OPNEXTBLK:
		o = *(int**)(o+1);
		if (*o != nOPRET || !lci->iego)
			goto top;
		o = lci->iego;
		lci->iego = 0;
		lci->plc = 0;
		if ((bt = lci->bt))
			lci->plc = &bt->lc[0];
		goto top;

	  default:
		goto bad;

	  }
 done:
	return rv;
	}

 static head *
Leval2(LCADJ_Info *lci, int *o)
{
	ASL *asl;
	Btest *bt;
	Lconstr *lc;
	Linform *lf, *lfd;
	head *rv, **w;
	int i, j, k, m, neg, nnum, nv, nz, **pop, *rp, *z, *zc;
	ograd *og, *og1, *og2, **pog;
	real *numvals, *r, t;

	asl = lci->asl;
	numvals = asl->i.numvals;
	nnum = lci->nnum;
	w = lci->w;
	nv = n_var;

	rv = 0;
 top:
	switch(*o) {

	  case OPRET:
		if ((bt = lci->bt) && bt->lc[0])
			rv = &bt->h;
		else {
			if (!(og = getog(lci, o[1])))
				return 0;
			lf = getlf(lci, og);
			rv = &lf->h;
			}
		goto done;

	  case OPMINUS0:
		neg = 1;
		goto plusminus0;

	  case OPPLUS0:
		neg = 0;
 plusminus0:
		rp = &o[1];
		i = o[2];
		j = o[3];
		o += 4;
		goto plusminus;

	  case OPMINUS01:
	  case OPMINUS10:
	  case OPMINUS2:
		neg = 1;
		goto plusminus1;

	  case OPPLUS01:
	  case OPPLUS10:
	  case OPPLUS2:
		neg = 0;
 plusminus1:
		rp = &o[2];
		i = o[3];
		j = o[4];
		o += 5;
 plusminus:
		if (!(og1 = getog2(lci, i))) {
 bad:
			return 0;
			}
		if (!(og2 = getog2(lci, j)))
			goto bad;
		if (neg) {
			for(og = og2; og; og = og->next)
				og->coef = -og->coef;
			}
		if (!(og = ogsum(lci, og1, og2)))
			og = new_og(lci, -1, 0.);
		lf = getlf(lci, og);
		w[*rp] = &lf->h;
		goto top;

	  case OPMULT0:
		rp = &o[1];
		i = o[2];
		j = o[3];
		o += 4;
		goto finishmult1;

	  case OPMULT01:
	  case OPMULT10:
	  case OPMULT2:
		rp = &o[2];
		i = o[3];
		j = o[4];
		o += 5;
 finishmult1:
		if (i < 0) {
			if ((i += nnum) < 0)
				goto bad;
			t = numvals[i];
			}
		else if (j < 0) {
			if ((j += nnum) < 0)
				goto bad;
			t = numvals[j];
			j = i;
			}
		else
			goto bad;
		if (!(og = getog2(lci, j)))
			goto bad;
 finishmult:
		lf = getlf(lci, og);
		while(og) {
			og->coef *= t;
			og = og->next;
			}
		w[*rp] = &lf->h;
		goto top;

	  case OPDIV0:
		rp = &o[1];
		i = o[2];
		j = o[3];
		o += 4;
		goto finishdiv;

	  case OPDIV10:
		rp = &o[2];
		i = o[3];
		j = o[4];
		o += 5;
 finishdiv:
		if (j >= 0 || (j += nnum) < 0)
			goto bad;
		if (!(t = numvals[j]))
			goto bad;
		if (!(og = getog2(lci, o[2])))
			goto bad;
		t =  1. / t;
		goto finishmult;

	  case OPUMINUS0:
		rp = &o[1];
		i = o[2];
		o += 3;
		goto finish_uminus;

	  case OPUMINUS1:
		rp = &o[2];
		i = o[3];
		o += 4;
 finish_uminus:
		if (!(og = getog2(lci, i)))
			goto bad;
		lf = getlf(lci, og);
		do og->coef = -og->coef; while((og = og->next));
		w[*rp] = &lf->h;
		goto top;

	  case n_OPLE:
		j = 0;
		goto finish_lincon;

	  case n_OPGE:
		j = 1;
		goto finish_lincon;

	  case n_OPNE:
		if (lci->bt)
			goto bad;
		if (!(bt = B_test2(lci, o[0], o[2], o[3])))
			goto bad;
 have_bt:
		lci->bt = bt;
		w[o[1]] = &bt->h;
		o += 4;
		goto top;

	  case OPEQ:
		j = 2;
 finish_lincon:
		if (!lci->bt && (bt = B_test2(lci, o[0], o[2], o[3])))
			goto have_bt;
		if (!(og1 = getog2(lci, o[2]))
		 || !(og2 = getog2(lci, o[3])))
			goto bad;
		if (!lci->plc)
			goto bad;
		t = 0.;
		if (og1->varno < 0) {
			t = -og1->coef;
			og1 = og_free1(lci, og1);
			}
		if (og2->varno < 0) {
			t += og2->coef;
			og2 = og_free1(lci, og2);
			}
		for(og = og2; og; og = og->next)
			og->coef = -og->coef;
		if (!og1)
			og = og2;
		else if (!og2)
			og = og1;
		else
			og = ogsum(lci, og1, og2);
		lc = (Lconstr*)tmem(lci, sizeof(Lconstr));
		lc->h.type = Contype;
		lc->h.n = j;
		lc->rhs = t;
		lc->og = og;
		*lci->plc = lc;
		*(lci->plc = &lc->next) = 0;
		w[o[1]] = &lc->h;
		o += 4;
		goto top;

	  case OPSUMLIST1:
		++o;
	  case OPSUMLIST0:
		k = o[2] + 3;
		t = 0.;
		if (!(zc = lci->zc)) {
			lci->zc = zc = (int*)new_mblk(htcl(nv*sizeof(int)));
			memset(zc, 0, nv*sizeof(int));
			}
		z = lci->z;
		r = lci->r;
		nz = 0;
		for(j = 3; j < k; ++j) {
			if ((i = o[j]) < 0) {
				if ((i += nnum) < 0)
					goto bad;
				t += numvals[i];
				}
			else if (i < nv) {
				if (!zc[i]) {
					z[nz] = i;
					zc[i] = ++nz;
					r[i] = 1.;
					}
				else if ((r[i] += 1.) == 0.) {
					z[--nz] = z[zc[i]-1];
					zc[i] = 0;
					}
				}
			else {
				lf = (Linform*)w[i];
				if (lf->h.type != Lintype) {
					while(nz > 0)
						zc[z[--nz]] = 0;
					goto bad;
					}
				lfd = lf;
				if (i < lci->mxv)
					lfd = 0;
				for(pog = &lf->og; (og = *pog); pog = &og->next) {
					if ((m = og->varno) < 0)
						t += og->coef;
					else if (!zc[m]) {
						z[nz] = m;
						zc[m] = ++nz;
						r[m] = og->coef;
						}
					else if ((r[m] += og->coef) == 0.) {
						z[--nz] = z[zc[m]-1];
						zc[m] = 0;
						}
					}
				if (lfd) {
					w[i] = 0;
					*pog = lci->freeog;
					lci->freeog = lf->og;
					lf->og = 0;
					*(Linform**)lf = lci->lffree;
					lci->lffree = lf;
					}
				}
			}
		if (nz > 1)
			qsortv(z, nz, sizeof(int), intcomp, NULL);
		og = 0;
		while(nz > 0) {
			zc[i = z[--nz]] = 0;
			(og1 = new_og(lci, i, r[i]))->next = og;
			og = og1;
			}
		if (t) {
			(og1 = new_og(lci, -1, t))->next = og;
			og = og1;
			}
		lf = getlf(lci, og);
		w[o[1]] = &lf->h;
		o += k;
		goto top;

	  case OPCOPY0:
		w[o[1]] = w[o[2]];
		o += 3;
		goto top;

	  case OPCOPY1:
	  case OPCOPY1a:
		w[o[2]] = w[o[3]];
		o += 4;
		goto top;

	  Opalign(OP_ORalign)
	  case n_OPOR:
		if (lci->plc)
			goto bad;
		bt = (Btest*)w[o[2]];
		if (bt != lci->bt || !bt)
			goto bad;
		lci->plc = &bt->lc[0];
		w[o[1]] = 0;
		pop = (int**)&o[3];
		o = (int*)&pop[1];
		goto top;

	  Opalign(OP_ANDalign)
	  case n_OPAND:
		if (!lci->plc)
			goto bad;
		pop = (int**)&o[3];
		o = (int*)&pop[1];
		goto top;

	  Opalign(OP_IMPELSE_align)
	  case n_OPIMPELSE:
		if (lci->plc || lci->iego)
			goto bad;
		bt = (Btest*)w[o[1]];
		if (bt != lci->bt || !bt)
			goto bad;
		lci->plc = &bt->lc[1];
		w[o[1]] = 0;
		pop = (int**)&o[2];
		o = pop[0];
		lci->iego = pop[1];
		goto top;

#ifdef X64_bit_pointers
	  case OP_GOTOalign:
	  case OPGOTO2align:
	  case OP_NEXTBLKalign:
	  case OPGOTOFalign:
	  case OPGOTOF2align:
	  case OPGOTOF2nalign:
		++o;
#endif
	  case OP_GOTO:
	  case OPGOTO2:
	  case OP_NEXTBLK:
	  case OPGOTOF:
	  case OPGOTOF2:
	  case OPGOTOF2n:
		o = *(int**)(o+1);
		if (*o != OPRET || !lci->iego)
			goto top;
		o = lci->iego;
		lci->iego = 0;
		lci->plc = 0;
		if ((bt = lci->bt))
			lci->plc = &bt->lc[0];
		goto top;

	  default:
		goto bad;

	  }
 done:
	return rv;
	}

 static void
chunkfree(ASL *asl, Tchunk **ptc)
{
	Tchunk *tc, *tc1;
	for(tc1 = *ptc; (tc = tc1); ) {
		tc1 = tc->prev;
		del_mblk(tc);
		}
	*ptc = 0;
	}

 static int
add_indicator(int ci, LCADJ_Info *lci, int *e)
{
	Btest *bt;
	Lconstr *lc;
	head *h;
	int Compl, i, j, k, vk, *z;
	ograd *og;
	real *r;

	if (lci->tchunks)
		chunkfree(lci->asl, &lci->tchunks);
	lci->freeog = 0;
	lci->lffree = 0;
	lci->tfree = lci->tfree0;
	lci->ntfree = Gulp;
	lci->bt = 0;
	lci->plc = 0;
	lci->iego = 0;
	if (!(h = lci->Leval(lci, e)) || h->type != BTtype) {
		lci->errinfo[0] = ci;
		return 1;
		}
	bt = (Btest*)h;
	vk = bt->h.n;
	Compl = bt->bval;
	r = lci->r;
	z = lci->z;
	for(j = 0; j < 2; ++j, Compl = 1 - Compl) {
		for(lc = bt->lc[j]; lc; lc = lc->next) {
			og = lc->og;
			for(k = 0; og; og = og->next, ++k) {
				r[k] = og->coef;
				z[k] = og->varno;
				}
			i = lci->add_indic(lci->v, vk, Compl, lc->h.n, k, z, r, lc->rhs);
			if (i) {
				lci->errinfo[0] = i;
				return 3;
				}
			}
		}
	return 0;
	}

 int
indicator_constrs_ASL(ASL *asl, void *v, Add_Indicator add_indic, int errinfo[2])
{
	LCADJ_Info lci;
	Linform *lf;
	cde *logc;
	cexp *ce;
	cexp1 *ce1;
	lincoef *lc, *lce;
	ograd *og, *og1;
	head *h;
	int i, n, nlogc, nbc, nbc1, nw, rc;
	linpart *lp;
	real chunk1[Gulp];

	if (!(nlogc = n_lcon))
		return 0;
	memset(&lci, 0, sizeof(lci));
	lci.Leval = asl->i.ASLtype <= ASL_read_fg ? Leval1 : Leval2;
	lci.v = v;
	lci.tfree0 = chunk1;
	n = n_var;
	nw = asl->i.wlen/sizeof(real) - n;
	i = htcl(n*(sizeof(int) + sizeof(ograd*) + sizeof(real)) + nw*sizeof(head*));
	lci.r = (real*)new_mblk(i);
	lci.w0 = (head**)(lci.r + n);
	lci.w = lci.w0 - n;
	lci.s = (ograd**)(lci.w0 + nw);
	lci.z = (int*)(lci.s + n);
	memset(lci.w0, 0, n*sizeof(ograd*) + nw*sizeof(head*));
	lci.n1og = (sizeof(ograd) + sizeof(real) - 1) / sizeof(real);
	lci.asl = asl;
	lci.mxv = asl->i.maxvar;
	lci.nnum = asl->i.numlen / sizeof(real);
	lci.nv = n_var;
	lci.numvals = asl->i.numvals;
	lci.nlv[1] = i = nlvb;
	lci.nlv[0] = i - nlvbi;
	lci.nlv[3] = i = nlvc;
	lci.nlv[2] = i - nlvci;
	if (nlvo > nlvc)
		i = nlvo;
	lci.nlv[5] = i;
	lci.nlv[4] = i - nlvoi;
	lci.errinfo = errinfo;
	lci.add_indic = add_indic;
	lci.zc = 0;
	nbc = combc;
	ce = ((ASL_fg*)asl)->I.cexps_;
	for(i = 0; i < nbc; ++i) {
		if ((h = lci.Leval(&lci, ce[i].o.e)) && h->type == Lintype && (lp = ce[i].lp)) {
			lc = lp->lc;
			lce = lc + lp->n;
			og = 0;
			while(lce > lc) {
				--lce;
				(og1 = new_og(&lci, lce->varno, lce->coef))->next = og;
				og = og1;
				}
			lf = (Linform*)h;
			lf->og = ogsum(&lci, lf->og, og);
			}
		lci.w0[i] = h;
		}
	nbc1 = comc1;
	ce1 = ((ASL_fg*)asl)->I.cexps1_;
	for(i = 0; i < nbc1; ++i) {
		if ((h = lci.Leval(&lci, ce1[i].o.e)) && h->type == Lintype && (lp = ce1[i].lp)) {
			lc = lp->lc;
			lce = lc + lp->n;
			og = 0;
			while(lce > lc) {
				--lce;
				(og1 = new_og(&lci, lce->varno, lce->coef))->next = og;
				og = og1;
				}
			lf = (Linform*)h;
			lf->og = ogsum(&lci, lf->og, og);
			}
		lci.w0[nbc+i] = h;
		}
	logc = ((ASL_fg*)asl)->I.lcon_de_;
	for(i = rc = 0; i < nlogc; ++i)
		if ((rc = add_indicator(i, &lci, logc[i].o.e)))
			break;
	if (lci.zc)
		del_mblk(lci.zc);
	if (lci.tchunks)
		chunkfree(asl, &lci.tchunks);
	if (lci.tchunks0)
		chunkfree(asl, &lci.tchunks0);
	del_mblk(lci.r);
	return rc;
	}
