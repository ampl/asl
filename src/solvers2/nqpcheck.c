/*******************************************************************
Copyright (C) 2017, 2018, 2019 AMPL Optimization, Inc.; written by David M. Gay.

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

#define GULP		200

 typedef struct
dyad {
	struct dyad *next;
	ograd *Lq, *Rq;
	} dyad;

 typedef struct
term {
	dyad	*Q, *Qe;
	ograd	*L, *Le;
	} term;

 typedef struct
Static {
	ASL_fg *asl;
	fint *_s_s, *_s_z;
	double *_s_x, *cv;
	term *_freeterm, *_term_block;
	ograd *_freeog, *_ograd_block;
	dyad *_freedyad, *_dyad_block;
	int _zerodiv;
	term **W, **_cterms;
	int _dyad_ntogo, nnum, nv1, nv2, _ograd_ntogo, _term_ntogo;
	Char **_M1state1, **_M1state2;
	term *(*comterm)(struct Static*, int);
	} Static;

#define M1state1	S->_M1state1
#define M1state2	S->_M1state2
#define cterms		S->_cterms
#define dyad_block	S->_dyad_block
#define dyad_ntogo	S->_dyad_ntogo
#define freedyad	S->_freedyad
#define freeog		S->_freeog
#define freeterm	S->_freeterm
#define ograd_block	S->_ograd_block
#define ograd_ntogo	S->_ograd_ntogo
#define s_s		S->_s_s
#define s_x		S->_s_x
#define s_z		S->_s_z
#define term_block	S->_term_block
#define term_ntogo	S->_term_ntogo
#define zerodiv		S->_zerodiv

 static void
free_blocks(Static *S)
{
	M1free_ASL(&S->asl->i, M1state1, M1state2);
	}

 static void
free_term(Static *S, term *t)
{
	t->Q = (dyad *)freeterm;
	freeterm = t;
	}

 static term *
new_term(Static *S, ograd *o)
{
	term *rv;

	if ((rv = freeterm))
		freeterm = (term *)rv->Q;
	else {
		if (!term_ntogo) {
			term_block = (term *)M1alloc_ASL(&S->asl->i,
				GULP*sizeof(term));
			term_ntogo = GULP;
			}
		rv = term_block++;
		--term_ntogo;
		}
	rv->L = rv->Le = o;
	rv->Q = rv->Qe = 0;
	return rv;
	}

 static void
free_og(Static *S, ograd *o)
{
	o->next = freeog;
	freeog = o;
	}

 static ograd *
new_og(Static *S, ograd *next, int i, real v)
{
	ograd *rv;

	if ((rv = freeog))
		freeog = rv->next;
	else {
		if (!ograd_ntogo) {
			ograd_block = (ograd *)M1alloc_ASL(&S->asl->i,
				GULP*sizeof(ograd));
			ograd_ntogo = GULP;
			}
		rv = ograd_block++;
		--ograd_ntogo;
		}
	rv->next = next;
	rv->varno = i;
	rv->coef = v;
	return rv;
	}

 static ograd *
ogdup(Static *S, ograd *og, ograd **oge)
{
	ograd *og0, *og1;

	og0 = og1 = new_og(S, 0, og->varno, og->coef);
	while((og = og->next))
		og1 = og1->next = new_og(S, 0, og->varno, og->coef);
	if (oge)
		*oge = og1;
	return og0;
	}

 static int
count(Static *S, ograd **ogp)
{
	int i, rv, nz;
	fint *s, *z;
	double t, *x;
	ograd *og, *og1;

	s = s_s;
	x = s_x;
	z = s_z;

	t = 0;
	nz = rv = 0;
	for(og = *ogp; og; og = og1) {
		og1 = og->next;
		if ((i = og->varno) < 0)
			t += og->coef;
		else if (!s[i]++)
			x[z[nz++] = i] = og->coef;
		else
			x[i] += og->coef;
		free_og(S, og);
		}
	while(nz > 0) {
		s[i = z[--nz]] = 0;
		if (x[i]) {
			og = new_og(S, og, i, x[i]);
			rv++;
			}
		}
	if (t)
		og = new_og(S, og, -1, t);
	*ogp = og;
	return rv;
	}

 static void
free_dyad(Static *S, dyad *t)
{
	t->next = freedyad;
	freedyad = t;
	}

 static dyad *
new_dyad(Static *S, dyad *next, ograd *L, ograd *R, int permute)
{
	dyad *rv;
	ograd *t;

	if (permute) {
		if (L == R) {
			count(S, &L);
			R = L;
			}
		else if (count(S, &L) > count(S, &R)) {
			t = L;
			L = R;
			R = t;
			}
		if (!L) /* e.g., 0*x*x from <<0;0,0>>x*x */
			/* with AMPL version < 20000216. */
			return next;
		}
	if ((rv = freedyad))
		freedyad = rv->next;
	else {
		if (!dyad_ntogo) {
			dyad_block = (dyad *)M1alloc_ASL(&S->asl->i,
				GULP*sizeof(dyad));
			dyad_ntogo = GULP;
			}
		rv = dyad_block++;
		--dyad_ntogo;
		}
	rv->next = next;
	rv->Lq = L;
	rv->Rq = R;
	return rv;
	}

 static term *
termsum(Static *S, term *L, term *R)
{
	if (!L || !R)
		return 0;
	if (L->Qe && (L->Qe->next = R->Q))
		L->Qe = R->Qe;
	else if (R->Q) {
		L->Q = R->Q;
		L->Qe = R->Qe;
		}
	if (L->Le && (L->Le->next = R->L))
		L->Le = R->Le;
	else if (R->L) {
		L->L = R->L;
		L->Le = R->Le;
		}
	free_term(S, R);
	return L;
	}

 static term *
scale(Static *S, term *T, register double t)
{
	register ograd *og;
	register dyad *d;

	if (T) {
		for(d = T->Q; d; d = d->next) {
			if (d->Lq == d->Rq)
				d->Rq = ogdup(S, d->Lq, 0);
			for(og = d->Lq; og; og = og->next)
				og->coef *= t;
			}
		for(og = T->L; og; og = og->next)
			og->coef *= t;
		}
	return T;
	}

 static term *ewalk1(Static*, int*), *ewalk2(Static*, int*);

 static term *
comterm1(Static *S, int i)
{
	ASL_fg* asl;
	cexp *c;
	cexp1 *c1;
	lincoef *L, *Le;
	linpart *lp;
	term *T;

	asl = S->asl;
	if (i < ncom0) {
		c = cexps + i;
		T = ewalk1(S, c->o.e);
		lp = c->lp;
		}
	else {
		c1 = cexps1 + (i - ncom0);
		T = ewalk1(S, c1->o.e);
		lp = c1->lp;
		}
	if (lp && T) {
		L = lp->lc;
		for(Le = L + lp->n; L < Le; L++)
			T = termsum(S, T, new_term(S, new_og(S, 0, L->varno, L->coef)));
		}
	return T;
	}

 static void
comeval1(Static *S, int i, int ie)
{
	term **Cterms;

	for(Cterms = cterms; i < ie; ++i)
		Cterms[i] = comterm1(S, i);
	}

 static term *
comterm2(Static *S, int i)
{
	ASL_pfgh* asl;
	cexp2 *c;
	int *e;
	lincoef *L, *Le;
	linpart *lp;
	term *T, *T1;

	asl = (ASL_pfgh*)S->asl;
	c = asl->I.cexps2_ + i;
	T = 0;
	if ((e = c->o.e) && !(T = ewalk2(S, e)))
		return 0;
	if ((lp = c->lp)) {
		L = lp->lc;
		for(Le = L + lp->n; L < Le; L++) {
			T1 = new_term(S, new_og(S, 0, L->varno, L->coef));
			if (!T)
				T = T1;
			else
				T = termsum(S, T, T1);
			}
		}
	return T;
	}

 static void
comeval2(Static *S, int i, int ie)
{
	ASL_pfgh *asl;
	cexp2 *c, *c1;
	int *dvsp0, j, k, nv1, *ndvsp;
	term **Cterms;

	asl = (ASL_pfgh*)S->asl;
	dvsp0 = asl->P.dvsp0;
	ndvsp = asl->P.ndvsp;
	c = asl->I.cexps2_;
	for(Cterms = cterms; i < ie; ++i) {
		if ((k = ndvsp[i])) {
			nv1 = S->nv1 >> 2;
			k += j = dvsp0[i];
			do {
				c1 = c + j;
				Cterms[c1->varno-nv1] = ewalk2(S, c1->o.e);
				} while(++j < k);
			}
		Cterms[i] = comterm2(S, i);
		}
	}

 static term *
termdup(Static *S, term *T)
{
	term *rv;
	ograd *og, *oge;
	dyad *Q, *Q1;

	if ((og = oge = T->L))
		og = ogdup(S, og, &oge);
	rv = new_term(S, og);
	rv->Le = oge;
	if (!(Q = T->Q))
		return rv;
	Q1 = rv->Qe = new_dyad(S, 0, ogdup(S, Q->Lq,0), ogdup(S, Q->Rq,0), 1);
	while((Q = Q->next))
		Q1 = new_dyad(S, Q1, ogdup(S, Q->Lq,0), ogdup(S, Q->Rq,0), 1);
	rv->Q = Q1;
	return rv;
	}

 static term *
tval(Static *S, int i)
{
	term *L;

	if (i >= S->nv2)
		return S->W[i];
	if (i < 0)
		return new_term(S, new_og(S, 0, -1 , S->cv[i + S->nnum]));
	if (i < S->nv1)
		return new_term(S, new_og(S, 0, i, 1.));
	i -= S->nv1;
	if (!(L = cterms[i]))
		return L;
	return termdup(S, L);
	}

 static term *
tval2(Static *S, int i)
{
	term *L;

	if (i >= S->nv2)
		return S->W[i];
	if (i < 0)
		return new_term(S, new_og(S, 0, -1 , S->cv[i + S->nnum]));
	if (i < S->nv1)
		return new_term(S, new_og(S, 0, i>>2, 1.));
	i -= S->nv1;
	i >>= 2;
	if (!(L = cterms[i]))
		return L;
	return termdup(S, L);
	}

 static term *
ewalk1(Static *S, int *e)
{
	int i, j, k, n;
	ograd *o, *oR;
	term *L, *R, *T, **w;

	w = S->W;

 top:
	switch(*e) {
	  case nOPRET:
		return tval(S, e[1]);

	  case nOPPLUS:
		if (!(L = tval(S, e[2])) || !(R = tval(S, e[3])))
			break;
		w[e[1]] = termsum(S, L, R);
		e += 4;
		goto top;

	  case nOPMINUS:
		if (!(L = tval(S, e[2])) || !(R = tval(S, e[3])))
			break;
		w[e[1]] = termsum(S, L, scale(S, R, -1.));
		e += 4;
		goto top;

	  case nOPUMINUS:
		if (!(L = tval(S, e[2])))
			break;
		w[e[1]] = scale(S, L, -1.);
		e += 3;
		goto top;

	  case nOPMULT:
		if (!(L = tval(S, e[2])) || !(R = tval(S, e[3])))
			break;
		i = e[1];
		e += 4;
		if (L->Q) {
			if (R->Q)
				break;
 qscale:
			o = R->L;
			if (o->next || o->varno >= 0)
				break;
			w[i] = scale(S, L, o->coef);
			free_og(S, o);
			free_term(S, R);
			goto top;
			}
		if (R->Q) {
			T = L;
			L = R;
			R = T;
			goto qscale;
			}
		o = L->L;
		oR = R->L;
		if (o->next || o->varno >= 0) {
			if (oR->next || oR->varno >= 0) {
				L->L = L->Le = 0;
				if (!(L->Q = L->Qe = new_dyad(S, 0,o,oR,1)))
					goto no_L;
				}
			else {
				scale(S, L, oR->coef);
				free_og(S, oR);
				}
			free_term(S, R);
			w[i] = L;
			goto top;
			}
 no_L:
		scale(S, R, o->coef);
		free_og(S, o);
		free_term(S, L);
		w[i] = R;
		goto top;

	  case nOPDIV2:
	  case nOPDIV3:
		/* only allow division by a constant */
		break;
	  case nOPDIV0:
	  case nOPDIV1:
		if (!(L = tval(S, e[2])) || !(R = tval(S, e[3])))
			break;
		i = e[1];
		e += 4;
		o = R->L;
		if (R->Q || o->next || o->varno >= 0)
			break;
		if (!o->coef) {
			zerodiv++;
			break;
			}
		w[i] = scale(S, L, 1./o->coef);
		free_og(S, o);
		free_term(S, R);
		goto top;

	  case nOPSUMLIST:
		n = e[2];
		j = 3;
		k = j + n;
		if (!(L = tval(S, e[j])))
			break;
		while(++j < k) {
			if (!(R = tval(S, e[j])))
				return 0;
			L = termsum(S, L, R);
			}
		w[e[1]] = L;
		e += k;
		goto top;

	  case OP2POW0:
	  case OP2POW1:
		i = e[1];
		if (!(L = tval(S, e[2])))
			break;
		if (!L || L->Q)
			break;
		e += 3;
		o = L->L;
		if (!o->next && o->varno < 0) {
			o->coef *= o->coef;
			w[i] = L;
			goto top;
			}
		L->Q = L->Qe = new_dyad(S, 0,o,o,1);
		L->L = L->Le = 0;
		w[i] = L;
		goto top;

	  case OPGOTO:
	  case OPNEXTBLK:
		e = *(int**)(e+1);	/* for chaining blocks of code */
		goto top;

#ifdef X64_bit_pointers
	  case OPGOTOalign:
	  case OPNEXTBLKalign:
		e = *(int**)(e+2);
		goto top;
#endif
		}
	return 0; /* nonlinear */
	}

 static term *
ewalk2(Static *S, int *e)
{
	int i, j, k, n, **pop;
	ograd *o, *oR;
	term *L, *R, *T, **w;

	w = S->W;

 top:
	switch(*e) {
	  case OPRET:
		return tval2(S, e[1]);

	  case OPPLUS01:
	  case OPPLUS10:
	  case OPPLUS2:
		++e;
	  case OPPLUS0:
		if (!(L = tval2(S, e[2])) || !(R = tval2(S, e[3])))
			break;
		w[e[1]] = termsum(S, L, R);
		e += 4;
		goto top;

	  case OPMINUS01:
	  case OPMINUS10:
	  case OPMINUS2:
		++e;
	  case OPMINUS0:
		if (!(L = tval2(S, e[2])) || !(R = tval2(S, e[3])))
			break;
		w[e[1]] = termsum(S, L, scale(S, R, -1.));
		e += 4;
		goto top;

	  case OPUMINUS1:
		++e;
	  case OPUMINUS0:
		if (!(L = tval2(S, e[2])))
			break;
		w[e[1]] = scale(S, L, -1.);
		e += 3;
		goto top;

	  case OPMULT01:
	  case OPMULT10:
	  case OPMULT2:
		++e;
	  case OPMULT0:
		if (!(L = tval2(S, e[2])) || !(R = tval2(S, e[3])))
			break;
		i = e[1];
		e += 4;
		if (L->Q) {
			if (R->Q)
				break;
 qscale:
			o = R->L;
			if (o->next || o->varno >= 0)
				break;
			w[i] = scale(S, L, o->coef);
			free_og(S, o);
			free_term(S, R);
			goto top;
			}
		if (R->Q) {
			T = L;
			L = R;
			R = T;
			goto qscale;
			}
		o = L->L;
		oR = R->L;
		if (o->next || o->varno >= 0) {
			if (oR->next || oR->varno >= 0) {
				L->L = L->Le = 0;
				if (!(L->Q = L->Qe = new_dyad(S, 0,o,oR,1)))
					goto no_L;
				}
			else {
				scale(S, L, oR->coef);
				free_og(S, oR);
				}
			free_term(S, R);
			w[i] = L;
			goto top;
			}
 no_L:
		scale(S, R, o->coef);
		free_og(S, o);
		free_term(S, L);
		w[i] = R;
		goto top;

	  case OPDIV01:
	  case OPDIV2:
		/* only allow division by a constant */
		break;
	  case OPDIV10:
		++e;
	  case OPDIV0:
		if (!(L = tval2(S, e[2])) || !(R = tval2(S, e[3])))
			break;
		i = e[1];
		e += 4;
		o = R->L;
		if (R->Q || o->next || o->varno >= 0)
			break;
		if (!o->coef) {
			zerodiv++;
			break;
			}
		w[i] = scale(S, L, 1./o->coef);
		free_og(S, o);
		free_term(S, R);
		goto top;

	  case OPSUMLIST1:
		++e;
	  case OPSUMLIST0:
		n = e[2];
		j = 3;
		k = j + n;
		if (!(L = tval2(S, e[j])))
			break;
		while(++j < k) {
			if (!(R = tval2(S, e[j])))
				return 0;
			L = termsum(S, L, R);
			}
		w[e[1]] = L;
		e += k;
		goto top;

	  case OP_2POW1:
		++e;
	  case OP_2POW0:
		i = e[1];
		if (!(L = tval2(S, e[2])))
			break;
		if (!L || L->Q)
			break;
		e += 3;
		o = L->L;
		if (!o->next && o->varno < 0) {
			o->coef *= o->coef;
			w[i] = L;
			goto top;
			}
		L->Q = L->Qe = new_dyad(S, 0,o,o,1);
		L->L = L->Le = 0;
		w[i] = L;
		goto top;

	  case OP_GOTO:
	  case OP_NEXTBLK:
	  case OPGOTO2:
		e = *(int**)(e+1);	/* for chaining blocks of code */
		goto top;

	  case OPGOTOF:
		pop = (int**)&e[1];
		e = (int*)&pop[1];
		goto top;

	  case OPGOTOF2:
	  case OPGOTOF2n:
		pop = (int**)&e[1];
		e = (int*)&pop[2];
		goto top;

	  case OPGOTOMM:
		e += 2;
		goto top;

	  case OPRETB:
		++e;
		goto top;

#ifdef X64_bit_pointers
	  case OP_GOTOalign:
	  case OP_NEXTBLKalign:
	  case OPGOTO2align:
		e = *(int**)(e+2);
		goto top;

	  case OPGOTOFalign:
		pop = (int**)&e[2];
		e = (int*)&pop[1];
		goto top;

	  case OPGOTOF2align:
	  case OPGOTOF2nalign:
		pop = (int**)&e[2];
		e = (int*)&pop[2];
		goto top;
#endif
		}
	return 0; /* nonlinear */
	}

#ifdef __cplusplus
extern "C" {
 static int comp(const void*, const void*, void*),
	    lcmp(const void*, const void*, void*);
}
#endif

 static int
comp(const void *a, const void *b, void *v)
{
	Not_Used(v);
	return (*(ograd **)a)->varno - (*(ograd **)b)->varno;
	}

 static ograd *
sortq(ograd *og0, ograd **q)
{
	ograd *og, **q1;
	int n;

	for(q1 = q, og = og0; og; og = og->next)
		*q1++ = og;
	if ((n = q1 - q) > 1) {
		qsortv(q, n, sizeof(ograd *), comp, NULL);
		og0 = 0;
		do {
			og = *--q1;
			og->next = og0;
			og0 = og;
			} while(q1 > q);
		}
	return og0;
	}

 static double
dsort(Static *S, term *T, ograd **q, cgrad **cgp, ograd **ogp, int arrays)
{
	cgrad *cg;
	ograd *og, *og1;
	double t, t1, rv, *x;
	dyad *Q, **pQ;

	x = s_x;

	rv = 0;
	count(S, &T->L);	/* accumulate */
	if (arrays) {
		if (ogp)
			for(og = *ogp; og; og = og->next)
				x[og->varno] = og->coef;
		else
			for(cg = *cgp; cg; cg = cg->next)
				x[cg->varno] = cg->coef;
		}
	if ((og = T->L) && og->varno < 0) {
		rv = og->coef;
		og = og->next;
		}
	for(; og; og = og->next)
		x[og->varno] += og->coef;

	pQ = &T->Q;
	while((Q = *pQ)) {
		og = Q->Lq;
		og1 = Q->Rq;
		t = t1 = 0;
		if (og->varno < 0) {
			t = og->coef;
			og = og->next;
			}
		if (og1->varno < 0) {
			t1 = og1->coef;
			Q->Rq = og1 = og1->next;
			rv += t*t1;
			}
		if (t)
			for(; og1; og1 = og1->next)
				x[og1->varno] += t*og1->coef;
		if (t1)
			for(og1 = og; og1; og1 = og1->next)
				x[og1->varno] += t1*og1->coef;
		if ((Q->Lq = sortq(og, q))
		 && (Q->Rq = og == Q->Rq ? Q->Lq : sortq(Q->Rq, q)))
			pQ = &Q->next;
		else {
			*pQ = Q->next;
			free_dyad(S, Q);
			}
		}
	if (arrays) {
		if (ogp)
			for(og = *ogp; og; og = og->next)
				og->coef = x[og->varno];
		else
			for(cg = *cgp; cg; cg = cg->next)
				cg->coef = x[cg->varno];
		}
	return rv;
	}

 static void
free_oglist(Static *S, ograd *og)
{
	ograd *og1;

	for(; og; og = og1) {
		og1 = og->next;
		free_og(S, og);
		}
	}

 static void
cterm_free(Static *S, term **cte)
{
	term **ct, *t;
	dyad *d, *d1;

	for(ct = cterms; ct < cte; ct++)
		if ((t = *ct)) {
			free_oglist(S, t->L);
			d1 = t->Q;
			while((d = d1)) {
				d1 = d->next;
				free_oglist(S, d->Lq);
				if (d->Rq != d->Lq)
					free_oglist(S, d->Rq);
				free_dyad(S, d);
				}
			}
	free(cterms);
	}

 static int
lcmp(const void *a, const void *b, void *v)
{
	Not_Used(v);
	return (int)(*(fint *)a - *(fint *)b);
	}

#ifndef Fint
#define Fint fint
#define Fints fint
#endif

 Fints
mqpcheck_ASL(ASL *a, int co, fint **rowqp, Fint **colqp, real **delsqp)
{
	typedef struct dispatch {
		struct dispatch *next;
		fint i, j, jend;
		} dispatch;
	ASL_fg *asl;
	ASL_pfgh *asl2;
	Fint  *colq, *colq1, nelq;
	Objrep *od, **pod;
	Static SS, *S;
	cde *c;
	cgrad *cg, **cgp, **cgq, *cq;
	dispatch *cd, *cd0, **cdisp, **cdisp0, *cdnext, **cdp;
	dyad *d, *d1, **q, **q1, **q2, **qe;
	fint *rowq, *rowq0, *rowq1, *s, *z;
	fint ftn, i, icol, j, ncom, nv, nz, nz1;
	int *C1, C10, akind, arrays, *cm, co0, dv0, dv1, *e, i1, nw, *ov, pass, *vmi;
	linarg *la;
	lincoef *lc, *lce;
	linpart *lp;
	ograd *og, *og1, *og2, *ogs, **ogp;
	ps_func *p;
	psb_elem *b, *be;
	psg_elem *g, *ge;
	real *L, *U, *delsq, *delsq0, *delsq1, objadj, *oc, *oc0, scal, t, *x;
	term *(*Comterm)(Static*, int), *T, *T1, **W;
	void (*Comeval)(Static*, int, int);

	e = 0;
	p = 0;
	if ((akind = a->i.ASLtype) != ASL_read_pfgh) {
		akind = ASL_read_fg;
		Comterm = comterm1;
		Comeval = comeval1;
		asl2 = 0;
		C10 = a->i.ncom0_;
		}
	else {
		Comterm = comterm2;
		Comeval = comeval2;
		asl2 = (ASL_pfgh*)a;
		C10 = 0;
		}
	ASL_CHECK(a, akind, "nqpcheck");
	asl = (ASL_fg*)a;
	if (co >= n_obj || co < -n_con)
		return -3L;
	od = 0;
	co0 = co;
	if (co >= 0) {
		if ((pod = asl->i.Or) && (od = pod[co])) {
			co = od->ico;
			if (!(cgp = asl->i.Cgrad0))
				cgp = Cgrad;
			goto use_Cgrad;
			}
		dv0 = combc;
		dv1 = dv0 + como;
		C1 = asl->i.o_cexp1st_;
		if (asl2)
			p = asl2->P.ops + co;
		else {
			c = obj_de + co;
			e = c->o.e;
			}
		ogp = Ograd + co;
		cgp = 0;
		}
	else {
		co = -1 - co;
		if ((cm = asl->i.cmap))
			co = cm[co];
		cgp = Cgrad;
 use_Cgrad:
		dv0 = comb;
		dv1 = combc;
		C1 = asl->i.c_cexp1st_;
		if (asl2)
			p = asl2->P.cps + co;
		else {
			c = con_de + co;
			e = c->o.e;
			}
		cgp += co;
		ogp = 0;
		}
	if (p) {
		if (p->pi.b >= p->pi.be && p->g >= p->ge)
			return 0;
		}
	else if (*e == nOPRET && e[1] < 0)
		return 0;

	memset(S = &SS, 0, sizeof(Static));
	SS.asl = asl;
	SS.comterm = Comterm;
	SS.nnum = asl->i.numlen /sizeof(real);
	SS.nv1 = asl->i.n_var0 + asl->i.nsufext[ASL_Sufkind_var];
	if (asl2) {
		ncom = asl2->P.max_var1_ - (asl->i.n_var0 + asl->i.nsufext[ASL_Sufkind_var]);
		SS.nv2 = 4*asl2->P.max_var1_;
		}
	else {
		ncom = ncom0 + ncom1;
		SS.nv2 = SS.nv1 + ncom;
		}
	nw = (int)(asl->i.wlen / sizeof(real)) - SS.nv2;
	if (asl->i.vmap && !asl->i.vminv)
		/* keep vminv from being lost in free_blocks(S) below */
		get_vminv_ASL(a);
	M1state1 = asl->i.Mbnext;
	M1state2 = asl->i.Mblast;
	nv = n_var;
	s_x = x = (double *)Malloc(nv*(sizeof(double)+2*sizeof(fint)) + nw*sizeof(term*));
	W = (term**)(x + nv);
	memset(W, 0, nw*sizeof(term*));
	SS.cv = asl->i.numvals;
	s_z = z = (fint *)(W + nw);
	s_s = s = z + nv;
	memset(s, 0, nv*sizeof(fint));
	SS.W = W -= SS.nv2;	/* we only access S.W[i] for i >= nv2 */
	ftn = Fortran;

	delsq = delsq0 = delsq1 = 0; /* silence buggy "not-initialized" warning */
	colq = colq1 = 0;				/* ditto */
	rowq = rowq0 = rowq1 = 0;			/* ditto */
	cd0 = 0;					/* ditto */
	cdisp = cdisp0 = 0;				/* ditto */

	if (ncom) {
		SS._cterms = (term **)Malloc(ncom*sizeof(term*));
		memset(SS._cterms, 0, ncom*sizeof(term*));
		}

	if (asl2) {
		for(la = asl2->P.lalist; la; la = la->lnext) {
			if ((i1 = la->u.v)) {
				oc0 = la->oc;
				oc = oc0 + la->nnz;
				ov = la->ov + la->nnz;
				og = 0;
				while(oc > oc0)
					og = new_og(S, og, *--ov, *--oc);
				SS._cterms[i1-SS.nv1] = new_term(S, og);
				}
			}
		SS.nv1 *= 4;
		}
	arrays = 1;
	if (rowqp)
		*rowqp = 0;
	else
		arrays = 0;
	if (colqp)
		*colqp = 0;
	else
		arrays = 0;
	if (delsqp)
		*delsqp = 0;
	else
		arrays = 0;

	zerodiv = 0;
	if (comb)
		Comeval(S, 0, comb);
	if (dv1 > dv0)
		Comeval(S, dv0, dv1);
	if (C1 && C1[co] < C1[co+1])
		Comeval(S, C10 + C1[co], C10 + C1[co+1]);
	if (p) {
		T = 0;
		ge = p->ge;
		for(g = p->g; g < ge; ++g) {
			if (*g->o != OP2POW_g
			 || g->o[3] != OPRET)
				goto done;
			T1 = 0;
			og = 0;
			if ((b = g->pi.b)) {
			 	if (b->o.f
				 || b->o.e[0] != OPRET
				 || (i1 = b->o.e[1]) < 0
				 || !(T1 = tval2(S, i1)))
					goto done;
				if (T1->Q) {
					T = 0;
					goto done;
					}
				og = T1->L;
				}
			if (g->g0)
				og = new_og(S, og, -1, g->g0);
			if ((lp = g->L)) {
				lc = lp->lc;
				for(lce = lc + lp->n; lc < lce; ++lc)
					og = new_og(S, og, lc->varno, lc->coef);;
				}
			count(S, &og);
			if (!og)
				continue;
			if (!T1)
				T1 = new_term(S, og);
			T1->L = T1->Le = 0;
			og1 = og;
			if ((scal = g->scale) != 1.) {
				if (og->varno < 0) {
					ogs = og;
					og2 = og->next;
					}
				else {
					ogs = 0;
					og2 = og;
					}
				for(og1 = 0; og2; og2 = og2->next)
					og1 = new_og(S, og1, og2->varno, scal*og2->coef);
				if (ogs)
					og1 = new_og(S, og1, ogs->varno, scal*ogs->coef);
				}
			if ((T1->Q = T1->Qe = new_dyad(S, 0, og, og1, 0))) {
				if (!T)
					T = T1;
				else
					T = termsum(S, T, T1);
				}
			}
		be = p->pi.be;
		for(b = p->pi.b; b < be; ++b) {
			if (!(T1 = ewalk2(S, b->o.e)))
				goto done;
			if (!T)
				T = T1;
			else
				T = termsum(S, T, T1);
			}
		if (zerodiv)
			goto done;
		}
	else if (!(T = ewalk1(S, e)) || zerodiv) {
 done:
		if (SS._cterms)
			free(SS._cterms);
		free_blocks(S);
		free(x);
		return T ? -2L : -1L;
		}

	if (cterms)
		cterm_free(S, cterms + ncom);
	if (od) {
		cgq = &od->cg;
		for(i = 0, cg = *cgp; cg; cg = cg->next) {
			if (cg->coef != 0.)
				++i;
			}
		if (i) {
			cq = (cgrad*)Malloc(i*sizeof(cgrad));
			for(cg = *cgp; cg; cg = cg->next) {
				if (cg->coef != 0.) {
					*cgq = cq;
					cgq = &cq->next;
					*cq = *cg;
					++cq;
					}
				}
			}
		*cgq = 0;
		}

	q = (dyad **)Malloc(nv*sizeof(dyad *));
	qe = q + nv;
	objadj = dsort(S, T, (ograd **)q, cgp, ogp, arrays);

	nelq = nz = nz1 = 0;
	/* In pass 0, we the count nonzeros in the lower triangle. */
	/* In pass 1, we compute the lower triangle and use column dispatch */
	/* (via the cdisp array) to copy the strict lower triangle to the */
	/* strict upper triangle.  This ensures symmetry. */
	for(pass = 0; pass < 2; pass++) {
		if (pass) {
			nelq += nelq - nz1;
			if (!nelq || !arrays)
				break;
			free(q);
			delsq1 = delsq = (double *)Malloc(nelq*sizeof(real));
			rowq1 = rowq = (fint *)Malloc(nelq*sizeof(fint));
			colq1 = colq = (Fint *)Malloc((nv+2)*sizeof(Fint));
			nelq = ftn;
			delsq0 = delsq - ftn;
			rowq0 = rowq - ftn;
			q = (dyad **)Malloc(nv*(sizeof(dyad*)
						+ sizeof(dispatch *)
						+ sizeof(dispatch)));
			qe = q + nv;
			cdisp = (dispatch**) qe;
			cdisp0 = cdisp - ftn;
			memset(cdisp, 0, nv*sizeof(dispatch*));
			cd0 = (dispatch *)(cdisp + nv);
			}
		memset(q, 0, nv*sizeof(dyad *));

		for(d = T->Q; d; d = d->next) {
			og = d->Rq;
			og1 = d->Lq;
			i = og->varno;
			while(og1 && og1->varno < i)
				og1 = og1->next;
			if (og1) {
				q1 = q + i;
				*q1 = new_dyad(S, *q1, og, og1, 0);
				}
			og1 = d->Lq;
			i = og1->varno;
			while(og && og->varno < i)
				og = og->next;
			if (og) {
				q1 = q + i;
				*q1 = new_dyad(S, *q1, og1, og, 0);
				}
			}
		vmi = asl->i.vmap ? get_vminv_ASL((ASL*)asl) : 0;
		for(icol = 0, q1 = q; q1 < qe; ++icol, ++q1) {
		    if (pass) {
			*colq++ = nelq;
			for(cd = cdisp[icol]; cd; cd = cdnext) {
				cdnext = cd->next;
				s[i = cd->i]++;
				x[z[nz++] = i] = delsq0[cd->j++];
				if (cd->j < cd->jend) {
					cdp = cdisp0 + rowq0[cd->j];
					cd->next = *cdp;
					*cdp = cd;
					}
				}
			}
		    if ((d = *q1))
			do {
				og = d->Lq;
				og1 = d->Rq;
				t = og->coef;
				for(; og1; og1 = og1->next) {
					if (!s[i = og1->varno]++)
						x[z[nz++] = i] =
							t*og1->coef;
					else
						x[i] += t*og1->coef;
					}
				if ((og1 = og->next)) {
				  og2 = d->Rq;
				  while (og2->varno < og1->varno)
				    if (!(og2 = og2->next)) {
					while((og1 = og->next))
						og = og1;
					break;
					}
				  d->Rq = og2;
				  }
				d1 = d->next;
				if ((og = og->next)) {
					i = og->varno;
					if (pass) {
						og1 = d->Rq;
						while(og1->varno < i)
							if (!(og1 = og1->next))
								goto d_del;
						d->Rq = og1;
						}
					d->Lq = og;
					q2 = q + i;
					d->next = *q2;
					*q2 = d;
					}
				else {
 d_del:
					free_dyad(S, d);
					}
				}
				while((d = d1));
		if (nz) {
			if (pass) {
				if (nz > 1)
					qsortv(z, nz, sizeof(fint), lcmp, NULL);
				for(i = nz1 = 0; i < nz; i++) {
					if ((t = x[j = z[i]])) {
						*delsq++ = t;
						if (vmi)
							j = vmi[j];
						*rowq++ = j + ftn;
						nelq++;
						z[nz1++] = j;
						}
					s[j] = 0;
					}
				for(i = 0; i < nz1; i++)
				    if ((j = z[i]) > icol && x[j]) {
					cd0->i = icol;
					cd0->j = colq[-1] + i;
					cd0->jend = nelq;
					cdp = cdisp + j;
					cd0->next = *cdp;
					*cdp = cd0++;
					break;
					}
				nz = 0;
				}
			else {
				while(nz > 0) {
					s[i = z[--nz]] = 0;
					if (x[i]) {
						++nelq;
						if (i == icol)
							++nz1;
						}
					}
				}
			}
		    }
		}
	free(q);
	free_blocks(S);
	free(x);
	if (od && od->cg)
		M1record(od->cg);
	if (nelq) {
		if (arrays) {
			/* allow one more for obj. adjustment */
			*colq = colq[1] = nelq;
			*rowqp = rowq1;
			*colqp = colq1;
			*delsqp = delsq1;
			}
		nelq -= ftn;
		}
	if (arrays) {
		if (od) {
			od->opify = 0; /*qp_opify_ASL;*/
			if ((t = od->c12) != 1.)
				for(i = 0; i < nelq; ++i)
					delsq1[i] *= t;
			objadj = t*objadj + od->c0a;
			for(i = 0, cg = *cgp; cg; cg = cg->next)
				++i;
			ogp = Ograd + co0;
			og2 = i ? (ograd*)M1alloc(i*sizeof(ograd)) : 0;
			for(cg = *cgp; cg; cg = cg->next) {
				*ogp = og = og2++;
				ogp = &og->next;
				og->varno = cg->varno;
				og->coef = t*cg->coef;
				}
			*ogp = 0;
			}
		else if (cgp && objadj != 0.) {
			if (Urhsx) {
				L = LUrhs + co;
				U = Urhsx + co;
				}
			else {
				L = LUrhs + 2*co;
				U = L + 1;
				}
			if (*L > negInfinity)
				*L -= objadj;
			if (*U < Infinity)
				*U -= objadj;
			objadj = 0;
			}
		if (co0 >= 0)
			asl->i.objconst[co0] = objadj;
		if (p)
			memset(p, 0, sizeof(ps_func));
		else {
			e[0] = nOPRET;
			e[1] = -1;
			}
		}
	return nelq;
	}


 Fints
nqpcheck_ASL(ASL *asl, int co, fint **rowqp, Fint **colqp, real **delsqp)
{
	Fints rv = mqpcheck_ASL(asl, co, rowqp, colqp, delsqp);
	if (rowqp && *rowqp) {
		M1record(*delsqp);
		M1record(*rowqp);
		M1record(*colqp);
		}
	return rv;
	}
