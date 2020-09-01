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

/* After fg_read() or pfgh_read(), degree(asl, co, pv) determines whether
   the indicated objective (co >= 0, 0 <= co < n_obj) or constraint
   (co < 0 indicates constraint 1 - co, 0 <= 1-co < n_con) is linear,
   quadratic, or higher-order.  In contrast to nqpcheck() and variants,
   this routine does not attempt to compute the number of nonzeros in
   the Hessian, so it can run much faster.  Possible return values:

	-1 ==> invalid co value
	 0 ==> constant
	 1 ==> linear
	 2 ==> quadratic
	 3 ==> more general nonlinear

   The pv argument (of type void **pv) should be NULL if degree_ASL(...)
   is to be called just once (e.g., to see if the objective is quadratic).
   If multiple calls are expected (e.g., for objectives and constraints),
   it may save time to use the pattern

	void *v = 0;
	for(...) { ... degree(asl, co, &v); ... }
	if (v) free(v);
*/

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

 static int
kind1(int *o, char *w, int lin)
{
	int i, j, k, m, n;

	if (!o)
		return 0;
 top:
	switch(*o) {

	  case OPRET:
		if ((i = o[1]) < 0 || !(i = w[i]))
			i = lin;
		return i;

	  case OPRETB:
		++o;
		goto top;

	  case nOPPLUS:
	  case nOPMINUS:
		if ((i = o[2]) < 0) {
			i = 0;
			if ((j = o[3]) >= 0)
				i = w[j];
			}
		else {
			i = w[i];
			if ((j = o[3]) >= 0 && (j = w[j]) > i)
				i = j;
			}
		w[o[1]] = i;
		o += 4;
		goto top;

	  case nOPUMINUS:
		w[o[1]] = w[o[2]];
		o += 3;
		goto top;

	  case nOPMULT:
		if ((i = o[2]) < 0)
			i = 0;
		else
			i = w[i];
		if ((j = o[3]) < 0)
			j = 0;
		else
			j = w[j];
		if ((i += j) > 3)
			i = 3;
		w[o[1]] = i;
		o += 4;
		goto top;

	  case OPDIV0:
		i = 0;
		goto nodiv;
	  case nOPDIV1:
	  case OPDIV2:
		if ((j = o[3]) >= 0)
			i = 3;
		else if ((i = o[2]) < 0)
			i = 0;
		else
			i = w[i];
 nodiv:
		w[o[1]] = i;
		o += 4;
		goto top;

	  case nOPSUMLIST:
		n = o[2];
		j = 3;
		m = 0;
		for(k = j + n; j < k; ++j) {
			if ((i = o[j]) >= 0) {
				i = w[i];
				if (m < i) {
					m = i;
					if (m >= 3)
						break;
					}
				}
			}
		w[o[1]] = m;
		o += k;
		goto top;

	  case OP2POW0:
		i = 0;
		goto nopow1;
	  case OP2POW1:
		if ((i = o[2]) < 0)
			i = 0;
		else {
			i = w[i];
			if ((i += i) > 3)
				i = 3;
			}
 nopow1:
		w[o[1]] = i;
		o += 3;
		goto top;

	  case OPGOTO:
	  case OPNEXTBLK:
		o = *(int**)(o+1);
		goto top;

#ifdef X64_bit_pointers
	  case OPGOTOalign:
	  case OPNEXTBLKalign:
		o = *(int**)(o+2);
		goto top;
#endif
	  default:
		break;
	  }
	return 3; /* nonlinear */
	}

 static int
kind2(int *o, char *w, int lin)
{
	int i, j, k, m, n;

	if (!o)
		return lin;
 top:
	switch(*o) {

	  case OPRET:
		if ((i = o[1]) < 0 || !(i = w[i]))
			i = lin;
		return i;

	  case OPPLUS0:
	  case OPMINUS0:
	  case OPMULT0:
	  case OPDIV0:
		w[o[1]] = 0;
		o += 4;
		goto top;

	  case OPPLUS10:
	  case OPPLUS01:
	  case OPPLUS2:
	  case OPMINUS10:
	  case OPMINUS01:
	  case OPMINUS2:
		if ((i = o[3]) < 0) {
			i = 0;
			if ((j = o[4]) >= 0)
				i = w[j];
			}
		else {
			i = w[i];
			if ((j = o[4]) >= 0 && (j = w[j]) > i)
				i = j;
			}
		w[o[2]] = i;
		o += 5;
		goto top;

	  case OPUMINUS0:
	  case OP_2POW0:
		w[o[1]] = 0;
		o += 3;
		goto top;

	  case OPUMINUS1:
		w[o[2]] = w[o[3]];
		o += 4;
		goto top;

	  case OPMULT10:
	  case OPMULT01:
	  case OPMULT2:
		if ((i = o[3]) < 0)
			i = 0;
		else
			i = w[i];
		if ((j = o[4]) < 0)
			j = 0;
		else
			j = w[j];
		if ((i += j) > 3)
			i = 3;
		w[o[2]] = i;
		o += 5;
		goto top;

	  case OPDIV10:
	  case OPDIV01:
	  case OPDIV2:
		if ((j = o[4]) >= 0)
			i = 3;
		else if ((i = o[3]) < 0)
			i = 0;
		else
			i = w[i];
		w[o[2]] = i;
		o += 5;
		goto top;

	  case OPSUMLIST0:
		m = 0;
		k = o[2] + 3;
		goto sum0;

	  case OPSUMLIST1:
		++o;
		n = o[2];
		j = 3;
		m = 0;
		for(k = j + n; j < k; ++j) {
			if ((i = o[j]) >= 0) {
				i = w[i];
				if (m < i) {
					m = i;
					if (m >= 3)
						break;
					}
				}
			}
 sum0:
		w[o[1]] = m;
		o += k;
		goto top;

	  case OP_2POW1:
		if ((i = o[3]) < 0)
			i = 0;
		else {
			i = w[i];
			if ((i += i) > 3)
				i = 3;
			}
		w[o[2]] = i;
		o += 4;
		goto top;

	  case OP_GOTO:
	  case OPGOTO2:
	  case OP_NEXTBLK:
	  case OPGOTOF:
	  case OPGOTOF2:
	  case OPGOTOF2n:
		o = *(int**)(o+1);
		goto top;

#ifdef X64_bit_pointers
	  case OP_GOTOalign:
	  case OPGOTO2align:
	  case OP_NEXTBLKalign:
	  case OPGOTOFalign:
	  case OPGOTOF2align:
	  case OPGOTOF2nalign:
		o = *(int**)(o+2);
		goto top;
#endif
	  default:
		break;
	  }
	return 3; /* nonlinear */
	}

 static void
cecomp(ASL *asl, int i, int n, char *w)
{
	cexp2 *ce, *ce1;
	char *wd;
	int j, k, *dvsp0, *ndvsp;

	wd = w + 4*asl->i.defvar0;
	ce = ((ASL_pfgh*)asl)->I.cexps2_;
	dvsp0 = ((ASL_pfgh*)asl)->P.dvsp0;
	ndvsp = ((ASL_pfgh*)asl)->P.ndvsp;
	for(; i < n; ++i) {
		if ((k = ndvsp[i])) {
			for(j = dvsp0[i], k += j; j < k; ++j)
				wd[4*j] = kind2(ce[j].o.e, w, 0);
			}
		ce1 = ce + i;
		wd[4*i] = kind2(ce1->o.e, w, ce1->lp ? 1 : 0);
		}
	}

 int
degree_ASL(ASL *asl, int co, void **pv)
{
	ASL_fg *asl1;
	ASL_pfgh *asl2;
	Objrep *od, **pod;
	cde  *cd1;
	cde2 *cd2;
	cexp  *ce0, *ce01;
	cexp1 *ce1, *ce11;
	cexp2 *ce2, *ce21;
	cgrad *cg;
	char *w, *wd;
	int ak, b1, *c1, *cm, i, nb, nbc, nbco, nc1, nd, nv, nv0, *o, rv, sf;
	linarg *la;
	ograd *og;
	ps_func *f;
	psb_elem *b, *be;
	psg_elem *g, *ge;
	size_t L;
	uint af1, afn, j, j0;

	ak = asl->i.ASLtype;
	if (ak < ASL_read_f || ak > ASL_read_pfgh)
		badasl_ASL(0,0,"degree_ASL");
	if (co >= n_obj || co < -n_con)
		return -1;
	sf = ak <= ASL_read_fg ? 1 : 4;
	nd = asl->i.defvar0;
	nb = comb;
	nbc = combc;
	nbco = nbc + como;
	nv = n_var;
	nv0 = asl->i.defvar0;
	b1 = nv0 + nbco;
	L = asl->i.wlen;
	if (pv && (w = *(char**)pv))
		wd = w + sf*nd;
	else {
		w = (char*)Malloc(L);
		memset(w, 1, sf*nv);
		wd = w + sf*nd;
		if (nb > 0)
			wd[0] = 15;
		if (nbc > nb)
			wd[sf*nb] = 15;
		if (nbco > nbc)
			wd[sf*nbc] = 15;
		if (pv)
			*pv = w;
		if (ak > ASL_read_fg && (la = ((ASL_pfgh*)asl)->P.lalist)) do
			w[4*la->u.v] = 1;
			while((la = la->lnext));
		}
	f = 0;
	od = 0;
	if (ak <= ASL_read_fg) {
		asl1 = (ASL_fg*)asl;
		asl2 = 0;
		ce0 = asl1->I.cexps_;
		ce1 = asl1->I.cexps1_;
		ce2 = 0;	/* avoid spurious warning */
		}
	else {
		asl1 = 0;
		asl2 = (ASL_pfgh*)asl;
		ce0 = 0; /*won't be used*/
		ce1 = 0; /*won't be used*/
		ce2 = asl2->I.cexps2_;
		}
	c1 = 0; /* not needed */
	if (co >= 0) {
		if ((pod = asl->i.Or) && (od = pod[co])) {
			co = od->ico;
			goto use_Cgrad;
			}
		if (ak <= ASL_read_fg)  {
			cd1= &((ASL_fg*)asl)->I.obj_de_[co];
			o = cd1->o.e;
			af1 = cd1->af1;
			afn = cd1->afn;
			if (nb > 0 && wd[0] == 15) {
				for(i = 0; i < nb; ++i) {
					ce01 = ce0 + i;
					wd[i] = kind1(ce01->o.e, w, ce01->lp ? 1 : 0);
					}
				}
			if (nbco > nbc && wd[nbc] == 15) {
				for(i = nbc; i < nbco; ++i) {
					ce01 = ce0 + i;
					wd[i] = kind1(ce01->o.e, w, ce01->lp ? 1 : 0);
					}
				}
			}
		else {
			c1 = asl2->i.o_cexp1st_;
			cd2 = &((ASL_pfgh*)asl)->I.obj2_de_[co];
			f = &((ASL_pfgh*)asl)->P.ops[co];
			o = cd2->o.e;
			af1 = cd2->af1;
			afn = cd2->afn;
			if (nb > 0 && wd[0] == 15)
				cecomp(asl, 0, nb, w);
			if (nbco > nbc && wd[4*nbc] == 15)
				cecomp(asl, nbc, nbco, w);
			}
		og = Ograd[co];
		cg = 0;
		}
	else {
		co = -1 - co;
		if ((cm = asl->i.cmap))
			co = cm[co];
 use_Cgrad:
		if (ak <= ASL_read_fg) {
			cd1 = &((ASL_fg*)asl)->I.con_de_[co];
			o = cd1->o.e;
			af1 = cd1->af1;
			afn = cd1->afn;
			if (nb > 0 && wd[0] == 15) {
				for(i = 0; i < nbc; ++i) {
					ce01 = ce0 + i;
					wd[i] = kind1(ce01->o.e, w, ce01->lp ? 1 : 0);
					}
				}
			else if (nbc > nb && wd[nb] == 15) {
				for(i = nb; i < nbc; ++i) {
					ce01 = ce0 + i;
					wd[i] = kind1(ce01->o.e, w, ce01->lp ? 1 : 0);
					}
				}
			}
		else {
			c1 = asl2->i.c_cexp1st_;
			cd2 = &((ASL_pfgh*)asl)->I.con2_de_[co];
			f = &((ASL_pfgh*)asl)->P.cps[co];
			o = cd2->o.e;
			af1 = cd2->af1;
			afn = cd2->afn;
			if (nb > 0 && wd[0] == 15)
				cecomp(asl, 0, nbc, w);
			else if (nbc > nb && wd[4*nb] == 15)
				cecomp(asl, nb, nbc, w);
			}
		cg = Cgrad[co];
		og = 0;
		}
	if (ak == ASL_read_f)
		goto linchk;
	else if (ak == ASL_read_fg) {
		if (afn) {
			afn += j = af1;
			do {
				ce11 = &ce1[j - b1];
				w[j] = kind1(ce11->o.e, w, ce11->lp ? 1 : 0);
				}
				while(++j < afn);
			}
		rv = kind1(o, w, 0);
		}
	else {
		if (c1 && (nc1 = c1[co+1] - (i = c1[co])) > 0)
			cecomp(asl, i, i + nc1, w);
		if (afn) {
			afn += j = af1;
			do {
				j0 = j - nv0;
				ce21 = &ce2[j0];
				w[4*j] = kind2(ce21->o.e, w, ce21->lp ? 1 : 0);
				}
				while(++j < afn);
			}
		rv = 0;
		if ((b = f->pi.b))
			for(be = f->pi.be; b < be; ++b) {
				i = kind2(b->o.e, w, 0);
				if (rv < i && (rv = i) >= 3)
					goto ret;
				}
		if ((g = f->g)) {
			rv = 0;
			for(ge = f->ge; g < ge; ++g) {
				o = g->o;
				if (*o != OP2POW_g)
					goto ret3;
				switch(o[3]) {
					case OPRET:
						break;
					case OP_GOTO:
					case OPGOTO2:
					case OP_NEXTBLK:
						if (**(int**)(o+4) != OPRET)
							goto ret3;
						break;
#ifdef X64_bit_pointers
					case OP_GOTOalign:
					case OPGOTO2align:
					case OP_NEXTBLKalign:
						if (**(int**)(o+5) != OPRET)
							goto ret3;
						break;
#endif
					default:
						goto ret3;
					}
				if ((b = g->pi.b)) {
					if (b->o.f || !(o = b->o.e))
						goto ret3;
					switch(o[0]) {
						case OPRET:
							break;
						case OP_GOTO:
						case OPGOTO2:
							if (**(int**)(o+1) != OPRET)
								goto ret3;
							break;
#ifdef X64_bit_pointers
						case OP_GOTOalign:
						case OPGOTO2align:
							if (**(int**)(o+2) != OPRET)
								goto ret3;
							break;
#endif
						default:
							goto ret3;
						}
					switch(w[o[1]]) {
					  case 1: rv = 2;
					  case 0: break;
					  default: goto ret3;
					  }
					}
				if (g->L)
					rv = 2;
				}
			}
		}
	if (rv > 3) {
 ret3:
		rv = 3;
		}
	else if (rv == 0) {
 linchk:
		while(og) {
			if (og->coef) {
				rv = 1;
				goto ret;
				}
			og = og->next;
			}
		while(cg) {
			if (cg->coef) {
				rv = 1;
				goto ret;
				}
			cg = cg->next;
			}
		}
 ret:
	if (!pv)
		free(w);
	return rv;
	}
