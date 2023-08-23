/*******************************************************************
Copyright (C) 2017, 2019 AMPL Optimization, Inc.; written by David M. Gay.

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
#define SKIP_NL2_DEFINES
#include "nlp2.h"
#include "opno.hd"
#include "opno2.h"
#include "opcode.hd"

#ifdef __cplusplus
extern "C" {
#endif

 struct
expr {
	int opno;
	int opcl;
	struct expr *L, *R;
	};

 typedef struct
exprf {
	int opno;
	int opcl;
	int fi;
	int n;
	expr *args[1];	/* really args[n] */
	} exprf;

 typedef struct
exprh {
	int opno;
	int opcl;
	int len;
	char sym[4]; /* really len + 1 */
	} exprh;

 typedef struct
expri {
	int opno;
	int opcl;
	expr *e, *T, *F;
	} expri;

 typedef struct
exprn {
	int opno;
	int opcl;
	real v;
	} exprn;

 typedef struct
exprp {
	int opno;
	int opcl;
	struct expr **pL, **pR;
	} exprp;

 typedef struct
exprv {
	int opno;
	int opcl;
	size_t varno;
	} exprv;

 typedef struct
Cexp {
	expr	*e;
	linpart	*L;
	} Cexp;

enum {	Unop = 1,
	Binop = 2,
	Varg = 3,
	Plterm = 4,
	EIf = 5,
	Esum = 6,
	Efunc = 7,
	Ehol = 8,
	Enum = 9,
	Evar = 10,
	nEblock = 4095
	};

 static exprn Two = { OPNUM, Enum, 2. };


 typedef struct
Eblock {
	union { struct Eblock *b;
		real align;
		} bnext;
	real Eb[nEblock];
	} Eblock;

typedef int Pf(FILE*, const char*, ...);

 typedef struct Staticfgw
{
	ASL *asl;
	FILE *nl_;
	Pf *pf_;
	Eblock *Ebusy;
	arglist *al;
	cexp1 *cexps1_;
	expr **w;
	real *Enext, *Elast;
	jmp_buf wjb;
	int com1off;
	int nv0;
	int keb;
	} Staticfgw;

 static expr*
Ealloc(Staticfgw *S, size_t len)
{
	ASL *asl;
	Eblock *eb;
	expr *rv;
	int k;
	size_t n;

	len = (len + sizeof(real) - 1) & ~(sizeof(real)-1);
	n = len / sizeof(real);
	if (S->Elast < S->Enext + n) {
		asl = S->asl;
		if (!(k = S->keb))
			S->keb = k = htcl(sizeof(Eblock));
		eb = 0;
		if (len > sizeof(eb->Eb))
			k = htcl(len + sizeof(eb->bnext));
		eb = (Eblock*)new_mblk(k);
		eb->bnext.b = S->Ebusy;
		S->Ebusy = eb;
		S->Enext = eb->Eb;
		S->Elast = (real*)((char*)eb + (sizeof(void*) << k));
		}
	rv = (expr*)S->Enext;
	S->Enext += n;
	return rv;
	}

 static void
Efree(Staticfgw *S)
{
	ASL *asl;
	Eblock *eb, *eb1;

	asl = S->asl;
	for(eb1 = S->Ebusy; (eb = eb1);) {
		eb1 = eb->bnext.b;
		Del_mblk_ASL(asl, eb);
		}
	S->Ebusy = 0;
	S->Elast = S->Enext = 0;
	}

 static int
aprintf(FILE *fd, const char *fmt, ...)
{
	char *s;
	char buf[32];
	va_list ap;
	int i, j;
	double x;
	int rc = 0;

	va_start(ap, fmt);
	if (*fmt != '%')
		rc++;
	for(;;) {
		for(;;) {
			switch(i = *fmt++) {
				case 0:	  goto done;
				case '%': break;
				default:  putc(i,fd);
					  continue;
				}
			break;
			}
		rc++;
		switch(*fmt++) {
			case 'c':
				i = va_arg(ap, int);
				putc(i,fd);
				continue;
			case 'd':
				i = va_arg(ap, int);
				if (i < 0) {
					putc('-',fd);
					i = -i;
					}
				s = buf;
				do {
					j = i / 10;
					*s++ = i - 10*j + '0';
					}
					while((i = j));
				do {
					i = *--s;
					putc(i,fd);
					}
					while(s > buf);
				continue;
			case '.':
				while(*fmt++ != 'g');
			case 'g':
				x = va_arg(ap, double);
				i = g_fmt(s = buf, x);
				goto have_s;
			case 's':
				s = va_arg(ap, char*);
 have_s:
				while((i = *s++))
					putc(i,fd);
				continue;
			default:
				fprintf(Stderr, "aprintf bug: unexpect fmt %s\n",
					fmt-1);
				exit(1);
			}
		}
 done:
	va_end(ap);
	return rc;
	}

#ifndef Int
#define Int Long
#endif

 static int
bprintf(FILE *fd, const char *fmt, ...)
{
	char *s;
	int i, rc;
	size_t len;
	union U { double x; short sh; Long L; Int i; char c; } u;
	va_list ap;

	va_start(ap, fmt);

	rc = 0;
	len = 0; /* silence buggy "not-initialized" warning */
	if ((i = *fmt) != '%') {
		fmt++;
#ifdef DMG
		if (i != 'i')
#endif
		{
		u.c = i;
		fwrite(&u.c, 1, 1, fd);
		rc++;
		}}

	for(;;) {
		while(*fmt == ' ')
			fmt++;
		if (*fmt++ != '%')
			break;
		switch(*fmt++) {
			case 'c':
				u.c = va_arg(ap, int);
				len = 1;
				break;
			case 'd':
				u.i = va_arg(ap, int);
				len = sizeof(Int);
				break;
			case '.':
				while(*fmt++ != 'g');
			case 'g':
				u.x = va_arg(ap, double);
				len = sizeof(double);
				break;
			case 'h':
				u.sh = va_arg(ap, int);
				len = sizeof(short);
				if (*fmt == 'd')
					fmt++;
				break;
			case 'l':
				u.L = (Long)va_arg(ap, long);
				len = sizeof(Long);
				if (*fmt == 'd')
					fmt++;
				break;
			case 's':
				s = va_arg(ap, char*);
				u.i = strlen(s);
				fwrite((char *)&u.i, sizeof(Int), 1, fd);
				fwrite(s, u.i, 1, fd);
				goto s_written;
			default:
				fprintf(Stderr, "bprintf bug: unexpect fmt %s\n",
					fmt-1);
				exit(1);
			}
		fwrite((char *)&u.L, len, 1, fd);
 s_written:
		rc++;
		}
	va_end(ap);
	return rc;
	}

 static expr*
ewalk1(Staticfgw *S, int *o, int *ostop)
{
	arglist *al;
	char *s;
	expr *e, **pe, *rv, **w;
	exprf *ef;
	exprh *eh;
	expri *eif;
	exprn *en;
	exprp *ep;
	exprv *v;
	int *at, i, j, k, n, *o1, *o2, op, **pop;
	plterm *p;
	real *bs, *rp;

	w = S->w;
 top:
	if (o == ostop)
		return 0;
	switch(*o) {
	  case nOPRET:
		rv = w[o[1]];
		goto done;
	  case nOPPLUS:
		op = OPPLUS;
 bin:
		e = Ealloc(S, sizeof(expr));
		e->opno = op;
		e->opcl = Binop;
		e->L = w[o[2]];
		e->R = w[o[3]];
		w[o[1]] = e;
		o += 4;
		goto top;
	  case nOPMINUS:
		op = OPMINUS;
		goto bin;
	  case nOPMULT:
		op = OPMULT;
		goto bin;
	  case nOPDIV0:
	  case nOPDIV1:
	  case nOPDIV2:
	  case nOPDIV3:
		op = OPDIV;
		goto bin;
	  case nOPREM0:
	  case nOPREM1:
		op = OPREM;
		goto bin;
	  case nOPPOW0:
	  case nOPPOW1:
	  case nOPPOW2:
	  case nOPPOW3:
	  case nOPPOW4:
		op = OPPOW;
		goto bin;
	  Opalign(OPLESSalign)
	  case nOPLESS0:
	  case nOPLESS1:
		op = OPLESS;
		goto bin;
	  Opalign(OPMINLIST1align)
	  case nOPMINLIST0:
	  case nOPMINLIST1:
		op = MINLIST;
 minmax:
		ep = (exprp*)Ealloc(S, sizeof(exprv) + o[2]*sizeof(expr*));
		ep->opno = op;
		ep->opcl = Varg;
		ep->pL = pe = (expr**)(ep + 1);
		k = o[2] + 3;
		i = 2;
		k = i + o[2];
		while(++i < k)
			*pe++ = w[o[i]];
		ep->pR = pe;
		*(exprp**)&w[o[1]] = ep;
		o += k;
		goto top;
	  Opalign(OPMAXLIST1align)
	  case nOPMAXLIST1:
	  case nOPMAXLIST0:
		op = MAXLIST;
		goto minmax;
	  case nFLOOR:
		op = FLOOR;
 un:
		e = Ealloc(S, sizeof(expr));
		e->opno = op;
		e->opcl = Unop;
		e->L = w[o[2]];
		e->R = 0;
		w[o[1]] = e;
		o += 3;
		goto top;
	  case nCEIL:
		op = CEIL;
		goto un;
	  case nOPABS0:
	  case nOPABS1:
		op = ABS;
		goto un;
	  case nOPUMINUS:
		op = OPUMINUS;
		goto un;
	  Opalign(OPORalign)
	  case nOPOR:
		op = OPOR;
 andor:
		e = Ealloc(S, sizeof(expr));
		e->opno = op;
		e->opcl = Binop;
		e->L = w[o[2]];
		pop = (int**)&o[3];
		if (*pop == o + 4) {
			e->R = e->L;
			o += 4;
			goto top;
			}
		ewalk1(S, pop[1], o = pop[0] - 3);
		e->R = w[o[2]];
		w[o[1]] = e;
		o += 3;
		goto top;
	  Opalign(OPANDalign)
	  case nOPAND:
		op = OPAND;
		goto andor;
	  case nOPLT:
		op = LT;
		goto bin;
	  case nOPNOTATMOST:
		op = OPNOTATMOST;
		goto bin;
	  case nOPLE:
		op = LE;
		goto bin;
	  case nOPATLEAST:
		op = OPATLEAST;
		goto bin;
	  case nOPEQ:
		op = EQ;
		goto bin;
	  case nOPEXACTLY:
		op = OPEXACTLY;
		goto bin;
	  case nOPGE:
		op = GE;
		goto bin;
	  case nOPATMOST:
		op = OPATMOST;
		goto bin;
	  case nOPGT:
		op = GT;
		goto bin;
	  case nOPNOTATLEAST:
		op = OPNOTATLEAST;
		goto bin;
	  case nOPNE:
		op = NE;
		goto bin;
	  case nOPNOTEXACTLY:
		op = OPNOTEXACTLY;
		goto bin;
#ifdef X64_bit_pointers
	  case OPIMPELSE_align:
		++o;
		goto more_impelse;
	  case OPIFnl0align:
		++o;
		goto more_if0;
#endif
	  case nOPIMPELSE:
 alignarg(more_impelse:)
		op = OPIMPELSE;
		pop = (int**)&o[5];
		goto more_opcond;
	  case nOPIFnl0:	/* and OPIFSYM */
 alignarg(more_if0:)
		pop = (int**)&o[5];
 more_if:
		op = OPIFnl;
 more_opcond:
		eif = (expri*)Ealloc(S, sizeof(expri));
		eif->opcl = EIf;
		eif->e = w[o[2]];
		ewalk1(S, pop[0], pop[2]);
		e = eif->T = w[o[3]];
		if (e->opcl == Ehol || (e->opcl == EIf && e->opno == OPIFSYM))
			op = OPIFSYM;
		eif->opno = op;
		ewalk1(S, pop[1], pop[3]);
		eif->F = w[o[4]];
		w[o[1]] = (expr*)eif;
		o = pop[4];
		goto top;
	  Opalign(OPIFnl1align)
	  case nOPIFnl1:
		pop = (int**)&o[6];
		goto more_if;
	  case nOP_tanh0:
	  case nOP_tanh1:
		op = OP_tanh;
		goto un;
	  case nOP_tan0:
	  case nOP_tan1:
		op = OP_tan;
		goto un;
	  case nOP_sqrt0:
	  case nOP_sqrt1:
		op = OP_sqrt;
		goto un;
	  case nOP_sinh0:
	  case nOP_sinh1:
		op = OP_sinh;
		goto un;
	  case nOP_sin0:
	  case nOP_sin1:
		op = OP_sin;
		goto un;
	  case nOP_log100:
	  case nOP_log101:
		op = OP_log10;
		goto un;
	  case nOP_log0:
	  case nOP_log1:
		op = OP_log;
		goto un;
	  case nOP_exp:
		op = OP_exp;
		goto un;
	  case nOP_cosh0:
	  case nOP_cosh1:
		op = OP_cosh;
		goto un;
	  case nOP_cos0:
	  case nOP_cos1:
		op = OP_cos;
		goto un;
	  case nOP_atanh0:
	  case nOP_atanh1:
		op = OP_atanh;
		goto un;
	  case nOP_atan20:
	  case nOP_atan21:
	  case nOP_atan22:
	  case nOP_atan23:
		op = OP_atan2;
		goto bin;
	  case nOP_atan0:
	  case nOP_atan1:
		op = OP_atan;
		goto un;
	  case nOP_asinh0:
	  case nOP_asinh1:
		op = OP_asinh;
		goto un;
	  case nOP_asin0:
	  case nOP_asin1:
		op = OP_asin;
		goto un;
	  case nOP_acosh0:
	  case nOP_acosh1:
		op = OP_acosh;
		goto un;
	  case nOP_acos0:
	  case nOP_acos1:
		op = OP_acos;
		goto un;
	  case nOPSUMLIST:
		op = OPSUMLIST;
 more_sumlist:
		n = o[2];
		ep = (exprp*)Ealloc(S, sizeof(exprp) + n*sizeof(expr*));
		ep->opno = op;
		ep->opcl = Esum;
		pe = ep->pL = (expr**)(ep + 1);
		ep->pR = pe + n;
		j = 3;
		for(k = j + n; j < k; ++j)
			*pe++ = w[o[j]];
		w[o[1]] = (expr*)ep;
		o += k;
		goto top;
	  case nOPintDIV:
		op = OPintDIV;
		goto bin;
	  case nOPprecision:
		op = OPprecision;
		goto bin;
	  case nOPround:
		op = OPround;
		goto bin;
	  case nOPtrunc:
		op = OPtrunc;
		goto bin;
	  case nOPCOUNT:
		op = OPCOUNT;
		goto more_sumlist;
	  case nOPNUMBEROF:
		op = OPNUMBEROF;
		goto more_sumlist;
	  case nOPNUMBEROFs:	/* TEMPORARY */
		op = OPNUMBEROFs;
		goto more_sumlist;
#ifdef X64_bit_pointers
	  case OPPLTERM0align:
	  case OPPLTERM1align:
		++o;
#endif
	  case nOPPLTERM0:
	  case nOPPLTERM1:
		i = o[1];
		v = (exprv*)w[o[2]];
		o += 3;
		p = (plterm*)o;
		n = p->n;
		bs = p->bs;
		o = (int*)&bs[2*n-1];
		w[i] = e = Ealloc(S, sizeof(expr));
		e->opno = OPPLTERM;
		e->opcl = Plterm;
		e->L = (expr*)p;
		e->R = (expr*)v;
		goto top;
	  case nOPANDLIST:
		op = ANDLIST;
		goto more_sumlist;
	  case nOPORLIST:
		op = ORLIST;
		goto more_sumlist;
	  case nOP_IFF:
		op = OP_IFF;
		goto bin;
	  case nOPALLDIFF:
		op = OPALLDIFF;
		goto more_sumlist;
	  case nOPSOMESAME:
		op = OPSOMESAME;
		goto more_sumlist;
	  case OP2POW0:
	  case OP2POW1:
		e = Ealloc(S, sizeof(expr));
		e->opno = OPPOW;
		e->opcl = Binop;
		e->L = w[o[2]];
		e->R = (expr*)&Two;
		w[o[1]] = e;
		o += 3;
		goto top;
#ifdef X64_bit_pointers
	  case OPCPOWalign:
		rv = w[o[1]];
		rp = (real*)&o[2];
		goto more_CPOW;
#endif
	  case nOPCPOW:
		rp = (real*)&o[1];
		rv = w[o[5]];
 alignarg(more_CPOW:)
		en = (exprn*)Ealloc(S, sizeof(exprn));
		en->opno = OPNUM;
		en->opcl = Enum;
		en->v = rp[0];
		w[o[1]] = e = Ealloc(S, sizeof(expr));
		e->opno = OPPOW;
		e->opcl = Binop;
		e->L = (expr*)en;
		e->R = rv;
		o += 7;
		goto top;

	  case OPFUNCALL0:
	  case OPFUNCALL1:
		k = o[1];
		al = S->al + o[2];
		n = al->n;
		at = al->at;
		ef = (exprf*)Ealloc(S, sizeof(exprf) + (n-1)*sizeof(expr*));
		ef->opno = OPFUNCALL;
		ef->opcl = Efunc;
		pe = ef->args;
		ef->n = n;
		ef->fi = o[3];
		o1 = o + 4;
		o2 = o1 + al->nr;
		for(i = 0; i < n; ++i) {
			if ((j = at[i]) < 0)
				e = w[o2[-(j + 1)]];
			else
				e = w[o1[j]];
			pe[i] = e;
			}
		w[k] = (expr*)ef;
		o = o1 + n;
		goto top;

	  case nOPHOL:
		s = (char*)(o+3);
		n = (int)strlen(s);
		eh = (exprh*)Ealloc(S, sizeof(exprh) + n - 3);
		eh->opno = OPHOL;
		eh->opcl = Ehol;
		eh->len = n;
		strcpy(eh->sym, s);
		w[o[1]] = (expr*)eh;
		o += o[2];
		goto top;

	  /*case nOPVARVAL: accessed directly */

	  case OPCOPY:
	  case OPCOPYSYM:
		w[o[1]] = w[o[2]];
		o += 3;
		goto top;

	  case OPGOTO:
	  case OPGOTO2:
	  case OP_NEXTBLK:
	  case OPGOTOF:
	  case OPGOTOF2:
	  case OPGOTOF2n:
		o = *(int**)(o+1);	/* for chaining blocks of code */
		goto top;

#ifdef X64_bit_pointers
	  case OPGOTOalign:
	  case OPGOTO2align:
	  case OP_NEXTBLKalign:
	  case OPGOTOFalign:
	  case OPGOTOF2align:
	  case OPGOTOF2nalign:
		o = *(int**)(o+2);
		goto top;
#endif

	  default:
		fprintf(Stderr, "\nUnexpected opno %d in eval1_ASL()\n", *o);
		fflush(Stderr);
		exit(1);
		rv = 0; /* not reached */
	  }
 done:
	return rv;
	}

#define pf S->pf_
#define nl S->nl_

 static void
eput(Staticfgw *S, expr *e)
{
	expr **ap, **ape;
	exprf *ef;
	exprh *eh;
	expri *eif;
	exprn *en;
	exprp *ep;
	exprv *v;
	exprp *va;
	int i, nop;
	plterm *p;
	real *r, *re;

 top:
	nop = e->opno;
	if ((i = e->opcl) < 7)
		(*pf)(nl, "o%d\n", nop);
	switch(i) {
	 case Unop:
		e = e->L;
		goto top;
	 case Binop:
		eput(S, e->L);
		e = e->R;
		goto top;
	 case Varg:
		va = (exprp *)e;
		ap = va->pL;
		ape = va->pR;
		(*pf)(nl, "%d\n", (int)(ape-ap));
		while(ap < ape)
			eput(S, *ap++);
		break;
	 case Plterm:
		p = (plterm*)e->L;
		(*pf)(nl, "%d\n", p->n);
		r = p->bs;
		re = r + 2*p->n - 1;
		while(r < re)
			(*pf)(nl, "n%g\n", *r++);
		e = e->R;
		goto top;
	 case EIf:
		eif = (expri*)e;
		eput(S, eif->e);
		eput(S, eif->T);
		e = eif->F;
		goto top;
	 case Esum:
		ep = (exprp*)e;
		ap = ep->pL;
		ape = ep->pR;
		(*pf)(nl, "%d\n", (int)(ape - ap));
		while(ap < ape)
			eput(S, *ap++);
		break;
	 case Efunc:
		ef = (exprf*)e;
		(*pf)(nl, "f%d %d\n", ef->fi, ef->n);
		ap = ef->args;
		ape = ap + ef->n;
		while(ap < ape)
			eput(S, *ap++);
		break;
	 case Ehol:
		eh = (exprh*)e;
		(*pf)(nl, "h%d:%s\n", eh->len, eh->sym);
		break;
	 case Enum:
		en = (exprn*)e;
		(*pf)(nl, "n%g\n", en->v);
		break;
	 case Evar:
		v = (exprv*)e;
		(*pf)(nl, "v%d\n", (int)v->varno);
		break;
	 default:
		fprintf(Stderr, "fg_write: unexpected type %d in eput.\n", i);
		longjmp(S->wjb, 1);
	 }
	}

 static void
Ewput(Staticfgw *S, int *o)
{
	expr *e;

	if ((e = ewalk1(S, o, 0))) {
		eput(S, e);
		Efree(S);
		}
	else
		(*pf)(nl, "n%g\n", 0.);
	}

/*#define offset_of(t,c) ((size_t)(char *)&((t*)0)->c)*/

 static void
coput(Staticfgw *S, int c, cde *de, int n, int *cexp1st, char *ot, int voff,
	int nn, real *oc, char *Not)
{
	cexp1 *ce;
	int i, i1, j, je, k;
	lincoef *Lc, *Lce;
	linpart *L;
	real t;

	if (cexp1st) {
		j = cexp1st[0];
		ce = S->cexps1_ + j;
		}
	else /* silence buggy "not-initialized" warnings */
		{ ce = 0; j = 0; }

	for(i = 0; i < n; i++) {
		if (cexp1st) {
			je = cexp1st[i1 = i + 1];
			i1 += voff;
			while(j < je) {
				k = S->com1off + j++;
				if ((L = ce->lp)) {
					(*pf)(nl, "V%d %d %d\n", k, L->n, i1);
					Lc = L->lc;
					for(Lce = Lc + L->n; Lc < Lce; ++Lc)
						(*pf)(nl, "%d %g\n", (int)Lc->varno, Lc->coef);
					}
				else
					(*pf)(nl, "V%d %d %d\n", k, 0, i1);
				Ewput(S, ce->o.e);
				ce++;
				}
			}
		if (ot)
			(*pf)(nl, "%c%d %d\n", c, i, ot[i]);
		else
			(*pf)(nl, "%c%d\n", c, i);
		Ewput(S, de[i].o.e);
		}
	t = 0.;
	for(n += nn; i < n; i++) {
		if (ot) {
			(*pf)(nl, "%c%d %d\n", c, i, Not ? *Not++ : 0);
			if (oc)
				t = *oc++;
			}
		else
			(*pf)(nl, "%c%d\n", c, i);
		(*pf)(nl, "n%g\n", t);
		}
	}

#undef pf
#undef nl

 static void
iguess(Pf *pf, FILE *nl, int c, real *x, char *havex, int n, int nn, real *y)
{
	int i, k;

	if (n + nn <= 0)
		return;
	i = k = 0;
	if (x) {
		if (havex) {
			while(i < n)
				if (havex[i++])
					k++;
			}
		else {
			while(i < n)
				if (x[i++])
					k++;
			}
		}
	if (y)
		for(i = 0; i < nn; i++)
			if (y[i])
				k++;
	if (!k)
		return;
	(*pf)(nl, "%c%d\n", c, k);
	if (x) {
		if (havex) {
			for(i = 0; i < n; i++)
				if (havex[i])
					(*pf)(nl, "%d %g\n", i, x[i]);
			}
		else {
			for(i = 0; i < n; i++)
				if (x[i])
					(*pf)(nl, "%d %g\n", i, x[i]);
			}
		}
	if (y) {
		for(i = 0; i < nn; i++)
			if (y[i])
				(*pf)(nl, "%d %g\n", i+n, y[i]);
		}
	}

 static void
br(Pf *pf, FILE *nl, int c, real *Lb, real *Ub, int n)
{
	int i;
	real L, U;

	if (n <= 0)
		return;
	if (c)
		(*pf)(nl, "%c\n", c);
	for(i = 0; i < n; i++) {
		L = *Lb++;
		U = Ub ? *Ub++ : *Lb++;
		if (L <= negInfinity)
			(*pf)(nl, U >= Infinity ? "3\n" : "1 %g\n", U);
		else
			(*pf)(nl, U >= Infinity ? "2 %g\n"
				: L == U ? "4 %g\n"
					 : "0 %g %g\n",
				L, U);
		}
	}

 static void
Gput(Pf *pf, FILE *nl, int c, int i, int n, ograd **ogp)
{
	ograd *og;
	int k;

	if (n <= 0)
		return;
	for(n += i; i < n; i++, ogp++) {
		if (!(og = *ogp))
			continue;
		k = 0;
		do k++;
			while((og = og->next));
		(*pf)(nl, "%c%d %d\n", c, i, k);
		for(og = *ogp; og; og = og->next)
			(*pf)(nl, "%d %g\n", og->varno, og->coef);
		}
	}

 static void
k2put(Pf *pf, FILE *nl, cgrad **cgp, int nc, int n, int k, int nnv,
	int nnc, ograd **ogp)
{
	cgrad *cg;
	ograd *og;
	int i, n1, *z;

	if (k) {
		n1 = n + nnv;
		z = (int*)Malloc(n1*sizeof(int));
		memset(z, 0, n1*sizeof(int));
		for(i = 0; i < nc; i++)
			for(cg = cgp[i]; cg; cg = cg->next)
				z[cg->varno]++;
		for(i = 0; i < nnc; i++)
			for(og = ogp[i]; og; og = og->next)
				z[og->varno]++;
		(*pf)(nl, "k%d\n", --n1);
		for(i = k = 0; i < n1; i++)
			(*pf)(nl, "%d\n", k += z[i]);
		free(z);
		}
	for(i = 0; i < nc; i++) {
		if (!(cg = cgp[i]))
			continue;
		k = 0;
		do k++;
			while((cg = cg->next));
		(*pf)(nl, "J%d %d\n", i, k);
		for(cg = cgp[i]; cg; cg = cg->next)
			(*pf)(nl, "%d %g\n", cg->varno, cg->coef);
		}
	Gput(pf, nl, 'J', nc, nnc, ogp);
	}

 static void
k1put(Pf *pf, FILE *nl, int *cs, real *a, int *rn, int nc, int n,
	int nnv, int nnc, ograd **ogp)
{
	int *cs1, i, j, j1, k, ftn, nz;
	cgrad *cg, *cg0, *cg1, **cgp, **cgq;
	ograd *og;

	ftn = cs[0];
	nz = cs[n] - ftn;
	k = n;
	if (nnc) {
		k += nnv;
		if (nz <= k)
			nz = k + 1;
		}
	cg0 = cg1 = (cgrad*)Malloc(nz*sizeof(cgrad) + nc*sizeof(cgrad*));
	cs1 = cs;
	if (nnc) {
		cs1 = (int*)cg1;
		for(i = 0; i < n; i++)
			cs1[i] = cs[i+1] - cs[i];
		while(i < k)
			cs1[i++] = 0;
		for(i = 0; i < nnc; i++)
			for(og = ogp[i]; og; og = og->next)
				cs1[og->varno]++;
		j = ftn;
		for(i = 0; i < k; i++) {
			j1 = j + cs1[i];
			cs1[i] = j;
			j = j1;
			}
		cs1[k] = j;
		}
	(*pf)(nl, "k%d\n", k - 1);
	for(i = 1; i < k; i++)
		(*pf)(nl, "%d\n", cs1[i] - ftn);
	memset(cgp = (cgrad**)(cg0 + nz), 0, nc*sizeof(cgrad*));
	k = cs[n] - ftn;
	for(i = n; --i >= 0;) {
		for(j = cs[i] - ftn; --k >= j;) {
			cg = cg1++;
			cg->coef = a[k];
			cgq = cgp + (cg->varno = rn[k] - ftn);
			cg->next = *cgq;
			*cgq = cg;
			}
		}
	k2put(pf, nl, cgp, nc, n, 0, nnv, nnc, ogp);
	free(cg0);
	}

 static int
LUcheck(int n, real *LU, real *u, int *nnep, int *nnrp)
{
	int i, nne, nnr;
	real L, U;

	nnr = nne = 0;
	if (!LU)
		return 1;
	for(i = 0; i < n; i++) {
		L = *LU++;
		U = u ? *u++ : *LU++;
		if (L < U) {
			if (L > negInfinity && U < Infinity)
				nnr++;
			}
		else if (U <= negInfinity
		 || L >= Infinity
		 || L > U
		 || L != L /* NaN */
		 || U != U)
			return 1;
		else
			nne++;
		}
	if (nnep) {
		*nnep = nne;
		*nnrp = nnr;
		}
	return 0;
	}

 static int
ogcheck(int n, int nn, ograd **ogp, int *nzp)
{
	int nz;
	ograd *og;

	if (!ogp)
		return 1;
	nz = 0;
	n += nn;
	while(nn--)
		for(og = *ogp++; og; og = og->next) {
			++nz;
			if (og->varno < 0
			 || og->varno >= n
			 || og->coef != og->coef
			 || og->coef == Infinity
			 || og->coef == negInfinity)
				return 1;
			}
	*nzp = nz;
	return 0;
	}

 static SufDesc*
reverse(SufDesc *sd)
{
	SufDesc *sn, *sp;
	sp = 0;
	while(sd) {
		sn = sd->next;
		sd->next = sp;
		sp = sd;
		sd = sn;
		}
	return sp;
	}

 int
fg_write_ASL(ASL *a, const char *stub, NewVCO *nu, int flags)
{
	ASL_fg *asl = (ASL_fg*)a;
	EvalWorkspace *ew;
	FILE *nl;
	Pf *pf;
	Staticfgw S;
	SufDesc *sd, *sd0;
	cexp *ce, *cee;
	char buf[256], *nbuf, *ts;
	const char *eol, *name, *obase, *s;
	expr **w;
	exprn *en, *en0;
	exprv *v;
	func_info *fi;
	int ak, c, i, j, *ip, *ipe, m, n, nnc, nne, nno, nnr, nnum, nnv, nnzc, nnzo;
	int nv, nx, oblen, rflag, rv;
	lincoef *Lc, *Lce;
	linpart *L;
	real *numv, *r, *re, t;
	static NewVCO nu0;

	ASL_CHECK(a, ASL_read_fg, "fg_write");
	if ((comc1 && !c_cexp1st) || (como1 && !o_cexp1st))
		return ASL_writeerr_badcexp1st;
	nnc = nne = nno = nnr = nnv = nnzc = nnzo = 0;
	if (!nu || (nu->nnv == 0 && nu->nnc == 0 && nu->nno == 0))
		nu = &nu0;
	else {
		nnc = nu->nnc;
		nno = nu->nno;
		nnv = nu->nnv;
		if ((nnv <= 0
		  || nnc < 0
		  || nno < 0
		  || nnc + nno <= 0
		  || nnc > 0) && !nu->LUnc)
			return ASL_writeerr_badNewVCO;
		if (LUcheck(nnv, nu->LUnv, nu->Unv, 0, 0))
			return ASL_writeerr_badNewVCO;
		n = n_var + nnv;
		if (nnc) {
			if (LUcheck(nnc, nu->LUnc, nu->Unc, &nnr, &nne))
				return ASL_writeerr_badNewVCO;
			if (ogcheck(n, nnc, nu->newc, &nnzc))
				return ASL_writeerr_badNewVCO;
			}
		if (nno) {
			if (ogcheck(n, nno, nu->newo, &nnzo))
				return ASL_writeerr_badNewVCO;
			if ((s = nu->ot))
			    for(i = 0; i < nno; i++)
				if (s[i] & ~1)
					return ASL_writeerr_badNewVCO;
			if ((r = nu->oc))
			    for(re = r + nno; r < re; r++) {
				if ((t = *r) <= negInfinity
				 || t >= Infinity
				 || t != t)
					return ASL_writeerr_badNewVCO;
				}
			}
		}

	s = name = obase = stub;
	while(*s)
	  switch(*s++) {
		case '/':
		case '\\':
			obase = s;
		}
	c = s - stub;
	nbuf = 0;
	oblen = s - obase;
	if (c <= 3 || strcmp(s - 3, ".nl")) {
		ts = buf;
		if (c + 4 > sizeof(buf))
			ts = nbuf = (char*)Malloc(c+4);
		memcpy(ts, stub, c);
		strcpy(ts+c, ".nl");
		name = ts;
		}
	else
		oblen -= 3;
	nl = fopen(name, "wb");
	if (nbuf)
		free(nbuf);
	if (!nl)
		return ASL_writeerr_openfail;
	memset(&S, 0, sizeof(S));
	S.asl = a;
	nv = n_var + nnv + ncom0 + ncom1;
	ew = asl->i.Ew0;
	S.al = ew->al;
	S.w = (expr**)ew->w;
	nnum = a->i.numlen / sizeof(real);
	en = en0 = (exprn*)new_mblk_ASL(a, htcl(nv*sizeof(exprv) + nnum*sizeof(exprn)));
	v = (exprv*)(en + nnum);
	numv = asl->i.numvals;
	w = S.w - nnum;
	while(w < S.w) {
		en->opno = -2;
		en->opcl = Enum;
		en->v = *numv++;
		*w++ = (expr*)en++;
		}
	for(i = 0; i < nv; ++i) {
		v->opno = -1;
		v->opcl = Evar;
		v->varno = i;
		*w++ = (expr*)v++;
		}
	rv = 0;
	i = setjmp(S.wjb);
	if (i) {
		rv = ASL_writeerr_badrops;
		goto ret;
		}
	if ((flags & ASL_write_ASCII) || (!(flags & ASL_write_binary) && !(binary_nl & 1))) {
		ak = 0;
		c = 'g';
		pf = aprintf;
		}
	else {
		ak = Arith_Kind_ASL;
		c = 'b';
		pf = bprintf;
		}
	S.nl_ = nl;
	S.pf_ = pf;
	eol = (char*)(flags & ASL_write_CR ? "\r\n" : "\n");
	fprintf(nl, "%c%d", c, n = ampl_options[0]);
	for(i = 1; i <= n; i++)
		fprintf(nl, " %d", ampl_options[i]);
	if (ampl_options[2] == 3)
		fprintf(nl, " %.g", ampl_vbtol);
	fprintf(nl, "\t# problem %.*s%s", oblen, obase, eol);
	fprintf(nl, " %d %d %d %d", n_var + nnv, n_con + nnc,
		n_obj + nno, nranges + nnr);
	s = "";
	if ((n = n_eqn + nne) >= 0) {
		fprintf(nl, " %d", n);
		s = ", eqns";
		}
	fprintf(nl, "\t# vars, constraints, objectives, ranges%s%s", s, eol);
	if (n_cc | nlcc)
		fprintf(nl, " %d %d %d %d%s%s", nlc, nlo, n_cc, nlcc,
		"\t# nonlinear constrs, objs; ccons: lin, nonlin", eol);
	else
		fprintf(nl, " %d %d\t# nonlinear constraints, objectives%s",
			nlc, nlo, eol);
	fprintf(nl, " %d %d\t# network constraints: nonlinear, linear%s",
		nlnc, lnc, eol);
	fprintf(nl, " %d %d %d%s%s", nlvc, nlvo, nlvb,
		"\t# nonlinear vars in constraints, objectives, both", eol);
	s = "";
	fprintf(nl, " %d %d", nwv, nfunc);
	if (ak | asl->i.flags) {
		fprintf(nl, " %d %d", ak, asl->i.flags);
		s = "; arith, flags";
		}
	fprintf(nl, "\t# linear network variables; functions%s%s", s, eol);
	fprintf(nl, " %d %d %d %d %d%s%s", nbv, niv, nlvbi, nlvci, nlvoi,
		"\t# discrete variables: binary, integer, nonlinear (b,c,o)",
		eol);
	fprintf(nl, " %d %d\t# nonzeros in Jacobian, gradients%s",
		nzc + nnzc, nzo + nnzo, eol);
	fprintf(nl, " 0 0\t# max name lengths: constraints, variables%s", eol);
	fprintf(nl, " %d %d %d %d %d\t# common exprs: b,c,o,c1,o1%s",
		comb, comc, como, comc1, como1, eol);

	for(i = 0; i < nfunc; i++) {
		fi = funcs[i];
		(*pf)(nl, "F%d %d %d %s\n", i, fi->ftype, fi->nargs, fi->name);
		}

	for(i = 0; i < 4; i++) {
		if (!(sd = asl->i.suffixes[i]))
			continue;
		nx = (&asl->i.n_var_)[i];
		for(sd = sd0 = reverse(sd); sd; sd = sd->next) {
			n = rflag = 0;
			if (sd->kind & ASL_Sufkind_real) {
				rflag = ASL_Sufkind_real;
				r = sd->u.r;
				re = r + nx;
				while(r < re)
					if (*r++)
						n++;
				}
			else {
				ip = sd->u.i;
				ipe = ip + nx;
				while(ip < ipe)
					if (*ip++)
						n++;
				}
			if (!n)
				continue;
			(*pf)(nl, "S%d %d %s\n", i | rflag, n, sd->sufname);
			j = 0;
			if (rflag) {
				r = sd->u.r;
				for(; j < nx; j++)
					if (r[j])
						(*pf)(nl, "%d %g\n", j, r[j]);
				}
			else {
				ip = sd->u.i;
				for(; j < nx; j++)
					if (ip[j])
						(*pf)(nl, "%d %d\n", j, ip[j]);
				}
			}
		reverse(sd0);
		}
	br(pf, nl, 'b', LUv, Uvx, n_var);
	br(pf, nl, 0, nu->LUnv, nu->Unv, nnv);
	if (!(flags & ASL_write_no_X0))
		iguess(pf, nl, 'x', X0, havex0, n_var, nnv, nu->x0);
	br(pf, nl, 'r', LUrhs, Urhsx, n_con);
	br(pf, nl, 0, nu->LUnc, nu->Unc, nnc);
	if (!(flags & ASL_write_no_pi0))
		iguess(pf, nl, 'd', pi0, havepi0, n_con, nnc, nu->d0);
	ce = cexps;
	n = n_var + nnv;
	for(cee = ce + comb + comc + como; ce < cee; ce++) {
		m = (L = ce->lp) ? (int)L->n : 0;
		(*pf)(nl, "V%d %d %d\n", n++, m, 0);
		if (L) {
			Lc = L->lc;
			for(Lce = Lc + m; Lc < Lce; ++Lc)
				(*pf)(nl, "%d %g\n", (int)Lc->varno, Lc->coef);
			}
		Ewput(&S, ce->o.e);
		}
	S.cexps1_ = asl->I.cexps1_;
	S.nv0 = n_var;
	S.com1off = S.nv0 + comb + comc + como;
	coput(&S, 'C', con_de, n_con, c_cexp1st, 0, 0, nnc, 0, 0);
	coput(&S, 'O', obj_de, n_obj, o_cexp1st, objtype, n_con,
		nno, nu->oc, nu->ot);
	if (A_vals)
		k1put(pf, nl, A_colstarts, A_vals, A_rownos, n_con, n_var,
			nnv, nnc, nu->newc);
	else
		k2put(pf, nl, Cgrad, n_con, n_var, 1, nnv, nnc, nu->newc);

	Gput(pf, nl, 'G', 0, n_obj, Ograd);
	Gput(pf, nl, 'G', n_obj, nno, nu->newo);

 ret:
	Del_mblk_ASL(a, en0);
	fclose(nl);
	return rv;
	}

 int
fg_wread_ASL(ASL *asl, FILE *f, int flags)
{
	want_xpi0 = 7;
	if (comc1)
		c_cexp1st = (int*)M1zapalloc((n_con + 1)*sizeof(int));
	if (como1)
		o_cexp1st = (int*)M1zapalloc((n_obj + 1)*sizeof(int));
	if (!(flags & ASL_keep_derivs)) {
		maxfwd = 0;
		want_derivs = 0;
		}
	if (!(flags & ASL_omit_all_suffixes))
		flags |= ASL_keep_all_suffixes;
	if (!(flags & ASL_forbid_missing_funcs))
		flags |= ASL_allow_missing_funcs;
	return fg_read_ASL(asl, f, flags);
	}

#ifdef __cplusplus
}
#endif
