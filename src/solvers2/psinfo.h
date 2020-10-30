/*******************************************************************
Copyright (C) 2017, 2020 AMPL Optimization, Inc.; written by David M. Gay.

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

#ifndef PSINFO_H0_included /*{*/

#ifndef EXPR_defined
#define EXPR_defined
 typedef struct expr expr;
#endif

 typedef struct
expr_nv {
	int op;
	int varno;
	} expr_nv;

 typedef union
pei {
	expr *e;
	expr **ep;
	int  *i;
	int  n;
	} pei;

#endif /*}*/
#ifdef PSHVREAD /*{{*/

#ifndef PSINFO_H2_included
#define PSINFO_H2_included
#ifndef NLP_H2_included
#undef PSINFO_H_included
#include "nlp2.h"
#endif
#endif /* PSINFO_H2_included */
#define cde cde2
#define la_ref la_ref2
#define range range2
#define rhead rhead2
#define psb_elem psb_elem2
#define psg_elem psg_elem2
#define ps_func ps_func2
#define dv_info dv_info2
#define split_ce split_ce2
#define ps_info ps_info2
#define psinfo psinfo2
#define Psbinfo Psbinfo2

 typedef struct
GOps {
	real O;
	real dL;
	real dL2;
	} GOps;

#else /* !PSHVREAD }{*/
#ifndef PSINFO_H1_included
#define PSINFO_H1_included
#ifndef NLP_H_included
#undef PSINFO_H_included
#include "nlp.h"
#endif
#endif
#endif /* PSHVREAD }}*/
#define PSINFO_H_included

 typedef struct la_ref la_ref;
 typedef struct range range;

 typedef struct
rhead {
	range *next, *prev;
	} rhead;

#ifndef PSINFO_H0_included

 typedef struct
Condptrs {
	derpblock *db;
	int *e;	/* for eval2() */
	int *f;	/* for hv_fwd  */
	int *b; /* for hv_back */
	int bder; /* >= 0 ==> hvg_back copies aO, adO here */
	} Condptrs;

 typedef struct
Minmaxptrs {
	derpblock *db;
	int *f;	/* for hv_fwd  */
	int *b; /* for hv_back */
	int  d; /* >= 0 ==> (Eresult*)&w[d] = this differentiable result */
	} Minmaxptrs;

#define MBLK_KMAX 30
#endif /* PSINFO_H0_included */

 typedef struct psb_elem psb_elem;

 struct
range {
	rhead	rlist;		/* list of all ranges */
	range	*hnext;		/* for hashing U */
	range	*hunext;	/* for hashing unit vectors */
	size_t	irange;		/* index of this range */
	size_t	uhlen;		/* sizeof(uHeswork) + (n-1)*sizeof(Ogptrs) */
	int	n;		/* rows in U */
	int	nv;		/* variables involved in U */
	int	nintv;		/* number of internal variables (non-unit */
				/* rows in U) */
	int	lasttermno;	/* termno of prev. use in this term */
				/* -1 ==> not yet used in this constr or obj. */
				/* Set to least variable (1st = 0) in this */
				/* range at the end of psedread. */
	int	lastgroupno;	/* groupno at last use of this term */
	unsigned int chksum;	/* for hashing */
	psb_elem *refs;		/* constraints and objectives with this range */
	int	*ui;		/* unit vectors defining this range */
				/* (for n >= nv) */
	linarg	**lap;		/* nonzeros in U: lap[i], 0 <= i < n */
	int	*cei;		/* common expressions: union over refs */
	int	hest;		/* nonzero ==> internal Hessian triangle */
				/* computed by hvpinit starts at ew->w + hest */
	};

#ifndef PSHVREAD
#define Ops Ops1
#endif

 struct
psb_elem {		/* basic element of partially-separable func */
	psb_elem *next;	/* for range.refs */
	range *U;
	int *ce;	/* common exprs if nonzero: ce[i], 1 <= i <= ce[0] */
	Ops o;		/* evaluation oplists: all, fwd, back */
	int conno;	/* constraint no. (if >= 0) or -2 - obj no. */
	int termno;
	int groupno;
	};
#undef Ops

 typedef struct
Psbinfo {
	derpblock **pdb;	/* pdb[0] = derps for the psb_elem's of this entity. */
				/* Subsequent entries are for cexps used, */
				/* followed by 0. */
	int *ce;		/* If not null, ce[i], 1 <= i <= ce[0], are */
				/* common exprs used (0 based), including split_ce */
				/* and linarg cells. */
	linarg **lap;		/* linargs used, followed by 0. */
	psb_elem *b, *be;	/* be - b = number of psb_elems. */
				/* When none, b == be == 0. */
	} Psbinfo;

 typedef struct
psg_elem {		/* group element details of partially-separable func */
	real	g0;	/* constant term */
	real	scale;
	Psbinfo pi;
	int	gm;	/* w[gm] = esum = result of summing g0, E and L */
			/* w[gm+1] = g1 = first derivative of g */
			/* w[gm+2] = g2 = 2nd derivative of g */
	int	nov;	/* length of ov array; see below */
	int	nu;	/* number of unary operators */
	int	*o;	/* for evaluating this element */
	int	*ov;	/* first deriv = g1 times (varno = ov[i], coef = w[gm+i+3]), */
			/* for 0 <= i < nov. */
	pei	g;	/* w[g->i[0]] = results of first unary operator (O, dO, dO2) */
			/* w[g->i[1]] = results of 2nd unary operator, etc. */
	pei	ge;	/* w[*ge.i] = results of last unary operator */
			/* There are ge->i - g->i  + 1 unary operators. */
	linpart *L;	/* the linear terms */
	} psg_elem;

 typedef struct
ps_func {
	Psbinfo pi;		/* the basic terms */
	psg_elem *g, *ge;	/* the group terms */
	} ps_func;

 typedef struct
dv_info {			/* defined variable info */
	ograd	*ll;		/* list of linear defined vars referenced */
	linarg	**nl;		/* nonlinear part, followed by 0 */
	real	scale;		/* scale factor for linear term */
	linarg	*lt;		/* linear term of nonlinear defined var */
	} dv_info;

 typedef struct
split_ce {
	range *r;
	int *ce;	/* common expressions (when r != 0) */
	linarg *la;	/* for new vars created by la_replace: r = 0, la != 0 */
	} split_ce;

#ifdef PSHVREAD

 struct
hes_fun {
	hes_fun *hfthread;
	cexp2	*c;
	ograd	*og;
	int	*vp;
	int	grdhes, n, nd;
	};

 typedef struct Hesoprod Hesoprod;
 struct
Hesoprod {
	Hesoprod *next;
	int *ov, *ove;	/* left vector */
	real *oc;	/* left vector */
	int *rv, *rve;	/* right vector */
	real *roc;	/* right vector */
	real coef;
	};

 typedef struct
Ogptrs {
	int *ov, *ove;
	real *oc;
	} Ogptrs;

 typedef struct uHeswork uHeswork;
 struct
uHeswork {
	uHeswork *next;
	range *r;
	int *ui, *uie;
	Ogptrs ogp[1];		/* scratch of length r->n */
	};

 typedef struct Ihinfo Ihinfo;
 struct
Ihinfo {
	Ihinfo *next;	/* for chaining ihinfo's with positive count */
	range *r;	/* list, on prev, of ranges with this ihd */
	int ihd;	/* internal Hessian dimension, min(n,nv) */
	int nr;		/* number of ranges with this ihd */
	};

#endif /* PSHVREAD */

 typedef struct
ps_info {
	Long merge;	/* for noadjust = 1 */
	ps_func	*cps;	/* constraints */
	ps_func	*ops;	/* objectives */
	dv_info	*dv;
	rhead rlist;
	linarg *lalist;	/* all linargs */
	int *dvsp0;	/* dvsp0[i] = subscript of first var into which */
			/* cexp i was split, 0 <= i < ncom */
	int *ndvsp;	/* cexp i was split into ndvsp[i] vars */
	size_t nran;	/* number of range structures (for sphes) */
	int nc1;	/* common expressions for just this function */
	int ns0;	/* initial number of elements */
	int ncom;	/* number of common expressions before splitting */
	int ndupdt;	/* duplicate linear terms in different terms */
	int ndupst;	/* duplicate linear terms in the same term */
	int nlttot;	/* total number of distinct linear terms */
	int ndvspcand;	/* # of defined variable candidates for splitting */
	int ndvsplit;	/* number of defined variables actually split */
	int ndvspin;	/* number of incoming terms from split defined vars */
	int ndvspout;	/* number of terms from split defined variables */
	int max_var1_;	/* used in psedread and pshvread */

#ifdef PSHVREAD
	/* Stuff for partially separable Hessian computations... */
	/* These arrays are allocated and zero-initialized by hes_setup, */
	/* which also supplies the cei field to ranges. */

	Ihinfo *ihi;
	Ihinfo *ihi1;	/* first with positive count */
	size_t zaplen;	/* for zeroing memory starting at &w[asl->P.rtodo] */
	int dOscratch, iOscratch, otodo, rtodo, utodo; /* subscripts in w for... */
	int nmax;	/* max{r in ranges} r->n */
	int ihdmax;	/* max possible ihd */
	int ihdmin;	/* min possible ihd > 0 and <= ihdmax, or 0 */
	int khesoprod;	/* used in new_Hesoprod in sputhes.c */
	int krnmax;	/* based on rnmax (below); set in hvpinit_nc_ASL() */
	int ndhmax;	/* Initial wh->ndhmax */
	int pshv_g1;	/* whether pshv_prod should multiply by g1 */
	int linmultr;	/* linear common terms used in more than one range */
	int linhesfun;	/* linear common terms in Hessian funnels */
	int nlmultr;	/* nonlin common terms used in more than one range */
	int nlhesfun;	/* nonlin common terms in Hessian funnels */
	int ncongroups;	/* # of groups in constraints */
	int nobjgroups;	/* # of groups in objectives */
	int rnmax;	/* max r->n for ranges r with r->n >= r->nv */
	int *zlsave;	/* for S->_zl */
	int *wkinit0, *wkinit2, *wkinitm1;	/* For initializing components of w */
						/* to 0, 2, -1, respectively. */
#endif /* PSHVREAD */
	split_ce *Split_ce;	/* for sphes_setup */
	} ps_info;

#ifdef PSHVREAD

 typedef struct
ASL_pfgh {
	Edagpars p;
	Edaginfo i;
	Edag2info I;
	ps_info2 P;
	} ASL_pfgh;

 typedef struct
Varval {
	real O;		/* 0: op value */
	real dO;	/* 1: deriv of op w.r.t. t in x + t*p */
	real aO;	/* 2: adjoint (in Hv computation) of op */
	real adO;	/* 3: adjoint (in Hv computation) of dO */
	} Varval;

 typedef struct
Eresult {
	real O;		/* 0: op value */
	real dO;	/* 1: deriv of op w.r.t. t in x + t*p */
	real aO;	/* 2: adjoint (in Hv computation) of O */
	real adO;	/* 3: adjoint (in Hv computation) of dO */
	real dL;	/* 4: deriv of op w.r.t. left operand L */
	real dL2;	/* 5: second partial w.r.t. L, L (or R,R for OPDIV) */
	real dR;	/* 6: deriv of op w.r.t. right operand R */
	real dLR;	/* 7: second partial w.r.t. L, R */
	real dR2;	/* 8: second partial w.r.t. R, R */
	} Eresult;

/* Unary ops and (and general binary ops with constant left operand) */
/* only use the first 2 fields. */
/* Abs and Less only use the first. */
/* Plus, Minus, Times do not use any. */
#else

 typedef struct
ASL_pfg {
	Edagpars p;
	Edaginfo i;
	Edag1info I;
	ps_info P;
	} ASL_pfg;

 typedef struct
Varval1 {
	real O;		/* 0: op value */
	} Varval1;

#endif /* PSHVREAD */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef PSINFO_H0_included
#define PSINFO_H0_included
typedef unsigned Long Ulong;

 struct
linarg {
	linarg *hnext;	/* for hashing */
	linarg *tnext;	/* next linear argument to this term */
	linarg *lnext;	/* for adjusting v->op */
	la_ref	*refs;	/* references */
	union {
		struct expr_vx *pv;
		int v;
		} u; /* variable (offset in w) that evaluates this linear term */
	int	*ov;	/* the nonzeros: varnos */
	real	*oc;	/* the nonzeros: coefs  */
	int	nnz;	/* number of nonzeros */
	int	termno;	/* helps tell whether new to this term */
	};

 typedef struct
tfinfo {
	char *wd;
	int *at;
	int *doff;
	int *fh;
	func_info *fi;
	int n;	/* total number of args */
	int nr;	/* number of numeric args */
	int nd; /* number differentiable numeric args */
	int ali; /* ew->al + ali = arglist for this invocation */
	} tfinfo;

#endif /* PSINFO_H0_included */
#ifdef PSHVREAD
 extern void duthes_ew_ASL(EvalWorkspace*, real *H, int nobj, real *ow, real *y);
 extern void duthese_ew_ASL(EvalWorkspace*, real *H, int nobj, real *ow, real *y, fint*);
 extern real eval2_ASL(int*, EvalWorkspace*);
 extern void fullhes_ew_ASL(EvalWorkspace*, real*H, fint LH, int nobj, real*ow, real *y);
 extern void fullhese_ew_ASL(EvalWorkspace*, real*H, fint LH, int nobj, real*ow, real *y, fint*);
 extern void hvpinit_ew_ASL(EvalWorkspace*, int hid_lim, int nobj, real *ow, real *y);
 extern void hvpinite_ew_ASL(EvalWorkspace*, int hid_lim, int nobj, real *ow, real *y, fint*);
 extern ASL_pfgh *pscheck_ASL(ASL*, const char*);
 extern void pshv_prod_ASL(EvalWorkspace*, range *r, int nobj, real *ow, real *y);
 extern fint sphes_setup_ew_ASL(EvalWorkspace*, SputInfo**, int nobj, int ow, int y, int ul);
 extern void sphes_ew_ASL(EvalWorkspace*, SputInfo**, real *H, int nobj, real *ow, real *y);
 extern void sphese_ew_ASL(EvalWorkspace*, SputInfo**, real *H, int nobj, real *ow, real *y, fint*);
 extern int xp_check_ASL(EvalWorkspace*, real*);
 extern void xpsg_check_ASL(EvalWorkspace*, int nobj, real *ow, real *y);
#else /* PSHVREAD */
 extern int xp1known_ew_ASL(EvalWorkspace*, real*, fint*);
#endif /* PSHVREAD */

#ifdef __cplusplus
	}

#define hvpinit(hx, no, ow, y) hvpinit_ew_ASL(asl->i.Ew0, hx, no, ow, y)

#endif
