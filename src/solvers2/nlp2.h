/*******************************************************************
Copyright (C) 2017 AMPL Optimization, Inc.; written by David M. Gay.

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

/* Variant of nlp.h for Hessian times vector computations. */

#ifndef NLP_H2_included
#define NLP_H2_included

#ifndef ASL_included
#include "asl.h"
#endif

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct cde2 cde2;
typedef struct cexp2 cexp2;
#ifndef EXPR_defined
#define EXPR_defined
typedef struct expr expr;
#endif
typedef struct hes_fun hes_fun;

 typedef struct
Ops {
	int *e;	/* all evaluations */
	int *f;	/* forward differentiable ops */
	int *b; /* reverse differentiable ops */
	} Ops;

 struct
cde2 {	/* same as cde (so far) */
	Ops o;
	linpart **c1lp; /* subset of this cde's cexp's that have linear parts */
	uint af1;	/* first defined var used only by this entity */
	uint afn;	/* number of such defined vars */
	};

 struct
cexp2 {
	Ops o;		/* evaluation oplist: all, fwd, back */
	derpblock *db;	/* basic derivative propagation (once per x) */
	derpblock *dbf;	/* funneled derivative prop */
	linpart   *lp;
	linarg	  *la;	/* corresponding to lp */
	int *vref;	/* vref != 0 ==> vref[i], 3 <= i < vr[0]+3 = vars referenced. */
			/* For vr[1]+3 <= i < vr[0]+3, the vars are defined vars. */
			/* vr[2] = funnelkind (0, 1 or 2). */
	hes_fun	*hfun;
	int varno;
	int splitvno;	/* first varno added to split this var (0 if none) */
	};

 typedef struct
Edag2info {
	cde2	*con2_de_;	/* constraint deriv. and expr. info */
	cde2	*lcon2_de_;	/* logical constraints */
	cde2	*obj2_de_;	/* objective  deriv. and expr. info */

			/* stuff for "defined" variables */
	cexp2	*cexps2_;
	cexp2	**dvfb;		/* funnels shared by constraints and objectives */
	cexp2	**dvfc;		/* funnels shared by constraints */
	cexp2	**dvfo;		/* funnels shared by objectives */
	hes_fun	*hesthread;
	char	*c_class;	/* class of each constraint: */
				/* 0 = constant */
				/* 1 = linear */
				/* 2 = quadratic */
				/* 3 = general nonlinear */
	char	*o_class;	/* class of each objective */
	char	*v_class;	/* class of each defined variable */
	int	c_class_max;	/* max of c_class values */
	int	o_class_max;	/* max of o_class values */
				/* The above are only computed if requested */
				/* by the ASL_find_c_class and */
				/* ASL_find_o_class bits of the flags arg */
				/* to pfgh_read() and pfg_read() */

	/* for parallel evaluations */

	size_t w_len;		/* length in bytes of EvalWorkspace.w */
	size_t w0_len;		/* w0 (below) is w0_len bytes long */
	size_t uhlen;
	real *w0;		/* w0 is copied to the start of EvalWorkspace.w */
	EvalWorkspace	*ew0;	/* used by evaluation routines called with NULL */
				/* for their EvalWorkspace argument */
	int gscrx;
	int ntseen;
	int nhop;
	} Edag2info;

 typedef struct
ASL_fgh {
	Edagpars  p;
	Edaginfo  i;
	Edag2info I;
	} ASL_fgh;

 extern void com21eval_ASL(EvalWorkspace*, int, int);
 extern void com2eval_ASL(EvalWorkspace*, int, int);
#ifdef __cplusplus
	}
#endif

#ifndef SKIP_NL2_DEFINES

#define cexps		asl->I.cexps2_
#define con_de		asl->I.con2_de_
#define f_b		asl->I.f2_b_
#define f_c		asl->I.f2_c_
#define f_o		asl->I.f2_o_
#define lcon_de		asl->I.lcon2_de_
#define obj_de		asl->I.obj2_de_

#define cde	cde2
#define cexp	cexp2
#define ei	ei2

#define com1eval	com21eval_ew_ASL
#define comeval		com2eval_ew_ASL
#undef  r_ops
#define r_ops		r2_ops_ASL

#ifndef PSHVREAD
#define f_OPIFSYM	f2_IFSYM_ASL
#define f_OPPLTERM	f2_PLTERM_ASL
#define f_OPFUNCALL	f2_FUNCALL_ASL
#define f_OP1POW	f2_1POW_ASL
#define f_OP2POW	f2_2POW_ASL
#define f_OPCPOW	f2_CPOW_ASL
#define f_OPPLUS	f2_PLUS_ASL
#define f_OPSUMLIST	f2_SUMLIST_ASL
#define f_OPHOL		f2_HOL_ASL
#define f_OPPOW		f2_POW_ASL
#define f_OPVARVAL	f2_VARVAL_ASL
#endif

/* treat if as vararg, minusL as plusL, binaryL as unary */

#endif /* SKIP_NL2_DEFINES */

#endif /* NLP_H2_included */
