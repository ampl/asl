/*******************************************************************
Copyright (C) 2016 AMPL Optimization, Inc.; written by David M. Gay.

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

#ifndef NLP_H_included
#define NLP_H_included

#ifndef ASL_included
#include "asl.h"
#endif

typedef struct cde cde;
typedef struct cexp cexp;
typedef struct cexp1 cexp1;

#define obj1val   obj1val_ew_ASL
#define obj1grd   obj1grd_ew_ASL
#define con1val   con1val_ew_ASL
#define jac1val   jac1val_ew_ASL
#define con1ival  con1ival_ew_ASL
#define con1grd   con1grd_ew_ASL
#define lcon1val  lcon1val_ew_ASL
#define x1known   x1known_ew_ASL
#ifndef ETYPE
#define ETYPE int
#endif

 typedef struct
Ops1 { int *e; } Ops1;

 struct
cde {
	Ops1	  o;
	derpblock *db;
	int	  *dvref;
	linpart   **c1lp; /* subset of this cde's cexp1's that have linear parts */
	uint	  af1;	  /* first defined var used only by this entity */
	uint	  afn;	  /* number of such defined vars */
	};

 struct
cexp1 {
	Ops1	 o;
	linpart *lp;
	};

 struct
cexp {
	Ops1	o;	/* nonlinear part */
	derpblock *db;	/* basic derivative propagation (once per x) */
	derpblock *dbf;	/* funneled derivative prop; == db if no funnel */
	linpart *lp;	/* linear part */
	int	*vref;	/* vref != 0 ==> vref[i], 3 <= i < vr[0]+3 = vars referenced. */
			/* For vr[1]+3 <= i < vr[0]+3, the vars are defined vars. */
			/* vr[2] = funnelkind (0, 1 or 2). */
	};

 typedef struct
Edag1info {
	cde	*con_de_;	/* constraint deriv. and expr. info */
	cde	*lcon_de_;	/* logical constraints */
	cde	*obj_de_;	/* objective  deriv. and expr. info */

			/* stuff for "defined" variables */
	cexp	*cexps_;
	cexp1	*cexps1_;
	cexp	**dvfb;		/* funnels shared by constraints and objectives */
	cexp	**dvfc;		/* funnels shared by constraints */
	cexp	**dvfo;		/* funnels shared by objectives */
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
	real *w0;		/* w0 is copied to the start of EvalWorkspace.w */
	EvalWorkspace	*ew0;	/* used by evaluation routines called with NULL */
				/* for their EvalWorkspace argument */
	} Edag1info;

 typedef struct
ASL_fg {
	Edagpars  p;
	Edaginfo  i;
	Edag1info I;
	} ASL_fg;

 typedef struct
Invd1 {
	func_info *fi;
	char *dig;
	int *at;
	int n, nr;
	} Invd1;

#ifdef __cplusplus
 extern "C" {
#endif
 extern void con1grd(EvalWorkspace*, int nc, real *X, real *G, fint *nerror);
 extern real con1ival(EvalWorkspace*,int nc, real *X, fint *ne);
 extern void con1val(EvalWorkspace*, real *X, real *F, fint *nerror);
 extern void jac1val(EvalWorkspace*, real *X, real *JAC, fint *nerror);
 extern int  lcon1val(EvalWorkspace*, int nc, real *X, fint *ne);
 extern void obj1grd(EvalWorkspace*, int nobj, real *X, real *G, fint *nerror);
 extern real obj1val(EvalWorkspace*, int nobj, real *X, fint *nerror);
 extern int x1known(EvalWorkspace*, real*, fint*);
#ifdef __cplusplus
	}
#endif

#define cexps	asl->I.cexps_
#define cexps1	asl->I.cexps1_
#define con_de	asl->I.con_de_
#define f_b	asl->I.f_b_
#define f_c	asl->I.f_c_
#define f_o	asl->I.f_o_
#define lcon_de	asl->I.lcon_de_
#define obj_de	asl->I.obj_de_

#endif /* NLP_H_included */
