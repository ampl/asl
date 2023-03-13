/* An interface from AMPL to Xpress-MP */
/* 1997-2008 Y. Colombani and others, Dash Optimization */
/* Requires Xpress-Optimizer libraries 18.00 or later */

/* Adjustments and additions to keywords and their descriptions, */
/* modifications for QPs, etc., by David M. Gay (2005). */

/* Remark: The options list MUST be in alphabetical order */

/* Some of the material in this file is derived from sample
   AMPL/solver interfaces that are available from netlib and
   which bear copyright notices of the following form... */

/****************************************************************
Copyright (C) 1997-2001 Lucent Technologies
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/

#include "xprs.h"
#if XPVERSION >= 30
#include <sys/types.h>
#include <sys/stat.h>
#endif
#include "nlp.h"
#include "getstub.h"
#include <signal.h>
#include <stdarg.h>

#undef  MININT
#define MININT (-XPRS_MAXINT - 1)

#ifndef XPRESS
#define XPRESS NULL
#endif

#if defined(X64_bit_pointers) && !defined(USE32BITOFFSETS)
#define ForceZ | ASL_use_Z
	typedef size_t CStype;
#undef nzc
#undef A_colstarts
#define nzc asl->i.nZc_
#define A_colstarts asl->i.A_colstartsZ_
#define XPRSloadglobal(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23) XPRSloadglobal64(a1,a2,a3,a4,a5,a6,a7,a8,(XPRSint64*)a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,(XPRSint64*)a21,a22,a23)
#define XPRSloadlp(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14) XPRSloadlp64(a1,a2,a3,a4,a5,a6,a7,a8,(XPRSint64*)a9,a10,a11,a12,a13,a14)
#define XPRSloadqcqp(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24) XPRSloadqcqp64(a1,a2,a3,a4,a5,a6,a7,a8,(XPRSint64*)a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,(XPRSint64*)a21,a22,a23,a24)
#define XPRSloadqcqpglobal(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33) XPRSloadqcqpglobal64(a1,a2,a3,a4,a5,a6,a7,a8,(XPRSint64*)a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,(XPRSint64*)a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,(XPRSint64*)a31,a32,a33)
#define XPRSloadqglobal(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27) XPRSloadqglobal64(a1,a2,a3,a4,a5,a6,a7,a8,(XPRSint64*)a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,(XPRSint64*)a25,a26,a27)
#define XPRSloadqp(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18) XPRSloadqp64(a1,a2,a3,a4,a5,a6,a7,a8,(XPRSint64*)a9,a10,a11,a12,a13,a14,a15,a16,a17,a18)
#else
	typedef int CStype;
#define ForceZ /*nothing*/
#endif


/* Basic status values returned by XPRESS function getbasis() */
#define XP_NBASLO 0  /* Vector non-basic at lower bound */
#define XP_BASIC  1  /* Vector basic */
#define XP_NBASUP 2  /* Vector non-basic at upper bound */


/* Define problem pointer */

XPRSprob prob;

typedef struct
 dims {
	double  *x;
	double  *y;
	int  *cstat;
	int  *rstat;
	SufDesc *csd;
	SufDesc *rsd;
	size_t nelq;
	int miqp;
  } dims;

 static int
	Ray = 0,
	Round = 1,
	advance = 1,
	lazy = 1,
	mipstart = 1,
	mipststat = 1,
	objrep = 2,
	sos = 1,
	sos2 = 1;

static char *set_known(Option_Info *oi, keyword *kw, char *v);
static char *set_int(Option_Info *oi, keyword *kw, char *v);
static char *set_dbl(Option_Info *oi, keyword *kw, char *v);
static char *set_fln(Option_Info *oi, keyword *kw, char *v);
static void nonlin(int n, char *what);
static void amplin(char *stub, char *argv[], dims*);
static void amplout(dims*);
static void show_times(void);
static void xperror(const char *where, ...);
static void killtempprob(void);
static void mip_priorities(void);

static ASL *asl;
struct LU_bounds {real lower, upper;};
static double objadj;

enum glstat_e {
      GLSTAT_NOPROB           = 0,
      GLSTAT_LP_UNFINISHED    = 1,
      GLSTAT_LP_FINISHED      = 2,
      GLSTAT_UNFINISHED_NOSOL = 3,
      GLSTAT_UNFINISHED_SOL   = 4,
      GLSTAT_FINISHED_NOSOL   = 5,
      GLSTAT_FINISHED_SOL     = 6,
      GLSTAT_MIP_UNBOUNDED    = 7
     };

enum lpstat_e {
      LPSTAT_OPTIMAL      = 1,
      LPSTAT_INFEASIBLE   = 2,
      LPSTAT_CUTOFF       = 3,
      LPSTAT_UNFINISHED   = 4,
      LPSTAT_UNBOUNDED    = 5,
      LPSTAT_CUTOFFINDUAL = 6,
      LPSTAT_UNSOLVED	  = 7,
      LPSTAT_NONCONVEX    = 8
     };

enum known_parameters {
	set_primal,        /* Choose algorithm */
#ifdef RWA_DEBUG
	set_debug,
#endif
	set_dual,
	set_barrier,
	set_maxim,         /* What to do */
	set_minim,
	set_relax,
	set_timing,
	set_iis,
	set_network,
	set_bestbound
     };

 typedef struct
PBuf { char *s, *se; } PBuf;

 static void
Bpf(PBuf *b, const char *fmt, ...)
{
	va_list ap;

	if (b->s < b->se) {
		va_start(ap, fmt);
		b->s += Vsnprintf(b->s, b->se - b->s, fmt, ap);
		va_end(ap);
		}
	}

 typedef struct Defer_setting Defer_setting;
 struct
Defer_setting {
	Defer_setting *next;
	int (*Setter)(Defer_setting*);
	keyword *kw;
	int ipar;
	union {
		int i;
		double d;
		} u;
	};

 static Defer_setting *Defer1, **Defer_next = &Defer1;

#if XPVERSION >= 39
#undef XPRSsetcbmessage
#define XPRSsetcbmessage(a,b,c) XPRSaddcbmessage(a,b,c,0)
#ifndef USE_FICO_MINMAXIM
 static int XPRS_CC
myXPRSmaxim(XPRSprob prob, const char *flags)
{
	if (XPRSchgobjsense(prob, XPRS_OBJ_MAXIMIZE))
		fprintf(Stderr, "XPRSchgobjsense(prob, XPRS_OBJ_MAXIMIZE) failed\n");
	return (flags[1] == 'g') ? XPRSmipoptimize(prob, flags) : XPRSlpoptimize(prob, flags);
	}
 static int XPRS_CC
myXPRSminim(XPRSprob prob, const char *flags)
{
	if (XPRSchgobjsense(prob, XPRS_OBJ_MINIMIZE))
		fprintf(Stderr, "XPRSchgobjsense(prob, XPRS_OBJ_MINIMIZE) failed\n");
	return (flags[1] == 'g') ? XPRSmipoptimize(prob, flags) : XPRSlpoptimize(prob, flags);
	}
#define XPRSmaxim myXPRSmaxim
#define XPRSminim myXPRSminim
#endif
#endif

 static Defer_setting *
new_Defer_setting(int (*Dset)(Defer_setting*), keyword *kw, int ipar)
{
	Defer_setting *ds = (Defer_setting*)M1alloc(sizeof(Defer_setting));
	*Defer_next = ds;
	Defer_next = &ds->next;
	ds->Setter = Dset;
	ds->kw = kw;
	ds->ipar = -ipar;
	return ds;
	}

 void
Do_Defer(VOID)
{
	Defer_setting *ds;
	int nbad;

	*Defer_next = 0;
	nbad = 0;
	for(ds = Defer1; ds; ds = ds->next)
		nbad += (*ds->Setter)(ds);
	if (nbad)
		exit(2);
	}

static char probname[L_tmpnam+4]; /* Use a temporary problem name */
static char *endbasis, *logfile, *startbasis, *writeprob;
static const char *wpflags;
static double Times[4];      /* Timing stats */
static int bestbound;	     /* return .bestbound suffix */
int timing=0;
int iis_find=0;              /* find IIS */
#ifdef RWA_DEBUG
char debugopt[]="                   ";
#endif

static int nobj=1;           /* Which objective we optimise */
                             /* Which optimization function */
static int (XPRS_CC *Optimise)(XPRSprob ,const char *);

static int prtmsg = 3;       /* Message output level */

char optimopt[4]={'X', 0,0,0};  /* Which algorithm */
/* [0] is b|d|X(meaning primal will be used - warning: non standard) */
/* [1] is l|g|blank */

#ifdef XPRS_MSP_SOLPRB_OBJ /*{*/

#include "xprs_mse_defaulthandler.h"

typedef int (XPRS_CC *mse_handler)(XPRSmipsolenum, XPRSprob prob, XPRSmipsolpool, void*,
				   int*, const double*, const int, const double, double*, int*, int*);
static int (XPRS_CC *MSEopt)(XPRSmipsolenum, XPRSprob, XPRSmipsolpool, mse_handler, void*, int*);
static XPRSmipsolenum mse;
static XPRSmipsolpool msp;
static int dualred, dupcol, nbest, npool;
static char *poolstub;

 extern int XPRS_CC XPRSgetcalist(const int **field_header_ids, const char ***field_names, const int **field_type_flags);

 static int
icompar(const void *a, const void *b, void *v)
{ return ((const int*)v)[*(int*)a] - ((const int*)v)[*(int*)b]; }

 static int
ncompar(const void *a, const void *b, void *v)
{ return strcmp(((const char**)v)[*(int*)a], ((const char**)v)[*(int*)b]); }

 static char *
set_par(Option_Info *oi, keyword *kw, char *v0)
{
	int c, d, f, i, id, j, k, l, m, nk, *p, q, slen, ti;
	char *b, *be, buf[1024], pbuf[64], *pe, *pn, *rv, *v;
	const char **kn, *what;
	double td;
	static const char **knames;
	static const int *kid, *ktype;
	static int *Pi, *Pn, n;

	enum {Ivar = 1, Dvar = 2, Svar = 4, Cvar = 0x20};

	if (!knames) {
		if ((l = XPRSgetcalist(&kid, &knames, &ktype))) {
			what = "XPRSgetcalist";
 getfail:
			printf("Bug: surprise %s failure (return %d).\n", what, l);
			exit(1);
			}
		for(kn = knames; *kn; ++kn);
		n = (int)(kn - knames);
		if (n <= 0) {
			printf("XPRSgetcalist provided no names.\n");
			exit(1);
			}
		Pi = (int*)M1alloc(n*(2*sizeof(int)));
		Pn = Pi + n;
		for(i = j = 0; i < n; ++i) {
			if (ktype[i] & Cvar) {
				Pi[j] = i;
				Pn[j++] = i;
				}
			}
		if (j > 1) {
			qsortv(Pi, j, sizeof(int), icompar, (void*)kid);
			qsortv(Pn, j, sizeof(int), ncompar, (void*)knames);
			}
		n = j;
		}

	c = *(unsigned char*)(rv = v0);
	what = kw->name;
	nk = 0;
	if (c >= '0' && c <= '9') {
		f = (int)strtol(v0, &rv, 10);
		if ((c = *rv) != '=') {
 noeq:
			m = rv - v0;
			printf("Expected '=' and a parameter value after "
				"%s%s%.*s, not \"%s\".\n", what, oi->eqsign, m, v0, rv);
 bad:
			badopt_ASL(oi);
			goto ret;
			}
		nk = 1;
		/* binary search */
		p = Pi;
		m = n;
		while(m > 0) {
			k = m >> 1;
			d = kid[p[k]] - f;
			if (!d) {
				j = p[k];
				goto have_f;
				}
			if (d > 0)
				m = k;
			else {
				p += ++k;
				m -= k;
				}
			}
		printf("control param %d not found.\n", f);
		goto bad;
		}
	if (c == '?') {
		m = 0;
		for(i = n; i > 0; ) {
			j = strlen(knames[Pn[--i]]);
			if (m < j)
				m = j;
			}
		printf("\n%d control parameters:\n%-*s  number  value\n", n, m, "name");
		for(i = 0; i < n; ++i) {
			j = Pn[i];
			k = ktype[j];
			id = kid[j];
			printf("%-*s %6d   ", m, knames[j], id);
			if (k & Ivar) {
				ti = 0;
				if ((l = XPRSgetintcontrol(prob, id, &ti))) {
					what = "XPRSgetintcontrol";
					goto getfail;
					}
				printf("%d\n", ti);
				}
			else if (k & Dvar) {
				td = 0.;
				if ((l = XPRSgetdblcontrol(prob, id, &td))) {
					what = "XPRSgetdblcontrol";
					goto getfail;
					}
				printf("%#.g\n", td);
				}
			else if (k & Svar) {
				l = XPRSgetstringcontrol(prob, id, buf, (int)sizeof(buf), &slen);
				if (l) {
					what = "XPRSgetstringcontrol";
					goto getfail;
					}
				printf("\"%s\"\n", buf);
				}
			else
				printf("unknown ktype #%x\n", k);
			}
		oi->option_echo &= ~ASL_OI_echothis;
		goto ret;
		}
	pn = pbuf;
	pe = pn + sizeof(pbuf) - 1;
	for(;;) {
		if (c <= ' ' || c == '=' || c == ',')
			break;
		if (pn >= pe) {
			for(v = rv; (c = *(unsigned char*)++v) > ' ' && c != '=' && c != ',';);
			printf("Oversize param name \"%.*s\".\n", (int)(v-v0), v0);
			goto bad;
			}
		if (c >= 'a' && c <= 'z')
			c += 'A' - 'a';
		*pn++ = c;
		c = *(unsigned char*)++rv;
		}
	if (c != '=')
		goto noeq;
	if (pn <= pbuf) {
		printf("Expected a param name after \"%s%s\".\n", what, oi->eqsign);
		goto bad;
		}
	*pn = 0;
	/* binary search */
	p = Pn;
	m = n;
	while(m > 0) {
		k = m >> 1;
		d = strcmp(knames[p[k]], pbuf);
		if (!d) {
			f = kid[j = p[k]];
			goto have_f;
			}
		if (d > 0)
			m = k;
		else {
			p += ++k;
			m -= k;
			}
		}
	printf("control param \"%s\" not found.\n", pbuf);
	goto bad;
 have_f:
	while(*(unsigned char*)++rv <= ' ') {
		if (!*rv) {
			printf("Expected a value after \"param%s%.*s\".\n",
				oi->eqsign, (int)(rv-v0), v0);
			goto bad;
			}
		}
	k = ktype[j];
	if ((c = *rv) == '?') {
		oi->option_echo &= ~ASL_OI_echothis;
		if (nk)
			printf("param=%d=", f);
		else
			printf("param=%s=", knames[j]);
		if (k & Ivar) {
			ti = 0;
			if ((l = XPRSgetintcontrol(prob, f, &ti))) {
				what = "XPRSgetintcontrol";
				goto getfail;
				}
			printf("%d\n", ti);
			}
		else if (k & Dvar) {
			td = 0.;
			if ((l = XPRSgetdblcontrol(prob, f, &td))) {
				what = "XPRSgetdblcontrol";
				goto getfail;
				}
			printf("%#.g\n", td);
			}
		else if (k & Svar) {
			l = XPRSgetstringcontrol(prob, f, buf, (int)sizeof(buf), &slen);
			if (l) {
				what = "XPRSgetstringcontrol";
				goto getfail;
				}
			printf("\"%s\"\n", buf);
			}
		else
			printf("unknown ktype #%x\n", k);
		goto ret;
		}
	m = rv - v0;
	if (k & (Ivar|Dvar)) {
		v = rv;
		td = strtod(v, &rv);
		if (rv <= v || ((c = *(unsigned char*)rv) > ' ' && c != ',')) {
			printf("Expected a numeric value after \"param=%.*s\", not \"%s\".\n",
				m, v0, v);
			goto bad;
			}
		if (k & Ivar) {
			ti = td;
			if (ti != td)
				printf("Assuming \"param%s%.*s%d\".\n", oi->eqsign, m, v0, ti);
			if ((l = XPRSsetintcontrol(prob, f, ti))) {
				what = "XPRSsetintcontrol";
				goto getfail;
				}
			}
		else if ((l = XPRSsetdblcontrol(prob, f, td))) {
			what = "XPRSsetdblcontrol";
			goto getfail;
			}
		oi->option_echo &= ~ASL_OI_echothis;
		goto ret;
		}
	b = buf;
	be = buf + sizeof(buf);
	if (c == '"' || c == '\'') {
		for(;;) {
			if ((q = *++rv) == c && *++rv != c)
				break;
			if (!q) {
				i = c == '"' ? '\'' : '"';
				printf("Missing closing %c%c%c in value for \"param%s%.*s\".\n",
					i, c, i, oi->eqsign, m, v0);
				badopt_ASL(oi);
				return rv;
				}
			if (b < be)
				*b++ = q;
			}
		}
	else while(c && c != ',') {
		if (b < be)
			*b++ = c;
		c = *++rv;
		}
	if (b >= be)  {
		printf("Oversize value for \"%s=", what);
		if (nk)
			printf("%d=\".\n", f);
		else
			printf("%s=\".\n", knames[j]);
		oi->option_echo &= ~ASL_OI_echothis;
		badopt_ASL(oi);
		}
	else {
		*b = 0;
		if ((l = XPRSsetstrcontrol(prob, f, buf))) {
			what = "XPRSsetstrcontrol";
			goto getfail;
			}
		}
 ret:
	while(*(unsigned char*)rv > ' ' && *rv++ != ',');
	return rv;
	}

 static int XPRS_CC
MSEsolve(XPRSprob prb, const char *op)
{ return MSEopt(mse, prb, msp, XPRS_mse_defaulthandler, 0, &nbest); }

#endif /*}*/

#if XPVERSION >= 30
 static char *tunerdir, *tunermethodfile;
 static int archconsistent;
#endif

static char
	advance_desc[]		= "whether to use an initial basis, if available:\n\
			0 = no, overriding mipstartstatus;\n\
			1 = yes (default), subject to mipstartstatus.\n\
		In an AMPL session, \"option send_statuses 0;\" is preferable\n\
		to \"option xpress_options '... advance=0 ...';\".",
#ifdef XPRS_ALGAFTERCROSSOVER
	algaftercrossover_desc[] = "algorithm for final cleanup after running the barrier\n\
		algorithm:\n\
			1 = automatic choice (default)\n\
			2 = dual simplex\n\
			3 = primal simplex\n\
			4 = concurrent",
#endif
#ifdef XPRS_ALGAFTERNETWORK
	algafternetwork_desc[] = "algorithm for final cleanup after the network simplex\n\
		algorithm:\n\
			1 = automatic choice (default)\n\
			2 = dual simplex\n\
			3 = primal simplex",
#endif
#if XPVERSION >= 30
	archconsistent_desc[]	= "whether to force the same execution path to be\n\
		independent of the platform architecture:\n\
			0 = no (default)\n\
			1 = yes",
#endif
#ifdef XPRS_AUTOCUTTING
	autocutting_desc[] =
		"whether to automatically decide if to generate cutting\n\
		planes at local nodes (overriden by cutfreq):\n\
			-1 = automatic (default)\n\
			 0 = disabled\n\
			 1 = enabled",
#endif
	autoperturb_desc[]	= "whether to introduce perturbations when the simplex method\n\
		encounters too many degenerate pivots:\n\
			1 = yes (default); 0 = no",
#ifdef XPRS_AUTOSCALING
	autoscaling_desc[]  = "whether the Optimizer should automatically select between\n\
		different scaling algorithms:\n\
		    -1 = automatic (default)\n\
		     0 = disabled\n\
		     1 = cautious strategy.  Non-standard scaling will only\n\
			    be selected if it appears to be clearly superior\n\
		     2 = moderate strategy\n\
		     3 = aggressive strategy.  Standard scaling will only be\n\
			    selected if it appears to be clearly superior",
#endif
	backtrack_desc[]	= "choice of next node when solving MIP problems:\n\
			-1 = automatic choice (default)\n\
			 1 = withdrawn; formerly choice 2 until a feasible\n\
				integer solution has been found, then\n\
				Forrest-Hirst-Tomlin choice\n\
			 2 = node with best estimated solution\n\
			 3 = node with best bound on the solution (default)\n\
			 4 = deepest node (depth-first search)\n\
			 5 = highest node (breadth-first search)\n\
			 6 = earliest-created node\n\
			 7 = most recently created node\n\
			 8 = random choice\n\
			 9 = node with fewest LP relaxation infeasibilities\n\
			10 = combination of 2 and 9\n\
			11 = combination of 2 and 4",
	backtracktie_desc[]	= "how to break ties for the next MIP node:  same choices as\n\
		for \"backtrack\"",
#ifdef XPRS_BARALG /*{*/
	baralg_desc[]		= "which barrier algorithm to use with \"barrier\":\n\
			-1 = automatic choice (default with just \"barrier\")\n\
			 1 = infeasible-start barrier algorithm\n\
			 2 = homogeneous self-dual barrier algorithm\n\
			 3 = start with 2 and maybe switch to 1 while solving",
#ifdef XPRS_BARCORES
	barcores_desc[]		= "if positive, number of CPU cores to assume present when\n\
		using the barrier algorithm.  Default = -1, which means\n\
		automatic choice.",
#endif
#endif /*}*/
	barcrash_desc[]		= "choice of crash procedure for crossover:\n\
			0 = no crash\n\
			1-6 = available strategies (default 4):\n\
			1 = most conservative, 6 = most aggressive",
	bardualstop_desc[]	= "barrier method convergence tolerance on\n\
		dual infeasibilities; default = 0 (automatic choice)",
	Bstop_desc[]	= "barrier method convergence tolerance on the relative\n\
		duality gap; default = 0",
#ifdef XPRS_BARGAPTARGET
	bargaptarget_desc[] =
		"barrier algorithm target tolerance for the relative duality\n\
		gap.  If not satisfied and no further progress is possible\n\
		but barstopgap is satisfied, then the current solution is\n\
		considered optimal.",
#endif
	barindeflimit_desc[]	=
		"maximum indefinite factorizations to allow in the barrier\n\
		algorithm for solving a QP: stop when the limit is hit;\n\
		default = 15",
	bariterlimit_desc[]	= "maximum number of Newton Barrier iterations; default = 500",
#ifdef XPRS_BARKERNEL
	barkernel_desc[]	= "how the barrier algorithm weights centrality:\n\
			>= +1.0 ==> more emphasis on centrality\n\
			<= -1.0 ==> each iteration, adaptively select a value\n\
				from [+1, -barkernel].\n\
		Default = 1.",
#endif
#ifdef XPRS_BAROBJPERTURB
	barobjperturb_desc[] =
		"defines how the barrier perturbs the objective (default\n\
		1e-6); values >0 let the optimizer decide if to perturb the\n\
		objective, values <0 force the perturbation:\n\
			n > 0 = automatic decison, scale n\n\
			    0 = turn off perturbation\n\
			n < 0 = force perturbation by abs(n)",
#endif
	barobjscale_desc[]	= "how the barrier algorithm scales the objective:\n\
			 -1 = automatic chocie (default)\n\
			  0 = scale by the geometric mean of the objective\n\
			      coefficients\n\
			> 0 = scale so the argest objective coefficient in\n\
			      absolute value is <= barobjscale.\n\
		When the objective is quadratic, the quadratic diagonal\n\
		is used in determining the scale.",
	barorder_desc[]		= "Cholesky factorization pivot order for barrier algorithm:\n\
			0 = automatic choice (default)\n\
			1 = minimum degree\n\
			2 = minimum local fill\n\
			3 = nested dissection",
#ifdef XPRS_BARORDERTHREADS
	barorderthreads_desc[]	= "number of threads to use when choosing a pivot order for\n\
		Cholesky factorization; default 0 ==> automatic choice.",
#endif
	barout_desc[]		= "amount of output for the barrier method:\n\
			0 = no output\n\
			1 = each iteration (default)",
	barpresolve_desc[]	= "level of barrier-specific presolve effort:\n\
			0 = use standard presolve (default)\n\
			1 = use more effort\n\
			2 = do full matrix eliminations for size reduction",
	barprimalstop_desc[]	= "barrier method convergence tolerance on\n\
		primal infeasibilities; default = 0 (automatic choice)",
#ifdef XPRS_BARREFITER
	barrefiter_desc[] =
		"maximum number of refinement iterations, helpful when the\n\
		the solution is near to the optimum using barrier or crossover:\n\
			    0 = default\n\
			n > 0 = perform n refinement iterations",
#endif
#ifdef XPRS_BARREGULARIZE
	barreg_desc[]		= "regularization to use with \"barrier\":\n\
			-1 = automatic choice (default with just \"barrier\")\n\
			Values >= 0 are the sum of:\n\
			1 = use \"standard\" regularization\n\
			2 = use \"reduced\" regularization: less perturbation\n\
				than \"standard\" regularization\n\
			4 = keep dependent rows in the KKT system\n\
			8 = keep degenerate rows in the KKT system",
#endif
	barrier_desc[]		= "[no assignment] use the Newton Barrier algorithm",
	barstart_desc[]		= "choice of starting point for barrier method:\n"
#if XPVERSION >= 30
		"\t\t\t-1 = use incoming solution for warm start\n"
#endif
		"\t\t	 0 = automatic choice (default)\n\
			 1 = heuristics based on magnitudes of matrix entries\n\
			 2 = use pseudoinverse of constraint matrix\n"
#if XPVERSION >= 30
		"\t\t\t 3 = unit starting point for homogeneous self-dual\n\
			    barrier algorithm."
#endif
	,
	barstepstop_desc[]	= "barrier method convergence tolerance: stop when\n\
		step size <= barstepstop; default = 1e-10",
	barthreads_desc[]	= "number of threads used in the Newton Barrier algorithm;\n\
		default = -1 (determined by \"threads\")",
	bestbound_desc[] = "[no assignment] return suffix .bestbound for the best known\n\
		bound on the objective value.  The suffix is on the problem\n\
		and objective and is +Infinity for minimization problems and\n\
		-Infinity for maximization problems if there are no integer\n\
		variables or if an integer feasible solution has not yet\n\
		been found.",
	bigmmethod_desc[]	= "0 = phase I/II, 1 = BigM method (default)",
	branchchoice_desc[]	= "whether to explore branch with min. or max. estimate first:\n\
			0 = explore branch with min. estimate first (default)\n\
			1 = explore branch with max. estimate first\n\
			2 = if an incumbent solution exists, first explore\n\
				the branch satisfied by the incumbent;\n\
				otherwise use choice 0 (min. est. first)"
#if XPVERSION >= 22
		"\n\
			3 (default) = explore the first branch that moves the\n\
				branching variable away from its value at the\n\
				root node; if the branching entity is not a\n\
				simple variable, assume branchchoice=0"
#endif
			,
	branchdisj_desc[]	= "whether to branch on general split disjunctions while\n\
		solving MIPs:\n\
			-1 = automatic choice (default)\n\
			 0 = disabled\n\
			 1 = cautious strategy: create branches only for\n\
				general integers with a wide range\n\
			 2 = moderate strategy\n\
			 3 = aggressive strategy:  create disjunctive branches\n\
				for both binary and integer variables",
	branchstruct_desc[]	= "whether to search for special structure during branch and\n\
		bound:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes",
	breadthfirst_desc[]	=
		"number of MIP nodes included in best-first search\n\
		(default 11) before switching to local-first search",
	cachesize_desc[]	= "cache size in Kbytes -- relevant to Newton Barrier:\n\
			-1 = determined automatically\n\
		default = system-dependent (-1 for Intel)",
	choleskyalg_desc[]	= "type of Cholesky factorization used for barrier: sum of\n\
			  1 ==> manual matrix blocking\n\
			  2 ==> single pass with manual blocking\n\
			  4 ==> nonseparable QP relaxation\n\
			  8 ==> manual corrector weight (honor \"16\" bit)\n\
			 16 ==> manual corrector weight \"on\"\n\
			 32 ==> manual refinement\n\
			 64 ==> use preconditioned conjugate gradients\n\
			128 ==> refine with QMR (quasi-minimal residual)\n\
			default = -1 (automatic choice)",
	choleskytol_desc[]	= "zero tolerance for Cholesky pivots in the\n\
		Newton Barrier algorithm; default = 1e-15",
#ifdef XPRS_CLAMPING
	clamping_desc[] =
		"control adjustements of the returned solution values\n\
		 such that they are always within bounds:\n\
			-1 ==> determined automatically\n\
			 0 ==> adjust primal solution to be within\n\
			       primal bounds (default)\n\
			 1 ==> adjust primal slack values to be within\n\
			       primal bounds\n\
			 2 ==> adjust dual solution to be within dual\n\
			       bounds\n\
			 3 ==> adjust reduced costs to be within dual\n\
			       bounds",
#endif
	convexitychk_desc[]	= "whether to check convexity before solving:\n\
			0 = no\n\
			1 = yes (default)",
#ifdef XPRS_PRECONEDECOMP
	conedecomp_desc[]	=
		"whether to decompose regular and rotated cone constraints\n\
		having more than two elements and to use the result in an\n\
		outer approximation:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes, unless the cone variable is fixed by XPRESS's\n\
				presolve\n\
			 2 = yes, even if the cone variable is fixed\n\
			 3 = yes, but only for outer approximations",
#endif
#ifdef XPRS_CORESPERCPU
	corespercpu_desc[]	=
		"number of cores to assume per cpu; default = -1 ==> number\n\
		detected; barrier cache = cachesize / corespercpu",
#endif
	covercuts_desc[]	=
		"for MIPS, the number of rounds of lifted-cover inequalities\n\
		at the top node; default = -1 ==> automatic choice",
#ifdef XPRS_CPUPLATFORM
	cpuplatform_desc[] = "whether the Newton Barrier method should use AVX or SSE2\n\
		instructions on platforms that offer both:\n"
#if XPVERSION >= 30
		"\t\t	-2 = highest supported [Generic, SSE2, AVX, or AVX2]\n\
			-1 = highest deterministic support (default; no AVX2)\n"
#else
		"\t\t	-1 = automatic choice (default)\n"
#endif
		"\t\t	 0 = use generic code: neither AVX nor SSE2\n\
			 1 = use SSE2\n\
			 2 = use AVX"
#if XPVERSION >= 29
			"\n\t\t\t 3 = use AVX2"
#endif
			,
#endif
	cputime_desc[]		= "which times to report when logfile is specified:\n\
			0 = elapsed time (default)\n\
			1 = CPU time\n\
			2 = process time\n\
		You may need to experiment to see how cputime=1 and\n\
		cputime=2 differ (if they do) on your system.",
	crash_desc[]		= "type of simplex crash:\n\
			0 = none\n\
			1 = one-pass search for singletons\n\
			2 = multi-pass search for singletons (default)\n\
			3 = multi-pass search including slacks\n\
			4 = at most 10 passes, only considering slacks\n\
			    at the end\n\
			n = (for n > 10) like 4, but at most n-10 passes",
	crossover_desc[]	= "whether to find a simplex basis after the barrier alg.:\n\
			-1 = automatic choice (default)\n\
			 0 = no crossover\n\
			 1 = primal crossover first\n\
			 2 = dual crossover first",
#ifdef XPRS_CROSSOVERITERLIMIT
	crossoveritlim_desc[] = "limit on crossover iterations after the barrier\n\
		algorithm; default = 2147483645",
#endif
#ifdef XPRS_CROSSOVEROPS
	crossoverops_desc[]	= "bit vector affecting crossover after the barrier\n\
		algorithm:  sum of\n\
			1 = return the barrier solution (rather than the last\n\
			    intermediate solution) when crossover stop early\n\
			2 = skip the second crossover stage\n\
			4 = skip pivots that are \"less numerically reliable\"\n\
			8 = do a slower but more numerically stable crossover",
#endif
#ifdef XPRS_CROSSOVERTHREADS
	crossoverthreads_desc[]	= "limit on threads used during crossover;\n\
		default not specified in the Release 8.2 documentation",
#endif
#ifdef XPRS_CROSSOVERACCURACYTOL
	crossovertol_desc[]	=
		"tolerance (default 1e-6) for deciding whether to adjust the\n\
		relative pivot tolerance during crossover when a new basis\n\
		factorization is necessary.  Errors in the recalculated\n\
		basic solution above this tolerance cause the pivot\n\
		tolerance to be adjusted.",
#endif
	cutdepth_desc[]		= "maximum MIP tree depth at which to generate cuts:\n\
			0  = no cuts\n\
			-1 = automatic choice (default)",
	cutfactor_desc[]	= "limit on number of cuts and cut coefficients added\n\
		while solving MIPs:\n\
			-1 = automatic choice (default)\n\
			 0 = do not add cuts\n\
			 > 0 ==> multiple of number of original constraints",
	cutfreq_desc[]		=
		"MIP cuts are only generated at tree depths that are integer\n\
		multiples of cutfreq; -1 = automatic choice (default)",
	cutselect_desc[]	= "detailed control of cuts at MIP root node:  sum of\n\
			    32 = clique cuts\n\
			    64 = mixed-integer founding (MIR) cuts\n\
			   128 = lifted cover cuts\n\
			  2048 = flow path cuts\n\
			  4096 = implication cuts\n\
			  8192 = automatic lift-and-project strategy\n\
			 16384 = disable cutting from cut rows\n\
			 32768 = lifted GUB cover cuts\n\
			 65536 = zero-half cuts\n\
			131072 = indicator-constraint cuts\n\
			    -1 = all available cuts (default)",
	cutstrategy_desc[]	=
		"how aggressively to generate MIP cuts; more ==> fewer nodes\n\
		but more time per node:\n\
			-1 = automatic choice (default)\n\
			 0 = no cuts\n\
			 1 = conservative strategy\n\
			 2 = moderate strategy\n\
			 3 = aggressive strategy",
	defaultalg_desc[]	=
		"algorithm to use when none of \"barrier\", \"dual\", or \"primal\"\n\
		is specified:\n\
			1 = automatic choice (default)\n\
			2 = dual simplex\n\
			3 = primal simplex\n\
			4 = Newton Barrier",
#ifdef XPRS_DEGRADEFACTOR
	degradefactor_desc[]	=
		"factor to multiply estimated degradations in MIP objective\n\
		value from exploring an unexplored node; default = 1.0",
#endif
	densecollimit_desc[]	=
		"number of nonzeros above which a column is treated as dense\n\
		in the barrier algorithm's Cholesky factorization:\n\
			0 = automatic choice (default)",
	deterministic_desc[]	= "whether a MIP search should be deterministic:\n\
			0 = no\n\
			1 = yes (default)"
#if XPVERSION >= 38
		"\n\t\t\t2 = yes, with opportunistic root LP solve"
#endif
				,
	dual_desc[]		= "[no assignment] use the dual simplex algorithm",
	dualgradient_desc[]	= "dual simplex pricing strategy:\n\
			-1 = automatic choice\n\
			 0 = Devex\n\
			 1 = steepest edge",
#ifdef XPRS_DUALIZE
	dualize_desc[]		= "whether to convert the primal problem to its dual and solve\n\
		the converted problem:\n\
			-1 = automatic choice (default)\n\
			 0 = no: solve primal problem\n\
			 1 = yes: solve dual problem",
#endif /*XPRS_DUALIZE*/
#ifdef XPRS_DUALIZEOPS
	dualizeops_desc[]	= "when solving the dual problem after deriving it from the\n\
		primal, whether to use primal simplex if dual simplex was\n\
		specified and vice versa:\n\
			0 = no\n\
			1 = yes (default)",
#endif
	dualstrategy_desc[]	= "how to remove infeasibilities when re-optimizing\n\
		with the dual algorithm during MIP solves:\n\
			0 = use primal algorithm\n\
			1 = use dual algorithm (default)",
#ifdef XPRS_DUALTHREADS
	dualthreads_desc[]	= "limit on number of threads used by parallel dual simplex,\n\
		overriding \"threads\"; default -1 ==> use \"threads\"",
#endif
	eigenvaltol_desc[]	= "regard the matrix in a quadratic form as indefinite if its\n\
		smallest eigvenalue is < -eigevnaltol; default = 1e-6",
#ifdef XPRS_ELIMFILLIN
	elimfillin_desc[]	= "maximum fillins allowed for a presolve elimination;\n\
		default = 10.\n",
#endif
	elimtol_desc[]		= "Markowitz tolerance for the elimination phase of\n\
		XPRESS's presolve; default = 0.001",
	etatol_desc[] =		"zero tolerance on eta elements; default varies with XPRESS\n\
		version; default = 1e-12 or 1e-13 with some versions.\n\
		Use etatol=? to see the current value.",
	feaspump_desc[]		= "whether to run the Feasibility Pump heuristic at the top\n\
		node during branch-and-bound:  one of\n\
			0 = no (default)\n\
			1 = yes\n\
			2 = only if other heurstics found no integer solution",
#ifdef XPRS_FEASTOLPERTURB
	feastol_perturb_desc[] = "how much a feasible primal basic solution is allowed to\n\
		be perturbed when performing basis changes.  The tolerance\n\
		specified by \"feastol\" is always considered as an upper\n\
		limit for the perturbations; default = 1.0E-06",
#endif
#ifdef XPRS_FEASTOLTARGET
	feastol_target_desc[]	= "feasibility tolerance on constraints for solution refiner\n\
		(see refineops):  if feastol_target > 0 is specified, it is\n\
		used instead of feastol",
#endif
#ifdef XPRS_MAXGLOBALFILESIZE
	globalfilemax_desc[]	= "maximum megabytes for temporary files storing the global\n\
		search tree:  a new file is started if globalfilemax\n\
		megabytes would be exceeded",
#endif
#ifdef XPRS_GLOBALFILELOGINTERVAL
	globalloginterval_desc[] = "seconds between additions to the logfile about, additions\n\
		to the \"global file\", a temporary file written during a\n\
		global search.  Default = 60.",
#endif
	gomcuts_desc[]		= "gomory cuts at root: -1 = automatic choice (default)",
	heureffort_desc[]	= "factor affecting how much work local search heuristics\n\
		should expend.  Default = 1; higher values cause more\n\
		local searches over larger neighborhoods.",


	hdive_rand_desc[]	= "value between 0 and 1 inclusive affecting randomization\n\
		in the diving heuristic:  0 (default) ==> none;\n\
			1 ==> full;\n\
			intermediate values ==> intermediate behavior",
#ifdef XPRS_HEURDIVESOFTROUNDING
	hdive_rounding_desc[]	= "whether to use soft rounding in the MIP diving heuristic\n\
		(to push variables to their bounds via the objective rather\n\
		than fixing them):\n\
			-1 = automatic choice (default)\n\
			 0 = no soft rounding\n\
			 1 = cautious soft rounding\n\
			 2 = aggressive soft rounding",
#endif

	hdive_speed_desc[]	= "controls tradeoff between speed and solution quality\n\
		in the diving heuristic:  an integer between -2 and 3:\n\
			-2 = automatic bias toward quality\n\
			-1 = automatic bias toward speed (default)\n\
			 0 = emphasize quality\n\
			 4 = emphasize speed\n\
			 1-3 = intermediate emphasis",
	hdive_strategy_desc[]	= "strategy for diving heuristic:  integer between -1 and 10:\n\
			-1 = automatic choice (default)\n\
			 0 = do not use the diving heursistic\n\
			1-10 = preset strategies for diving",
	heurdepth_desc[]	=
#if XPVERSION >= 22
		"deprecated:  no longer has any effect:\n\t\t"
#endif
		"maximum depth of branch-and-bound tree search at which to\n\
		apply heuristics; 0 = no heuristics; default = -1",

#ifdef XPRS_HEUREMPHASIS
	heuremphasis_desc[] =
		"epmphasis for the heuristic search for branch and\n\
		bound.  Setting it to 1 gets a gap quicker at the\n\
		expense of time to optimality:\n\
			-1 = default strategy\n\
			 0 = disable heuristics\n\
			 1 = focus on reducing the gap early\n\
			 2 = extremely aggressive heuristics",
#endif
#ifdef XPRS_HEURFORCESPECIALOBJ
	heurforcespecobj_desc[] = "whether to use special objective heuristics on large\n\
		problems and even if an incumbant exists:\n\
			0 = no (default)\n\
			1 = yes.",
#endif
	heurfreq_desc[]		=
		"during branch and bound, heuristics are applied at nodes\n\
		whose depth from the root is zero modulo heurfreq;\n\
		default = -1 (automatic choice)",
	heurmaxsol_desc[]	=
#if XPVERSION >= 22
		"deprecated:  no longer has any effect:\n\t\t"
#endif
		"maximum number of heuristic solutions to find during branch-\n\
		and-bound tree search; default = -1 (automatic choice)",
	heurnodes_desc[]	=
#if XPVERSION >= 22
		"deprecated:  no longer has any effect:\n\t\t"
#endif
		"maximum nodes at which to use heuristics during\n\
		branch-and-bound tree search; default = 1000",
	heurroot_desc[]		= "bit vector controlling local search heuristics to apply at\n\
		the root node:  sum of\n\
			  1 = large-neighborhood search: may be slow, but may\n\
				find solutions far from the incumbent\n\
			  2 = small-neighborhood search about node LP solution\n\
			  4 = small-neighborhood search about integer solutions\n\
			  8 = local search near multiple integer solutions\n\
			 16 = no effect\n"
#if XPVERSION >= 29
		"\t\t\t 32 = local search without an objective; may only be\n\
				done when no feasible solution is available\n"
#if XPVERSION >= 30
		"\t\t\t 64 = local search with an auxiliary objective; may\n\
				be done when no feasible solution is available\n\
		default = 117",
#else
		"\t\tdefault = 53",
#endif
#else
		"\t\tdefault = 17",
#endif
#ifdef XPRS_HEURSEARCHROOTCUTFREQ
	heurrootcutfreq_desc[]	= "how often to run the local search heuristic while\n\
		cutting at the root node:\n\
			-1 ==> automatic choice (default)\n\
			 0 ==> never\n\
			 n > 0 ==> do n cutting rounds between runs of the\n\
				local search heuristic",
#endif
	heursearch_desc[]	= "how often the local search heurstic should be run during\n\
		branch-and-bound:\n\
			-1 = automatic choice (default)\n\
			 0 = never\n\
			 n > 0 ==> every n nodes",
#ifdef XPRS_HEURSTRATEGY
	heurstrategy_desc[] =
#ifdef XPRS_HEUREMPHASIS
	"deprecated, use heuremphasis:\n\t\t"
#endif
	"heuristic strategy for branch and bound: one of\n\
			-1 = automatic choice (default)\n\
			 0 = no heuristics\n\
			 1 = basic heuristics\n\
			 2 = enhanced heuristics\n\
			 3 = extensive heuristics",
#endif
	heurthreads_desc[]	= "number of threads for the root node of\n\
		branch-and-bound:\n\
			-1 = determined from \"threads\" keyword\n\
			 0 = no separate threads (default)\n\
			 n > 0 ==> use n threads",
	heurtree_desc[]		= "heuristics to apply during tree search:  sum of\n\
		the same values as for heurroot; default 17",
	iis_desc[]	=
		"[no assignment] if the problem is infeasible, find an\n\
		Irreducible Independent Set of infeasible constraints\n\
		and return it in suffix .iis.  If changing the bounds\n\
		on just one constraint or variable could remove the\n\
		infeasibility, return suffix .iso with value 1 for\n\
		each such constraint or variable.",
	indlinbigm_desc[]	= "largest \"big M\" value to use in converting indicator\n\
		constraints to regular constraints; default = 1e5.",
	indprelinbigm_desc[]	= "largest \"big M\" value to use in converting indicator\n\
		constraints to regular constraints during XPRESS\n\
		presolve; default = 100.0",
#ifdef XPRS_INPUTTOL
		inputtol_desc[] =
		"tolerance on input elements (default 0.0); any value v where\n\
		abs(v) <= inputtol is treated as 0",
#endif
	invertfreq_desc[]	=
		"maximum simplex iterations before refactoring the basis:\n\
		-1 = automatic choice (default)",
	invertmin_desc[]	=
		"minimum simplex iterations before refactoring the basis:\n\
		default = 3",
	keepbasis_desc[]	= "basis choice for the next LP iteration:\n\
			0 = ignore previous basis\n\
			1 = use previous basis (default)\n\
			2 = use previous basis only if the number of basic\n\
				variables == number of constraints",
	keepnrows_desc[]	=
		"1 (default) if unconstrained rows are to be kept, else 0",
	lazy_desc[]		=
		"whether to regard constraints with nonzero .lazy suffix\n\
		values as lazy (i.e., delayed) constraints if the problem\n\
		is a MIP:\n\
			0 = no\n\
			1 = yes (default)",
	lnpbest_desc[]		=
		"number of global infeasible entities for which to create\n\
		lift-and-project cuts during each round of Gomory cuts\n\
		at the top node; default = 50",
	lnpiterlimit_desc[]	=
		"maximum iterations for each lift-and-project cut;\n\
		default = -1 (automatic choice)",
#ifdef XPRS_LPREFINEITERLIMIT
	lpref_itlim_desc[] = "limit on simplex iterations used by the solution refiner\n\
		(see refineops); default = -1 ==> automatic choice",
#endif
	localchoice_desc[]	= "when to backtrack between two child nodes\n\
		during a \"dive\":\n\
			1 = never backtrack from the first child unless it\n\
				is dropped (i.e., is infeasible or cut off)\n\
			2 = always solve both nodes first\n\
			3 = automatic choice (default)",
	logfile_desc[]		= "name of log file; default = no log file",
#ifdef XPRS_LPFOLDING
	lpfolding_desc[]	= "whether to attempt exploiting symmetries by \"LP Folding\":\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes.\n",
#endif
	lpthreads_desc[]	= "number of threads in concurrent LP solves:\n\
			-1 = determine from \"threads\" keyword (default)\n\
			n > 0 ==> use n threads",
	lpiterlimit_desc[]	=
		"simplex iteration limit; default = 2147483647 = 2^31 - 1",
	lplog_desc[]		=
		"frequency of printing simplex iteration log; default = 100",
	markowitztol_desc[]	=
		"Markowitz tolerance used when factoring the basis matrix\n\
		default = 0.01",
	matrixtol_desc[]	= "zero tolerance on matrix elements; default = 1e-9",
	maxcuttime_desc[]	= "maximum time (CPU seconds) to spend generating cuts\n\
		and reoptimizing; default = 0 ==> no limit",
	maxiis_desc[]		= "maximum number of Irreducible Infeasible Sets to find:\n\
			-1 = no limit (default)\n\
			 0 = none",
	maxim_desc[]		= "[no assignment] force maximization of the objective",
#ifdef XPRS_MAXIMPLIEDBOUND
	maximpliedbound_desc[]	=
		"when preprocessing MIP problems, only use computed bounds\n\
		at most maximpliedbound (default 1e8) in absolute value",
#endif
	maxlocalbt_desc[]	= "max height above current node to look for a local backtrack\n\
		candidate node; default = 1",
	maxmipsol_desc[]	= "maximum number of integer solutions to find:\n\
			0 = no limit (default)",
	maxnode_desc[]		=
		"maximum number of MIP nodes to explore; default = 2147483647",
	maxpagelines_desc[]	= "maximum output lines between page breaks in logfile;\n\
		default = 23",
	maxlogscale_desc[]	= "max log2 of factors used in scaling; must be >= 0 and\n\
		<= 64; default 64",
#if XPVERSION <= 34
#define XPRS_MAXMEMORYSOFT XPRS_MAXMEMORY
#endif
#ifdef XPRS_MAXMEMORYHARD
	maxmemoryhard_desc[]	= "hard limit (integer number of megabytes) on memory\n\
		allocated, causing early termination if exceeded\n\
			0 (default) = no limit",
#endif
#ifdef XPRS_MAXMEMORYSOFT
	maxmemory_desc[]	= "limit (integer number of megabytes) on memory used:\n\
			-1 = automatic choice (default)\n\
			>0 = target megabytes of memory to use",
#endif
#ifdef XPRS_MAXMIPTASKS
	maxmiptasks_desc[]	= "maximum tasks to run in parallel during a MIP solve:\n\
			   -1 ==> use mipthreads\n\
			n > 0 ==> at most n tasks running at once\n\
		For maxmiptasks > 0, branch-and-bound nodes are solved in a\n\
		deterministic way, but the barrier algorithm (if used) may\n\
		cause a nondeterministic MIP solve unless barthreads = 1.",
#endif
#ifdef XPRS_MAXSLAVE
	maxslaves_desc[]	= "how many processors to use in parallel when solving\n\
		mixed-integer problems:  must be at most the number of\n\
		processors available and licensed",
#endif

#ifdef XPRS_MAXSTALLTIME
	maxstalltime_desc[] =
		"maximum time in seconds that the Optimizer will continue to\n\
		search for improving solution after finding a new incumbent:\n\
			    0 ==> no limit (default)\n\
			n > 0 ==> stop after n seconds without a\n\
			          new incumbent (no effet before\n\
			          the first has been found",
#endif
	maxtime_desc[]		= "limit on solution time:  for maxtime=n (an integer),\n\
			n < 0 ==> stop LP or MIP search after -n seconds\n\
			n = 0 ==> no time limit (default)\n\
			n > 0 ==> for MIP problems, stop after n seconds\n\
				  if a feasible solution has been found;\n\
				  otherwise continue until a feasible\n\
				  solution has been found.",
	minim_desc[]		= "[no assignment] force minimization of the objective",
	mipabscutoff_desc[]	=
		"initial MIP cutoff:  ignore MIP nodes with objective values\n\
		worse than mipabscutoff; default = 1e40 for minimization,\n\
		-1e40 for maximization",
	mipabsstop_desc[]	=
		"stop MIP search if abs(MIPOBJVAL - BESTBOUND) <= mipabsstop\n\
		default = 0",
	mipaddcutoff_desc[]	=
		"amount to add to the objective function of the best integer\n\
		solution found to give the new MIP cutoff; default -1e-5",
#ifdef XPRS_MIPCOMPONENTS
	mipcomponents_desc[] =
		"determines whether disconnected components in a MIP should\n\
		be solved as separate MIPs:\n\
			-1 ==> automatic (default)\n\
			 0 ==> disable\n\
			 1 ==> enable",
#endif

#ifdef XPRS_MIPCONCURRENTNODES
	mipconcurrentnodes_desc[] =
		"node limit to choose the winning solve when concurrent\n\
		solves are enabled:\n\
			   -1 ==> automatic (default)\n\
			n > 0 ==> number of nodes to complete",
	mipconcurrentsolves_desc[] =
		"select the number of concurrent solves to start for a MIP:\n\
			   -1 ==> enabled, the number of concurrent solves\n\
				  depends on mipthreads\n\
			 0, 1 ==> disabled (default)\n\
			n > 1 ==> number of concurrent solves = n",
#endif
#ifdef XPRS_MIPDUALREDUCTIONS
	mipdualreductions_desc[] = "kinds of dual reductions allowed during branch and bound:\n\
			0 ==> none\n\
			1 ==> all (default)\n\
			2 ==> restrict dual reductions to continuous variables.\n\
		If poolnbest > 1 is specified, specifying\n\
		mipdualreductions = 2 might be prudent.",
#endif
#ifdef XPRS_MIPKAPPAFREQ
	mipkappafreq_desc[] = "during branch-and-bound, how often to compute\n\
		basis condition numbers:\n\
			0 ==> never (default)\n\
			1 ==> every node\n\
			n > 1 ==> once per node at level n of the\n\
				branch-and-bound tree.\n\
		When mipkappafreq > 0, a final summary shows the number of\n\
		sampled nodes that are\n\
			\"stable\": kappa < 10^7\n\
			\"suspicious\": 10^7 <= kappa < 10^10\n\
			\"unstable\": 10^10 <= kappa < 10^13\n\
			\"ill-posed\": 10^13 <= kappa.\n\
		A \"Kappa attention level\" between 0 and 1 is also reported.\n\
		Condition numbers use the Frobenius norms of the basis\n\
		and its inverse.",
#endif
	miplog_desc[]		= "MIP printing level to logfile (default -100):\n\
			-n = print summary line every n MIP nodes\n\
			 0 = no MIP summary lines\n\
			 1 = only print a summary at the end\n\
			 2 = log each solution found\n\
			 3 = log each node",
	mipops_desc[]		= "MIP solver options:  one of\n\
			0 = traditional primal first phase (default)\n\
			1 = Big M primal first phase\n\
			2 = traditional dual first\n\
			3 = Big M dual first\n\
			4 = always use artificial bounds in dual\n\
			5 = use original basis only when warmstarting\n\
			6 = skip primal bound flips for ranged primals\n\
			7 = also do single-pivot crash\n\
			8 = suppress aggressive dual perturbations",
	mippresolve_desc[]	= "MIP presolve done at each node: sum of\n\
		         1 = reduced-cost fixing\n\
		         2 = logical preprocessing of binary variables\n\
		         4 = ignored; replaced by \"preprobing\"\n\
		         8 = allow changing continuous-variable bounds\n\
		        16 = allow dual reductions\n\
		        32 = allow global tightening of the problem\n\
		        64 = use objective function\n\
		       128 = allow restarting\n\
		       256 = allow use of symmetry\n\
		default = -1 (automatic choice)",
#ifdef XPRS_MIPRAMPUP
	miprampup_desc[]	= "whether to limit the number of parallel tasks\n\
		during the ramp-up phase of the parallel MIP algorithm:\n\
			-1 = automatic choice (default)\n\
			 0 = no: use as many tasks as possible\n\
			 1 = yes, until finished with initial dives",
#endif
#ifdef XPRS_MIPREFINEITERLIMIT
	miprefiterlim_desc[]	=
		"max. simplex iterations per reoptimization in MIP refiner\n\
		when refineops is 2 or 3; default -1 ==> automatic choice",
#endif
	miprelcutoff_desc[]	=
		"fraction of best integer solution found to add to MIP\n\
		cutoff; default 1e-4",
	miprelstop_desc[]	= "stop MIP search if\n\
		abs(MIPOBJVAL - BESTBOUND) < miprelstop * abs(BESTBOUND);\n\
		default = 0.0001",
#ifdef XPRS_MIPRESTART
	miprestart_desc[] =
		"MIP: control strategy for in-tree restarts:\n\
		-1 = determined automatically (default)\n\
		 0 = disable in-tree restarts\n\
		 1 = normal aggressiveness\n\
		 2 = higher aggressiveness",
#endif
#ifdef XPRS_MIPRESTARTGAPTHRESHOLD
	miprestartgapthreshold_desc[]=
		"MIP: initial gap threshold to delay in-tree restart;\n\
		the restart is delayed if the relative gap is below the\n\
		threshold (default 0.02)",
#endif
#ifdef XPRS_MIPRESTARTFACTOR
	miprestartfactor_desc[] =
	"MIP: fine tune initial conditions to trigger an in-tree\n\
		restart; values > 1 increase the aggressiveness, < 1\n\
		decrease it (default 1.0)",
#endif
	mipstartstatus_desc[]	= "whether to use incoming statuses on MIP problems;\n\
		default 1 ==> yes",
	mipstartvalue_desc[]	= "whether to use the specified initial guess (if supplied)\n\
		when solving a MIP problem:\n\
			0 = no\n\
			1 = yes (default)",
#ifdef XPRS_MIPTERMINATIONMETHOD
	mipstop_desc[]		= "how to stop a MIP solve when a time or node limit is\n\
		reached:\n\
			0 = stop tasks as soon as possible (default)\n\
			1 = let currently running tasks finish, but do not\n\
				start new ones",
#endif
#ifdef XPRS_MIPTARGET
	miptarget_desc[]	=
		"initial MIP target objective value; default = 1e40,\n\
		used in some node-selection rules and updated as MIP\n\
		solutions are found",
#endif
	mipthreads_desc[]	= "number of threads to use solving mixed-integer\n\
		programming problems:\n\
			-1 = use \"threads\" keyword (default)\n\
			n > 0 ==> use n threads",
	miptol_desc[]		= "integer feasibility tolerance; default = 5e-6",
#ifdef XPRS_MIPTOLTARGET
	miptoltarget_desc[]	= "value of miptol used for refining equalities on MIP\n\
		problems when refineops is 2 or 3; default = 0",
#endif
#ifdef XPRS_MIQCPALG
	miqcpalg_desc[]		=
		"algorithm for solving mixed-integer problems with quadratic\n\
		or second-order cone constraints:\n\
			-1 = automatic choice (default)\n\
			 0 = barrier algorithm during branch and bound\n\
			 1 = outer approximations during branch and bound",
#endif
#ifdef XPRS_NETSTALLLIMIT
	netstalllimmit_desc[] =
		"limit the number of degenerate pivots of the network\n\
		simplex algorithmm before switching to primal or dual:\n\
			   -1 ==> automatic\n\
			    0 ==> no limit\n\
			n > 0 ==> limit to n iterations",
#endif
	network_desc[]		= "[no assignment] try to find and exploit an embedded network",
	nodefilebias_desc[]	=
#if (XPVERSION) >= 34
		"deprecated and ignored:\n"
#endif
		"a value between 0 and 1 (inclusive) that influences\n\
		operations when \"treememlimit\" (on how much of the\n\
		branch-and-bound tree should be kept in memory) has\n\
		been exceeded:\n\
			  0 ==> compress every node before writing anything to\n\
				the \"nodefile\";\n\
			  1 ==> write nodes to the \"nodefile\" immediately;\n\
		values between 0 and 1 give intermediate behavior;\n\
		default = 0.5",
#ifdef XPRS_NODEPROBINGEFFORT
	nodeprobingeffort_desc[] =
		"effort put into probing during branch and bound; the\n\
		number is used as a multiplier on the default amount of\n\
		work.  Set to 0 to disable node probing; default 1.",
#endif
	nodeselection_desc[]	= "next MIP node control:\n\
			1 = local first:  choose among descendant and sibling\n\
			    nodes if available, else from all outstanding nodes\n\
			2 = best first of all outstanding nodes\n\
			3 = local depth first:  choose among descendant and\n\
			    sibling nodes if available, else from deepest nodes\n\
			4 = best first for breadthfirst nodes, then local first\n\
			5 = pure depth first:  choose among deepest nodes.\n\
		The default is determined from matrix characteristics.",
	objrep_desc[]		= "Whether to replace\n\
			minimize obj: v;\n\
		with\n\
			minimize obj: f(x)\n\
		when variable v appears linearly in exactly one\n\
		constraint of the form\n\
			s.t. c: v >= f(x);\n\
		or\n\
			s.t. c: v == f(x);\n\
		Possible objrep values:\n\
			0 = no\n\
			1 = yes for v >= f(x)\n\
			2 = yes for v == f(x) (default)\n\
			3 = yes in both cases\n\
		For a maximization problem, \"<=\" replaces \">=\".",
#ifdef XPRS_OBJSCALEFACTOR
	objscalefactor_desc[] = "Power of 2 (default 0) by which the objective is scaled.\n\
		Nonzero objscalfactor values override automatic global\n\
		objective scaling.",
#endif
	optimalitytol_desc[]	= "tolerance on reduced cost; default = 1e-6",
#ifdef XPRS_OPTIMALITYTOLTARGET
	opttol_target_desc[]	= "feasibility tolerance on reduced costs for solution refiner\n\
		(see refineops):  default = 0; if opttol_target > 0 is\n\
		specified, it is used instead of optimalitytol.",
#endif
	outlev_desc[]		= "message level:\n\
			1 = all\n\
			2 = information\n\
			3 = warnings & errors only (default)\n\
			4 = errors\n\
			5 = none",
	outputtol_desc[]	= "zero tolerance on print values; default 1e-5",
	penalty_desc[]		= "minimum absolute penalty variable coefficient;\n\
	default = automatic choice",
	param_desc[] =
		"Used with syntax \"param=name=value\" (no spaces), where\n\
		\"name\" is the name of an XPRESS control parameter and\n\
		\"value\" is to be assigned to that parameter.  If value is\n\
		?, report the current value of the parameter.  If name is\n\
		a string control, value can be a quoted string or a\n\
		sequence of nonblank characters other than comma.  This\n\
		facility provides a way to modify control parameters,\n\
		identified by name or number, that have not (yet) been\n\
		assigned a keyword.  As a special case, \"param=?\" requests\n\
		a list of all control parameters and their current values.",
#ifdef XPRS_PREPERMUTESEED
	permuteseed_desc[]	= "seed for the random-number generator used by prepermute;\n\
		default = 1",
#endif
#ifdef XPRS_PERTURB
	perturb_desc[]		= "perturbation factor if autoperturb is set to 1;\n\
		0 = default = automatic choice."
#if XPVERSION >= 33

		"\nDeprecated; overridden by primalperturb and dualperturb,\n"
		"which should be used instead of perturb.",
#endif
#endif
	primalperturb_desc[]	=
		"Factor by which to possibly perturb the problem in the\n"
		"dual primal algorithm.  If >= 0, overrides \"perturb\".\n"
		"Default -1 ==> automatic choice; 0 ==> no perturbatation.",
	dualperturb_desc[]	=
		"Factor by which to possibly perturb the problem in the\n"
		"dual simplex algorithm.  If >= 0, overrides \"perturb\".\n"

#endif
		,		"Default -1 ==> automatic choice; 0 ==> no perturbatation."
	pivottol_desc[]		= "zero tolerance for pivots; default = 1e-9",
#ifdef XPRS_MSP_SOLPRB_OBJ /*{*/
/*		"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"	*/
	pooldualred_desc[] =
		"Whether to suppress removal of dominated solutions (via\n"
		"\"dual reductions\") when poolstub is specified:\n\
			0 = yes (default, which can be expensive)\n\
			1 = no\n\
			2 = honor presolveops bit 3 (2^3 = 8)",
	pooldupcol_desc[] =
		"Whether to suppress duplicate variable removal when\n"
		"poolstub is specified:\n\
			0 = yes (default, which can be expensive)\n\
			1 = no\n\
			2 = honor presolveops bit 5 (2^5 = 32)",
	pooldups_desc[] =
		"How poolstub should handle duplicate solutions:\n\
			0 = retain all duplicates\n\
			1 = discard exact matches\n\
			2 = discard exact matches of continuous variables\n\
				and matches of rounded values of discrete\n\
				variables\n\
			3 = default: discard matches of rounded values of\n\
				discrete variables\n"
		"Rounding of discrete variables is affected by poolmiptol\n"
		"and poolfeastol.",
	poolfeastol_desc[] =
		"Zero tolerance for discrete variables in the solution\n"
		"pool (see poolstub); default = 1e-6.",
	poolmiptol_desc[] =
		"Error (nonintegrality) allowed in discrete variables\n"
		"in the solution pool (see poolstub); default = 5e-6.",
	poolnbest_desc[] =
		"Whether the solution pool (see poolstub) should contain\n"
		"inferior solutions.  When poolnbest = n > 1, the\n"
		"solution pool is allowed to keep the n best solutions.",
	poolstub_desc[] =
		"Stub for solution files in the MIP solution pool.\n"
		"Ignored unless some variables are integer or binary.\n"
		"A pool of alternate MIP solutions is computed if\n"
		"poolstub is specified, and the solutions in this pool\n"
		"are written to files\n\n\
		   (poolstub & '1') ... (poolstub & |solution pool|),\n\n"
		"where |solution pool| is the number of solutions in the\n"
		"solution pool.  That is, file names are obtained by\n"
		"appending 1, 2, ... |solution pool| to poolstub.  The\n"
		"value of |solution pool| is returned in suffix npool\n"
		"on the objective and problem.",
#endif /*}*/
	ppfactor_desc[]		= "partial-pricing candidate-list size factor; default = 1.0",
#ifdef XPRS_PREANALYTICCENTER
	preanalyticcenter_desc[] = "whether to compute and use analytic centers while solving\n\
		MIP problems:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes, but only for variable fixing\n\
			 2 = yes, but only for computing reduced costs\n\
			 3 = yes, for both variable fixing and reduced costs.",
#endif
#ifdef XPRS_PREBASISRED
	prebasisred_desc[]	= "whether XPRESS's presolve should try to use a lattice basis\n\
		reduction algorithm:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes.",
#endif
#ifdef XPRS_PREBNDREDCONE
	prebndredcone_desc[]	= "for MIP problems, whether to use cone constraints to\n\
		reduce bounds on variables:\n\
			 0 = no\n\
			 1 = yes\n\
			-1 = default (undocumented)",
#endif
#ifdef XPRS_PREBNDREDQUAD
	prebndredquad_desc[]	= "for MIP problems, whether to use convex quadratic\n\
		constraints to reduce bounds on variables:\n\
			 0 = no\n\
			 1 = yes\n\
			-1 = default (undocumented)",
#endif
	precoefelim_desc[]	= "whether XPRESS's presolve should recombine constraints:\n\
			0 = no,\n\
			1 = yes, as many as possible\n\
			2 = yes, cautiously (default)",
#ifdef XPRS_PRECOMPONENTS
	precomponents_desc[]	= "whether XPRESS's presolve should detect and separately\n\
		solve independent MIP subproblems:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes",
#endif
#ifdef XPRS_PRECONVERTSEPARABLE
	preconvertsep_desc[] = "How to reformulate problems with nondiagonal quadratic\n\
		objectives or constraints:\n\
			-1 = automatic choice (default)\n\
			 0 = no reformulation\n\
			 1 = reformulate to diagonal constraints\n\
			 2 = also allow reduction to second-order cones\n\
			 3 = also convert the objective to a constraint.",
#endif
	predomcol_desc[]	= "whether XPRESS's presolve should remove variables\n\
		when solving MIP problems:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes, cautiously\n\
			 2 = yes, check all candidates",
	predomrow_desc[]	= "whether XPRESS's presolve should remove constraints\n\
		when solving MIP problems:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes, cautiously\n\
			 2 = yes, medium strategy\n\
			 3 = yes, check all candidates",
#ifdef XPRS_PREDUPROW
	preduprow_desc[] = "how XPRESS's presolve should deal with duplicate rows\n\
		in MIP problems:\n\
			-1 = automatic choice (default)\n\
			 0 = do not remove duplicate rows (constraints)\n\
			 1 = remove duplicate rows identical in all variables\n\
			 2 = like 1 but allowing simple penalty variables\n\
			 3 = like 1 but allowing more complex penalty variables",
#endif
#ifdef XPRS_PREFOLDING
	prefolding_desc[] =
		"choose if folding aggregate continuous column in an\n\
		equitable partition:\n\
			-1 = automatic choiche (default)\n\
			 0 = disabled\n\
			 1 = enabled",
#endif
#ifdef XPRS_PREIMPLICATIONS
	preimplications_desc[] = "whether XPRESS's presolve should use implication\n\
		structures to remove redundant rows:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes",
#endif
#ifdef XPRS_PRELINDEP
	prelindep_desc[]	= "whether to check for and remove linearly dependent\n\
		equality constraints:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes",
#endif
#ifdef XPRS_PREOBJCUTDETECT
	preobjcutdetect_desc[]	= "on MIP problems, whether to check for constraints\n\
		that are (nearly) parallel to a linear objective function\n\
		and can be removed safely:\n\
			0 = no\n\
			1 = yes (default)",
#endif
#ifdef XPRS_PREPERMUTE
	prepermute_desc[]	= "whether to randomly permute variables or constraints before\n\
		applying XPRESS's presolve:  sum of\n\
			1 ==> permute constraints\n\
			2 ==> permute variables\n\
			4 ==> permute global MIP information\n\
		default = 0; see permuteseed",
#endif
	preprobing_desc[]	= "how much probing on binary variables to do during XPRESS's\n\
		presolve:\n\
			-1 = automatic choice (default)\n\
			 0 = none\n\
			 1 = light probing\n\
			 2 = full probing\n\
			 3 = repeated full probing",
	presolve_desc[]		= "whether to use XPRESS's presolver:\n\
			0 = no\n\
			1 = yes, removing redundant bounds (default)\n\
			2 = yes, retaining redundant bounds",
#ifdef XPRS_PRESOLVEMAXGROW
	presolvemaxgrow_desc[] = "factor by which the number of nonzero coefficients\n\
		may grow during XPRESS's presolve; default = 0.1",
#endif
	presolveops_desc[]	= "reductions to use in XPRESS's presolve:  sum of\n\
			    1 = 2^0  = remove singleton columns\n\
			    2 = 2^1  = remove singleton constraints (rows)\n\
			    4 = 2^2  = forcing row removal (whatever that is)\n\
			    8 = 2^3  = dual reductions\n\
			   16 = 2^4  = redundant constraint (row) removal\n\
			   32 = 2^5  = duplicate variable removal\n\
			   64 = 2^6  = duplicate constraint removal\n\
			  128 = 2^7  = strong dual reductions\n\
			  256 = 2^8  = variable eliminations\n\
			  512 = 2^9  = no IP reductions\n\
			 1024 = 2^10 = no semicontinuous variable detection\n\
			 2048 = 2^11 = no advanced IP reductions\n\
			16384 = 2^14 = remove linearly dependent constraints\n\
			32768 = 2^15 = no integer variable and SOS detection\n\
		default = 511 (bits 0-8 set).",
#ifdef XPRS_PRESOLVEPASSES
	presolvepasses_desc[] = "Number of rounds to use in the XPRESS presolve algorithm;\n"
			"default = 1.",
#endif
	pricingalg_desc[]	= "primal simplex pricing method:\n\
			-1 = partial pricing\n\
			 0 = automatic choice (default)\n\
			 1 = Devex pricing",
	primal_desc[]		= "[no assignment] use the primal simplex algorithm",
#ifdef XPRS_PRIMALOPS
	primalops_desc[]	= "primal simplex options:  sum of\n\
			1 = 2^0 = aggressive dj scaling\n\
			2 = 2^1 = conventional dj scaling\n\
			4 = 2^2 = reluctant switching back to partial pricing\n\
			8 = 2^3 = dynamic switching between cheap and expensive pricing\n\
			default = all of the above; if bits 0 and 1 are the same (both on or\n\
			both off), choose dj scaling automatically",
#endif
	primalunshift_desc[]	= "whether the primal alg. calls the dual to unshift:\n\
			0 = yes (default)\n\
			1 = no",
	pseudocost_desc[]	=
		"default pseudo-cost assumed for forcing an integer variable\n\
		to an integer value; default = 0.01",
	pseudocost_ud_desc[] = 	"how to update pseudocosts during branch-and-bound:\n\
			-1 = automatic choice (default)\n\
			 0 = no updates\n\
			 1 = use only regular branches\n\
			 2 = use regular and strong branch results\n\
			 3 = use results from all nodes",
#ifdef XPRS_QCCUTS
	qccuts_desc[] = "when using miqcpalg=1 to solve a mixed-integer problem that\n\
		has quadratic constraints or second-order cone constraints,\n\
		the number of rounds of outer approximation cuts at the top\n\
		node:  default = -1 means automatic choice.",
	qcrootalg_desc[] = "when using miqcpalg=1 to solve a mixed-integer problem that\n\
		has quadratic constraints or second-order cone constraints,\n\
		the algorithm for solving the root node:\n\
			-1 = automatic choice (default)\n\
			 0 = barrier algorithm\n\
			 1 = dual simplex on outer approximations",
#endif
	quadunshift_desc[]	= "whether quadratic simplex should do an extra\n\
			purification after finding a solution:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes",
	ray_desc[]		=
		"whether to return a ray of unboundedness in suffix .unbdd:\n\
			0 ==> no (default)\n\
			1 ==> yes, after suppressing XPRESS's presolve\n\
			2 ==> yes, without suppressing XPRESS's presolve\n\
		The last setting (ray=2) may give wrong results when\n\
		XPRESS's presolve detects infeasibility.  Both ray=1 and\n\
		ray=2 cause reoptimization with primal simplex if some other\n\
		algorithm was used.  No ray is returned for MIP problems.",
#ifdef XPRS_REFINEOPS
	refineops_desc[]	= "whether refine equalities -- to reduce infeasibilities\n\
		in constraints that should hold as equalities: sum of\n\
			 1 ==> refine LP solutions\n\
			 2 ==> refine MIP solutions;\n"
#if XPVERSION <= 33
		"\t\tdefault = 3 (do both)",
#else
		"\t\t	 4 ==> refine the final MIP solution found\n\
			 8 ==> refine each node of the search tree\n\
			16 ==> refine non-global solutions\n"
#if XPVERSION > 34
		"\t\t	32 ==> refine all solutions\n\
			64 ==> use higher precision during iterative refinement\n"
#endif
#if XPVERSION >= 36
		"\t    128 ==> use the primal simplex algorithm for refining\n\
	    256 ==> use the dual simplex algorithm for refining\n\
	    512 ==> refine MIP solutions such that rounding them\n\
				keeps the problem feasible when reoptimized\n\
	   1024 ==> attempt to refine MIP solutions such that\n\
				rounding them keeps the problem feasible when\n\
				reoptimized, but accept integers solutions\n\
				even if refinement fails.\n"
#endif
		"\t\tdefault = 1 + 2 + 16 = 19.",
#endif
#endif
	relaxtreemem_desc[]	= "fraction of memory limit by which to relax \"treememlimit\"\n\
		when too much structural data appears; default 0.1",
	relpivottol_desc[]	= "relative pivot tolerance default = 1e-6",
	repairindefq_desc[]	= "whether to repair indefinite quadratic forms:\n\
			0 = yes\n\
			1 = no (default)",
#ifdef XPRS_RESOURCESTRATEGY
	resourcestrategy_desc[] = "whether to allow nondeterministic decisions to cope with\n\
		low memory (affected by maxmemory and maxmemoryhard):\n\
			0 = no (default)\n\
			1 = yes",
#endif
	rootpresolve_desc[]	= "whether to presolve after root cutting and heuristics:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes",
	round_desc[]		=
		"whether to round integer variables to integral values before\n\
		returning the solution, and whether to report that XPRESS\n\
		returned noninteger values for integer values:  sum of\n\
			 1 ==> round nonintegral integer variables\n\
			 2 ==> do not modify solve_result\n\
			 4 ==> do not modify solve_message\n\
			 8 ==> report modifications even if maxerr < 1e-9\n\
		Modifications take place only if XPRESS assigned nonintegral\n\
		values to one or more integer variables, and (for round < 8)\n\
		are reported if the maximum deviation from integrality\n\
		exceeded 1e-9.  Default = 1.",
	sbbest_desc[]		= "For MIP problems, the number of infeasible\n\
		global entities on which to perform strong branching;\n\
		default -1 ==> automatic choice.",
	sbeffort_desc[]		= "multiplier on strong-branching controls that\n\
		are set to \"automatic\"; default = 1.0",
	sbestimate_desc[]	= "how to compute pseudo costs from the local node\n\
		when selecting an infeasible entity to branch on:\n\
			-1 = automatic choice (default)\n\
			1-6 = particular strategies (not described)",
	sbiterlimit_desc[]	=
		"Number of dual iterations to perform the strong branching;\n\
		0 ==> none; default = -1 (automatic choice)",
	sbselect_desc[]		= "size of candidate list for strong branching:\n\
			-2 = low-effort automatic choice (default)\n\
			-1 = high-effort automatic choice\n\
			n >= 0 ==> include max(n, sbbest) candidates",
	scaling_desc[]		=
		"how to scale the constraint matrix before optimizing: sum of\n\
			    1 = 2^0 = row scaling\n\
			    2 = 2^1 = column scaling\n\
			    4 = 2^2 = row scaling again\n\
			    8 = 2^3 = maximum scaling\n\
			   16 = 2^4 = Curtis-Reid\n\
			   32 = 2^5 = scale by maximum element (rather\n\
					than by geometric mean)\n\
			   64 = 2^6 = no special handing for big-M constraints\n\
			  128 = 2^7 = objective-function scaling\n\
			  256 = 2^8 = excluding quadratic part of constraint\n\
					when calculating scaling factors\n\
			  512 = 2^9 = scale before presolve\n\
			 1024 = 2^10 = do not scale constraints (rows) up\n\
			 2048 = 2^11 = do not scale variables down"
#if XPVERSION >= 23
				"\n\
			 4096 = 2^12 = do global objective function scaling\n\
			 8192 = 2^13 = do right-hand side scaling\n\
			16384 = 2^14 = disable aggressive quadratic scaling\n\
			32768 = 2^15 = disable explicit slack scaling."
#endif
				"\n\
		Default = 163.",
#ifdef XPRS_SIFTPASSES
	siftpasses_desc[] =
		"how quickly we allow the worker problems to grow during the\n\
		sifting algorithm; large values might reduce the number of\n\
		iterations but increase the solve time for each.  Default 4.",
#endif
#ifdef XPRS_SIFTING
	sifting_desc[] = "when using dual simplex, whether to enable sifting,\n\
		which can speed up the solve when there are many more\n\
		variables than constraints:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes",
#endif
#ifdef XPRS_SIFTPRESOLVEOPS
	siftpresolveops_desc[] =
		"presolve operations for solving the subproblems during\n\
		sifting:\n\
			 -1 = use presolveops value (default)\n\
			> 0 = use this value",
#endif
#ifdef XPRS_SIFTSWITCH
			siftswitch_desc[] =
			"determines which algoorithm to use during sifting\n\
			   -1 ==> dual simplex\n\
			    0 ==> barrier\n\
			n > 0 ==> barrier if the number of dual\n\
				  infeasibilities > n else dual simplex",
#endif
#ifdef XPRS_SLEEPONTHREADWAIT
	sleeponthreadwait_desc[]	=
		"whether threads should sleep while awaiting work:\n\
			0 = no (busy-wait)\n\
			1 = yes (sleep; may add overhead)\n\
		default = -1 (automatic choice)",
#endif
	sos_desc[]		= "whether to use explicit SOS information; default 1 ==> yes",
	sos2_desc[]		= "whether to tell XPRESS about SOS2 constraints for\n\
		nonconvex piecewise-linear terms; default 1 ==> yes",
	sosreftol_desc[]	= "minimum relative gap between reference row entries;\n\
		default = 1e-6",
#ifdef XPRS_SYMMETRY
	symmetry_desc[]		= "amount of effort to detect symmetry in MIP problems:\n\
			0 = none: do not attempt symmetry detection\n\
			1 = modest effort (default)\n\
			2 = aggressive effort",
#endif
#ifdef XPRS_TEMPBOUNDS
	tempbounds_desc[]	= "whether dual simplex should put temporary bounds on\n\
		unbounded variables:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes",
#endif
	threads_desc[]		= "default number of threads to use:\n\
			-1 = automatic choice (based on hardware)\n\
			 n > 0 ==> use n threads",
#if XPVERSION >= 30 /*{*/
	tunerdir_desc[]		= "directory for tuner results; specifying tunerdir causes\n\
		the XPRESS tuner to solve the problem several times\n\
		to find good settings for solving similar problems.\n\
		Results are stored in tunerdir and its subdirectories.",
	tunerhistory_desc[]	= "when tunerdir is specified, whether to reuse previous\n\
		tuner results and/or to augment them:\n\
			0 = discard previous tuner results\n\
			1 = ignore previous tuner results,\n\
				but add new results to them\n\
			2 = reuse previous tuner results and add\n\
				new results to them (default).",
	tunermaxtime_desc[]	= "maximum seconds to run the tuner when tunerdir is\n\
		specified.  Default 0 ==> no limit.  Use \"maxtime\" to limit\n\
		the time the tuner uses for each problem solved.",
	tunermethod_desc[]	= "method for tuning when tunerdir is specified:\n\
			-1 = automatic choice (default)\n\
			 0 = default LP tuner\n\
			 1 = default MIP tuner\n\
			 2 = more elaborate MIP tuner\n\
			 3 = root-focused MIP tuner\n\
			 4 = tree-focused MIP tuner\n\
			 5 = simple MIP tuner\n\
			 6 = default SLP tuner\n\
			 7 = default MISLP tuner"
#if XPVERSION >= 34
			"\n\
			 8 = MIP tuner using primal heuristics."
#endif
			,
	tunermethodfile_desc[]	= "name of a file that can be read to specify the\n\
		method for tuning (overriding tunermethod) when tunerdir\n\
		is specified.",
	tunerpermute_desc[]	= "when running the XPRESS tuner and tunerpermute = n > 0,\n\
		solve the original problem and n permutations thereof.",
	tunertarget_desc[]	= "what to measure to compare two problem solutions\n\
		when running the XPRESS tuner (what to measure):\n\
			-1 = automatic choice (default)\n\
			 0 = solution time, then integrality gap\n\
			 1 = solution time, then best bound\n\
			 2 = solution time, then best integer solution\n\
			 3 = the \"primal dual integral\", whatever that is\n\
			 4 = just solution time (default for LPs)\n\
			 5 = just objective value\n\
			 6 = validation number (probably not relevant)"
#if XPVERSION >= 32
			"\n\
			 7 = integrality gap only\n\
			 8 = best bound only\n\
			 9 = best integer solution only"
#endif
			".",
	tunerthreads_desc[] =	"number of tuner threads to run in parallel:\n\
		default -1 ==> automatic choice.\n\
		\"threads\" controls the number of threads for each solve.\n\
		The product of threads and tunerthreads should not exceed\n\
		the number of threads the system can run in parallel.",
#endif /*}*/
	trace_desc[]		= "whether to explain infeasibility:\n\
			0 = no (default)\n\
			1 = yes",
	treecovercuts_desc[]	=
		"number of rounds of lifted-cover inequalities at MIP nodes\n\
		other than the top node (cf covercuts);\n\
		default = -1 (automatic choice)",
	treecompress_desc[]	= "level of effort at data compression when branch-and-bound\n\
		memory exceeds \"treememlimit\":  higher ==> greater effort\n\
		(taking more time); default = 2",
	treecuts_desc[] 	= "cuts to generate at nodes during tree search:  sum of\n\
			    32 = 2^5  = clique cuts\n\
			    64 = 2^6  = mixed-integer rounding (MIR) cuts\n\
			    64 = 2^7  = lifted-cover cuts\n\
			  2048 = 2^11 = flow-path cuts\n\
			  4096 = 2^12 = implication cuts\n\
			  8192 = 2^13 = lift-and-project cuts\n\
			 16384 = 2^14 = disable cutting from row cuts\n\
			 32768 = 2^15 = lifted GUB cover cuts\n\
			 65536 = 2^16 = zero-half cuts\n\
			131072 = 2^17 = indicator cuts.\n\
		Default = 259839 (same effect as -2305).",
	treegomcuts_desc[]	=
		"number of rounds of Gomory cuts to generate at MIP nodes\n\
		other than the top node (cf covercuts);\n\
		default = -1 (automatic choice)",
	treeoutlev_desc[]	= "how much to report about branch-and-bound trees\n\
		(if allowed by outlev):  sum of\n\
			1 = regular summaries\n\
			2 = report tree compression and output to nodefile\n\
		default = 3",
#ifdef XPRS_TREEPRESOLVE
	treepresolve_desc[] = "how much presolving to apply to nodes of the MIP\n\
		branch-and-bound tree:\n\
			-1 = automatic choice (default)\n\
			 0 = none\n\
			 1 = cautious\n\
			 2 = moderate\n\
			 3 = aggressive",
#endif
	treememlimit_desc[]	= "an integer: soft limit in megabytes on memory to use for\n\
		branch-and-bound trees.  Default = 0 ==> automatic choice.",
	treememtarget_desc[] 	= "fraction of \"treememlimit\" to try to recover by compression\n\
		or writing to nodefile when  \"treememlimit\" is exceeded.\n\
		Default = 0.2",
	varselection_desc[]	=
		"how to score the integer variables at a MIP node, for\n\
			branching on a variable with minimum score:\n\
			-1 = automatic choice (default)\n\
			 1 = minimum of the 'up' and 'down' pseudo-costs\n\
			 2 = 'up' pseudo-cost + 'down' pseudo-cost\n\
			 3 = maximum of the 'up' and 'down' pseudo-costs plus\n\
			     twice their minimum\n\
			 4 = maximum of the 'up' and 'down' pseudo-costs\n\
			 5 = the 'down' pseudo-cost\n\
			 6 = the 'up' pseudo-cost",
	version_desc[] =
		"Report version details before solving the problem.  This is\n\
		a single-word \"phrase\" that does not accept a value\n\
		assignment.",
	writeprob_desc[] = "Name of file to which the problem is written\n\
		in a format determined by the name's suffix:\n\
			.mps = MPS file;\n\
			.lp = LP file.",
	/* WS_desc = variant of solvers/ws_desc.c without "=..." and with extra tabs */
	WS_desc[] = "solution report without -AMPL: sum of\n\
			1 = write .sol file\n\
			2 = print primal variable values\n\
			4 = print dual variable values\n\
			8 = do not print solution message";

#ifdef XPRS_CONCURRENTTHREADS
#undef XPRS_LPTHREADS
#define XPRS_LPTHREADS XPRS_CONCURRENTTHREADS
#endif
        /* The list of Xpress-MP options */
        /******MUST BE IN ALPHABETIC ORDER!****/

static keyword keywds[]={
  KW("advance",		I_val,	 &advance,		advance_desc),
#ifdef XPRS_ALGAFTERCROSSOVER
  KW("algaftercrossover", set_int, XPRS_ALGAFTERCROSSOVER, algaftercrossover_desc),
#endif
#ifdef XPRS_ALGAFTERNETWORK
  KW("algafternetwork",	set_int, XPRS_ALGAFTERNETWORK, algafternetwork_desc),
#endif
#if XPVERSION >= 30
  KW("archconsistent",	I_val,	&archconsistent,	archconsistent_desc),
#endif
#ifdef XPRS_AUTOCUTTING
  KW("autocutting",	set_int, XPRS_AUTOCUTTING,	autocutting_desc),
#endif
  KW("autoperturb",	set_int, XPRS_AUTOPERTURB,	autoperturb_desc),
#ifdef XPRS_AUTOSCALING
	KW("autoscaling", set_int, XPRS_AUTOSCALING, autoscaling_desc),
#endif
  KW("backtrack",	set_int, XPRS_BACKTRACK,	backtrack_desc),
  KW("backtracktie",	set_int, XPRS_BACKTRACKTIE,	backtracktie_desc),
#ifdef XPRS_BARALG /*{*/
  KW("baralg",		set_int, XPRS_BARALG,		baralg_desc),
#ifdef XPRS_BARCORES
  KW("barcores",	set_int, XPRS_BARCORES,		barcores_desc),
#endif
#endif /*}*/
  KW("barcrash",	set_int, XPRS_BARCRASH,		barcrash_desc),
  KW("bardualstop",	set_dbl, XPRS_BARDUALSTOP,	bardualstop_desc),
  KW("bargapstop",	set_dbl, XPRS_BARGAPSTOP,	Bstop_desc),
#ifdef XPRS_BARGAPTARGET
  KW("bargaptarget",	set_dbl, XPRS_BARGAPTARGET,	bargaptarget_desc),
#endif
  KW("barindeflimit",	set_int, XPRS_BARINDEFLIMIT,	barindeflimit_desc),
  KW("bariterlimit",	set_int, XPRS_BARITERLIMIT,	bariterlimit_desc),
#ifdef XPRS_BARKERNEL
  KW("barkernel",	set_dbl, XPRS_BARKERNEL,	barkernel_desc),
#endif
#ifdef XPRS_BAROBJPERTURB
  KW("barobjperturb", set_dbl, XPRS_BAROBJPERTURB, barobjperturb_desc),
#endif
  KW("barobjscale",	set_dbl, XPRS_BAROBJSCALE,	barobjscale_desc),
  KW("barorder",	set_int, XPRS_BARORDER,		barorder_desc),
#ifdef XPRS_BARORDERTHREADS
  KW("barorderthreads",	set_int, XPRS_BARORDERTHREADS,	barorderthreads_desc),
#endif
  KW("baroutput",	set_int, XPRS_BAROUTPUT,	barout_desc),
  KW("barpresolve",	set_int, XPRS_BARPRESOLVEOPS,	barpresolve_desc),
  KW("barprimalstop",	set_dbl, XPRS_BARPRIMALSTOP,	barprimalstop_desc),
#ifdef XPRS_BARREFITER
  KW("barrefiter",	set_int, XPRS_BARREFITER,	barrefiter_desc),
#endif
#ifdef XPRS_BARREGULARIZE
  KW("barreg",		set_int, XPRS_BARREGULARIZE,	barreg_desc),
#endif
  KW("barrier",		set_known, set_barrier,		barrier_desc),
  KW("barstart",	set_int, XPRS_BARSTART,		barstart_desc),
  KW("barstepstop",	set_dbl, XPRS_BARSTEPSTOP,	barstepstop_desc),
  KW("barthreads",	set_int, XPRS_BARTHREADS,	barthreads_desc),
  KW("basisin",		set_fln, &startbasis,		"load initial basis from specified file"),
  KW("basisout",	set_fln,  &endbasis,		"save final basis to specified file"),
  KW("bestbound",	set_known, set_bestbound,	bestbound_desc),
  KW("bigm",		set_dbl, XPRS_BIGM,		"infeasibility penalty; default = 1024"),
  KW("bigmmethod",	set_int, XPRS_BIGMMETHOD,	bigmmethod_desc),
  KW("branchchoice",	set_int, XPRS_BRANCHCHOICE,	branchchoice_desc),
  KW("branchdisj",	set_int, XPRS_BRANCHDISJ,	branchdisj_desc),
  KW("branchstruct",	set_int, XPRS_BRANCHSTRUCTURAL, branchstruct_desc),
  KW("breadthfirst",	set_int, XPRS_BREADTHFIRST,	breadthfirst_desc),
  KW("cachesize",	set_int, XPRS_CACHESIZE,	cachesize_desc),
  KW("choleskyalg",	set_int, XPRS_CHOLESKYALG,	choleskyalg_desc),
  KW("choleskytol",	set_dbl, XPRS_CHOLESKYTOL,	choleskytol_desc),
#ifdef XPRS_CLAMPING
  KW("clamping",	set_int, XPRS_CLAMPING,	clamping_desc),
#endif
  KW("concurrentthreads", set_int, XPRS_LPTHREADS,	"synonym for lpthreads"),
#ifdef XPRS_PRECONEDECOMP
  KW("conedecomp",	set_int, XPRS_PRECONEDECOMP,	conedecomp_desc),
#endif
  KW("convexitychk",	set_int, XPRS_IFCHECKCONVEXITY,	convexitychk_desc),
#ifdef XPRS_CORESPERCPU
  KW("corespercpu",	set_int, XPRS_CORESPERCPU,	corespercpu_desc),
#endif
  KW("covercuts",	set_int, XPRS_COVERCUTS,	covercuts_desc),
#ifdef XPRS_CPUPLATFORM
  KW("cpuplatform",	set_int, XPRS_CPUPLATFORM,	cpuplatform_desc),
#endif
  KW("cputime",		set_int, XPRS_CPUTIME,		cputime_desc),
  KW("crash",		set_int, XPRS_CRASH,		crash_desc),
  KW("crossover",	set_int, XPRS_CROSSOVER,	crossover_desc),
#ifdef XPRS_CROSSOVERITERLIMIT
  KW("crossoveritlim", set_int, XPRS_CROSSOVERITERLIMIT, crossoveritlim_desc),
#endif
#ifdef XPRS_CROSSOVEROPS
  KW("crossoverops",	set_int, XPRS_CROSSOVEROPS,	crossoverops_desc),
#endif
#ifdef XPRS_CROSSOVERTHREADS
  KW("crossoverthreads", set_int, XPRS_CROSSOVERTHREADS, crossoverthreads_desc),
#endif
#ifdef XPRS_CROSSOVERACCURACYTOL
  KW("crossovertol",	set_dbl, XPRS_CROSSOVERACCURACYTOL, crossovertol_desc),
#endif
  KW("cutdepth",	set_int, XPRS_CUTDEPTH,		cutdepth_desc),
  KW("cutfactor",	set_dbl, XPRS_CUTFACTOR,	cutfactor_desc),
  KW("cutfreq",		set_int, XPRS_CUTFREQ,		cutfreq_desc),
  KW("cutselect",	set_int, XPRS_CUTSELECT,	cutselect_desc),
  KW("cutstrategy",	set_int, XPRS_CUTSTRATEGY,	cutstrategy_desc),
#ifdef RWA_DEBUG
  KW("debug",		set_known, set_debug,		"RWA's debug switch [no assignment]"),
#endif
  KW("defaultalg",	set_int, XPRS_DEFAULTALG,	defaultalg_desc),
#ifdef XPRS_DEGRADEFACTOR
  KW("degradefactor",	set_dbl, XPRS_DEGRADEFACTOR,	degradefactor_desc),
#endif
  KW("densecollimit",	set_int, XPRS_DENSECOLLIMIT,	densecollimit_desc),
  KW("deterministic",	set_int, XPRS_DETERMINISTIC,	deterministic_desc),
  KW("dual",		set_known, set_dual,		dual_desc),
  KW("dualgradient",	set_int, XPRS_DUALGRADIENT,	dualgradient_desc),
#ifdef XPRS_DUALIZE
  KW("dualize",		set_int, XPRS_DUALIZE,		dualize_desc),
#endif
#ifdef XPRS_DUALIZEOPS
  KW("dualizeops",	set_int, XPRS_DUALIZEOPS,	dualizeops_desc),
#endif
#if XPVERSION >= 33
  KW("dualperturb",	set_dbl, XPRS_DUALPERTURB,	dualperturb_desc),
#endif
  KW("dualstrategy",	set_int, XPRS_DUALSTRATEGY,	dualstrategy_desc),
#ifdef XPRS_DUALTHREADS
  KW("dualthreads",	set_int, XPRS_DUALTHREADS,	dualthreads_desc),
#endif
  KW("eigenvaltol",	set_dbl, XPRS_EIGENVALUETOL,	eigenvaltol_desc),
#ifdef XPRS_ELIMFILLIN
  KW("elimfillin",	set_int, XPRS_ELIMFILLIN,	elimfillin_desc),
#endif
  KW("elimtol",		set_dbl, XPRS_ELIMTOL,		elimtol_desc),
  KW("etatol",		set_dbl, XPRS_ETATOL,		etatol_desc),
  KW("feaspump",	set_int, XPRS_FEASIBILITYPUMP,	feaspump_desc),
  KW("feastol",		set_dbl, XPRS_FEASTOL,		"zero tolerance on RHS; default = 1e-6"),
#ifdef XPRS_FEASTOLPERTURB
  KW("feastol_perturb", set_dbl, XPRS_FEASTOLPERTURB, feastol_perturb_desc),
#endif
#ifdef XPRS_FEASTOLTARGET
  KW("feastol_target",	set_dbl, XPRS_FEASTOLTARGET,	feastol_target_desc),
#endif
#ifdef XPRS_MAXGLOBALFILESIZE
  KW("globalfilemax",	set_int, XPRS_MAXGLOBALFILESIZE, globalfilemax_desc),
#endif
#ifdef XPRS_GLOBALFILELOGINTERVAL
  KW("globalloginterval", set_int, XPRS_GLOBALFILELOGINTERVAL, globalloginterval_desc),
#endif
  KW("gomcuts",		set_int, XPRS_GOMCUTS,		gomcuts_desc),
  KW("hdive_rand",	set_dbl, XPRS_HEURDIVERANDOMIZE, hdive_rand_desc),
#ifdef XPRS_HEURDIVESOFTROUNDING
  KW("hdive_rounding",	set_int, XPRS_HEURDIVESOFTROUNDING, hdive_rounding_desc),
#endif
  KW("hdive_speed",	set_int, XPRS_HEURDIVESPEEDUP,	hdive_speed_desc),
  KW("hdive_strategy",	set_int, XPRS_HEURDIVESTRATEGY,	hdive_strategy_desc),
  KW("heurdepth",	set_int, XPRS_HEURDEPTH,	heurdepth_desc),
  KW("heureffort",	set_dbl, XPRS_HEURSEARCHEFFORT,	heureffort_desc),
#ifdef XPRS_HEUREMPHASIS
	  KW("heuremphasis", set_int, XPRS_HEUREMPHASIS, heuremphasis_desc),
#endif
#ifdef XPRS_HEURFORCESPECIALOBJ
  KW("heurforcespecobj", set_int, XPRS_HEURFORCESPECIALOBJ, heurforcespecobj_desc),
#endif
  KW("heurfreq",	set_int, XPRS_HEURFREQ,		heurfreq_desc),
  KW("heurmaxsol",	set_int, XPRS_HEURMAXSOL,	heurmaxsol_desc),
  KW("heurnodes",	set_int, XPRS_HEURNODES,	heurnodes_desc),
  KW("heurroot",	set_int, XPRS_HEURSEARCHROOTSELECT, heurroot_desc),
#ifdef XPRS_HEURSEARCHROOTCUTFREQ
  KW("heurrootcutfreq",	set_int, XPRS_HEURSEARCHROOTCUTFREQ, heurrootcutfreq_desc),
#endif
  KW("heursearch",	set_int, XPRS_HEURSEARCHFREQ,	heursearch_desc),
#ifdef XPRS_HEURSTRATEGY
  KW("heurstrategy",	set_int, XPRS_HEURSTRATEGY,	heurstrategy_desc),
#endif
  KW("heurthreads",	set_int, XPRS_HEURTHREADS,	heurthreads_desc),
  KW("heurtree",	set_int, XPRS_HEURSEARCHTREESELECT, heurtree_desc),
  KW("iis",		set_known, set_iis,		iis_desc),
  KW("indlinbigm",	set_dbl, XPRS_INDLINBIGM,	indlinbigm_desc),
  KW("indprelinbigm",	set_dbl, XPRS_INDPRELINBIGM,	indprelinbigm_desc),
#ifdef XPRS_INPUTTOL
  KW("inputtol",	set_dbl, XPRS_INPUTTOL, inputtol_desc),
#endif
  KW("invertfreq",	set_int, XPRS_INVERTFREQ,	invertfreq_desc),
  KW("invertmin",	set_int, XPRS_INVERTMIN,	invertmin_desc),
  KW("keepbasis",	set_int, XPRS_KEEPBASIS,	keepbasis_desc),
  KW("keepnrows", 	set_int, XPRS_KEEPNROWS,	keepnrows_desc),
  KW("lazy",		I_val, &lazy,			lazy_desc),
  KW("lnpbest",		set_int, XPRS_LNPBEST,		lnpbest_desc),
  KW("lnpiterlimit",	set_int, XPRS_LNPITERLIMIT,	lnpiterlimit_desc),
  KW("localchoice",	set_int, XPRS_LOCALCHOICE,	localchoice_desc),
  KW("logfile",		set_fln, &logfile,		logfile_desc),
#ifdef XPRS_LPFOLDING
  KW("lpfolding",	set_int, XPRS_LPFOLDING,	lpfolding_desc),
#endif
  KW("lpiterlimit",	set_int, XPRS_LPITERLIMIT,	lpiterlimit_desc),
  KW("lplog",		set_int, XPRS_LPLOG,		lplog_desc),
#ifdef XPRS_LPREFINEITERLIMIT
  KW("lpref_itlim",	set_int, XPRS_LPREFINEITERLIMIT, lpref_itlim_desc),
#endif
  KW("lpthreads",	set_int, XPRS_LPTHREADS,	lpthreads_desc),
  KW("markowitztol",	set_dbl, XPRS_MARKOWITZTOL, 	markowitztol_desc),
  KW("matrixtol",	set_dbl, XPRS_MATRIXTOL,	matrixtol_desc),
  KW("maxcuttime",	set_int, XPRS_MAXCUTTIME,	maxcuttime_desc),
  KW("maxiis",		set_int, XPRS_MAXIIS,		maxiis_desc),
  KW("maxim",		set_known, set_maxim,		maxim_desc),
  KW("maximise",	set_known, set_maxim,		maxim_desc),
  KW("maximize",	set_known, set_maxim,		maxim_desc),
#ifdef XPRS_MAXIMPLIEDBOUND
  KW("maximpliedbound",	set_dbl, XPRS_MAXIMPLIEDBOUND,	maximpliedbound_desc),
#endif
  KW("maxlocalbt",	set_int, XPRS_MAXLOCALBACKTRACK, maxlocalbt_desc),
  KW("maxlogscale",	set_int, XPRS_MAXSCALEFACTOR,	maxlogscale_desc),
#ifdef XPRS_MAXMEMORYSOFT
  KW("maxmemory",	set_int, XPRS_MAXMEMORYSOFT,	maxmemory_desc),
#endif
#ifdef XPRS_MAXMEMORYHARD
  KW("maxmemoryhard",	set_int, XPRS_MAXMEMORYHARD,	maxmemoryhard_desc),
#endif
  KW("maxmipsol",	set_int, XPRS_MAXMIPSOL,	maxmipsol_desc),
#ifdef XPRS_MAXMIPTASKS
  KW("maxmiptasks",	set_int, XPRS_MAXMIPTASKS,	maxmiptasks_desc),
#endif
  KW("maxnode",		set_int, XPRS_MAXNODE,		maxnode_desc),
  KW("maxpagelines",	set_int, XPRS_MAXPAGELINES,	maxpagelines_desc),
#ifdef XPRS_MAXSLAVE
  KW("maxslaves",	set_int, XPRS_MAXSLAVE,		maxslaves_desc),
#endif
#ifdef XPRS_MAXSTALLTIME
  KW("maxstalltime", set_int, XPRS_MAXSTALLTIME, maxstalltime_desc),
#endif
  KW("maxtime",		set_int, XPRS_MAXTIME,		maxtime_desc),
  KW("minim",		set_known, set_minim,		minim_desc),
  KW("minimise",	set_known, set_minim,		minim_desc),
  KW("minimize",	set_known, set_minim,		minim_desc),
  KW("mipabscutoff",	set_dbl, XPRS_MIPABSCUTOFF,	mipabscutoff_desc),
  KW("mipabsstop",	set_dbl, XPRS_MIPABSSTOP,	mipabsstop_desc),
  KW("mipaddcutoff",	set_dbl, XPRS_MIPADDCUTOFF,	mipaddcutoff_desc),

#ifdef XPRS_MIPCOMPONENTS
  KW("mipcomponents", set_int, XPRS_MIPCOMPONENTS, mipcomponents_desc),
#endif
#ifdef XPRS_MIPCONCURRENTNODES
  KW("mipconcurnodes", set_int, XPRS_MIPCONCURRENTNODES, mipconcurrentnodes_desc),
  KW("mipconcursolves", set_int, XPRS_MIPCONCURRENTSOLVES, mipconcurrentsolves_desc),
#endif
#ifdef XPRS_MIPDUALREDUCTIONS
  KW("mipdualreductions", set_int, XPRS_MIPDUALREDUCTIONS, mipdualreductions_desc),
#endif
#ifdef XPRS_MIPKAPPAFREQ
  KW("mipkappafreq",	set_int, XPRS_MIPKAPPAFREQ,	mipkappafreq_desc),
#endif
  KW("miplog",		set_int, XPRS_MIPLOG,		miplog_desc),
  KW("mipops",		set_int, XPRS_QSIMPLEXOPS,	mipops_desc),
  KW("mippresolve",	set_int, XPRS_MIPPRESOLVE,	mippresolve_desc),
#ifdef XPRS_MIPRAMPUP
  KW("miprampup",	set_int, XPRS_MIPRAMPUP,	miprampup_desc),
#endif
#ifdef XPRS_MIPREFINEITERLIMIT
  KW("miprefiterlim",	set_int, XPRS_MIPREFINEITERLIMIT, miprefiterlim_desc),
#endif
  KW("miprelcutoff",	set_dbl, XPRS_MIPRELCUTOFF,	miprelcutoff_desc),
  KW("miprelstop",	set_dbl, XPRS_MIPRELSTOP,	miprelstop_desc),
#ifdef XPRS_MIPRESTART
  KW("miprestart", set_int, XPRS_MIPRESTART, miprestart_desc),
#endif
#ifdef XPRS_MIPRESTARTFACTOR
  KW("miprestartfactor", set_dbl, XPRS_MIPRESTARTFACTOR, miprestartfactor_desc),
#endif
#ifdef XPRS_MIPRESTARTGAPTHRESHOLD
  KW("miprestartgaptol", set_dbl, XPRS_MIPRESTARTGAPTHRESHOLD, miprestartgapthreshold_desc),
#endif
  KW("mipstart",	I_val, &mipstart,		"synonym for mipstartvalue"),
  KW("mipstartstatus",	I_val, &mipststat,		mipstartstatus_desc),
  KW("mipstartvalue",	I_val, &mipstart,		mipstartvalue_desc),
#ifdef XPRS_MIPTERMINATIONMETHOD
  KW("mipstop",		set_int, XPRS_MIPTERMINATIONMETHOD, mipstop_desc),
#endif
#ifdef XPRS_MIPTARGET
  KW("miptarget",	set_dbl, XPRS_MIPTARGET,	miptarget_desc),
#endif
  KW("mipthreads",	set_int, XPRS_MIPTHREADS,	mipthreads_desc),
  KW("miptol",		set_dbl, XPRS_MIPTOL,		miptol_desc),
#ifdef XPRS_MIPTOLTARGET
  KW("miptoltarget",	set_dbl, XPRS_MIPTOLTARGET,	miptoltarget_desc),
#endif
#ifdef XPRS_MIQCPALG
  KW("miqcpalg",	set_int, XPRS_MIQCPALG,		miqcpalg_desc),
#endif
#ifdef XPRS_NETSTALLLIMIT
  KW("netstalllimit", set_int, XPRS_NETSTALLLIMIT, netstalllimmit_desc),
#endif
  KW("network",		set_known, set_network,		network_desc),
#ifdef XPRS_GLOBALFILEBIAS
  KW("nodefilebias",	set_dbl, XPRS_GLOBALFILEBIAS,	nodefilebias_desc),
#endif
#ifdef XPRS_NODEPROBINGEFFORT
  KW("nodeprobingeffort", set_dbl, XPRS_NODEPROBINGEFFORT, nodeprobingeffort_desc),
#endif
  KW("nodeselection",	set_int, XPRS_NODESELECTION,	nodeselection_desc),
  KW("objno",		I_val, &nobj,			"objective number (0=none, 1=first...)"),
  KW("objrep",		I_val, &objrep,			objrep_desc),
#ifdef XPRS_OBJSCALEFACTOR
  KW("objscalefactor",	set_dbl, XPRS_OBJSCALEFACTOR,	objscalefactor_desc),
#endif
  KW("optimalitytol",	set_dbl, XPRS_OPTIMALITYTOL,	optimalitytol_desc),
#ifdef XPRS_OPTIMALITYTOLTARGET
  KW("opttol_target", set_dbl, XPRS_OPTIMALITYTOLTARGET, opttol_target_desc),
#endif
  KW("outlev",		I_val,	 &prtmsg,		outlev_desc),
  KW("outputtol",	set_dbl, XPRS_OUTPUTTOL,	outputtol_desc),
  KW("param",		set_par, 0,			param_desc),
  KW("penalty",		set_dbl, XPRS_PENALTY,		penalty_desc),
#ifdef XPRS_PREPERMUTESEED
  KW("permuteseed",	set_int, XPRS_PREPERMUTESEED,	permuteseed_desc),
#endif
#ifdef XPRS_PERTURB
  KW("perturb",		set_dbl, XPRS_PERTURB,		perturb_desc),
#endif
  KW("pivottol",	set_dbl, XPRS_PIVOTTOL,		pivottol_desc),
#ifdef XPRS_MSP_SOLPRB_OBJ /*{*/
  KW("pooldualred",	I_val, &dualred,		pooldualred_desc),
  KW("pooldupcol",	I_val, &dupcol,			pooldupcol_desc),
  KW("pooldups",	set_int, XPRS_MSP_DUPLICATESOLUTIONSPOLICY, pooldups_desc),
  KW("poolfeastol",	set_dbl, XPRS_MSP_SOL_FEASTOL,	poolfeastol_desc),
  KW("poolmiptol",	set_dbl, XPRS_MSP_SOL_MIPTOL,	poolmiptol_desc),
  KW("poolnbest",	I_val, &nbest,			poolnbest_desc),
  KW("poolstub",	set_fln, &poolstub,		poolstub_desc),
#endif /*}*/
  KW("ppfactor",	set_dbl, XPRS_PPFACTOR,		ppfactor_desc),
#ifdef XPRS_PREANALYTICCENTER
  KW("preanalyticcenter", set_int, XPRS_PREANALYTICCENTER, preanalyticcenter_desc),
#endif
#ifdef XPRS_PREBASISRED
  KW("prebasisred",	set_int, XPRS_PREBASISRED,	prebasisred_desc),
#endif
#ifdef XPRS_PREBNDREDCONE
  KW("prebndredcone",	set_int, XPRS_PREBNDREDCONE,	prebndredcone_desc),
#endif
#ifdef XPRS_PREBNDREDQUAD
  KW("prebndredquad",	set_int, XPRS_PREBNDREDQUAD,	prebndredquad_desc),
#endif
  KW("precoefelim",	set_int, XPRS_PRECOEFELIM,	precoefelim_desc),
#ifdef XPRS_PRECOMPONENTS
  KW("precomponents",	set_int, XPRS_PRECOMPONENTS,	precomponents_desc),
#endif
#ifdef XPRS_PRECONVERTSEPARABLE
  KW("preconvertsep",	set_int, XPRS_PRECONVERTSEPARABLE,	preconvertsep_desc),
#endif
  KW("predomcol",	set_int, XPRS_PREDOMCOL,	predomcol_desc),
  KW("predomrow",	set_int, XPRS_PREDOMROW,	predomrow_desc),
#ifdef XPRS_PREDUPROW
  KW("preduprow",	set_int, XPRS_PREDUPROW,	preduprow_desc),
#endif
#ifdef XPRS_PREFOLDING
  KW("prefolding", set_int, XPRS_PREFOLDING, prefolding_desc),
#endif
#ifdef XPRS_PREIMPLICATIONS
  KW("preimplications",	set_int, XPRS_PREIMPLICATIONS,	preimplications_desc),
#endif
#ifdef XPRS_PRELINDEP
  KW("prelindep",	set_int, XPRS_PRELINDEP,	prelindep_desc),
#endif
#ifdef XPRS_PREOBJCUTDETECT
  KW("preobjcutdetect",	set_int, XPRS_PREOBJCUTDETECT,	preobjcutdetect_desc),
#endif
#ifdef XPRS_PREPERMUTE
  KW("prepermute",	set_int, XPRS_PREPERMUTE,	prepermute_desc),
#endif
  KW("preprobing",	set_int, XPRS_PREPROBING,	preprobing_desc),
  KW("presolve",	set_int, XPRS_PRESOLVE,		presolve_desc),
#ifdef XPRS_PRESOLVEMAXGROW
  KW("presolvemaxgrow",	set_dbl, XPRS_PRESOLVEMAXGROW,	presolvemaxgrow_desc),
#endif
  KW("presolveops",	set_int, XPRS_PRESOLVEOPS,	presolveops_desc),
#ifdef XPRS_PRESOLVEPASSES
  KW("presolvepasses",	set_int, XPRS_PRESOLVEPASSES,	presolvepasses_desc),
#endif
  KW("pricingalg",	set_int, XPRS_PRICINGALG,	pricingalg_desc),
  KW("primal",		set_known, set_primal,		primal_desc),
#ifdef XPRS_PRIMALOPS
  KW("primalops",	set_int, XPRS_PRIMALOPS,	primalops_desc),
#endif
#if XPVERSION >= 33
  KW("primalperturb",	set_dbl, XPRS_PRIMALPERTURB,	primalperturb_desc),
#endif
  KW("primalunshift",	set_int, XPRS_PRIMALUNSHIFT,	primalunshift_desc),
  KW("pseudocost",	set_dbl, XPRS_PSEUDOCOST,	pseudocost_desc),
  KW("pseudocost_ud",	set_int, XPRS_HISTORYCOSTS,	pseudocost_ud_desc),
#ifdef XPRS_QCCUTS
  KW("qccuts",		set_int, XPRS_QCCUTS,		qccuts_desc),
  KW("qcrootalg",	set_int, XPRS_QCROOTALG,	qcrootalg_desc),
#endif
  KW("quadunshift",	set_int, XPRS_QUADRATICUNSHIFT,	quadunshift_desc),
  KW("ray",		I_val, &Ray,			ray_desc),
#ifdef XPRS_REFINEOPS
  KW("refineops",	set_int, XPRS_REFINEOPS,	refineops_desc),
#endif
  KW("relax",		set_known, set_relax,		"[no assignment] ignore integrality"),
  KW("relaxtreemem",	set_dbl, XPRS_RELAXTREEMEMORYLIMIT, relaxtreemem_desc),
  KW("relpivottol",	set_dbl, XPRS_RELPIVOTTOL,	relpivottol_desc),
  KW("repairindefq",	set_int, XPRS_REPAIRINDEFINITEQ, repairindefq_desc),
#ifdef XPRS_RESOURCESTRATEGY
  KW("resourcestrategy", set_int, XPRS_RESOURCESTRATEGY, resourcestrategy_desc),
#endif
  KW("rootpresolve",	set_int, XPRS_ROOTPRESOLVE,	rootpresolve_desc),
  KW("round",		I_val, &Round,			round_desc),
  KW("sbbest",		set_int, XPRS_SBBEST,		sbbest_desc),
  KW("sbeffort",	set_dbl, XPRS_SBEFFORT,		sbeffort_desc),
  KW("sbestimate",	set_int, XPRS_SBESTIMATE,	sbestimate_desc),
  KW("sbiterlimit",	set_int, XPRS_SBITERLIMIT,	sbiterlimit_desc),
  KW("sbselect",	set_int, XPRS_SBSELECT,		sbselect_desc),
  KW("scaling",		set_int, XPRS_SCALING,		scaling_desc),
#ifdef XPRS_SIFTING
  KW("sifting",		set_int, XPRS_SIFTING,		sifting_desc),
#endif
#ifdef XPRS_SIFTPASSES
  KW("siftpasses", set_int, XPRS_SIFTPASSES, siftpasses_desc),
#endif
#ifdef XPRS_SIFTPRESOLVEOPS
  KW("siftpresolveops", set_int, XPRS_SIFTPRESOLVEOPS, siftpresolveops_desc),
#endif
#ifdef XPRS_SIFTSWITCH
  KW("siftswitch", set_int, XPRS_SIFTSWITCH, siftswitch_desc),
#endif
#ifdef XPRS_SLEEPONTHREADWAIT
  KW("sleeponthreadwait", set_int, XPRS_SLEEPONTHREADWAIT, sleeponthreadwait_desc),
#endif
  KW("sos",		I_val, &sos,			sos_desc),
  KW("sos2",		I_val, &sos2,			sos2_desc),
  KW("sosreftol",	set_dbl, XPRS_SOSREFTOL,	sosreftol_desc),
#ifdef XPRS_SYMMETRY
  KW("symmetry",	set_int, XPRS_SYMMETRY,		symmetry_desc),
#endif
#ifdef XPRS_TEMPBOUNDS
  KW("tempbounds",	set_int, XPRS_TEMPBOUNDS,	tempbounds_desc),
#endif
  KW("threads",		set_int, XPRS_THREADS,		threads_desc),
  KW("timing",		set_known, set_timing,		"[no assignment] give timing statistics"),
  KW("trace",		set_int, XPRS_TRACE,		trace_desc),
  KW("treecompress",	set_int, XPRS_TREECOMPRESSION,	treecompress_desc),
  KW("treecovercuts",	set_int, XPRS_TREECOVERCUTS,	treecovercuts_desc),
  KW("treecuts",	set_int, XPRS_TREECUTSELECT,	treecuts_desc),
  KW("treegomcuts",	set_int, XPRS_TREEGOMCUTS,	treegomcuts_desc),
  KW("treememlimit",	set_int, XPRS_TREEMEMORYLIMIT,	treememlimit_desc),
  KW("treememtarget",	set_dbl, XPRS_TREEMEMORYSAVINGTARGET, treememtarget_desc),
  KW("treeoutlev",	set_int, XPRS_TREEDIAGNOSTICS,	treeoutlev_desc),
#ifdef XPRS_TREEPRESOLVE
  KW("treepresolve",	set_int, XPRS_TREEPRESOLVE,	treepresolve_desc),
#endif
#if XPVERSION >= 30 /*{*/
  KW("tunerdir",	set_fln, &tunerdir,		tunerdir_desc),
  KW("tunerhistory",	set_int, XPRS_TUNERHISTORY,	tunerhistory_desc),
  KW("tunermaxtime",	set_int, XPRS_TUNERMAXTIME,	tunermaxtime_desc),
  KW("tunermethod",	set_int, XPRS_TUNERMETHOD,	tunermethod_desc),
  KW("tunermethodfile",	set_fln, &tunermethodfile,	tunermethodfile_desc),
  KW("tunerpermute",	set_int, XPRS_TUNERPERMUTE,	tunerpermute_desc),
  KW("tunertarget",	set_int, XPRS_TUNERTARGET,	tunertarget_desc),
  KW("tunerthreads",	set_int, XPRS_TUNERTHREADS,	tunerthreads_desc),
#endif /*}*/
  KW("varselection",	set_int, XPRS_VARSELECTION,	varselection_desc),
  KW("version",		Ver_val, 0,			version_desc),
  KW("wantsol",		WS_val, 0,			WS_desc),
  KW("writeprob",	set_fln, &writeprob,		writeprob_desc)
     };

static Option_Info Oinfo = { "xpress", NULL, "xpress_options",
		keywds,nkeywds,0,"XPRESS", 0,0,0,0,0, 20220112, 0,0,0,0,0,0,0,
		ASL_OI_tabexpand | ASL_OI_addnewline };

 static char *
strcpy1(char *t, const char *s, char *te)
{
	while(t < te && (*t = *s++))
		t++;
	return t;
	}

 static void
adjust_version(void)
{
	char *t, *te;
	static char vbuf[64];
#ifdef XPRS_XPRESSVERSION
	char rbuf[16];
	int vlen;

	vlen = 0;
	if (XPRSgetstringattrib(prob, XPRS_XPRESSVERSION, rbuf, (int)sizeof(rbuf), &vlen))
		*rbuf = 0;
#endif
	te = vbuf + sizeof(vbuf) - 1;
	*te = 0;
	t = strcpy1(Oinfo.bsname = vbuf, "XPRESS ", te);
#if defined(EXPRESS_RELEASE)
	t = strcpy1(t, EXPRESS_RELEASE, te);
#elif defined(XPRS_XPRESSVERSION)
	t = strcpy1(t, rbuf, te);
#endif
	t = strcpy1(t, "(", te);
	if (!XPRSgetversion(t))
		t += strlen(t);
	*t++ = ')';
	Oinfo.version = ++t;
	t = strcpy1(t, "AMPL/", te);
	t = strcpy1(t, Oinfo.bsname, te); /* "t =..." so we can check t-vbuf in a debugger */
	}

 static int breaking;
 static jmp_buf Jb;

 static void
intcatch(int n)
{
	printf("\n<BREAK> (xpress)\n");
	fflush(stdout);
	if (!breaking++)
		XPRSinterrupt(prob, XPRS_STOP_CTRLC);
	else if (breaking > 3)
		longjmp(Jb, 2);
	}

/* Xpress-MP callback in case the user wants some output from Optimizer */
void XPRS_CC xpdisplay(XPRSprob prob, void *data, const char *ch, int n, int msglvl)
{
  /*
   msglvl gives the message level as follows:
   * 1 dialogue
   * 2 information
   * 3 warnings
   * 4 errors
   * a negative value indicates the XPRESS is about to finish and
   * buffers should be flushed.
   */

  /* You could divert the messages to your own log file if you wanted */

  if (msglvl < 0)
    fflush(NULL);
  else if (msglvl >= prtmsg && (msglvl != 4 || strncmp(ch, "?899 ", 5)))
    printf("%s\n", ch);
}

/***********************/
/* Print abort message */
/***********************/
static void xperror(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  fprintf(stderr, progname ? "%s: Error " : "Error ", progname);
  vfprintf(stderr, fmt, ap);
  fprintf(stderr, ".\n");
  va_end(ap);
  exit(1);
}

#ifndef LICERROR
 static void
licerror(int iret)
{
	char buf[512];

	buf[0] = 0;
	XPRSgetlicerrmsg(buf, sizeof(buf));
	xperror("initializing Xpress-MP (return code %d):\n%s\n",
		iret, buf);
	}
#endif

/******************************/
/* Delete .sol and .glb files */
/******************************/
static void killtempprob(void)
{
  int len = strlen(probname);
#if 0 /* would be needed if a call on XPRSwritebinsol were added */
  strcat(probname,".sol");
  remove(probname);
  probname[len]='\0';
#endif
  strcat(probname,".glb");
  remove(probname);
  probname[len]='\0';
}

#if XPVERSION >= 21
 static char iis_table[] = "\n\
0	non	not in the iis\n\
1	low	at lower bound\n\
2	fix	fixed\n\
3	upp	at upper bound\n\
4	mem	member\n\
5	pmem	possible member\n\
6	plow	possibly at lower bound\n\
7	pupp	possibly at upper bound\n\
8	bug\n";
#endif

static SufDecl
suftab[] = {
  { "bestbound", 0, ASL_Sufkind_obj | ASL_Sufkind_outonly},
  { "bestbound", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly},
  { "direction", 0, ASL_Sufkind_var },
#if XPVERSION >= 21
  { "iis", iis_table, ASL_Sufkind_var | ASL_Sufkind_outonly },
  { "iis", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
  { "iso", 0, ASL_Sufkind_var | ASL_Sufkind_outonly },
  { "iso", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
#endif
  { "lazy", 0, ASL_Sufkind_con },
#ifdef XPRS_MSP_SOLPRB_OBJ
  { "npool", 0, ASL_Sufkind_obj | ASL_Sufkind_outonly},
  { "npool", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly},
#endif
  { "priority", 0, ASL_Sufkind_var },
  { "ref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
  { "sos", 0, ASL_Sufkind_var },
  { "sos", 0, ASL_Sufkind_con },
  { "sosno", 0, ASL_Sufkind_var | ASL_Sufkind_real },
  { "sosref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
  { "sstatus", 0, ASL_Sufkind_var, 1 },
  { "sstatus", 0, ASL_Sufkind_con, 1 },
  { "unbdd", 0, ASL_Sufkind_var | ASL_Sufkind_real },
  };


#if XPVERSION >= 30
 static int
tune(void) {
#if !defined(S_IFREG) && defined(_S_IFREG) /* MS Windows */
#define S_IFREG _S_IFREG
#define S_IFDIR _S_IFDIR
#endif
#ifdef _WIN32
#define Mkdir(x) _mkdir(x)
#else
#define Mkdir(x) mkdir(x, 0755)
#endif
	char *b, errbuf[512], *s;
	const char *fmt, *name;
	int c, j, k;
	size_t L;
	struct stat sb;

	b = 0;
	for(s = tunerdir; ; ++s) {
		if (!(c = *s)
		|| c == '/'
#ifdef _WIN32
		|| c == '\\' || c == ':'
#endif
			) {
			*s = 0;
			if (stat(tunerdir, &sb)) {
				if (Mkdir(tunerdir))
					xperror("Failed to create directory \"%s\"", tunerdir);
				}
			else if (!(sb.st_mode & S_IFDIR)) {
				fmt = "Invalid tunerdir=\"%s\".";
				solve_result_num = 530;
 badtune:
				L = strlen(tunerdir) + 40;
				s = (char*)M1alloc(L);
				snprintf(s, L, fmt, tunerdir);
 badret:
				write_sol(s, 0, 0, &Oinfo);
				return 1;
				}
			if (!(*s = c) || !s[1])
				break;
			b = s;
			}
		}
	if (b) {
		c = *b;
		*b = 0;
		name = tunerdir;
		}
	else
		name = ".";
	if (XPRSsetstrcontrol(prob, XPRS_TUNEROUTPUTPATH, name))
		xperror("Setting tuneroutputpath to \"%s\" failed", name);
	if (b)
		*b++ = c;
	else
		b = tunerdir;
	if (XPRSsetstrcontrol(prob, XPRS_TUNERSESSIONNAME, b))
		xperror("Setting tunersessionname to \"%s\" failed", b);
	if (tunermethodfile) {
		if (stat(name = tunermethodfile, &sb)
		|| !(sb.st_mode & S_IFREG)
		|| XPRSsetstrcontrol(prob, XPRS_TUNERMETHODFILE, tunermethodfile)) {
			fmt = "Bad tunermethodfile=\"%s\".";
			solve_result_num = 532;
			goto badtune;
			}
		}
	if ((k = XPRStune(prob, optimopt))) {
		memset(errbuf, 0, sizeof(errbuf));
		XPRSgetlasterror(prob, errbuf);
		for(j = 0; j < sizeof(errbuf) && errbuf[j]; ++j);
		s = (char*)M1alloc(L = j + 80);
		snprintf(s,L,"Surprise return %d from XPRStune(): \"%.*s\"", k, j, errbuf);
		solve_result_num = 531;
		goto badret;
		}
	return 0;
	}
#endif

/**********************/
/* The main procedure */
/**********************/

int main(int argc, char *argv[])
{
 char *stub;
 int iret;
 dims d;
 void (*oic)(int);

 Times[0] = xectim_();

 /* for debugging, allow -=, -v, -? to work */
 if (!(iret = XPRSinit(XPRESS)) && (iret = XPRScreateprob(&prob)))
	xperror("Creating Xpress-MP problem(return code %d)\n%s", iret, "Have you set XPRESS?");
 adjust_version();
 if (argc == 2 && *(stub = argv[1]) == '-' && stub[1] && !stub[2]) {
	asl = (ASL*)ASL_alloc(ASL_read_fg);
	getstub(&argv, &Oinfo);
	return 0;
	}
 if (iret)
   licerror(iret);

#ifdef UNIX
  XPRSsetlogfile(prob,"/dev/null");  /* Kill output from XPRSinit() */
#endif

 XPRSsetcbmessage(prob, xpdisplay, NULL);

 asl = (ASL*)ASL_alloc(ASL_read_fg);		/* Allocate a structure */
 if(!(stub = getstub(&argv, &Oinfo)))		/* Get the 'stub' name */
   usage_ASL(&Oinfo,1);

 if(tmpnam(probname)==NULL) xperror("tmpnam() failed: cannot obtain a problem name");
                                                 /* Get temporary problem name */

 suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));

 amplin(stub,argv,&d);        /* Read and load the problem */

#if XPVERSION >= 31
 if (archconsistent)
	XPRS_ge_setarchconsistency(1);
#endif

#ifdef XPRS_SOLUTIONFILE
 if(optimopt[1]!='g')
   XPRSsetintcontrol(prob,XPRS_SOLUTIONFILE,0);  /* Don't save solution */
#endif


 if(startbasis!=NULL)                            /* Load an initial basis */
   if(XPRSreadbasis(prob,startbasis,""))
     xperror("loading an initial basis");

 Times[1] = xectim_();

 breaking = 0;
 iret = setjmp(Jb);
 oic = signal(SIGINT, intcatch);
 if (iret) {
	signal(SIGINT, oic);
	--iret;
	Times[2] = xectim_();
	if (amplflag | (Oinfo.wantsol & 1)) {
		if (solve_result_num <= 0)
			solve_result_num = 520;
		iret = 0;
		write_sol("Break!", 0, 0, &Oinfo);
		}
	else
		printf("Break received.\n");
	}
#if XPVERSION >= 30
 else if (tunerdir && tune())
	; /* write_sol already called */
#endif
 else {
	iret = Optimise(prob,optimopt);                /* Optimise the problem */
	if (iret < 0 || (iret > 8 && iret != 32))
		xperror("optimising the problem:  surprise return %d", iret);
	Times[2] = xectim_();

	if((endbasis!=NULL)&&(optimopt[1]!='g')        /* Save the final basis */
	  && XPRSwritebasis(prob,endbasis,""))
		xperror("saving the current basis");

	amplout(&d);                                    /* Export the solution */
	}
 show_times();
 XPRSfree();
 return iret;
}

/*************************/
/* Set a known parameter */
/*************************/
 static char *
set_known(Option_Info *oi, keyword *kw, char *v)
{
#ifdef RWA_DEBUG
 char *d;
#endif
 switch((size_t) kw->info)
 {
  case set_primal:   optimopt[0]='p';  break;
#ifdef RWA_DEBUG
  case set_debug:  d=debugopt;
                        while( *v != '\0' && *v != ' ' ) *d++ = *v++;
      *d++ = '\0';
                        break;
#endif
  case set_dual:	optimopt[0] = 'd';	break;
  case set_barrier:	optimopt[0] = 'b';	break;
  case set_relax:	optimopt[1] = 'l';	break;
  case set_maxim:	Optimise = XPRSmaxim;	break;
  case set_minim:	Optimise = XPRSminim;	break;
  case set_timing:	timing = 1;		break;
  case set_iis:		iis_find = 1;		break;
  case set_network:	optimopt[2] = 'n';	break;
  case set_bestbound:	bestbound = 1;		break;
  default:          printf("Unknown option %s\n",kw->name);
                    badopt_ASL(oi);
 }
 return v;
}

/**************************************/
/* Set an Xpress-MP integer parameter */
/**************************************/

 static int
IntDset(Defer_setting *ds)
{
	if (XPRSsetintcontrol(prob,ds->ipar,ds->u.i)) {
		printf("The value %d could not be assigned belatedly to %s\n",
			ds->u.i, ds->kw->name);
		return 1;
		}
	return 0;
	}

static char *set_int(Option_Info *oi, keyword *kw, char *v)
{
  Defer_setting *ds;
  int ipar,isval;
  char *rv;

  ipar = (int) (unsigned long) kw->info;

  if ((*v=='?') && (v[1]<=' ')) {
    if (ipar < 0)
	ipar = -ipar;
    XPRSgetintcontrol(prob,ipar,&isval);
    printf("%s=%d\n",kw->name,isval);
    oi->option_echo &= ~ASL_OI_echothis;
    return v+1;
  }
  isval = (int)strtol(v,&rv,10);
  if (v==rv) {
    printf("Expected a numeric value for %s, not \"%s\"\n",kw->name,v);
    badopt_ASL(oi);
  } else if (ipar < 0) {
	ds = new_Defer_setting(IntDset, kw, ipar);
	ds->u.i = isval;
  } else if (XPRSsetintcontrol(prob,ipar,isval)) {
    printf("The value %d is not allowed for %s\n",isval,kw->name);
    badopt_ASL(oi);
  }

  return rv;
}

/*************************************/
/* Set an Xpress-MP double parameter */
/*************************************/

 static int
DblDset(Defer_setting *ds)
{
	if (XPRSsetdblcontrol(prob,ds->ipar,ds->u.d)) {
		printf("The value %g could not be assigned belatedly to %s\n",
			ds->u.d, ds->kw->name);
		return 1;
		}
	return 0;
	}

static char *set_dbl(Option_Info *oi, keyword *kw, char *v)
{
 Defer_setting *ds;
 int ipar;
 double dgval;
 char *rv;

 ipar= (int) (unsigned long) kw->info;
 if((*v=='?') && (v[1]<=' '))
 {
  if (ipar < 0)
	ipar = -ipar;
  XPRSgetdblcontrol(prob,ipar,&dgval);
  printf("%s=%g\n",kw->name,dgval);
  oi->option_echo &= ~ASL_OI_echothis;
  return v+1;
 }
 dgval= strtod(v,&rv);
 if(v==rv)
 {
  printf("Expected a numeric value for %s, not \"%s\"\n",kw->name,v);
  badopt_ASL(oi);
 }
 else if (ipar < 0) {
	ds = new_Defer_setting(DblDset, kw, ipar);
	ds->u.d = dgval;
	}
 else
  if(XPRSsetdblcontrol(prob,ipar,dgval))
  {
   printf("The value %g is not allowed for %s\n",dgval,kw->name);
   badopt_ASL(oi);
  }
 return rv;
}

/***********************************/
/* Set a filename option/parameter */
/***********************************/
static char *set_fln(Option_Info *oi, keyword *kw, char *v)
{
 char *rv, *t,**f;
 int n, q;

 if(!(q = *v))
 {
  printf("rejecting %s: no following file name\n", kw->name);
  badopt_ASL(oi);
  return v;
 }
 f = (char**)kw->info;
 if (q == '?' && v[1] <= ' ')
 {
  if (!(t = *f))
	t = "";
  printf("%s=\"%s\"\n", kw->name, t);
  oi->option_echo &= ~ASL_OI_echothis;
  return v + 1;
 }
 /* Allow optional quoting of file name with ' or " */
 /* Quoted file names may contain  blanks. */
 if (q == '"' || q == '\'')
	for(rv = ++v; *rv != q && *rv; ++rv);
 else
	for(rv = v; *(unsigned char*)++rv > ' ';);
 if (!(n = rv - v))
	t = 0;
 else {
	t = M1alloc(n + 1);
	strncpy(t, v, n);
	t[n] = '\0';
	}
 *f = t;
 if (q && *rv == q)
	rv++;
 return rv;
}

/*********************************/
/* Exit if problem is non linear */
/*********************************/

 static int
nlify(int k, const char *msg)
{
#ifndef NLIFY_LL
#define NLIFY_LL 72
#endif
	/* break error message into lines <= NLIFY_LL characters long */
	const char *s;
	int j;

	for(;;) {
		while(*msg <= ' ')
			if (!*msg++)
				return k;
		for(s = msg; *++s > ' ';);
		j = (int)(s - msg);
		if (j + k > NLIFY_LL)
			k = fprintf(Stderr, "\n%.*s", j, msg) - 1;
		else
			k += fprintf(Stderr, " %.*s", j, msg);
		msg += j;
		}
	}

 static void
nonlin(int n, char *what)
{
	if (n) {
		n = fprintf(Stderr, "Sorry, %s", Oinfo.bsname);
		n = nlify(n, "cannot handle");
		nlify(n, what);
		putc('\n', Stderr);
		exit(4);
		}
	}

/*****************************************/
/* Read priorities for integral entities */
/*****************************************/

 static int
nzeros(int *p, int n)
{
  int *pe = p + n;
  n = 0;
  while(p < pe)
    if (*p++)
      n++;
  return n;
}


 static int
priority_suf(void)
{
  int baddir, badpri, i, j, k, listsize, nnames;
  int *colindex, *dir, *p, *priority;
  char *direction;
  SufDesc *dp, *dd;

  baddir = badpri = i = listsize = 0;
  dd = suf_get("direction", ASL_Sufkind_var);
  dir = dd->u.i;
  dp = suf_get("priority", ASL_Sufkind_var);
  if (!(p = dp->u.i) && !dir)
    return 0;
  direction = 0;
  nnames = n_var;
  k = p ? 2 : 1;
  if (p && dir) {
    for(; i < nnames; i++)
      if (p[i] || dir[i])
        listsize++;
    }
  else if (p)
    listsize = nzeros(p, nnames);
  else
    listsize = nzeros(dir, nnames);
  if (!listsize)
    return 1;
  priority = 0;
  colindex = (int*)M1alloc(k*listsize*sizeof(int));
  if (dir)
    direction = (char*)M1alloc(listsize);
  i = k = 0;
  if (dir && p) {
    priority = colindex + listsize;
    for(; i < nnames; i++)
      if (p[i] || dir[i]) {
        colindex[k] = i;
        if ((j = p[i]) < 0) {
          badpri++;
          j = 0;
          }
        priority[k] = j;
        switch(dir[i]) {
          case -1:
          j = 'D';
          break;
          case 1:
          j= 'U';
          break;
          default:
          baddir++;
          /* no break */
          case 0:
          j = 'N';
          }
        direction[k++] = j;
        }
    }
  else if (p) {
    priority = colindex + listsize;
    for(; i < nnames; i++)
      if (p[i]) {
        colindex[k] = i;
        if ((j = p[i]) < 0) {
          badpri++;
          j = 0;
          }
        priority[k++] = j;
        }
    }
  else {
    for(; i < nnames; i++)
      if (dir[i]) {
        switch(dir[i]) {
          case -1:  j = 'D';      break;
          case 1:   j = 'U';      break;
          default:  baddir++;   /* no break */
          case 0:   j = 'N';
          }
        direction[k] = j;
        colindex[k++] = i;
        }
    }
  listsize = k;
  if (baddir)
    fprintf(Stderr,
   "Treating %d .direction values outside [-1, 1] as 0.\n",
      baddir);
  if (badpri)
    fprintf(Stderr,
      "Treating %d negative .priority values as 0\n",
      badpri);
  if (XPRSloaddirs(prob,listsize,colindex,priority,direction,NULL,NULL))
    xperror("loading priorities");
  return 1;
  }

 static void
 mip_priorities(void)
{
 int nbpri, *start, *pri, *num;
 int ndir, *mcols, *mpri;
 int p, c, curr;

 if (priority_suf())
  return;
 if((nbpri=mip_pri(&start, &num, &pri, 2147483647))>0)
 {
  ndir=0;
  for(p=0;p<nbpri;p++) ndir+=num[p];  /* Number of priorities */
  mcols=(int *)M1alloc(ndir*sizeof(int));
  mpri=(int *)M1alloc(ndir*sizeof(int));

  curr=0;
  for(p=0;p<nbpri;p++)        /* Create mcols & mpri... */
   for(c=0;c<num[p];c++)
   {
    mcols[curr]=start[p]+c;
    mpri[curr]=pri[p];
    curr++;
   }
  if(XPRSloaddirs(prob,ndir,mcols,mpri,NULL,NULL,NULL))
    xperror("loading priorities");
 }
}

 static void
stat_map(int *stat, int n, int *map, int mx, char *what)
{
  int bad, i, i1=0, j, j1=0;
  static char badfmt[] = "XPRESS driver: %s[%d] = %d\n";

  for(i = bad = 0; i < n; i++) {
    if ((j = stat[i]) >= 0 && j <= mx)
      stat[i] = map[j];
    else {
      stat[i] = 0;
      i1 = i;
      j1 = j;
      if (!bad++)
        fprintf(Stderr, badfmt, what, i, j);
      }
    }
  if (bad > 1) {
    if (bad == 2)
      fprintf(Stderr, badfmt, what, i1, j1);
    else
      fprintf(Stderr,
    "Xpress-MP driver: %d messages about bad %s values suppressed.\n",
        bad-1, what);
    }
  }

 static void
get_statuses(dims *d)
{
  static int map[] = {0, 1, 3, 0, 2, 0, 2};

  if (!advance)
	return;
  if ((!mipststat && niv + nbv)
   || (!(d->csd->kind & ASL_Sufkind_input)
   && !(d->rsd->kind & ASL_Sufkind_input)))
    return;
  stat_map(d->cstat, n_var, map, 7, "incoming cstat");
  stat_map(d->rstat, n_con, map, 7, "incoming rstat");
  if (XPRSloadbasis(prob,d->rstat, d->cstat))
    xperror("loading statuses");
  }

 static void
qcadj(CStype **pnelq, int **pqcrows, int **qcol1, int **qcol2, real **pqv, void **v)
{
	/* Adjust for quadratic constraints. */

	QPinfo *qpi;
	CStype j, k, *ka, *ka1, *nelq, nnz;
	cgrad **cgp, *cg, *cga, *cg1, **cgt;
	char errbuf[256];
	double *a, *a1, *qv, *x, *y;
	int i, i1, i2, j1, je, kl, nl, m, n, nc, nqc, nqv;
	int *cn, *col0, *col1, *col2, *ia, *ia1, *qcrows, *rowq;
	size_t L, *cq, js, ks, nqcnl;
	ssize_t nq;

	n = n_var;
	nqc = nlc;
	nqv = nlvc;
	ia = A_rownos;
	ka = A_colstarts;
	a = A_vals;

	j = nzc;
	for(i = kl = nl = 0; i < j; i++) {
		if (ia[i] < nqc)
			++nl;
		else
			++kl;
		}

	cg = cga = (cgrad*)Malloc((nqc+n)*sizeof(cgrad*) + nl*sizeof(cgrad));
	Cgrad = cgp = (cgrad**)(cga + nl);
	cgt = cgp + nqc;
	memset(cgp, 0, (nqc+n)*sizeof(cgrad*));

	j = 0;
	for(i = 0; i < n; i++) {
		for(k = ka[i+1]; j < k; j++)
			if ((i1 = ia[j]) < nqc) {
				cg->varno = i;
				cg->coef = a[j];
				cg->next = cgp[i1];
				cgp[i1] = cg++;
				}
		}
	x = LUrhs;
	for(nqcnl = i = 0; i < nqc; i++) {
		i2 = i << 1;
		if (x[i2] > negInfinity && x[i2+1] < Infinity) {
			snprintf(errbuf, sizeof(errbuf),
			 "constraint %s, which is not convex quadratic since it is %s constraint.",
				con_name(i), x[i2] == x[i2+1] ? "an equality"
						: "a two-sided");
			nonlin(1, errbuf);
			}
		nq = mqpcheckv(-(i+1), 0, v);
		if (nq < 0) {
			nonlin(nq == -2,
			 "a quadratic constraint involving division by 0");
			nonlin(1, "a nonquadratic nonlinear constraint");
			}
		nonlin(nq == 0,
		 "a driver bug: no quadratic terms in a \"nonlinear\" constraint");
		nqcnl += nq;
		}

	x = *pqv = qv = (double*)Malloc(nqc*sizeof(CStype)
					+ nqcnl*(2*sizeof(int) + sizeof(double)));
	*pnelq = nelq = (CStype*)(qv + nqcnl);
	*qcol1 = col1 = (int*)(nelq + nqc);
	*qcol2 = col2 = col1 + nqcnl;

	for(i = m = nl = 0; i < nqc; i++) {
		nq = mqpcheckv(-(i+1), &qpi, v);
		rowq = qpi->rowno;
		cq = qpi->colbeg;
		cn = qpi->colno;
		nc = qpi->nc;
		y = qpi->delsq;
		col0 = col1;
		/* Discard elements below the diagonal, and account for XPRESS's */
		/* expected scaling of quadratic constraints. */
		js = cq[0];
		for(j = 0; j < nc; ) {
			k = cn[j++];
			for(ks = cq[j]; js < ks; ++js) {
				if (rowq[js] <= k) {
					*x++ = 0.5 * y[js];
					*col1++ = rowq[js];
					*col2++ = k;
					}
				}
			}
		nelq[i] = col1 - col0;
		free(qpi);
		}

	nl = 0;
	while(i > 0) {
		for(cg = cgp[--i]; cg; cg = cg1) {
			cg1 = cg->next;
			if (cg->coef != 0.) {
				++nl;
				cg->next = cgt[j = cg->varno];
				cgt[j] = cg;
				cg->varno = i;
				}
			}
		}

	nnz = kl + nl + ka[n] - ka[nqv] + 1;
	L = nnz*sizeof(real) + (nqc + nnz)*sizeof(int) + (n + 2)*sizeof(CStype);
	A_vals = a1 = (real*)Malloc(L);
	A_colstarts = ka1 = (CStype*)(a1 + nnz);
	A_rownos = ia1 = (int*)(ka1 + n + 2);
	*pqcrows = qcrows = (int*)(ia1 + nnz);
	*ka1 = 0;
	for(i = i1 = j = j1 = 0; i < n; ++i) {
		for(cg = cgt[i]; cg; cg = cg->next) {
			ia1[j1] = cg->varno;
			a1[j1++] = cg->coef;
			}
		for(je = ka[i+1]; j < je; ++j) {
			if (ia[j] >= nqc) {
				ia1[j1] = ia[j];
				a1[j1++] = a[j];
				}
			}
		*++ka1 = j1;
		}
	for(i = 0; i < nqc; ++i)
		qcrows[i] = i;
	free(cga);
	free(a);
	}

 typedef struct
IndicInfo {
	int nic, nr;
	int *rn;
	int *indv;
	int *comps;
	} IndicInfo;

 static int
add_indic(void *v, int iv, int compl, int sense, int nz, int *ig, real *g, real rhs)
{
	IndicInfo *II = (IndicInfo*)v;
	int i, j, k, mstart[2];

	mstart[0] = 0;
	mstart[1] = nz;
	if (sense <= 1) {
		i = XPRSaddrows(prob, 1, nz, "LGE" + sense, &rhs, NULL, mstart, ig, g);
		k = 1;
		}
	else {
		/* Since XPRESS disallows equalities here, we supply two inequalities. */
		i = XPRSaddrows(prob, 1, nz, "L", &rhs, NULL, mstart, ig, g)
		    ||
		    XPRSaddrows(prob, 1, nz, "G", &rhs, NULL, mstart, ig, g);
		k = 2;
		}
	if (i)
		return i;
	for(j = 0; j < k; ++j) {
		i = II->nic++;
		II->rn[i] = II->nr++;
		II->indv[i] = iv;
		II->comps[i] = compl ? -1 : 1;
		}
	return 0;
	}

 static void
indicator_constrs(void)
{
	IndicInfo II;
	int errinfo[2], i, nlogc;

	nlogc = 4*n_lcon;	/* Allow space for "else" constraints and two */
				/* inequalities instead of an equality. */
	II.nic = 0;
	if (XPRSgetintattrib(prob, XPRS_ROWS, &II.nr))
		xperror("XPRSgetintattrib(XPRS_ROWS)");
	II.rn = (int*)Malloc(3*nlogc*sizeof(int));
	II.indv = II.rn + nlogc;
	II.comps  = II.indv + nlogc;
	if ((i = indicator_constrs_ASL(asl, &II, add_indic, errinfo))) {
		switch(i) {
		  case 1:
			xperror("logical constraint %s is not an indicator constraint.\n",
				lcon_name(errinfo[0]));
			break;
		  case 2:
			xperror("logical constraint %s is not an indicator constraint\n\
	due to bad comparison with %s.\n", lcon_name(errinfo[0]), var_name(errinfo[1]));
			break;
		  case 3:
			xperror("adding indicator row");
			break;
		  default:
			xperror("indicator_constrs_ASL");
		  }
		}
	if (XPRSsetindicators(prob, II.nic, II.rn, II.indv, II.comps))
		xperror("XPRSsetindicators");
	free(II.rn);
	}

 static void
lazy_adj(int *z)
{
	const char *s;
	int i, j, nc, nlin, nq, nqc;

	nc = n_con;
	nqc = nlc;
	for(i = nq = 0; i < nqc; ++i)
		if (z[i])
			++nq;
	if (nq) {
		s = nq == 1 ? "" : "s";
		fprintf(Stderr, "Ignoring .lazy suffix%s on %d quadratic constraint%s.\n",
			s, nq, s);
		}
	for(nlin = 0; i < nc; ++i)
		if (z[i])
			++nlin;
	if (!nlin)
		return;
	j = 0;
	for(i = nqc; i < nc; ++i)
		if (z[i])
			z[j++] = i;
	if (XPRSloaddelayedrows(prob, j, z))
		xperror("load delayed rows");
	}

/************************************************/
/* Read the matrix and call loadprob/loadglobal */
/************************************************/
 static void
amplin(char *stub, char *argv[], dims *d)
{
 CStype *ka, *msstart, nz, *qmn;
 FILE *nl;
 QPinfo *qpi;
 SufDesc *lzd;
 char *qgtype, *qrtype, *qstype, *s;
 double *L, *U, *a, *dref, *obj, *q, *q0, *q1, *qe, *qv, *rhs, *rng;
 int *cn, *ia, *mgcols, *mqc1=NULL, *mqc2=NULL, *mscols;
 int *ims0, *qcol1, *qcol2, *qcrows, *rq, *rq1, *vmi;
 int i, iret, j, k, m, m1, n, n0, n_bv, ncol, ngents, nqc, nsets, row;
 ograd *og;
 size_t *colq, nelq, nq;
 ssize_t nelqs;
 struct LU_bounds *rhs_bounds;
 void *v;
#ifndef NO_XPR_CON_INDICATOR
#define ALLOW_CLP ASL_allow_CLP
 int nlogc;
#else
#define ALLOW_CLP 0
#define nlogc 0
#endif
 static int repmap[4] = { 0, ASL_obj_replace_ineq,  ASL_obj_replace_eq,
			     ASL_obj_replace_ineq | ASL_obj_replace_eq };

 nl = jac0dim(stub, (fint)strlen(stub));

 /* allow elbow room for objadj */
 m = n_con;
#ifndef NO_XPR_CON_INDICATOR
 nlogc = n_lcon;
#endif
 m1 = m + 4*nlogc + 1;	/* "+ 1" due to nextra == 1 for sstatus */
 n = n_var + 1;
 if (!m)
	m = 1;	/* we'll add a constraint to bypass a defect in XPRESS */
 nz = nzc + 1;
 d->cstat = (int*)M1zapalloc((m1 + n)*sizeof(int));
 d->rstat = d->cstat + n;
 d->csd = suf_iput("sstatus", ASL_Sufkind_var, d->cstat);
 d->rsd = suf_iput("sstatus", ASL_Sufkind_con, d->rstat);
 d->miqp = d->nelq = 0;

 obj = (real*)Malloc((5*m1+4*n)*sizeof(real) + m*sizeof(char));
 rng = obj + n;
 memset(obj, 0, (n+m1)*sizeof(real));
 LUv = rng + m1;
 Uvx = LUv + n;
 LUrhs = Uvx + n;
 rhs = LUrhs + 2*m1;
 d->x = rhs + m1;
 d->y = d->x + n; memset(d->y, 0, sizeof(real) * m1);
 qrtype = (char*)(d->y + m1);

 A_vals = a = (real*)Malloc(nz*(sizeof(real)+sizeof(int))+(n+1)*sizeof(CStype));
 A_colstarts = ka = (CStype*)(a + nz);
 A_rownos = (int*)(ka + n + 1);

 if (!n_obj)
  nobj = 0;
 if (getopts(argv, &Oinfo))     /* Set options */
 {
	#ifdef XPRESSLIB
	 solve_result_num = 503;
	 asl->i.uinfo = "Error in $xpress_options.";
	 return;
	#endif
	 exit(1);
 }

 if (writeprob) { /* check validity */
	if (!(s = strrchr(writeprob, '.'))) {
 badname:
		fprintf(Stderr,
			"Expected \"writeprob=...\" to specify a filename ending in \".lp\"\n"
			"or \".mps\"; got \"%s\".\n", writeprob);
		exit(1);
		}
	if (!strcmp(s, ".mps"))
		wpflags = "";
	else if (!strcmp(s, ".lp"))
		wpflags = "l";
	else
		goto badname;
	}
 want_derivs = 0;
 ngents = niv /*+ nbv*/ + nlvbi + nlvci + nlvoi;
 if (mipstart && ngents + nbv)
	want_xpi0 = 5;
#if XPVERSION >= 30
  else if (optimopt[0] == 'b')
	want_xpi0 = 7;
#endif
 if (objrep > 3)
	objrep = 3;
 else if (objrep < 0)
	objrep = 0;
 qp_read(nl,ALLOW_CLP | repmap[objrep] ForceZ);
 if(logfile != NULL)
 {
  if(XPRSsetlogfile(prob,logfile))
    xperror("opening the logfile");
  XPRSsetintcontrol(prob,XPRS_OUTPUTLOG,1);        /* Chat mode */
 }
 qstype=NULL;
 msstart=NULL;
 mscols=NULL;
 qv = dref = NULL;

 i = sos ? 0 : ASL_suf_sos_ignore_sosno;
 if (!sos2)
  i |= ASL_suf_sos_ignore_amplsos;
 nsets = suf_sos(i, 0, &qstype, 0,0, &ims0, &mscols, &dref);
 ngents += nbv; /* suf_sos may adjust nbv */
 if (nsets) {
	/* It is very unlikely that ims0 would really need to be of type size_t*, */
	/* but we need to satisfy the msstart argument type in calls below. */
	msstart = (CStype*)M1alloc((nsets+1)*sizeof(CStype));
	for(i = 0; i <= nsets; ++i) {
		msstart[i] = j = ims0[i];
		if (j < 0) /* very unlikely */ {
			fprintf(Stderr,"msstart adjustment botch!\n");
			exit(1);
			}
		}
	}
 if ((ngents>0) && (optimopt[1]=='l')) {
	printf("Ignoring integrality of %d variable%s.\n",
	    ngents, ngents > 1 ? "s" : "");
	nlvoi = nbv = niv = ngents = 0;
	}
 m = n_con;
 n = n_var;

                        /* Preparing the objective function */
 obj_no = --nobj;
 og = 0;
 nq = 0;
 nelq = 0;
 q0 = 0;
 v = 0;
 qpi = 0;
 if(nobj < n_obj && nobj >= 0) {
	if (nlo && (nelqs = mqpcheckv(nobj, &qpi, &v))) {
		if (nelqs < 0) {
			nonlin(nelqs == -2,
			 "a quadratic objective involving division by 0");
			nonlin(1, "a non-quadratic nonlinear objective");
			}
		d->nelq = nelq = nelqs;
		d->miqp = ngents;
		/* discard quadratic terms below the diagonal */
		q = q1 = q0 = qpi->delsq;
		rq = rq1 = mqc1 = qpi->rowno;
		colq = qpi->colbeg;
		ncol = qpi->nc;
		cn = qpi->colno;
		for(j = 0; j < ncol;) {
			i = cn[j];
			for(qe = q0 + colq[++j]; q < qe; q++)
				if ((*rq1 = *rq++) <= i) {
					rq1++;
					*q1++ = *q;
					}
			colq[j] = q1 - q0;
			}
		nq = colq[ncol];
		mqc2 = (int*)Malloc(nq*sizeof(int));
		for(rq = mqc1, j = k = 0; j < ncol; ) {
			i = cn[j];
			for(rq1 = mqc1 + colq[++j]; rq < rq1; rq++)
				mqc2[k++] = i;
			}
		}
	og = Ograd[nobj];
	objadj = objconst(nobj);
	if(!Optimise)
		Optimise = objtype[nobj] == 0 ? XPRSminim : XPRSmaxim;
	}
  else if (nobj != -1) {
	fprintf(Stderr,"Objective %d does not exist.\n",nobj+1);
	exit(1);
	}
  else if (!Optimise)
	Optimise = XPRSminim;

  qmn = 0;
  qcrows = qcol1 = qcol2 = 0;
  if ((nqc = nlc)) {
	qcadj(&qmn, &qcrows, &qcol1, &qcol2, &qv, &v);
	++d->nelq;
	}
  a = A_vals; /* qcadj may have changed A_vals, A_rownos and A_colstarts */
  ia = A_rownos;
  ka = A_colstarts;

  if (og) {
	if (asl->i.vmap) {
		vmi = get_vminv_ASL(asl);
		do obj[vmi[og->varno]] = og->coef;
		   while((og = og->next));
		}
	else {
		do obj[og->varno] = og->coef;
		   while((og = og->next));
		}
	}

  n0 = n;
  if(objadj) {      /* Getting adjustment value */
                    /* If necessary, add a new variable */
	LUv[n]=objadj;
	Uvx[n]=objadj;
	obj[n]=1.0;
	n++;                    /* pretend 1 more column */
	A_colstarts[n]=A_colstarts[n0];
	}

      /* Make RHS, QRTYPE and RNG */

 rhs_bounds=(struct LU_bounds *)LUrhs;

 for(row=0;row<m;row++)
  if(rhs_bounds[row].upper==Infinity)
  {
   if(rhs_bounds[row].lower==negInfinity)
   {
    qrtype[row]='N';        /* Non-binding constraint */
    rhs[row]=0;
   }
   else
    {
     qrtype[row]='G';       /* >= constraint */
     rhs[row]=rhs_bounds[row].lower;
    }
  }
  else
   if(rhs_bounds[row].lower==rhs_bounds[row].upper)
    {
     qrtype[row]='E';       /* == constraint */
     rhs[row]=rhs_bounds[row].lower;
    }
   else
    if(rhs_bounds[row].lower==negInfinity)
     {
      qrtype[row]='L';      /* <= constraint */
      rhs[row]=rhs_bounds[row].upper;
     }
    else
     {
      qrtype[row]='R';      /* Range constraint */
      rhs[row]=rhs_bounds[row].upper;
      rng[row]=rhs_bounds[row].upper-rhs_bounds[row].lower;
     }

 mgcols = 0;
 if(ngents + nsets == 0)    /* Is it just LP, or MIP ? */
 {
  iret = nqc
	? XPRSloadqcqp(prob,probname,n,m,qrtype,rhs,rng,obj,
		ka,NULL,ia,a,LUv,Uvx,
		nq, mqc1, mqc2, q0,
		nqc,qcrows,qmn,qcol1,qcol2,qv)
	: nq
	? XPRSloadqp(prob,probname,n,m,qrtype,rhs,rng,obj,
		ka,NULL,ia,a,LUv,Uvx,
		nq, mqc1, mqc2, q0)
	: XPRSloadlp(prob,probname,n,m,qrtype,rhs,rng,obj,
		ka,NULL,ia,a,LUv,Uvx);
  if(iret)
    xperror("loading the problem");
 }
 else
 {
  qgtype = 0;
  if (ngents) {
	mgcols = (int *) Malloc(ngents*(sizeof(int)+1));
	qgtype = (char *)(mgcols + ngents);

	L = LUv;
	U = Uvx;
	j = 0;
	n_bv = nlvb;
	for(i = nlvb - nlvbi; i < n_bv; ++i, ++j) {	/* nonlinear integer or binary variables */
						/* in both constraints and objectives */
		qgtype[j] = L[i] == 0. && U[i] == 1. ? 'B' : 'I';
		mgcols[j] = i;
		if (L[i] < MININT)
			L[i] = MININT;
		if (U[i] > XPRS_MAXINT)
			U[i] = XPRS_MAXINT;
		}
	n_bv = nlvc;
	for(i = n_bv - nlvci; i < n_bv; ++i, ++j) {	/* nonlinear integer or binary variables */
							/* just in constraints */
		qgtype[j] = L[i] == 0. && U[i] == 1. ? 'B' : 'I';
		mgcols[j] = i;
		if (L[i] < MININT)
			L[i] = MININT;
		if (U[i] > XPRS_MAXINT)
			U[i] = XPRS_MAXINT;
		}
	n_bv += nlvo - nlvc;
	for(i = n_bv - nlvoi; i < n_bv; ++i, ++j) {	/* nonlinear integer or binary variables */
							/* just in objectives */
		qgtype[j] = L[i] == 0. && U[i] == 1. ? 'B' : 'I';
		mgcols[j] = i;
		if (L[i] < MININT)
			L[i] = MININT;
		if (U[i] > XPRS_MAXINT)
			U[i] = XPRS_MAXINT;
		}

	n_bv = n_var;
	for(i = n_bv - (nbv + niv); i < n_bv; ++i, ++j) { /* linear integer or binary variables */
		qgtype[j] = L[i] == 0. && U[i] == 1. ? 'B' : 'I';
		mgcols[j] = i;
		if (L[i] < MININT)
			L[i] = MININT;
		if (U[i] > XPRS_MAXINT)
			U[i] = XPRS_MAXINT;
		}
	}

  if( !(m) ) { /* code around XPRESS defect: add a nonbinding constraint */
	m = 1;
	A_vals[0] = 1.;
	qrtype[0] = 'N';
	rhs[0] = 0.;
	}

  iret = nqc
	? XPRSloadqcqpglobal(prob,probname,n,m,qrtype,rhs,rng,obj,
		ka,NULL,ia,a,LUv,Uvx,
		nq, mqc1, mqc2, q0,
		nqc,qcrows,qmn,qcol1,qcol2,qv,
		ngents,nsets,qgtype,mgcols,NULL/*mplim*/,
		qstype,msstart,mscols,dref)
	: nq
	? XPRSloadqglobal(prob,probname,n,m,qrtype,rhs,rng,obj,
         	ka,NULL,ia,a,LUv,Uvx,
		nq, mqc1, mqc2, q0,
		ngents,nsets,qgtype,mgcols,NULL/*mplim*/,
		qstype,msstart,mscols,dref)
	: XPRSloadglobal(prob,probname,n,m,qrtype,rhs,rng,obj,
		ka,NULL,ia,a,LUv,Uvx,
		ngents,nsets,qgtype,mgcols,NULL/*mplim*/,
		qstype,msstart,mscols,dref);
  if(iret)
	xperror("loading the problem");
#ifndef NO_XPR_CON_INDICATOR
	if (nlogc)
		indicator_constrs();
#endif /*NO_XPR_CON_INDICATOR*/

  if(optimopt[1]==0)
  {
   optimopt[1]='g';       /* Search will be global */
   mip_priorities();      /* using provided priorities */
   if (lazy && (lzd = suf_get("lazy", ASL_Sufkind_con)) && lzd->u.i)
	lazy_adj(lzd->u.i);
  }
  if (mipstart && optimopt[1] == 'g' && X0)
	XPRSloadmipsol(prob, X0, &iret);
 }
 if (mqc2)
	free(mqc2);
 if (qpi)
	free(qpi);
 if (v)
	mqpcheckv_free(&v);
 if (mgcols)
	free(mgcols);
 if (qv)
	free(qv);
 free(a);

 Do_Defer();

#if XPVERSION >= 30
 if ((X0 || pi0) && optimopt[0] == 'b') {
	i = 0;
	XPRSgetintcontrol(prob, XPRS_BARSTART, &i);
	if (i == -1)
		XPRSloadlpsol(prob, X0, 0, pi0, 0, &i);
	}
#endif

#ifdef RWA_DEBUG
 if( strstr(debugopt,"save") ) XPRSsave(prob); /* save matrix to internal file */
#endif
 atexit(killtempprob);    /* Ensure temp files are removed on exit */
 get_statuses(d);
 if (!optimopt[1] && optimopt[2] == 'n') {
	optimopt[1] = optimopt[2];
	optimopt[2] = 0;
	}
#ifdef XPRS_MSP_SOLPRB_OBJ /*{*/
	if (ngents && poolstub) {
		XPRS_msp_create(&msp);
		XPRS_msp_probattach(msp, prob);
		if (nbest > 1) {
			MSEopt = Optimise == XPRSmaxim ? XPRS_mse_maxim : XPRS_mse_minim;
			Optimise = MSEsolve;
			XPRS_mse_create(&mse);
			}
		XPRSsetintcontrol(prob, XPRS_HEUREMPHASIS, 0);
		if (dualred != 2 || dupcol != 2) {
			j = XPRSgetintcontrol(prob, XPRS_PRESOLVEOPS, &i);
			j = k = 0;
			if (dualred != 2) {
				j = 1 << 3;
				if (dualred == 1)
					k = j;
				}
			if (dupcol != 2) {
				j |= 1 << 5;
				if (dupcol == 1)
					k |= 1 << 5;
				}
			i = (i & ~j) | k;
			XPRSsetintcontrol(prob, XPRS_BARPRESOLVEOPS, i);
			}
		}
#endif /*}*/
 if (writeprob && (i = XPRSwriteprob(prob, writeprob, wpflags))) {
	fprintf(Stderr, "Could not write \"%s\".\n",
		writeprob, i);
	exit(1);
	}
}

 static int
send_statuses(dims *d)
{
  int *cstat, *rstat;
  static int map[] = {3, 1, 4, 2};

  if (!(asl->i.flags & 1) && !amplflag)
    return 0;
  if (d->miqp)
	return 1; /* avoid unsuppressable message from XPRSgetbasis */
  cstat = d->cstat;
  rstat = d->rstat;
  memset(cstat, 0, n_var*sizeof(int));
  memset(rstat, 0, n_con*sizeof(int));
  if (XPRSgetbasis(prob,rstat, cstat))
	return 1;
  stat_map(cstat, n_var, map, 3, "outgoing cstat");
  stat_map(rstat, n_con, map, 3, "outgoing rstat");
  return 0;
}

 static int
getvec(int *mrow, double *dmat, int mxelt, int *pnelt, int jvec)
{
  int nrow;
  char qrtype;          /* Row type for slack variable */
  static int mbeg[2];

  XPRSgetintattrib(prob,XPRS_ROWS, &nrow);
  if (jvec < nrow) {   /* Vector is a slack/surplus */
    *pnelt = 1;
    if (mxelt < 1) return 0;
  if (XPRSgetrowtype(prob,&qrtype, jvec, jvec)) return 1;

    *mrow = jvec; *dmat = (qrtype == 'G') ? -1.0 : 1.0;
  }
  else {                /* Vector is a structural */
    if (XPRSgetcols(prob,mbeg, mrow, dmat, mxelt, pnelt, jvec - nrow, jvec - nrow))
      return 1;
  }

  return 0;
}

 void
unpack(int *mind, double *dnz, int size, int nnz, double *dvec)
{
  double *d = dvec;

  while (size-- > 0)
    *(d++) = 0.0;
  while (nnz-- > 0)
    dvec[*(mind++)] = *(dnz++);
}

 static int
send_ray(dims *d, PBuf *B)
{
  int *cstat, *mrow, *pivrow, *rstat, i, j, junb, nb, ncol, nelt, nrow, nrseq, ns;
  double *dmat, *dvec, *unbdd, dscale;

  int pstat,lstat;

  if (optimopt[1] == 'g') {
	if (!Ray)
		Ray = 3; /* don't say "not requested" */
	return 1;
	}
  switch(Ray) {
	default: return 1;
	case 2:
		if (optimopt[0] != 'p')
			goto use_primal;
		break;
	case 1:
		XPRSgetintcontrol(prob, XPRS_PRESOLVE, &i);
		if (!i)
			break;
		XPRSsetintcontrol(prob, XPRS_PRESOLVE, 0);
 use_primal:
		j = optimopt[0];
		optimopt[0] = 'p';
		if((*Optimise)(prob,optimopt))
			xperror("optimising the problem in send_ray");
		optimopt[0] = j;
		XPRSgetintattrib(prob, XPRS_LPSTATUS, &i);
		if (i != LPSTAT_UNBOUNDED)
			Bpf(B, "\nSurprise LPSTATUS = %d computing .unbdd", i);
		XPRSgetintattrib(prob, XPRS_BARITER, &nb);
		XPRSgetintattrib(prob, XPRS_SIMPLEXITER, &ns);
		if (nb)
			Bpf(B, "\n%d extra barrier iteations computing .unbdd", nb);
		if (ns)
			Bpf(B, "\n%d extra simplex iteations computing .unbdd", ns);
	}
  XPRSgetintattrib(prob,XPRS_PRESOLVESTATE, &pstat);
  XPRSgetintattrib(prob,XPRS_LPSTATUS, &lstat);

  XPRSgetintattrib(prob,XPRS_ROWS, &nrow);
  XPRSgetintattrib(prob,XPRS_COLS, &ncol);
  XPRSgetintattrib(prob,XPRS_SPAREROWS, &nrseq);
  nrseq += nrow;

  dmat  = (double*) M1alloc((nrseq+ncol+2*nrow)*sizeof(int) + 2*nrow*sizeof(double));
  dvec  = dmat + nrow;
  rstat = (int*)(dvec + nrow);
  cstat = rstat + nrseq;
  pivrow= cstat + ncol;
  mrow  = pivrow + nrow;

  i = (n_var > ncol) ? n_var : ncol;
  unbdd = (double*)M1zapalloc(i*sizeof(double));

  if (XPRSgetunbvec(prob,&junb))
    return 1;

  if (getvec(mrow, dmat, nrow, &nelt, junb))
    return 1;

  unpack(mrow, dmat, nrow, nelt, dvec);

  if (XPRSftran(prob,dvec))
    return 1;

  if (XPRSgetbasis(prob,rstat, cstat))
    return 1;

  if (XPRSgetpivotorder(prob,pivrow))
    return 1;

  dscale = (rstat[junb] == XP_NBASUP) ? -1 : 1;

  if (junb >= nrow) unbdd[junb-nrow] = dscale;

  dscale = -dscale;

  for (i = 0; i < nrow; i++){
    j = pivrow[i];
    if (j < nrseq)
  continue; /* it's a row that's basic */
    j -=nrseq;           /* get col seq number starting at 0 */
    if (cstat[j] == XP_BASIC)
  unbdd[j] = dscale * dvec[i];
  }

  suf_rput("unbdd", ASL_Sufkind_var, unbdd);

  return 0;
}


 static int
xround(double *x, fint n, int assign, double *w)
{
	double d, dx, *xe, y;
	int m = 0;

	dx = *w;
	for(xe = x + n; x < xe; x++) {
		y = floor(*x + 0.5);
		if ((d = *x - y) != 0.) {
			if (d < 0)
				d = -d;
			if (dx < d)
				dx = d;
			m++;
			if (assign)
				*x = y;
			}
		}
	*w = dx;
	return m;
	}

#ifdef XPRS_MSP_SOLPRB_OBJ /*{*/
 static void
poolwrite(PBuf *B)
{
	PBuf B1;
	char buf[32], buf1[1200], *fname;
	const char *who;
	int i, id, j, k, n, nret, nsols, nx, *sid;
	real obj, *x;
	size_t L;

	if (mse)
		XPRS_mse_getintattrib(mse, XPRS_MSE_SOLUTIONS, &npool);
	else
		XPRS_msp_getintattrib(msp, XPRS_MSP_SOLUTIONS, &npool);
	suf_iput("npool", ASL_Sufkind_obj, &npool);
	suf_iput("npool", ASL_Sufkind_prob, &npool);
	if (npool <= 0) {
		Bpf(B, "\n%d solutions in solution pool.", npool);
		return;
		}
	B1.s = buf1;
	B1.se = buf1 + sizeof(buf1);
	n = n_var;
	L = strlen(poolstub);
	x = (real*)Malloc(npool*sizeof(int) + n*sizeof(real) + L + 32);
	sid = (int*)(x + n);
	fname = (char*)(sid + npool);
	strcpy(fname, poolstub);
	nret = nsols = 0;
	if (mse) {
		who = "XPRS_mse_getsollist";
		XPRS_mse_getsollist(mse, XPRS_MSE_METRIC_MIPOBJECT, 1, npool, sid, &nret, &nsols);
		}
	else {
		who = "XPRS_msp_getsollist";
		XPRS_msp_getsollist(msp, 0, 0, 0, 1, npool, sid, &nret, &nsols);
		}
	if (nret != npool) {
		Bpf(B, "\n%s reports %d rather than %d solutions in pool.",
				who, nret, npool);
			npool = nret;
			}
	for(i = k = 0; i < npool;) {
		id = sid[i++];
		j = nx = 0;
		XPRS_msp_getsol(msp, id, &j, x, 0, n-1, &nx);
		if (nx <= 0)
			break;
		while(nx < n)
			x[nx++] = 0.;
		if (mse)
			XPRS_mse_getsolmetric(mse, id, &j, XPRS_MSE_METRIC_MIPOBJECT, &obj);
		else
			nx = XPRS_msp_getdblattribprobsol(msp, prob, id, &j, XPRS_MSP_SOLPRB_OBJ, &obj);
		g_fmtop(buf, obj);
		Bpf(&B1, "Solution pool member %d (of %d); objective %s", i, npool, buf);
		sprintf(fname+L,"%d.sol", i);
		if (write_solf_ASL(asl, buf1, x, 0, 0, fname))
			break;
		++k;
		B1.s = buf1;
		}
	Bpf(B, "\nWrote %d solution%s from the solution pool (poolstub=\"%s\").",
		k, k == 1 ? "" : "s", poolstub);
	free(x);
	}
#endif /*}*/

/***************************************************************/
/* Check the Xpress-MP status parameter and write the solution */
/***************************************************************/
 static void
amplout(dims *d)
{
 PBuf B;
 char hbuf[1024], buf[32], *wb;
 double objvalue, w;
 int *cstatus, *rstatus;
 int didbarrier, i, ipstat, lpstat=-1, n;
 int nbit = 0, nbs, ncol, nint, nround, nrow, nsit = 0;
 real bb, *pbb, *x, *x1, *y;
 typedef struct { char *msg; int code; } Sol_info;
 static Sol_info report[]={
    { "Problem has not been loaded", 500 },
    { "Optimal solution found", 000 },
    { "Infeasible problem", 200 },
    { "Objective is worse than cutoff", 100 },
    { "Unfinished optimization", 400 },
    { "Unbounded problem", 300 },
    { "Cutoff in dual", 101 },
    { "Problem unsolved", 502 }, /* should not happen */
    { "Problem is not convex", 510 }	/*8*/
    };
 static Sol_info repglb[]={
  { "Problem has not been loaded", 500},
  { "LP has not been optimized (probably LP Infeasible)", 501},
  { "LP has been optimized", 001},
  { "Global search incomplete - no integer solution found", 401},
  { "Global search incomplete", 102},
  { "Global search complete - infeasible (no integer solution found)", 201},
  { "Global search complete", 002},
  { "Unbounded problem with some integer variables", 301}
  };

 B.s = hbuf;
 B.se = hbuf + sizeof(hbuf);
 n = n_var;
 x = d->x;
 y = d->y;
 if (nobj >= 0)
	bb = objtype[nobj] ? negInfinity : Infinity;
 Bpf(&B, "%s: ", Oinfo.bsname);
 didbarrier = (optimopt[0]=='b');
 if(optimopt[1]=='g')    /* We did a Global search */
 {
  XPRSgetintattrib(prob,XPRS_MIPSTATUS, &ipstat);
  if (ipstat <= 2)  /* ...but never started global, so .sol file not created */
  {

#ifdef XPRS_SOLUTIONFILE
   XPRSsetintcontrol(prob,XPRS_SOLUTIONFILE,0);/*  In case we try to get the soln anyway:
                             prevent looking for absent .sol file */
               /* changed from seticv(N_IFMEM,1|4|8|16)*/
#endif
  }
  if(ipstat > GLSTAT_LP_FINISHED)  /* We have a valid IP sol on the .sol file */
 {
   XPRSgetdblattrib(prob,XPRS_BESTBOUND,&bb);
   Bpf(&B, "%s", repglb[ipstat].msg);
   solve_result_num = repglb[ipstat].code;
   switch (ipstat)
   {
    case GLSTAT_UNFINISHED_NOSOL:
    case GLSTAT_FINISHED_NOSOL:
            g_fmtop(buf,bb);
            Bpf(&B, "\nBest bound determined so far %s", buf);
	   /* no break */
    case GLSTAT_MIP_UNBOUNDED:
            x=y=NULL;
            break;
    case GLSTAT_UNFINISHED_SOL:
    case GLSTAT_FINISHED_SOL:
            XPRSgetdblattrib(prob,XPRS_MIPOBJVAL,&objvalue);
            g_fmtop(buf,objvalue);
            Bpf(&B, "\nBest integer solution found %s", buf);
            XPRSgetintattrib(prob,XPRS_MIPSOLS,&nbs);
            if(nbs>1)      /* At least 2 solutions here */
               Bpf(&B, "\n%d integer solutions have been found",nbs);
            if(XPRSgetmipsol(prob,x,NULL)) /* there are no dual variables for mip solutions */
            xperror("preparing solution file");
            break;
    default:
            xperror("unrecognised global status");
   }
   XPRSgetintattrib(prob,XPRS_NODES,&nbit);
   Bpf(&B, "\n%d branch and bound node%s", nbit, (nbit!=1) ? "s" : "");
  }}
 else        /* Just LP minim or maxim - even if prob was global */
 {
  XPRSgetintattrib(prob,XPRS_LPSTATUS, &lpstat);

  Bpf(&B, "%s",report[lpstat].msg);
  solve_result_num = report[lpstat].code;
  if(lpstat==LPSTAT_OPTIMAL)
  {
   XPRSgetdblattrib(prob,XPRS_LPOBJVAL, &objvalue);
   g_fmtop(buf,objvalue);
   Bpf(&B, "\nObjective %s", buf);
  }
  XPRSgetintattrib(prob,XPRS_BARITER,&nbit);
  XPRSgetintattrib(prob,XPRS_SIMPLEXITER,&nsit);

  if (didbarrier || nbit > 0) {
	i = -1;
	XPRSgetintcontrol(prob,XPRS_BARTHREADS,&i);
	if (i > 1)
		Bpf(&B, "\n%d processors used.", i);
	}

  /* Get IIS */
  if(lpstat==LPSTAT_INFEASIBLE && iis_find==1) {
#if XPVERSION >= 21
	char *ciso, *ct, *viso, *vt;
	const char *what;
	int *ci, *cis, *cn, j, k, nc, nciis, ni, nisoc, nisov, nviis, *vi, *vis, *vn;

	i = -1;
	XPRSiisfirst(prob,1,&i);
	if (i == 0) {
		nciis = nviis = -1;
		XPRSgetiisdata(prob, 1, &nciis, &nviis, 0,0,0,0,0,0,0,0);
		ni = nciis + nviis;
		if (nciis >= 0 && nviis >= 0 && ni > 0) {
			nc = n_con;
			cn = (int*)M1zapalloc(2*(n+nc)*sizeof(int) + ni*(sizeof(int)+2));
			vn = cn + nciis;
			ci = vn + nviis;
			vi = ci + nc;
			cis = vi + n;
			vis = cis + nc;
			ct = (char*)(vis + n);
			vt = ct + nciis;
			ciso = vt + nviis;
			viso = ciso + nciis;
			XPRSgetiisdata(prob, 1, &nciis, &nviis, cn, vn, ct, vt, 0, 0, ciso, viso);
			nisoc = nisov = 0;
			for(i = 0; i < nciis; ++i) {
				k = cn[i];
				switch(ct[i]) {
				  case 'L': j = 3; break;
				  case 'G': j = 1; break;
				  case '1':
				  case '2':
				  case 'E': j = 5; break;
				  default:  j = 4;
				  }
				ci[k] = j;
				if (ciso[i]) {
					cis[k] = 1;
					++nisoc;
					}
				}
			for(i = 0; i < nviis; ++i) {
				k = vn[i];
				switch(vt[i]) {
				  case 'U': j = 3; break;
				  case 'L': j = 1; break;
				  case 'F': j = 5; break;
				  default:  j = 4;
				  }
				vi[k] = j;
				if (viso[i]) {
					vis[k] = 1;
					++nisov;
					}
				}
			suf_iput("iis", ASL_Sufkind_var, vi);
			suf_iput("iis", ASL_Sufkind_con, ci);
			Bpf(&B, "\nReturning an IIS involving %d variables and %d constraints.",
				nviis, nciis);
			if (nisoc + nisov) {
				suf_iput("iso", ASL_Sufkind_var, vis);
				suf_iput("iso", ASL_Sufkind_con, cis);
				what = "constraints and variables";
				if (!nisoc)
					what = "variables";
				else if (!nisov)
					what = "constraints";
				Bpf(&B, "\nNonzero values of suffix iso identify single %s\n"
					"whose bounds might cause the infeasibility.", what);
				}
			}
		}
#else
    XPRSiis(prob, "");
#endif
  }


/* If (lpstat==LPSTAT_INFEASIBLE or lpstat==LPSTAT_UNBOUNDED) and 'its'==0
   then no solution is available (presolve has proven infeasibility or
   unboundedness). */

  if (lpstat==LPSTAT_NONCONVEX
  || (lpstat==LPSTAT_INFEASIBLE && iis_find==1)
  ||  lpstat==LPSTAT_UNSOLVED)
	x = y = NULL;
  else if (nbit > 0
   || nsit > 0
   || !(lpstat==LPSTAT_INFEASIBLE
   || lpstat==LPSTAT_UNBOUNDED
   || lpstat==LPSTAT_NONCONVEX)) {
    if(XPRSgetlpsol(prob,x,NULL,y,NULL))
      xperror("preparing solution file");
  }
  else {

    /* there will be no basis, so set an all slack one */

    rstatus = (int *) d->y;
    cstatus = (int *) d->x;
    XPRSgetintattrib(prob,XPRS_ROWS,&nrow);
    XPRSgetintattrib(prob,XPRS_COLS,&ncol);
    for(i=0;i<nrow;i++) rstatus[i]=1;
    for(i=0;i<ncol;i++) cstatus[i]=0;
    XPRSloadbasis(prob,rstatus,cstatus);
    x = y = NULL;
  }
  if (nsit > 0)
	Bpf(&B, "\n%d simplex iteration%s", nsit, nsit == 1 ? "" : "s");
  if (nbit > 0)
	Bpf(&B, "\n%d barrier iteration%s", nbit, nbit == 1 ? "" : "s");
 }
 if (bestbound && nobj >= 0) {
	if (n_obj > 1) {
		pbb = (real*)M1zapalloc(n_obj*sizeof(real));
		pbb[nobj] = bb;
		}
	else
		pbb = &bb;
	suf_rput("bestbound", ASL_Sufkind_obj, pbb);
	suf_rput("bestbound", ASL_Sufkind_prob, &bb);
	}
 if ((nbit > 0 && !nsit) || send_statuses(d)) {
	d->csd->kind &= ~ASL_Sufkind_output;
	d->rsd->kind &= ~ASL_Sufkind_output;
	Bpf(&B, "\nNo basis.");
	}
 if (lpstat == LPSTAT_UNBOUNDED && send_ray(d,&B))
  Bpf(&B, "\nNo unbounded vector%s.", Ray ? "" : " requested");

#ifdef XPRS_CONES
 if (d->nelq) {
	i = 0;
	XPRSgetintattrib(prob,XPRS_CONES,&i);
	if (i > 0.)
		Bpf(&B, "\n%d cones detected.", i);
	}
#endif

 /* make sure integer variables near to integer values are integers */

 if (niv + nbv + nlvoi && x && optimopt[1] != 'l' && Round >= 0) {
	nround = 0;
	w = 0;
	if ((nint = niv + nbv)) {
		x1 = x + n - nint;
		nround = xround(x1, nint, Round & 1, &w);
		}
	if ((nint = nlvoi)) {
		x1 = x + (nlvo - nint);
		nround += xround(x1, nint, Round & 1, &w);
		}
	if (w < 1e-9 && !(Round & 8))
		nround = 0;
	else if (nround) {
		if (solve_result_num < 200 && !(Round & 2))
			solve_result_num += 10;
		if (Round & 4)
			nround = 0;
		else if (!(Round & 1))
			nround = -nround;
		}
	if (nround) {
		wb = "";
		if (nround < 0) {
			nround = -nround;
			wb = "would be ";
			}
		Bpf(&B, "\n%d integer variables %srounded (maxerr = %g).\n%s\n",
			nround, wb, w,
			"Reducing miptol to something < maxerr might help.");
		}
	}
#ifdef XPRS_MSP_SOLPRB_OBJ /*{*/
	if (msp)
		poolwrite(&B);
#endif /*}*/

 write_sol(hbuf,x,y,&Oinfo);
}

/***************************/
/* Some timing information */
/***************************/
static void show_times(void)
{
 Times[3] = xectim_();
 if(timing)
  printf("\nTimes (seconds):\nInput =  %g\nSolve =  %g\nOutput = %g\n",
      Times[1] - Times[0], Times[2] - Times[1],
      Times[3] - Times[2]);
}
