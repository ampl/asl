/*******************************************************************
Copyright (C) 2020 AMPL Optimization, Inc.; written by David M. Gay.

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
#include "gurobi_c.h"
#if GRB_VERSION_MAJOR >= 6
#include "nlp.h"
#include "r_qp.hd"
#endif
#include "getstub.h"
#include "signal.h"
#include "time.h"

#ifndef Sig_ret_type
#define Sig_ret_type void
#define SigRet /*nothing*/
#endif

#ifndef SigRet
#define SigRet /*nothing*/
#endif

#ifndef Sig_func_type
typedef void sig_func_type(int);
#endif

#undef ALLOW_GUROBI_SERVER
#if GRB_VERSION_MAJOR > 5 || (GRB_VERSION_MAJOR == 5 && GRB_VERSION_MINOR >= 5)
#ifndef DISALLOW_GUROBI_SERVER
#define ALLOW_GUROBI_SERVER
#endif
#endif

#if GRB_VERSION_MAJOR < 7 || (GRB_VERSION_MAJOR == 7 && GRB_VERSION_MINOR < 5)
#undef NO_MOkwf
#define NO_MOkwf
#endif

#if GRB_VERSION_MAJOR >= 6 && defined(X64_bit_pointers)
#define USE_Z  ASL_use_Z
typedef size_t CS_type;
#define GRBloadmodel GRBXloadmodel
#undef A_colstarts
#undef nzc
#define A_colstarts asl->i.A_colstartsZ_
#define nzc asl->i.nZc_
#else
#define USE_Z 0
typedef int CS_type;
#endif

 typedef struct
Dims {
	double	*c;
	double	*x;
	double	*y;
	double	*y0;
	int	*cstat;
	int	*rstat;
	SufDesc	*csd;
	SufDesc *rsd;
	char	*mb, *mbend;
	double	*xams;
	char	*fname;
	char	*fname_end;
	Option_Info *oi;
	int	kiv;
	int	missing;
	int	nc0, nv0;
	int	objprec;
	int	objsense;
	} Dims;

 typedef struct
Filename {
	struct Filename *next;
	char *name;
	} Filename;

 typedef struct
Ext_info {
	char *ext;
	int aftersol;
	} Ext_info;

 typedef struct
mint_values {
	int L;
	int U;
	int val;
	} mint_values;

 enum { /* sf_mint f values */
	set_iis		= 0,
	set_relax	= 1,
	set_mipstval	= 2,
	set_objno	= 3,
	set_sos		= 4,
	set_sos2	= 5,
	set_timing	= 6,
	set_basis	= 7,
	set_intstart	= 8,
	set_outlev	= 9,
	set_bestbound	= 10,
	set_solnsens	= 11,
	set_retmipgap	= 12,
	set_rays	= 13,
	set_priorities	= 14,
	set_feasrelax	= 15,
	set_warmstart	= 16,
	set_objrep	= 17,
	set_basisdebug	= 18,
	set_lazy	= 19,
	set_multiobj	= 20,
	set_round	= 21,
	set_kappa = 22,
	set_concurrentwinmethod = 23,
	set_maxvio = 24,
	set_work = 25
	};

 static mint_values
mint_val[26] = {
	/* set_iis */		{0, 1, 0},
	/* set_relax */		{0, 1, 0},
	/* set_mipstval */	{0, 1, 1},
	/* set_objno */		{0, 0/*n_obj*/,	1},
	/* set_sos */		{0, 1, 1},
	/* set_sos2 */		{0, 1, 1},
	/* set_timing */	{0, 3, 0},
	/* set_basis */		{0, 3, 3},
	/* set_intstart */	{0, 1, 1},
	/* set_outlev */	{0, 1, 0},
	/* set_bestbound */	{0, 1, 0},
	/* set_solnsens */	{0, 1, 0},
	/* set_retmipgap */	{0, 7, 0},
	/* set_rays */		{0, 3, 3},
	/* set_priorities */	{0, 1, 1},
	/* set_feasrelax */	{0, 6, 0},
#ifdef GRB_DBL_ATTR_VARHINTVAL
	/* set_warmstart */	{0, 4, 1},
#else
	/* set_warmstart */	{0, 2, 1},
#endif
	/* set_objrep */	{0, 3, 2},
	/* set_basisdebug */	{0, 2, 0},
	/* set_lazy */		{0, 1, 1},
	/* set_multiobj */	{0, 1, 0},
	/* set_round */		{0, 7, 7},
	/* set_kappa */   {0, 3, 0},
	/* set_concurrentwinmethod */ {0, 1, 0},
	/* set_maxvio */ {0, 1, 0},
	/* set_work */ {0, 1, 0}
	};

#define want_iis	mint_val[0].val
#define relax		mint_val[1].val
#define mipstval	mint_val[2].val
#define nobjno		mint_val[3].U
#define objno		mint_val[3].val
#define sos		mint_val[4].val
#define sos2		mint_val[5].val
#define time_flag	mint_val[6].val
#define basis		mint_val[7].val
#define intstart	mint_val[8].val
#define outlev		mint_val[9].val
#define bestbound	mint_val[10].val
#define solnsens	mint_val[11].val
#define retmipgap	mint_val[12].val
#define rays		mint_val[13].val
#define priorities	mint_val[14].val
#define feasrelax	mint_val[15].val
#define warmstart	mint_val[16].val
#define objrep		mint_val[17].val
#define basisdebug	mint_val[18].val
#define lazy		mint_val[19].val
#define multiobj	mint_val[20].val
#define Round		mint_val[21].val
#define kappa		mint_val[22].val
#define want_concurrentwinmethod mint_val[23].val
#define want_maxvio mint_val[24].val
#define want_work mint_val[25].val

 static int fixedmethod = -12345;
 static real round_reptol = 1e-9;

#ifdef ALLOW_GUROBI_SERVER /*{*/
 enum {
	ReplayInt = 0,
	ReplayDbl = 1,
	ReplayStr = 2,
	ReplayGulp = 15
	};

 typedef struct
ReplayItem {
	int kind;
	const char *code;
	union {
		double d;
		int i;
		const char *s;
		} u;
	} ReplayItem;

 typedef struct
ReplayBlock {
	struct ReplayBlock *next;
	int nitems;
	ReplayItem RI[ReplayGulp];
	} ReplayBlock;

 static ReplayBlock RB0;
 static ReplayBlock *RB = &RB0;
#endif /*}*/
 static Filename *Wflist[3];
 static GRBmodel *grbmodel;
 static char *logfile, verbuf[64];
 static double Times[5];
 static int breaking;
 static jmp_buf Jb;
#ifdef GRB_STR_PAR_CLOUDHOST
 static char *cloudhost;
 static char cloudhost_desc[] = "Host for Gurobi Instant Cloud.";
#endif

#ifdef ALLOW_GUROBI_SERVER /*{*/
#if GRB_VERSION_MAJOR >= 7
 static char *cloudid, *cloudkey, *cloudpool;
 static char
	cloudid_desc[] = "Use Gurobi Instant Cloud with this \"accessID\".",
	cloudkey_desc[] = "Use Gurobi Instant Cloud with this \"secretKey\".\n\
		Both cloudid and cloudkey are required.",
	cloudpool_desc[] = "Optional \"machine pool\" to use with Gurobi Instant Cloud.";
#endif
#if GRB_VERSION_MAJOR >= 8
 static int cloudpriority;
 static char cloudpriority_desc[] = "Priority of Cloud job, an integer >= -100 and <= 100.\n\
		Default 0.  Jobs with priority 100 run immediately -- use\n\
		caution when specifying this value.";
#endif
 static char *server, *server_passwd, *serverlic;
#if GRB_VERSION_MAJOR < 8
 static int server_port = DEFAULT_CS_PORT;
#else
 static char *server_group, *server_router;
 static int server_insecure;
#endif
 static int server_priority = DEFAULT_CS_PRIORITY;
#if GRB_VERSION_MAJOR > 9 || (GRB_VERSION_MAJOR == 9 && GRB_VERSION_MINOR >= 5)
 static int server_timeout = -1;
#else
 static double server_timeout = -1.;
#endif
 static char
	server_desc[] = "Comma-separated list of Gurobi compute servers, specified\n\
		either by name or by IP address.  Default: run Gurobi locally\n\
		(i.e., do not use a remote Gurobi server).",

#if GRB_VERSION_MAJOR >= 8
	server_group_desc[] = "Name of Compute Server Group",
	server_insecure_desc[] = "Whether to use \"insecure mode\" with the Gurobi Compute\n\
		Server.  Should be left at default value (0) unless an\n\
		administrator specifies another value.",
#else
	server_port_desc[] = "IP port to use for Gurobi compute server(s);\n\
		-1 ==> use default.",
#endif
	server_passwd_desc[] = "Password (if needed) for specified Guruobi compute server(s).",

	server_priority_desc[] = "Priority for Gurobi compute server(s).  Default = 0.\n\
		Highest priority = 100.",

#if GRB_VERSION_MAJOR >= 8
	server_router_desc[] = "Name or IP address of router for Compute Server;\n\
		can be \"\" if not used.",
#endif
	server_timeout_desc[] = "Report job as rejected by Gurobi compute server if the\n\
		job is not started within server_timeout seconds.\n\
		Default = -1 (no limit).",

	serverlic_desc[] = "Name of file containing \"server = ...\" and possibly\n\
		values for server_password, server_port, and server_timeout.\n\
		Synonyms for server: computeserver, servers.\n\
		Synonym for server_password: password.";

 static void
RPRecord(Option_Info *oi, keyword *kw, int kind, void *val)
{
	ASL *asl;
	ReplayBlock *rb, *rb1;
	ReplayItem *ri;

	rb = RB;
	if (rb->nitems >= ReplayGulp) {
		asl = oi->asl;
		rb->next = rb1 = (ReplayBlock *) M1alloc(sizeof(ReplayBlock));
		rb1->nitems = 0;
		RB = rb = rb1;
		}
	ri = &rb->RI[rb->nitems++];
	ri->code = (const char*)kw->info;
	switch(ri->kind = kind) {
	  case ReplayInt: ri->u.i = *(int*)val; break;
	  case ReplayDbl: ri->u.d = *(double*)val; break;
	  case ReplayStr: ri->u.s = (const char*)val;
	  }
	}

 static void
Replay(GRBenv *env)
{
	ReplayBlock *rb;
	ReplayItem *ri, *rie;

	for(rb = &RB0; rb; rb = rb->next) {
		rie = &rb->RI[rb->nitems];
		for(ri = rb->RI; ri < rie; ++ri) {
			switch(ri->kind) {
			  case ReplayInt:
				GRBsetintparam(env, ri->code, ri->u.i);
				break;
			  case ReplayDbl:
				GRBsetdblparam(env, ri->code, ri->u.d);
				break;
			  case ReplayStr:
				GRBsetstrparam(env, ri->code, ri->u.s);
			  }
			}
		}
	}

 static int
badserverlic(ASL *asl, const char *what, int rn)
{
	char *s;
	size_t L;

	L = strlen(serverlic) + strlen(what) + 32;
	asl->i.uinfo = s = M1alloc(L);
	snprintf(s, L, "%s serverlic file \"%s\".", what, serverlic);
	solve_result_num = rn;
	return 1;
	}

 static int
server_licread(ASL *asl)
{
	FILE *f;
	char buf[4096], *s;

	static keyword lrkeywds[] = {
	{ "computeserver", C_val, &server, "Synonym for server." },
	{ "password", C_val, &server_passwd, "Synonym for server_password." },
	{ "server", C_val, &server, server_desc },
#if GRB_VERSION_MAJOR >= 8
	{ "server_group", C_val, &server_group, server_group_desc },
#endif
	{ "server_password", C_val, &server_passwd, server_passwd_desc },
#if GRB_VERSION_MAJOR < 8
	{ "server_port", I_val, &server_port, server_port_desc },
#endif
	{ "server_priority", I_val, &server_priority, server_priority_desc },
#if GRB_VERSION_MAJOR >= 8
	{ "server_router", C_val, &server_router, server_router_desc },
#endif
#if GRB_VERSION_MAJOR > 9 || (GRB_VERSION_MAJOR == 9 && GRB_VERSION_MINOR >= 5)
	{ "server_timeout", I_val, &server_timeout, server_timeout_desc },
#else
	{ "server_timeout", D_val, &server_timeout, server_timeout_desc },
#endif
	{ "servers",  C_val, &server, "Synonym for server." } };
	static Option_Info lrinfo = { 0,0,0, lrkeywds,
		(int)(sizeof(lrkeywds)/sizeof(keyword)), ASL_OI_never_echo };

	if (!(f = fopen(serverlic, "r")))
		return badserverlic(asl, "Cannot open", 530);
	lrinfo.asl = asl;
 nextline:
	while(lrinfo.n_badopts == 0 && fgets(s = buf, sizeof(buf), f)) {
		while(*s && lrinfo.n_badopts == 0) {
			while(*s <= ' ') {
				if (!*s++)
					goto nextline;
				}
			if (*s == '#')
				goto nextline;
			s = get_opt_ASL(&lrinfo, s);
			}
		}
	fclose(f);
	if (lrinfo.n_badopts)
		return badserverlic(asl, "Bad assignment in", 531);
	return 0;
	}
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3
 static double ams_eps, ams_epsabs;
 static int ams_limit, ams_mode, ams_modeseen;
 static char *ams_stub;
#endif

#if GRB_VERSION_MAJOR >= 5 /*{*/
 static real lbpen = 1., rhspen = 1., ubpen = 1.;
#ifdef GRB_INT_PAR_TUNEOUTPUT /*{*/
  static char *tunebase;

 static int
tunerun(ASL *asl, GRBenv *env, GRBmodel *mdl, char **tunemsg)
{
	char *b, *s, *tbuf, tbuf0[4096];
	const char *fmt, *fmt2, *em, *what;
	int i, j, k, m, n;
	size_t L;
	static char prm[] = ".prm";

	what = "GRBtunemodel";
	fmt2 = "Return %d from %s: %s().";
	L = 64;
	em = 0;
	if ((i = GRBtunemodel(mdl))) {
 trouble:
		fmt = "Return %d from %s().";
		if ((em = GRBgeterrormsg(env))) {
			L += strlen(em);
			fmt = fmt2;
			}
 trouble2:
		*tunemsg = s = (char*)M1alloc(L);
		snprintf(s, L, fmt, i, what, em);
 badret:
		solve_result_num = 532;
		return 1;
		}
	n = -1;
	fmt = "No tuning results available: return %d from %s().";
	fmt2 = "No tuning results available: return %d from %s():\n%s.";
	L = 80;
	if ((i = GRBgetintattr(mdl, GRB_INT_ATTR_TUNE_RESULTCOUNT, &n))) {
		what = "GRBgetintattr";
		goto trouble;
		}
	if (n <= 0) {
		L = 32;
		fmt = "No tuning results found.";
		goto trouble2;
		}
	for(s = tunebase, b = 0; *s; ++s) {
		if (*s == '.' && s[1] == 'p' && s[2] == 'r' && s[3] == 'm') {
			b = s;
			s += 3;
			}
		}
	if (b)
		j = b - tunebase;
	else {
		j = s - tunebase;
		b = prm;
		}
	tbuf = tbuf0;
	if ((L = s - tunebase + 16) > sizeof(tbuf0))
		tbuf = (char*)M1alloc(L);
	if (j > 0)
		memcpy(tbuf, tunebase, j);
	for(k = m = 0; n > 0; ++k) {
		if ((i = GRBgettuneresult(mdl, --n))) {
			what = "GRBgettuneresult";
			if (!k)
				goto trouble;
			L = L + 80;
			*tunemsg = s = (char*)M1alloc(L);
			snprintf(s, L, "Surprise return %d from %s() after writing"
				" %d %.*s*%s files.", i, what, k, j, tunebase, b);
			goto badret;
			}
		m = snprintf(tbuf+j, L-j, "%d%s", n+1, b);
		what = "GRBwriteparams";
		if ((i = GRBwriteparams(env, tbuf)))
			goto trouble2;
		}
	*tunemsg = s = (char*)M1alloc(L = 2*m + 64);
	switch(k) {
	 case 1:
		snprintf(s, L, "Wrote tuning parameter file \"%s\".", tbuf);
		break;
	 case 2:
		snprintf(s, L,
			"Wrote tuning parameter files \"%s\" and \"%.*s2%s\".",
			tbuf, j, tunebase, b);
		break;
	 default:
		snprintf(s, L,
			"Wrote %d tuning parameter files \"%s\" ... \"%.*s%d%s\".",
			k, tbuf, j, tunebase, k, b);
	 }
	return 0;
	}
#endif /*}*/
#if GRB_VERSION_MAJOR >= 6 /*{*/

 static real pl_bigm = 1e6;
 static char pl_bigm_desc[] =
		"When some variables appear in piecewise-linear terms in the\n\
		objective and AMPL's \"option pl_lineraize 0\" is specified,\n\
		lower bounds of -pl_bigm are assumed for such variables that\n\
		are not bounded below and upper bounds of +pl_bigm are\n\
		assumed for such variables that are not bounded above.\n\
		Default = 1e6.";

 typedef struct
PLterm {
	int vno;
	int n;
	real *x, *y;
	} PLterm;

 static PLterm *
new_plterm(ASL_fg *asl, PLterm *p, expr *e, real **wp)
{
	int n, needadj, vno;
	plterm *pt;
	real L, U, adj, b, *bs, *bse, t, *x, *y;

	p->vno = vno = (expr_v*)e->R.e - var_e;
	pt = e->L.p;
	p->n = n = pt->n + 1;
	bs = pt->bs;
	bse = bs + 2*n - 4;
	p->x = x = *wp;
	p->y = y = x + n;
	*wp = y + n;
	L = LUv[vno];
	U = Uvx[vno];
	if (L <= negInfinity)
		L = -pl_bigm;
	if (L > bs[1])
		L = bs[1];
	if (L > 0.)
		L = 0.;
	if (U >= Infinity)
		U = pl_bigm;
	if (U < bse[-1])
		U = bse[-1];
	if (U < 0.)
		U = 0.;
	*x = L;
	t = adj = 0.;
	if (L <= 0. && bs[1] >= 0) {
		*y = L*bs[0];
		L = 0.;
		needadj = 0;
		}
	else {
		needadj = 1;
		*y = (L - bs[1])*bs[1];
		L = bs[1];
		}
	while(bs < bse) {
		*++x = b = bs[1];
		if (needadj && b >= 0.) {
			adj = t - L*bs[0];
			needadj = 0;
			}
		*++y = t += (b - L)*bs[0];
		L = b;
		bs += 2;
		}
	if (needadj)
		adj = t - L*bs[0];
	*++x = U;
	*++y = t + (U - L)*bs[0];
	if (adj != 0.)
		for(x = p->y; x <= y; ++x)
			*x -= adj;
	return p + 1;
	}

 static PLterm *
get_plterms(ASL_fg *asl, int objnum)
{
	PLterm *p, *rv;
	cde *o;
	expr *e, *e1, **ep, **ep1, **epe;
	int np, nt;
	real *w;

	o = &obj_de[objnum];
	e = o->e;
	np = nt = 0;
	switch(Intcast e->op) {
	  case OPPLUS:
		e1 = e->L.e;
		if (Intcast e1->op == OPPLTERM) {
			nt = 1;
			np = e1->L.p->n;
			}
		e1 = e->R.e;
		if (Intcast e1->op == OPPLTERM) {
			++nt;
			np += e1->L.p->n;
			}
		break;
	  case OPSUMLIST:
		for(ep = e->L.ep, epe = e->R.ep; ep < epe; ++ep) {
			e1 = *ep;
			if (Intcast e1->op == OPPLTERM) {
				++nt;
				np += e1->L.p->n;
				}
			}
		break;
	  case OPPLTERM:
		nt = 1;
		np = e->L.p->n;
	  }
	if (!nt)
		return 0;
	np += nt++;
	np <<= 1;
	rv = p = (PLterm*)Malloc(nt*sizeof(PLterm) + np*sizeof(real));
	w = (real*)(p + nt);
	switch(Intcast e->op) {
	  case OPPLUS:
		e1 = e->L.e;
		if (Intcast e1->op == OPPLTERM) {
			p = new_plterm(asl, p, e1, &w);
			if (nt == 1) {
				o->e = e->R.e;
				break;
				}
			}
		e1 = e->R.e;
		if (Intcast e1->op == OPPLTERM) {
			p = new_plterm(asl, p, e1, &w);
			if (nt == 2)
				goto nowzero;
			else
				o->e = e->L.e;
			}
		break;
	  case OPSUMLIST:
		for(ep = ep1 = e->L.ep, epe = e->R.ep; ep < epe; ++ep) {
			e1 = *ep;
			if (Intcast e1->op == OPPLTERM)
				p = new_plterm(asl, p, e1, &w);
			else
				*ep1++ = e1;
			}
		if (ep1 == e->L.ep)
			goto nowzero;
		e->R.ep = ep1;
		break;
	  case OPPLTERM:
		p = new_plterm(asl, p, e, &w);
nowzero:
		e->op = f_OPNUM;
		((expr_n*)e)->v = 0.;
	  }
	memset(p, 0, sizeof(PLterm));
	return rv;
	}

#endif /*}*/
#endif /*}*/

 static void
badretfmt(int rc, const char *fmt, ...)
{
	ASL *asl = cur_ASL;
	va_list ap;
	char buf[8192], *s;
	int k;

	va_start(ap, fmt);
	k = Vsnprintf(buf, sizeof(buf)-1, fmt, ap) + 1;
	if (rc) {
		solve_result_num = rc;
		memcpy(s = (char*)M1alloc(k), buf, k);
		asl->i.uinfo = s;
		}
	else
		fprintf(Stderr, "%s\n", buf);
	va_end(ap);
	}

 static void
failed(GRBenv *env, const char *what)
{
	const char *s = GRBgeterrormsg(env);
	if (s)
		badretfmt(501, "%s failed:\n\t%s.\n", what, s);
	else
		badretfmt(501, "%s failed.\n", what);
	longjmp(Jb,1);
	}

 static void enamefailed(GRBenv *env, const char *what, const char *name);

 static void
namefailed(const char *what, const char *name)
{
	badretfmt(506, "%s(\"%s\") failed.", what, name);
	longjmp(Jb,1);
	}

 static void
badival(Option_Info *oi, keyword *kw, int t, int L, int U)
{
	printf("rejecting %s %d; must be between %d and %d\n",
		kw->name, t, L, U);
	badopt_ASL(oi);
	}

 static char *
sf_mint(Option_Info *oi, keyword *kw, char *v)
{
	int t;
	char *rv;
	int i = Intcast kw->info;
	mint_values *m = mint_val + i;

	if (*v == '?' && v[1] <= ' ') {
		printf("%s=%d\n", kw->name, m->val);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = (int)strtol(v, &rv, 10);
	if (rv == v) {
		printf("Expected an integer value for %s, not \"%s\"\n",
			kw->name, v);
		badopt_ASL(oi);
		return v;
		}
	if (t < m->L || t > m->U) {
		badival(oi,kw,t,m->L,m->U);
		return rv;
		}
	m->val = t;
	return rv;
	}

 static char *
sf_dpar(Option_Info *oi, keyword *kw, char *v)
{
	GRBenv *env;
	double p[4], t;
	char *parname, *rv;

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;

	if (*v == '?' && v[1] <= ' ') {
		if (GRBgetdblparam(env, parname, &t))
			namefailed("GRBgetdblparam", parname);
		printf("%s=%.g\n", kw->name, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = strtod(v, &rv);
	if (rv == v) {
		printf("Expected a numeric value for %s, not \"%s\"\n",
			kw->name, v);
		badopt_ASL(oi);
		return v;
		}
	if (GRBsetdblparam(env, parname, t)) {
		if (GRBgetdblparaminfo(env, parname, p, p+1, p+2, p+3))
			namefailed("GRBsetdblparam", parname);
		badretfmt(506, "%s must be >= %.g and <= %.g.", kw->name, p[1], p[2]);
		badopt_ASL(oi);
		}
#ifdef ALLOW_GUROBI_SERVER
	else
		RPRecord(oi, kw, ReplayDbl, &t);
#endif
	return rv;
	}

 static void
int_rangerr(Option_Info *oi, keyword *kw)
{
	GRBenv *env;
	char *fmt, *parname;
	int p[4];

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;
	if (GRBgetintparaminfo(env, parname, p, p+1, p+2, p+3))
		namefailed("GRBsetintparam", parname);
	fmt = p[2] == p[1] + 1
		? "%s must be %d or %d."
		: "%s must be >= %d and <= %d.";
	badretfmt(506, fmt, kw->name, p[1], p[2]);
	badopt_ASL(oi);
	}

 static char *
sf_ipar(Option_Info *oi, keyword *kw, char *v)
{
	GRBenv *env;
	int t;
	char *parname, *rv;

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;

	if (*v == '?' && v[1] <= ' ') {
		if (GRBgetintparam(env, parname, &t))
			namefailed("GRBgetintparam", parname);
		printf("%s=%d\n", kw->name, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = (int)strtol(v, &rv, 10);
	if (rv == v) {
		printf("Expected an integer value for %s, not \"%s\"\n",
			kw->name, v);
		badopt_ASL(oi);
		return v;
		}
	if (GRBsetintparam(env, parname, t))
		int_rangerr(oi, kw);
#ifdef ALLOW_GUROBI_SERVER
	else
		RPRecord(oi, kw, ReplayInt, &t);
#endif
	return rv;
	}

 static char*
sf_iparlog(Option_Info *oi, keyword *kw, char *v)
{
	GRBenv *env;
	int t;
	char *parname, *rv;

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;

	if (*v == '?' && v[1] <= ' ') {
		if (GRBgetintparam(env, parname, &t))
			namefailed("GRBgetintparam", parname);
		printf("%s=%d\n", kw->name, t);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	t = (int)strtol(v, &rv, 10);
	if (rv == v) {
		printf("Expected an integer value for %s, not \"%s\"\n",
			kw->name, v);
		badopt_ASL(oi);
		return v;
		}
	if (GRBsetintparam(env, parname, t))
		int_rangerr(oi, kw);
#ifdef ALLOW_GUROBI_SERVER
	else
		RPRecord(oi, kw, ReplayInt, &t);
#endif
	return rv;
	}

#if GRB_VERSION_MAJOR >= 3
 static char*
Im_val(Option_Info *oi, keyword *kw, char *v)
{
	char *rv = sf_ipar(oi, kw, v);
	GRBgetintparam((GRBenv*)oi->uinfo, "PoolSearchMode", &ams_mode);
	ams_modeseen = 1;
	return rv;
	}
#endif

 static Filename *
fn_value(char **pv, const char *what, keyword *kw)
{
	ASL *asl = cur_ASL;
	Filename *f;
	char *s, *t, *v;
	int c, q;
	size_t L;

	v = *pv;
	q = *v;
	if (q == '"' || q == '\'') {
		s = ++v;
		for(;;) {
			if (!(c = *s)) {
				printf("Bad %s \"%s\" for %s\n", what, v, kw->name);
				*pv = v;
				return 0;
				}
			++s;
			if (c == q && *s != q)
				break;
			}
		}
	else {
		q = 0;
		for(s = v; *s > ' '; ++s);
		if (s == v) {
			printf("Missing %s for %s\n", what, kw->name);
			return 0;
			}
		}
	L = s - v;
	f = M1alloc(sizeof(Filename) + L + 1);
	f->name = t = (char*)(f + 1);
	if (q) {
		for(s = v;; ++s, ++t) {
			if ((*t = *s) == q) {
				if (*++s == q)
					continue;
				break;
				}
			}
		}
	else {
		memcpy(t, v, L);
		t += L;
		}
	*t = 0;
	*pv = s;
	return f;
	}

#if GRB_VERSION_MAJOR > 1 /*{*/
#define GRB_MAJ2(x) x

 static char *
sf_spar(Option_Info *oi, keyword *kw, char *v)
{
	Filename *f;
	GRBenv *env;
	char *parname, tbuf[GRB_MAX_STRLEN + 8];

	env = (GRBenv*)oi->uinfo;
	parname = (char*)kw->info;

	if (*v == '?' && v[1] <= ' ') {
		memset(tbuf, 0, sizeof(tbuf));
		if (GRBgetstrparam(env, parname, tbuf))
			namefailed("GRBgetstrparam", parname);
		printf("%s=\"%s\"\n", kw->name, tbuf);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	if ((f = fn_value(&v, "value", kw))
	 && GRBsetstrparam(env, parname, f->name))
		enamefailed(env, "GRBsetstrparam", parname);
#ifdef ALLOW_GUROBI_SERVER
	else
		RPRecord(oi, kw, ReplayStr, f->name);
#endif
	return v;
	}
#else
#define GRB_MAJ2(x) /*nothing*/
#endif /*}*/

 static char *
sf_wfile(Option_Info *oi, keyword *kw, char *v)
{
	Ext_info *e;
	Filename *f, **pf;
	char *dot, *t;
	int i, q;
	static Ext_info W_ext[] = {
		{"bas",1},
		{"fix_lp",2},
		{"fix_mps",2},
		{"lp",0},
		{"mps",0},
		{"prm",0},
		{"rew",0},
		{"rlp",0},
		{"sol",1},
		{0,0}};

	q = *v;
	if (q == '?' && v[1] <= ' ') {
		for(i = 0; i < 3; ++i)
			for(f = Wflist[i]; f; f = f->next)
				printf("%s=\"%s\"\n", kw->name, f->name);
		oi->option_echo &= ~ASL_OI_echothis;
		return v + 1;
		}
	if (!(f = fn_value(&v, "file name", kw))) {
		printf("Bad file name \"%s\" for %s\n", v, kw->name);
		badopt_ASL(oi);
		return v;
		}
	dot = 0;
	for(t = f->name; *t; )
		if (*t++ == '.')
			dot = t;
	if (dot)
		for(e = W_ext; e->ext; ++e)
			if (!strcmp(e->ext, dot))
				goto good_ext;
	printf("File name for %s must end in one of\n", kw->name);
	for(e = W_ext; e->ext; ++e)
		printf("\t.%s\n", e->ext);
	badopt_ASL(oi);
	goto ret;
 good_ext:
	pf = &Wflist[e->aftersol];
	f->next = *pf;
	*pf = f;
	if (e->aftersol == 2)	/* replace _ by . */
		dot[3] = '.';
 ret:
	return v;
	}

#if GRB_VERSION_MAJOR >= 3 /*{*/

 static char *
sf_pf(Option_Info *oi, keyword *kw, char *v)
{
	FILE *f;
	Filename *fn;
	GRBenv *env;
	char buf[512], *fname, *s, *s1, *se;
	int lineno;
	static char extra[]  = "Line %d of paramfile \"%s\":\n\t"
		"expected a name and value, but got \"%s\".";
	static char failed[] = "Line %d of paramfile \"%s\":\n\t"
		"GRBsetstrparam(\"Dummy\", \"%s\") failed:\n\t%s.";
	static char missing[] = "Missing value in line %d of paramfile \"%s\".";

	env = (GRBenv*)oi->uinfo;

	if ((fn = fn_value(&v, "value", kw))) {
		fname = fn->name;
		if (!(f = fopen(fname, "r"))) {
			badretfmt(511, "Cannot open paramfile \"%s\".", fname);
			longjmp(Jb,1);
			}
		lineno = 0;
 nextline:
		while(fgets(buf, sizeof(buf), f)) {
			++lineno;
			for(s = buf; *s <= ' '; ++s)
				if (!*s)
					goto nextline;
			if (*s == '#')
				goto nextline;
			for(s1 = s; *++s1 > ' '; );
			while(*s1 <= ' ')
				if (!*s1++) {
					badretfmt(512, missing, lineno, fname);
					longjmp(Jb,1);
					}
			for(se = s1; *++se; );
			while(--se > s1 && *se <= ' ');
			se[1] = 0;
			while(*++s1 > ' ');
			if (*s1) {
				for(se = s1; *++se; ) {
					if (*se > ' ') {
						badretfmt(513, extra, lineno, fname, s);
						longjmp(Jb,1);
						}
					}
				*s1 = 0;
				}
			if (GRBsetstrparam(env, "Dummy", s)) {
				badretfmt(514, failed, lineno, fname, s, GRBgeterrormsg(env));
				longjmp(Jb,1);
				}
			}
		fclose(f);
		}
	return v;
	}

 static char aggfill_desc[] = "Amount of fill allowed during aggregation during\n\
			gurobi's presolve "
#if GRB_VERSION_MAJOR > 5 || (GRB_VERSION_MAJOR == 5 && GRB_VERSION_MINOR >= 1)
	"(default -1).";
#else
	"(default 10).";
#endif
#endif /*}*/

 static char aggregate_desc[] = "Whether to use aggregation during Gurobi presolve:\n\
			0 = no (sometimes reduces numerical errors)\n\
			1 = yes (default).";

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char ams_eps_desc[] = "Relative tolerance for reporting alternate MIP solutions\n\
			(default 0 = no limit).";

 static char ams_epsabs_desc[] = "Absolute tolerance for reporting alternate MIP solutions\n\
			(default 0 = no limit).";

 static char ams_limit_desc[] = "Limit on number of alternate MIP solutions written,\n\
			with no limit when ams_limit = 0 (default).";

#ifdef GRB_INT_PAR_POOLSEARCHMODE
 static char ams_mode_desc[] = "Search mode for MIP solutions when ams_stub is specified\n\
			to request finding several alternative solutions:\n\
			0 = ignore ams_stub; seek just one optimal solution\n\
			1 = make some effort at finding additional solutions\n\
			2 = seek \"ams_limit\" best solutions (default).";
#endif

 static char ams_stub_desc[] = "Stub for alternative MIP solutions, written to files with\n\
			names obtained by appending \"1.sol\", \"2.sol\", etc., to\n\
			ams_stub.  The number of such files written is affected\n\
			by four keywords:\n\
			  ams_mode specifies how much effort to expend;\n\
			  ams_limit gives the maximum number of files written,\n\
				with no limit when ams_limit = 0 (default);\n\
			  ams_eps gives a relative tolerance on the objective\n\
				values of alternative solutions; and\n\
			  ams_epsabs gives an absolute tolerance on how much\n\
				worse the objectives can be.\n\
			The number of alternative MIP solution files written is\n\
			returned in suffix npool on the objective and problem.\n\
			Suffix poolignore can be used to keep only one solution\n\
			with the best objective value among solutions that\n\
			differ only in variables where the suffix is 1.";

 static char barconvtol_desc[] = "Tolerance on the relative difference between the\n\
			primal and dual objectives for stopping the barrier\n\
			algorithm (default 1e-8).";

 static char barcorrectors_desc[] = "Limit on the number of central corrections done in\n\
			each barrier iteration (default -1 = automatic choice)";

#ifdef GRB_INT_PAR_BARHOMOGENEOUS /*{*/
 static char barhomogeneous_desc[] = "Whether to use the homogeneous barrier algorithm\n\
		(e.g., when method=2 is specified):\n\
			-1 = only when solving a MIP node relaxation (default)\n\
			 0 = never\n\
			 1 = always.\n\
		The homogeneous barrier algorithm can detect infeasibility or\n\
		unboundedness directly, without crossover, but is a bit slower\n\
		than the nonhomogeneous barrier algorithm.";
#endif /*}*/

 static char bariterlimit_desc[] = "Limit on the number of barrier iterations (default 1000).";

 static char barorder_desc[] = "Ordering used to reduce fill in sparse-matrix factorizations\n\
				during the barrier algorithm:\n\
			-1 = automatic choice (default)\n\
			 0 = approximate minimum degree\n\
			 1 = nested dissection.";
#endif /*}*/

#ifdef GRB_DBL_PAR_BARQCPCONVTOL
 static char barqcptol_desc[] = "Convergence tolerance on the relative difference between\n\
		primal and dual objective values for barrier algorithms when\n\
		solving problems with quadratic constraints (default 1e-6).";
#endif

 static char basis_desc[] = "Whether to use or return a basis:\n\
			0 = no\n\
			1 = use incoming basis (if provided)\n\
			2 = return final basis\n\
			3 = both (1 + 2 = default)."
#if GRB_VERSION_MAJOR >= 5
			"\n\
		For problems with integer variables and quadratic constraints,\n\
		basis = 0 is assumed quietly unless qcpdual=1 is specified."
#endif
;
 static char basisdebug_desc[] = "Whether to honor basis and solnsens when an "
		"optimal solution\n\
		was not found:\n\
			0 = only if a feasible solution was found (default)\n\
			1 = yes\n\
			2 = no.";

#ifdef GRB_DBL_PAR_BESTBDSTOP
 static char bestbndstop_desc[] = "Stop once the best bound on the objective value\n\
		is at least as good as this value.";
#endif

 static char bestbound_desc[] = "Whether to return suffix .bestbound for the "
		"best known bound\n\
		on the objective value:\n\
			0 = no (default)\n\
			1 = yes\n\
		The suffix is on the objective and problem and is +Infinity\n\
		for minimization problems and -Infinity for maximization\n\
		problems if there are no integer variables or if an integer\n\
		feasible solution has not yet been found.";

#ifdef GRB_DBL_PAR_BESTOBJSTOP
 static char bestobjstop_desc[] = "Stop after a feasible solution with objective value\n\
		at least as good as this value has been found.";
#endif

#ifdef GRB_INT_PAR_BQPCUTS
 static char bqpcuts_desc[] = "Whether to enable Boolean Quadric Polytope cut generation:\n\
			-1 = automatic choice (default)\n\
			 0 = disallow BQP cuts\n\
			 1 = enable moderate BQP cut generation\n\
			 2 = enable aggressive BQP cut generation.\n\
		Values 1 and 2 override the \"cuts\" keyword.";
#endif

#ifdef GRB_INT_PAR_CONCURRENTMIP
 static char concurrentmip_desc[] =
		"How many independent MIP solves to allow at once when multiple\n\
		threads are available.  The available threads are divided as\n\
		evenly as possible among the concurrent solves.  Default = 1.";
#endif

#ifdef GRB_INT_ATTR_CONCURRENTWINMETHOD
 static char concurrentwinmethod_desc[] =
		"Whether to return the winning method after a continuous\n\
		problem has been solved with concurrent optimization:\n\
			0 = do not return (default)\n\
			1 = return as problem suffix \"concurrentwinmethod\".\n\
		See option \"method\" for a description of the returned values.";
#endif

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char crossover_desc[] = "How to transform a barrier solution to a basic one:\n\
			-1 = automatic choice (default)\n\
			 0 = none: return an interior solution\n\
			 1 = push dual vars first, finish with primal simplex\n\
			 2 = push dual vars first, finish with dual simplex\n\
			 3 = push primal vars first, finish with primal simplex\n\
			 4 = push primal vars first, finish with dual simplex";
 static char crossoverbasis_desc[] = "strategy for initial basis construction during crossover:\n\
			0 = favor speed (default)\n\
			1 = favor numerical stability.";
#endif /*}*/

#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
 static char branchdir_desc[] = "Which child node to explore first when branching:\n\
			-1 = explore \"down\" branch first\n\
			 0 = explore \"most promising\" branch first (default)\n\
			 1 = explore \"up\" branch first.";
#endif /*}*/
 static char cutagg_desc[] = "Maximum number of constraint aggregation passes\n\
		during cut generation (-1 = default = no limit);\n\
		overrides \"cuts\".";

 static char cutoff_desc[] = "If the optimal objective value is worse than cutoff,\n\
		report \"objective cutoff\" and do not return a solution.\n\
		Default: +Infinity for minimizing, -Infinity for maximizing.";

#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
 static char cutpasses_desc[] = "Maximum number of cutting-plane passes to do\n\
		during root-cut generation; default = -1 ==> automatic choice.";
#endif /*}*/

#ifdef GRB_INT_PAR_DEGENMOVES
 static char degenmoves_desc[] = "Limit on the number of degenerate simplex moves -- for use\n\
		when too much time is taken after solving the initial root\n\
		relaxation of a MIP problem and before cut generation or root\n\
		heuristics have started.  Default -1 ==> automatic choice.";
#endif

#ifdef GRB_INT_PAR_DISCONNECTED
 static char disconnected_desc[] =
		"Whether to exploit independent MIP sub-models:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = use moderate effort\n\
			 2 = use aggressive effort.";
#endif

#ifdef GRB_INT_PAR_DUALREDUCTIONS  /*{*/
 static char dualreductions_desc[] = "Whether Gurobi's presolve should use dual reductions, which\n\
		may be useful on a well-posed problem but can prevent\n\
		distinguishing whether a problem is infeasible or unbounded:\n\
			0 = no\n\
			1 = yes (default).";
#endif /*}*/

 static char cuts_desc[] = "Global cut generation control, valid unless overridden\n\
		by individual cut-type controls:\n\
			-1 = automatic choice (default)\n\
			 0 = no cuts\n\
			 1 = conservative cut generation\n\
			 2 = aggressive cut generation"
			GRB_MAJ2("\n\t\t\t 3 = very aggressive cut generation.")
		;

#if GRB_VERSION_MAJOR >= 5 /*{*/
 static char feasrelax_desc[] = "Whether to modify the problem into a feasibility\n\
		relaxation problem:\n\
			0 = no (default)\n\
			1 = yes, minimizing the weighted sum of violations\n\
			2 = yes, minimizing the weighted count of violations\n\
			3 = yes, minimizing the sum of squared violations\n\
			4-6 = same objective as 1-3, but also optimize the\n\
				original objective, subject to the violation\n\
				objective being minimized\n\
		Weights are given by suffixes .lbpen and .ubpen on variables\n\
		and .rhspen on constraints (when positive), else by keywords\n\
		lbpen, ubpen, and rhspen, respectively (default values = 1).\n\
		Weights <= 0 are treated as Infinity, allowing no violation.";
#ifdef GRB_DBL_PAR_FEASRELAXBIGM
 static char feasrelaxbigm_desc[] =
		"Value of \"big-M\" sometimes used with constraints when doing\n\
		a feasibility relaxation.  Default = 1e6.";
#endif
#endif /*}*/

 static char feastol_desc[] = "Primal feasibility tolerance (default 1e-6).";

 static char fixedmethod_desc[] =
		"Value of \"method\" to use when seeking a basis for MIP problems\n\
		when \"basis=2\" or (the default) \"basis=3\" has been specified.\n\
		Default: if \"method\" is 0 or 1 then \"method\" else 1.";

 static char gomory_desc[] = "Maximum number of Gomory cut passes during cut generation\n\
		(-1 = default = no limit); overrides \"cuts\".";

 static char heurfrac_desc[] = "Fraction of time to spend in MIP heuristics (default 0.05)";

 static char iisfind_desc[] = "Whether to return an IIS (irreducible infeasible set of\n\
		constraints and variable bounds, via suffix .iis on\n\
		constraints and variables) when the problem is infeasible:\n\
			0 = no (default)\n\
			1 ==> yes.\n\
		When iisfind=1 and the problem is infeasible, suffixes\n\
		iisforce, iisforcelb, and iisforceub can be used to influence\n\
		the IIS found, either forcing an entity to be in or not in the\n\
		IIS or letting the algorithm decide:\n\
			0 = algorithm decides (default)\n\
			1 = force to be excluded from the IIS\n\
			2 = force to be included in the IIS.";

#if GRB_VERSION_MAJOR > 1 /*{*/
 static char iismethod_desc[] = "Which method to use when finding an IIS (irreducible\n\
		infeasible set of constraints, including variable bounds):\n\
			-1 = automatic choice (default)\n\
			 0 = often faster than method 1\n\
			 1 = can find a smaller IIS than method 0.";
#endif /*}*/

#ifdef GRB_INT_PAR_INFPROOFCUTS
 static char infproofcuts_desc[] = "Whether to generate infeasibility proof cuts:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = moderate cut generation\n\
			 2 = aggressive cut generation.";
#endif
#ifdef GRB_INT_PAR_INTEGRALITYFOCUS
 static char integralityfocus_desc[] = "Setting this parameter to 1 requests the solver to work\n\
		harder at finding solutions that are still (nearly) feasible\n\
		when all integer variables are rounded to exact integral\n\
		values:\n\
			0 = no (default)\n\
			1 = yes.";
#endif
#ifdef GRB_DBL_PAR_IMPROVESTARTGAP
 static char isg_desc[] = "Optimality gap below which the MIP solver switches from\n\
		trying to improve the best bound to trying to find better\n\
		feasible solutions (default 0).";
#endif

#ifdef GRB_DBL_PAR_IMPROVESTARTNODES
 static char isn_desc[] = "Number of MIP nodes after which the solution strategy\n\
		will change from improving the best bound to finding better\n\
		feasible solutions (default Infinity).";
#endif

#ifdef GRB_DBL_PAR_IMPROVESTARTTIME
 static char ist_desc[] = "Execution seconds after which the MIP solver switches from\n\
		trying to improve the best bound to trying to find better\n\
		feasible solutions (default Infinity).";
#endif

 static char intfeastol_desc[] = "Feasibility tolerance for integer variables (default 1e-05).";

 static char intstart_desc[] = "When there are integer variables, whether to use\n\
		an initial guess (if available):\n\
			0 = no\n\
			1 = yes (default).";
#if GRB_VERSION_MAJOR >= 8
 static char kappa_desc[] = "Whether to return the estimated condition number (kappa) of\n\
		the optimal basis (default 0): sum of\n\
			1 = report kappa in Gurobi's result message;\n\
			2 = return kappa in the solver-defined suffix kappa on\n\
			    the objective and problem.\n\
	The request is ignored when there is no optimal basis.";
#endif
#if GRB_VERSION_MAJOR >= 6
 static char lazy_desc[] = "Whether to honor suffix .lazy on linear constraints in\n\
		problems with binary or integer variables:\n\
			0 = no (ignore .lazy)\n\
			1 = yes (default).\n\
		Lazy constraints are indicated with .lazy values of "
#if (GRB_VERSION_MAJOR == 9 && GRB_VERSION_MINOR >= 1) || GRB_VERSION_MAJOR > 9
		"-1, "
#endif
		"1, 2,\n\
		or 3 and are ignored until a solution feasible to the\n\
		remaining constraints is found.  What happens next depends\n\
		on the value of .lazy:\n\
		      -1 ==> treat the constraint as a user cut; the\n\
			      constraint must be redundant with respect to the\n\
			      rest of the model, although it can cut off LP\n\
			      solutions;\n\
			1 ==> the constraint may still be ignored if another\n\
			      lazy constraint cuts off the current solution;\n\
			2 ==> the constraint will henceforth be enforced if it\n\
			      is violated by the current solution;\n\
			3 ==> the constraint will henceforth be enforced.";
#endif
#ifdef GRB_INT_PAR_LIFTPROJECTCUTS
 static char liftprojectcuts_desc[] = "Whether to generate lift-and-project cuts:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = moderate cut generation\n\
			 2 = aggressive cut generation.";
#endif
 static char logfile_desc[] = "Name of file to receive log lines (default: none)"
		GRB_MAJ2(";\n\t\t\timplies outlev = 1.")
		;

 static char logfreq_desc[] = "Interval in seconds between log lines (default 5).";
#ifdef GRB_INT_PAR_LPWARMSTART
 static char lpwarmstart_desc[] = "Controls how gurobi uses warm-start in LP optimization:\n\
			0 = ignore information\n\
			1 = use information to solve the original problem\n\
			2 = use the (crushed) information to solve the\n\
			    presolved problem.\n\
		Note that if presolve is disabled, 1 prioritizes basis\n\
		statuses, 2 prioritizes start vectors.  Default = 1.";
#endif
 static char maxmipsub_desc[] = "Maximum number of nodes for RIMS heuristic to explore\n\
		on MIP problems (default 500).";

#ifdef GRB_DBL_ATTR_MAX_VIO
 static char maxvio_desc[] = "Whether to return the maximum of all (unscaled) violations to\n\
		the current problem:\n\
			0 = do not return (default)\n\
			1 = return in the problem suffix \"maxvio\".";
#endif

#ifdef GRB_DBL_PAR_MEMLIMIT
 static char memlimit_desc[] = "Maximum amount of memory available to Gurobi (in GB, default\n\
		no limit). The solution will fail if more memory is needed.";
#endif

#ifdef GRB_DBL_PAR_SOFTMEMLIMIT
 static char softmemlimit_desc[] = "Maximum amount of memory available to Gurobi (in GB; default\n\
		= no limit). The solution is returned even if more memory\n\
		could be used.";
#endif

#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
 static char minrelnodes_desc[] = "Number of nodes for the Minimum Relaxation heuristic to\n\
		explore at the MIP root node when a feasible solution has not\n\
		been found by any other heuristic; default -1 ==> automatic\n\
		choice.";
#endif /*}*/
#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char mipfocus_desc[]  = "MIP solution strategy:\n\
			0 = balance finding good feasible solutions and\n\
			    proving optimality (default)\n\
			1 = favor finding feasible solutions\n\
			2 = favor proving optimality\n\
			3 = focus on improving the best objective bound.";
#endif /*}*/

 static char mipstart_desc[] = "Whether to use initial guesses in problems with\n\
		integer variables:\n\
			0 = no\n\
			1 = yes (default).";

#ifdef GRB_INT_PAR_MIQCPMETHOD
 static char miqcpmethod_desc[] =
		"Method for solving mixed-integer quadratically constrained\n\
		(MIQCP) problems:\n\
			-1 = automatic choice (default)\n\
			 0 = solve continuous QCP relaxations at each node\n\
			 1 = use linearized outer approximations.";
#endif

#if GRB_VERSION_MAJOR >= 7 /*{*/
 static char multiobj_desc[] = "Whether to do multi-objective optimization:\n\
			0 = no (default)\n\
			1 = yes\n\
		When multiobj = 1 and several objectives are present, suffixes\n\
		.objpriority, .objweight, .objreltol, and .objabstol on the\n\
		objectives are relevant.  Objectives with greater .objpriority\n\
		values (integer values) have higher priority.  Objectives with\n\
		the same .objpriority are weighted by .objweight.  Objectives\n\
		with positive .objabstol or .objreltol are allowed to be\n\
		degraded by lower priority objectives by amounts not exceeding\n\
		the .objabstol (absolute) and .objreltol (relative) limits.\n\
		The objective must all be linear."
#ifndef NO_MOkwf
						"  Objective-specific\n\
		convergence tolerances and method values may be assigned via\n\
		keywords of the form obj_n_name, such as obj_1_method for the\n\
		first objective."
#endif
		;

 static char multiobjmethod_desc[] = "Choice of optimization algorithm for lower-priority\n\
		objectives:\n\
			-1 = automatic choice (default)\n\
			 0 = primal simplex\n\
			 1 = dual simplex\n\
			 2 = ignore warm-start information; use the algorithm\n\
				specified by the method keyword.\n\
		The method keyword determines the algorithm to use for the\n\
		highest priority objective.";

 static char multiobjpre_desc[] = "How to apply Gurobi's presolve when doing\n\
		multi-objective optimization:\n\
			-1 = automatic choice (default)\n\
			 0 = do not use Gurobi's presolve\n\
			 1 = conservative presolve\n\
			 2 = aggressive presolve, which may degrade lower-\n\
				priority objectives.";
#endif /*}*/

#ifdef GRB_INT_PAR_NETWORKALG /* 10.0 */
 static char networkalg_desc[] = "Controls whether to use network simplex, if an LP is a\n\
		network problem:\n\
			-1 = automatic choice (default)\n\
			 0 = do not use network simplex\n\
			 1 = use network sinmplex.";
#endif

#ifdef GRB_INT_PAR_NLPHEUR /* 9.5 */
 static char nlpheur_desc[] = "Controls the NLP heuristic, affecting non-convex quadratic\n\
		problems:\n\
			 0 = Do not apply heuristic\n\
			 1 = Apply heuristic (default).";
#endif

#if GRB_VERSION_MAJOR >= 4 /*{*/
 static char nodemethod_desc[] = "Algorithm used to solve relaxed MIP node problems:\n\
			-1 = automatic choice (default)\n\
			 0 = primal simplex\n\
			 1 = dual simplex (default)\n\
			 2 = barrier.";
#endif /*}*/

#ifdef GRB_INT_PAR_NONCONVEX
 static char nonconvex_desc[] = "How to handle non-convex quadratic objectives and constraints:\n\
			-1 = default choice (currently the same as 1)\n\
			 0 = complain about nonquadratic terms\n\
			 1 = complain if Gurobi's presolve cannot discard or\n\
				eliminate nonquadratic terms\n\
			 2 = translate quadratic forms to bilinear form and use\n\
				spatial branching.";
#endif

#ifdef GRB_DBL_PAR_NORELHEURTIME
 static char norelheurtime_desc[] = "Limits the amount of time spent in the NoRel heuristic;\n\
		see the description of norelheurwork for details.  This\n\
		parameter will introduce non determinism; use norelheurwork\n\
		for deterministic results.  Default 0.";
 static char norelheurwork_desc[] = "Limits the amount of work spent in the NoRel heuristic.\n\
		This heuristics searches for high-quality feasible solutions\n\
		before solving the root relaxation.  The work metrix is hard\n\
		to define precisely, as it depends on the machine.  Default 0.";
#endif

#ifdef GRB_INT_PAR_NUMERICFOCUS
 static char numericfocus_desc[] = "How much to try detecting and managing numerical issues:\n\
			0 = automatic choice (default)\n\
			1-3 = increasing focus on more stable computations.";
#endif

#ifdef GRB_INT_PAR_OBBT
 static char obbt_desc[] = "Controls aggressiveness of Optimality-Based Bound Tightening:\n\
			-1 = automatic choice (default)\n\
			 0 = do not use OBBT\n\
			 1 = low aggressiveness\n\
			 2 = moderate  aggressiveness\n\
			 3 = high aggressiveness.";
#endif

 static char objno_desc[] = "Objective to optimize:\n\
			0 = none\n\
			1 = first (default, if available),\n\
			2 = second (if available), etc.";

#if GRB_VERSION_MAJOR >= 4 /*{*/
 static char objrep_desc[] = "Whether to replace\n\
			minimize obj: v;\n\
		with\n\
			minimize obj: f(x)\n\
		when variable v appears linearly in exactly one constraint\n\
		of the form\n\
			s.t. c: v >= f(x);\n\
		or\n\
			s.t. c: v == f(x);\n\
		Possible objrep values:\n\
			0 = no\n\
			1 = yes for v >= f(x)\n\
			2 = yes for v == f(x) (default)\n\
			3 = yes in both cases\n\
		For maximization problems, \">= f(x)\" is changed to \"<= f(x)\"\n\
		in the description above.";
#endif /*}*/
#if GRB_VERSION_MAJOR > 1 /*{*/
 static char objscale_desc[] = "How to scale the objective:\n\
			0 ==> automatic choice (default)\n\
			negative >= -1 ==> divide by max abs. coefficient\n\
					   raised to this power\n\
			positive ==> divide by this value.";
#endif /*}*/

 static char opttol_desc[] = "Optimality tolerance on reduced costs (default 1e-6).";

 static char outlev_desc[] = "Whether to write Gurobi log lines (chatter) to stdout:\n\
			0 = no (default)\n\
			1 = yes (see logfreq).";

#define Overrides_cuts "Overrides \"cuts\"; choices as for \"cuts\"."
 static char overrides_cuts[] = Overrides_cuts;

 static char perturb_desc[] = "Magnitude of simplex perturbation (when needed; default 2e-4).";

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char param_desc[] = "General way to specify values of both documented and\n\
		undocumented Gurobi parameters; value should be a quoted\n\
		string (delimited by ' or \") containing a parameter name, a\n\
		space, and the value to be assigned to the parameter.  Can\n\
		appear more than once.  Cannot be used to query current\n\
		parameter values.";

 static char paramfile_desc[] = "Name of file (surrounded by 'single' or \"double\" quotes if the\n\
		name contains blanks) of parameter names and values for them.\n\
		Lines that start with # are ignored.  Otherwise, each nonempty\n\
		line should contain a name and a value, separated by a space.";

 static char predeprow_desc[] = "Whether Gurobi's presolve should remove linearly\n\
		dependent constraint-matrix rows:\n\
			-1 = only for continuous models (default)\n\
			 0 = never\n\
			 1 = for all models.";

 static char predual_desc[] = "Whether gurobi's presolve should form the dual of a\n\
				continuous model:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes\n\
			 2 = form both primal and dual and use two threads to\n\
			     choose heuristically between them.";
#endif /*}*/

#ifdef GRB_INT_PAR_PARTITIONPLACE
 static char partitionplace_desc[] =
		"Whether and how to use the .partition suffix on variables\n\
		in the partition heuristic for MIP problems: sum of\n\
			 1 ==> when the branch-and-cut search ends\n\
			 2 ==> at nodes in the branch-and-cut search\n\
			 4 ==> after the root-cut loop\n\
			 8 ==> at the start of the root-cut loop\n\
			16 ==> before solving the root relaxation.\n\
		Default = 15.  Values of .parition determine how variables\n\
		participate in the partition heuristic.  Variables with\n\
			.partition = -1 are always held fixed;\n\
			.partition = 0 can vary in all sub-MIP models;\n\
			.partition > 0 can vary only in in that sub-MIP model.\n\
		The partition heuristic is only run when partition_place is\n\
		between 1 and 31 and some variables have suitable .partition\n\
		suffix values.";
#endif

#ifdef GRB_INT_PAR_PREMIQCPFORM
 static char premiqcpform_desc[] =
		"For mixed-integer quadratically constrained (MIQCP) problems,\n\
		how Gurobi should transform quadratic constraints:\n\
			-1 = automatic choice (default)\n\
			 0 = retain MIQCP form\n\
			 1 = transform to second-order cone contraints\n\
			 2 = transform to rotated cone constraints\n\
		Choices 0 and 1 work with general quadratic constraints.\n\
		Choices 1 and 2 only work with constraints of suitable forms.";
#endif

#ifdef GRB_DBL_PAR_PRESOS1BIGM
 static char
	presos1bigm_desc[] = "Big-M for converting SOS1 constraints to binary form:\n\
			-1 = automatic choice (default)\n\
			 0 = no conversion\n\
		Large values (e.g., 1e8) may cause numeric trouble.",
	presos2bigm_desc[] = "Big-M for converting SOS2 constraints to binary form:\n\
			-1 = automatic choice\n\
			 0 = no conversion (default)\n\
		Large values (e.g., 1e8) may cause numeric trouble.";
#endif

#ifdef GRB_INT_PAR_PRESOS1ENCODING
 static char
	 presos1enc_desc[] = "Encoding used for SOS1 reformulation:\n\
			-1 = automatic choice (default)\n\
			 0 = multiple choice model, produces an LP relaxation\n\
			     that is easy to solve\n\
			 1 = incremental model, reduces the amount of branching\n\
			 2 = formulation whose LP relaxation is easier to solve\n\
			 3 = formulation with better branching behavior,\n\
			     requires sum of the variables == 1.\n\
		Options 0 and 1 produce reformulations that are linear in\n\
		size; options 2 and 3 use reformulation logarithmic in size.\n\
		Option 2 and 3 apply only when all the variables are >=0.",
	 presos2enc_desc[] = "Encoding used for SOS2 reformulation, see presos1enc.";
#endif

#ifdef GRB_STR_ATTR_SERVER /*{*/
 static char
#if GRB_VERSION_MAJOR >= 6
	pool_distmip_desc[] = "Number of machines in the server pool (if specified by\n\
		pool_servers) to use for solving each MIP instance.",
#endif
	pool_mip_desc[] =
		"Number of independent MIP jobs (default 0) to generate and\n\
		solve using the server pool (if specified by pool_servers).\n\
		Gurobi automatically chooses different algorithm parameter\n\
		values for each job.",
	pool_password_desc[] = "Password for the server pool (if needed).",
	pool_servers_desc[] =
		"Comma-separated list of server names or IP addresses\n\
		of machines in the server pool (default \"\" = none).",
	pool_tunejobs_desc[] =
		"Number of parallel tuning jobs (default 0) to run on the\n\
		server (if specified by pool_servers).  Tuning results are not\n\
		normalized by server performance, so tuning is most effective\n\
		when all the servers in the server pool have similar\n\
		performance characteristics.";
#endif /*}*/

 static char presolve_desc[] = "Whether to use Gurobi's presolve:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = conservative presolve\n\
			 2 = aggressive presolve.";

 static char pricing_desc[] = "Pricing strategy:\n\
			-1 = automatic choice (default)\n\
			 0 = partial pricing\n\
			 1 = steepest edge\n\
			 2 = Devex\n\
			 3 = quick-start steepest edge.";

#ifdef GRB_INT_ATTR_BRANCHPRIORITY
 static char priorities_desc[] =
		"Whether to use the variable.priority suffix with MIP problems.\n\
		When several branching candidates are available, a variable\n\
		with the highest .priority is chosen for the next branch.\n\
		Priorities are nonnegative integers (default 0).\n\
		Possible values for \"priorities\":\n\
			0 = ignore .priority; assume priority 0 for all vars\n\
			1 = use .priority if present (default).";
#endif

#if GRB_VERSION_MAJOR >= 4 /*{*/
 static char psdtol_desc[] = "Maximum diagonal perturbation to correct indefiniteness\n\
		in quadratic objectives (default 1e-6).";
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char pumppasses_desc[] = "Number of feasibility-pump passes to do after the MIP root\n\
		when no other root heuristoc found a feasible solution.\n\
		Default -1 = automatic choice.";
#endif /*}*/

#ifdef GRB_INT_PAR_QCPDUAL
 static char qcpdual_desc[] = "Whether to compute dual variables when the problem\n\
		has quadratic constraints (which can be expensive):\n\
			0 = no (default)\n\
			1 = yes.";
#endif

#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
 static char rays_desc[] = "Whether to return suffix .unbdd if the objective is unbounded\n\
		or suffix .dunbdd if the constraints are infeasible:\n\
			0 = neither\n\
			1 = just .unbdd\n\
			2 = just .dunbdd\n\
			3 = both (default).";
#endif /*}*/

 static char relax_desc[] = "Whether to relax integrality:\n\
			0 = no (default)\n\
			1 = yes: treat integer and binary variables\n\
				as continuous.";

#ifdef GRB_INT_PAR_RELAXLIFTCUTS
 static char relaxliftcuts_desc[] = "Whether to enable relax-and-lift cut generation:\n\
			-1 = automatic choice (default)\n\
			 0 = disallow relax-and-lift cuts\n\
			 1 = enable moderate relax-and-lift cut generation\n\
			 2 = enable aggressive relax-and-lift cut generation.\n\
		Values 1 and 2 override the \"cuts\" keyword.";
#endif

 static char return_mipgap_desc[] =
		"Whether to return mipgap suffixes or include mipgap values\n\
		(|objectve - best_bound|) in the solve_message:  sum of\n\
			1 = return relmipgap suffix (relative to |obj|);\n\
			2 = return absmipgap suffix (absolute mipgap);\n\
			4 = suppress mipgap values in solve_message.\n\
		Default = 0.  The suffixes are on the objective and problem.\n\
		Returned suffix values are +Infinity if no integer-feasible\n\
		solution has been found, in which case no mipgap values are\n\
		reported in the solve_message.";

#ifdef GRB_INT_PAR_RLTCUTS
 static char rltcuts_desc[] = "Whether to enable generation of cuts by the Relaxation\n\
		Linearization Technique (RLT):\n\
			-1 = automatic choice (default)\n\
			 0 = disallow RLT cuts\n\
			 1 = enable moderate RLT cut generation\n\
			 2 = enable aggressive RLT cut generation.\n\
		Values 1 and 2 override the \"cuts\" keyword.";
#endif

 static char round_desc[] =
		"Whether to round integer variables to integral values before\n\
		returning the solution, and whether to report that GUROBI\n\
		returned noninteger values for integer values:  sum of\n\
			 1 ==> round nonintegral integer variables\n\
			 2 ==> modify solve_result\n\
			 4 ==> modify solve_message\n\
		Default = 7.  Modifications that were or would be made are\n\
		reported in solve_result and solve_message only if the maximum\n\
		deviation from integrality exceeded round_reptol.";

 static char round_reptol_desc[] =
		"Tolerance for reporting rounding of integer variables to\n\
		integer values; see \"round\".  Default = 1e-9.";

#ifdef GRB_INT_PAR_SEED
 static char seed_desc[] =
		"Random number seed (default 0), affecting perturbations that\n\
		may influence the solution path.";
#endif

#ifdef GRB_INT_PAR_SIFTING /*{ new in 4.6 */
 static char sifting_desc[] =
		"Whether to use sifting within the dual simplex algorithm,\n\
		which can be useful when there are many more variables than\n\
		constraints:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes, moderate sifting\n\
			 2 = yes, aggressive sifting.";

  static char siftmethod_desc[] =
		"Algorithm to use for sifting with the dual simplex method:\n\
			-1 = automatic choice (default)\n\
			 0 = primal simplex\n\
			 1 = dual simplex\n\
			 2 = barrier.";
#endif /*}*/

 static char scale_desc[] = "Whether to scale the problem:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes"
#if GRB_VERSION_MAJOR >= 6
			"\n\
			 2 = yes, more aggressively"
#if GRB_VERSION_MAJOR >= 8
			"\n\
			 3 = yes, even more aggressively"
#endif
#endif
			"."
			;

 static char simplex_desc[] =
#if GRB_VERSION_MAJOR < 4
			"Which algorithm to use:"
#else
			"Which algorithm to use for non-MIP problems or for the root\n\
		node of MIP problems:"
#endif
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
			"\n\
			-1 automatic (default): 3 for LP, 2 for QP, 1 for MIP\n\
				root node\n\
			 0 = primal simplex\n\
			 1 = dual simplex"
#else
			"\n\
			0 = primal simplex\n\
			1 = dual simplex (default)"
#endif /*}*/
#if GRB_VERSION_MAJOR >= 3 /*{*/
			"\n\
			 2 = barrier"
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
			"\n\
			 3 = nondeterministic concurrent (several solves in\n\
				parallel)\n\
			 4 = deterministic concurrent"
#endif
#if GRB_VERSION_MAJOR >= 8
			"\n\
			 5 = deterministic concurrent simplex."
#endif
#endif /*}*/
			;

 static char solnsens_desc[] = "Whether to return suffixes for solution sensitivities, i.e.,\n\
		ranges of values for which the optimal basis remains optimal:\n\
			0 = no (default)\n\
			1 = yes:  suffixes return on variables are\n\
				.sensobjlo = smallest objective coefficient\n\
				.sensobjhi = greatest objective coefficient\n\
				.senslblo = smallest variable lower bound\n\
				.senslbhi = greatest variable lower bound\n\
				.sensublo = smallest variable upper bound\n\
				.sensubhi = greatest variable upper bound\n\
			suffixes for constraints are\n\
				.sensrhslo = smallest right-hand side value\n\
				.sensrhshi = greatest right-hand side value.\n\
		For range constraints Lconst <= expr <= Uconst or equivalently\n\
		Uconst >= expr >= Lconst (where both Lconst and Uconst are\n\
		constants and expr is an expression involving variables,\n\
		.sensrhslo and .sensrhshi apply to Lconst and\n\
				.sensrhslo2 = smallest Uconst\n\
				.sensrhslo2 = greatest Uconst.\n\
		Note that AMPL converts constraints with a single relational\n\
		operator into the form expr relop Const.  If you write\n\
		Const relop expr, AMPL converts it to -(expr) relop -(Const).\n\
		You need to take this into account when examining the\n\
		.sensrhslo and .sensrhshi values."
#if GRB_VERSION_MAJOR >= 5
			"\n\
		For problems with integer variables and quadratic constraints,\n\
		solnsens = 0 is quietly assumed."
#endif
;

 static char sos_desc[] = "whether to honor declared suffixes .sosno and .ref describing\n\
		SOS sets:\n\
			0 = no\n\
			1 = yes (default):  each distinct nonzero .sosno\n\
				value designates an SOS set, of type 1 for\n\
				positive .sosno values and of type 2 for\n\
				negative values.  The .ref suffix contains\n\
				corresponding reference values.";

 static char sos2_desc[] = "Whether to tell Gurobi about SOS2 constraints for nonconvex\n\
		piecewise-linear terms:\n\
			0 = no\n\
			1 = yes (default), using suffixes .sos and .sosref\n\
				provided by AMPL.";

#ifdef GRB_INT_PAR_STARTNODELIMIT
 static char startnodelimit_desc[] = "Limit on how many branch-and-bound nodes to explore when\n\
		doing a partial MIP start:\n\
			-2 = suppress MIP start processing\n\
			-1 = use submipnodes (default)\n\
			>= 0 ==> specific node limit.";
#endif
#ifdef GRB_INT_PAR_SUBMIPNODES
 static char submipnodes_desc[] = "Limit on nodes explored by MIP-based heuristics, e.g., RINS.\n\
		Default = 500.";
#endif

 static char threads_desc[] = "How many threads to use when using the barrier algorithm\n"
			"or solving MIP problems; default 0 ==> automatic choice.";

 static char timing_desc[] = "Whether to report timing:\n\
			0 = no (default)\n\
			1 = report times on stdout\n\
			2 = report times on stderr.";

#ifdef GRB_INT_PAR_TUNEOUTPUT /*{*/

 static char tunebase_desc[] = "Base name for results of running Gurobi's search for better\n\
		parameter settings.  The search is run only when tunebase\n\
		is specified.  Results are written to files with names derived\n\
		from tunebase by appending \".prm\" if \".prm\" does not occur in\n\
		tunebase and inserting 1, 2, ... (for the first, second,\n\
		... set of parameter settings) before the right-most \".prm\".\n\
		The file with \"1\" inserted is the best set and the solve\n\
		results returned are for this set.  In a subsequent \"solve;\",\n\
		you can use paramfile=... to apply the settings in results\n\
		file ... .";

 static char tuneoutput_desc[] = "Amount of tuning output when tunebase is specified:\n\
			0 = none\n\
			1 = summarize each new best parameter set\n\
			2 = summarize each set tried (default)\n\
			3 = summary plus detailed solver output for each trial.";

 static char tuneresults_desc[] = "Limit on the number of tuning result files to write\n\
		when tunerbase is specified.  The default (-1) is to write\n\
		results for all parameter sets on the efficient frontier.";

 static char tunetimelim_desc[] = "Time limit (in seconds) on tuning when tunebase\n\
		is specified.  Default -1 ==> automatic choice of time limit.";

 static char tunetrials_desc[] = "Number of trials for each parameter set when tunebase\n\
		is specified, each with a different random seed value.\n\
		Default = 3.";

#endif /*}*/

 static char varbranch_desc[] = "MIP branch variable selection strategy:\n\
			-1 = automatic choice (default)\n\
			 0 = pseudo reduced-cost branching\n\
			 1 = pseudo shadow-price branching\n\
			 2 = maximum infeasibility branching\n\
			 3 = strong branching.";

 static char version_desc[] = "Report version details before solving the problem.  This is a\n\
		single-word \"phrase\" that does not accept a value assignment.";
 static char work_desc[] = "Whether to report the amount of (deterministic) work spent on\n\
		the latest optimization:\n\
			0 = do not report (default)\n\
			1 = report in the problem suffix \"work\".";
 static char writeprob_desc[] = "Name of a GUROBI-format file to be written (for debugging);\n\
		must end in one of \".bas\", \".lp\", \".mps\", \".prm\", \".rew\",\n\
		\".rlp\", \".sol\", or for the \"fixed\" model used to recover a\n\
		basis or dual values for problems with integer variables or\n\
		quadratic constraints, \".fix_lp\" or \".fix_mps\"; the '_' will\n\
		be replaced by '.' in the name of the file written for\n\
		\".fix_lp\" or \".fix_mps\".  Can appear more than once with\n\
		different filenames.";

 /* WS_desc_ASL = modified solvers/ws_desc.c: extra initial tab; must not be static. */
 char WS_desc_ASL[] = "=... solution report without -AMPL: sum of\n\
			1 ==> write .sol file\n\
			2 ==> print primal variable values\n\
			4 ==> print dual variable values\n\
			8 ==> do not print solution message.";

#if GRB_VERSION_MAJOR > 1 /*{*/

 static char multprice_norm_desc[] = "Choice of norm used in multiple pricing:\n\
			-1 = automatic choice (default)\n\
			 0, 1, 2, 3 = specific choices:  hard to describe,\n\
				but sometimes a specific choice will perform\n\
				much better than the automatic choice.";

 static char nodefiledir_desc[] = "Directory where MIP tree nodes are written after memory\n\
		for them exceeds nodefilestart; default \".\"";

 static char nodefilestart_desc[] = "Gigabytes of memory to use for MIP tree nodes;\n\
		default = Infinity (no limit, i.e., no node files written).";

#ifdef GRB_INT_PAR_PREMIQPMETHOD /*{*/
 static char premiqpmethod_desc[] =
#ifdef GRB_INT_PAR_PREQLINEARIZE
			"Deprecated; replaced by preqlinearize.\n\
		Same possible values as preqlinearize.";
#else
			"How Gurobi's presolve should treat MIQP problems:\n\
			-1 = automatic choice (default)\n\
			 0 = leave the problem as an MIQP\n\
			 1 = try to transform an MIQP to an MILP.";
#endif
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char prepasses_desc[] = "Limit on the number of Gurobi presolve passes:\n\
			-1 = automatic choice (automatic)\n\
			 n >= 0: at most n passes.";
#endif /*}*/

#ifdef GRB_INT_PAR_PREQLINEARIZE /*{*/
 static char preqlinearize_desc[] = "How Gurobi's presolve should treat quadratic problems:\n\
			-1 = automatic choice (default)\n\
			 0 = do not modify the quadratic part(s)\n"
#if GRB_VERSION_MAJOR < 8 || (GRB_VERSION_MAJOR == 8 && GRB_VERSION_MINOR < 1)
			"\t\t\t 1 = try to linearize quadratic parts";
#else
			"\t\t\t 1 or 2 = try to linearize quadratic parts:\n\
				1 = focus on a strong LP relaxation\n\
				2 = focus on a compact relaxation.";
#endif
#endif /*}*/

#ifdef GRB_INT_PAR_PRESPARSIFY
 static char presparsify_desc[] =
		"Whether Gurobi's presolve should use its \"sparsify reduction\",\n\
		which sometimes gives significant problem-size reductions:\n\
			-1 = automatic choice\n\
			 0 = no\n\
			 1 = yes.";
#endif

 static char quad_desc[] = "Whether simplex should use quad-precision:\n\
			-1 = automatic choice (default)\n\
			 0 = no\n\
			 1 = yes.";

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char resultfile_desc[] = "Name of a file of extra information written after\n\
		completion of optimization.  The name's suffix determines what\n\
		is written:\n\
			.sol  solution vector\n\
			.bas  simplex basis\n\
			.mst  integer variable solution vector.";

 static char rins_desc[] = "How often to apply the RINS heuristic for MIP problems:\n\
			-1 = automatic choice (default)\n\
			 0 = never\n\
			 n > 0: every n-th node.";
#endif /*}*/

#if GRB_VERSION_MAJOR < 4 /*{*/
 static char rootmethod_desc[] = "Algorithm for MIP root relaxation:\n\
			0 = primal simplex\n\
			1 = dual simplex (default)"
#if GRB_VERSION_MAJOR >= 3
			"\n\
			2 = barrier."
#endif
			;
#endif /*}*/
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char symmetry_desc[] = "MIP symmetry detection:\n\
			-1 = automatic choice (default)\n\
			 0 = none\n\
			 1 = conservative\n\
			 2 = agressive.";
#endif /*}*/

#if GRB_VERSION_MAJOR >= 3 /*{*/
 static char warmstart_desc[] = "Whether to use incoming primal and dual variable values\n\
		(if both are available) in a simplex warm start:\n\
			0 = no;\n\
			1 = yes if there is no incoming basis (default);\n\
			2 = yes, ignoring the incoming basis (if any);\n"
#ifdef GRB_DBL_ATTR_VARHINTVAL
"			3 = no, but on MIP problems, use the incoming primal\n\
			    values as hints, ignoring the .hintpri suffix;\n\
			4 = similar to 3, but use the .hintpri suffix on\n\
			    variables:  larger (integer) values give greater\n\
			    priority to the initial value of the associated\n\
			    variable.\n"
#endif
"		Note that specifying basis=0 or basis=2 causes there to be\n\
		no incoming basis."
#ifdef GRB_DBL_ATTR_VARHINTVAL
		"  This is relevant to warmstart values\n\
		1, 3, and 4.  For continuous problems, warmstart values >= 2\n\
		are treated as 1."
#endif
#ifdef GRB_INT_PAR_LPWARMSTART
	 "\n\t\tNormally an incoming solution vector disables Gurobi's\n\
		LP presolve; to enable it set lpwarmstart to 2."
#endif
		;
#endif /*}*/

#ifdef GRB_DBL_PAR_WORKLIMIT
 static char worklimit_desc[] =
	 "Maximum amount of work expended (in work units); in contrast\n\
		to timelim, work limits are deterministic (default no limit).\n";
#endif

#ifdef GRB_INT_PAR_ZEROOBJNODES
 static char zeroobjnodes_desc[] =
		"Number of nodes to explore at the root MIP node if no other\n\
		heuristic has found a feasible solution.  Default = -1\n\
		(automatic choice).";
#endif

#define VP (void*)

#if GRB_VERSION_MAJOR >= 4
#define Method "Method"
#else
#define Method "LPMethod"
#endif

 static keyword
keywds[] = {	/* must be in alphabetical order */

#if GRB_VERSION_MAJOR >= 3
	{ "aggfill", sf_ipar, "AggFill", aggfill_desc },
#endif
	{ "aggregate", sf_ipar, "Aggregate", aggregate_desc },
#if GRB_VERSION_MAJOR >= 3 /*{*/
	{ "ams_eps", D_val, &ams_eps, ams_eps_desc },
	{ "ams_epsabs", D_val, &ams_epsabs, ams_epsabs_desc },
	{ "ams_limit", I_val, &ams_limit, ams_limit_desc },
#ifdef GRB_INT_PAR_POOLSEARCHMODE
	{ "ams_mode", Im_val, "PoolSearchMode", ams_mode_desc },
#endif
	{ "ams_stub", C_val, &ams_stub, ams_stub_desc },
	{ "barconvtol", sf_dpar, "BarConvTol", barconvtol_desc },
	{ "barcorrectors", sf_ipar, "BarCorrectors", barcorrectors_desc },
#ifdef GRB_INT_PAR_BARHOMOGENEOUS
	{ "barhomogeneous", sf_ipar, "BarHomogeneous", barhomogeneous_desc },
#endif
	{ "bariterlim",  sf_ipar, "BarIterLimit", bariterlimit_desc },
	{ "barorder", sf_ipar, "BarOrder", barorder_desc },
#endif /*}*/
#ifdef GRB_DBL_PAR_BARQCPCONVTOL
	{ "barqcptol", sf_dpar, "BarQCPConvTol", barqcptol_desc },
#endif
	{ "basis", sf_mint, VP set_basis, basis_desc },
	{ "basisdebug", sf_mint, VP set_basisdebug, basisdebug_desc },
#ifdef GRB_DBL_PAR_BESTBDSTOP
	{ "bestbndstop", sf_dpar, GRB_DBL_PAR_BESTBDSTOP, bestbndstop_desc },
#endif
	{ "bestbound", sf_mint, VP set_bestbound, bestbound_desc },
#ifdef GRB_DBL_PAR_BESTOBJSTOP
	{ "bestobjstop", sf_dpar, GRB_DBL_PAR_BESTOBJSTOP, bestobjstop_desc },
#endif
#ifdef GRB_INT_PAR_BQPCUTS
	{ "bqpcuts", sf_ipar, GRB_INT_PAR_BQPCUTS, bqpcuts_desc },
#endif
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	{ "branchdir", sf_ipar, "BranchDir", branchdir_desc },
#endif /*}*/
	{ "cliquecuts", sf_ipar, "CliqueCuts", overrides_cuts },
#ifdef GRB_STR_PAR_CLOUDHOST
	{ "cloudhost", C_val, &cloudhost, cloudhost_desc },
#endif
#if GRB_VERSION_MAJOR >= 7 && defined(ALLOW_GUROBI_SERVER) /*{*/
	{ "cloudid", C_val, &cloudid, cloudid_desc },
	{ "cloudkey", C_val, &cloudkey, cloudkey_desc },
	{ "cloudpool", C_val, &cloudpool, cloudpool_desc },
#if GRB_VERSION_MAJOR >= 8
	{ "cloudpriority", I_val, &cloudpriority, cloudpriority_desc },
#endif
#endif /*}*/
#ifdef GRB_INT_PAR_CONCURRENTMIP
	{ "concurrentmip", sf_ipar, GRB_INT_PAR_CONCURRENTMIP, concurrentmip_desc },
#endif
#ifdef GRB_INT_ATTR_CONCURRENTWINMETHOD
	{ "concurrentwin", sf_mint, VP set_concurrentwinmethod, concurrentwinmethod_desc},
#endif
	{ "covercuts", sf_ipar, "CoverCuts", overrides_cuts },
#if GRB_VERSION_MAJOR >= 3 /*{*/
	{ "crossover", sf_ipar, "Crossover", crossover_desc },
	{ "crossoverbasis", sf_ipar, "CrossoverBasis", crossoverbasis_desc },
#endif /*}*/
	{ "cutagg", sf_ipar, "CutAggPasses", cutagg_desc },
	{ "cutoff", sf_dpar, "Cutoff", cutoff_desc },
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	{ "cutpasses", sf_ipar, "CutPasses", cutpasses_desc },
#endif /*}*/
	{ "cuts", sf_ipar, "Cuts", cuts_desc },
#ifdef GRB_INT_PAR_DEGENMOVES
	{ "degenmoves", sf_ipar, GRB_INT_PAR_DEGENMOVES, degenmoves_desc },
#endif
#ifdef GRB_INT_PAR_DISCONNECTED
	{ "disconnected", sf_ipar, GRB_INT_PAR_DISCONNECTED, disconnected_desc },
#endif
#ifdef GRB_INT_PAR_DUALREDUCTIONS
	{ "dualreductions", sf_ipar, "DualReductions", dualreductions_desc },
#endif
#if GRB_VERSION_MAJOR >= 5 /*{*/
	{ "feasrelax", sf_mint, VP set_feasrelax, feasrelax_desc },
#ifdef GRB_DBL_PAR_FEASRELAXBIGM
	{ "feasrelaxbigm", sf_dpar, GRB_DBL_PAR_FEASRELAXBIGM, feasrelaxbigm_desc },
#endif
#endif /*}*/
	{ "feastol", sf_dpar, "FeasibilityTol", feastol_desc },
	{ "fixedmethod", I_val, &fixedmethod, fixedmethod_desc },
	{ "flowcover", sf_ipar, "FlowCoverCuts", "flowcover cuts:  " Overrides_cuts },
	{ "flowpath", sf_ipar, "FlowPathCuts", "flowpath cuts:  " Overrides_cuts },
	{ "gomory", sf_ipar, "GomoryPasses", gomory_desc },
	{ "gubcover", sf_ipar, "GUBCoverCuts", "gubcover cuts:  " Overrides_cuts },
	{ "heurfrac", sf_dpar, "Heuristics", heurfrac_desc },
	{ "iisfind", sf_mint, VP set_iis, iisfind_desc },
#if GRB_VERSION_MAJOR > 1 /*{*/
	{ "iismethod", sf_ipar, "IISMethod", iismethod_desc },
#endif /*}*/
	{ "implied", sf_ipar, "ImpliedCuts", "implied cuts:  " Overrides_cuts },
#ifdef GRB_DBL_PAR_IMPROVESTARTGAP
	{ "improvegap", sf_dpar, GRB_DBL_PAR_IMPROVESTARTGAP, isg_desc },
#endif
#ifdef GRB_DBL_PAR_IMPROVESTARTTIME
	{ "improvetime", sf_dpar, GRB_DBL_PAR_IMPROVESTARTTIME, ist_desc },
#endif
#ifdef GRB_DBL_PAR_IMPROVESTARTNODES
	{ "impstartnodes", sf_dpar, GRB_DBL_PAR_IMPROVESTARTNODES, isn_desc },
#endif
#ifdef GRB_INT_PAR_INFPROOFCUTS
	{ "infproofcuts", sf_ipar, GRB_INT_PAR_INFPROOFCUTS, infproofcuts_desc },
#endif
#ifdef GRB_INT_PAR_INTEGRALITYFOCUS
	{ "integrality", sf_ipar, GRB_INT_PAR_INTEGRALITYFOCUS, integralityfocus_desc },
#endif
	{ "intfeastol", sf_dpar, "IntFeasTol", intfeastol_desc },
	{ "intstart", sf_mint, VP set_intstart, intstart_desc },
	{ "iterlim", sf_dpar, "IterationLimit", "iteration limit (default: no limit)" },
#if GRB_VERSION_MAJOR >= 8
	{ "kappa", sf_mint, VP set_kappa, kappa_desc },
#endif
#if GRB_VERSION_MAJOR >= 6
	{ "lazy", sf_mint, VP set_lazy, lazy_desc },
#endif
#if GRB_VERSION_MAJOR >= 5
	{ "lbpen", D_val, &lbpen, "See feasrelax." },
#endif
#ifdef GRB_INT_PAR_LIFTPROJECTCUTS
	{ "liftprojectcuts", sf_ipar, GRB_INT_PAR_LIFTPROJECTCUTS, liftprojectcuts_desc },
#endif
	{ "logfile", C_val, &logfile, logfile_desc },
	{ "logfreq", sf_iparlog, "DisplayInterval", logfreq_desc },
	{ "lpmethod", sf_ipar, Method, "Synonym for method." },
#ifdef GRB_INT_PAR_LPWARMSTART
	{ "lpwarmstart", sf_ipar, GRB_INT_PAR_LPWARMSTART, lpwarmstart_desc },
#endif
	{ "maxmipsub", sf_ipar, "SubMIPNodes", maxmipsub_desc },
#ifdef GRB_DBL_ATTR_MAX_VIO
	{ "maxvio", sf_mint, VP set_maxvio, maxvio_desc},
#endif
#ifdef GRB_DBL_PAR_MEMLIMIT
	{ "memlimit", sf_dpar, GRB_DBL_PAR_MEMLIMIT, memlimit_desc},
#endif
	{ "method", sf_ipar, Method, simplex_desc },
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	{ "minrelnodes", sf_ipar, "MinRelNodes", minrelnodes_desc },
#endif /*}*/
#if GRB_VERSION_MAJOR >= 3
	{ "mipfocus", sf_ipar, "MIPFocus", mipfocus_desc },
#endif
	{ "mipgap", sf_dpar, "MipGap", "max. relative MIP optimality gap (default 1e-4)" },
#if GRB_VERSION_MAJOR >= 3
	{ "mipgapabs", sf_dpar, "MipGapAbs", "absolute MIP optimality gap (default 1e-10)" },
#endif
	{ "mipsep", sf_ipar, "MIPSepCuts", "MIPsep cuts:  " Overrides_cuts },
	{ "mipstart", sf_mint, VP set_mipstval, mipstart_desc },
#ifdef GRB_INT_PAR_MIQCPMETHOD
	{ "miqcpmethod", sf_ipar, GRB_INT_PAR_MIQCPMETHOD, miqcpmethod_desc },
#endif
	{ "mircuts", sf_ipar, "MIRCuts", "MIR cuts:  " Overrides_cuts },
#if GRB_VERSION_MAJOR >= 4
	{ "modkcuts", sf_ipar, "ModKCuts", "mod-k cuts:  " Overrides_cuts },
#endif
#if GRB_VERSION_MAJOR >= 7
	{ "multiobj", sf_mint, VP set_multiobj, multiobj_desc },
	{ "multiobjmethod", sf_ipar, "MultiObjMethod", multiobjmethod_desc },
	{ "multiobjpre", sf_ipar, "MultiObjPre", multiobjpre_desc },
#endif
#if GRB_VERSION_MAJOR > 1 /*{*/
	{"multprice_norm", sf_ipar, "NormAdjust", multprice_norm_desc},
#ifdef GRB_INT_PAR_NETWORKALG
	{ "networkalg", sf_ipar, GRB_INT_PAR_NETWORKALG, networkalg_desc},
#endif
#if GRB_VERSION_MAJOR >= 3
	{ "networkcuts", sf_ipar, "NetworkCuts", "Network cuts:  " Overrides_cuts },
#endif
#ifdef GRB_INT_PAR_NLPHEUR
	{ "nlpheur", sf_ipar, GRB_INT_PAR_NLPHEUR, nlpheur_desc },
#endif
	{"nodefiledir", sf_spar, "NodefileDir", nodefiledir_desc},
	{"nodefilestart", sf_dpar, "NodefileStart", nodefilestart_desc},
#endif /*}*/
	{ "nodelim", sf_dpar, "NodeLimit", "maximum MIP nodes to explore (default: no limit)" },
#if GRB_VERSION_MAJOR >= 4 /*{*/
	{ "nodemethod", sf_ipar, "NodeMethod", nodemethod_desc },
#endif /*}*/
#ifdef GRB_INT_PAR_NONCONVEX
	{ "nonconvex", sf_ipar, GRB_INT_PAR_NONCONVEX, nonconvex_desc },
#endif
#ifdef GRB_DBL_PAR_NORELHEURTIME
	{"norelheurtime",  sf_dpar, GRB_DBL_PAR_NORELHEURTIME, norelheurtime_desc },
	{"norelheurwork",  sf_dpar, GRB_DBL_PAR_NORELHEURWORK, norelheurwork_desc },
#endif
#if GRB_VERSION_MAJOR > 1 /*{*/
	{ "normadjust", sf_ipar, "NormAdjust", "Synonym for multprice_norm." },
#endif /*}*/
#ifdef GRB_INT_PAR_NUMERICFOCUS
	{ "numericfocus", sf_ipar, GRB_INT_PAR_NUMERICFOCUS, numericfocus_desc },
#endif
#ifdef GRB_INT_PAR_OBBT
	{ "obbt", sf_ipar, GRB_INT_PAR_OBBT, obbt_desc },
#endif
	{ "objno", sf_mint, VP set_objno, objno_desc },
#if GRB_VERSION_MAJOR >= 4
	{ "objrep", sf_mint, VP set_objrep, objrep_desc },
#endif
#if GRB_VERSION_MAJOR > 1 /*{*/
	{ "objscale", sf_dpar, "ObjScale", objscale_desc },
#endif /*}*/
	{ "opttol", sf_dpar, "OptimalityTol", opttol_desc },
	{ "outlev", sf_mint, VP set_outlev, outlev_desc },
#if GRB_VERSION_MAJOR >= 3
	{ "param", sf_spar, "Dummy", param_desc },
	{ "paramfile", sf_pf, 0, paramfile_desc },
#endif
#ifdef GRB_INT_PAR_PARTITIONPLACE
	{ "partitionplace", sf_ipar, GRB_INT_PAR_PARTITIONPLACE, partitionplace_desc },
#endif
	{ "perturb", sf_dpar, "PerturbValue", perturb_desc },
	{ "pivtol", sf_dpar, "MarkowitzTol", "Markowitz pivot tolerance (default 7.8125e-3)" },
	/*GRB_MAJ2(({"precrush", sf_ipar, "PreCrush", precrush_desc},))*/
#ifdef GRB_STR_ATTR_SERVER /*{*/
#if GRB_VERSION_MAJOR >= 6
	{ "pl_bigm", D_val, &pl_bigm, pl_bigm_desc },
	{ "pool_distmip", sf_ipar, GRB_INT_PAR_DISTRIBUTEDMIPJOBS, pool_distmip_desc },
	{ "pool_mip", sf_ipar, GRB_INT_PAR_CONCURRENTJOBS, pool_mip_desc },
	{ "pool_password", sf_spar, GRB_STR_PAR_WORKERPASSWORD, pool_password_desc },
	{ "pool_servers", sf_spar, GRB_STR_PAR_WORKERPOOL, pool_servers_desc },
#else
	{ "pool_mip", sf_ipar, GRB_INT_PAR_CONCURRENTMIPJOBS, pool_mip_desc },
	{ "pool_password", sf_spar, GRB_STR_PAR_SERVERPASSWORD, pool_password_desc },
	{ "pool_servers", sf_spar, GRB_STR_PAR_SERVERPOOL, pool_servers_desc },
#endif
	{ "pool_tunejobs", sf_ipar, GRB_INT_PAR_TUNEJOBS, pool_tunejobs_desc },
#endif /*}*/
#if GRB_VERSION_MAJOR >= 3
	{ "poolgap", D_val, &ams_eps, "Synonym for ams_eps." },
	{ "poolgapabs", D_val, &ams_epsabs, "Synonym for ams_epsabs." },
#endif
#ifdef GRB_INT_PAR_POOLSEARCHMODE
	{ "poolsearchmode", sf_ipar, "PoolSearchMode", "Synonym for ams_mode." },
#endif
#ifdef GRB_INT_PAR_POOLSOLUTIONS
	{ "poolsolutions", I_val, &ams_limit, "Synonym for ams_limit." },
#endif
#if GRB_VERSION_MAJOR >= 3
	{ "poolstub", C_val, &ams_stub, "Synonym for ams_stub." },
	{ "predeprow", sf_ipar, "PreDepRow", predeprow_desc },
	{ "predual", sf_ipar, "PreDual", predual_desc },
#ifdef GRB_INT_PAR_PREMIQCPFORM
	{ "premiqcpform", sf_ipar, GRB_INT_PAR_PREMIQCPFORM, premiqcpform_desc },
#endif
#ifdef GRB_INT_PAR_PREMIQPMETHOD
	{ "premiqpmethod", sf_ipar, "PreMIQPMethod", premiqpmethod_desc },
#endif
	{ "prepases", sf_ipar, "PrePasses", "deprecated synonym for \"prepasses\"" },
	{ "prepasses", sf_ipar, "PrePasses", prepasses_desc },
#endif
#ifdef GRB_INT_PAR_PREQLINEARIZE
	{ "preqlinearize", sf_ipar, GRB_INT_PAR_PREQLINEARIZE, preqlinearize_desc },
#endif
	{ "presolve", sf_ipar, "Presolve", presolve_desc },
#ifdef GRB_DBL_PAR_PRESOS1BIGM
	{ "presos1bigm", sf_dpar, GRB_DBL_PAR_PRESOS1BIGM, presos1bigm_desc },
	{ "presos2bigm", sf_dpar, GRB_DBL_PAR_PRESOS2BIGM, presos2bigm_desc },
#endif
#ifdef GRB_INT_PAR_PRESOS1ENCODING /* new in 9.5 */
	{ "presos1enc", sf_ipar, GRB_INT_PAR_PRESOS1ENCODING, presos1enc_desc },
	{ "presos2enc", sf_ipar, GRB_INT_PAR_PRESOS2ENCODING, presos2enc_desc },
#endif
#ifdef GRB_INT_PAR_PRESPARSIFY /* new in 4.6 */
	{ "presparsify", sf_ipar, GRB_INT_PAR_PRESPARSIFY, presparsify_desc },
#endif
	{ "pricing", sf_ipar, "SimplexPricing", pricing_desc },
#ifdef GRB_INT_ATTR_BRANCHPRIORITY
	{ "priorities", sf_mint, VP set_priorities, priorities_desc },
#endif
#if GRB_VERSION_MAJOR >= 4
	{ "psdtol", sf_dpar, "PSDTol", psdtol_desc },
#endif
#if GRB_VERSION_MAJOR >= 3
	{ "pumppasses", sf_ipar, "PumpPasses", pumppasses_desc },
#endif
#ifdef GRB_INT_PAR_QCPDUAL
	{ "qcpdual", sf_ipar, "QCPDual", qcpdual_desc },
#endif
#if GRB_VERSION_MAJOR > 1 /*{*/
	{ "quad", sf_ipar, "Quad", quad_desc},
#endif /*}*/
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	{ "rays", sf_mint, VP set_rays, rays_desc },
#endif /*}*/
	{ "relax", sf_mint, VP set_relax, relax_desc },
#ifdef GRB_INT_PAR_RELAXLIFTCUTS
	{ "relaxliftcuts", sf_ipar, GRB_INT_PAR_RELAXLIFTCUTS, relaxliftcuts_desc },
#endif
#if GRB_VERSION_MAJOR >= 3
	{ "resultfile", sf_spar, "ResultFile", resultfile_desc },
#endif
	{ "return_mipgap", sf_mint, VP set_retmipgap, return_mipgap_desc },
#if GRB_VERSION_MAJOR >= 5
	{ "rhspen", D_val, &rhspen, "See feasrelax." },
#endif
#if GRB_VERSION_MAJOR >= 3
	{ "rins", sf_ipar, "RINS", rins_desc },
#endif
#ifdef GRB_INT_PAR_RLTCUTS
	{ "rltcuts", sf_ipar, GRB_INT_PAR_RLTCUTS, rltcuts_desc },
#endif
#if GRB_VERSION_MAJOR > 1 && GRB_VERSION_MAJOR < 4 /*{*/
	{"rootmethod", sf_ipar, "RootMethod", rootmethod_desc},
#endif /*}*/
	{ "round", sf_mint, VP set_round, round_desc },
	{ "round_reptol", D_val, &round_reptol, round_reptol_desc },
	{ "scale", sf_ipar, "ScaleFlag", scale_desc },
#ifdef GRB_INT_PAR_SEED
	{ "seed", sf_ipar, GRB_INT_PAR_SEED, seed_desc },
#endif
#ifdef ALLOW_GUROBI_SERVER
	{ "server", C_val, &server, server_desc },
#if GRB_VERSION_MAJOR >= 8
	{ "server_insecure", I_val, &server_insecure, server_insecure_desc },
#endif
	{ "server_password", C_val, &server_passwd, server_passwd_desc },
#if GRB_VERSION_MAJOR < 8
	{ "server_port", I_val, &server_port, server_port_desc },
#endif
	{ "server_priority", I_val, &server_priority, server_priority_desc },
#if GRB_VERSION_MAJOR > 9 || (GRB_VERSION_MAJOR == 9 && GRB_VERSION_MINOR >= 5)
	{ "server_timeout", I_val, &server_timeout, server_timeout_desc },
#else
	{ "server_timeout", D_val, &server_timeout, server_timeout_desc },
#endif
	{ "serverlic", C_val, &serverlic, serverlic_desc },
	{ "servers",  C_val, &server, "Synonym for server." },
#endif
#ifdef GRB_INT_PAR_SIFTING /* new in 4.6 */
	{ "sifting", sf_ipar, GRB_INT_PAR_SIFTING, sifting_desc },
	{ "siftmethod", sf_ipar, GRB_INT_PAR_SIFTMETHOD, siftmethod_desc },
#endif
	{ "simplex", sf_ipar, Method, "Synonym for method." },
#ifdef GRB_DBL_PAR_SOFTMEMLIMIT
	{ "softmemlimit", sf_dpar, GRB_DBL_PAR_SOFTMEMLIMIT, softmemlimit_desc },
#endif
	{ "solnlimit", sf_ipar, "SolutionLimit", "maximum MIP solutions to find (default 2e9)" },
	{ "solnsens", sf_mint, VP set_solnsens, solnsens_desc },
	{ "sos", sf_mint, VP set_sos, sos_desc },
	{ "sos2", sf_mint, VP set_sos2, sos2_desc },
#ifdef GRB_INT_PAR_STARTNODELIMIT
	{ "startnodelimit", sf_ipar, GRB_INT_PAR_STARTNODELIMIT, startnodelimit_desc },
#endif
#if GRB_VERSION_MAJOR >= 3
	{ "submipcuts", sf_ipar, "SubMIPCuts", "sub-MIP cuts:  " Overrides_cuts },
#ifdef GRB_INT_PAR_SUBMIPNODES
	{ "submipnodes", sf_ipar, GRB_INT_PAR_SUBMIPNODES, submipnodes_desc },
#endif
	{ "symmetry", sf_ipar, "Symmetry", symmetry_desc },
#endif
	{ "threads", sf_ipar, "Threads", threads_desc },
	{ "timelim", sf_dpar, "TimeLimit", "limit on solve time (in seconds; default: no limit)" },
	{ "timing", sf_mint, VP set_timing, timing_desc },
#ifdef GRB_INT_PAR_TUNEOUTPUT
	{ "tunebase", C_val, &tunebase, tunebase_desc },
	{ "tuneoutput", sf_ipar, GRB_INT_PAR_TUNEOUTPUT, tuneoutput_desc },
	{ "tuneresults", sf_ipar, GRB_INT_PAR_TUNERESULTS, tuneresults_desc },
	{ "tunetimelimit", sf_dpar, GRB_DBL_PAR_TUNETIMELIMIT, tunetimelim_desc },
	{ "tunetrials", sf_ipar, GRB_INT_PAR_TUNETRIALS, tunetrials_desc },
#endif
#if GRB_VERSION_MAJOR >= 5
	{ "ubpen", D_val, &ubpen, "See feasrelax." },
#endif
	{ "varbranch", sf_ipar, "VarBranch", varbranch_desc },
	{ "version", Ver_val, 0, version_desc },
	{ "wantsol", WS_val, 0, WS_desc_ASL+5 },
#if GRB_VERSION_MAJOR >= 5
	{ "warmstart", sf_mint, VP set_warmstart, warmstart_desc },
#endif
#ifdef GRB_DBL_ATTR_WORK
	{ "work", sf_mint, VP set_work, work_desc},
#endif
#ifdef GRB_DBL_PAR_WORKLIMIT
	{ "worklimit", sf_dpar,GRB_DBL_PAR_WORKLIMIT, worklimit_desc },
#endif
	{ "writeprob", sf_wfile, 0, writeprob_desc }
#if GRB_VERSION_MAJOR > 1 /*{*/
	,{"zerohalfcuts", sf_ipar, "ZeroHalfCuts", "zero-half cuts:  " Overrides_cuts }
#endif /*}*/
#ifdef GRB_INT_PAR_ZEROOBJNODES /* new in 4.6 */
	,{"zeroobjnodes", sf_ipar, GRB_INT_PAR_ZEROOBJNODES, zeroobjnodes_desc }
#endif
	};

#ifdef NO_MOkwf
#define MOkwf 0
#else
 static fint MOkwf(char *key, fint klen);
#endif

 static Option_Info
Oinfo = { "gurobi", verbuf, "gurobi_options", keywds, nkeywds, ASL_OI_keep_underscores, verbuf,
	  0, MOkwf,0,0,0, 20230310, 0,0,0,0,0,0,0, ASL_OI_tabexpand | ASL_OI_addnewline };

 static void
enamefailed(GRBenv *env, const char *what, const char *name)
{
	fprintf(Stderr, "%s(\"%s\") failed:\n\t%s.\n", what, name, GRBgeterrormsg(env));
	++Oinfo.n_badopts;
	}

#if 0
 static char iisforce_table[] = "\n\
0	auto (default)\n\
1	force out\n\
2	force in\n";
#endif

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

 static SufDecl
suftab[] = {
	{ "absmipgap", 0, ASL_Sufkind_obj   | ASL_Sufkind_outonly },
	{ "absmipgap", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
	{ "bestbound", 0, ASL_Sufkind_obj   | ASL_Sufkind_outonly },
	{ "bestbound", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
#ifdef GRB_INT_ATTR_CONCURRENTWINMETHOD
	{ "concurrentwinmethod", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly },
#endif
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5
	{ "dunbdd", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },
#endif
#ifdef GRB_INT_ATTR_VARHINTPRI
	{ "hintpri", 0, ASL_Sufkind_var },
#endif
	{ "iis", iis_table, ASL_Sufkind_var | ASL_Sufkind_outonly },
	{ "iis", 0, ASL_Sufkind_con | ASL_Sufkind_outonly },

#ifdef GRB_INT_ATTR_IIS_LBFORCE
	{ "iisforce", 0, ASL_Sufkind_con | ASL_Sufkind_input},
	{ "iisforcelb", 0, ASL_Sufkind_var | ASL_Sufkind_input},
	{ "iisforceub", 0, ASL_Sufkind_var | ASL_Sufkind_input},
#endif
#if GRB_VERSION_MAJOR >= 8
	{ "kappa", 0, ASL_Sufkind_obj | ASL_Sufkind_outonly | ASL_Sufkind_real},
	{ "kappa", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly | ASL_Sufkind_real},
#endif
#if GRB_VERSION_MAJOR >= 6
	{ "lazy", 0, ASL_Sufkind_con },
#endif
#if GRB_VERSION_MAJOR >= 5
	{ "lbpen", 0, ASL_Sufkind_var | ASL_Sufkind_real },
#endif
#ifdef GRB_DBL_ATTR_MAX_VIO
	{ "maxvio", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly },
#endif
#if GRB_VERSION_MAJOR >= 3
	{ "npool", 0, ASL_Sufkind_obj   | ASL_Sufkind_outonly },
	{ "npool", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
#endif
#ifdef GRB_INT_PAR_POOLSEARCHMODE
	{"poolignore", 0, ASL_Sufkind_var | ASL_Sufkind_input},
#endif
#if GRB_VERSION_MAJOR >= 7
	{ "objabstol", 0, ASL_Sufkind_obj | ASL_Sufkind_real },
	{ "objpriority", 0, ASL_Sufkind_obj },
	{ "objreltol", 0, ASL_Sufkind_obj | ASL_Sufkind_real },
	{ "objweight", 0, ASL_Sufkind_obj | ASL_Sufkind_real },
#endif
#if GRB_VERSION_MAJOR >= 8
	{ "partition", 0, ASL_Sufkind_var },
#endif
#ifdef GRB_INT_ATTR_BRANCHPRIORITY
	{ "priority", 0, ASL_Sufkind_var },
#endif
	{ "ref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "relmipgap", 0, ASL_Sufkind_obj   | ASL_Sufkind_outonly },
	{ "relmipgap", 0, ASL_Sufkind_prob  | ASL_Sufkind_outonly },
#if GRB_VERSION_MAJOR >= 5
	{ "rhspen", 0, ASL_Sufkind_con | ASL_Sufkind_real },
#endif
	{ "senslbhi",  0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "senslblo",  0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensobjhi", 0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensobjlo", 0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensrhshi", 0, ASL_Sufkind_con | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensrhshi2", 0, ASL_Sufkind_con | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensrhslo", 0, ASL_Sufkind_con | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensrhslo2", 0, ASL_Sufkind_con | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensubhi",  0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sensublo",  0, ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_outonly },
	{ "sos", 0, ASL_Sufkind_var },
	{ "sos", 0, ASL_Sufkind_con },
	{ "sosno", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "sosref", 0, ASL_Sufkind_var | ASL_Sufkind_real },
	{ "sstatus", 0, ASL_Sufkind_var, 1 },
	{ "sstatus", 0, ASL_Sufkind_con, 1 }
#if GRB_VERSION_MAJOR >= 5
	,{ "ubpen", 0, ASL_Sufkind_var | ASL_Sufkind_real }
#endif
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	,{ "unbdd", 0, ASL_Sufkind_var | ASL_Sufkind_outonly }
#endif
#ifdef GRB_DBL_ATTR_WORK
	,{ "work", 0, ASL_Sufkind_prob | ASL_Sufkind_outonly }
#endif
	};

#ifdef GRB_INT_ATTR_BRANCHPRIORITY
 static void
add_priorities(ASL *asl, GRBenv *env, GRBmodel *mdl)
{
	SufDesc *dp;
	int *p;

	if ((dp = suf_get("priority", ASL_Sufkind_var))
	 && (p = dp->u.i)
	 && GRBsetintattrarray(mdl, GRB_INT_ATTR_BRANCHPRIORITY, 0, n_var, p))
		failed(env, "GRBsetintattrarray(\"BranchPriority\")");
	}
#endif

#ifdef GRB_INT_ATTR_IIS_LBFORCE
 static void
force_iis(ASL* asl, GRBenv* env, GRBmodel* mdl)
{
	int i, j, *z;
	int map[3] = { -1, 0, 1 };
	SufDesc* dp;
	int m, n, nq;

	n = n_var;
	m = n_con;
	if((dp = suf_get("iisforce", ASL_Sufkind_con)) && (z = dp->u.i)) {
		for (i = 0; i < m; i++)
			z[i] = map[(j = z[i]) >= 0 && j <= 2 ? j : 0];
		nq = nlc;
		if (nq)
			GRBsetintattrarray(mdl, GRB_INT_ATTR_IIS_QCONSTRFORCE, 0, nq, z);
		if (m - nq)
			GRBsetintattrarray(mdl, GRB_INT_ATTR_IIS_CONSTRFORCE, 0, m - nq, z + nq);
		}
	if ((dp = suf_get("iisforcelb", ASL_Sufkind_var)) && (z = dp->u.i)) {
		for (i = 0; i < n; i++)
			z[i] = map[(j = z[i]) >= 0 && j <= 2 ? j : 0];
		GRBsetintattrarray(mdl, GRB_INT_ATTR_IIS_LBFORCE, 0, n, dp->u.i);
		}
	if ((dp = suf_get("iisforceub", ASL_Sufkind_var)) && (z = dp->u.i)) {
		z = dp->u.i;
		for (i = 0; i < n; i++)
			z[i] = map[(j = z[i]) >= 0 && j <= 2 ? j : 0];
		GRBsetintattrarray(mdl, GRB_INT_ATTR_IIS_UBFORCE, 0, n, dp->u.i);
		}
	}
#endif

 static void
show_times(void)
{
	FILE *f;
	int i;

	Times[3] = xectim_();
	Times[4] = time(0) - Times[4];
	for(i = 1; i <= 2; i++)
	    if (time_flag & i) {
		f = i == 1 ? stdout : Stderr;
		fprintf(f, "\nTimes (seconds):\nInput =  %g"
			"\nSolve =  %g (summed over threads)"
			"\nOutput = %g\nElapsed ",
			Times[1] - Times[0], Times[2] - Times[1],
			Times[3] - Times[2]);
		fprintf(f, Times[4] < 1. ? "< 1\n" : "= %g\n", Times[4]);
		}
	}

 Sig_ret_type
intcatch(int n)
{
	printf("\n<BREAK> (gurobi)\n", n);
	fflush(stdout);
	if (++breaking > 3)
		longjmp(Jb, 2);
	signal(SIGINT, intcatch);
	if (grbmodel)
		GRBterminate(grbmodel);
	SigRet;
	}

 static const char*
retiis(ASL *asl, GRBenv *env, GRBmodel *mdl, Dims *d, const char *what, int *srp)
{
	SufDesc *sd;
	char buf[128], *rv;
	int *c, i, j, k, kv, m, n, nr, *s, *v;

	m = n_con;
	n = n_var;
	nr = n + nranges;
	s = d->rstat;
	if (GRBgetintattrarray(mdl, "IISConstr", 0, m-nlc, s+nlc))
		failed(env, "GRBgetintattrarray(\"IISConstr\")");
	c = v = 0;
	for(i = k = kv = 0; i < m; ++i) {
		if (s[i]) {
			c = (int*)M1zapalloc(m*sizeof(int));
			for(; i < m; ++i)
				if (s[i]) {
					c[i] = 4;
					++k;
					}
			break;
			}
		}
	if (GRBgetintattrarray(mdl, "IISLB", 0, nr, s = d->cstat))
		failed(env, "GRBgetintattrarray(\"IISLB\")");
	for(i = 0; i < n; ++i) {
		if (s[i]) {
			v = (int*)M1zapalloc(n*sizeof(int));
			for(; i < n; ++i)
				if (s[i]) {
					v[i] = 1;
					++kv;
					}
			break;
			}
		}
	for(i = n; i < nr; ++i) {
		if (s[i]) {
			if (!c)
				c = (int*)M1zapalloc(m*sizeof(int));
			for(; i < nr; ++i)
				if (s[i] && !c[i]) {
					c[i] = 4;
					++k;
					}
			break;
			}
		}
	if (GRBgetintattrarray(mdl, "IISUB", 0, nr, s = d->cstat))
		failed(env, "GRBgetintattrarray(\"IISUB\")");
	for(i = 0; i < n; ++i) {
		if (s[i]) {
			if (!v)
				v = (int*)M1zapalloc(n*sizeof(int));
			for(; i < n; ++i)
				if (s[i] && !v[i]) {
					v[i] = 3;
					++kv;
					}
			break;
			}
		}
	for(i = n; i < nr; ++i) {
		if (s[i]) {
			if (!c)
				c = (int*)M1zapalloc(m*sizeof(int));
			for(; i < nr; ++i)
				if (s[i] && !c[i]) {
					c[i] = 4;
					++k;
					}
			break;
			}
		}
	sd = 0;
	if (c)
		sd = suf_iput("iis", ASL_Sufkind_con, c);
	if (v)
		suf_iput("iis", ASL_Sufkind_var, v);
	else if (sd)
		sd->table = iis_table;
	*srp = 201;
	if (k) {
		j = Snprintf(buf, sizeof(buf), "%s\nReturning an IIS of %d constraints",
			what, k);
		j += Snprintf(buf+j, sizeof(buf)-j, kv ? " and %d variables." : ".", kv);
		}
	else if (kv)
		j = Snprintf(buf, sizeof(buf), "%s\nReturning an IIS of %d variables.",
			what, kv);
	else {
		j = Snprintf(buf, sizeof(buf), "%s; empty IIS!");
		*srp = 202;
		}
	rv = (char*)M1alloc(++j);
	memcpy(rv, buf, j);
	return rv;
	}

 static void
dpf(Dims  *d, const char *fmt, ...)
{
	va_list ap;

	if (d->mbend > d->mb) {
		va_start(ap, fmt);
		d->mb += Vsnprintf(d->mb, d->mbend - d->mb, fmt, ap);
		va_end(ap);
		}
	}

 static void
missing_msg(Dims *d)
{
	static const char *missing[3] = {
		"primal",
		"dual",
		"primal or dual"
		};
	dpf(d, "\nNo %s variables returned.", missing[d->missing - 1]);
	}

#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
#ifdef RAYDEBUG
#define Debug(x) x
#else
#define Debug(x) /*nothing*/
#endif

 static int
send_ray(ASL *asl, Dims *d, GRBenv *env, GRBmodel *mdl)
{
	int i, n, n0;
	real *c, t, *y;

	n = n_var;
	n0 = d->nv0;
	y = (real *)M1zapalloc(n0 * 2*sizeof(real));
	c = y + n;
	if (n > n0)
		n = n0;
	if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_UNBDRAY, 0, n, y)) {
		Debug(printf("Get UnbdRay failed: %s\n", GRBgeterrormsg(env)));
		return 0;
		}
	if (!GRBgetdblattrarray(mdl, GRB_DBL_ATTR_OBJ, 0, n, c)) {
		t = 0.;
		for(i = 0; i < n; ++i)
			t += c[i]*y[i];
		if (d->objsense * t > 0.)
			for(i = 0; i < n; ++i)
				y[i] = -y[i];
		}
	suf_rput("unbdd", ASL_Sufkind_var, y);
	return 1;
	}

 static int
send_dray(ASL *asl, Dims *d, GRBenv *env, GRBmodel *mdl)
{
	char *sense;
	int i, n, n0;
	real *rhs, t, t1, *y;

	n = n_con;
	n0 = d->nc0;
	y = (real *)M1zapalloc(n*sizeof(real) + n0*(sizeof(real) + 1));
	rhs = y + n;
	sense = (char*)(rhs + n0);
	if (n > n0)
		n = n0;
	if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_FARKASDUAL, 0, n, y)) {
		Debug(printf("Get FarkasDual failed: %s\n", GRBgeterrormsg(env)));
		return 0;
		}
	if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_RHS, 0, n, rhs))
		Debug(printf("Get RHS failed: %s\n", GRBgeterrormsg(env)));
	else if (GRBgetcharattrarray(mdl, GRB_CHAR_ATTR_SENSE, 0, n, sense))
		Debug(printf("Get SENSE failed: %s\n", GRBgeterrormsg(env)));
	else {
		t = 0.;
		for(i = 0; i < n; ++i) {
			t1 = y[i]*rhs[i];
			if (sense[i] == '>')
				t1 = -t1;
			t += t1;
			}
		if (d->objsense * t < 0.)
			for (i = 0; i < n; ++i)
				y[i] = -y[i];
		}
	suf_rput("dunbdd", ASL_Sufkind_con, y);
	return 1;
	}

#undef Debug
#endif /*}*/

#if GRB_VERSION_MAJOR >= 5 /*{*/
 static void
pen_set(real *x, int n, SufDesc *d, real pen)
{
	real *s, t, *xe;

	if (pen <= 0.)
		pen = GRB_INFINITY;
	xe = x + n;
	if (d && (s = d->u.r)) {
		while(x < xe)
			*x++ = (t = *s++) > 0. ? t : pen;
		return;
		}
	while(x < xe)
		*x++ = pen;
	}

 static real *
do_feasrelax(ASL *asl, GRBenv *env, GRBmodel *mdl, const char **objqual, real *fto)
{
	SufDesc *dlb, *drhs, *dub;
	int mr, nc, nlb, nr, nrhs, nub, nv, nvr, rt;
	real pen, t;
	real *lbp, *lbr, *lu, *lue, *p, *rhp, *rv, *ubp, *ubr, *x;
	size_t L;

	nlb = lbpen > 0.;
	nub = ubpen > 0.;
	nrhs = rhspen > 0.;
	if ((dlb = suf_get("lbpen", ASL_Sufkind_var | ASL_Sufkind_input)))
		++nlb;
	if ((drhs = suf_get("rhspen", ASL_Sufkind_con | ASL_Sufkind_input)))
		++nrhs;
	if ((dub = suf_get("ubpen", ASL_Sufkind_var | ASL_Sufkind_input)))
		++nub;
	if (!(nlb + nub + nrhs))
		return 0;
	nc = n_con;
	nv = n_var;
	if ((nr = nranges) && nrhs) {
		++nlb;
		++nub;
		}
	nvr = nv + nr;
	L = 0;
	if (nlb)
		L += nvr;
	if (nub)
		L += nvr;
	if (nrhs)
		L += nc;
	lbp = rhp = ubp = 0;
	x = (real*)Malloc(L*sizeof(real));
	if (nlb) {
		pen_set(lbp = x, nv, dlb, lbpen);
		x += nvr;
		}
	if (nub) {
		pen_set(ubp = x, nv, dub, ubpen);
		x += nvr;
		}
	if (nrhs) {
		pen_set(rhp = x, nc, drhs, rhspen);
		if (nr) {
			lbr = lbp + nv;
			ubr = ubp + nv;
			if ((pen = rhspen) <= 0.)
				pen = Infinity;
			lu = LUrhs;
			lue = lu + 2*nc;
			if (drhs && (p = drhs->u.r)) do {
				t = *p++;
				if (lu[0] > negInfinity
				 && lu[1] < Infinity
				 && lu[0] < lu[1]) {
					if (t <= 0.)
						t = pen;
					*lbr++ = *ubr++ = t;
					}
				} while((lu += 2) < lue);
			else do
				*lbr++ = *ubr++ = pen;
				while(--nr > 0);
			}
		}
	if ((rt = feasrelax - 1) >= 3) {
		rt -= 3;
		mr = 1;
		rv = fto;
		}
	else {
		mr = 0;
		rv = 0;
		*objqual = "feastol ";
		}
	if (GRBfeasrelax(mdl, rt, mr, lbp, ubp, rhp, fto))
		failed(env, "GRBfeasrelax");
	return rv;
	}
#endif /*}*/

 static int
xround(real *x, int n, int assign, real *emax)
{
	int m;
	real d, e, *xe, y;

	e = *emax;
	m = 0;
	for(xe = x + n; x < xe; x++) {
		y = floor(*x + 0.5);
		if ((d = *x - y) != 0.) {
			if (d < 0)
				d = -d;
			if (e < d)
				e = d;
			m++;
			if (assign)
				*x = y;
			}
		}
	*emax = e;
	return m;
	}

 static int
solround(ASL *asl, real *x, real *emax)
{
	int assign, nint, nround;
	real *x2;

	assign = Round & 1;
	*emax = 0.;
	nround = 0;
	if ((nint = niv + nbv)) {
		x2 = x + n_var - nint;
		nround = xround(x2, nint, assign, emax);
		}
	if ((nint = nlvbi)) {
		x2 = x + (nlvb - nint);
		nround += xround(x2, nint, assign, emax);
		}
	if ((nint = nlvci)) {
		x2 = x + (nlvc - nint);
		nround += xround(x2, nint, assign, emax);
		}
	if ((nint = nlvoi)) {
		x2 = x + (nlvo - nint);
		nround += xround(x2, nint, assign, emax);
		}
	if (*emax <= round_reptol)
		nround = 0;
	return nround;
	}

 static const char*
statmsg(ASL *asl, GRBenv *env, GRBmodel *mdl, int i, Dims *d, int *wantobj)
{
	char buf[64], *rv1;
	const char *rv, *rva;
	int m, n, nc, ns, nv, nvr, objwant, sr, srd, srp;
	real *x, *y;
	size_t L;
#if GRB_VERSION_MAJOR >= 7
	real b[4], obj, objdiff, t;
#endif

	*wantobj = 0;
	if (i) {
		sr = 502;
		d->x = d->y = 0;
		switch(i) {
		  case GRB_ERROR_OUT_OF_MEMORY:
			rv = "ran out of memory";
			break;
		  case GRB_ERROR_NO_LICENSE:
			if (!(rv = GRBgeterrormsg(env)) || !*rv)
				rv = "invalid license";
			break;
		  case GRB_ERROR_SIZE_LIMIT_EXCEEDED:
			rv = "problem size limit exceeded";
			break;
		  case GRB_ERROR_IIS_NOT_INFEASIBLE:
			rv = "bug: IIS problem is infeasible";
			break;
#ifdef GRB_ERROR_NETWORK
		  case GRB_ERROR_NETWORK:
			rv = "Gurobi Compute Server not reached";
			sr = 540;
			break;
#endif
#ifdef GRB_ERROR_JOB_REJECTED
		  case GRB_ERROR_JOB_REJECTED:
			rv = "Rejected by Gurobi Compute Server -- perhaps the\n"
				"queue was too full or queueing time was exceeded.";
			sr = 541;
			break;
#endif
#ifdef GRB_ERROR_NOT_SUPPORTED
		  case GRB_ERROR_NOT_SUPPORTED:
			if (server) {
				rv = "Feature not supported by Gurobi Compute Server";
				sr = 542;
				}
			else {
				rv = "Feature not supported";
				sr = 543;
				}
			break;
#endif
#if GRB_VERSION_MAJOR >= 4 /*{*/
		  case GRB_ERROR_Q_NOT_PSD:
			rv = "quadratic objective is not positive definite";
#if GRB_VERSION_MAJOR >= 5
			if (nlc)
				rv = objno > 0 && objno <= nlo
					? "quadratic objective or constraint is not positive definite"
					: "quadratic constraint is not positive definite";
#endif
			sr = 524;
			break;
#endif /*}*/
		  default:
			Snprintf(buf, sizeof(buf), "surprise return %d from GRBoptimize", i);
			rv = rv1 = M1alloc(strlen(buf)+1);
			strcpy(rv1, buf);
		  }
		return rv;
		}
	if (GRBgetintattr(mdl, GRB_INT_ATTR_STATUS, &i))
		failed(env, "GRBgetintattr(STATUS)");
	nc = n_con;
	nv = n_var;
	nvr = nv + nranges;
	objwant = 1;
	sr = 0;
	switch(i) {
	  case GRB_OPTIMAL:
		rv = "optimal solution";
		break;
	  case GRB_INFEASIBLE:
		objwant = 0;
		nc = srd = 0;
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
		if (rays & 2)
			srd = send_dray(asl, d, env, mdl);
#endif /*}*/
		if (!want_iis) {
			sr = 200;
			rv = "infeasible";
			}
		else if (GRBcomputeIIS(mdl)) {
			sr = 202;
			rv = "infeasible; no IIS found";
			}
		else
			rv = retiis(asl, env, mdl, d, "infeasible", &sr);
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
		if (srd) {
 have_srd:
			sr += 3;
			rv1 = (char*)M1alloc(L = strlen(rv) + 40);
			snprintf(rv1, L, "%s; constraint.dunbdd returned.", rv);
			rv = rv1;
			}
#endif /*}*/
		break;
	  case GRB_INF_OR_UNBD:
		objwant = 0;
		nc = nv = 0;
#if GRB_VERSION_MAJOR > 4 || (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) /*{*/
		srd = srp = 0;
		if (rays & 1)
			srp = send_ray(asl, d, env, mdl);
		if (rays & 2 && srp != 1)
			srd = send_dray(asl, d, env, mdl);
#endif /*}*/
		if (!want_iis) {
			sr = 300;
			rv = "infeasible or unbounded";
			}
		else if (GRBcomputeIIS(mdl)) {
			sr = 301;
			rv = "infeasible or unbounded; no IIS";
			srp = srd = 0;
			}
		else
			rv = retiis(asl, env, mdl, d, "infeasible or unbounded", &sr);
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
		if (srd)
			goto have_srd;
		if (srp)
			goto have_srp;
#endif /*}*/
		break;
	  case GRB_UNBOUNDED:
		rv = "unbounded";
		objwant = 0;
		nc = nv = srp = 0;
		sr = 300;
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
		if (rays & 1) {
			srp = send_ray(asl, d, env, mdl);
 have_srp:
			sr += 2;
			rv1 = (char*)M1alloc(L = strlen(rv) + 30);
			snprintf(rv1, L, "%s; variable.unbdd returned.", rv);
			rv = rv1;
			}
#endif /*}*/
		break;
	  case GRB_CUTOFF:
		rv = "objective cutoff";
		sr = 400;
		break;
	  case GRB_ITERATION_LIMIT:
		rv = "iteration limit with a feasible solution";
		rva = "interation limit without a feasible solution";
		sr = 401;
		goto solcheck;
	  case GRB_NODE_LIMIT:
		rv = "node limit with a feasible solution";
		rva = "node limit without a feasible solution";
		sr = 402;
		goto solcheck;
	  case GRB_TIME_LIMIT:
		rv = "time limit with a feasible solution";
		rva = "time limit without a feasible solution";
		sr = 403;
		goto solcheck;
	  case GRB_SOLUTION_LIMIT:
		rv = "solution limit";
		sr = 404;
		break;
	  case GRB_INTERRUPTED:
		rv = "interrupted with a feasible solution";
		rva = "interrupted without a feasible solution";
		sr = 405;
 solcheck:
		ns = 0;
		if (GRBgetintattr(mdl, GRB_INT_ATTR_SOLCOUNT, &ns) || ns <= 0) {
			rv = rva;
			sr += 10;
			}
		break;
	  case GRB_NUMERIC:
		rv = "numeric error";
		sr = 520;
		break;
#ifdef GRB_SUBOPTIMAL
	  case GRB_SUBOPTIMAL:
		rv = "suboptimal";
		sr = 100;
		break;
#endif
#if GRB_VERSION_MAJOR >= 7
	  case GRB_USER_OBJ_LIMIT:
		rv = "bestobjstop or bestbndstop reached";
		sr = 103;
		if (!GRBgetdblattr(mdl, GRB_DBL_ATTR_OBJVAL, &obj)) {
			objdiff = Infinity;
			if (!GRBgetdblparaminfo(env, GRB_DBL_PAR_BESTBDSTOP, b,b+1,b+2,b+3)
			 && (t = d->objsense*(b[0] - obj)) >= 0.) {
				objdiff = t;
				rv = "bestbndstop reached";
				sr = 101;
				}
			if (!GRBgetdblparaminfo(env, GRB_DBL_PAR_BESTOBJSTOP, b,b+1,b+2,b+3)
			 && (t = d->objsense*(b[0] - obj)) >= 0. && t <= objdiff) {
				rv = "bestobjstop reached";
				sr = 102;
				}
			}
		break;
#endif
	  default:
		Snprintf(buf, sizeof(buf), "surprise status %d after GRBoptimize", i);
		rv = rv1 = (char*)M1alloc(strlen(buf)+1);
		sr = 530;
		strcpy(rv1, buf);
	  }
	x = y = 0;
	if ((n = nvr + nc + 2*n_lcon) > 0) {
		x = (real*)M1alloc(n*sizeof(real));
		d->y0 = y = x + nvr;
		}
	m = 0;
	if (!nv)
		x = 0;
	else if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_X, 0, nvr, x)) {
		x = 0;
		m = 1;
		if (sr == 103) {
			sr = 104;
			rv = "bestobjstop or bestbndstop reached with no solution available";
			}
		else if (sr < 200) {
			sr = 570;
			rv = "solution found but not available (Gurobi bug?)";
			}
		}
	solve_result_num = sr;
	if (!nc)
		y = 0;
#if GRB_VERSION_MAJOR >= 5
	else if ((nlc && GRBgetdblattrarray(mdl, GRB_DBL_ATTR_QCPI, 0, nlc, y))
		|| (nc > nlc && GRBgetdblattrarray(mdl, GRB_DBL_ATTR_PI, 0, nc-nlc, y+nlc)))
#else
	else if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_PI, 0, nc, y))
#endif
		{
		y = 0;
		m += 2;
		}
	d->missing = m;
	d->x = x;
	d->y = y;
	if (!x)
		objwant = 0;
	*wantobj = objwant;
	return rv;
	}

typedef char *(*what_name)(ASL*,int);

 static int
stat_map(ASL *asl, int *stat, int n, int *map, int mx, const char *what, what_name wn)
{
	int bad, i, i1, j, j1;
	static char badfmt[] = "gurobi driver: %s.sstatus = %d\n";

	bad = i1 = j1 = 0;
	for(i = 0; i < n; i++) {
		if ((j = stat[i]) >= 0 && j <= mx)
			stat[i] = map[j];
		else {
			stat[i] = 0;
			i1 = i;
			j1 = j;
			if (!bad++)
				fprintf(Stderr, badfmt, wn(asl,i), j);
			}
		}
	if (bad > 1) {
		if (bad == 2)
			fprintf(Stderr, badfmt, what, i1, j1);
		else
			fprintf(Stderr,
		"gurobi driver: %d messages about bad %s.sstatus values suppressed.\n",
				bad-1, what);
		}
	return bad;
	}

#ifdef DO_BASIS_CHECK
#define Basis_check(x) x
#else
#define Basis_check(x)
#endif

 static int
get_input_statuses(ASL *asl, GRBenv *env, GRBmodel *mdl, Dims *d, real *resid)
{
	int bad, *cs, i, m, n, nvr, *rs, *rsta;
#ifdef DO_BASIS_CHECK
	int nb;
#endif
	real *lu;
#ifdef OLD_CMAP_IN
	static int cmap[] = {-1, 0, -1, -1, -1, -1, -1};
#else
	static int cmap[] = {0, 0, -1, -1, -1, -1, -1};
#endif
	static int vmap[] = {-3, 0, -3, -1, -2, -3, -3};

	if (!(d->csd->kind & ASL_Sufkind_input) && !(d->rsd->kind & ASL_Sufkind_input))
		return 0;
	n = n_var;
	m = n_con;
	cs = d->cstat;
	rs = d->rstat;
	nvr = n + nranges;
	Basis_check(nb = 0;)	/* number of basic entities */
	if (nvr > n) {
		rsta = cs + n;
		lu = LUrhs;
		for(i = 0; i < m; ++i, lu += 2) {
			if (lu[0] > negInfinity && lu[0] < lu[1] && lu[1] < Infinity) {
				/* a range constraint */
				if ((resid ? (resid[i] >= lu[1] || resid[i] <= lu[0])
					   : (0. >= lu[1] || 0. <= lu[0]))) {
					rs[i] = 1;
					*rsta = 0;
					Basis_check(++nb;) /* basic range constraint */
					}
				else {
					rs[i] = 0;
					*rsta = 1;
					}
				++rsta;
				}
#ifdef DO_BASIS_CHECK
			else if (rs[i] == 1)
				++nb;
#endif
			}
		}
#ifdef DO_BASIS_CHECK
	else {
		for(i = 0; i < m; ++i) {
			if (rs[i] == 1)
				++nb;
		}	}
	for(i = 0; i < n; ++i) {
		if (cs[i] == 1)
			++nb;
		}
	if (nb > m) {
		/* too big a basis */
		for(i = 0; i < m; ++i) { /* consider constraints first */
			if (rs[i] == 1) {
				rs[i] = 0;
				if (--nb <= m)
					goto have_basis;
				}
			}
		for(i = 1; i < n; ++i) { /* now try variables */
			if (cs[i] == 1) {
				cs[i] = 0;
				if (--nb <= m)
					goto have_basis;
				}
			}
		}
	else if (nb && nb < m) { /* too small a basis; nb should be > 0 */
		for(i = 0; i < m; ++i) { /* first try to make new constraints basic */
			if (rs[i] == 0) {
				rs[i] = 1;
				if (++nb >= m)
					goto have_basis;
				}
			}
		for(i = 0; i < m; ++i) {
			if (rs[i] != 1) {
				rs[i] = 1;
				if (++nb >= m)
					goto have_basis;
				}
			}
		for(i = 0; i < n; ++i) {
			if (cs[i] != 1) {
				cs[i] = 1;
				if (++nb >= m)
					goto have_basis;
				}
			}
		}
 have_basis:
#endif
	bad = stat_map(asl, cs, n, vmap, 6, "_svar", var_name_ASL);
	bad += stat_map(asl, rs, m, cmap, 6, "_scon", con_name_ASL);
	if (bad)
		return 0;
	if (GRBsetintattrarray(mdl, GRB_INT_ATTR_VBASIS, 0, nvr, cs))
		failed(env, "GRBsetintattrarray(\"VBasis\")");
	if (GRBsetintattrarray(mdl, GRB_INT_ATTR_CBASIS, 0, m, rs))
		failed(env, "GRBsetintattrarray(\"CBasis\")");
	return 1;
	}

 static void
intbasis_fail(Dims *d, const char *call)
{
	dpf(d, "\nintbasis trouble: GRB%s failed.", call);
	}

 static GRBmodel*
fixed_model(ASL *asl, GRBmodel *mdl0, Dims *d)
{
	Filename *fn;
	GRBenv *env;
	GRBmodel *mdl;
	double f, *y;
	int i, k;
#if GRB_VERSION_MAJOR >= 5
	int m, m1, nqc;
#endif
	static char *statusname[] = {
		"infeasible",
		"infeasible or unbounded",
		"unbounded",
		"cutoff",
		"iteration limit",
		"node limit",
		"time limit",
		"solution limit",
		"interrupted",
		"numeric difficulty",
		"suboptimal"
		};

	if (!(mdl = GRBfixedmodel(mdl0)))
		return 0;
	if (!(env = GRBgetenv(mdl))) {
		dpf(d, "\nGRBgetenv failed in fixed_model().");
 badret:
		GRBfreemodel(mdl);
		return 0;
		}
	if (GRBsetintparam(env, "Presolve", 0)) {
		intbasis_fail(d, "setintparam(\"Presolve\")");
		goto badret;
		}
	k = -12345;
	GRBgetintparam(env, Method, &k);
	if (fixedmethod < -1 || fixedmethod >= 5) { /* not specified or invalid */
		if (k >= 2 || k < 0)
			fixedmethod = 1;
		else
			fixedmethod = k;
		}
	if (k != fixedmethod)
		GRBsetintparam(env, Method, fixedmethod);
	if ((fn = Wflist[2])) {
		GRBupdatemodel(mdl);
		do {
			if (GRBwrite(mdl, fn->name))
				enamefailed(env, "GRBwrite", fn->name);
			} while((fn = fn->next));
		}
	if (GRBoptimize(mdl)) {
		intbasis_fail(d, "optimize()");
		goto badret;
		}
	if (GRBgetintattr(mdl, GRB_INT_ATTR_STATUS, &i)) {
		intbasis_fail(d, "getintattr()");
		goto badret;
		}
	if (i != GRB_OPTIMAL) {
		if (i >= GRB_INFEASIBLE && i <= GRB_SUBOPTIMAL)
			dpf(d, "\nGRBoptimize of fixed model: %s.",
				statusname[i-GRB_INFEASIBLE]);
		else
			dpf(d, "\nSurprise status %d after GRBoptimize of fixed model.",
				i);
		goto badret;
		}
	if (d->missing & 2 && (y = d->y0)) {
#if GRB_VERSION_MAJOR < 5
		if (!GRBgetdblattrarray(mdl, GRB_DBL_ATTR_PI, 0, n_con, y))
#else
		m = n_con;
		nqc = nlc;
		k = 0;
		if (nqc > 0)
			k = GRBgetdblattrarray(mdl, GRB_DBL_ATTR_QCPI, 0, nqc, y);
		if ((m1 = m - nqc) > 0 && !k)
			k = GRBgetdblattrarray(mdl, GRB_DBL_ATTR_PI, 0, m1, y+nqc);
		if (!k)
#endif
			{
			d->y = y;
			d->missing &= ~2;
			}
		}
	if (!GRBgetdblattr(mdl, GRB_DBL_ATTR_ITERCOUNT, &f)) {
		if (f > 0.)
			dpf(d, "\nplus %.0f simplex iteration%s for intbasis",
				f, "s" + (f == 1.));
		}
	return mdl;
	}
#undef Method

 static void
get_output_statuses(ASL *asl, GRBmodel *mdl, Dims *d)
{
	int i, j, m, n, nr, rv, *s;
#ifndef OLD_CMAP_OUT
	int k, *rsta;
	real *lu;
	static int cmap[6] = {2, 1, 2, 4, 3, 1};
#endif
	static int vmap[4] = {2, 4, 3, 1};

	m = n_con;
	n = n_var;
	nr = n + nranges;
	rv = 1;
	if (GRBgetintattrarray(mdl, GRB_INT_ATTR_VBASIS, 0, nr, d->cstat)
	 || GRBgetintattrarray(mdl, GRB_INT_ATTR_CBASIS, 0, m,  d->rstat)) {
		/*failed(env, "GRBgetintattrarray(\"VBasis\")");*/
		d->csd->kind &= ~ASL_Sufkind_output;
		d->rsd->kind &= ~ASL_Sufkind_output;
		goto ret;
		}
	rv = 0;
	s = d->cstat;
	for(i = 0; i < n; ++i) {
		if ((j = s[i] + 3) >= 0 && j <= 3)
			s[i] = vmap[j];
		else {
			badretfmt(504, "Surprise VBasis[%d] = %d.", i, j-3);
			goto ret;
			}
		}
#ifndef OLD_CMAP_OUT
	rsta = s + n;
	lu = LUrhs;
#endif
	s = d->rstat;
	for(i = 0; i < m; ++i
#ifndef OLD_CMAP_OUT
			      , lu += 2
#endif
					) {
		j = s[i] + 1;
		if (j < 0 || j > 1) {
			badretfmt(505, "Surprise CBasis[%d] = %d.", i, j-1);
			goto ret;
			}
#ifdef OLD_CMAP_OUT
		s[i] = j;
#else
		if (lu[0] > negInfinity && lu[0] < lu[1] && lu[1] < Infinity) {
			if ((k = *rsta + 5) < 2 || k > 5) {
				badretfmt(504, "Surprise VBasis[%d] = %d.", rsta - d->cstat, k-5);
				goto ret;
				}
			++rsta;
			if (j != 1)
				j = k;
			}
		s[i] = cmap[j];
#endif
		}
 ret:
	if (rv)
		dpf(d, "\nNo basis.");
	}

 static void
nl_iv_adj(ASL *asl, int j, int k, char *vtype, real *x)
{
	/* This will be needed once gurobi can handle nonlinear discrete variables, */
	/* e.g., in QPs. */

	int i, i0;
	real *L, *U;

	L = LUv;
	U = Uvx;
	i0 = k - j;
	if (vtype)
		for(i = i0; i < k; ++i)
			vtype[i] = L[i] == 0. && U[i] == 1. ? 'B' : 'I';
	if (x)
		for(i = i0; i < k; ++i) {
			if (x[i] < L[i])
				x[i] = L[i];
			else if (x[i] > U[i])
				x[i] = U[i];
			else
				x[i] = floor(x[i] + .5);
			}
	}

 typedef struct
Sensname {
	char *aname, *gname;
	int iscon;
	} Sensname;

 static void
put_sens(ASL *asl, GRBmodel *mdl)
{
	Sensname *sn, *sne;
	int i, len[3], nc, nr, nv, nvr;
	real *a, *a0, *lu, *rlo2, *rhi2, *ubhi, *ublo;
	size_t L;
	static Sensname Snames[] = {
		{ "sensobjlo", "SAObjLow", 0 },
		{ "sensobjhi", "SAObjUp", 0 },
		{ "sensrhslo", "SARHSLow", 1 },
		{ "sensrhshi", "SARHSUp", 1 },
		{ "senslblo",  "SALBLow", 2 },
		{ "senslbhi",  "SALBUp", 2 },
		{ "sensublo",  "SAUBLow", 2 },
		{ "sensubhi",  "SAUBUp", 2 }};
	static int ak[2] = {ASL_Sufkind_var, ASL_Sufkind_con};

	nr = nranges;
	len[0] = nv = n_var;
	len[1] = nc = n_con;
	len[2] = nvr = nv + nr;
	L = (6*nv + 2*nc)*sizeof(real);
	if (nr)
		L += (2*nc + 4*nr)*sizeof(real);
	a0 = a = (real*)M1alloc(L);
	for(sn = Snames, sne = sn + sizeof(Snames)/sizeof(Sensname); sn < sne; ++sn) {
		if (GRBgetdblattrarray(mdl, sn->gname, 0, len[sn->iscon], a))
  	  		namefailed("GRBgetdblattrarray", sn->gname);
		a += len[sn->iscon];
		}
	if (nr) {
		memset(rlo2 = a, 0, 2*nc*sizeof(real));
		rhi2 = rlo2 + nc;
		lu = LUrhs;
		ublo = a0 + 3*nv + 2*(nc+nvr);
		ubhi = ublo + nvr;
		for(i = 0; i < nc; ++i, lu += 2) {
			if (lu[0] > negInfinity && lu[1] < Infinity && lu[0] < lu[1]) {
				rlo2[i] = ublo[i] + lu[0];
				rhi2[i] = ubhi[i] + lu[0];
				}
			}
		}
	for(a = a0, sn = Snames; sn < sne; ++sn) {
		suf_rput(sn->aname, ak[sn->iscon & 1], a);
		a += len[sn->iscon];
		}
	if (nr) {
		suf_rput("sensrhslo2", ASL_Sufkind_con, rlo2);
		suf_rput("sensrhshi2", ASL_Sufkind_con, rhi2);
		}
	}

#if GRB_VERSION_MAJOR >= 3 /*{*/

 enum {fname_endlen = 32};

 static int
ams_write(ASL *asl, GRBenv *env, GRBmodel *mdl, Dims *dims, int nsols, real bestobj, int mip1)
{
	Option_Info *oi;
	char *fname, *fname_end, msg[200];
	int havetol, i, j, k, nvr;
	real ba, emax, *c, obj, t, ta, *x;
	size_t L;

	if (nsols > ams_limit && ams_limit > 0)
		nsols = ams_limit;
	nvr = n_var + nranges;
	L = strlen(ams_stub);
	dims->xams = x = (real*)Malloc(nvr*sizeof(real) + L + fname_endlen + sizeof(Option_Info));
	dims->oi = oi = (Option_Info*)(x + nvr);
	memset(oi, 0, sizeof(*oi));
	oi->wantsol = 9;
	dims->fname = fname = (char*)(oi+1);
	dims->fname_end = fname_end = fname + L;
	memcpy(fname, ams_stub, L);
	havetol = 0;
	ba = 0.; /* silence erroreous warning */
	if (ams_eps > 0. || ams_epsabs > 0.) {
		havetol = 1;
		if ((ba = bestobj) < 0.)
			ba = -ba;
		}
	c = dims->c;
	for(i = 1; i <= nsols; ++i) {
		if (mip1) {
			if (GRBsetintparam(env, "SolutionNumber", i-1))
				namefailed("GRBsetintparam", "SolutionNumber");
			if (GRBgetdblattrarray(mdl, GRB_DBL_ATTR_Xn, 0, nvr, x))
				namefailed("GRBgetdblattrarray", GRB_DBL_ATTR_Xn);
			if (!GRBgetdblattr(mdl, GRB_DBL_ATTR_OBJVAL, &obj)) {
				obj = 0.;
				for(j = 0; j < nvr; ++j)
					obj += c[j]*x[j];
				}
			}
		else {
			if (nsols > 1) {
				badretfmt(571, "Expected just one solution when problem "
					"has no integer variables, but found %d.\n", nsols);
				longjmp(Jb, 1);
				}
			memcpy(x, dims->x, nvr*sizeof(real));
			obj = bestobj;
			}
		if (havetol) {
			t = dims->objsense*(obj - bestobj);
			if (ams_epsabs > 0. && t > ams_epsabs)
				break;
			if (ams_eps > 0. && t > 0.) {
				if ((ta = obj) < 0.)
					ta = -ta;
				if (ta < ba)
					ta = ba;
				if (t/ta > ams_eps)
					break;
				}
			}
		j = Snprintf(msg, sizeof(msg),
			"Alternative MIP solution %d, objective = %.*g",
			i, dims->objprec, obj);
		if ((Round & 5) && (k = solround(asl, x, &emax)) && (Round & 4))
			Snprintf(msg+j, sizeof(msg)-j,
				"\n%d integer variables %srounded to integers;"
				" maxerr = %g", k, Round & 1 ? "" : "would be ", emax);
		Snprintf(fname_end, fname_endlen, "%d.sol", i);
		if (write_solf_ASL(asl, msg, x, 0, oi, fname))
			break;
		dpf(dims, "\n%s", msg);
		}
	if (GRBsetintparam(env, "SolutionNumber", 0)) /* restore */
		namefailed("GRBsetintparam", "SolutionNumber");
	return i-1;
	}
#if GRB_VERSION_MAJOR >= 7 /*{*/

 static int
add_indic(void *v, int iv, int compl, int sense, int nz, int *ig, real *g, real rhs)
{
	static char gsense[3] = { GRB_LESS_EQUAL, GRB_GREATER_EQUAL, GRB_EQUAL };

	return GRBaddgenconstrIndicator((GRBmodel*)v, NULL, iv, 1-compl, nz, ig, g, gsense[sense], rhs);
	}

#ifndef NO_MOkwf /*{*/
 typedef struct
MOkwval {
	struct MOkwval *next;
	keyword *kw;
	char *val;
	int Objno;
	} MOkwval;

 static MOkwval *first_MOkw, **pnext_MOkw = &first_MOkw;

 static fint
MOkwf(char *key, fint klen) {
	ASL *asl;
	MOkwval *v;
	char *key0;
	int c, on;
	keyword *kw;

	key0 = key;
	if (klen < 7 || strncmp(key, "obj_", 4)) {
 badkey:
		printf("Unrecognized keyword \"%s\"\n", key0);
		return 1;
		}
	on = 0;
	for(key += 4;;) {
		c = *key++;
		if (c == '_')
			break;
		if (c < '0' || c > '9')
			return 1;
		on = 10*on + (c - '0');
		}
	asl = Oinfo.asl;
	if (on <= 0 || on > n_obj) {
		printf("Rejecting obj_%d; obj_n must have 1 <= n <= %d\n", on, n_obj);
		return 1;
		}
	if (!(kw = (keyword*)b_search_ASL(Oinfo.keywds, (int)sizeof(keyword),
					Oinfo.n_keywds, &key, &Oinfo.eqsign)))
		goto badkey;
	*pnext_MOkw = v = (MOkwval*)mem(sizeof(MOkwval) + strlen(key) + 1);
	pnext_MOkw = &v->next;
	v->next = 0;
	v->kw = kw;
	strcpy(v->val = (char*)(v+1), key);
	v->Objno = on - 1;
	return 0;
	}
#endif /*}*/

#endif /*}*/
#endif /*}*/

enum { MBL = 8192 };

 static GRBenv *env0;

 FILE *GRB_nl_ASL;	/* Foil gcc optimization bug: with -O, the "nl = 0" */
			/* assignment below does not happen, and after longjmp */
			/* we get an erroneous attempt to fclose(nl). */

 static void
setSuffixKappa(ASL* asl, GRBmodel* mdl, double kvalue)
{
	double* buffer = (real*)M1zapalloc(nobjno * sizeof(real));
	buffer[objno - 1] = kvalue;
	suf_rput("kappa", ASL_Sufkind_obj, buffer);
	suf_rput("kappa", ASL_Sufkind_prob, buffer);
}

 static double
reportDblProblemSuffix(ASL* asl, GRBmodel* mdl,
	 const char* grbAttr, const char* suffixName)
{
	double val;
	int error;
	error = GRBgetdblattr(mdl, grbAttr, &val);
	if(error)
		 failed(env0, "GRBgetdblattr");
	suf_rput(suffixName, ASL_Sufkind_prob, &val);
	return val;
 	}

 static int
reportIntProblemSuffix(ASL* asl, GRBmodel* mdl,
	 const char* grbAttr, const char* suffixName)
{
	int val, error;
	error = GRBgetintattr(mdl, grbAttr, &val);
	if (error)
		 failed(env0, "GRBgetintattr");
	suf_iput(suffixName, ASL_Sufkind_prob, &val);
	return val;
	}

 int
main(int argc, char **argv)
{
	ASL *asl;
	CS_type *cs, *csr, j, j1, k, k1, nz, nzcr;
	Dims dims;
	/*FILE *nl;*/ #define nl GRB_nl_ASL
	Filename *fn;
	GRBenv *env;
	GRBmodel *mdl, *mdl1;
	char mbuf[MBL];
	char **cnames, **cnam0, *hx0, *sense, *smsg, *sostype, *stub, **vnames, *vtype;
	const char *sc, *solmsg;
#if GRB_VERSION_MAJOR >= 8
	SufDesc *partition;
	const char *group, *router;
	int *part;
#endif
	int i, i1, lvi, mip, mip1, nc, nfree, nlvi, no, nqc, nsosnz, nrange, nsos;
	int nv, nvr, objno1, rc, wantobj;
	int *cmsave, *rnr, *rsta, *sosbeg, *sosind, *sostypes;
	int vinfo[3], *vlen, *vlenr, *vmi;
	ograd *og;
	real absmipgap, f, obj, oc, relmipgap, t;
	real *A, *Ar, *bb, *lu, *lxr, *rhs, *sosref, *uxr, *x, *y;
	sig_func_type *oic;
	size_t L;
	static int sos_types[2] = { GRB_SOS_TYPE1, GRB_SOS_TYPE2 };
#if GRB_VERSION_MAJOR >= 4 /*{*/
	QPinfo *qpif = 0;
	void *v;
	int *colno, *colq, i2, nelq, nqcol, *rowq;
	real *qmat;
	size_t *colbeg, is, js;
	ssize_t nelqf = 0;
#if GRB_VERSION_MAJOR >= 5 /*{*/
	QPinfo *qpi;
	cgrad *cg, **cgp;
	char qsense;
	const char *objqual = "";
	int *ia, nlnz, *rn, *zc, *zr;
	real *a1, fto, *pfto, qrhs, *resid;
	ssize_t nelqc;
#ifdef GRB_INT_PAR_TUNEOUTPUT
	char *tunemsg;
#endif
#if GRB_VERSION_MAJOR >= 6 /*{*/
	PLterm *plt = 0, *plt0;
	SufDesc *lazyd;
	int *lz, rv;
#endif /*}*/
#endif /*}*/
#if GRB_VERSION_MAJOR >= 7 /*{*/
	cde *ocd;
	expr *en;
	int errinfo[2], nlogc, *objpri, *operm;
	real *objabs, *objrel, *objwt;
#ifndef NO_MOkwf
	GRBenv *env1;
	MOkwval *mov, **pmov;
	keyword *kw;
#endif
#define ALLOW_CLP ASL_allow_CLP
#else
#define nlogc 0
#define ALLOW_CLP 0
#endif /*}*/
	static int repmap[4] = { 0, ASL_obj_replace_ineq,  ASL_obj_replace_eq,
				    ASL_obj_replace_ineq | ASL_obj_replace_eq };
#endif /*}*/

	Times[0] = xectim_();
	Times[4] = time(0);
#ifdef LICENSE_FILE
	if (!(stub = getenv("GRB_LICENSE_FILE")) || !*stub)
		putenv("GRB_LICENSE_FILE=" LICENSE_FILE);
#endif

	oic = 0;
	env0 = 0;
	mdl = 0;
	memset(&dims, 0, sizeof(dims));
	vinfo[0] = vinfo[1] = vinfo[2] = 0;
	GRBversion(&vinfo[0], &vinfo[1], &vinfo[2]);
	Snprintf(verbuf, sizeof(verbuf), "Gurobi %d.%d.%d", vinfo[0], vinfo[1], vinfo[2]);
	Lic_info_add_ASL = "Portions Copyright Gurobi Optimization, LLC, 2020.";

#ifdef main
	if (!(asl = asl1))
#endif
	asl = ASL_alloc(ASL_read_fg);
	/*static initialization: Oinfo.option_echo = ASL_OI_tabexpand | ASL_OI_addnewline;*/
	if (!(stub = getstub(&argv, &Oinfo)))
		usage_ASL(&Oinfo, 1);
	nl = jac0dim(stub, 0);
	nqc = nlc;
#if GRB_VERSION_MAJOR < 5
	if (nqc > 0) {
		asl->i.uinfo = "Gurobi can't handle nonlinear constraints.";
		solve_result_num = 522;
		rc = 1;
		goto bailout;
		}
#endif
	if (n_cc > 0) {
		asl->i.uinfo = "Gurobi can't handle complementarity constraints.";
		solve_result_num = 567;
		rc = 1;
		goto bailout;
		}
	if (!(nobjno = no = n_obj))
		objno = 0;
	rc = 1;
	/* Passing 0 for file logfile; would make logfile an option if we could	*/
	/* specify it after processing $gurobi_options, but we cannot do so,	*/
	/* and we need to have env before reading $gurobi_options, as values	*/
	/* in env may be changed by $gurobi_options. */
	if ((!env0 && (i1 = GRBloadenv(&env0,0))) || !env0) {
		solve_result_num = 500;
		solmsg = GRBgeterrormsg(env0);
		goto ws_now;
		}
	Oinfo.uinfo = (char*)env0;
	nlvi = nlvbi + nlvci + nlvoi;
	lvi = nbv + niv;
	dims.kiv = nlvi + lvi;
	GRBsetintparam(env0, "OutputFlag", 0); /* possibly changed by getopts... */
	rc = setjmp(Jb);
	if (rc) {
 bailout:
		if (nl)
			fclose(nl);
		--rc;
		if (solve_result_num > 0 && asl->i.uinfo){
			if (amplflag | (Oinfo.wantsol & 1))
				rc = 0;
			L = strlen(Oinfo.bsname) + strlen(asl->i.uinfo);
			solmsg = smsg = (char*)M1alloc(L+3);
			sprintf(smsg, "%s: %s", Oinfo.bsname, asl->i.uinfo);
 ws_now:
			write_sol(solmsg, 0, 0, &Oinfo);
			}
		goto done;
		}
	if (getopts(argv, &Oinfo)) {
		if ((solmsg = asl->i.uinfo))
			goto ws_now;
		solve_result_num = 503;
		solmsg = "Bad $gurobi_options.";
		goto ws_now;
		}
#if GRB_VERSION_MAJOR >= 7
	if (ams_limit > 0)
		GRBsetintparam(env0, "PoolSolutions", ams_limit);
#endif
#ifdef ALLOW_GUROBI_SERVER /*{*/
	if (serverlic && server_licread(asl)) {
		solmsg = asl->i.uinfo;
		goto ws_now;
		}
	if (server) {
		if (env0) {
			GRBfreeenv(env0);
			env0 = 0;
			}
#if GRB_VERSION_MAJOR >= 8
		if (!(group = server_group))
			group = "";
		if (!(router = server_router))
			router = "";
#endif
#if GRB_VERSION_MAJOR > 9 || (GRB_VERSION_MAJOR == 9 && GRB_VERSION_MINOR >= 5)
		if (server_timeout < 1)
			server_timeout = 2000000000;
		if ((i = GRBemptyenv(&env0)) || (i = GRBstartenv(env0))
		 || (i = GRBsetstrparam(env0, GRB_STR_PAR_LOGFILE, logfile))
		 || (i = GRBsetstrparam(env0, GRB_STR_PAR_COMPUTESERVER, server))
		 || (i = GRBsetstrparam(env0, GRB_STR_PAR_CSROUTER, router))
		 || (i = GRBsetstrparam(env0, GRB_STR_PAR_SERVERPASSWORD, server_passwd))
		 || (i = GRBsetstrparam(env0, GRB_STR_PAR_CSGROUP, group))
		 || (i = GRBsetintparam(env0, GRB_INT_PAR_CSTLSINSECURE, server_insecure))
		 || (i = GRBsetintparam(env0, GRB_INT_PAR_CSPRIORITY, server_priority))
		 || (i = GRBsetintparam(env0, GRB_INT_PAR_SERVERTIMEOUT, server_timeout)))
#else
		if ((i = GRBloadclientenv(&env0, logfile, server,
#if GRB_VERSION_MAJOR < 8
				server_port, server_passwd,
#else
				router, server_passwd, group, server_insecure,
#endif
				server_priority, server_timeout)))
#endif
		{
			switch(i) {
			 case GRB_ERROR_NETWORK:
				solmsg = "Could not talk to Gurobi Compute Server.";
				solve_result_num = 601;
				break;
			 case GRB_ERROR_JOB_REJECTED:
				solmsg = "Job rejected by Gurobi Compute Server.";
				solve_result_num = 602;
				break;
			 case GRB_ERROR_NO_LICENSE:
				solmsg = "No license for specified Gurobi Compute Server.";
				solve_result_num = 603;
				break;
			 default:
				solmsg = mbuf;
				snprintf(mbuf, sizeof(mbuf),
					"Surprise return %d from GRBloadclientenv().", i);
				solve_result_num = 604;
			 }
			goto ws_now;
			}
#if GRB_VERSION_MAJOR >= 7
 new_env:
#endif
		Oinfo.uinfo = (char*)env0;
		logfile = 0;
		GRBsetintparam(env0, "OutputFlag", 0);
		Replay(env0);
		}
#if GRB_VERSION_MAJOR >= 7 /*{*/
	else if (cloudid && cloudkey) {
		if (env0) {
			GRBfreeenv(env0);
			env0 = 0;
			}
#if GRB_VERSION_MAJOR > 9 || (GRB_VERSION_MAJOR == 9 && GRB_VERSION_MINOR >= 5)
		if ((i = GRBemptyenv(&env0))
		 || (i = GRBsetstrparam(env0, GRB_STR_PAR_CLOUDHOST, cloudhost))
		 || (i = GRBsetstrparam(env0, GRB_STR_PAR_CLOUDACCESSID, cloudid))
		 || (i = GRBsetstrparam(env0, GRB_STR_PAR_CLOUDSECRETKEY, cloudkey))
		 || (i = GRBsetstrparam(env0, GRB_STR_PAR_CLOUDPOOL, cloudpool))
		 || (i = GRBsetintparam(env0, GRB_INT_PAR_CSPRIORITY, server_priority))
		 || (i = GRBstartenv(env0))
		 || (i = GRBsetstrparam(env0, GRB_STR_PAR_LOGFILE, logfile)))
#else
		if ((i = GRBloadcloudenv(&env0, logfile, cloudid, cloudkey, cloudpool
#if GRB_VERSION_MAJOR >= 8
				, cloudpriority
#endif
			)))
#endif
			{
			switch(i) {
			 case GRB_ERROR_NETWORK:
				solmsg = "Could not talk to Gurobi Instant Cloud.";
				solve_result_num = 601;
				break;
			 case GRB_ERROR_JOB_REJECTED:
				solmsg = "Job rejected by Gurobi Instant Cloud.";
				solve_result_num = 602;
				break;
			 case GRB_ERROR_NO_LICENSE:
				solmsg = "No license for specified Gurobi Instant Cloud.";
				solve_result_num = 603;
				break;
			 case GRB_ERROR_CLOUD:
				solmsg = "Bad value for cloudid or cloudkey, or Gurobi Cloud out of reach.";
				solve_result_num=605;
				break;
			 default:
				solmsg = mbuf;
				snprintf(mbuf, sizeof(mbuf),
					"Surprise return %d from GRBloadcloudenv().", i);
				solve_result_num = 604;
			 }
			goto ws_now;
			}
		goto new_env;
		}
#endif /*}*/
#endif /*}*/
	if (relax)
		dims.kiv = 0;
	if (!(mip = dims.kiv))
		Round = 0;
	breaking = 3;
	oic = signal(SIGINT, intcatch);
	nrange = nranges;
	dims.nv0 = nv = n_var;
	nvr = nv + nrange;
	dims.nc0 = nc = n_con;
	nz = nzc;
	nzcr = nz + nrange;
	L = (2*(nvr+nc)+nzcr+nrange) * sizeof(real)
		+ (nvr + 1)*sizeof(CS_type)
		+ (nvr + nzcr)*sizeof(int) + nc + nvr + nv;
	A_vals = A = (real*)Malloc(L);
	Ar = A + nz;
	LUv = lxr = A + nzcr;
	Uvx = uxr = lxr + nvr;
	y = uxr + nvr;
	rhs = y + nc;
	A_colstarts = cs = (CS_type*)(rhs + nc);
	vlen = vlenr = (int*)(cs + nvr + 1);
	A_rownos = rnr = vlen + nvr;
	sense = (char*)(rnr + nzcr);
	vtype = sense + nc;
	havex0 = hx0 = vtype + nvr;
	suf_declare(suftab, sizeof(suftab)/sizeof(SufDecl));
#if GRB_VERSION_MAJOR >= 7
	nlogc = 2*n_lcon;
#endif
	dims.cstat = (int*)M1zapalloc((nvr+nc+nlogc+2)*sizeof(int));
	dims.rstat = dims.cstat + nvr + 1;
	rsta = dims.cstat + nv;
	dims.csd = suf_iput("sstatus", ASL_Sufkind_var, dims.cstat);
	dims.rsd = suf_iput("sstatus", ASL_Sufkind_con, dims.rstat);

#if GRB_VERSION_MAJOR >= 5
	resid = 0;
	if (warmstart == 2)
		basis &= ~1;
	if (mip) {
		want_xpi0 = 1;
		if (basis & 1
		 && dims.csd->kind & ASL_Sufkind_input
		 && dims.rsd->kind & ASL_Sufkind_input
		 && warmstart == 1)
			warmstart = 0;
		}
	else if (warmstart >= 2)
		warmstart = 1;
	if (warmstart) {
		want_xpi0 = 3;
		if (nrange) {
			pi0 = 0;
			resid = y;
			asl->i.nsufext[ASL_Sufkind_var] += nrange;
			/* include space for range slacks in X0 */
			}
		}
#else
	if (mip)
		want_xpi0 = 1;
#endif
#if GRB_VERSION_MAJOR >= 4
	qp_read(nl, ALLOW_CLP | USE_Z | repmap[objrep]);
	nv = n_var;
	nc = n_con;
	nqc = nlc;
	nz = nzc;
	nvr = nv + nrange;
	v = 0;
#else
	fg_read(nl,0);
#endif
	nl = 0;	/* was closed by qp_read */
	nfree = nsos = 0;
	dims.c = 0;
	dims.objsense = 1;
	oc = 0.;
	qpif = 0;
	if ((objno1 = objno - 1) >= 0 && objno1 < n_obj) {
#if GRB_VERSION_MAJOR < 4
		if (objno1 >= no - nlo) {
			asl->i.uinfo = "Gurobi cannot handle nonlinear objectives.";
			solve_result_num = ASL_readerr_nonlin;
			rc = 1;
			goto bailout;
			}
#else
#if GRB_VERSION_MAJOR >= 6 /*{*/
		plt = get_plterms((ASL_fg*)asl, objno1);
#endif /*}*/
		/* must call mqpcheck() after qp_read so objconst() will work right */
		if (nlo && (nelqf = mqpcheckv(objno1, &qpif, &v)) < 0) {
			if (nelqf == -2) {
				solve_result_num = 523;
				asl->i.uinfo = "Cannot handle a quadratic objective involving division by 0";
				}
			else {
				solve_result_num = 521;
				asl->i.uinfo = "Gurobi cannot handle general nonlinear objectives.";
				}
			rc = 1;
			goto bailout;
			}
#endif
		if (objtype[objno1])
			dims.objsense = -1;
		oc = objconst(objno1);
		dims.c = M1zapalloc(nvr * sizeof(real));
		if (asl->i.vmap) {
			vmi = get_vminv_ASL(asl);
			for(og = Ograd[objno1]; og; og = og->next)
				dims.c[vmi[og->varno]] = og->coef;
			}
		else
			for(og = Ograd[objno1]; og; og = og->next)
				dims.c[og->varno] = og->coef;
		}
	else
		objno1 = -1;
	obj_no = objno1;	/* for write_sol() */
	if (mip || (sos
	    && suf_get("sosno", ASL_Sufkind_var | ASL_Sufkind_input)
	    && suf_get("ref", ASL_Sufkind_var | ASL_Sufkind_input))) {
		i = ASL_suf_sos_explict_free;
		if (!sos)
			i |= ASL_suf_sos_ignore_sosno;
		if (!sos2)
			i |= ASL_suf_sos_ignore_amplsos;
		if ((nsos = suf_sos(i, &nsosnz, &sostype, 0, 0,
				&sosbeg, &sosind, &sosref))) {
			nv = n_var;
			nrange = nranges;
			nvr = nv + nrange;
			nc = n_con;
			nz = nzc;
			nzcr = nz + nrange;
			lvi = nbv + niv;
			}
		}
	Ar = A + nz;
	rnr = A_rownos + nz;
	lxr = LUv + nv;
	uxr = Uvx + nv;
	rsta = dims.cstat + nv;
#if GRB_VERSION_MAJOR >= 5 /*{*/
	if (!pi0) {
		if (nc) {
			if (warmstart <= 2)
				warmstart = 0;
			memset(y, 0, nc*sizeof(real));
			}
		}
	else if (warmstart && nrange && (x = X0)) {
		memset(resid, 0, nc*sizeof(real));
		rn = A_rownos;
		for(j = i = 0; i < nv; ) {
			t = x[i];
			for(k = cs[++i]; j < k; ++j)
				resid[rn[j]] += t*A[j];
			}
		}
#endif /*}*/
	for(i = 0; i < nv; ++i)
		*vlenr++ = cs[i+1] - cs[i];
	csr = cs + nv;
	for(lu = LUrhs, i = 0; i < nc; ++i, lu += 2) {
		rhs[i] = lu[0];
		if (lu[0] <= negInfinity) {
			if (lu[1] >= Infinity)
				++nfree;
			else {
				sense[i] = '<';
				rhs[i] = lu[1];
				}
			}
		else if (lu[1] >= Infinity)
			sense[i] = '>';
		else {
			sense[i] = '=';
			if (lu[1] > lu[0]) {
				*rnr++ = i;
				*csr++ = Ar - A;
				*Ar++ = -1.;
				*vlenr++ = 1;
				*lxr++ = 0.;
				*uxr++ = lu[1] - lu[0];
				*rsta++ = y[i] > 0. ? -1 : y[i] < 0. ? -2 : 0;
				}
			}
		}
	if (nfree) {
		fprintf(stderr, "The problem has %d free row%s, which gurobi"
			" cannot handle.\nYou need to let AMPL's presolve "
			"remove free rows.\n", nfree, nfree > 1 ? "s" : "");
		return 1;
		}
	if (nvr)
		*csr = Ar - A;

	memset(vtype, 'C', nvr);
	if (mip) {
		if (nlvi) {
			k = nlvb;
			if ((j = nlvbi))
				nl_iv_adj(asl, j, k, vtype, X0);
			k = nlvc;
			if ((j = nlvci))
				nl_iv_adj(asl, j, k, vtype, X0);
			k += nlvo - nlvc;
			if ((j = nlvoi))
				nl_iv_adj(asl, j, k, vtype, X0);
			}
		k = nv - lvi;
		if ((j = nbv)) {
			memset(vtype+k, 'B', j);
			k += j;
			if (X0)
				nl_iv_adj(asl, j, k, 0, X0);
			}
		if ((j = niv)) {
			memset(vtype+k, 'I', j);
			k += j;
			if (X0)
				nl_iv_adj(asl, j, k, 0, X0);
			}
		}
	cnam0 = cnames = vnames = 0;
	if (maxcolnamelen + maxrownamelen > 0
	 && (Wflist[0] || Wflist[1] || Wflist[2])) {
		if (maxcolnamelen) {
			vnames = (char**)mem(nvr*sizeof(char*));
			for(i = 0; i < nv; ++i)
				vnames[i] = var_name(i);
			rnr -= nrange;
			for(i = 0; i < nrange; ++i) {
				i1 = snprintf(mbuf, sizeof(mbuf), "%s$range",
						con_name(rnr[i])) + 1;
				memcpy(vnames[nv+i] = (char*)mem(i1), mbuf, i1);
				}
			}
		if (maxrownamelen) {
			cnam0 = cnames = (char**)mem(nc*sizeof(char*));
			for(i = 0; i < nc; ++i)
				cnames[i] = con_name(i);
			cnames += nqc;
			}
		}
#if GRB_VERSION_MAJOR >= 5 /*{*/
	if (nqc) {
		if (mip && (GRBgetintparam(env0, "QCPDual", &i) || i == 0))
			basis = solnsens = 0;
		ia = A_rownos;
		j = nzc;
		for(i = nlnz = 0; i < j; i++)
			if (ia[i] < nqc)
				++nlnz;
		for(i = k = 0; i < nqc; ++i) {
			nelqc = mqpcheckv(-(i+1), 0, &v);
			if (k < nelqc)
				k = nelqc;
			}
		a1 = (real*)Malloc(nqc*sizeof(cgrad*) + nlnz*sizeof(cgrad)
				+ nv*(sizeof(int) + sizeof(real)) + k*(2*sizeof(int)));
		cg = (cgrad*)(a1 + nv);
		Cgrad = cgp = (cgrad**)(cg + nlnz);
		zc = (int*)(cgp + nqc);
		zr = zc + k;
		rn = zr + k;
		memset(cgp, 0, nqc*sizeof(cgrad*));
		k = cs[i = nvr];
		while(--i >= 0) {
			j = cs[i];
			j1 = k1 = k;
			while(k > j) {
				if ((i1 = ia[--k]) < nqc) {
					cg->varno = i;
					cg->coef = A[k];
					cg->next = cgp[i1];
					cgp[i1] = cg++;
					}
				else  {
					A[--j1] = A[k];
					ia[j1] = i1 - nqc;
					}
				}
			cs[i] = j1;
			vlen[i] = k1 - j1;
			}
		if (GRBloadmodel(env0, &mdl, "foo",
				nvr, nc-nqc, dims.objsense, oc, dims.c,
				sense+nqc, rhs+nqc, cs, vlen, ia, A,
				LUv, Uvx, vtype,  vnames, cnames) || !mdl)
			failed(env0, "GRBloadmodel");
		lu = LUrhs;
		/* We need cmap = 0 in mqpcheckv since we just recomputed Cgrad. */
		cmsave = asl->i.cmap;
		asl->i.cmap = 0;
		for(i = 0; i < nqc; ++i, lu += 2) {
			nelqc = mqpcheckv(-(i+1), &qpi, &v);
			if (nelqc < 0) {
				switch(nelqc) {
				  case -2:
					solve_result_num = 525;
					asl->i.uinfo =
					 "a quadratic constraint involving division by 0.";
					break;
				  default:
					solve_result_num = 522;
					asl->i.uinfo =
					 "Gurobi can't handle nonquadratic nonlinear constraints.";
					}
				rc = 1;
				goto bailout;
				}
			if (lu[0] <= negInfinity) {
				qsense = '<';
				qrhs = lu[1];
				}
			else if (lu[1] >= Infinity) {
				qsense = '>';
				qrhs = lu[0];
				}
			else {
				if (lu[0] == lu[1]) {
					qsense = '=';
					qrhs = lu[0];
					if (!qpi)
						goto allow_eq;
#if GRB_VERSION_MAJOR >= 9
					if (!GRBgetintparam(env0, GRB_INT_PAR_NONCONVEX, &i1)
					    && i1 >= 2)
						goto allow_eq;
#endif
					solve_result_num = 526;
					asl->i.uinfo =
					 "Gurobi cannot handle quadratic equality constraints"
#if GRB_VERSION_MAJOR >= 9
						"\nunless \"nonconvex=2\" is specified"
#endif
						".";
					}
				else {
					solve_result_num = 527;
					asl->i.uinfo =
					 "Gurobi cannot handle quadratic range constraints.";
					}
				rc = 1;
				goto bailout;
				}
 allow_eq:
			for(j = 0, cg = cgp[i]; cg; cg = cg->next)
				if (cg->coef != 0.) {
					rn[j] = cg->varno;
					a1[j++] = cg->coef;
					}
			if (!qpi) {
				if (GRBaddconstr(mdl, j, rn, a1, qsense, qrhs, 0))
					failed(env0, "GRBaddconstr");
				continue;
				}
			colbeg = qpi->colbeg;
			colno = qpi->colno;
			rowq = qpi->rowno;
			qmat = qpi->delsq;
			nqcol = qpi->nc;
			is = colbeg[0];
			for(i1 = k= 0; i1 < nqcol; ++i1) {
				j1 = colno[i1];
				for(js = colbeg[i1+1]; is < js; ++is) {
					k1 = rowq[is];
					if (k1 <= j1) {
						zc[k] = j1;
						zr[k] = k1;
						if (j1 == k1)
							qmat[is] *= 0.5;
						qmat[k++] = qmat[is];
						}
					}
				}
			if (GRBaddqconstr(mdl, j, rn, a1, k, zr, zc, qmat, qsense,
						qrhs, cnam0 ? cnam0[i] : 0))
				failed(env0, "GRBaddqconstr");
			free(qpi);
			}
		asl->i.cmap = cmsave;
		free(a1);
		GRBupdatemodel(mdl);
		}
	else
#endif /*}*/
	if (GRBloadmodel(env0, &mdl, "foo",
			nvr, nc, dims.objsense, oc, dims.c,
			sense, rhs, cs, vlen, A_rownos, A_vals,
			LUv, Uvx, vtype,  vnames, cnames) || !mdl)
		failed(env0, "GRBloadmodel");
	x = 0;
	if (X0 && (!mip || mipstval)) {
		x = X0;
		for(i = 0; i < nv; ++i)
			if (!hx0[i]) {
				x[i] = GRB_UNDEFINED;
				warmstart = 0;
				}
		}

	if (!(env = GRBgetenv(mdl)))
		failed(env0, "GRBgetenv");
	Oinfo.uinfo = (char*)env;

#if GRB_VERSION_MAJOR >= 4 /*{*/
	if (nelqf) {
		colbeg = qpif->colbeg;
		colno = qpif->colno;
		rowq = qpif->rowno;
		qmat = qpif->delsq;
		nqcol = qpif->nc;
		is = colbeg[0];
		for(i1 = nelq = 0; i1 < nqcol; ++i1) {
			i = colno[i1];
			for(js = colbeg[i1+1]; is < js; ++is)
				if (rowq[is] <= i)
					++nelq;
			}
		colq = (int*)Malloc(nelq*sizeof(int));
		is = colbeg[0];
		for(j1 = i1 = 0; i1 < nqcol; ++i1) {
			i = colno[i1];
			for(js = colbeg[i1+1]; is < js; ++is) {
				if ((i2 = rowq[is]) <= i) {
					colq[j1] = i;
					rowq[j1] = i2;
					qmat[j1] = qmat[is];
					if (i2 == i)
						qmat[j1] *= 0.5;
					++j1;
					}
				}
			}
		if (GRBaddqpterms(mdl, nelq, rowq, colq, qmat))
			failed(env, "GRBaddqpterms");
		free(qpif);
		}
#if GRB_VERSION_MAJOR >= 6 /*{*/
	if ((plt0 = plt)) {
		for(; plt->n; ++plt)
			GRBsetpwlobj(mdl, plt->vno, plt->n, plt->x, plt->y);
		free(plt0);
		}
	if (lvi && lazy && (lazyd = suf_get("lazy", ASL_Sufkind_con | ASL_Sufkind_input))) {
		i1 = nc - nqc;
		lz = lazyd->u.i + nqc;
		for(i = 0; i < i1; ++i) {

#if((GRB_VERSION_MAJOR == 9 && GRB_VERSION_MINOR == 1) || GRB_VERSION_MAJOR > 9) /*{{*/
			if ((i2 = lz[i]) >= -1) {
#if (GRB_VERSION_MAJOR == 9 && GRB_VERSION_MINOR == 1) /*{{{*/
				if (i2 == -1)
				{
					if (sense[i2] == '>')
					{
						fprintf(Stderr, "Cannot specify user cuts that are >= inequalities, see:\n"
							"https://support.gurobi.com/hc/en-us/articles/360052709391\n");
						exit(1);
					}
					if (!Wflist[0]) // If the user did not request a model write, fail, as gurobi would
						// crash due to a bug (reported on 19/11/2020)
					{
						fprintf(Stderr, "Gurobi version 9.1 has a bug that prevents model solution if writeprob is not\n"
							"specified. If the solver crashes, see the 'writeprob' option.\n");
						fflush(Stderr); // flush otherwise the message is not shown
					}
				}
#endif /*}}}*/
				if (i2 == 0)
					continue;
#else /*{{*/
			if ((i2 = lz[i]) > 0) {
#endif /*}}*/
				if (i2 > 3)
					i2 = 3;
				rv = GRBsetintattrelement(mdl, GRB_INT_ATTR_LAZY, i, i2);
#if 1  /*{{*/
				if (rv) {
					fprintf(Stderr, "Surprise return %d from "
						"GRBsetintattrelement(mdl, \"Lazy\", %d, %d)\n",
						rv, i, i2);
					exit(1);
					}
#endif  /*}}*/
				}
			}
		}
#endif /*}*/
	if (v)
		mqpcheckv_free(&v);
#endif /*}*/
	if (nsos) {
		sostypes = (int*)Malloc(nsos*sizeof(int));
		for(i = 0; i < nsos; ++i)
			sostypes[i] = sos_types[sostype[i] - '1'];
		if (GRBaddsos(mdl, nsos, nsosnz, sostypes, sosbeg, sosind, sosref))
			failed(env, "GRBaddsos");
		free(sostypes);
  	  	}
#if GRB_VERSION_MAJOR >= 7 /*{*/
	if (nlogc) {
		if ((i = indicator_constrs_ASL(asl, mdl, add_indic, errinfo))) {
			solmsg = "indicator_constrs_ASL";
			switch(i) {
			  case 1:
				badretfmt(563,
					"logical constraint %s is not an indicator constraint.\n",
					lcon_name(errinfo[0]));
				break;
			  case 2:
				badretfmt(563,
					"logical constraint %s is not an indicator constraint\n\
		due to bad comparison with %s.\n", lcon_name(errinfo[0]), var_name(errinfo[1]));
				break;
			  case 3:
				solmsg = "GRBaddgenconstrIndicator";
				/* no break; */
			  default:
				badretfmt(565, "%s failed; error code %d.", solmsg, i);
			  }
			rc = 1;
			goto bailout;
			}
		if (cnames) {
			nlogc >>= 1;
			GRBupdatemodel(mdl);
			for(i = 0; i < nlogc; ++i)
				GRBsetstrattrelement(mdl, GRB_STR_ATTR_GENCONSTRNAME,
					i, lcon_name_ASL(asl, i));
			}
		}
	operm = 0;
	if (no > 1 && multiobj && objno1 >= 0) {
		ocd = ((ASL_fg*)asl)->I.obj_de_;
		operm = (int*)M1alloc(no*sizeof(int));
		operm[0] = objno1;
		for(i = 0; i < objno1; ++i)
			operm[i+1] = i;
		while(++i < no)
			operm[i] = i;
		objpri = 0;
		objabs = objrel = objwt = 0;
		if ((lazyd = suf_get("objweight", ASL_Sufkind_obj | ASL_Sufkind_input))) {
			objwt = lazyd->u.r;
			for(i = i1 = 0; i < no; ++i) {
				if ((t = objwt[operm[i]]) < 0.) {
					badretfmt(564, "%s.objweight = %.g is negative.\n",
						obj_name(i), objwt[i]);
					rc = 1;
					goto bailout;
					}
				if (t > 0.) {
					en = ocd[operm[i1++] = operm[i]].e;
					if (en && en->op != (efunc*)OPNUM) {
						badretfmt(564, "%s is nonlinear.\n",
							obj_name(operm[i]));
						rc = 1;
						goto bailout;
						}
					}
				}
			if (i1 < 2) {
				badretfmt(564, "Expected at least 2 positive .objweight values; found %d.\n", i1);
				rc = 1;
				goto bailout;
				}
			if (objwt[objno1] <= 0.) {
				badretfmt(564, "Expected a positive weight for objective objno = %d.\n", objno);
				rc = 1;
				goto bailout;
				}
			no = i1;
			}
		if ((lazyd = suf_get("objpriority", ASL_Sufkind_obj | ASL_Sufkind_input))) {
			objpri = lazyd->u.i;
			for(i = 0; i < no; ++i) {
				if (objpri[i1 = operm[i]] < 0) {
					badretfmt(564, "%s.objpriority = %d should be nonnegative.\n",
						obj_name(i1), objpri[i1]);
					rc = 1;
					goto bailout;
					}
				}
			}
		if ((lazyd = suf_get("objabstol", ASL_Sufkind_obj | ASL_Sufkind_input))) {
			objabs = lazyd->u.r;
			for(i = 0; i < no; ++i) {
				if (objabs[i1 = operm[i]] < 0.) {
					badretfmt(564, "%s.objabstol = %.g should be nonnegative.\n",
						obj_name(i1), objabs[i1]);
					rc = 1;
					goto bailout;
					}
				}
			}
		if ((lazyd = suf_get("objreltol", ASL_Sufkind_obj | ASL_Sufkind_input))) {
			objrel = lazyd->u.r;
			for(i = 0; i < no; ++i) {
				if (objrel[i1 = operm[i]] < 0.) {
					badretfmt(564, "%s.objreltol = %.g should be nonnegative.\n",
						obj_name(i1), objrel[i1]);
					rc = 1;
					goto bailout;
					}
				}
			}
		if (GRBsetintattr(mdl, GRB_INT_ATTR_NUMOBJ, no))
			failed(env, "GRBsetintattr(\"NumObj\")");
		for(i = 0; i < no; ++i) {
			if (GRBsetintparam(env, GRB_INT_PAR_OBJNUMBER, i))
				failed(env, "GRBsetintparam(\"ObjNumber\")");
			i1 = operm[i];
			for(og = Ograd[i1]; og; og = og->next) {
				if (GRBsetdblattrelement(mdl, "ObjN", og->varno, og->coef))
					failed(env, "GRBsetdblattrelement(\"ObjN\")");
				}
			if ((t = objconst(i1)) && GRBsetdblattr(mdl, "ObjNCon", t))
				failed(env, "GRBsetdblattr(\"ObjNCon\")");
			}
		if (objpri || objabs || objrel || objwt) {
			for(i = 0; i < no; ++i) {
				if (GRBsetintparam(env, GRB_INT_PAR_OBJNUMBER, i))
					failed(env, "GRBsetintparam(\"ObjNumber\")");
				i1 = operm[i];
				if (objpri && GRBsetintattr(mdl, "ObjNPriority", objpri[i1]))
					failed(env, "GRBsetintattr(\"ObjNPriority\")");
				if (objwt && GRBsetdblattr(mdl, "ObjNWeight", objwt[i1]))
					failed(env, "GRBsetdbkattr(\"ObjNWeight\")");
				if (objabs && GRBsetdblattr(mdl, "ObjNAbsTol", objabs[i1]))
					failed(env, "GRBsetdbkattr(\"ObjNAbsTol\")");
				if (objrel && GRBsetdblattr(mdl, "ObjNRelTol", objrel[i1]))
					failed(env, "GRBsetdbkattr(\"ObjNRelTol\")");
				}
			}
#ifndef NO_MOkwf
		for(pmov = &first_MOkw; (mov = *pmov); pmov = &mov->next) {
			if (!(env1 = GRBgetmultiobjenv(mdl, mov->Objno))) {
				snprintf(mbuf, sizeof(mbuf),
					 "GRBgetmultiobjenv(mdl, %d)", mov->Objno);
				failed(env, mbuf);
				}
			Oinfo.uinfo = (char*)env1;
			kw = mov->kw;
			kw->kf(&Oinfo, kw, mov->val);
			if (Oinfo.n_badopts) {
				badretfmt(501, "\"obj_%d_%s %s\" rejected.\n",
					mov->Objno+1, kw->name, mov->val);
				longjmp(Jb,1);
				}
			}
		if (first_MOkw)
			Oinfo.uinfo = (char*)env;
#endif
		}
#endif /*}*/
#if GRB_VERSION_MAJOR >= 3 /*{*/
	if (mip && ams_stub && (!ams_modeseen || ams_mode == 1 || ams_mode == 2))
	{
		if(!ams_modeseen)
			GRBsetintparam(env, "PoolSearchMode", 2);
#ifdef GRB_INT_ATTR_POOLIGNORE
		if ((lazyd = suf_get("poolignore", ASL_Sufkind_var | ASL_Sufkind_input))
		 && GRBsetintattrarray(mdl, GRB_INT_ATTR_POOLIGNORE, 0, nv, lazyd->u.i))
			failed(env, "GRBsetintattrarray(GRB_INT_ATTR_POOLIGNORE)");
#endif
	}
#endif /*}*/
#if GRB_VERSION_MAJOR >= 5 /*{{*/
	if (warmstart && x) {
		if (basis & 1) {
			basis &= ~1;
			if (get_input_statuses(asl, env, mdl, &dims, resid) && warmstart <= 2)
				warmstart = 0;
			}
		if (warmstart) {
			if (nrange) {
				lu = LUrhs;
				j = nv;
				for(i = 0; i < nc; ++i, lu += 2) {
					if (lu[0] > negInfinity
					 && lu[1] < Infinity
					 && lu[0] < lu[1])
						x[j++] = resid[i] - rhs[i];
					}
				}
			if (warmstart <= 2) {
				if (GRBsetdblattrarray(mdl, GRB_DBL_ATTR_PSTART, 0, nvr, x))
					failed(env, "GRBsetdblattrarray(PStart)");
				if (pi0
				 && nc > nqc
				 && GRBsetdblattrarray(mdl, GRB_DBL_ATTR_DSTART, 0, nc-nqc, pi0))
					failed(env, "GRBsetdblattrarray(DStart)");
				}
#ifdef GRB_DBL_ATTR_VARHINTVAL
			else {
				if (GRBsetdblattrarray(mdl, GRB_DBL_ATTR_VARHINTVAL, 0, nvr, x))
					failed(env, "GRBsetdblattrarray(Varhint)");
				else if (warmstart == 4
					&& (lazyd = suf_get("hintpri", ASL_Sufkind_var | ASL_Sufkind_input))
					&& GRBsetintattrarray(mdl, GRB_INT_ATTR_VARHINTPRI, 0, nv, lazyd->u.i))
					failed(env, "GRBsetdblattrarray(Hintpri)");
				}
#endif
			}
		}
#else /*}{*/
	if (basis & 1)
		get_input_statuses(asl, env, mdl, &dims, resid);
#endif /*}}*/
	free(A);
	if (x && GRBsetdblattrarray(mdl, GRB_DBL_ATTR_START, 0, nv, x))
		failed(env, "GRBsetdblattrarray(START)");
	if (logfile) {
		if (GRBsetintparam(env, "OutputFlag", 1))
			namefailed("GRBsetintparam", "OutputFlag");
#if GRB_VERSION_MAJOR == 1
		if (GRBsetlogfilename(env, logfile))
			failed(env, "GRBsetlogfilename");
#elif GRB_VERSION_MAJOR == 2
		if (GRBsetstrparam(env, "LogfileName", logfile))
			namefailed("GRBsetstrparam", "LogfileName");
#else
		if (GRBsetstrparam(env, GRB_STR_PAR_LOGFILE, logfile))
			namefailed("GRBsetstrparam", GRB_STR_PAR_LOGFILE);
#endif
		}
	else if (outlev && GRBsetintparam(env, "OutputFlag", 1))
		namefailed("GRBsetintparam", "OutputFlag");
	if ((fn = Wflist[0])) {
		GRBupdatemodel(mdl);
		do {
			if (GRBwrite(mdl, fn->name))
				enamefailed(env, "GRBwrite", fn->name);
			} while((fn = fn->next));
		}
#if (GRB_VERSION_MAJOR == 4 && GRB_VERSION_MINOR >= 5) || GRB_VERSION_MAJOR >= 5 /*{*/
	if (rays)
		GRBsetintparam(env, "InfUnbdInfo", 1);
#endif /*}*/
#ifdef GRB_INT_ATTR_IIS_LBFORCE
	if (want_iis)
		force_iis(asl, env, mdl);
#endif
#ifdef GRB_INT_ATTR_BRANCHPRIORITY
	if (mip && priorities)
		add_priorities(asl, env, mdl);
#endif
#if GRB_VERSION_MAJOR >= 8
	if (mip
	&& (partition = suf_get("partition", ASL_Sufkind_var))
	&& (part = partition->u.i)
	&& GRBsetintattrarray(mdl, GRB_INT_ATTR_PARTITION, 0, n_var, part))
		failed(env, "GRBsetintattrarray(\"Partition\")");
#endif
#if GRB_VERSION_MAJOR >= 5
	pfto = 0;
	if (feasrelax)
		pfto = do_feasrelax(asl, env, mdl, &objqual, &fto);
#endif
	breaking = 1;
	Times[1] = xectim_();
#ifdef GRB_INT_PAR_TUNEOUTPUT
	tunemsg = 0;
	if (tunebase && tunerun(asl, env, mdl, &tunemsg)) {
		write_sol(tunemsg, 0, 0, &Oinfo);
		goto done;
		}
#endif
	grbmodel = mdl;
	i = GRBoptimize(mdl);
	grbmodel = 0;
	Times[2] = xectim_();
	// Get the solver message - as a side effect, sets solve_result_num
	solmsg = statmsg(asl, env, mdl, i, &dims, &wantobj);
	dims.mb = mbuf;
	dims.mbend = mbuf + sizeof(mbuf);
	dpf(&dims, "%s: %s", verbuf, solmsg);
	absmipgap = relmipgap = Infinity;
	dims.objprec = 0;
	if (wantobj && !GRBgetdblattr(mdl, GRB_DBL_ATTR_OBJVAL, &obj)) {
		dims.objprec = obj_prec();
#if GRB_VERSION_MAJOR >= 5
		dpf(&dims, "; %sobjective %.*g", objqual, dims.objprec, obj);
		if (pfto)
			dpf(&dims, "\nfeasrelax objective = %.*g", dims.objprec, *pfto);
#else
		dpf(&dims, "; objective %.*g", dims.objprec, obj);
#endif
#if GRB_VERSION_MAJOR >= 7
		if (multiobj && no > 1) {
			dpf(&dims, "\nIndividual objective values:");
			for(i = 0; i < no; ++i) {
				if (GRBsetintparam(env, GRB_INT_PAR_OBJNUMBER, i))
					failed(env, "GRBsetintparam(\"ObjNumber\")");
				t = Infinity;
				if (GRBgetdblattr(mdl, GRB_DBL_ATTR_OBJNVAL, &t))
					failed(env, "GRBgetdblattr(GRB_DBL_ATTR_OBJVAL)");
				dpf(&dims, "\n\t%s = %.*g", obj_name(operm[i]), dims.objprec, t);
				}
			}
#endif
		}
	else
		solnsens = 0;
	if ((bestbound | (retmipgap&3)) && objno > 0) {
		if (GRBgetdblattr(mdl, GRB_DBL_ATTR_OBJBOUND, &f))
			f = Infinity * dims.objsense;
		else if (wantobj) {
			if ((absmipgap = obj - f) < 0.)
				absmipgap = -absmipgap;
			if ((t = obj) < 0.)
				t = -t;
			relmipgap = absmipgap / (1e-10 + t);
			}
		if (retmipgap & 1) {
			bb = (real*)M1zapalloc(nobjno*sizeof(real));
			bb[objno - 1] = relmipgap;
			suf_rput("relmipgap", ASL_Sufkind_obj, bb);
			suf_rput("relmipgap", ASL_Sufkind_prob, bb);
			}
		if (retmipgap & 2) {
			bb = (real*)M1zapalloc(nobjno*sizeof(real));
			bb[objno - 1] = absmipgap;
			suf_rput("absmipgap", ASL_Sufkind_obj, bb);
			suf_rput("absmipgap", ASL_Sufkind_prob, bb);
			}
		if (bestbound) {
			bb = (real*)M1zapalloc(nobjno*sizeof(real));
			bb[objno - 1] = f;
			suf_rput("bestbound", ASL_Sufkind_obj, bb);
			suf_rput("bestbound", ASL_Sufkind_prob, bb);
			}
		}

#if GRB_VERSION_MAJOR >= 8
	if (kappa && !solve_result_num)
	{
		double kvalue;
		if(GRBgetdblattr(mdl, GRB_DBL_ATTR_KAPPA, &kvalue))
			failed(env, "GRBgetdblattr(GRB_DBL_ATTR_KAPPA)");
		if (kappa & 2)
			setSuffixKappa(asl, mdl, kvalue);
		if (kappa & 1)
			dpf(&dims, "\nkappa value: %g", kvalue);
	}
#endif
#if GRB_VERSION_MAJOR >= 3
	if (!GRBgetintattr(mdl, GRB_INT_ATTR_BARITERCOUNT, &i) && i > 0)
		dpf(&dims, "\n%d barrier iterations", i);
#endif
	if (!GRBgetdblattr(mdl, GRB_DBL_ATTR_ITERCOUNT, &f) && f > 0.)
		dpf(&dims, "\n%.0f simplex iterations", f);
	if (mip && !GRBgetdblattr(mdl, GRB_DBL_ATTR_NODECOUNT, &f) && f > 0.)
		dpf(&dims, "\n%.0f branch-and-cut nodes", f);
	mdl1 = mdl;
	mip1 = mip;
	if ((i = solve_result_num)) {
		switch(basisdebug) {
		  case 0:
			if (i < 200)
				break;
		  case 2:
			mip = 0;
		  }
		}
	if (mip && ((basis & 2) | solnsens))
		mdl1 = fixed_model(asl, mdl, &dims);
	if (absmipgap > 0. && absmipgap < Infinity && !(retmipgap & 4))
		dpf(&dims, "\nabsmipgap = %.3g, relmipgap = %.3g", absmipgap, relmipgap);
	if (mdl1) {
		if (solnsens)
			put_sens(asl, mdl1);
		if (basis & 2)
			get_output_statuses(asl, mdl1, &dims);
		if (mdl1 != mdl)
			GRBfreemodel(mdl1);
		}
	if (dims.missing)
		missing_msg(&dims);
#if GRB_VERSION_MAJOR >= 3 /*{*/
	if (ams_stub && wantobj) {
		if (!GRBgetintattr(mdl, GRB_INT_ATTR_SOLCOUNT, &i) && i > 0
		 && (i1 = ams_write(asl, env, mdl, &dims, i, obj, mip1)) > 0) {
			dpf(&dims, "\n%d alternative MIP solutions written to \"%s1.sol\"\n"
				"%s \"%s%d.sol\".", i1, ams_stub,
				i1 == 2 ? "and" : "...", ams_stub, i1);
			if (i1 < i)
				dpf(&dims, "\nIgnoring %d other inferior alternative MIP solutions.",
					i - i1);
			dpf(&dims, "\nAlternative solutions do not include dual variable values.");
			dpf(&dims, "\nBest solution is available in \"%s1.sol\".", ams_stub);
			free(dims.xams);
			}
		else
			i1 = 0;
		suf_iput("npool", ASL_Sufkind_obj, &i1);
		suf_iput("npool", ASL_Sufkind_prob, &i1);
		}
#ifdef GRB_INT_PAR_TUNEOUTPUT
	if (tunemsg)
		dpf(&dims, "\n%s", tunemsg);
#endif
#endif /*}*/
	if (Round && dims.x && (k = solround(asl, dims.x, &t))) {
		if (Round & 2)
			solve_result_num = 3 - (Round & 1);
		if (Round & 4) {
			sc = k > 1 ? "s" : "";
			dpf(&dims, "\n%d integer variable%s %srounded to integer%s;"
			" maxerr = %g", k, sc, Round & 1 ? "" : "would be ", sc, t);
			}
		}
#ifdef GRB_DBL_ATTR_MAX_VIO
	if (want_maxvio) {
		t = reportDblProblemSuffix(asl, mdl, GRB_DBL_ATTR_MAX_VIO, "maxvio");
		dpf(&dims, "\nmax violation = %g\n", t);
		}
#endif
#ifdef GRB_INT_ATTR_CONCURRENTWINMETHOD
	if (want_concurrentwinmethod) {
		i = reportIntProblemSuffix(asl, mdl, GRB_INT_ATTR_CONCURRENTWINMETHOD, "concurrentwinmethod");
		dpf(&dims, "\nconcurrent winner = %d", i);
		}
#endif
#ifdef GRB_DBL_ATTR_WORK
	if (want_work) {
		t = reportDblProblemSuffix(asl, mdl, GRB_DBL_ATTR_WORK, "work");
		dpf(&dims, "\nwork = %g\n", t);
		}
#endif
	write_sol(mbuf, dims.x, dims.y, &Oinfo);
	for(fn = Wflist[1]; fn; fn = fn->next)
		if (GRBwrite(mdl, fn->name))
			enamefailed(env, "GRBwrite", fn->name);
 done:
	if (oic)
		signal(SIGINT, oic);
	if (mdl)
		GRBfreemodel(mdl);
	if (env0)
		GRBfreeenv(env0);
	if (!amplflag && solve_result_num >= 500)
		rc = 1;
	ASL_free(&asl);
	show_times();
	return rc;
	}
