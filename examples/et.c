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

#include "asl_pfgh.h"
#ifdef _ASL_EW_
#define EW(x) x
 typedef real (*Ffunc)(EvalWorkspace*,int nobj,real*X,fint*nerror);
 typedef void (*Gfunc)(EvalWorkspace*,int nobj,real*X,real*G,fint*ne);
#define FirstX ew->x0kind |= ASL_first_x
#undef sputinfo
#define sputinfo ew->Sputinfo
#else
#define EW(x)
 typedef real (*Ffunc)(ASL*,int nobj,real*X,fint*nerror);
 typedef void (*Gfunc)(ASL*,int nobj,real*X,real*G,fint*ne);
#define FirstX x0kind = ASL_first_x
#endif

 typedef unsigned long ULong;

 typedef struct
Pinfo {
	ASL *asl;
	real *p, *v, *w;
	int nx, nxg;
	} Pinfo;

 typedef struct Tcg Tcg;
 struct
Tcg {
	Tcg *next;
	int varno;
	int goff;
	};

typedef const char cchar;
typedef unsigned char uchar;
typedef char *(*Namef)(ASL*,int);

static ASL *asl;
static Pinfo *P, *Pavail, *Pfree, *Plast;
static int Rtime = 1, funcshow = 1, npinfo = 8, nreps = 10000;
static Jmp_buf exit_jb, *mejb;
static real aeps = 0, eps = 1e-8, haeps = 0, heps = 4.8e-6;
static real aextol = -1., extol = 1e-5;

#undef psinfo
#define psinfo ((ASL_pfgh*)asl)->P

 static char *mname[5] = { "f", "fg", "fgh", "pfg", "pfgh" };

#define MaxTemp 8

 static int n_to_free;
 static real *x_to_free[MaxTemp];

 static real *
tempvec(int n)
{
	int n1;

	if (n_to_free >= MaxTemp) {
		fprintf(Stderr, "Bug! MaxTemp exceeded!\n");
		exit(1);
		}
	if (n <= 1)
		n1 = 1;
	else
		for(n1 = 1, --n; n; n >>= 1)
			n1 <<= 1;
	return x_to_free[n_to_free++] = (real *)Malloc(n1*sizeof(real));
	}

 static void
freeall(void)
{
	while(n_to_free > 0)
		free(x_to_free[--n_to_free]);
	}

 void
mainexit_ASL(int n)
{
#ifdef _ASL_EW_
	EvalWorkspace *ew;
	freeall();
	if (asl && (ew = asl->i.Ew0))
		ew->x0kind = ASL_first_x;
#else
	freeall();
	if (asl)
		x0kind = ASL_first_x;
#endif
	longjmp(mejb->jb, n);
	}

 static char *
skipblank(char *s)
{
	while(*(uchar*)s <= ' ' && *s)
		s++;
	return s;
	}

 static char *
get_range(char *s, int *lb, int *ub, int mx)
{
	int i;

	*lb = 0;
	*ub = mx;
	if (!s)
		return s;
	s = skipblank(s);
	if (s[0] >= '0' && s[0] <= '9') {
		i = (int)strtol(s,&s,10);
		if (i > mx)
			i = mx;
		*lb = i;
		if (s[0] == '-')
			goto have_dash;
		if (i < mx)
			*ub = i + 1;
		return s;
		}
	if (s[0] == '-') {
 have_dash:
		i = *++s;
		if (i >= '0' && i <= '9') {
			i = (int)strtol(s,&s,10);
			if (i < mx)
				*ub = i + 1;
			}
		}
	return s;
	}

 static char *
get_prange(char *s0, int *lb, int *ub, int mx, int *p, int *r)
{
	char *s;

	if (s0)
		s0 = skipblank(s0);
	s = get_range(s0, lb, ub, mx);
	*r = s > s0;
	if (s && *(s = skipblank(s)) == 'p') {
		*p = 1;
		s++;
		}
	return s;
	}

 static int
curchk(void)
{
	if (asl)
		return 0;
	printf("# No current problem!\n");
	return 1;
	}

 static char *
cline(char *s)
{
	int i, n;
	real *c, *c0, *ce, *lu, t;

	if (curchk())
		return 0;
	s = get_range(s, &i, &n, n_con);
	if (i >= n)
		return s;
	lu = LUrhs + 2*i;
	if (i == 0 && n == n_con) {
		c = c0 = tempvec(n_con);
		conval(X0,c,0);
		ce = c + n_con;
		if (maxrownamelen)
		    for(; c < ce; c++, lu += 2)
			printf("# %s = %g, lslack = %g, uslack = %g\n",
				con_name((int)(c-c0)),
				*c, *c - lu[0], lu[1] - *c);
		else
		    for(; c < ce; c++, lu += 2)
			printf("# c[%d] = %g, lslack = %g, uslack = %g\n",
				c-c0, *c, *c - lu[0], lu[1] - *c);
		freeall();
		return s;
		}
	xknown(X0);
	if (maxrownamelen)
	    for(; i < n; i++, lu += 2) {
		t = conival(i, X0, 0);
		printf("# %s = %g, lslack = %g, uslack = %g\n",
			con_name(i), t, t - lu[0], lu[1] - t);
		}
	else
	    for(; i < n; i++, lu += 2) {
		t = conival(i, X0, 0);
		printf("# c[%d] = %g, lslack = %g, uslack = %g\n",
			i, t, t - lu[0], lu[1] - t);
		}
	xunknown();
	return s;
	}

#ifdef CUR_GRADIENTS	/* should not be needed after 20010531 */
 static void
cur_gradients(void)
{
	int i, k;
	real *h;

	if (P->nx = P->nxg)
		return;
	k = n_var;
	if (k < nzjac)
		k = nzjac;
	h = tempvec(k);
	for(i = 0; i < n_obj; i++)
		objgrd(i, X0, h, 0);
	if (n_con)
		jacval(X0, h, 0);
	freeall();
	P->nxg = P->nx;
	}
#else
#define cur_gradients()
#endif

 static real
reldif(real a, real b)
{
	real d;

	d = a - b;
	if (d < 0)
		d = -d;
	if (a < 0)
		a = -a;
	if (b < 0)
		b = -b;
	if (a < b)
		a = b;
	if (a > 0)
		d /= a;
	return d;
	}

 static real
fdstep(real x)
{
	if (x > 0)
		return eps * x + aeps;
	if (x < 0)
		return eps * x - aeps;
	return eps + aeps;
	}

 static void
fdgrads(int i, int n, int p, int isobj, Ffunc ffunc, Gfunc gfunc)
{
#ifdef _ASL_EW_
	EvalWorkspace *ew = asl->i.Ew0;
#else
#define ew asl
#endif
	int j, k;
	real dx, ea, ex, f, fdg, fdk, *g, gi, x0;
	static char *what[2] = { "Constraint", "Objective" };

	g = tempvec(n_var);
	for(; i < n; i++) {
		if (maxrownamelen)
			printf("# %s", isobj ? obj_name(i) : con_name(i));
		else
			printf("# %s %d", what[isobj], i);
		printf(p ? ":\n" : ": ");
		f = ffunc(ew, i, X0, 0);
		gfunc(ew, i, X0, g, 0);
		k = -1;
		ex = -1.;
		fdk = 0.;
		for(j = 0; j < n_var; j++) {
			dx = fdstep(x0 = X0[j]);
			X0[j] = x0 + dx;
			fdg = (ffunc(ew, i, X0, 0) - f) / dx;
			X0[j] = x0;
			ea = reldif(fdg, gi = g[j]);
			if (p && ea > extol)
			    printf("# g[%ld] = %g, fd = %g, rel. error = %.2g\n",
				j, gi, fdg, ea);
			if (ex < ea) {
				ex = ea;
				k = j;
				fdk = fdg;
				}
			}
		printf("# max. rel. error = %.2g", ex);
		printf(ex > extol ? " (g[%d] = %g, fd = %g)\n" : "\n",
			k, g[k], fdk);
		}
	freeall();
#undef ew
	}

 static char *
dgline(int p, char *s)
{
	int i, n, r;

	s = get_prange(s, &i, &n, n_obj, &p, &r);
	fdgrads(i, n, p, 1, asl->p.Objval, asl->p.Objgrd);
	return s;
	}

 static real
wf(real *x)
{
	int i;
	real f, t;
	real *w = P->w;

	f = 0;
	for(i = 0; i < n_obj; i++)
		if ((t = *w++))
			f += t*objval(i, x, 0);
	w = pi0;
	for(i = 0; i < n_con; i++)
		if ((t = *w++))
			f += t*conival(i, x, 0);
	return f;
	}

 static void
wg(real *x, real *g, real *gt)
{
	cgrad *cg, **cgp;
	int i, j, k, m;
	real *J, *g1, *gi, *ge;
	real t;
	real *w = P->w;

	memset(g, 0, n_var*sizeof(real));
	ge = g + n_var;
	for(i = 0; i < n_obj; i++)
		if ((t = *w++)) {
			objgrd(i, x, g1 = gt, 0);
			for(gi = g; gi < ge; )
				*gi++ += t * *g1++;
			}
	w = pi0;
	m = n_con;
#if 0
	for(i = 0; i < m; i++) {
		if ((t = *w++)) {
			congrd(i, x, g1 = gt, 0);
			for(gi = g; gi < ge; )
				*gi++ += t * *g1++;
			}
		}
#else
	/* exploit sparsity ==> much faster on large problems */
	j = -1;
	for(i = k = 0; i < m; ++i) {
		if (w[i]) {
			if (++k > 1) {
				J = (real*)Malloc(nzjac*sizeof(real));
				jacval(x, J, 0);
				cgp = Cgrad;
				for(; j < m; ++j) {
					if ((t = w[j])) {
						for(cg = cgp[j]; cg; cg = cg->next)
							g[cg->varno] += t*J[cg->goff];
						}
					}
				free(J);
				return;
				}
			j = i;
			}
		}
	if (j >= 0) {
		t = w[j];
		congrd(j, x, g1 = gt, 0);
		for(gi = g; gi < ge; )
			*gi++ += t * *g1++;
		}
#endif
	}

 static void
fdhstep(real *x, real *fx, real *dx)
{
	int i, n;
	real t, ta;

	n = n_var;
	for(i = 0; i < n; i++) {
		ta = haeps;
		if ((t = x[i]) < 0)
			ta = -ta;
		x[i] += dx[i] = t == 0. ? heps + ta : t*heps + ta;
		fx[i] = wf(x);
		x[i] = t;
		}
	}

 typedef struct
FDHinfo {
	real f, xi, dxi, fxi;
	real *dx, *fx, *x;
	int i;
	} FDHinfo;

 static real
fdh(FDHinfo *f, int j)
{
	real h, rv, xj, *x;

	x = f->x;
	h = f->dx[j];
	xj = x[j];
	if (f->i == j) {
		x[j] = f->xi - h;
		rv = (f->fx[j] + wf(x) - 2*f->f) / (h*h);
		}
	else {
		x[j] = xj + h;
		rv = (wf(x) - f->fxi - f->fx[j] + f->f)/(h*f->dxi);
		}
	x[j] = xj;
	return rv;
	}

#ifdef _ASL_EW_ /*{*/

static void notyet(const char *s)
{ printf("\n%s not yet implemented.\n", s); }

 int
fgh_read_ASL(ASL *a, FILE *f, int flags)
{
	notyet("fgh_read");
	return 0;
	}

 int
pfg_read_ASL(ASL *a, FILE *f, int flags)
{
	notyet("pfg_read");
	return 0;
	}
#endif /*}*/

 static void
dHdline(int p)
{
	FDHinfo F;
	int i, ix, j, jx, n;
	real *H, *Hi;
	real Hx, ea, ex, t, t1, tx;

	n = n_var;
	F.dx = tempvec(((n*(n+1))>>1) + 2*n);
	F.fx = F.dx + n;
	Hi = H = F.fx + n;
	F.f = wf(F.x = X0);
	wg(F.x, F.fx, F.dx);
	duthes(H, -1, P->w, pi0);
	fdhstep(F.x, F.fx, F.dx);
	ex = -1.;
	Hx = tx = 0.;
	ix = jx = -1;
	for(i = 0; i < n; i++) {
		F.i = i;
		F.dxi = F.dx[i];
		F.fxi = F.fx[i];
		F.x[i] = (F.xi = F.x[i]) + F.dxi;
		for(j = 0; j <= i; j++) {
			t = fdh(&F,j);
			t1 = *Hi++;
			if (p)
				printf("# H[%d,%d] = %g, Fd = %g\n", j,i,t1,t);
			ea = reldif(t, t1);
			if (ex < ea) {
				ex = ea;
				ix = j;
				jx = i;
				Hx = t1;
				tx = t;
				}
			}
		F.x[i] = F.xi;
		}
	printf("# duthes: max rel error %.2g at (%d,%d)\n", ex, ix,jx);
	if (ex > extol)
		printf("# Hes(%d,%d) = %g, fd = %g\n", ix,jx,Hx,tx);
	freeall();
	}

 static void
dhdline(int p)
{
	real *g, *gt, *g0, *H, *Hi, *Hi1, *x;
	real Hx, dx, ea, ex, t, t1, tx, xi;
	int i, ix, j, jx, n;

	n = n_var;
	g0 = tempvec(((n*(n+1))>>1) + 3*n);
	g = g0 + n;
	gt = g + n;
	Hi = H = gt + n;
	wg(x = X0, g0, gt);
	duthes(H, -1, P->w, pi0);
	ex = -1.;
	ix = jx = -1;
	Hx = tx = 0.;
	for(i = 0; i < n; i++) {
		x[i] += dx = fdstep(xi = x[i]);
		wg(x,g,gt);
		x[i] = xi;
		for(j = 0; j <= i; j++) {
			t = (g[j] - g0[j]) / dx;
			t1 = *Hi++;
			if (p)
				printf("# H[%d,%d] = %g, fd = %g\n", j,i,t1,t);
			ea = reldif(t, t1);
			if (ex < ea) {
				ex = ea;
				ix = j;
				jx = i;
				Hx = t1;
				tx = t;
				}
			}
		for(Hi1 = Hi-1; j < n; j++) {
			t = (g[j] - g0[j]) / dx;
			t1 = *(Hi1 += j);
			if (p)
				printf("# H[%d,%d] = %g, fd = %g\n", j,i,t1,t);
			ea = reldif(t, t1);
			if (ex < ea) {
				ex = ea;
				ix = j;
				jx = i;
				Hx = t1;
				tx = t;
				}
			}
		}
	printf("# duthes: max rel error %.2g at (%d,%d)\n", ex, ix,jx);
	if (ex > extol)
		printf("# Hes(%d,%d) = %g, fd = %g\n", ix,jx,Hx,tx);
	freeall();
	}

 static void
dHfline(int p)
{
	FDHinfo F;
	int i, ix, j, jx, n;
	fint LH;
	real *H, *Hi;
	real Hx, ea, ex, t, t1, tx;

	LH = n = n_var;

	F.dx = tempvec(n*(n + 2));
	F.fx = F.dx + n;
	Hi = H = F.fx + n;
	F.f = wf(F.x = X0);
	wg(F.x, F.fx, F.dx);
	fullhes(H, LH, -1, P->w, pi0);
	fdhstep(F.x, F.fx, F.dx);
	ex = -1.;
	Hx = tx = 0.;
	ix = jx = -1;
	for(i = 0; i < n; i++) {
		F.i = i;
		F.dxi = F.dx[i];
		F.fxi = F.fx[i];
		F.x[i] = (F.xi = F.x[i]) + F.dxi;
		for(j = 0; j < n; j++) {
			t = fdh(&F,j);
			t1 = *Hi++;
			if (p)
				printf("# H[%d,%d] = %g, Fd = %g\n", j,i,t1,t);
			ea = reldif(t, t1);
			if (ex < ea) {
				ex = ea;
				ix = j;
				jx = i;
				Hx = t1;
				tx = t;
				}
			}
		F.x[i] = F.xi;
		}
	printf("# fullhes: max rel error %.2g at (%d,%d)\n", ex, ix,jx);
	if (ex > extol)
		printf("# Hes(%d,%d) = %g, fd = %g\n", ix,jx,Hx,tx);
	freeall();
	}

 static void
dhfline(int p)
{
	real *g, *gt, *g0, *H, *Hi, *x;
	real Hx, dx, ea, ex, t, t1, tx, xi;
	int i, ix, j, jx, n;
	fint LH;

	LH = n = n_var;

	g0 = tempvec(n*(n + 3));
	g = g0 + n;
	gt = g + n;
	Hi = H = gt + n;
	wg(x = X0, g0, gt);
	fullhes(H, LH, -1, P->w, pi0);
	ex = -1.;
	ix = jx = -1;
	Hx = tx = 0.;
	for(i = 0; i < n; i++) {
		x[i] += dx = fdstep(xi = x[i]);
		wg(x,g,gt);
		x[i] = xi;
		for(j = 0; j < n; j++) {
			t = (g[j] - g0[j]) / dx;
			t1 = *Hi++;
			if (p)
				printf("# H[%d,%d] = %g, fd = %g\n", j,i,t1,t);
			ea = reldif(t, t1);
			if (ex < ea) {
				ex = ea;
				ix = j;
				jx = i;
				Hx = t1;
				tx = t;
				}
			}
		}
	printf("# fullhes: max rel error %.2g at (%d,%d)\n", ex, ix,jx);
	if (ex > extol)
		printf("# Hes(%d,%d) = %g, fd = %g\n", ix,jx,Hx,tx);
	freeall();
	}

 static void
dHsline(int p, int uptri)
{
	FDHinfo F;
	int i, ib, ix, ixa, j, jx, jxa, n;
	fint *hcs, *hr, *hr0, *hre, nsph;
	real *H, *Hi;
	real Hx, Hxa, Hxd, ea, ex, exa, t, t1, ta, tx;
	EW(EvalWorkspace *ew = asl->i.Ew0;)

	n = n_var;
	nsph = sphsetup(-1,1,1,uptri);
	hr = hr0 = sputinfo->hrownos;
	hcs = sputinfo->hcolstarts;
	F.dx = tempvec(nsph + 2*n);
	F.fx = F.dx + n;
	Hi = H = F.fx + n;
	F.f = wf(F.x = X0);
	wg(F.x, F.fx, F.dx);
	sphes(H, -1, P->w, pi0);
	fdhstep(F.x, F.fx, F.dx);
	ex = exa = -1.;
	ix = ixa = jx = jxa = -1;
	ib = n-1;
	Hxa = Hxd = Hx = tx = 0.;
	for(i = j = 0; i < n; i++) {
		F.i = i;
		F.dxi = F.dx[i];
		F.fxi = F.fx[i];
		F.x[i] = (F.xi = F.x[i]) + F.dxi;
		hre = hr0 + *++hcs;
		if (uptri)
			ib = i;
		for(j = 0; j <= ib; j++) {
			t1 = 0;
			if (hr < hre && *hr == j) {
				t1 = *Hi++;
				hr++;
				}
			t = fdh(&F,j);
			if (p)
				printf("# H[%d,%d] = %g, Fd = %g\n", j,i,t1,t);
			ea = reldif(t, t1);
			if (ex < ea) {
				ex = ea;
				ix = j;
				jx = i;
				Hx = t1;
				tx = t;
				}
			if ((ta = t - t1) < 0)
				ta = -ta;
			if (exa < ta) {
				exa = ta;
				ixa = j;
				jxa = i;
				Hxd = t;
				Hxa = t1;
				}
			}
		F.x[i] = F.xi;
		}
	printf("# sputhes: max rel error %.2g at (%d,%d)\n", ex, ix,jx);
	if (ex > extol) {
		printf("# Hes(%d,%d) = %g, Fd = %g\n", ix,jx,Hx,tx);
		if (exa > aextol && aextol >= 0 && (ix != ixa || jx != jxa))
			printf("# Hes(%d,%d) = %g, fd = %g, diff = %g\n",
				ixa, jxa, Hxa, Hxd, Hxa-Hxd);
		}
	freeall();
	}

 static void
dhsline(int p, int uptri)
{
	EW(EvalWorkspace *ew = asl->i.Ew0;)
	real *g, *gt, *g0, *H, *Hi, *x;
	real Hx, Hxa, Hxd, dx, ea, ex, exa, t, t1, ta, tx, xi;
	int i, ib, ix, ixa, j, jx, jxa, n;
	fint *hcs, *hr, *hr0, *hre, nsph;

	n = n_var;
	nsph = sphsetup(-1,1,1,uptri);
	hr = hr0 = sputinfo->hrownos;
	hcs = sputinfo->hcolstarts;
	g0 = tempvec(nsph + 3*n);
	g = g0 + n;
	gt = g + n;
	Hi = H = gt + n;
	wg(x = X0, g0, gt);
	sphes(H, -1, P->w, pi0);
	ex = exa = -1.;
	ix = ixa = jx = jxa = -1;
	ib = n-1;
	Hx = Hxa = Hxd = tx = 0.;
	for(i = j = 0; i < n; i++) {
		x[i] += dx = fdstep(xi = x[i]);
		wg(x,g,gt);
		x[i] = xi;
		hre = hr0 + *++hcs;
		if (uptri)
			ib = i;
		for(j = 0; j <= ib; j++) {
			t1 = 0;
			if (hr < hre && *hr == j) {
				t1 = *Hi++;
				hr++;
				}
			t = (g[j] - g0[j]) / dx;
			if (p)
				printf("# H[%d,%d] = %g, fd = %g\n", j,i,t1,t);
			ea = reldif(t, t1);
			if (ex < ea) {
				ex = ea;
				ix = j;
				jx = i;
				Hx = t1;
				tx = t;
				}
			if ((ta = t - t1) < 0)
				ta = -ta;
			if (exa < ta) {
				exa = ta;
				ixa = j;
				jxa = i;
				Hxd = t;
				Hxa = t1;
				}
			}
		}
	printf("# sputhes: max rel error %.2g at (%d,%d)\n", ex, ix,jx);
	if (ex > extol) {
		printf("# Hes(%d,%d) = %g, fd = %g\n", ix,jx,Hx,tx);
		if (exa > aextol && aextol >= 0 && (ix != ixa || jx != jxa))
			printf("# Hes(%d,%d) = %g, fd = %g, diff = %g\n",
				ixa, jxa, Hxa, Hxd, Hxa-Hxd);
		}
	freeall();
	}

 static void
dHvline(int p)
{
	FDHinfo F;
	int i, ix, j, jx, n;
	real *g, *gt, *hv, *v;
	real Hx, ea, ex, t, t1, tx;

	n = n_var;

	F.dx = tempvec(6*n);
	F.fx = F.dx + n;
	v = F.fx + n;
	hv = v + n;
	g = hv + n;
	gt = g + n;
	memset(v, 0, n_var*sizeof(real));
	fdhstep(F.x = X0, F.fx, F.dx);
	F.f = wf(F.x);
	ex = -1.;
	tx = Hx = 0.;
	ix = jx = -1;
	for(i = 0; i < n; i++) {
		v[i] = 1.;
		wg(F.x, g, gt);
		hvcomp(hv, v, -1, P->w, pi0);
		v[i] = 0.;
		F.i = i;
		F.dxi = F.dx[i];
		F.fxi = F.fx[i];
		F.x[i] = (F.xi = F.x[i]) + F.dxi;
		for(j = 0; j < n; j++) {
			t = fdh(&F,j);
			t1 = hv[j];
			if (p)
				printf("# H[%d,%d] = %g, Fd = %g\n", j,i,t1,t);
			ea = reldif(t, t1);
			if (ex < ea) {
				ex = ea;
				ix = j;
				jx = i;
				Hx = t1;
				tx = t;
				}
			}
		F.x[i] = F.xi;
		}
	printf("# Hvcomp: max rel error %.2g at (%d,%d)\n", ex, ix,jx);
	if (ex > extol)
		printf("# Hes(%d,%d) = %g, fd = %g\n", ix,jx,Hx,tx);
	freeall();
	}

 static void
dhvline(int p)
{
	real *g, *gt, *g0, *hv, *v, *x;
	real Hx, dx, ea, ex, t, t1, tx, xi;
	int i, ix, j, jx, n;

	n = n_var;

	g0 = tempvec(5*n);
	g = g0 + n;
	gt = g + n;
	v = gt + n;
	hv = v + n;
	memset(v, 0, n_var*sizeof(real));
	x = X0;
	ex = -1.;
	ix = jx = -1;
	Hx = tx = 0.;
	for(i = 0; i < n; i++) {
		v[i] = 1.;
		wg(x, g0, gt);
		hvcomp(hv, v, -1, P->w, pi0);
		v[i] = 0.;
		x[i] += dx = fdstep(xi = x[i]);
		wg(x,g,gt);
		x[i] = xi;
		for(j = 0; j < n; j++) {
			t = (g[j] - g0[j]) / dx;
			t1 = hv[j];
			if (p)
				printf("# H[%d,%d] = %g, fd = %g\n", j,i,t1,t);
			ea = reldif(t, t1);
			if (ex < ea) {
				ex = ea;
				ix = j;
				jx = i;
				Hx = t1;
				tx = t;
				}
			}
		}
	printf("# Hvcomp: max rel error %.2g at (%d,%d)\n", ex, ix,jx);
	if (ex > extol)
		printf("# Hes(%d,%d) = %g, fd = %g\n", ix,jx,Hx,tx);
	freeall();
	}

 typedef struct
Dhl {
	void (*DHD)(int);
	void (*DHF)(int);
	void (*DHS)(int,int);
	void (*DHV)(int);
	}Dhl;

 static char *
dhline(int p, char *s, Dhl *dhl) {
	int k = 2;
	int uptri = 1;
	if (s) {
		while(*s <= ' ' && *s)
			++s;
		switch(*s++) {
		 case 'b': k = 2; uptri = 0; break;
		 case 'd': k = 0; break;
		 case 'f': k = 1; break;
		 case 's': k = 2; break;
		 case 'v': k = 3; break;
		 default: --s;
		 }
		s = skipblank(s);
		if (*s == 'p') {
			p = 1;
			s = skipblank(s+1);
			}
		}
	switch(k) {
	 case 0: dhl->DHD(p);		break;
	 case 1: dhl->DHF(p);		break;
	 case 2: dhl->DHS(p,uptri);	break;
	 case 3: dhl->DHV(p);
	 }
	return s;
	}

 static void
transpose(int m, cgrad **cgp, int n, Tcg **tcp, Tcg *tcfree)
{
	cgrad *cg, **cgpi;
	Tcg *tcg, **tcpi;

	tcpi = tcp + n;
	while(tcpi > tcp)
		*--tcpi = 0;
	cgpi = cgp + m;
	while(cgpi > cgp)
		for(--m, cg = *--cgpi; cg; cg = cg->next) {
			tcpi = tcp + cg->varno;
			tcg = tcfree++;
			tcg->varno = m;
			tcg->goff = cg->goff;
			tcg->next = *tcpi;
			*tcpi = tcg;
			}
	}

 static void
fdjac(int p)
{
	int i, j, k, L, m, n;
	real ag, dx, ea, ex, fdg, gi, *J, *R, *R1, t, x0;
	Tcg *tcfree, *tcg, **tcgp;

	m = n_con;
	n = n_var;
	R = tempvec(2*m + nzjac
		+ (nzjac*sizeof(Tcg) + n*sizeof(Tcg*))/sizeof(real) + 1);
	R1 = R + m;
	J = R1 + m;
	tcfree = (Tcg*)(J + nzjac);
	tcgp = (Tcg**)(tcfree + nzjac);

	conval(X0, R, 0);
	jacval(X0, J, 0);
	transpose(m, Cgrad, n, tcgp, tcfree);
	i = 0;
	j = k = -1;
	ex = -1.;
	ag = fdg = 0.;
	for(i = 0; i < n; i++) {
		dx = fdstep(x0 = X0[i]);
		X0[i] = x0 + dx;
		conval(X0, R1, 0);
		X0[i] = x0;
		tcg = *tcgp++;
		for(L = 0; L < m; L++) {
			t = (R1[L] - R[L]) / dx;
			gi = 0.;
			if (tcg && tcg->varno == L) {
				gi = J[tcg->goff];
				tcg = tcg->next;
				}
			ea = reldif(t, gi);
			if (p && ea > extol)
		    		printf("# J[%d,%d] = %g, fd = %g, rel. error = %.2g\n",
					L, i, gi, t, ea);
			if (ex < ea) {
				ex = ea;
				j = L;
				k = i;
				fdg = t;
				ag = gi;
				}
			}
		}
	printf("# max. rel. Jacobian error = %.2g", ex);
	printf(ex > extol ? " (Jacobian = %g, fd = %g at J[%d,%d])\n" : "\n",
		ag, fdg, j, k);
	freeall();
	}

 static char *
djline(int p, char *s)
{
	int i, n, r;

	s = get_prange(s, &i, &n, n_con, &p, &r);
	if (r)
		fdgrads(i, n, p, 0, asl->p.Conival, asl->p.Congrd);
	else if (n_con > 0)
		fdjac(p);
	return s;
	}

 static char *
dline(char *s) {
	static Dhl dhfd = { dHdline, dHfline, dHsline, dHvline };
	static Dhl dhgd = { dhdline, dhfline, dhsline, dhvline };
	int i, p = 0;

	if (curchk())
		return 0;
	if (*s == 'p') {
		s = skipblank(s+1);
		p = 1;
		}
	switch(*s++) {
	  case 'H': s = dhline(p,s,&dhfd);	break;
	  case 'g': s = dgline(p,s);		break;
	  case 'h': s = dhline(p,s,&dhgd);	break;
	  case 'j': s = djline(p,s);		break;
	  default:
		--s;
		i = asl->i.ASLtype;
		if (i >= ASL_read_fg) {
			dgline(p,0);
			djline(p,0);
			if (i == ASL_read_pfgh)
				dhline(p,0,&dhgd);
			}
	  }
	return s;
	}

 static int
addind(void *V, int iv, int compl, int sense, int nz, int *ig, real *g, real rhs)
{
	int i;
	printf("\nIndicator constraint %d: iv = %d, compl = %d, sense = %d, rhs = %.g, nz = %d:\n",
		++*(int*)V, iv, compl, sense, rhs, nz);
	for(i = 0; i < nz; ++i)
		printf("\tg[%d] = %.g\n", ig[i], g[i]);
	return 0;
	}

 static void
report_indicators(void)
{
	int ei[2], i, ni, nl;

	nl = n_lcon;
	if (nl <= 0) {
		printf("No indicator constraints\n");
		return;
		}
	ei[0] = ei[1] = 0;
	ni = 0;
	if ((i = indicator_constrs_ASL(asl, &ni, addind, ei)))
		printf("return %d from indicator_constrs_ASL(); ei = %d, %d\n", i, ei[0], ei[1]);
	}

 static char *
Lline(char *s) {
	int i, n;

	if (curchk())
		return 0;
	s = skipblank(s);
	if (*s == 'i') {
		report_indicators();
		return s + 1;
		}
	s = get_range(s, &i, &n, n_lcon);
	if (i >= n)
		return s;
	if (maxrownamelen)
		for(; i < n; i++)
			printf("# %s = %s\n", lcon_name(i), lconval(i, X0, 0) ? "true" : "false");
	else
		for(; i < n; i++)
			printf("# L[%d] = %s\n", i, lconval(i, X0, 0) ? "true" : "false");
	return s;
	}

 static char *
fline(char *s) {
	int i, n;

	if (curchk())
		return 0;
	s = get_range(s, &i, &n, n_obj);
	if (i >= n)
		return s;
	if (maxrownamelen)
		for(; i < n; i++)
			printf("# %s = %g\n", obj_name(i), objval(i, X0, 0));
	else
		for(; i < n; i++)
			printf("# f[%d] = %g\n", i, objval(i, X0, 0));
	return s;
	}

 static char *
gline(char *s) {
	int i, n;
	real *g, *ge, *gi;

	if (curchk())
		return 0;
	s = get_range(s, &i, &n, n_obj);
	if (i >= n)
		return s;
	g = tempvec(n_var);
	ge = g + n_var;
	for(; i < n; i++) {
		printf("\n# Gradient of objective %d:\n", i);
		objgrd(i, X0, g, 0);
		gi = g;
		if (maxcolnamelen)
			for(; gi < ge; gi++)
				printf("# g[%s] = %g\n",
					var_name((int)(gi-g)), *gi);
		else
			for(; gi < ge; gi++)
				printf("# g[%d] = %g\n", gi-g, *gi);
		}
	freeall();
	return s;
	}

 static void
v_init(void)
{
	int nv;
	real *v, *ve;

	nv = n_var;
	P->v = v = (real*)M1alloc(nv*sizeof(real));
	ve = v + nv;
	while(v < ve)
		*v++ = 1.;
	}

 static char *
hline(char *s) {
	EW(EvalWorkspace *ew;)
	int i, j, k, nz, uptri;
	real *h, *h1, *he, t;
	fint N, *hcs, *hr, nsph;

	if (curchk())
		return 0;
	cur_gradients();
	EW(if (!(ew = asl->i.Ew0))
		asl->i.Ew0 = ew = ewalloc2_ASL(asl);)
	k = n_var;
	switch(*s++) {
	  case 'd':
		duthes(h = tempvec((k*(k+1))>>1), -1, P->w, pi0);
		j = 0;
		if (maxcolnamelen)
		    for(; j < k; j++) {
			for(i = 0; i <= j; i++)
			    if ((t = *h++))
				printf("# h[%s,%s] = %g\n", var_name(i),
					var_name(j), t);
			}
		else
		    for(; j < k; j++)
			for(i = 0; i <= j; i++)
			    if ((t = *h++))
				printf("# h[%d,%d] = %g\n", i, j, t);
		break;

	  case 'f':
		N = k;
		fullhes(h = tempvec(k*k), N, -1, P->w, pi0);
		j = 0;
		if (maxcolnamelen)
		    for(; j < k; j++) {
			for(i = 0; i < k; i++)
			    if ((t = *h++))
				printf("# h[%s,%s] = %g\n",
					var_name(i), var_name(j), t);
			}
		else
		    for(; j < k; j++)
			for(i = 0; i < k; i++)
			    if ((t = *h++))
				printf("# h[%d,%d] = %g\n", i, j, t);
		break;

	  case 'v':
		if (!P->v)
			v_init();
		hvcomp(h = tempvec(k), P->v, -1, P->w, pi0);
		printf("# h = H*v:\n");
		for(i = 0; i < k; ++i) {
			if (h[i])
				printf("# h[%d] = %g\n", i, h[i]);
			}
		break;

	  default:
		--s;
		/* no break; */
	  case 's':
		uptri = 1;
		goto spcase;
	  case 'b':
		uptri = 0;
		goto spcase;
	  case 'z':
		nz = 1;
		goto spcase1;
	  case 'l':
		uptri = 2;
	  spcase:
		nz = 0;
		if (*s == 'z') {
			nz = 1;
			++s;
			}
 spcase1:
		nsph = sphsetup(-1,1,1,uptri);
		sphes(h = tempvec(nsph), -1, P->w, pi0);
		h1 = h;
		hcs = sputinfo->hcolstarts;
		hr = sputinfo->hrownos;
		i = 0;
		if (maxcolnamelen)
		    for(; i < k; i++) {
			he = h + *++hcs;
			while(h1 < he) {
				t = *h1++;
				if (!nz || t != 0.)
					printf("# h[%s,%s] = %g\n",
						var_name((int)*hr++),
						var_name(i), t);
				}
			}
		else
		    for(; i < k; i++) {
			he = h + *++hcs;
			while(h1 < he) {
				t = *h1++;
				if (!nz || t != 0.)
					printf("# h[%ld,%d] = %g\n", *hr++, i, t);
				}
			}
	  }
	freeall();
	return s;
	}

 static char *
jline(char *s) {
	int i, n;
	real *J, *J1, *Je;
	cgrad *cg, **cgp;

	if (curchk())
		return 0;
	s = get_range(s, &i, &n, n_con);
	if (i >= n)
		return s;
	if (i == 0 && n == n_con) {
		J = tempvec(nzjac);
		jacval(X0, J, 0);
		cgp = Cgrad;
		i = 0;
		if (maxcolnamelen && maxrownamelen)
		    for(; i < n_con; i++)
			for(cg = *cgp++; cg; cg = cg->next)
				printf("# J[%s,%s] = %g\n", con_name(i),
					var_name(cg->varno), J[cg->goff]);
		else
		    for(; i < n_con; i++)
			for(cg = *cgp++; cg; cg = cg->next)
				printf("# J[%d,%d] = %g\n", i, cg->varno,
					J[cg->goff]);
		goto done;
		}
	J = tempvec(n_var);
	Je = J + n_var;
	xknown(X0);
	if (maxcolnamelen)
	    for(; i < n; i++) {
		congrd(i, X0, J, 0);
		for(J1 = J; J1 < Je; J1++)
			if (*J1)
				printf("# J[%s,%s] = %g\n", con_name(i),
					var_name((int)(J1-J)), *J1);
		}
	else
	    for(; i < n; i++) {
		congrd(i, X0, J, 0);
		for(J1 = J; J1 < Je; J1++)
			if (*J1)
				printf("# J[%d,%d] = %g\n", i, (int)(J1-J), *J1);
		}
	xunknown();
 done:
	freeall();
	return s;
	}

 static void
show_curno(void)
{
	if (asl)
		printf("# Current problem number = %d, mode %d (%s), file %.*s.nl\n",
			(int)(P - Pavail), asl->i.ASLtype,
			mname[asl->i.ASLtype-1], (int)(stub_end - filename), filename);
	else
		printf("# No current problem\n");
	}

 static char *
getprobno(char *s)
{
	int c1 = *s, i;

	if (c1 >= '0' && c1 <= '9') {
		i = (int)strtol(s,&s,10);
		if (i < npinfo && Pavail[i].asl) {
			P = Pavail + i;
			asl = P->asl;
			}
		}
	return s;
	}

 static char *
pline(char *s)
{
	s = getprobno(s);
	show_curno();
	return s;
	}

 static char *
sline(char *s)
{
	if (s)
		s = getprobno(s);
	show_curno();
	if (asl) {
		if (n_lcon)
			printf("# %d variables, %d constraints (including %d logical), %d objectives\n",
				n_var, nclcon, n_lcon, n_obj);
		else
			printf("# %d variables, %d constraints, %d objectives\n",
				n_var, n_con, n_obj);
		}
	return s;
	}

 static char *
rline(char *s)
{
	char *stub;
	fint stublen;
	int b, flags, mode;
#ifdef _ASL_EW_
#define rmode mode
#define smax '5'
#else
#define smax '6'
	int rmode;
#endif
	FILE *f;
	real t, t0, *w, *we;
	ASL *asl1;

	if (!*s)
		return s;
	for(stub = s; *++s > ' '; );
	stublen = s - stub;
	mode = 5;
	flags = ASL_allow_CLP;
	if (*s) {
		*s++ = 0;
		s = skipblank(s);
		if (s[0] >= '1' && s[0] <= smax) {
			mode = *s++ - '0';
			s = skipblank(s);
			}
		if (mode >= 4 && mode != 6)
			flags |= ASL_findgroups;
		if (*s >= '1' && *s <= '9')
			flags = (int)strtol(s,&s,10);
		else if (*s == '#')
			flags = (int)strtol(s+1,&s,16);
		else if (*s == '0') {
			b = 8;
			if (s[1] == 'x' || s[1] == 'X') {
				s += 2;
				b = 16;
				}
			flags = (int)strtol(s,&s,b);
			}
		}
#ifndef _ASL_EW_
	if ((rmode = mode) == 6)
		mode = 2;
#endif
	asl1 = ASL_alloc(mode);
	/* Use asl1 rather than asl in case of longjmp!!! */
#define asl asl1
	return_nofile = 1;
	f = jac0dim(stub, stublen);
	if (!f) {
		printf("# Can't open %s\n", stub);
		return 0;
		}
	want_xpi0 = 27;
	t0 = xectim_();
	switch(rmode) {
	 case 1: f_read_ASL   (asl, f, flags);	break;
	 case 2: fg_read_ASL  (asl, f, flags);	break;
	 case 3: fgh_read_ASL (asl, f, flags);	break;
	 case 4: pfg_read_ASL (asl, f, flags);	break;
	 case 5: pfgh_read_ASL(asl, f, flags);  break;
#ifndef _ASL_EW_
	 case 6: qp_read_ASL  (asl, f, flags);
#endif
	 }
	if (!X0)
		X0 = (real*)M1zapalloc(n_var*sizeof(real));
	t = xectim_();
	if (mode > 1 && funcshow) {
		show_funcs();
		funcshow = 0;
		}
#undef asl
	if (Pfree == Plast) {
		Pavail = Realloc(Pavail, 2*sizeof(Pinfo)*npinfo);
		memset(Pavail + npinfo, 0, npinfo*sizeof(Pinfo));
		npinfo <<= 1;
		Plast = Pavail + npinfo;
		}
	P = Pfree++;
	while(Pfree < Plast && Pfree->asl)
		Pfree++;
	P->asl = asl = asl1;
	P->nx = 1;
	P->nxg = 0;
	P->v = 0;
	P->w = 0;
	if (n_obj > 0) {
		w = P->w = (real*)M1alloc(n_obj*sizeof(real));
		for(we = w + n_obj; w < we; w++)
			*w = 1.;
		}
	P->p = 0;
	sline(0);	/* show stats */
	if (Rtime)
		printf("## Read time %g seconds\n", t-t0);
	return s;
	}

 static char *
uline(char *s) {
	Pinfo *P0 = P;
	s = getprobno(s);
	if (asl) {
		if (P->p) {
			free(P->p);
			P->p = 0;
			}
		ASL_free(&P->asl);
		if (Pfree > P)
			Pfree = P;
		if (P != P0) {
			P = P0;
			asl = P->asl;
			}
		else
			for(P = Pavail; ; P++) {
				if (P == Plast) {
					P = 0;
					break;
					}
				if ((asl = P->asl))
					break;
				}
		}
	sline(0);
	return s;
	}

#ifndef NO_YLINES
 static int iran(int);
 static int yseed = 1;
#include <time.h>
#endif

 static int slen(char *s)
{
	char *s1;
	for(s1 = s; *(uchar*)s1 > ' '; ++s1);
	return (int)(s1 - s);
	}

 static char *
wxyline(char *s, real *x, cchar *name, int nx, int *ncp, Namef nf, real *scale)
{
	char *s1;
	int i, j, j0, j1, m, mx, n, nc, seed;
	real L, U, t, t1;

	nc = 0;
	j = j1 = -1;
	i = *s;
	switch(i) {
	  case '*':
		m = n = nx;
		i = j0 = 0;
		s1 = s;
		goto get_t;
	  case 'p':
		++s;
		break;
	  case ',':
		++s;
	  case 0:
		i = 'p';
		break;
#ifndef NO_YLINES /*{*/
	  case 'R':
 Rval:
		L = 0.;
		U = 1.;
		if ((s = skipblank(s+1)) && *s) {
			L = strtod(s, &s1);
			if (s1 <= s || *(uchar*)s1 > ' ') {
				printf("Bad lower bound \"%.*s\" for random numbers.\n",
					slen(s), s);
				return 0;
				}
			if (!(s = skipblank(s1))) {
				printf("Missing upper bound for random numbers.\n");
				return 0;
				}
			U = strtod(s, &s1);
			if (s1 <= s || *(uchar*)s1 > ' ') {
				printf("Bad upper bound \"%.*s\" for random numbers.\n",
					slen(s), s);
				return 0;
				}
			if (U <= L) {
				printf("Bad upper bound %.g for random numbers:\n"
					"\tmust be greater than lower bound %.g.\n",
					U, L);
				return 0;
				}
			U -= L;
			s = s1;
			}
		seed = yseed;
		for(i = 0; i < nx; ++i)
			x[i] = L + U*(seed = iran(seed)) / 2147483647.;
		if (scale)
			for(i = 0; i < nx; ++i)
				x[i] /= scale[i];
		return s;
#endif /*}*/
		}
	if (i == 'p') {
		s = get_range(skipblank(s), &i, &n, nx);
		if (scale) {
			if (maxrownamelen)
			    for(; i < n; i++) {
				t = x[i];
				U = scale[i];
				printf("# %s[%d] = \"%s\" = %g/%g = %g\n",
					name, i, nf(asl,i), t*U, U, t);
				}
			else
			    for(; i < n; i++) {
				t = x[i];
				U = scale[i];
				printf("# %s[%d] = %g/%g = %g\n", name, i, t*U, U, t);
				}
			goto done;
			}
		if (maxrownamelen)
			for(; i < n; i++)
				printf("# %s[%d] = \"%s\" = %g\n", name, i, nf(asl,i), x[i]);
		else
			for(; i < n; i++)
				printf("# %s[%d] = %g\n", name, i, x[i]);
		goto done;
		}
	s = get_range(skipblank(s), &i, &n, nx);
	j0 = i;
	while(i < n) {
		t = strtod(s,&s1);
		if (s == s1) {
			if (*(s1 = skipblank(s)) != '*')
				break;
			m = n - i;
			goto get_t;
			}
		m = 1;
		if (*s1 == '*') {
			m = (int)t;
			if (m > (mx = n-i))
				m = mx;
 get_t:
			t = strtod(s = ++s1, &s1);
			if (s == s1) {
				while(*s && *s <= ' ')
					++s;
				if (*s == 'R')
					goto Rval;
				break;
				}
			}
		s = s1;
		s = s1;
		if (*s == ',')
			++s;
		while(m-- > 0) {
			if (j >= 0) {
				if (j >= j1)
					j = j0;
				t = x[j++];
				}
			U = scale ? t / scale[i] : t;
			if (x[i] != U) {
				x[i] = U;
				nc++;
				}
			i++;
			if (m > 1) {
				t1 = strtod(s, &s1);
				if (s != s1) {
					t = t1;
					j1 = i + 1;
					if (*(s = s1) == ',')
						++s;
					}
				}
			else if (j1 > 0 && j < 0)
				j = j0;
			}
		}
	if (nc && ncp)
		++*ncp;
 done:
	return s;
	}

 static char *
vline(char *s) {
	if (curchk())
		return 0;
	if (!P->v)
		v_init();
	return wxyline(s, P->v, "v", n_var, 0, var_name_ASL, 0);
	}

 static char *
wline(char *s) {
	if (curchk())
		return 0;
	return wxyline(s, P->w, "w", n_obj, 0, obj_name_ASL, 0);
	}

 static char *
xline(char *s) {
	if (curchk())
		return 0;
	return wxyline(s, X0, "x", n_var, &P->nx, var_name_ASL, asl->i.vscale);
	}

 static char *
yline(char *s)
{
	if (curchk())
		return 0;
	return wxyline(s, pi0, "y", n_con, 0, con_name_ASL, asl->i.cscale);
	}

 static void
scalevec(void(*f)(ASL*,int,real,fint*), int k, int n)
{
	real a, b;
	int i = 0;

	switch(k) {
	  case 0:
		a = 0;
		b = 1.;
		while(i < n) {
			b = -b;
			a += 10;
			(*f)(asl, i++, a * b, 0);
			}
		break;
	  case 1:
		while(i < n)
			(*f)(asl, i++, 10., 0);
		break;
	  case 2:
		while(i < n)
			(*f)(asl, i++, -10., 0);
	  }
	}

 static char *
Sline(char *s)
{
	int i;
	if (curchk())
		return 0;
	i = (int)strtol(s,&s,10);
	if (n_con > 0)
	    switch(i) {
		case 1:	scalevec(conscale_ASL,0,n_con); break;
		case 2: scalevec(varscale_ASL,0,n_var); break;
		case 3:	scalevec(conscale_ASL,1,n_con); break;
		case 4: scalevec(varscale_ASL,1,n_var); break;
		case 5:	scalevec(conscale_ASL,2,n_con); break;
		case 6: scalevec(varscale_ASL,2,n_var); break;
		}
	return s;
	}

 static char *
Vline(char *s)
{
	int b, c, i, j, n;
	real bj, *lu, *lue, t, tj, tm, *x, xi;

	if (curchk())
		return s;
	b = j = 0;
	tm = Infinity;
	x = X0;
	if (*s == 'v') {
		c = 'v';
		s = get_range(s+1, &i, &n, n_var);
		lu = LUv;
		lue = lu + 2*n;
		x += i;
		for(lu += 2*i; lu < lue; lu += 2, ++i) {
			xi = *x++;
			t = xi - lu[0];
			if (tm > t) {
				tm = t;
				tj = xi;
				bj = lu[0];
				j = i;
				b = 'L';
				}
			else {
				t = lu[1] - xi;
				if (tm > t) {
					tm = t;
					tj = xi;
					bj = lu[1];
					j = i;
					b = 'U';
					}
				}
			}
		}
	else {
		i = 0;
		n = n_con;
		c = 'c';
		if (*s) {
			if (*s == 'c')
				++s;
			s = get_range(s, &i, &n, n);
			}
		lu = LUrhs;
		lue = lu + 2*n;
		xknown(x);
		for(lu += 2*i; lu < lue; lu += 2, ++i) {
			xi = conival(i, x, 0);
			t = xi - lu[0];
			if (tm > t) {
				tm = t;
				tj = xi;
				bj = lu[0];
				j = i;
				b = 'L';
				}
			else {
				t = lu[1] - xi;
				if (tm > t) {
					tm = t;
					tj = xi;
					bj = lu[1];
					j = i;
					b = 'U';
					}
				}
			}
		xunknown();
		}
	printf("Minimum %c slack = %g", c, tm);
	if (b)
		printf(": %c[%d] = %g, %c[%d] = %g\n", c,j,tj,b,j,bj);
	else
		putchar('\n');
	return s;
	}

 static char *
Xline(char *s)
{
	char buf[1024], *solmsg, *t, *te;
	int m;
	real *x, *y;

	if (curchk())
		return s;
	if (*(uchar*)s > ' ') {
		t = buf;
		for(te = t + sizeof(buf) - 1; t < te && (*t = *(uchar*)s) > ' '; ++s, ++t);
		*t = 0;
		solmsg = fread_soln(buf, &x, &y);
		s = skipblank(s);
		}
	else
		solmsg = read_soln(&x, &y);
	if (solmsg) {
		printf("\n%s\n", solmsg);
		if (x) {
			memcpy(X0, x, n_var*sizeof(real));
			free(x);
			}
		else
			printf("No primal solution found in \"%s\"\n", buf);
		if ((m = n_con)) {
			if (y) {
				memcpy(pi0, y, m*sizeof(real));
				free(y);
				}
			else
				printf("No dual solution found in \"%s\"\n", buf);
			}
		}
	return s;
	}

 static char *
lline(char *s)
{
	lagscale(strtod(s,&s), 0);
	return s;
	}

 static char *
oline(char *s)
{
	int i, n;

	if (curchk())
		return 0;
	s = get_range(s, &i, &n, n_obj);
	for(; i < n; ++i)
		printf("# objconst(%d) = %.g\n", i, objconst(i));
	return s;
	}

 typedef struct
Tolinfo {
	real *tol;
	char *desc;
	} Tolinfo;

 static void
tshow(Tolinfo *t)
{
	printf("%c %-6g\t# %s\n", *t->desc, *t->tol, t->desc+1);
	}

 static char *
adjtol(char *s, Tolinfo *t)
{
	char *se;
	real r;

	s = skipblank(s);
	if (*s == '?') {
		tshow(t);
		return s+1;
		}
	r = strtod(s,&se);
	if (se > s) {
		*t->tol = r;
		return se;
		}
	printf("Bad value \"%.*s\" for tolerance key '%c'\n", slen(s), s, *t->desc);
	return "";
	}

 static char *
tline(char *s)
{
	Tolinfo *t;
	char c;
	int acted;
	static Tolinfo tolinfo[] = {
		{ &aextol,	"aabs. error tol." },
		{ &eps,		"erel. f.d. step factor" },
		{ &aeps,	"fabs. f.d. step increment" },
		{ &heps,	"hrel. f.d. step factor for Hessian from func diffs" },
		{ &haeps,	"kabs. f.d. step incr. for Hessian from func diffs" },
		{ &extol,	"rrel. error tol." },
		{ 0, 0 }
		};

	acted = 0;
 top:
	for(;;) {
		s = skipblank(s);
		switch((c = *s++)) {
		  case '?':
			for(t = tolinfo; t->tol; t++)
				tshow(t);
			acted = 1;
			continue;
		  case 0:
		  case '#':
			--s;
		  case ',':
		  case ';':
			break;
		  default:
			for(t = tolinfo; t->tol; t++)
				if (c == *t->desc) {
					s = adjtol(s, t);
					acted = 1;
					goto top;
					}
			printf("Unknown tolerance key '%c'; %s\n", c,
				"use t? to see keys and current settings.");
			return "";
		  }
		break;
		}
	if (!acted)
		for(t = tolinfo; t->tol; t++)
			tshow(t);
	return s;
	}

#ifndef NO_YLINES /*{*/

 static real ytol;

#define i15 (1<<15)
#define i16 (1<<16)
#define mult 16807
#define modul 2147483647

 static int
iran(int jran)
{
	/* Based on "A More Portable Fortran Random Number Generator"
	 * by Linus Schrage [ACM Trans. Math. Software 5#2(1979), 132-138].
	 * Derived from the variant in Klingman's netgen program.
	 */
	int ixhi, ixlo, ixahi, ixalo, irthi, iover,
		irtlo, ifulhi, leftlo;

	ixhi = jran / i16;
	ixlo = jran - ixhi * i16;
	ixalo = ixlo * mult;
	leftlo = ixalo / i16;
	ixahi = ixhi * mult;
	ifulhi = ixahi + leftlo;
	irtlo = ixalo - leftlo * i16;
	iover = ifulhi / i15;
	irthi = ifulhi - iover * i15;
	jran = irtlo - modul + irthi * i16 + iover;
	if (jran & 0x80000000L)
		jran += modul;
	return jran;
	}

 static char *
Yhelp(char *s)
{
	printf("\nY... commands:\n\n\
Yhelp			print this message\n\
Yseed			print current seed for generating p (0 if none)\n\
Yseed nnn		generate p from seed nnn\n\
Ytest			compare hvcomp and hvcomp[ds] results for all i\n\
Ytest nnn		compare hvcomp and hvcomp[ds] just for i = nnn\n\
Ytime [k [m [nnn]]	run timing loops: k iterations for hvcomp[ds], m for hvcomp\n\
			-- if nnn is present, just for i = nnn; else for all i\n\
Ytol [eps]		have \"Ytest\" complain only about errors > eps\n\
P range values		assign p (for Ytest; clears yseed)\n\
P p [range]		print current p (generate if not yet set)\n\
P			print whole p (for Ytest; generate if not set)\n");

	if (s) {
		if (!strncmp(s,"help",4)) {
			s += 4;
			if (*(uchar*)s > ' ')
				s = 0;
			}
		else if (*s == '?') {
			if (*(uchar*)++s > ' ')
				s = 0;
			}
		else
			s = 0;
		}
	return s;
	}

 static real *
Pgen(void)
{
	int seed;
	real *p, *pe;

	P->p = p = (real*)Malloc(n_var*sizeof(real));
	seed = yseed;
	pe = p + n_var;
	while(p < pe)
		*p++ = (seed = iran(seed)) / 2147483647.;
	return P->p;
	}

 static char *
Rline(char *s, int * retnow)
{
	char *se;
	int seed;

	s = skipblank(s);
	if (!*s) {
		printf("random seed %d\n", yseed);
 ret0:
		if (retnow)
			*retnow = 1;
		return 0;
		}
	seed = (int)strtol(s, &se, 10);
	if (seed < 0 && *(uchar*)se > ' ') {
		printf("Bad seed value: %.*s\n", slen(s), s);
		goto ret0;
		}
	if (!seed) {
		seed = (int)time(0);
		if (seed < 0)
			seed = -seed;
		if (seed <= 0)
			seed = 1;
		printf("random seed %d\n", seed);
		}
	yseed = seed;
	if (retnow)
		*retnow = 0;
	return se;
	}

 static char *
Yseed(char *s)
{
	char *se;
	int retnow;

	if (strncmp(s, "seed", 4))
		return Yhelp(0);
	s += 4;
	if (*(uchar*)s > ' ')
		return Yhelp(0);
	se = Rline(s, &retnow);
	if (!retnow)
		Pgen();
	return se;
	}

 static char *
Ytest(char *s)
{
	cchar *name;
	char *se;
	int co, i, ie, j, jx, n, nbad, nc, needhead;
	real emax, t, vmax;
	real *p, *q, *r, *y;
	varno_t k, nz, rk, *z;

	if (strncmp(s,"test",4))
		return Yhelp(0);
	s += 4;
	if (*(uchar*)s > ' ')
		return Yhelp(0);
	s = skipblank(s);
	nc = n_con;
	ie = nlc;
	i = -n_obj;
	if (*s) {
		co = (int)strtol(s,&se,10);
		if (*(uchar*)se > ' ')
			return Yhelp(0);
		if (co < i || co >= nc) {
			printf("\n\"Ytest nnn\" must have nnn in [%d, %d)\n", i, nc);
			if (ie < nc)
				printf("Constraints starting with %d are linear\n", ie);
			return 0;
			}
		i = co;
		ie = co + 1;
		}
	n = n_var;
	if (!(p = P->p))
		p = Pgen();
	q = tempvec(3*n + nc);
	r = q + n;
	y = r + n;
	z = (varno_t*)(y + nc);
	if (nc)
		memset(y, 0, nc*sizeof(real));
	printf("\n");
	nbad = needhead = 0;
	for(;i < ie; ++i) {
		if (i >= 0) {
			name = con_name(i);
			y[i] = 1.;
			hvcomp(q, p, -1, 0, y);
			y[i] = 0.;
			}
		else {
			j = -1 - i;
			name = obj_name(j);
			hvcomp(q, p, j, 0, 0);
			}
		hvcompd(r, p, i);
		jx = -1;
		emax = -1.;
		vmax = 0.;
		for(j = 0; j < n; ++j) {
			t = q[j] - r[j];
			if (t < 0.)
				t = -t;
			if (emax < t) {
				emax = t;
				jx = j;
				}
			if ((t = q[j]) < 0.)
				t = -t;
			if (vmax < t)
				vmax = t;
			if ((t = r[j]) < 0.)
				t = -t;
			if (vmax < t)
				vmax = t;
			}
		needhead = 1;
		if (emax > ytol) {
			vmax = emax / vmax;
			printf("%d: %s: max abs error %.3g (rel %.3g): "
				"q[%d] = %g, r[%d] = %g\n",
				i, name, emax, vmax, jx, q[jx], jx, r[jx]);
			needhead = 0;
			}
		else if (ytol == 0.) {
			printf("%d: %s: max abs error 0\n", i, name);
			needhead = 0;
			}
		nz = hvcomps(q, p, i, n, z);
		if (ytol == 0.) {
			printf("hvcomps(q,p) returns %ld\n", (long)nz);
			needhead = 0;
			}
		for(k = 0; k < nz; ++k) {
			if (q[k] != r[rk = z[k]]) {
				if (needhead) {
					printf("hvcomps(q,p) returns %ld\n", (long)nz);
					needhead = 0;
					++nbad;
					}
				printf("**BOTCH** q[%ld] = %.3g but r[%ld] = %.3g; diff = %.3g\n",
					(long)k, q[k], (long)z[k], rk, q[k] - r[rk]);
				}
			r[rk] = 0.;
			}
		for(j = 0; j < n; ++j) {
			if (r[j]) {
				if (needhead) {
					printf("hvcomps(q,p) returns %ld\n", (long)nz);
					++nbad;
					}
				printf("***BOTCH*** r[%d] = %.3g != 0\n", j, r[j]);
				break;
				}
			}
		}
	if (!nbad && needhead)
		printf("No errors > ytol = %.3g\n", ytol);
	freeall();
	return 0;
	}

 static char *
Ytime(char *s)
{
	cchar *name;
	char *se;
	int co, i, i0, ie, k, k1, m, n, nc;
	real t, t0;
	real *p, *q, *r, *y;
	varno_t *z;

	if (strncmp(s,"time",4))
		return Yhelp(0);
	s += 4;
	if (*(uchar*)s > ' ')
		return Yhelp(0);
	k = m = 1;
	s = skipblank(s);
	i = -n_obj;
	nc = n_con;
	ie = nlc;
	name = 0;
	if (*s) {
		k = (int)strtol(s,&se,10);
		if (k < 0 || *(uchar*)se > ' ')
			return Yhelp(0);
		s = skipblank(se);
		if (!*s)
			goto have_args;
		m = (int)strtol(s,&se,10);
		if (m < 0 || *(uchar*)se > ' ')
			return Yhelp(0);
		s = skipblank(se);
		if (!*s)
			goto have_args;
		co = (int)strtol(s,&se,10);
		if (*(uchar*)se > ' ')
			return Yhelp(0);
		if (co < i || co >= nc) {
			printf("\n\"Ytest nnn\" must have nnn in [%d, %d)\n", i, nc);
			if (ie < nc)
				printf("Constraints starting with %d are linear\n", ie);
			return 0;
			}
		i = co;
		ie = co + 1;
		if (i >= 0)
			name = con_name(i);
		else
			name = obj_name(-1 - i);
		}
 have_args:
	if (k+m == 0)
		return Yhelp(0);
	n = n_var;
	if (!(p = P->p))
		p = Pgen();
	q = tempvec(3*n + nc);
	r = q + n;
	y = r + n;
	z = (varno_t*)(y + nc);
	if (nc)
		memset(y, 0, nc*sizeof(real));
	i0 = i;
	if (k > 0) {
		t0 = xectim_();
		for(k1 = 0; k1 < k; ++k1) {
			for(i = i0; i < ie; ++i)
				hvcompd(r, p, i);
			}
		t = xectim_() - t0;
		if (name)
			printf("\n%d: %s: %.3g sec per hvcompd call\n", i0, name, t/k);
		else
			printf("\nrange %d .. %d: %.3g sec per cycle of hvcompd calls\n",
				i0, ie-1, t/k);
		t0 = xectim_();
		for(k1 = 0; k1 < k; ++k1) {
			for(i = i0; i < ie; ++i)
				hvcomps(r, p, i, n, z);
			}
		t = xectim_() - t0;
		if (name)
			printf("\n%d: %s: %.3g sec per hvcomps call\n", i0, name, t/k);
		else
			printf("\nrange %d .. %d: %.3g sec per cycle of hvcomps calls\n",
				i0, ie-1, t/k);
		}
	if (m > 0) {
		t0 = xectim_();
		for(k1 = 0; k1 < m; ++k1) {
			for(i = i0; i < ie; ++i) {
				if (i >= 0) {
					y[i] = 1.;
					hvcomp(q, p, -1, 0, y);
					y[i] = 0.;
					}
				else
					hvcomp(q, p, -1 - i, 0, 0);
				}
			}
		t = xectim_() - t0;
		if (name)
			printf("\n%d: %s: %.3g sec per hvcomp call\n", i0, name, t/k);
		else
			printf("\nrange %d .. %d: %.3g sec per cycle of hvcomp calls\n",
				i0, ie-1, t/k);
		}
	printf("\n");
	freeall();
	return 0;
	}

 static char *
Ytol(char *s)
{
	char *se;
	real t;

	if (strncmp(s,"tol",3))
		goto bad;
	s += 3;
	if (*(uchar*)s > ' ')
		goto bad;
	s = skipblank(s);
	if (!*s)
		printf("Ytol %g\n", ytol);
	else {
		t = strtod(s, &se);
		if (t < 0. || *(uchar*)se > ' ') {
 bad:
			return Yhelp(0);
			}
		ytol = t;
		s = se;
		}
	return s;
	}

 static char *
Yline(char *s)
{
	if (curchk())
		return 0;
	switch(*s) {
	  case 'h': return Yhelp(s);
	  case 's': return Yseed(s);
	  case 't': switch(s[1]) {
			case 'e': return Ytest(s);
			case 'i': return Ytime(s);
			case 'o': return Ytol(s);
			}
	  }
	return Yhelp(0);
	}

 static char *
Pline(char *s)
{
	real *p;

	if (curchk())
		return 0;
	if (!(p = P->p))
		p = Pgen();
	return wxyline(s, p, "P", n_var, 0, var_name_ASL, asl->i.vscale);
	}

#endif /*}*/

 static char *
Mline(char *s)
{
	ASL *asl0;

	if (curchk())
		return 0;
	asl0 = asl;
	s = getprobno(s);
	if (asl != asl0)
		show_curno();
	if (asl) {
		printf("%s:\n\t%lu bytes temporarily used during .nl reading\n"
			"\t%lu bytes retained after .nl reading\n",
			filename, (ULong)asl->i.temp_rd_bytes, (ULong)asl->i.rd_M1z_bytes);
#ifdef _ASL_EW_
		if (!asl->i.Ew0)
			printf("\t*** No Ew0 allocated.\n");
		else
			printf("\t%lu bytes without Ew0:\n\t%lu bytes per EvalWorkspace\n",
				(ULong)(asl->i.rd_M1z_bytes - asl->i.ew_bytes),
				(ULong)asl->i.ew_bytes);
#endif
		if (asl->i.rd_M1z_bytes != asl->i.tot_M1z_bytes)
			printf("\t%lu bytes currently in use\n", asl->i.tot_M1z_bytes);
		}
	return s;
	}

 static void
treport(char *s0, char *s, real t)
{
	printf("%g = %g sec / %d rep:\tT %.*s\n", t/nreps, t, nreps, (int)(s-s0), s0);
	}

 static char *
Tcline(char *s, char *s0, int wd)
{
	EW(EvalWorkspace *ew = asl->i.Ew0;)
	int i, i0, j, n;
	real *c, t, t0;

	s = get_range(s, &i0, &n, n_con);
	if (i0 >= n)
		return s;
	if (!wd)
		want_derivs = 0;
	if (i0 == 0 && n == n_con) {
		c = tempvec(n_con);
		t0 = xectim_();
		for(j = 0; j < nreps; ++j) {
			FirstX;
			conval(X0,c,0);
			}
		t = xectim_() - t0;
		freeall();
		}
	else {
		t0 = xectim_();
		for(j = 0; j < nreps; ++j) {
			FirstX;
			xknown(X0);
			for(i = i0; i < n; ++i)
				conival(i, X0, 0);
			}
		t = xectim_() - t0;
		xunknown();
		}
	if (!wd)
		want_derivs = 1;
	treport(s0, s, t);
	return s;
	}

 static char *
Tfline(char *s, char *s0, int wd)
{
	EW(EvalWorkspace *ew = asl->i.Ew0;)
	int i, i0, j, n;
	real t, t0;

	s = get_range(s, &i0, &n, n_obj);
	if (i0 >= n)
		return s;
	if (!wd)
		want_derivs = 0;
	t0 = xectim_();
	for(j = 0; j < nreps; ++j) {
		FirstX;
		for(i = i0; i < n; ++i)
			objval(i, X0, 0);
		}
	t = xectim_() - t0;
	if (!wd)
		want_derivs = 1;
	treport(s0, s, t);
	return s;
	}

 static char *
Tgline(char *s, char *s0)
{
	EW(EvalWorkspace *ew = asl->i.Ew0;)
	int i, i0, j, n;
	real *g, t, t0;

	s = get_range(s, &i0, &n, n_obj);
	if (i0 >= n)
		return s;
	g = tempvec(n_var);
	t0 = xectim_();
	for(j = 0; j < nreps; ++j) {
		FirstX;
		for(i = i0; i < n; ++i)
			objgrd(i, X0, g, 0);
		}
	t = xectim_() - t0;
	freeall();
	treport(s0, s, t);
	return s;
	}

 static char *
Thline(char *s, char *s0)
{
	int ir, k, uptri;
	real *h, t, t0;
	fint N, nsph;

	k = n_var;
	switch(*s++) {
	  case 'd':
		h = tempvec((k*(k+1))>>1);
		duthes(h, -1, P->w, pi0);
		t0 = xectim_();
		for(ir = 0; ir < nreps; ++ir)
			duthes(h, -1, P->w, pi0);
		break;

	  case 'f':
		N = k;
		h = tempvec(k*k);
		fullhes(h, N, -1, P->w, pi0);
		t0 = xectim_();
		for(ir = 0; ir < nreps; ++ir)
			fullhes(h, N, -1, P->w, pi0);
		break;

	  case 'v':
		if (!P->v)
			v_init();
		h = tempvec(k);
		hvcomp(h, P->v, -1, P->w, pi0);
		t0 = xectim_();
		for(ir = 0; ir < nreps; ++ir)
			hvcomp(h, P->v, -1, P->w, pi0);
		break;

	  default:
		--s;
		/* no break; */
	  case 's':
		uptri = 1;
		goto spcase;
	  case 'b':
		uptri = 0;
		goto spcase;
	  case 'l':
		uptri = 2;
	  spcase:
		nsph = sphsetup(-1,1,1,uptri);
		h = tempvec(nsph);
		sphes(h, -1, P->w, pi0);
		t0 = xectim_();
		for(ir = 0; ir < nreps; ++ir)
			sphes(h, -1, P->w, pi0);
		break;
	  }
	t = xectim_() - t0;
	freeall();
	treport(s0, s, t);
	return s;
	}

 static char *
Tjline(char *s, char *s0)
{
	EW(EvalWorkspace *ew = asl->i.Ew0;)
	int i, i0, ir, ms, n;
	real *J, t, t0;

	s = get_range(s, &i0, &n, n_con);
	if (i0 >= n)
		return s;
	if (i0 == 0 && n == n_con) {
		J = tempvec(nzjac);
		t0 = xectim_();
		for(ir = 0; ir < nreps; ++ir) {
			FirstX;
			jacval(X0, J, 0);
			}
		}
	else {
		J = tempvec(n_var);
		ms = asl->i.congrd_mode;
		asl->i.congrd_mode = 1;
		t0 = xectim_();
		for(ir = 0; ir < nreps; ++ir) {
			FirstX;
			xknown(X0);
			for(i = i0; i < n; ++i)
				congrd(i, X0, J, 0);
			}
		asl->i.congrd_mode = ms;
		xunknown();
		}
	t = xectim_() - t0;
	freeall();
	treport(s0, s, t);
	return s;
	}

 static char *
Tline(char *s)
{
	char *s0, *se;
	int c, n;

	c = *(s0 = s);
	if (!c) {
		printf("Current T repetitions = %d\n", nreps);
		return s;
		}
	if (c >= '1' && c <= '9') {
		n = (int)strtod(s,&se);
		if (n > 0 && *se <= ' ') {
			nreps = n;
			s = skipblank(se);
			c = *s;
			}
		else {
			if (*se > ' ')
				++se;
			printf("Bad repetition count \"%.*s\"\n", (int)(se-s), s);
 skiprest:
			while(*se)
				++se;
			return se;
			}
		}

	if (c) {
		s = skipblank(s+1);
		if (curchk())
			return 0;
		}
	switch(c) {
	  case 0:
		return s;
	  case 'b':
		s = Tcline(s,s0,0); break;
	  case 'c':
		s = Tcline(s,s0,1); break;
	  case 'e':
		s = Tfline(s,s0,0); break;
	  case 'f':
		s = Tfline(s,s0,1); break;
	  case 'g':
		s = Tgline(s,s0); break;
	  case 'h':
		s = Thline(s,s0); break;
	  case 'j':
		s = Tjline(s,s0); break;
	  default:
		printf("Bad \"T\" operation '%c'\n", *s);
		goto skiprest;
	  }
	return s;
	}

 static char *Iline(char*);

 static void
usage(void)
{
	printf("%s\n%s\n%s\n%s\n%s\n%s\n",
	"Allow multiple problems in memory at once, numbered from 0.\n\
Single letter commands, followed by zero or more operands:\n\n\
p n			make n the current problem;\n\
			show the current problem number and name\n\n\
r filename mode	[flags]	load problem; reply with\n\
			problem number and statistics\n\n\
	mode values:	1 = linear functions only (f_read)\n\
			2 = functions and gradients (fg_read)\n\
			3 = func, grad, and (Lagr.) Hessian (fgh_read)\n\
			4 = partially sep. func & grad (pfg_read)", "\
			5 = partially sep. func, grad, Hes (pfgh_read)\n"
#ifndef _ASL_EW_
			"\t\t\t6 = use qp_read (for use with \"Li\")\n"
#endif
	"\n\
	flags: reader flags in asl.h plus 0x1000000 to show read time\n\n\
u n			unload problem number n; if n was current,\n\
			the smallest available problem number\n\
			(if any) becomes current\n\n\
x range values		assign values to x components in range\n\
			range = lower[-upper]\n\
			first subscript = 0\n\
			n*x = n copies of x\n\
			n*x1,x2[,x3...] = n values x1,x2[,x3...],x1,x2,...\n\
			*x = fill rest of range with x\n\
			*x1,x2[,x3...] = fill rest as for n*x1,x2...."
#ifndef NO_YLINES
			"\n\
			R = fill with random values in [0, 1]\n\
			R L U = fill with random values in [L, U]"
#endif
			"\n\n\
x p [range]		print x\n\n\
y range values		assign y (negative of Langrange multipliers)\n\
y p [range]\n\
y r			set y to least-squares residuals\n\
w range values		set objective weights", "\
w p [range]		print objective weights\n\
v range values		set v (for Hessian * v products)\n\
v p [range]		print v; initially v = all ones\n\
\n\
f [range]		print objective value\n\
g [range]		print gradient\n\
c [range]		print constraint values\n\
j [range]		print jacobian\n\
h [kind]		print Hessian of Lagrangian:","\
			Hessian of w*objectives + y*constraints\n\
	kind values:	d = duthes (dense upper triangle)\n\
			f = fullhes (full dense Hessian)\n\
			s = sparse upper triangle (default)\n\
			l = sparse lower triangle\n\
			b = sparse, both triangles","\
			v = Hessian-vector product: H * v\n\
			z = show only nonzero values (may be appended to s,l,b;\n\
				otherwise assume l)\n\
d [p] g [range]		check objective gradients by finite-differences\n\
d [p] j [range]		check constraint gradients by f.d.\n", "\
d [p] h	[kind]		check Hessian by f.d.; kind 'v' = hvcomp\n\
d [p] H [kind]		check Hessian by f.d. of function values\n\
L [range]		print logical constraints\n\
L i			report indicator constraints\n\
M			report memory statistics"
#ifndef NO_YLINES
			"\n\
R seed			seed for random number generator (default 1)"
#endif
			"\n\
S n			n = 1,3,5 = scale constraints, 2,4,6 = scale variables\n\
T  [n] [bcefghj] ...	do n repetitions of timing e, f, g, b, c, j, or h evals\n\
				... = options for plain [cfghj]\n\
				b = c with want_derivs = 0\n\
				e = f with want_derivs = 0\n\
				g includes f, j includes c\n\
				h includes f and c (if implied by w and y)\n\
				default n = previous n (initially 10000)\n\
V [cv] [range]		report minimum constraint or variable-bound slack\n\
			(with negative ==> violation):\n\
				c ==> constraint (default)\n\
				v ==> variable\n\
X [filename]			read x and y from filename (default stub.sol)"
#ifndef NO_YLINES
			"\n\
Y...			commands to test hvcompi(); Yhelp for details"
#endif
			"\n\
l scale			call lagscale(scale)\n\
o [range]		print objective constant (0 unless linear)\n\
#			ignore remaining text on this line\n\
*			print and otherwise ignore this line\n\
q			quit (this file)\n\
t key=val key=val...	set tolerance key; key = ? ==> list settings\n\
			val = ? ==> list value for this key\n\
< filename		read commands from filename\n\
,			command separator (optional)\n");
	}

 static int
process(FILE *f, int badquit)
{
	char buf[4096], c, *s, *s0;
	int rc = 0;

	mejb = &exit_jb;
	if (setjmp(exit_jb.jb)) {
		rc++;
		if (badquit)
			return rc;
		}
	for(; fgets(s = buf,sizeof(buf),f); fflush(stdout))
	    for(;;) {
		s = skipblank(s0 = s);
		if (!(c = *s++))
			break;
		s = skipblank(s);
		switch(c) {
		  case '*': printf("%s", s0-1);
			    /* no break */
		  case '#': s = ""; continue;
		  case ',': break;
		  case ';': break;
		  case '<': if (!(s = Iline(s)))
				return 1;
			    break;
		  case 'L': s = Lline(s); break;
		  case 'S': s = Sline(s); break;
		  case 'V': s = Vline(s); break;
		  case 'X': s = Xline(s); break;
		  case 'c': s = cline(s); break;
		  case 'd': s = dline(s); break;
		  case 'f': s = fline(s); break;
		  case 'g': s = gline(s); break;
		  case 'h': s = hline(s); break;
		  case 'j': s = jline(s); break;
		  case 'l': s = lline(s); break;
		  case 'o': s = oline(s); break;
		  case 'p': s = pline(s); break;
		  case 'q': return rc;
		  case 'r': s = rline(s); break;
		  case 's': s = sline(s); break;
		  case 't': s = tline(s); break;
		  case 'u': s = uline(s); break;
		  case 'v': s = vline(s); break;
		  case 'w': s = wline(s); break;
		  case 'x': s = xline(s); break;
		  case 'y': s = yline(s); break;
#ifndef NO_YLINES
		  case 'R': s = Rline(s,0); break;
		  case 'P': s = Pline(s); break;
		  case 'Y': s = Yline(s); break;
#endif
		  case 'M': s = Mline(s); break;
		  case 'T': s = Tline(s); break;
		  default:
			printf("# Bad code '%c'\n", c);
			rc++;
		  case '?': usage(); s = ""; continue;
		  }
		if (!s) {
			if (badquit)
				return rc+1;
			break;
			}
		}
	return rc;
	}

 static char *
Iline(char *s)
{
	FILE *f;
	char buf[256];
	int c, i, q;
	unsigned char *su;

	if (!(q = *s))
		return s;
	su = (unsigned char *)(s+1);
	i = 0;
	if (q == '"' || q == '\'') {
		while((c = su[i]) != q) {
			if (c < ' ') {
				printf("Bad character 0x%02x in filename %.*s...\"\n",
					c, i+2, s);
				return 0;
				}
			if (i >= sizeof(buf)-1)
				goto bufoverflow;
			buf[i++] = c;
			}
		if (!i) {
			printf("Emtpy filename %c%c\n", q,q);
			return 0;
			}
		if (s[i += 2] > ' ') {
			printf("Unexpected character '%c' after closing %c\n", s[i], q);
			return 0;
			}
		}
	else {
		while((c = s[i]) > ' ') {
			if (i >= sizeof(buf) - 1) {
 bufoverflow:
				buf[i] = 0;
				printf("Oversize filename in \"< %s...\"\n", buf);
				return 0;
				}
			buf[i++] = c;
			}
		}
	buf[i] = 0;
	if (!(f = fopen(buf, "r"))) {
		printf("Cannot open \"%s\"\n", buf);
			return 0;
		}
	c = process(f,1);
	fclose(f);
	if (c)
		return 0;
	return s + i;
	}

 static int
ix_usage(void)
{
	const char **o = ix_details_ASL, *s;

	printf("-i options:\n");
	while((s = *o++))
		printf("\t%s\n", s);
	return 0;
	}

 int
main(int argc, char **argv)
{
	FILE *f;
	char *s;
	int rc, u;

	Stderr_init_ASL();
	progname = argv[0];
	Pavail = Pfree = (Pinfo*)Malloc(npinfo*sizeof(Pinfo));
	Plast = Pavail + npinfo;
	memset(Pavail, 0, npinfo*sizeof(Pinfo));

	rc = 1;
	u = 0;
	while ((s = *++argv) && *s == '-') {
 more:
		switch(s[1]) {
		 case 'f':
			funcshow = 0;
			if (s[2]) {
				++s;
				goto more;
				}
			continue;
		 case 'i':
			if (s[2]) {
				if (s[2] == '?' && !s[3])
					return ix_usage();
				i_option_ASL = s + 2;
				continue;
				}
			if (!(s = *++argv))
				goto Usage;
			if (*s == '?' && !s[1])
				return ix_usage();
			i_option_ASL = s;
			continue;
		 case 'u':
			u = 1;
			if (!s[2])
				rc = 0;
			goto Usage;
		 case 't':
			Rtime = 0;
			if (s[2]) {
				++s;
				goto more;
				}
			continue;
		 case '-':
			if (!s[2])
				break;
			if (!strcmp(s+2,"help"))
				rc = 0;
			goto Usage;
		 case '?':
		 case 'h':
			if (!s[2])
				rc = 0;
		 default:
 Usage:
		 f = rc ? Stderr : stdout;
		 fprintf(f,"Usage: %s [option] [file [file...]]\n"
			"to test nonlinear amplsolver function and derivative routines.\n"
			"Options:\n\t-f  ==> suppress listing available imported functions\n"
			"\t-ix ==> import user-defined functions from x; -i? gives details\n"
			"\t-t  ==> suppress report of read time\n"
			"\t-u  ==> show detailed usage and exit\n", progname);
		 if (u)
			usage();
		 else if (!rc)
			printf("Invoke and type \"?\" for command summary.\n");
		 return rc;
		 }
		}
	if (!s)
		return process(stdin,0);
	rc = 0;
	do {
		f = fopen(s,"r");
		if (!f) {
			fprintf(Stderr, "Can't open %s\n", s);
			return 1;
			}
		rc |= process(f,1);
		fclose(f);
		} while((s = *++argv));
	return rc;
	}
