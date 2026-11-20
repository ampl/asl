/* Derivative consistency check for the ASL readers (solvers and solvers2).
 *
 * usage: derivcheck stub reader seed order [dumpfile [reffile]]
 *   reader:   fg | fgh | pfg | pfgh
 *   seed:     0 = start point of the .nl file, k > 0 = k-th perturbed point
 *   order:    v = values before derivatives, g = derivatives first,
 *             r = read only, allowing logical constraints (reader tests),
 *             e = expected evaluation errors: evaluate exactly at the
 *                 start point (seed must be 0), values first, then again
 *                 derivatives first after an evaluation at a nearby point
 *                 (so that nothing cached for the start point is reused);
 *                 the routines listed in "stub.err" must report an error
 *                 (nonzero *nerror), all others must succeed
 *   dumpfile: receives all computed values, one "key i j value" per line
 *   reffile:  dump of another reader/library at the same point; function
 *             values, gradients, Jacobian entries and the evaluation
 *             errors must agree with it, record for record
 *
 * Gradients and Jacobians are compared with central differences of the
 * function values, Hessians (hvcomp, sphes, fullhes) with differences of
 * the Lagrangian gradient and with each other.  An evaluation error skips
 * the comparisons of the routine that reported it; any error skips the
 * Hessian checks.
 * "stub.err": one "routine index" line per expected error, routine one of
 * objval, objgrd (index = objective), conval, jacval (index 0); "#" starts
 * a comment.
 * Exit status: 0 = consistent, 3 = mismatch, 2 = usage error.
 * Compile with -DDERIVCHECK_ASL2 for solvers2.
 */
#include "asl.h"
#include <string.h>

static ASL *asl;
static FILE *dump;
static int nbad, nerr;
static double worst;

/* reference values from another reader */
typedef struct { char key[16]; int i, j; double v; int seen; } RefVal;
static RefVal *ref;
static int nref;

/* evaluation errors: expected ("stub.err", order e) and reported */
typedef struct { char key[16]; int i; } ErrRec;
#define MAXERR 64
static ErrRec experr[MAXERR], goterr[MAXERR];
static int nexperr, ngoterr, checkerr;

static unsigned long long rs = 12345;
static double rnd(void) {	/* uniform in (0,1) */
	rs = rs * 6364136223846793005ULL + 1442695040888963407ULL;
	return ((rs >> 11) + 0.5) / 9007199254740992.0;
}

static void cmp(const char *key, int i, int j, double a, double fd);

static void load_ref(const char *name) {
	FILE *f = fopen(name, "r");
	RefVal r;
	int cap = 0;
	if (!f) {
		fprintf(stderr, "cannot open reference %s\n", name);
		exit(2);
	}
	r.seen = 0;
	while (fscanf(f, "%15s %d %d %lf", r.key, &r.i, &r.j, &r.v) == 4) {
		if (nref >= cap)
			ref = (RefVal*)realloc(ref, (cap = 2*cap + 64) * sizeof(RefVal));
		ref[nref++] = r;
	}
	fclose(f);
}

static void mismatch(const char *what, const char *key, int i, int j) {
	if (++nbad <= 5)
		fprintf(stderr, "  MISMATCH %s: %s[%d,%d]\n", what, key, i, j);
	worst = 1.;
}

/* Every record (value, derivative or "*_err" evaluation error) must have a
 * counterpart in the other reader's dump, and vice versa (see ref_unseen).
 * Hessian entries are not compared across readers: fg has none. */
static void cmp_ref(const char *key, int i, int j, double v) {
	int k, k1 = -1;
	if (key[0] == 'H')
		return;
	for (k = 0; k < nref; ++k)	/* a repeated record: next unseen copy */
		if (ref[k].i == i && ref[k].j == j && !strcmp(ref[k].key, key)) {
			k1 = k;
			if (!ref[k].seen)
				break;
		}
	if (k1 < 0) {
		mismatch("not from the other reader", key, i, j);
		return;
	}
	ref[k1].seen = 1;
	cmp("other-reader", i, j, v, ref[k1].v);
}

/* records of the other reader that this one did not produce */
static void ref_unseen(void) {
	int k;
	for (k = 0; k < nref; ++k)
		if (ref[k].key[0] != 'H' && !ref[k].seen)
			mismatch("only from the other reader", ref[k].key, ref[k].i,
				ref[k].j);
}

static void load_experr(const char *stub) {
	char name[1024], line[256], routine[16];
	FILE *f;
	int i;
	snprintf(name, sizeof(name), "%s.err", stub);
	if (!(f = fopen(name, "r"))) {
		fprintf(stderr, "cannot open %s\n", name);
		exit(2);
	}
	while (fgets(line, sizeof(line), f))
		if (sscanf(line, "%11s %d", routine, &i) == 2 && routine[0] != '#'
		 && nexperr < MAXERR) {
			snprintf(experr[nexperr].key, 16, "%s_err", routine);
			experr[nexperr++].i = i;
		}
	fclose(f);
	checkerr = 1;
}

static int has_err(const ErrRec *e, int n, const char *key, int i) {
	int k;
	for (k = 0; k < n; ++k)
		if (e[k].i == i && !strcmp(e[k].key, key))
			return 1;
	return 0;
}

/* compare the errors of one evaluation pass with the expected ones */
static void check_errs(const char *pass) {
	int k;
	for (k = 0; k < nexperr; ++k)
		if (!has_err(goterr, ngoterr, experr[k].key, experr[k].i)) {
			fprintf(stderr, "  (%s)", pass);
			mismatch("expected error not reported", experr[k].key,
				experr[k].i, 0);
		}
	for (k = 0; k < ngoterr; ++k)
		if (!has_err(experr, nexperr, goterr[k].key, goterr[k].i)) {
			fprintf(stderr, "  (%s)", pass);
			mismatch("unexpected error", goterr[k].key, goterr[k].i, 0);
		}
	ngoterr = 0;
}

static void put(const char *key, int i, int j, double v) {
	if (dump)
		fprintf(dump, "%s %d %d %.10g\n", key, i, j, v);
	if (ref)
		cmp_ref(key, i, j, v);
	if (strstr(key, "_err") && ngoterr < MAXERR) {
		strcpy(goterr[ngoterr].key, key);
		goterr[ngoterr++].i = i;
	}
}

static void cmp(const char *key, int i, int j, double a, double fd) {
	double e = fabs(a - fd) / (1. + fabs(a) + fabs(fd));
	if (!(e <= 1e-4)) {	/* also catches NaN */
		if (++nbad <= 5)
			fprintf(stderr, "  MISMATCH %s[%d,%d]: analytic %.8g  fd %.8g\n",
				key, i, j, a, fd);
		if (e > worst || e != e)
			worst = e;
	}
}

static int reader_code(const char *r) {
	if (!strcmp(r, "fg"))   return ASL_read_fg;
	if (!strcmp(r, "fgh"))  return ASL_read_fgh;
	if (!strcmp(r, "pfg"))  return ASL_read_pfg;
	if (!strcmp(r, "pfgh")) return ASL_read_pfgh;
	fprintf(stderr, "unknown reader %s\n", r);
	exit(2);
	return 0;
}

static void make_point(real *x, int seed) {
	int j, k;
	real lb, ub, w;
	for (k = 0; k < 17 * seed; ++k)
		rnd();
	for (j = 0; j < n_var; ++j) {
		lb = LUv[2*j]; ub = LUv[2*j+1];
		x[j] = X0[j];
		if (seed)
			x[j] += (rnd() - .3) * .4 * (1. + fabs(x[j]));
		if (x[j] <= lb || x[j] >= ub) {	/* move strictly inside */
			w = seed ? rnd() : .5;
			if (lb > negInfinity && ub < Infinity)
				x[j] = lb + (.1 + .8*w) * (ub - lb);
			else if (lb > negInfinity)
				x[j] = lb + .1 + w;
			else if (ub < Infinity)
				x[j] = ub - .1 - w;
		}
	}
}

/* gradient of the Lagrangian: objgrd(nobj) + J' y; returns 0 on error */
static int lag_grad(int nobj, real *x, real *y, real *g, real *J) {
	fint ne = 0;
	cgrad *cg;
	int i;
	memset(g, 0, n_var * sizeof(real));
	if (nobj >= 0) {
		objgrd(nobj, x, g, &ne);
		if (ne) return 0;
	}
	if (n_con) {
		jacval(x, J, &ne);
		if (ne) return 0;
		for (i = 0; i < n_con; ++i)
			for (cg = Cgrad[i]; cg; cg = cg->next)
				g[cg->varno] += y[i] * J[cg->goff];
	}
	return 1;
}

/* values at x; fok[i], *cok: objective i and the constraints evaluated */
static void values(real *x, real *f, real *c, int *fok, int *cok) {
	fint ne;
	int i;
	for (i = 0; i < n_obj; ++i) {
		ne = 0;
		f[i] = objval(i, x, &ne);
		fok[i] = !ne;
		if (ne) { ++nerr; put("objval_err", i, 0, 1); }
		else put("f", i, 0, f[i]);
	}
	*cok = 0;
	if (n_con) {
		ne = 0;
		conval(x, c, &ne);
		*cok = !ne;
		if (ne) { ++nerr; put("conval_err", 0, 0, 1); }
		else for (i = 0; i < n_con; ++i) put("c", i, 0, c[i]);
	}
}

/* objective values and constraint bodies at a difference point */
static void values_at(real *x, real *f, real *c, int *fok, int *cok) {
	fint ne;
	int i;
	*fok = 1;
	for (i = 0; i < n_obj; ++i) {
		ne = 0;
		f[i] = objval(i, x, &ne);
		if (ne) *fok = 0;
	}
	*cok = 1;
	if (n_con) {
		ne = 0;
		conval(x, c, &ne);
		if (ne) *cok = 0;
	}
}

/* Evaluate at a point near x, so that the next evaluation at x starts
 * afresh instead of reusing the values ASL keeps for x. */
static void leave_point(real *x, real *xp, real *f, real *c) {
	int j, fok, cok;
	for (j = 0; j < n_var; ++j)
		xp[j] = x[j] + 1e-3 * (1. + fabs(x[j]));
	values_at(xp, f, c, &fok, &cok);
}

static void derivs(real *x, real *G, real *J, int *gok, int *jok) {
	fint ne;
	cgrad *cg;
	int i, j;
	for (i = 0; i < n_obj; ++i) {
		ne = 0;
		objgrd(i, x, G + i*n_var, &ne);
		gok[i] = !ne;
		if (ne) { ++nerr; put("objgrd_err", i, 0, 1); }
		else for (j = 0; j < n_var; ++j) put("g", i, j, G[i*n_var + j]);
	}
	*jok = 0;
	if (n_con) {
		ne = 0;
		jacval(x, J, &ne);
		*jok = !ne;
		if (ne) { ++nerr; put("jacval_err", 0, 0, 1); }
		else for (i = 0; i < n_con; ++i)
			for (cg = Cgrad[i]; cg; cg = cg->next)
				put("J", i, cg->varno, J[cg->goff]);
	}
}

int main(int argc, char **argv) {
	FILE *nl;
	fint ne;
	cgrad *cg;
	int hes, i, j, k, nobj, rd, seed, jok, *gok, *fok, cok, haveh, nz;
	real *x, *xp, *f, *c, *fp, *cp, *fm, *cm, *G, *J, *y, *H, *Hs, *hv, *p;
	real *gp, *gm, *J2, h, t;
	char order;

	if (argc < 5) {
		fprintf(stderr, "usage: %s stub reader seed order [dump [ref]]\n", argv[0]);
		return 2;
	}
	rd = reader_code(argv[2]);
	seed = atoi(argv[3]);
	order = argv[4][0];
	dump = argc > 5 ? fopen(argv[5], "w") : 0;
	if (argc > 6)
		load_ref(argv[6]);
	asl = ASL_alloc(rd);
	nl = jac0dim(argv[1], (fint)strlen(argv[1]));
	X0 = (real*)Malloc(n_var * sizeof(real));
	k = order == 'r' ? ASL_allow_CLP : 0;
	switch (rd) {
	  case ASL_read_fg:   fg_read(nl, k);   break;
#ifndef DERIVCHECK_ASL2	/* solvers2 declares but does not define these two */
	  case ASL_read_fgh:  fgh_read(nl, k);  break;
	  case ASL_read_pfg:  pfg_read(nl, k);  break;
#endif
	  case ASL_read_pfgh: pfgh_read(nl, k); break;
	}
	if (order == 'r') {
		printf("read n=%d m=%d\n", n_var, n_con);
		return 0;
	}
	hes = rd == ASL_read_fgh || rd == ASL_read_pfgh;
	k = n_var * sizeof(real);
	x = Malloc(2*k); xp = x + n_var;
	f = Malloc(3 * (n_obj + 1) * sizeof(real)); fp = f + n_obj + 1; fm = fp + n_obj + 1;
	c = Malloc(3 * (n_con + 1) * sizeof(real)); cp = c + n_con + 1; cm = cp + n_con + 1;
	G = Malloc((n_obj + 1) * k);
	J = Malloc(2 * (nzc + 1) * sizeof(real)); J2 = J + nzc + 1;
	y = Malloc((n_con + 1) * sizeof(real));
	gok = Malloc(2 * (n_obj + 1) * sizeof(int)); fok = gok + n_obj + 1;
	H = Malloc((size_t)n_var * k); Hs = Malloc((size_t)n_var * k);
	hv = Malloc(4*k); p = hv + n_var; gp = p + n_var; gm = gp + n_var;
	make_point(x, seed);
	if (order == 'e') {	/* exactly the start point, maybe on a bound */
		if (seed) {
			fprintf(stderr, "order e needs seed 0\n");
			return 2;
		}
		memcpy(x, X0, k);
		load_experr(argv[1]);
	}
	for (i = 0; i < n_con; ++i)
		y[i] = .5 + rnd();
	nobj = n_obj > 0 ? 0 : -1;

	if (order == 'g') {
		derivs(x, G, J, gok, &jok);
		values(x, f, c, fok, &cok);
	} else {
		values(x, f, c, fok, &cok);
		derivs(x, G, J, gok, &jok);
	}
	if (checkerr) {	/* order e: both call orders must give the errors */
		check_errs("values first");
		leave_point(x, xp, fp, cp);
		derivs(x, G, J, gok, &jok);
		values(x, f, c, fok, &cok);
		check_errs("derivatives first");
	}

	/* Hessian of the Lagrangian, column by column, then sparse and dense */
	haveh = 0;
	if (hes && n_var <= 400 && !nerr) {
		haveh = 1;
		for (j = 0; j < n_var; ++j) {
			memset(p, 0, k);
			p[j] = 1.;
			hvcomp(hv, p, nobj, 0, y);
			for (i = 0; i < n_var; ++i) {
				H[i + (size_t)j*n_var] = hv[i];
				if (hv[i]) put("H", i, j, hv[i]);
			}
		}
		for (j = 0; j < n_var; ++j)	/* symmetry */
			for (i = 0; i < j; ++i)
				cmp("Hsym", i, j, H[i + (size_t)j*n_var], H[j + (size_t)i*n_var]);
#ifndef DERIVCHECK_ASL2
		if (rd == ASL_read_pfgh)
#endif
		{
			real *Hsp;
			fint *hcs, *hrn;
			nz = (int)sphsetup(nobj, 0, 1, 1);
			Hsp = Malloc((nz + 1) * sizeof(real));
			sphes(Hsp, nobj, 0, y);
			memset(Hs, 0, (size_t)n_var * k);
			hcs = sputinfo->hcolstarts; hrn = sputinfo->hrownos;
			for (j = 0; j < n_var; ++j)
				for (i = hcs[j]; i < hcs[j+1]; ++i)
					Hs[hrn[i] + (size_t)j*n_var] = Hsp[i];
			for (j = 0; j < n_var; ++j)
				for (i = 0; i <= j; ++i)
					cmp("sphes-vs-hvcomp", i, j, Hs[i + (size_t)j*n_var],
						H[i + (size_t)j*n_var]);
			fullhes(Hs, n_var, nobj, 0, y);
			for (j = 0; j < n_var; ++j)
				for (i = 0; i < n_var; ++i)
					cmp("fullhes-vs-hvcomp", i, j, Hs[i + (size_t)j*n_var],
						H[i + (size_t)j*n_var]);
			free(Hsp);
		}
	}

	/* finite differences, for each routine that evaluated without error */
	if (n_var <= 400)
	for (j = 0; j < n_var; ++j) {
		int okfp, okfm, okcp, okcm, okgp = 0, okgm = 0;
		h = 1e-5 * (1. + fabs(x[j]));
		memcpy(xp, x, k);
		xp[j] = x[j] + h;
		values_at(xp, fp, cp, &okfp, &okcp);
		if (haveh && okfp && okcp) okgp = lag_grad(nobj, xp, y, gp, J2);
		xp[j] = x[j] - h;
		values_at(xp, fm, cm, &okfm, &okcm);
		if (haveh && okfm && okcm) okgm = lag_grad(nobj, xp, y, gm, J2);
		if (okfp && okfm)
			for (i = 0; i < n_obj; ++i)
				if (gok[i])
					cmp("objgrd", i, j, G[i*n_var + j], (fp[i] - fm[i]) / (2*h));
		if (jok && okcp && okcm) {
			for (i = 0; i < n_con; ++i) {
				t = 0.;
				for (cg = Cgrad[i]; cg; cg = cg->next)
					if (cg->varno == j) t = J[cg->goff];
				cmp("jacval", i, j, t, (cp[i] - cm[i]) / (2*h));
			}
		}
		if (okgp && okgm)
			for (i = 0; i < n_var; ++i)
				cmp("hvcomp", i, j, H[i + (size_t)j*n_var], (gp[i] - gm[i]) / (2*h));
	}

	/* back at x: values and Hessian must be reproduced (stale-state check) */
	if (n_var <= 400) {
		for (i = 0; i < n_obj; ++i) {
			ne = 0;
			t = objval(i, x, &ne);
			if (!ne && fok[i]) cmp("objval-again", i, 0, t, f[i]);
			ne = 0;
			objgrd(i, x, gp, &ne);
			if (!ne && gok[i])
				for (j = 0; j < n_var; ++j) cmp("objgrd-again", i, j, gp[j], G[i*n_var + j]);
		}
		if (n_con) {
			ne = 0;
			conval(x, cp, &ne);
			if (!ne && cok) for (i = 0; i < n_con; ++i) cmp("conval-again", i, 0, cp[i], c[i]);
			ne = 0;
			jacval(x, J2, &ne);
			if (!ne && jok) for (i = 0; i < nzc; ++i) cmp("jacval-again", i, 0, J2[i], J[i]);
		}
		if (haveh) {
			for (j = 0; j < n_var; ++j) p[j] = .5 + .1*(j % 7);
			hvcomp(hv, p, nobj, 0, y);
			for (i = 0; i < n_var; ++i) {
				t = 0.;
				for (j = 0; j < n_var; ++j) t += H[i + (size_t)j*n_var] * p[j];
				cmp("hvcomp-again", i, 0, hv[i], t);
			}
		}
	}
	if (ref)
		ref_unseen();
	if (dump)
		fclose(dump);
	printf("n=%d m=%d nlc=%d nlo=%d evalerrs=%d mismatches=%d worst=%.3g\n",
		n_var, n_con, nlc, nlo, nerr, nbad, worst);
	return nbad ? 3 : 0;
}
