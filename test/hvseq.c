/* Hessian-vector products over a sequence of points, one product per point.
 *
 * usage: hvseq stub [npoints [mode [zeromask]]]
 * Pass 1 goes through npoints (default 12) random points near the start
 * point; at each it evaluates functions and gradients and computes one
 * product of the Hessian of the Lagrangian.  Pass 2 compares every product
 * with central differences of the gradient of the Lagrangian.  A product
 * computed with the Hessian of an earlier point shows up as a mismatch.
 * mode (default c) selects the routine:
 *   c  hvcomp;
 *   z  hvcomp with ihd_limit = 0 (no Hessian kept between products);
 *   d  hvcompd, once for the objective and once per constraint, summed
 *      with the constraint weights;
 *   s  hvcomps, likewise; its sparse result is scattered by the indices;
 *   D, S  as d and s, but the products are the first calls at each point:
 *      hvcompd and hvcomps evaluate the functions themselves;
 *   m  hvcompd for the objective as the first call at each point, then the
 *      constraints are evaluated and hvcomp gives the product;
 *   o  the objective is evaluated and hvcomp gives its product alone, then
 *      hvcompd adds the constraints as the first calls for them;
 *   h  the objective is evaluated and hvcomp gives its product alone, then
 *      the constraints are evaluated and hvcomp gives the product of the
 *      Lagrangian at the same point (funnels must be recomputed);
 *   F, G  as D, but at the first point the products are the first
 *      evaluation calls of the run, before the pass that finds the points
 *      to skip: F calls hvcompd for the constraints first, then the
 *      objective, G the objective first.  Shared defined variables of the
 *      functions not yet evaluated have no values then;
 *   q  the objective is evaluated and hvcomp gives its product alone, then
 *      hvcomp gives the product of the Lagrangian without evaluating the
 *      constraints first (hvcomp evaluates them);
 *   Q  the mirror: the constraints are evaluated and hvcomp gives their
 *      product alone, then that of the Lagrangian;
 *   r, R  as q, Q with ihd_limit = 0;
 *   y  hvcomp for the Lagrangian with nothing evaluated before (in solvers2
 *      through hvcompe with a null nerror, which must act as hvcomp);
 *   Y  as y with ihd_limit = 0.
 * In modes D, S, m, o, h, F and G, a first pass evaluates the functions at
 * every point to find the points to skip; consecutive points differ, so
 * each product starts at a new x (with npoints = 1 it does not).  Modes q,
 * Q, r, R, y, Y have no such pass: every call at a point is a first call,
 * and an evaluation error ends the run.
 * zeromask (default 0): bit i set gives constraint i weight 0 (i < 30).
 * Works with solvers and solvers2 (pfgh_read).
 * Exit status: 3 on mismatch, 4 when fewer than two points (one if
 * npoints is 1) could be compared, so a run that skips every point fails.
 */
#include "asl_pfgh.h"

static ASL *asl;
static real *Y;
static int mode = 'c';
static varno_t *Z;	/* indices of the nonzeros from hvcomps */
static unsigned long long seed = 12345;

/* uniform in [0,1] */
static real rnd(void) {
	seed = seed*6364136223846793005ULL + 1442695040888963407ULL;
	return (real)((seed >> 33) % 2001)/2000.;
	}

/* NaN fails v - v == 0 as well */
static int finite_value(real v) { return v - v == 0.; }

/* A point near the start point, strictly inside the variable bounds:
 * functions that are only defined inside the bounds stay defined. */
static void make_point(real *x) {
	int j; real lb, ub, w;
	for (j = 0; j < n_var; j++) {
		lb = LUv[2*j]; ub = LUv[2*j+1];
		x[j] = X0[j] + .3*(2.*rnd() - 1.)*(1. + fabs(X0[j]));
		if (x[j] > lb && x[j] < ub) continue;
		w = rnd();
		if (lb > negInfinity && ub < Infinity) x[j] = lb + (.1 + .8*w)*(ub - lb);
		else if (lb > negInfinity) x[j] = lb + .1 + w;
		else if (ub < Infinity) x[j] = ub - .1 - w;
		}
	}

/* Direction scaled with the point: x +/- h*p then changes every component
 * by a representable amount, whatever the magnitude of x. */
static void make_direction(const real *x, real *p) {
	int j;
	for (j = 0; j < n_var; j++)	/* a fixed variable does not move */
		p[j] = LUv[2*j] < LUv[2*j+1] ? (1. - .125*(j % 7))*(1. + fabs(x[j])) : 0.;
	}

/* Finite-difference step: 1e-6, reduced so that x +/- h*p stays inside
 * the variable bounds (at most half the distance to the nearest bound). */
static real fd_step(const real *x, const real *p) {
	int j; real h = 1e-6, room;
	for (j = 0; j < n_var; j++) {
		if (!p[j]) continue;
		room = x[j] - LUv[2*j];
		if (room > LUv[2*j+1] - x[j]) room = LUv[2*j+1] - x[j];
		if (h*fabs(p[j]) > .5*room) h = .5*room/fabs(p[j]);
		}
	return h;
	}

static int lagrgrad(real *x, real *gl, real *tmp) {
	fint ne = 0; int i, j;
	for (j = 0; j < n_var; j++) gl[j] = 0.;
	if (n_obj) { objgrd(0, x, gl, &ne); if (ne) return 1; }
	for (i = 0; i < n_con; i++) {
		congrd(i, x, tmp, &ne); if (ne) return 1;
		for (j = 0; j < n_var; j++) gl[j] += Y[i]*tmp[j];
		}
	return 0;
	}

/* Hessian of one function (co >= 0: constraint co, -1: objective 0)
 * times p, added with weight t to hv; w: n_var reals. */
static void add_function_product(int co, real t, real *p, real *hv, real *w) {
	varno_t k, nz;
	int j;
	if (mode == 'd' || mode == 'D' || mode == 'F' || mode == 'G') {
		hvcompd(w, p, co);
		for (j = 0; j < n_var; j++) hv[j] += t*w[j];
		return;
		}
	nz = hvcomps(w, p, co, n_var, Z);	/* returns the number of nonzeros */
	for (k = 0; k < nz && k < n_var; k++) hv[Z[k]] += t*w[k];
	}

/* The product of mode 'd' or 's': objective plus weighted constraints. */
static void product_by_function(real *p, real *hv, real *w) {
	int i, j;
	for (j = 0; j < n_var; j++) hv[j] = 0.;
	if (n_obj) add_function_product(-1, 1., p, hv, w);
	for (i = 0; i < n_con; i++) add_function_product(i, Y[i], p, hv, w);
	/* hvcompd and hvcomps leave x "known" (both libraries, unlike hvcomp):
	 * later evaluations at other points would return the values at this
	 * point.  Calling programs must call xunknown() afterwards. */
	xunknown();
	}

/* Values and gradients at x; returns 1 on an evaluation error.  With
 * no_obj, only the constraints.  w: n_var reals, c: n_con reals. */
static int evaluate_at(real *x, real *w, real *c, int no_obj) {
	fint ne = 0; int i;
	if (n_obj && !no_obj) objval(0, x, &ne);
	if (n_obj && !no_obj && !ne) objgrd(0, x, w, &ne);
	if (n_con && !ne) conval(x, c, &ne);
	for (i = 0; i < n_con && !ne; i++) congrd(i, x, w, &ne);
	return ne != 0;
	}

/* Modes D, S, m, o, h, q, Q, r, R, y, Y: the product is the first call at
 * x; in modes D, S, m, o and h its functions
 * were evaluated without error at an earlier pass. */
static void first_product_at(real *x, real *p, real *hv, real *w, real *c) {
	fint ne = 0; int i, j;
	xknown(x);	/* hvcompd works at the most recent x */
	if (mode == 'q' || mode == 'r') {
		if (n_obj) { objval(0, x, &ne); objgrd(0, x, w, &ne); }
		hvcomp(hv, p, n_obj ? 0 : -1, 0, 0);
		hvcomp(hv, p, n_obj ? 0 : -1, 0, n_con ? Y : 0);
		xunknown();
		return;
		}
	if (mode == 'Q' || mode == 'R') {
		evaluate_at(x, w, c, 1);
		hvcomp(hv, p, -1, 0, n_con ? Y : 0);
		hvcomp(hv, p, n_obj ? 0 : -1, 0, n_con ? Y : 0);
		xunknown();
		return;
		}
	if (mode == 'y' || mode == 'Y') {
#ifdef _ASL_EW_	/* hvcompe with nerror == 0: as hvcomp */
		asl->p.Hvcompe(asl->i.Ew0, hv, p, n_obj ? 0 : -1, 0, n_con ? Y : 0, (fint*)0);
#else
		hvcomp(hv, p, n_obj ? 0 : -1, 0, n_con ? Y : 0);
#endif
		xunknown();
		return;
		}
	if (mode == 'h') {
		if (n_obj) { objval(0, x, &ne); objgrd(0, x, w, &ne); }
		hvcomp(hv, p, n_obj ? 0 : -1, 0, 0);
		xunknown();
		evaluate_at(x, w, c, 1);
		hvcomp(hv, p, n_obj ? 0 : -1, 0, n_con ? Y : 0);
		return;
		}
	if (mode == 'o') {
		if (n_obj) { objval(0, x, &ne); objgrd(0, x, w, &ne); }
		hvcomp(hv, p, n_obj ? 0 : -1, 0, 0);
		for (i = 0; i < n_con; i++) {
			hvcompd(w, p, i);	/* evaluates constraint i at x */
			for (j = 0; j < n_var; j++) hv[j] += Y[i]*w[j];
			}
		xunknown();
		return;
		}
	if (mode != 'm') {
		product_by_function(p, hv, w);
		return;
		}
	if (n_obj) hvcompd(w, p, -1);	/* evaluates the objective at x */
	xunknown();
	evaluate_at(x, w, c, 1);	/* hvcomp needs the constraint values */
	hvcomp(hv, p, n_obj ? 0 : -1, 0, n_con ? Y : 0);
	}

/* Modes F and G at the first point: hvcompd for each function as the
 * first evaluation calls of the run; returns 1 on an evaluation error. */
static int product_first_calls(real *x, real *p, real *hv, real *w) {
#ifdef _ASL_EW_	/* solvers2: hvcompde reports errors */
	fint ne = 0; int i, j, k, n = n_con + (n_obj > 0);
	for (j = 0; j < n_var; j++) hv[j] = 0.;
	xknown(x);
	for (k = 0; k < n && !ne; k++) {
		if (mode == 'F')	/* constraints first */
			i = k < n_con ? k : -1;
		else			/* the objective first */
			i = n_obj ? k - 1 : k;
		asl->p.Hvcompde(asl->i.Ew0, w, p, i, &ne);
		for (j = 0; j < n_var && !ne; j++) hv[j] += (i < 0 ? 1. : Y[i])*w[j];
		}
	xunknown();
	return ne != 0;
#else
	return 1;
#endif
	}

/* Values and gradients at x, then one product; returns 1 on an evaluation
 * error.  w: n_var reals, c: n_con reals. */
static int product_at(real *x, real *p, real *hv, real *w, real *c) {
	if (evaluate_at(x, w, c, 0)) return 1;
	if (mode == 'd' || mode == 's')
		product_by_function(p, hv, w);
	else	/* no weight vector: it would need one weight per objective */
		hvcomp(hv, p, n_obj ? 0 : -1, 0, n_con ? Y : 0);
	return 0;
	}

/* Compares hv with central differences along p; w: 4*n_var reals.
 * Returns 0 if equal, 1 on a mismatch (reported), 2 if not evaluable. */
static int compare_at(int k, real *x, real *p, real *hv, real *w, int report) {
	int j, n = n_var; real h = fd_step(x, p), d, sc;
	for (j = 0; j < n; j++) w[j] = x[j] + h*p[j];
	if (lagrgrad(w, w + n, w + 3*n)) return 2;
	for (j = 0; j < n; j++) w[j] = x[j] - h*p[j];
	if (lagrgrad(w, w + 2*n, w + 3*n)) return 2;
	for (j = 0; j < n; j++) {
		d = (w[n+j] - w[2*n+j])/(2*h); sc = 1. + fabs(d) + fabs(hv[j]);
		/* an infinite value would pass the tolerance test (inf <= inf) */
		if (finite_value(d) && finite_value(hv[j]) && fabs(d - hv[j]) <= 1e-4*sc)
			continue;
		if (report)
			printf("  MISMATCH point %d, hv[%d]: mode %c %.10g  fd %.10g\n", k, j, mode, hv[j], d);
		return 1;
		}
	return 0;
	}

int main(int argc, char **argv) {
	FILE *nl; int K = 12, n, m, i, k, r;
	int bad = 0, skipped = 0, checked = 0;
	real *X, *P, *HV, *w, *c; char *err;
	unsigned long zmask = 0;
	if (argc < 2) { fprintf(stderr, "usage: hvseq stub [npoints [mode [zeromask]]]\n"); return 2; }
	if (argc > 2) {	/* npoints: a positive integer */
		long lk = strtol(argv[2], &err, 10);
		if (*argv[2] == 0 || *err || lk < 1 || lk > 1000000) {
			fprintf(stderr, "usage: hvseq stub [npoints [mode]]; npoints"
				" must be a positive integer, not '%s'\n", argv[2]);
			return 2;
			}
		K = (int)lk;
		}
	if (argc > 3) mode = argv[3][0];
	if (!strchr("czdsDSmohFGqQrRyY", mode) || !mode) { fprintf(stderr, "unknown mode %s\n", argv[3]); return 2; }
	if (argc > 4) {	/* zeromask: constraints with weight 0 */
		zmask = strtoul(argv[4], &err, 10);
		if (*argv[4] == 0 || *err || zmask >= 1UL << 30) {
			fprintf(stderr, "zeromask must be an integer < 2^30, not '%s'\n", argv[4]);
			return 2;
			}
		}
	asl = ASL_alloc(ASL_read_pfgh);
	if (strchr("zrRY", mode)) ihd_limit = 0;
	nl = jac0dim(argv[1], (fint)strlen(argv[1]));
	n = n_var; m = n_con;
	X0 = (real*)Malloc(n*sizeof(real));
	pfgh_read(nl, 0);
	X = (real*)Malloc((size_t)K*n*sizeof(real)); P = (real*)Malloc((size_t)K*n*sizeof(real));
	HV = (real*)Malloc((size_t)K*n*sizeof(real)); Y = (real*)Malloc((m+1)*sizeof(real));
	w = (real*)Malloc(4*n*sizeof(real)); c = (real*)Malloc((m+1)*sizeof(real));
	err = (char*)Malloc(K);
	Z = (varno_t*)Malloc((n+1)*sizeof(varno_t));
	for (i = 0; i < m; i++)
		Y[i] = i < 30 && zmask >> i & 1 ? 0. : 1. + .25*(i % 5);
	for (k = 0; k < K; k++) {
		make_point(X + (size_t)k*n);
		make_direction(X + (size_t)k*n, P + (size_t)k*n);
		}
	/* pass 1: one product per point */
	if (strchr("qQrRyY", mode))	/* first calls, no pass to find errors */
		for (k = 0; k < K; k++) {
			err[k] = 0;
			first_product_at(X + (size_t)k*n, P + (size_t)k*n, HV + (size_t)k*n, w, c);
			}
	else if (strchr("DSmohFG", mode)) {
		k = 0;
		if (mode == 'F' || mode == 'G')	/* before any other evaluation */
			err[k++] = (char)product_first_calls(X, P, HV, w);
		for (i = k; i < K; i++)
			err[i] = (char)evaluate_at(X + (size_t)i*n, w, c, 0);
		for (; k < K; k++)
			if (!err[k])
				first_product_at(X + (size_t)k*n, P + (size_t)k*n, HV + (size_t)k*n, w, c);
		}
	else
		for (k = 0; k < K; k++)
			err[k] = (char)product_at(X + (size_t)k*n, P + (size_t)k*n, HV + (size_t)k*n, w, c);
	/* pass 2: finite differences */
	for (k = 0; k < K; k++) {
		r = err[k] ? 2 : compare_at(k, X + (size_t)k*n, P + (size_t)k*n, HV + (size_t)k*n, w, bad < 5);
		if (r == 2) skipped++;
		else { checked++; bad += r; }
		}
	printf("mode=%c n=%d m=%d points=%d skipped=%d mismatching points=%d\n", mode, n, m, K, skipped, bad);
	if (bad) return 3;
	if (checked < (K < 2 ? K : 2)) {
		printf("  too few points compared: %d\n", checked);
		return 4;
		}
	return 0;
	}
