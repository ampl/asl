/* Hessian-vector products over a sequence of points, one product per point.
 *
 * usage: hvseq stub [npoints]
 * Pass 1 goes through npoints (default 12) random points near the start
 * point; at each it evaluates functions and gradients and calls hvcomp
 * once.  Pass 2 compares every product with central differences of the
 * gradient of the Lagrangian.  A product computed with the Hessian of an
 * earlier point shows up as a mismatch.
 * Works with solvers and solvers2 (pfgh_read).
 * Exit status: 3 on mismatch, 4 when fewer than two points (one if
 * npoints is 1) could be compared, so a run that skips every point fails.
 */
#include "asl_pfgh.h"

static ASL *asl;
static real *Y;
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

/* Values and gradients at x, then one product; returns 1 on an evaluation
 * error.  w: n_var reals, c: n_con reals. */
static int product_at(real *x, real *p, real *hv, real *w, real *c) {
	fint ne = 0; int i;
	if (n_obj) objval(0, x, &ne);
	if (n_obj && !ne) objgrd(0, x, w, &ne);
	if (n_con && !ne) conval(x, c, &ne);
	for (i = 0; i < n_con && !ne; i++) congrd(i, x, w, &ne);
	if (ne) return 1;
	/* no weight vector: it would need one weight per objective */
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
			printf("  MISMATCH point %d, hv[%d]: hvcomp %.10g  fd %.10g\n", k, j, hv[j], d);
		return 1;
		}
	return 0;
	}

int main(int argc, char **argv) {
	FILE *nl; int K = argc > 2 ? atoi(argv[2]) : 12, n, m, i, k, r;
	int bad = 0, skipped = 0, checked = 0;
	real *X, *P, *HV, *w, *c; char *err;
	if (argc < 2) { fprintf(stderr, "usage: hvseq stub [npoints]\n"); return 2; }
	asl = ASL_alloc(ASL_read_pfgh);
	nl = jac0dim(argv[1], (fint)strlen(argv[1]));
	n = n_var; m = n_con;
	X0 = (real*)Malloc(n*sizeof(real));
	pfgh_read(nl, 0);
	X = (real*)Malloc((size_t)K*n*sizeof(real)); P = (real*)Malloc((size_t)K*n*sizeof(real));
	HV = (real*)Malloc((size_t)K*n*sizeof(real)); Y = (real*)Malloc((m+1)*sizeof(real));
	w = (real*)Malloc(4*n*sizeof(real)); c = (real*)Malloc((m+1)*sizeof(real));
	err = (char*)Malloc(K);
	for (i = 0; i < m; i++) Y[i] = 1. + .25*(i % 5);
	for (k = 0; k < K; k++) {
		make_point(X + (size_t)k*n);
		make_direction(X + (size_t)k*n, P + (size_t)k*n);
		}
	/* pass 1: one product per point */
	for (k = 0; k < K; k++)
		err[k] = (char)product_at(X + (size_t)k*n, P + (size_t)k*n, HV + (size_t)k*n, w, c);
	/* pass 2: finite differences */
	for (k = 0; k < K; k++) {
		r = err[k] ? 2 : compare_at(k, X + (size_t)k*n, P + (size_t)k*n, HV + (size_t)k*n, w, bad < 5);
		if (r == 2) skipped++;
		else { checked++; bad += r; }
		}
	printf("n=%d m=%d points=%d skipped=%d mismatching points=%d\n", n, m, K, skipped, bad);
	if (bad) return 3;
	if (checked < (K < 2 ? K : 2)) {
		printf("  too few points compared: %d\n", checked);
		return 4;
		}
	return 0;
	}
