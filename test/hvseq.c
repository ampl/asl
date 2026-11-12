/* Hessian-vector products over a sequence of points, one product per point.
 *
 * usage: hvseq stub [npoints]
 * Pass 1 goes through npoints (default 12) random points near the start
 * point; at each it evaluates functions and gradients and calls hvcomp
 * once.  Pass 2 compares every product with central differences of the
 * gradient of the Lagrangian.  A product computed with the Hessian of an
 * earlier point shows up as a mismatch.
 * Works with solvers and solvers2 (pfgh_read).  Exit status 3 on mismatch.
 */
#include "asl_pfgh.h"

static ASL *asl;
static real *Y;

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

int main(int argc, char **argv) {
	FILE *nl; int K = argc > 2 ? atoi(argv[2]) : 12, n, m, i, j, k, bad = 0, skipped = 0;
	unsigned long long s = 12345;
	real *X, *HV, *p, *w, *x, *hv, ow = 1., h = 1e-6, d, sc, u;
	char *err; fint ne;
	if (argc < 2) { fprintf(stderr, "usage: hvseq stub [npoints]\n"); return 2; }
	asl = ASL_alloc(ASL_read_pfgh);
	nl = jac0dim(argv[1], (fint)strlen(argv[1]));
	n = n_var; m = n_con;
	X0 = (real*)Malloc(n*sizeof(real));
	pfgh_read(nl, 0);
	X = (real*)Malloc((size_t)K*n*sizeof(real)); HV = (real*)Malloc((size_t)K*n*sizeof(real));
	Y = (real*)Malloc((m+1)*sizeof(real)); p = (real*)Malloc(n*sizeof(real));
	w = (real*)Malloc(4*n*sizeof(real)); err = (char*)Malloc(K);
	for (i = 0; i < m; i++) Y[i] = 1. + .25*(i % 5);
	for (j = 0; j < n; j++) p[j] = 1. - .125*(j % 7);
	for (k = 0; k < K; k++)
		for (j = 0; j < n; j++) {
			s = s*6364136223846793005ULL + 1442695040888963407ULL;
			u = (real)((s >> 33) % 2001)/1000. - 1.;	/* in [-1,1] */
			X[(size_t)k*n + j] = X0[j] + .3*u*(1. + fabs(X0[j]));
			}
	/* pass 1: one product per point */
	for (k = 0; k < K; k++) {
		x = X + (size_t)k*n; hv = HV + (size_t)k*n; ne = 0;
		if (n_obj) objval(0, x, &ne);
		if (n_obj && !ne) objgrd(0, x, w, &ne);
		if (m && !ne) conval(x, w + n, &ne);
		for (i = 0; i < m && !ne; i++) congrd(i, x, w, &ne);
		if ((err[k] = ne != 0)) { skipped++; continue; }
		hvcomp(hv, p, n_obj ? 0 : -1, n_obj ? &ow : 0, m ? Y : 0);
		}
	/* pass 2: finite differences */
	for (k = 0; k < K; k++) {
		if (err[k]) continue;
		x = X + (size_t)k*n; hv = HV + (size_t)k*n;
		for (j = 0; j < n; j++) w[j] = x[j] + h*p[j];
		if (lagrgrad(w, w + n, w + 3*n)) { skipped++; continue; }
		for (j = 0; j < n; j++) w[j] = x[j] - h*p[j];
		if (lagrgrad(w, w + 2*n, w + 3*n)) { skipped++; continue; }
		for (j = 0; j < n; j++) {
			d = (w[n+j] - w[2*n+j])/(2*h); sc = 1. + fabs(d) + fabs(hv[j]);
			if (!(fabs(d - hv[j]) <= 1e-4*sc)) {
				if (++bad <= 5)
					printf("  MISMATCH point %d, hv[%d]: hvcomp %.10g  fd %.10g\n", k, j, hv[j], d);
				break;
				}
			}
		}
	printf("n=%d m=%d points=%d skipped=%d mismatching points=%d\n", n, m, K, skipped, bad);
	return bad ? 3 : 0;
	}
