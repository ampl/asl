/* sphes when its setup runs inside the call (solvers2, pfgh_read).
 *
 * usage: sphcheck stub [e]
 * Evaluates the functions at the start point, then calls sphes without
 * sphsetup, so that sphes sets up its structure itself, and then calls
 * sphes again at the same point: both results must be equal.  The setup
 * uses cells that the evaluations filled (e.g. the branch an if took), so
 * the first result was wrong when the setup came after the evaluations.
 * With "e", the objective is expected to have an evaluation error at the
 * start point, and sphsetup after it must return (it used to end the run).
 * With "s", "d" or "f", a second derivative is undefined at the start
 * point: after the function evaluations, sphsetup and sphes ("s"), duthes
 * ("d") or fullhes ("f") must report the evaluation error (they used to
 * return a Hessian without that term, or with Inf or stale values); built
 * with WANT_INFNAN, they report no error and the first diagonal element
 * is +Infinity.
 * Exit status: 3 if the results differ, 4 if no error was reported where
 * one was expected, 0 otherwise (also when an evaluation error leaves
 * nothing to compare).
 */
#include <string.h>
#include "asl_pfgh.h"

int main(int argc, char **argv) {
	ASL *asl;
	FILE *nl;
	fint ne = 0, nz;
	real *c, *g, *H1, *H2, d, t, *x, *y;
	int i, j, nobj;
	if (argc < 2) { fprintf(stderr, "usage: sphcheck stub [e|s|d|f]\n"); return 2; }
	asl = ASL_alloc(ASL_read_pfgh);
	nl = jac0dim(argv[1], (fint)strlen(argv[1]));
	X0 = (real*)Malloc(n_var*sizeof(real));
	pfgh_read(nl, 0);
	x = X0;
	g = (real*)Malloc((n_var + n_con + 1)*sizeof(real));
	c = g + n_var;
	y = (real*)Malloc((n_con + 1)*sizeof(real));
	for (i = 0; i < n_con; i++) y[i] = 1. + .25*(i % 5);
	nobj = n_obj > 0 ? 0 : -1;
	if (argc > 2 && *argv[2] == 'e') {
		if (n_obj) objval(0, x, &ne);
		nz = sphsetup(nobj, 0, n_con > 0, 0);
		printf("evaluation error %d, then sphsetup: nnz=%d\n", (int)ne, (int)nz);
		return 0;
		}
	if (argc > 2 && strchr("sdf", *argv[2])) {	/* Hessian error expected */
		if (n_obj) objval(0, x, &ne);
		if (n_con) conval(x, c, &ne);
		H1 = (real*)Malloc(((size_t)n_var*n_var + 1)*sizeof(real));
		if (*argv[2] == 's') {
			sphsetup(nobj, 0, n_con > 0, 1);
			sphes(H1, nobj, 0, n_con ? y : 0);
			}
		else if (*argv[2] == 'd')
			duthes(H1, nobj, 0, n_con ? y : 0);
		else
			fullhes(H1, n_var, nobj, 0, n_con ? y : 0);
		/* With WANT_INFNAN no error: H[0] (the first diagonal */
		/* element in all three layouts) is then +Infinity */
		printf("no evaluation error reported; H[0] = %g: %s\n", H1[0],
			H1[0] > 0. && H1[0] * 0. != 0. ? "H0_PLUS_INF" : "H0_OTHER");
		return 4;
		}
	if (n_obj) objval(0, x, &ne);
	if (n_obj && !ne) objgrd(0, x, g, &ne);
	if (n_con && !ne) conval(x, c, &ne);
	for (i = 0; i < n_con && !ne; i++) congrd(i, x, g, &ne);
	if (ne) { printf("evaluation error at the start point: not compared\n"); return 0; }
	/* at most the full matrix: no sphsetup, so both triangles */
	H1 = (real*)Malloc(2*((size_t)n_var*n_var + 1)*sizeof(real));
	H2 = H1 + (size_t)n_var*n_var + 1;
	sphes(H1, nobj, 0, n_con ? y : 0);	/* no sphsetup: sets up inside */
	nz = sputinfo->hcolstarts[n_var] - sputinfo->hcolstarts[0];
	sphes(H2, nobj, 0, n_con ? y : 0);
	for (j = 0, d = 0.; j < nz; j++) {
		t = fabs(H1[j] - H2[j]) / (1. + fabs(H2[j]));
		if (t > d) d = t;
		}
	printf("n=%d nnz=%d max relative difference of the first sphes: %.3g\n",
		n_var, (int)nz, d);
	return d > 1e-12 ? 3 : 0;
	}
