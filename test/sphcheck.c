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
 * Exit status: 3 if the results differ, 0 otherwise (also when an
 * evaluation error leaves nothing to compare).
 */
#include "asl_pfgh.h"

int main(int argc, char **argv) {
	ASL *asl;
	FILE *nl;
	fint ne = 0, nz;
	real *c, *g, *H1, *H2, d, t, *x, *y;
	int i, j, nobj;
	if (argc < 2) { fprintf(stderr, "usage: sphcheck stub [e]\n"); return 2; }
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
