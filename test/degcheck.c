/* degree_ASL against mqpcheckv: solvers2 with pfgh_read; with DEGCHECK_FG
 * the solvers library with qp_read, as its QP drivers read (its
 * degree_ASL needs the operator codes, which fg_read replaces).
 *
 * usage: degcheck stub
 * For every objective and every nonlinear constraint, a function that
 * degree_ASL calls quadratic (2) must be accepted by mqpcheckv, which
 * drivers then use to pass the quadratic terms to the solver (e.g.
 * Knitro); a rejection would stop such a driver.
 * Exit status: 3 if degree_ASL and mqpcheckv disagree.
 */
#ifdef DEGCHECK_FG
#include "asl.h"
#else
#include "asl_pfgh.h"
#endif

int main(int argc, char **argv) {
	ASL *asl;
	FILE *nl;
	QPinfo *qpi;
	void *vd = 0, *vp = 0;
	int co, d, bad = 0, nq = 0;
	if (argc < 2) { fprintf(stderr, "usage: degcheck stub\n"); return 2; }
	nl = 0;
#ifdef DEGCHECK_FG
	asl = ASL_alloc(ASL_read_fg);
	nl = jac0dim(argv[1], (fint)strlen(argv[1]));
	qp_read(nl, 0);
#else
	asl = ASL_alloc(ASL_read_pfgh);
	nl = jac0dim(argv[1], (fint)strlen(argv[1]));
	pfgh_read(nl, ASL_cc_simplify | ASL_sep_U_arrays);
#endif
	for (co = -nlc; co < n_obj; co++) {	/* co < 0: constraint -1 - co */
		d = degree_ASL(asl, co, &vd);
		if (d != 2)
			continue;
		++nq;
		qpi = 0;
		if (mqpcheckv(co, &qpi, &vp) < 0) {
			printf("  %s %d: degree_ASL 2, rejected by mqpcheckv\n",
				co < 0 ? "constraint" : "objective", co < 0 ? -1 - co : co);
			++bad;
			}
		if (qpi)
			free(qpi);
		}
	printf("n_obj=%d nlc=%d quadratic=%d disagreements=%d\n", n_obj, nlc, nq, bad);
	return bad ? 3 : 0;
	}
