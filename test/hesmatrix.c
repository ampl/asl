/* Hessian entry points x prior evaluations x weights (both libraries).
 *
 * usage: hesmatrix stub WHO
 * WHO names the functions with an undefined derivative at the start point:
 * "o" (objective 0), "c" (the constraints) or "oc" (both); with "n" the
 * derivative becomes Inf or NaN with WANT_INFNAN (pow"; not sqrt').  For every
 * combination of
 *   - entry point: sphes without sphsetup, sphsetup + sphes, sphsetup +
 *     sphese, hvcomp, sphsetup + hvcomp, duthes, duthese, fullhes, fullhese;
 *   - prior evaluation at x: none (the Hessian call evaluates first),
 *     objval, conval, both, xknown;
 *   - weights: (nobj, ow, y) = (0, -, 1), (0, -, -), (0, {0}, 1),
 *     (-1, {1}, -), (-1, -, 1), (-1, {1}, 0),
 * the stub is read again (a fresh state), the calls are made, and an
 * evaluation error must be reported exactly when a function named by WHO
 * has a nonzero weight.  Built with WANT_INFNAN and "n", no error is
 * expected.  Without an error, the result must agree with fullhes for the
 * same weights (sphes: every element of the upper triangle, so missing
 * structure shows; duthes: the upper triangle; hvcomp: column 0).
 * Exit status: 3 if an outcome differs, 0 otherwise.
 */
#include <math.h>
#include <setjmp.h>
#include <string.h>
#include "asl_pfgh.h"

enum { ENTRIES = 9, PRIORS = 5, WEIGHTS = 6 };
static const char *entry_name[ENTRIES] = { "sphes", "setup+sphes",
	"setup+sphese", "hvcomp", "setup+hvcomp", "duthes", "duthese",
	"fullhes", "fullhese" };
static const char *prior_name[PRIORS] = { "none", "objval", "conval",
	"both", "xknown" };

/* The weights of combination w: nobj, ow (0 or a pointer), y (0, ones or
   zeros); the objective and the constraints are used when their weight
   is nonzero */
static void weights(int w, int *nobj, real **ow, real **y, real *zero1,
	real *one1, real *ones, real *zeros, int *uo, int *uc)
{
	static const int nob[WEIGHTS] = { 0, 0, 0, -1, -1, -1 };
	*nobj = nob[w];
	*ow = w == 2 ? zero1 : w == 3 || w == 5 ? one1 : 0;
	*y = w == 0 || w == 2 || w == 4 ? ones : w == 5 ? zeros : 0;
	*uo = w == 0 || w == 1 || w == 3 || w == 5;
	*uc = w == 0 || w == 2 || w == 4;
	}

/* |a - b| small relative to their size (both finite) */
 static int
differ(real a, real b)
{
	if (!isfinite(a) || !isfinite(b))
		return 0;	/* WANT_INFNAN: not compared */
	return fabs(a - b) > 1e-10*(1. + fabs(a) + fabs(b));
	}

/* Compare the result of entry e (H, upper triangle by the layout of e)
   with the dense fullhes F; returns 1 if they differ */
 static int
compare(ASL *asl, int e, real *H, real *F)
{
	int i, j, k, n = n_var;
	real *D;

	if (e <= 2) {	/* sphes: sparse upper triangle, all of it */
		D = (real*)Malloc(n*n*sizeof(real));
		memset(D, 0, n*n*sizeof(real));
		for(j = 0; j < n; ++j)
			for(k = sputinfo->hcolstarts[j];
			    k < sputinfo->hcolstarts[j+1]; ++k)
				D[sputinfo->hrownos[k] + j*n] = H[k];
		for(j = 0, k = 0; j < n && !k; ++j)
			for(i = 0; i <= j; ++i)
				if (differ(D[i + j*n], F[i + j*n]))
					k = 1;
		free(D);
		return k;
		}
	if (e <= 4) {	/* hvcomp: column 0 */
		for(i = 0; i < n; ++i)
			if (differ(H[i], F[i]))
				return 1;
		return 0;
		}
	if (e <= 6) {	/* duthes: packed upper triangle by columns */
		for(j = 0; j < n; ++j)
			for(i = 0; i <= j; ++i)
				if (differ(H[i + j*(j+1)/2], F[i + j*n]))
					return 1;
		}
	return 0;	/* fullhes itself */
	}

/* One combination; returns 1 if an evaluation error was reported, sets
   *val to 1 if the result differs from fullhes */
static int run(const char *stub, int e, int p, int w, int *val)
{
	ASL *asl;
	FILE *nl;
	Jmp_buf jb;
	fint ne = 0;
	real *H, *c, *ones, *zeros, *ow, *y, *v, zero1 = 0., one1 = 1.;
	int i, nobj, uo, uc, err;

	asl = ASL_alloc(ASL_read_pfgh);
	nl = jac0dim((char*)stub, (fint)strlen(stub));
	X0 = (real*)Malloc(n_var*sizeof(real));
	pfgh_read(nl, 0);
	H = (real*)Malloc((n_var*(n_var + 2) + 2*n_con + 2)*sizeof(real));
	v = H + n_var*n_var;
	c = v + n_var + 1;
	ones = c + n_con;
	zeros = (real*)Malloc((n_con + 1)*sizeof(real));
	for(i = 0; i < n_con; ++i) {
		ones[i] = 1.;
		zeros[i] = 0.;
		}
	weights(w, &nobj, &ow, &y, &zero1, &one1, ones, zeros, &uo, &uc);
	err = 0;
	err_jmp = &jb;	/* an error longjmps here silently */
	if (setjmp(jb.jb))
		err = 1;
	else {
		if (p == 1 || p == 3)
			objval(0, X0, &ne);
		if ((p == 2 || p == 3) && n_con)
			conval(X0, c, &ne);
		if (p == 4)
			xknown(X0);
		switch(e) {
		  case 0: sphes(H, nobj, ow, y); break;
		  case 1:
		  case 2:
			sphsetup(nobj, ow != 0, y != 0, 1);
			if (e == 1)
				sphes(H, nobj, ow, y);
			else if (sphese(H, nobj, ow, y, &ne), ne)
				err = 1;
			break;
		  case 3:
		  case 4:
			if (e == 4)
				sphsetup(nobj, ow != 0, y != 0, 1);
			for(i = 0; i < n_var; ++i)
				v[i] = i == 0;
			hvcomp(H, v, nobj, ow, y);
			break;
		  case 5: duthes(H, nobj, ow, y); break;
		  case 6: if (duthese(H, nobj, ow, y, &ne), ne) err = 1; break;
		  case 7: fullhes(H, n_var, nobj, ow, y); break;
		  default: if (fullhese(H, n_var, nobj, ow, y, &ne), ne) err = 1;
		  }
		}
	*val = 0;
	if (!err) {	/* fullhes in the same state for the same weights */
		real *F = (real*)Malloc(n_var*n_var*sizeof(real));
		err_jmp = &jb;
		if (!setjmp(jb.jb)) {
			fullhes(F, n_var, nobj, ow, y);
			*val = compare(asl, e, H, F);
			}
		free(F);
		}
	err_jmp = 0;
	free(zeros);
	free(H);
	ASL_free(&asl);
	return err;
	}

int main(int argc, char **argv)
{
	int bad, e, err, expect, n, o, c, p, w, uo, uc, nobj, val;
	real *ow, *y, z = 0, *zp = &z;

	if (argc < 3) {
		fprintf(stderr, "usage: hesmatrix stub o|c|oc[n]\n");
		return 2;
		}
	o = strchr(argv[2], 'o') != 0;
	c = strchr(argv[2], 'c') != 0;
	bad = n = 0;
	for(e = 0; e < ENTRIES; ++e)
	 for(p = 0; p < PRIORS; ++p)
	  for(w = 0; w < WEIGHTS; ++w) {
		weights(w, &nobj, &ow, &y, zp, zp, zp, zp, &uo, &uc);
		expect = (o && uo) || (c && uc);
#ifdef WANT_INFNAN
		if (strchr(argv[2], 'n'))
			expect = 0;
#endif
		err = run(argv[1], e, p, w, &val);
		++n;
		if (val) {
			++bad;
			printf("VALUES %s, prior %s, weights %d: differ from fullhes\n",
				entry_name[e], prior_name[p], w);
			}
		else if (err != expect) {
			++bad;
			printf("MISMATCH %s, prior %s, weights %d: %s\n",
				entry_name[e], prior_name[p], w,
				err ? "error reported" : "no error");
			}
		}
	printf("hesmatrix: %d combinations, %d mismatches\n", n, bad);
	return bad ? 3 : 0;
	}
