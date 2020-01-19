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

/* Test program for nqpcheck:  read one or more problems given as  */
/* command-line arguments and apply nqpcheck to all objectives and */
/* constraints: nq >= 0 means a quadratic form with nq nonzeros in */
/* whole Q matrix (not just its upper triangle), and nq = -1 means */
/* a nonlinear objective or constraint.  This also illustrates use */
/* of ASL_free() to dispose of an ASL pointer.			   */

#include "asl.h"

 static real
qterm(ASL *asl, fint *colq, fint *rowq, real *delsq)
{
	real t, t1, *x, *x0, *xe;
	fint *rq0, *rqe;

	t = 0.;
	x0 = x = X0;
	xe = x + n_var;
	rq0 = rowq;
	while(x < xe) {
		t1 = *x++;
		rqe = rq0 + *++colq;
		while(rowq < rqe)
			t += t1*x0[*rowq++]**delsq++;
		}
	return 0.5 * t;
	}

 int
main(int argc, char **argv)
{
	char *s;
	FILE *nl;
	int i, j, nz;
	fint **colqp, ne, nq, **rowqp, *z, *z0;
	real *LU, **delsqp, t;
	ASL *asl;

	if (argc < 2) {
		printf("usage: %s stub [stub ...]\n\tto apply nqpcheck to all objectives and constraints.\n", argv[0]);
		return 0;
		}

	while(s = *++argv) {
		if (argc > 3)
			printf("\n%s:\n", s);
		asl = ASL_alloc(ASL_read_fg);
		nl = jac0dim(s, (long)strlen(s));
		X0 = (real*)M1alloc(n_var*sizeof(real));
		qp_read(nl,0);
		nz = n_obj + n_con;
		colqp = (fint**)Malloc(nz*(sizeof(fint)
					   + 2*sizeof(fint*)
					   + sizeof(real)));
		rowqp = colqp + nz;
		delsqp = (real **)(rowqp + nz);
		z = z0 = (fint*)(delsqp + nz);
		for(i = 0; i < n_obj; i++)
			*z++ = nqpcheck(i, rowqp+i, colqp+i, delsqp+i);
		for(j = 0; j < n_con; i++, j++)
			*z++ = nqpcheck(-1-j, rowqp+i, colqp+i, delsqp+i);
		qp_opify();
		for(i = 0; i < n_obj; i++) {
			ne = 0;
			t = objval(i, X0, &ne);
			if (ne) {
				printf("Error evaluating objective %d\n", i);
				continue;
				}
			nq = z0[i];
			if (nq > 0)
				t += qterm(asl, colqp[i], rowqp[i], delsqp[i]);
			printf("Objective %d: nq = %d, value = %.g\n",
					i, (int)nq, t);
			}
		LU = LUrhs;
		for(j = 0; j < n_con; i++, j++) {
			ne = 0;
			t = conival(j, X0, &ne);
			if (ne) {
				printf("Error evaluating constraint %d\n", j);
				continue;
				}
			nq = z0[i];
			if (nq > 0)
				t += qterm(asl, colqp[i], rowqp[i], delsqp[i]);
			printf("Constraint body %d: nq = %d, value = %.g\n",
				j, (int)nq, t);
			printf("\t\tLslack = %g, Uslack = %g\n",
				t - LU[0], LU[1] - t);
			LU += 2;
			}
		free(colqp);
		ASL_free((ASL**)&asl);
		}
	return 0;
	}
