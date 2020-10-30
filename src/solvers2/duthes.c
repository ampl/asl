/*******************************************************************
Copyright (C) 2016, 2020 AMPL Optimization, Inc.; written by David M. Gay.

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

 static void
add_op(real *w, real *H, psg_elem *g, real t)
{
	int i, j, k, nov, *ov;
	real *Hj, *oc, t1;

	ov = g->ov;
	oc = w + g->gm + 3;
	nov = g->nov;
	for(i = 0; i < nov; ++i) {
		if ((t1 = t * oc[i])) {
			j = ov[i];
			Hj = H + ((j*(j+1))>>1);
			for(k = 0;; ++k) {
				Hj[ov[k]] += t1*oc[k];
				if (k == i)
					break;
				}
			}
		}
	}

 void
duthes_ew_ASL(EvalWorkspace *ew, real *H, int nobj, real *ow, real *y)
{
	/* dense upper triangle of Hessian */
	ASL *a;
	ASL_pfgh *asl;
	Varval *V, *v;
	int i, j, n, no, noe, *ov, *ov1, *ov1e, *ove;
	linarg *la, **lap, **lap1, **lape;
	ps_func *p, *pe;
	psg_elem *g, *ge;
	range *r, *r0;
	real *Hj, *cscale, g2, *oc, *oc1, *owi, *s, *si, t, t1, *w;

	asl = (ASL_pfgh*)(a = ew->asl);
	ASL_CHECK(a, ASL_read_pfgh, "duthes");
	xpsg_check_ASL(ew, nobj, ow, y);
	if (nobj >= 0 && nobj < n_obj) {
		no = nobj;
		noe = no + 1;
		owi = ow ? ow + no : &edag_one_ASL;
		}
	else {
		nobj = -1;
		no = noe = 0;
		if ((owi = ow))
			noe = n_obj;
		}

	n = c_vars >= o_vars ? c_vars : o_vars;
	if (n <= 0)
		return;
	++ew->stats.hesmat;
	V = (Varval*)(w = ew->w);
	s = w + asl->P.dOscratch;
	memset(H, 0, (n*(n+1) >> 1) * sizeof(real));

	r0 = (range*)&asl->P.rlist;
	for(r = asl->P.rlist.next; r != r0; r = r->rlist.next) {
		if ((j = r->n) <= 0)
			continue;
		lap = r->lap;
		lape = lap + j;
		si = s;
		while(lap < lape) {
			*si = 1;
			pshv_prod_ASL(ew, r, nobj, ow, y);
			*si++ = 0;
			la = *lap++;
			ov = la->ov;
			oc = la->oc;
			ove = ov + la->nnz;
			while(ov < ove) {
				t = *oc++;
				j = *ov++;
				Hj = H + ((j*(j+1))>>1);
				for(lap1 = r->lap; lap1 < lape; ) {
					la = *lap1++;
					v = V + la->u.v;
					if (!(t1 = t * v->aO))
						continue;
					ov1 = la->ov;
					ov1e = ov1 + la->nnz;
					for(oc1 = la->oc;ov1 < ov1e && (i = *ov1) <= j; ++ov1)
						Hj[i] += t1**oc1++;
					}
				}
			}
		}
	if (asl->P.nobjgroups)
		for(; no < noe; no++)
			if ((t = *owi++)) {
				p = asl->P.ops + no;
				g = p->g;
				for(ge = p->ge; g < ge; g++)
					if ((g2 = w[g->gm + 2]))
						add_op(w, H, g, t*g2);
				}
	if (asl->P.ncongroups && y) {
		cscale = asl->i.lscale;
		p = asl->P.cps;
		for(pe = p + n_con; p < pe; p++, y++)
			if ((t = cscale ? *cscale++ * *y : *y))
				for(g = p->g, ge = p->ge; g < ge; g++)
					if ((g2 = w[g->gm + 2]))
						add_op(w, H, g, t*g2);
		}
	if ((s = asl->i.vscale))
		for(i = 0; i < n; i++) {
			t = s[i];
			for(j = 0; j <= i; j++)
				*H++ *= t * s[j];
			}
	}

/* Variant of duthes that has a final nerror argument, working
   similarly to the final nerror argument to objval_(), etc. */

 void
duthese_ew_ASL(EvalWorkspace *ew, real *H, int nobj, real *ow, real *y, fint *nerror)
{
	Jmp_buf **Jp, *Jsave, b;

	Jp = !nerror || *nerror >= 0 ? &ew->err_jmpw : &ew->err_jmpw1;
	Jsave = *Jp;
	*Jp = &b;
	*nerror = 0;
	if (setjmp(b.jb))
		*nerror = 1;
	else
		duthes_ew_ASL(ew, H, nobj, ow, y);
	*Jp = Jsave;
	}
