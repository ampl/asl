/*******************************************************************
Copyright (C) 2016 AMPL Optimization, Inc.; written by David M. Gay.

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
add_op(real *w, real *H, psg_elem *g, real t, fint n)
{
	int i, j, k, nov, *ov;
	real *Hj, *oc, t1;

	ov = g->ov;
	oc = w + g->gm + 3;
	nov = g->nov;
	for(i = 0; i < nov; ++i) {
		if ((t1 = t * oc[i])) {
			j = ov[i];
			Hj = H + n*j;
			for(k = 0;; ++k) {
				Hj[ov[k]] += t1*oc[k];
				if (k == i)
					break;
				}
			}
		}
	}

 void
fullhes_ew_ASL(EvalWorkspace *ew, real *H, fint LH, int nobj, real *ow, real *y)
{
	/* full Hessian */
	ASL *a;
	ASL_pfgh *asl;
	Varval *V, *v;
	int i, j, n, no, noe, *ov, *ov1, *ov1e, *ove;
	linarg *la, **lap, **lap1, **lape;
	ps_func *p, *pe;
	psg_elem *g, *ge;
	range *r, *r0;
	real *Hi, *Hj, *Hje, *cscale, g2, *oc, *oc1, *owi, *s, *si, t, t1, *w;

	asl = (ASL_pfgh*)(a = ew->asl);
	ASL_CHECK(a, ASL_read_pfgh, "fullhes");
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

	/* First compute upper triangle */

	for(Hj = H, i = 1; i <= n; Hj += LH - i++)
		for(Hje = Hj + i; Hj < Hje; )
			*Hj++ = 0;

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
			ove = ov + la->nnz;
			oc = la->oc;
			while(ov < ove) {
				j = *ov++;
				t = *oc++;
				Hj = H + LH*j;
				for(lap1 = r->lap; lap1 < lape; ) {
					la = *lap1++;
					v = V + la->u.v;
					if (!(t1 = t * v->aO))
						continue;
					oc1 = la->oc;
					ov1 = la->ov;
					ov1e = ov1 + la->nnz;
					for(; ov1 < ov1e && (i = *ov1) <= j; ++ov1)
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
						add_op(w, H, g, t*g2, LH);
				}
	if (asl->P.ncongroups && y) {
		cscale = asl->i.lscale;
		p = asl->P.cps;
		for(pe = p + n_con; p < pe; p++, y++)
			if ((t = cscale ? *cscale++ * *y : *y))
				for(g = p->g, ge = p->ge; g < ge; g++)
					if ((g2 = w[g->gm + 2]))
						add_op(w, H, g, t*g2, LH);
		}

	if ((s = asl->i.vscale)) {
		Hi = H;
		for(i = 0; i < n; i++, Hi += LH) {
			t = s[i];
			for(j = 0; j <= i; j++)
				Hi[j] *= t * s[j];
			}
		}

	/* Now copy upper triangle to lower triangle. */

	for(Hj = H, i = 1; i < n; i++) {
		Hi = H + i;
		Hj = H + i*LH;
		Hje = Hj + i;
		while(Hj < Hje) {
			*Hi = *Hj++;
			Hi += LH;
			}
		}
	}

/* Variant of fullhes() that has a final nerror argument, working
   similarly to the final nerror argument to objval_(), etc. */

 void
fullhese_ew_ASL(EvalWorkspace *ew, real *H, fint LH, int nobj, real *ow, real *y, fint *nerror)
{
	Jmp_buf **Jp, *Jsave, b;

	Jp = !nerror || *nerror >= 0 ? &ew->err_jmpw : &ew->err_jmpw1;
	Jsave = *Jp;
	*Jp = &b;
	*nerror = 0;
	if (setjmp(b.jb))
		*nerror = 1;
	else
		fullhes_ew_ASL(ew, H, LH, nobj, ow, y);
	*Jp = Jsave;
	}
