/*******************************************************************
Copyright (C) 2017, 2020 AMPL Optimization, Inc.; written by David M. Gay.

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

#include "jacpdim.h"
#include "opno2.h" /* for Varval */

#ifdef __cplusplus
extern "C" {
#endif

extern void dv_comp_ASL(EvalWorkspace*, int, int);

 int
xp_check_ASL(EvalWorkspace *ew, real *x)
{
	ASL_pfgh *asl;
	Varval *V;
	int *ov, *ove, *vm;
	linarg *la;
	real *oc, t, *vscale, *w, *xe;

	asl = (ASL_pfgh*)ew->asl;
	if (!(ew->x0kind & ASL_first_x) && !memcmp(ew->Lastx, x, x0len)) {
		++ew->stats.oldx;
		return 0;
		}
	++ew->stats.newx;
	if (ew->Derrs)
		deriv_errclear_ASL(ew);
	ew->wantderiv = want_derivs;
	memcpy(ew->Lastx, x, x0len);
	ew->nxval++;
	ew->ihdcur = 0;
	ew->x0kind = asl->i.x0kindinit;
	xe = (real*)((char*)x + x0len);
	V = (Varval*)(w = ew->w);
	vscale = asl->i.vscale;
	if ((vm = asl->i.vmap)) {
		if (vscale)
			while(x < xe)
				V[*vm++].O = *vscale++ * *x++;
		else
			while(x < xe)
				V[*vm++].O = *x++;
		}
	else {
		if (vscale)
			while(x < xe)
				(V++)->O = *vscale++ * *x++;
		else
			while(x < xe)
				(V++)->O = *x++;
		V = (Varval*)w;
		}

	for(la = asl->P.lalist; la; la = la->lnext) {
		oc = la->oc;
		ov = la->ov;
		ove = ov + la->nnz;
		t = V[*ov].O**oc;
		while(++ov < ove)
			t += V[*ov].O**++oc;
		V[la->u.v].O = t;
		}
	errno = 0;
	if (comb)
		dv_comp_ASL(ew, 0, comb);
	return 1;
	}

 int
xp2known_ew_ASL(EvalWorkspace *ew, real *X, fint *nerror)
{
	ASL *a;
	Jmp_buf err_jmp0;
	int ij, rc;

	a = ew->asl;
	ASL_CHECK(a, ASL_read_pfgh, "xp2known");
	rc = 1;
	if (ew->x0kind & ASL_xknown_ignore)
		goto ret;
	if (nerror && *nerror >= 0) {
		ew->err_jmpw = &err_jmp0;
		ij = setjmp(err_jmp0.jb);
		if ((*nerror = ij))
			goto done;
		}
	errno = 0;	/* in case f77 set errno opening files */
	rc = xp_check_ASL(ew, X);
	ew->x0kind |= ASL_x_known;
 done:
	ew->err_jmpw = 0;
 ret:
	return rc;
	}

 typedef struct Umultinfo Umultinfo;
 struct
Umultinfo {
	Umultinfo *next;
	int *ov, *ov1, *ove;
	real *oc, *oc1;
	ograd *og, *og0;
	int i, v;
	};

 static void
bigUmult(EvalWorkspace *ew, real *h, range *r, int nobj, real *ow, Umultinfo *u, real *y)
{
	ASL_pfgh *asl;
	Umultinfo *u0, *u1, *ue, **utodo, **utodoi;
	Varval *V, *v;
	int *imap, *iv;
	int i, j, n, nv, *ov, *ove;
	linarg *la, **lap;
	real *oc, *s, t, *w;

	asl = (ASL_pfgh*)ew->asl;
	w = ew->w;
	V = (Varval*)w;
	s = w + asl->P.dOscratch;
	utodo = utodoi = (Umultinfo**)(w + asl->P.utodo);
	n = r->n;
	u0 = u;
	imap = (int*)(u + n);
	iv = r->ui;
	nv = r->nv;
	for(i = 0; i < nv; i++) {
		imap[j = *iv++] = i;
		utodo[j] = 0;
		}
	lap = r->lap;
	for(i = 0; i < n; i++) {
		la = *lap++;
		u->v = la->u.v;
		u->i = i;
		u->ove = (u->ov = u->ov1 = la->ov) + la->nnz;
		u->oc = u->oc1 = la->oc;
		utodoi = utodo + *la->ov;
		u->next = *utodoi;
		*utodoi = u++;
		}
	ue = u;
	iv = r->ui;
	for(i = 0; i < nv; i++) {
		utodoi = utodo + *iv++;
		u1 = *utodoi;
		*utodoi = 0;
		for(u = u1; u; u = u->next)
			s[u->i] = *u->oc1;
		v = V + i;
		v->dO = 1.;
		v->aO = v->adO = 0.;
		pshv_prod_ASL(ew, r, nobj, ow, y);
		v->dO = 0.;
		h += i;
		for(j = 0; j <= i; j++)
			h[j] = 0.;
		while((u = u1)) {
			u1 = u->next;
			s[u->i] = 0.;
			if ((ov = u->ov1 + 1) < u->ove) {
				u->ov1 = ov;
				++u->oc1;
				utodoi = utodo + *ov;
				u->next = *utodoi;
				*utodoi = u;
				}
			}
		for(u = u0; u < ue; u++) {
			if ((t = V[u->v].aO)) {
				ov = u->ov;
				oc = u->oc;
				ove = u->ove;
				for(; ov < ove && (j = imap[*ov]) <= i; ++ov)
					h[j] += t**oc++;
				}
			}
		}
	}

 void
hvpinit_nc_ASL(EvalWorkspace *ew, int ndhmax, int nobj, real *ow, real *y)
{
	ASL_pfgh *asl;
	Ihinfo *ihi;
	Umultinfo *u;
	Varval *V, *v;
	int i, ihc, n1, *ov, *ove, *ui, *uie;
	linarg *la, **lap, **lap1, **lape;
	range *r;
	real *h, *oc, *s, *si, *w;

	asl = (ASL_pfgh*)ew->asl;
	ew->nhvprod = 0;
	ihc = 0;
	if (ndhmax > asl->P.ihdmax)
		ndhmax = asl->P.ihdmax;
	if ((ew->ndhmax = ndhmax) <= 0)
		goto done;
	if (!(ihi = asl->P.ihi1) || ndhmax < asl->P.ihdmin)
		return;
	if (nobj < 0 || nobj >= n_obj)
		nobj = -1;
	w = ew->w;
	V = (Varval*)w;
	s = w + asl->P.dOscratch;
	u = 0;
	for(ihc = 0; ihi->ihd <= ndhmax; ihi = ihi->next) {
		ihc = ihi->ihd;
		for(r = ihi->r; r; r = r->rlist.prev) {
			for(ui = r->ui, uie = ui + r->nv; ui < uie; ++ui) {
				v = V + *ui;
				v->aO = v->dO = v->adO = 0.;
				}
			h = w + r->hest;
			if ((n1 = r->n) < r->nv) {
				si = s;
				lape = lap = r->lap;
				for(i = 0; i < n1; i++) {
					*si = 1.;
					la = *lape++;
					oc = la->oc;
					ov = la->ov;
					for(ove = ov + la->nnz; ov < ove; ++ov) {
						v = V + *ov;
						v->dO = *oc++;
						v->aO = v->adO = 0.;
						}
					pshv_prod_ASL(ew, r, nobj, ow, y);
					for(ov = la->ov; ov < ove; ++ov)
						V[*ov].dO = 0.;
					*si++ = 0;
					lap1 = lap;
					do *h++ = V[(*lap1++)->u.v].aO;
					   while(lap1 < lape);
					}
				}
			else {
				if (!u) {
					if (!(i = asl->P.krnmax))
						asl->P.krnmax = i =
							htcl(asl->P.rnmax*sizeof(Umultinfo)
								+ n_var*sizeof(int));
					u = (Umultinfo*)new_mblk(i);
					}
				bigUmult(ew, h, r, nobj, ow, u, y);
				}
			}
		}
	if (u)
		del_mblk(u);
 done:
	ew->ihdcur = ihc;
	}

 void
hvpinit_ew_ASL(EvalWorkspace *ew, int ndhmax, int nobj, real *ow, real *y)
{
	ASL_pfgh *asl = (ASL_pfgh*)ew->asl;
	ASL_CHECK(((ASL*)asl), ASL_read_pfgh, "hvpinit");
	xpsg_check_ASL(ew, nobj, ow, y);
	hvpinit_nc_ASL(ew, ndhmax, nobj, ow, y);
	}

/* Variant of hvpinit_ew() that has a final nerror argument, working
   similarly to the final nerror argument to objval_(), etc. */

 void
hvpinite_ew_ASL(EvalWorkspace *ew, int ndhmax, int nobj, real *ow, real *y, fint *nerror)
{
	ASL_pfgh *asl;
	Jmp_buf **Jp, *Jsave, b;

	asl = (ASL_pfgh*)ew->asl;
	ASL_CHECK(((ASL*)asl), ASL_read_pfgh, "hvpinite");
	Jp = !nerror || *nerror >= 0 ? &ew->err_jmpw : &ew->err_jmpw1;
	Jsave = *Jp;
	*Jp = &b;
	*nerror = 0;
	if (setjmp(b.jb))
		*nerror = 1;
	else
		hvpinit_nc_ASL(ew, ndhmax, nobj, ow, y);
	*Jp = Jsave;
	}

#ifdef __cplusplus
}
#endif
