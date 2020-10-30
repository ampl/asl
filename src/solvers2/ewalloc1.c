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

#include "nlp.h"

 EvalWorkspace*
ewalloc1_ASL(ASL *asl)
{
	EvalWorkspace *rv;
	Invd1 *invd, **pinvd;
	arglist *al;
	char *s;
	const char **sa;
	func_info *fi;
	int i, nc, nf, nv, nv0, nvx;
	real *d, *lastx, *ra, *w, *z;
	size_t L, Le, Lx, Lxv;

	switch (asl->i.ASLtype) {
	  case ASL_read_f:
	  case ASL_read_fg:
		break;
	  default:
		fprintf(Stderr, "\nUnexpected ASLtype = %d in ewalloc1_AS(()\n", asl->i.ASLtype);
		fflush(Stderr);
		exit(1);
		}
	Lx = Lxv = 0;
	nf = asl->i.nfinv;
	Le = (sizeof(EvalWorkspace) + sizeof(real) - 1) & ~(sizeof(real)-1);
	L = Le + asl->i.derplen + asl->i.wlen + asl->i.numlen + nf*sizeof(arglist)
		+ asl->i.ra_max*sizeof(real) + asl->i.sa_max*sizeof(char*);
	nc = nclcon;
	if (asl->i.wlen)
		L += Lxv = (n_obj + nc + n_lcon)*sizeof(size_t);
	if (asl->i.vscale || (asl->i.vmap && asl->i.vmap[ASL_Sufkind_var]))
		L += Lx = x0len;
	if (asl->i.ew_bytes < L)
		asl->i.ew_bytes = L;
	rv = (EvalWorkspace*)M1alloc(L);
	ACQUIRE_MBLK_LOCK(&asl->i, MemLock);
	++asl->i.n_ew0;
	*asl->i.pewthread = rv;
	asl->i.pewthread = &rv->ewthread;
	FREE_MBLK_LOCK(&asl->i, MemLock);
	s = (char*)rv + Le;
	memset(rv, 0, sizeof(EvalWorkspace));
	rv->asl = asl;
	rv->x0kind = ASL_first_x;
	lastx = 0;
	if (Lx) {
		lastx = (real*)s;
		s += Lx;
		}
	memcpy(s, asl->i.numvals, asl->i.numlen);
	if (!asl->i.wlen) /* ASL_read_f */
		return rv;
	rv->w = w = (real*)(s += asl->i.numlen);
	if (!(rv->Lastx = lastx))
		rv->Lastx = w;
	nv = asl->i.n_var_;
	nv0 = asl->i.n_var0;
	nvx = nv0 + asl->i.nsufext[ASL_Sufkind_var];
	rv->dv = z = w + nvx;
	rv->dv1 = z + combc + como;
	s += asl->i.wlen;
	ra = (real*)s;
	s = (char*)(ra + asl->i.ra_max);
	rv->derps = d = (real*)s;
	rv->noxval = (size_t*)(s += asl->i.derplen);
	memset(rv->noxval, 0, Lxv);
	rv->ncxval = rv->noxval + n_obj;
	rv->nlxval = rv->ncxval + nc;
	sa = (const char**)(rv->nlxval + n_lcon);
	s = (char*)(sa + asl->i.sa_max);
	rv->wantderiv = want_derivs;
	if (asl->i.X0_) {
		memcpy(w, asl->i.X0_, nv0*sizeof(real));
		if (nv > nv0)
			memset(w + nv0, 0, (nv-nv0)*sizeof(real));
		}
	else
		memset(w, 0, nv*sizeof(real));
	if (nf) {
		rv->al = al = (arglist*)s;
		pinvd = (Invd1**)asl->i.invd;
		memset(al, 0, nf*sizeof(arglist));
		for(i = 0; i < nf; ++i, ++al) {
			invd = pinvd[i];
			al->n = al->nin = invd->n;
			if ((al->nr = invd->nr))
				al->ra = ra;
			al->at = invd->at;
			al->dig = invd->dig;
			if ((al->nsin = invd->n - invd->nr))
				al->sa = sa;
			fi = invd->fi;
			al->f = (function*)fi;
			al->funcinfo = fi->funcinfo;
			al->AE = asl->i.ae;
			}
		}
	return rv;
	}
