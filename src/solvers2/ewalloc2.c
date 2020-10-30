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

#define PSHVREAD
#include "asl.h"
#include "psinfo.h"
#ifdef NANDEBUG
#include "signal.h"
#include "fenv.h"
#include "fpu_control.h"
#endif

 static void
wk_init(real *w, int *z, real t)
{
	int i, n;

	n = *z;
	for(i = 1; i <= n; ++i)
		w[z[i]] = t;
	}

 EvalWorkspace*
ewalloc2_ASL(ASL *a)
{
	ASL_pfgh *asl;
	EvalWorkspace *rv;
	Varval *v, *ve;
	arglist *al;
	char *s;
	const char **sa;
	func_info *fi;
	int i, nc, nf, nlogc, nlv, no, nv, nv0, nvx;
	real *d, *lastx, *ra, *w;
	size_t L, La, Le, Lh, Lh0, Lo, Lu, Luw, Lx, Lxc, *nxc;
	tfinfo **ptfi, *tfi;
#ifdef NANDEBUG
	typedef union {real r; unsigned int x[2]; } U;
	U u0, *u, *ue;
#endif

	switch (a->i.ASLtype) {
	  /* case ASL_read_fgh: */
	  /* case ASL_read_pfg: */
	  case ASL_read_pfgh:
		break;
	  default:
		fprintf(Stderr, "\nUnexpected ASLtype = %d in ewalloc2_AS(()\n", a->i.ASLtype);
		fflush(Stderr);
		exit(1);
		}
	asl = (ASL_pfgh*)a;
	nc = nclcon;
	nlogc = n_lcon;
	no = n_obj;
	nf = asl->i.nfinv;
	La = nf * sizeof(arglist);
	Lo = (nc + no)*sizeof(real);
	Lu = asl->I.gscrx * sizeof(real);
	Lx = x0len;
	Le = (sizeof(EvalWorkspace) + sizeof(real) - 1) & ~(sizeof(real)-1);
	Lxc = ((2*(nc + no)+ nlogc)*sizeof(size_t) + sizeof(real) - 1) & ~(sizeof(real)-1);
	Lh =  (asl->I.nhop * sizeof(Hesoprod) + sizeof(real) - 1) & ~(sizeof(real)-1);
	Luw = (asl->I.uhlen + sizeof(real) - 1) & ~(sizeof(real)-1);
	if ((nlv = nlvc) < nlvo)
		nlv = nlvo;
	Lh0 = (nlv*sizeof(real) + asl->P.nran*sizeof(range*) + sizeof(real) - 1) & ~(sizeof(real)-1);
	L = La + Le + Lh + Lh0 + Lo + Lu + Luw + Lx + asl->i.derplen
		+ asl->i.wlen + asl->i.numlen + Lxc
		+ asl->i.ra_max*sizeof(real) + asl->i.sa_max*sizeof(char*);
	if (asl->i.ew_bytes < L)
		asl->i.ew_bytes = L;
	rv = (EvalWorkspace*)M1alloc(L);
	ACQUIRE_MBLK_LOCK(&asl->i, MemLock);
	++asl->i.n_ew0;
	*asl->i.pewthread = rv;
	asl->i.pewthread = &rv->ewthread;
	FREE_MBLK_LOCK(&asl->i, MemLock);
	s = (char*)rv + Le;
	memset(rv, 0, Le + Lxc);
	rv->asl = a;
	rv->nlv = nlv;
	rv->nxval = 1;
	rv->x0kind = ASL_first_x;
	nxc = (size_t*)s;
	rv->ncxval = nxc;
	rv->ncxgval = nxc += nc;
	rv->nlxval = nxc += nc;
	rv->noxval = nxc += nlogc;
	rv->noxgval = nxc += no;
	s += Lxc;
	rv->hop_free = 0;
	rv->hop_free0 = rv->hop_free = (Hesoprod*)s;
	rv->hop_free_end = rv->hop_free0 + asl->I.nhop;
	s += Lh;
	rv->uhw_free = 0;
	rv->uhw_free0 = (uHeswork*)s;
	rv->uhw_free_end = (uHeswork*)((char*)rv->uhw_free0 + asl->I.uhlen);
	s += Luw;
	rv->H0 = (real*)s;
	s += Lh0;
	lastx = (real*)s;
	s += Lx;
	rv->unopscr = (real*)s;
	s += Lu;
	rv->oyow0 = (real*)s;
	s += Lo;
	memcpy(s, asl->i.numvals, asl->i.numlen);
	rv->w = w = (real*)(s += asl->i.numlen);
#ifdef NANDEBUG
	u0.x[0] = 0x12345678;
	u0.x[1] = 0xfff7abcd;	/* signalling NaN */
	u = (U*)w;
	ue = (U*)(w + asl->P.rtodo);
	while(u < ue)
		*u++ = u0;
	__fpu_control = _FPU_IEEE - _FPU_EXTENDED + _FPU_DOUBLE - _FPU_MASK_IM;
	_FPU_SETCW(__fpu_control);
	signal(SIGFPE, fpecatch_ASL);
#endif
	nv = asl->i.n_var_;

	/* The following loop only matters under obscure conditions when setting */
	/* up Hessian compuations involving defined variables having linear parts */
	/* and used linearly in a constraint or objective. */
	/* With default AMPL settings, this situation does not arise. */
	v = (Varval*)w;
	for(ve = v + nv; v < ve; ++v)
		v->dO = v->aO = 0.;

	memset(w + asl->P.rtodo, 0, asl->P.zaplen);
	rv->Lastx = lastx;
	nv0 = asl->i.n_var0;
	nvx = nv0 + asl->i.nsufext[ASL_Sufkind_var];
	rv->dv = (Varval*)w + nvx;
	rv->dv1 = (Varval*)rv->dv + combc + como;
	s += asl->i.wlen;
	rv->derps = d = (real*)s;
	if ((i = asl->i.maxvar - asl->i.defvar0) > 0)
		memset(d + asl->i.defvar0, 0, i*sizeof(real));
	s += asl->i.derplen;
/*	rv->derpzap = d + nvx + combc + como + comc1 + como1; */
	rv->wantderiv = want_derivs;
	rv->ndhmax = asl->P.ndhmax;
	if (asl->i.X0_) {
		memcpy(w, asl->i.X0_, nv0*sizeof(real));
		if (nv > nv0)
			memset(w + nv0, 0, (nv-nv0)*sizeof(real));
		}
	else
		memset(w, 0, nv*sizeof(real));
	if (nf) {
		ra = (real*)s;
		sa = (const char**)(ra + asl->i.ra_max);
		rv->al = al = (arglist*)(sa + asl->i.sa_max);
		ptfi = (tfinfo**)asl->i.invd;
		memset(al, 0, nf*sizeof(arglist));
		for(i = 0; i < nf; ++i, ++al) {
			tfi = ptfi[i];
			al->n = al->nin = tfi->n;
			if ((al->nr = tfi->nr))
				al->ra = ra;
			al->at = tfi->at;
			al->dig = tfi->wd;
			if ((al->nsin = tfi->n - tfi->nr))
				al->sa = sa;
			fi = tfi->fi;
			al->f = (function*)fi;
			al->funcinfo = fi->funcinfo;
			al->AE = asl->i.ae;
			}
		}
	if (asl->P.wkinit0)
		wk_init(w, asl->P.wkinit0, 0.);
	if (asl->P.wkinit2)
		wk_init(w, asl->P.wkinit2, 2.);
	if (asl->P.wkinitm1)
		wk_init(w, asl->P.wkinitm1, -1.);
	return rv;
	}
