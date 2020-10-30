/*******************************************************************
Copyright (C) 2017, 2018, 2019, 2020 AMPL Optimization, Inc.; written by David M. Gay.

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
#include "opcode.hd"
#include "opno.hd"

#define oddnum(op,n) (((size_t)&op[n]) & 7)
#ifdef X64_bit_pointers
#define isodd(op,n) (((size_t)&op[n]) & (sizeof(void*)-1))
#define opalign(op,n,a) if (((size_t)&op[n]) & (sizeof(void*)-1)) *op++ = a;
#define alignarg(x) x
#define alignnum(a,b) (b)
#else
#define isodd(op,n) 0
#define opalign(op,n,a) /*nothing*/
#define alignarg(x)   /*nothing*/
#define alignnum(a,b) (a)
#endif

#ifdef Just_Linear
#define who "f_read"
#define fg_read_ASL f_read_ASL
#else
#define who "fg_read"
#ifdef __cplusplus
extern "C" {
 static real Missing_func(arglist*);
 static int compar(const void*, const void*, void*);
	}
#endif /* __cplusplus */
#endif /* Just_Linear */

#undef nzc
#undef r_ops

 typedef struct
Numhash { struct Numhash *nhnext, *nhthread; real *nhval; } Numhash;

 typedef union { real r; uint ui[2]; } UIR;

 typedef struct
Static {
	int _nv0;
	ASL *a;
	ASL_fg *asl;
	Invd1 **pinvd, **pinvde;
	Numhash **numht, *numhtf, *numhtfend, *numhtfirst, **numhtpth;
	int *opfirst, *opnext0, *_opnext, *_oplast;
	void **numhtf0;
	real *htvals, *htvals_end;
	int _wlast, k_seen;
	uint numht_mask, numht_n;
#ifndef Just_Linear
	derpblock curdb, **firstdb0, **firstdb1, **firstdbe;
	derp *first_d;
	int (*_holread) ANSI((EdRead*));
	int *_zc, *_zci, *c1s;
	int _amax1, _last_cex, _lastj, _max_var, _nv01, _nv011;
	int _nv0b, _nv0c, com1, k_afree;
	int kfirstdb, kinvd, negone, nvar0, nvinc, one;
	uint afirst, _alast, _nzc, ncond, nderp;
	uint *afree, *afree0, *afree1, *c1a, *oval;
#endif /* Just_Linear */
#if 1 /*TEMP DEBUG*/
	EdRead *R;
#endif
	} Static;

#define alast		S->_alast
#define amax1		S->_amax1
#define holread		S->_holread
#define if2list		S->_if2list
#define if2list_end	S->_if2list_end
#define iflist		S->_iflist
#define last_cex	S->_last_cex
#define lastj		S->_lastj
#define max_var		S->_max_var
#define nv0		S->_nv0
#define nv01		S->_nv01
#define nv011		S->_nv011
#define nv0b		S->_nv0b
#define nv0c		S->_nv0c
#define nzc		S->_nzc
#define oplast		S->_oplast
#define opnext		S->_opnext
#define varg2list	S->_varg2list
#define varg2list_end	S->_varg2list_end
#define varglist	S->_varglist
#define wlast		S->_wlast
#define zc		S->_zc
#define zci		S->_zci

#ifdef Just_Linear /*{{*/

 static void
sorry_nonlin(EdRead *R)
{
	fprintf(Stderr,
		"Sorry, %s can't handle nonlinearities.\n",
		progname ? progname : "");
	exit_ASL(R,ASL_readerr_nonlin);
	}
#else /*}{ Just_Linear */

 static uint
new_a(Static *S)
{
	if (S->afree > S->afree0)
		return *--S->afree;
	return ++alast;
	}

 static void
free_a(Static *S, uint a)
{
	ASL_fg *asl;
	size_t n;
	uint *t;

	if (S->afree >= S->afree1) {
		asl = S->asl;
		if (!S->k_afree) {
			S->afree = S->afree0 = (uint*)new_mblk(S->k_afree = 3);
			S->afree1 = (uint*)((char*)S->afree + (sizeof(void*) << S->k_afree));
			}
		else {
			n = S->afree1 - S->afree0;
			t = (uint*)new_mblk(++S->k_afree);
			memcpy(t, S->afree0, n*sizeof(uint));
			del_mblk(S->afree0);
			S->afree0 = t;
			S->afree = t + n;
			S->afree1= S->afree + n;
			}
		}
	*S->afree++ = a;
	}

#include "r_opn.hd"

 static void
db_reset(Static *S)
{
	S->curdb.de = S->curdb.d0;
	S->curdb.next = 0;
	S->curdb.nxt = 0;
	}
#endif /*}} Just_Linear */

 static int
numind(Static *S, real x)
{
	ASL *a;
	Numhash **nht, *p, **p0, **p00, *p1, **p2;
	real *r;
	size_t L;
	uint h, m;
	void **vp;

	h = (((UIR*)&x)->ui[0] ^ ((UIR*)&x)->ui[1]) % S->numht_mask;
	p0 = p00 = &S->numht[h];
	while((p = *p0)) {
		if (x == *p->nhval)
			goto ret;
		p0 = &p->nhnext;
		}
	if (S->numht_n > S->numht_mask) {
		m = S->numht_n;
		a = S->a;
		a->i.temp_rd_bytes += L = m*sizeof(real);
		r = (real*)Malloc(2*L);
		memcpy(r + m, S->htvals, L);
		free(S->htvals);
		S->htvals = r;
		S->htvals_end = r += 2*m;
		a->i.temp_rd_bytes += L = sizeof(Numhash*)*S->numht_n;
		nht = (Numhash**)Malloc(L *= 2);
		memset(nht, 0, L);
		m = S->numht_mask = 2*S->numht_mask + 1;
		*S->numhtpth = 0;
		for(p = S->numhtfirst; p; p = p->nhthread) {
			p->nhval = --r;
			h = (((UIR*)r)->ui[0] ^ ((UIR*)r)->ui[1]) % m;
			p2 = &nht[h];
			while((p1 = *p2))
				p2 = &p1->nhnext;
			*p2 = p;
			p->nhnext = 0;
			}
		h = (((UIR*)&x)->ui[0] ^ ((UIR*)&x)->ui[1]) % m;
		p00 = &nht[h];
		free(S->numht);
		S->numht = nht;
		}
	if (S->numhtf >= S->numhtfend) {
		vp = (void**)Malloc(L = 1024*sizeof(void*));
		S->a->i.temp_rd_bytes += L;
		*vp = S->numhtf0;
		S->numhtf0 = vp;
		S->numhtf = p1 = (Numhash*)(vp+1);
		S->numhtfend = p1 + (1023*sizeof(void*)/sizeof(Numhash));
		}
	*S->numhtpth = p = S->numhtf++;
	*(p->nhval = S->htvals_end - ++S->numht_n) = x;
	S->numhtpth = &p->nhthread;
	p->nhnext = *p00;
	*p00 = p;
 ret:
	return p->nhval - S->htvals_end;
	}

 static void
sorry_CLP(EdRead *R, const char *what)
{
	fprintf(Stderr,
		"Sorry, %s cannot handle %s.\n",
		progname ? progname : "", what);
	exit_ASL(R,ASL_readerr_CLP);
	}

 static Static *
ed_reset(Static *S, ASL *a)
{
	memset(S, 0, sizeof(Static));
	S->asl = (ASL_fg*)a;
	S->a = a;
	a->i.memLast = a->i.memNext = 0;
	return S;
	}

#ifndef Just_Linear

 static void
fscream(EdRead *R, const char *name, int nargs, const char *kind)
{
	scream(R, ASL_readerr_argerr,
		"line %ld: attempt to call %s with %d %sargs\n",
		R->Line, name, nargs, kind);
	}

 static void
note_firstdb(Static *S, derpblock *db)
{
	derpblock **fdb1;
	int n;

	if (S->firstdb1 >= S->firstdbe) {
		n = S->firstdbe - S->firstdb0;
		fdb1 = (derpblock**)new_mblk_ASL(S->a, ++S->kfirstdb);
		memcpy(fdb1, S->firstdb0, n*sizeof(derpblock*));
		Del_mblk_ASL(S->a, S->firstdb0);
		S->firstdb0 = fdb1;
		S->firstdb1 = fdb1 += n;
		S->firstdbe = fdb1 + n;
		}
	*S->firstdb1++ = db;
	}

/*DEBUG*/ derp *derpzork2;

 static void
new_derp(Static *S, uint a, uint b, int c)
{
	ASL *asl;
	derp *d;
	derpblock *db, *db1;

	++S->nderp;
	db = &S->curdb;
	if (db->d0 <= S->first_d) {
		asl = S->a;
		if (db->d0 < db->de) {
			db1 = (derpblock*)mem(sizeof(derpblock));
			*db1 = *db;
			if (!db1->nxt && !db1->next)
				note_firstdb(S, db1);
			db->next = db1;
			db->nxt = 1;
			}
		S->first_d = d = (derp*)M1alloc(8190*sizeof(derp));
		db->d0 = db->de = d + 8190;
		}
	d = --db->d0;
	/*DEBUG*/ if (d == derpzork2)
	/*DEBUG*/	printf("");
	d->a = a - 1;
	d->b = b - 1;
	d->c = c;
	}

#endif /* Just_Linear */

 static int opblock_gulp = 16384;

 static int*
new_opblock(Static *S, int need)
{
	ASL *asl;
	int n, *op, *opnew;

	asl = S->a;
	if (!S->opfirst && (op = S->opnext0)) {
		if (op + need <= oplast)
			return S->opfirst = op;
		opnext = op;
		}
	need += 4;
	if ((n = opblock_gulp) < need)
		n = need;
	opnew = (int*)M1alloc(n*sizeof(int));
	if ((op = opnext)) {
		if (isodd(op,1))
			*op++ = OPNEXTBLKalign;
		*op++ = OPNEXTBLK;
		*(int**)op = opnew;
		}
	if (!S->opfirst)
		S->opfirst = opnew;
	oplast = opnew + n - alignnum(2,4);
	return opnew;
	}

 static int*
nextop(Static *S, int n)
{
	int *op, *op1;

	op = opnext;
	op1 = op + n;
	if (op1 > oplast) {
		op = new_opblock(S, n);
		op1 = op + n;
		}
	opnext = op1;
	return op;
	}

 static int
eread(EdRead *R, uint *deriv, uint atop)
{
	ASL_fg *asl;
	Static *S;
	fint L1;
	int (*Xscanf)(EdRead*, const char*, ...);
	real r;
#ifndef Just_Linear
	Invd1 *invd, **pinvd;
	char *dig;
	derpblock dbsave, *db, *db1, *db2, **pdb;
	func_info *fi;
	int a0, a1, *at, i, ica, j, jca, k, k1, k2, k3, kd, numargs;
	int *op, *ope[2], *opg, *opg1, *opg2, *opg3, **pop, rv, symargs;
	plterm *p;
	real *b, t;
	size_t L;
	uint af, atop2, ia, ja, ka, kaf, *pia, *pja;
#else  /* Just_Linear */
	Not_Used(deriv);
#endif /* Just_Linear */

	S = (Static *)R->S;
	asl = S->asl;
	Xscanf = xscanf;
#ifndef Just_Linear
	ica = jca = 0;
#endif
	switch(edag_peek(R)) {
#ifdef Just_Linear
		case 'f':
		case 'h':
		case 'o':
		case 'v':
			sorry_nonlin(R);
#else
		case 'f':
			if (Xscanf(R, "%d %d", &i, &j) != 2)
				badline(R);
			fi = funcs[i];
			if (fi->nargs >= 0) {
				if (fi->nargs != j) {
 bad_nargs:
					fscream(R, fi->name, j, "");
					}
				}
			else if (-(1+fi->nargs) > j)
				goto bad_nargs;
			opg = (int*)new_mblk(htcl(j*(2*sizeof(int) + sizeof(uint) + 1)));
			at = opg + j;
			pja = (uint*)(at + j);
			dig = (char*)(pja + j);
			kd = symargs = numargs = 0;
			pia = 0;
			if (deriv)
				pia = &ia;
			ka = 0;
			af = S->afirst;
			k2 = k3 = -1;
			for(i = 0; i < j; ++i) {
				ia = ja = 0;
				S->oval = &ja;
				opg[i] = eread(R, pia, 0);
				pja[i] = ia;
				if (ia > af) {
					if (ia == atop || (ka > ia && ka != atop) || !ka) {
						if (ka)
							free_a(S, ka);
						ka = ia;
						k2 = i;
						k3 = numargs;
						}
					}
				if (ja == OPHOL || ja == OPIFSYM)
					at[i] = -(++symargs);
				else {
					dig[numargs] = 1;
					if (ia) {
						dig[numargs] = 0;
						++kd;
						}
					at[i] = numargs++;
					}
				}
			if (symargs && !(fi->ftype & 1))
				fscream(R, fi->name, symargs, "symbolic ");
			a1 = 0;
			if (kd) {
				if (!ka && !(ka = atop))
					ka = new_a(S);
				*deriv = ka;
				k = OPFUNCALL1;
				if (kd < numargs)
					a1 = (numargs + sizeof(int) - 1) / sizeof(int);
				}
			else {
				k = OPFUNCALL0;
				}
			op = nextop(S, j + 3);
			op[0] = k;
			op[1] = rv = wlast;
			i = rv + 1;
			if (kd)
				i += numargs;
			wlast = i;
			if (asl->i.ra_max < numargs)
				asl->i.ra_max = numargs;
			if (asl->i.sa_max < symargs)
				asl->i.sa_max = symargs;
			op[2] = i = asl->i.nfinv++;
			op += 3;
			pinvd = S->pinvd + i;
			if (pinvd >= S->pinvde) {
				if (!S->pinvde) {
					S->kinvd = k = 9;
					pinvd = S->pinvd = (Invd1**)new_mblk(k);
					S->pinvde = pinvd + (1 << k);
					}
				else {
					L = sizeof(void*) << S->kinvd;
					pinvd = (Invd1**)new_mblk(++S->kinvd);
					memcpy(pinvd, S->pinvd, L);
					del_mblk(S->pinvd);
					S->pinvd = pinvd;
					S->pinvde = pinvd + (1 << S->kinvd);
					pinvd += i;
					}
				}
			k = numargs + symargs;
			*pinvd = invd = (Invd1*)mem(sizeof(Invd1) + k*sizeof(int) + a1);
			invd->fi = fi;
			invd->nr = numargs;
			invd->at = 0;
			invd->dig = 0;
			if ((invd->n = k)) {
				memcpy(invd->at = (int*)(invd+1), at, k*sizeof(int));
				if (a1)
					memcpy(invd->dig = (char*)(invd->at + k), dig, numargs);
				}
			if (numargs) {
				for(i = 0; i < j; ++i) {
					if (at[i] >= 0)
						*op++ = opg[i];
					}
				if (kd) {
					k1 = rv + 1;
					if (k2 >= 0) {
						ia = pja[k2];
						pja[k2] = 0;
						new_derp(S, ia, ka, k1+k3);
						}
					while(i > 0) {
						if ((k2 = at[--i]) >= 0 && (ia = pja[i])) {
							new_derp(S, ia, ka, k1+k2);
							if (ia > ka)
								free_a(S, ia);
							}
						}
					}
				}
			if (symargs) {
				for(i = 0; i < j; ++i)
					if (at[i] < 0)
						*op++ = opg[i];
				}
			opnext = op;
			del_mblk(opg);
			return rv;

		case 'h':
			if ((pia = S->oval)) {
				*pia = OPHOL;
				S->oval = 0;
				}
			return holread(R);
#endif /* Just_Linear */
		case 's':
			if (Xscanf(R, "%hd", &L1) != 1)
				badline(R);
			r = L1;
			goto have_r;

		case 'l':
			if (Xscanf(R, "%ld", &L1) != 1)
				badline(R);
			r = L1;
			goto have_r;

		case 'n':
			if (Xscanf(R, "%lf", &r) != 1)
				badline(R);
 have_r:
			return numind(S, r);
#ifndef Just_Linear

		case 'o':
			break;

		case 'v':
			if (Xscanf(R,"%d",&k) != 1 || k < 0)
				badline(R);
			if (k >= S->nvar0)
				k += S->nvinc;
			if (k >= max_var)
				badline(R);
			if (k >= nv01) {
				if (deriv)
					*deriv = S->c1a[k];
				}
			else if (deriv) {
				*deriv = k + 1;
				if (!zc[k]) {
					zci[nzc++] = k;
					zc[k] = nzc;
					}
				}
			return k;

#endif /* Just_Linear */
		default:
			badline(R);
		}

#ifndef Just_Linear

	if (Xscanf(R, asl->i.opfmt, &k) != 1 || k < 0 || k >= N_OPS)
		badline(R);
	if ((pia = S->oval)) {
		*pia = k;
		S->oval = 0;
		}
	switch(k) {
	  case FLOOR:
		k = nFLOOR;
		goto ceil1;
	  case CEIL:
		k = nCEIL;
 ceil1:
		i = eread(R, 0, 0);
		rv = wlast++;
		goto uret;
	  case ABS:
		i = eread(R, deriv, 0);
		rv = wlast++;
		if (deriv && (ia = *deriv)) {
			if (ia > S->afirst)
				free_a(S, ia);
			k = nOPABS1;
			if (!(ja = atop))
				ja = new_a(S);
			new_derp(S, ia, *deriv = ja, wlast++);
			}
		else
			k = nOPABS0;
		goto uret;
	  case OPUMINUS:
		k = nOPUMINUS;
		i = eread(R, deriv, 0);
		rv = wlast++;
		if (deriv && (ia = *deriv)) {
			if (ia > S->afirst)
				free_a(S, ia);
			if (!S->negone)
				S->negone = numind(S, -1.);
			if (!(ja = atop))
				ja = new_a(S);
			new_derp(S, ia, *deriv = ja, S->negone);
			}
		goto uret;
	  case OP_tanh:
		k = nOP_tanh0;
		k1 = nOP_tanh1;
		goto unop;
	  case OP_tan:
		k = nOP_tan0;
		k1 = nOP_tan1;
		goto unop;
	  case OP_sqrt:
		k = nOP_sqrt0;
		k1 = nOP_sqrt1;
		goto unop;
	  case OP_sinh:
		k = nOP_sinh0;
		k1 = nOP_sinh1;
		goto unop;
	  case OP_sin:
		k = nOP_sin0;
		k1 = nOP_sin1;
		goto unop;
	  case OP_log10:
		k = nOP_log100;
		k1 = nOP_log101;
		goto unop;
	  case OP_log:
		k = nOP_log0;
		k1 = nOP_log1;
		goto unop;
	  case OP_exp:
		k = nOP_exp;
		i = eread(R, deriv, 0);
		rv = wlast++;
		if (deriv && (ia = *deriv)) {
			if (ia > S->afirst)
				free_a(S, ia);
			if (!(ja = atop))
				ja = new_a(S);
			new_derp(S, ia, *deriv = ja, rv);
			}
		goto uret;
	  case OP_cosh:
		k = nOP_cosh0;
		k1 = nOP_cosh1;
		goto unop;
	  case OP_cos:
		k = nOP_cos0;
		k1 = nOP_cos1;
		goto unop;
	  case OP_atanh:
		k = nOP_atanh0;
		k1 = nOP_atanh1;
		goto unop;
	  case OP_atan:
		k = nOP_atan0;
		k1 = nOP_atan1;
		goto unop;
	  case OP_asinh:
		k = nOP_asinh0;
		k1 = nOP_asinh1;
		goto unop;
	  case OP_asin:
		k = nOP_asin0;
		k1 = nOP_asin1;
		goto unop;
	  case OP_acosh:
		k = nOP_acosh0;
		k1 = nOP_acosh1;
		goto unop;
	  case OP_acos:
		k = nOP_acos0;
		k1 = nOP_acos1;
 unop:
		i = eread(R, deriv, 0);
		rv = wlast++;
		if (deriv && (ia = *deriv)) {
			if (ia > S->afirst)
				free_a(S, ia);
			if (!(ja = atop))
				ja = new_a(S);
			new_derp(S, ia, *deriv = ja, wlast++);
			k = k1;
			}
 uret:
		op = nextop(S, 3);
		op[0] = k;
		op[1] = rv;
		op[2] = i;
		return rv;
	  case OPNOT:
		if (deriv) {
			*deriv = 0;
			deriv = 0;
			}
		k = nOPNOT;
		i = eread(R, 0, 0);
		rv = wlast++;
		goto uret;
	  case OPPLUS:
		k = nOPPLUS;
		jca = S->one;
		atop2 = atop;
 bret2:
		rv = wlast++;
		if (deriv) {
			ia = ja = 0;
			*deriv = 0;
			if (atop)
				kaf = 0;
			else
				atop = kaf = new_a(S);
			i = eread(R, &ia, atop);
			j = eread(R, &ja, ia == atop2 ? 0 : atop2);
			if (ia | ja) {
				af = S->afirst;
				if (ja > af && ja != atop)
					free_a(S, ja);
				if (ia > af && ia != atop)
					free_a(S, ia);
				if (!(ka = atop))
					ka = new_a(S);
				*deriv = ka;
				if (ka == ja) {
					if (!atop2)
						new_derp(S, ja, ka, jca);
					ja = 0;
					}
				if (ia && ia != ka)
					new_derp(S, ia, ka, S->one);
				if (ja)
					new_derp(S, ja, ka, jca);
				}
			else if (kaf)
				free_a(S, atop);
			}
		else {
			i = eread(R, 0, 0);
			j = eread(R, 0, 0);
			}
 bret:
		op = nextop(S, 4);
		op[0] = k;
		op[2] = i;
		op[3] = j;
		return op[1] = rv;
	  case OPMINUS:
		if (deriv && !(jca = S->negone))
			S->negone = jca = numind(S, -1.);
		k = nOPMINUS;
		atop2 = 0;
		goto bret2;
	  case OPMULT:
		k = nOPMULT;
		rv = wlast++;
		if (deriv) {
			ia = ja = 0;
			*deriv = 0;
			i = eread(R, &ia, 0);
			j = eread(R, &ja, 0);
			if (ia | ja) {
				ica = j;
				jca = i;
 binderp:
				af = S->afirst;
				if (ja > af)
					free_a(S, ja);
				if (ia > af)
					free_a(S, ia);
				if (!(ka = atop))
					ka = new_a(S);
				*deriv = ka;
				if (ka == ja) {
					new_derp(S, ja, ka, jca);
					ja = 0;
					}
				if (ia)
					new_derp(S, ia, ka, ica);
				if (ja)
					new_derp(S, ja, ka, jca);
				}
			}
		else {
			i = eread(R, 0, 0);
			j = eread(R, 0, 0);
			}
		goto bret;
	  case OPDIV:
		k = nOPDIV0;
		if (deriv) {
			ia = ja = 0;
			*deriv = 0;
			i = eread(R, &ia, 0);
			j = eread(R, &ja, 0);
			rv = wlast++;
			if (ia | ja) {
				if (ia) {
					ica = rv + 1;
					if (ja) {
						k = nOPDIV3;
						wlast = rv + 3;
						jca = rv + 2;
						}
					else {
						k = nOPDIV1;
						wlast = rv + 2;
						}
					}
				else {
					k = nOPDIV2;
					wlast = rv + 2;
					jca = rv + 1;
					}
				goto binderp;
				}
			}
		else {
			i = eread(R, 0, 0);
			j = eread(R, 0, 0);
			rv = wlast++;
			}
		goto bret;
	  case OPREM:
		k = nOPREM0;
		ia = ja = 0;
		if (deriv) {
			*deriv = 0;
			i = eread(R, &ia, 0);
			j = eread(R, &ja, 0);
			}
		else {
			i = eread(R, 0, 0);
			j = eread(R, 0, 0);
			}
		rv = wlast++;
		if (ia | ja) {
			af = S->afirst;
			if (ja > af)
				free_a(S, ja);
			if (ia > af)
				free_a(S, ia);
			if (!(ka = atop))
				ka = new_a(S);
			*deriv = ka;
			if (ja) {
				wlast = rv + 2;
				k = nOPREM1;
				new_derp(S, ja, ka, rv+1);
				}
			if (ia && ia != ka)
				new_derp(S, ia, ka, S->one);
			}
		goto bret;
	  case OPPOW:
		k = nOPPOW0;
		ia = ja = 0;
		if (deriv) {
			*deriv = 0;
			i = eread(R, &ia, 0);
			j = eread(R, &ja, 0);
			}
		else {
			i = eread(R, 0, 0);
			j = eread(R, 0, 0);
			}
		rv = wlast++;
		if (ia | ja) {
			af = S->afirst;
			if (ja > af)
				free_a(S, ja);
			if (ia > af)
				free_a(S, ia);
			if (!(ka = atop))
				ka = new_a(S);
			*deriv = ka;
			if (ia) {
				if (ja) {
					wlast = rv + 3;
					k = nOPPOW3;
					new_derp(S, ia, ka, rv+1);
					new_derp(S, ja, ka, rv+2);
					}
				else {
					wlast = rv + 2;
					k = nOPPOW1;
					new_derp(S, ia, ka, rv+1);
					if (j < 0) {
						if ((t = S->htvals_end[j]) == 2.) {
							k = OP2POW1;
 pow2:
							op = nextop(S, 3);
							op[0] = k;
							op[2] = i;
							return op[1] = rv;
							}
						if (t == floor(t))
							k = nOPPOW4;
						}
					}
				}
			else {
				wlast = rv + 2;
				k = nOPPOW2;
				new_derp(S, ja, ka, rv+1);
				if (i < 0) {
					t = S->htvals_end[i];
					if (t > 0.) {
						op = nextop(S, 7);
						if (isodd(op,1)) {
							op[0] = OPCPOWalign;
							op[1] = j;
							b = (real*)&op[2];
							}
						else {
							op[0] = nOPCPOW;
							b = (real*)&op[1];
							op[5] = j;
							}
						b[0] = t;
						b[1] = log(t);
						opnext = op + 7;
						return op[6] = rv;
						}
					}
				}
			}
		else if (j < 0 && S->htvals_end[j] == 2.) {
			k = OP2POW0;
			goto pow2;
			}
		goto bret;
	  case OPLESS:
		k = nOPLESS0;
		if (deriv) {
			ia = ja = 0;
			*deriv = 0;
			dbsave = S->curdb;
			i = eread(R, &ia, 0);
			j = eread(R, &ja, 0);
			rv = wlast++;
			if (ia | ja) {
				++S->ncond;
				af = S->afirst;
				if (ja > af)
					free_a(S, ja);
				if (ia > af)
					free_a(S, ia);
				if (!(ka = atop))
					ka = new_a(S);
				op = nextop(S, alignnum(6,8));
				if (isodd(op,1)) {
					op[0] = OPLESSalign;
					op[1] = rv;
					pdb = (derpblock**)&op[2];
					op = (int*)&pdb[2];
					}
				else {
					op[0] = nOPLESS1;
					pdb = (derpblock**)&op[1];
					op = (int*)&pdb[2];
					*op++ = rv;
					}
				op[0] = i;
				op[1] = j;
				*deriv = ka;
				if (ia)
					new_derp(S, ia, ka, S->one);
				if (ja) {
					if (!S->negone)
						S->negone = numind(S, -1.);
					new_derp(S, ja, ka, S->negone);
					}
				db = pdb[0] = (derpblock*)mem(sizeof(derpblock));
				*db = S->curdb;
				S->curdb.next = 0;
				S->curdb.nxt = rv + 1;
				db = pdb[1] = (derpblock*)mem(sizeof(derpblock));
				*db = dbsave;
				if (!dbsave.nxt && !dbsave.next)
					note_firstdb(S, db);
				wlast = rv + 2;
				opnext = op + 2;
				return rv;
				}
			}
		else {
			i = eread(R, 0, 0);
			j = eread(R, 0, 0);
			rv = wlast++;
			}
		goto bret;
	  case OPOR:
		k = nOPOR;
		alignarg(k1 = OPORalign;)
		goto andor;
	  case OPAND:
		k = nOPAND;
		alignarg(k1 = OPANDalign;)
 andor:
		i = eread(R, 0, 0);
		op = nextop(S, alignnum(4,6));
		opalign(op, 3, k1)
		op[0] = k;
		op[1] = rv = wlast++;
		op[2] = i;
		pop = (int**)&op[3];
		opnext = (int*)&pop[1];
		j = eread(R, 0, 0);
		if (rv != j) {
			opg = nextop(S, 3);
			opg[0] = OPCOPY;
			opg[1] = rv;
			opg[2] = j;
			}
		pop[0] = opnext;
		return rv;
	  case LT:
		k = nOPLT;
 bnoderiv:
		if (deriv) {
			*deriv = 0;
			deriv = 0;
			}
		i = eread(R, 0, 0);
		j = eread(R, 0, 0);
		rv = wlast++;
		goto bret;
	  case LE:
		k = nOPLE;
		goto bnoderiv;
	  case EQ:
		k = nOPEQ;
		goto bnoderiv;
	  case GE:
		k = nOPGE;
		goto bnoderiv;
	  case GT:
		k = nOPGT;
		goto bnoderiv;
	  case NE:
		k = nOPNE;
		goto bnoderiv;
	  case OP_atan2:
		k = nOP_atan20;
		if (deriv) {
			*deriv = 0;
			ia = ja = 0;
			i = eread(R, &ia, 0);
			j = eread(R, &ja, 0);
			rv = wlast++;
			if (ia | ja) {
				wlast = rv + 2;
				ica = jca = rv + 1;
				if (ia) {
					k = nOP_atan21;
					if (ja) {
						k = nOP_atan23;
						jca = rv + 2;
						wlast = rv + 3;
						}
					}
				else
					k = nOP_atan22;
				goto binderp;
				}
			}
		else {
			i = eread(R, 0, 0);
			j = eread(R, 0, 0);
			rv = wlast++;
			}
		goto bret;
	  case OPintDIV:
		k = nOPintDIV;
		goto bnoderiv;
	  case OPprecision:
		k = nOPprecision;
		goto bnoderiv;
	  case OPround:
		k = nOPround;
		goto bnoderiv;
	  case OPtrunc:
		k = nOPtrunc;
		goto bnoderiv;
	  case OPATLEAST:
		k = nOPATLEAST;
		goto bnoderiv;
	  case OPATMOST:
		k = nOPATMOST;
		goto bnoderiv;
	  case OPEXACTLY:
		k = nOPEXACTLY;
		goto bnoderiv;
	  case OPNOTATLEAST:
		k = nOPNOTATLEAST;
		goto bnoderiv;
	  case OPNOTATMOST:
		k = nOPNOTATMOST;
		goto bnoderiv;
	  case OPNOTEXACTLY:
		k = nOPNOTEXACTLY;
		goto bnoderiv;
	  case OP_IFF:
		k = nOP_IFF;
		goto bnoderiv;
	  case MINLIST:
		k = nOPMINLIST0;
		k1 = nOPMINLIST1;
		alignarg(k2 = OPMINLIST1align;)
		goto maxminlist;
	  case MAXLIST:
		k = nOPMAXLIST0;
		k1 = nOPMAXLIST1;
		alignarg(k2 = OPMAXLIST1align;)
 maxminlist:
		i = -1;
		Xscanf(R, "%d", &i);
		if (i <= 0)
			badline(R);
		j = htcl(i*(sizeof(derpblock*) + sizeof(int)));
		pdb = (derpblock**)new_mblk(j);
		opg = (int*)(pdb + i);
		if (deriv) {
			*deriv = 0;
			dbsave = S->curdb;
			S->curdb.de = dbsave.d0;
			db = db2 = 0;
			ka = 0;
			af = S->afirst;
			for(j = 0; j < i; ++j) {
				ia = 0;
				opg[j] = eread(R, &ia, 0);
				if (ia) {
					if (ia > af)
						free_a(S, ia);
					if (!db) {
						if (dbsave.d0) {
							db2 = db = (derpblock*)mem(sizeof(derpblock));
							*db = dbsave;
							if (!dbsave.nxt && !dbsave.next)
								note_firstdb(S, db);
							for(k = 0; k < j; ++k)
								pdb[k] = db;
							S->curdb.nxt = 1;
							S->curdb.next = db;
							}
						else
							db = &dbsave;
						if (!(ka = atop))
							ka = new_a(S);
						*deriv = ka;
						free_a(S, ka);
						}
					if (ia != ka)
						new_derp(S, ia, ka, S->one);
					pdb[j] = db1 = (derpblock*)mem(sizeof(derpblock));
					*db1 = S->curdb;
					S->curdb.de = S->curdb.d0;
					S->curdb.nxt = (S->curdb.next = db2) ? 1 : 0;
					}
				else
					pdb[j] = db2;
				}
			if (db) {
				++S->ncond;
				j = alignnum(3,4) + i*(1 + sizeof(derpblock*)/sizeof(int));
				op = nextop(S, j);
				opalign(op,i+3,k2)
				op[0] = k1;
				op[1] = rv = wlast;
				op[2] = i;
				memcpy(op += 3, opg, i*sizeof(int));
				memcpy(op += i, pdb, i*sizeof(derpblock*));
				opnext = (int*)&((derpblock**)op)[i];
				wlast = rv + 2;
				del_mblk(pdb);
				S->curdb.de = S->curdb.d0;
				S->curdb.next = 0;
				S->curdb.nxt = rv + 1;
				if (ka != atop) {
					pia = S->afree0;
					pja = --S->afree;
					for(;; --pja) {
						if (pja < pia)
							scream(R, 1, "***Bad afree in minmaxlist\n");
						if (*pja == ka) {
							*pja = *S->afree;
							break;
							}
						}
					}
				return rv;
				}
			}
		else
			for(j = 0; j < i; ++j)
				opg[j] = eread(R, 0, 0);
		op = nextop(S, i + 3);
		op[0] = k;
		op[1] = rv = wlast++;
		op[2] = i;
		memcpy(op + 3, opg, i*sizeof(int));
		del_mblk(pdb);
		return rv;

	  case OPPLTERM: /* piece-wise linear */
		i = -1;
		Xscanf(R, "%d", &i);
		if (i <= 1)
			badline(R);
		if (deriv) {
			k = nOPPLTERM1;
			alignarg(k1 = OPPLTERM1align;)
			}
		else {
			k = nOPPLTERM0;
			alignarg(k1 = OPPLTERM0align;)
			}
		plterms++;
		j = 2*i - 1;
		k2 = (sizeof(plterm) + (j-1)*sizeof(real))/sizeof(int);
		op = nextop(S, k2 + alignnum(3,4));
		opalign(op,3,k1)
		op[0] = k;
		op[1] = rv = wlast;
		wlast = deriv ? rv + 2 : rv + 1;
		p = (plterm *)&op[3];
		p->n = i;
		b = p->bs;
		do {
			switch(edag_peek(R)) {
				case 's':
					if (Xscanf(R,"%hd",&L1) != 1)
						badline(R);
					r = L1;
					break;
				case 'l':
					if (Xscanf(R,"%ld",&L1) != 1)
						badline(R);
					r = L1;
					break;
				case 'n':
					if (Xscanf(R,"%lf",&r) == 1)
						break;
				default:
					badline(R);
				}
			*b++ = r;
			}
			while(--j > 0);
		opnext = (int*)b;
		if (b[-2] <= 0.)
			p->z = 2*i - 2;
		else {
			b = p->bs + 1;
			while(*b <= 0.)
				b += 2;
			p->z = (b - p->bs) - 1;
			}
		if (deriv) {
			ia = 0;
			op[2] = eread(R, deriv, 0);
			if ((ia = *deriv)) {
				if (!(ka = atop) && (ka = ia) <= S->afirst)
					ka = new_a(S);
				*deriv = ka;
				new_derp(S, ia, ka, rv+1);
				}
			}
		else
			op[2] = eread(R, 0, 0);
		return rv;
	  case OPIFnl:
		k = eread(R, 0, 0);
		op = nextop(S, alignnum(13,21));
		rv = a0 = wlast;
		ia = ja = ka = 0;
		pia = pja = 0;
		db = db1 = 0;
		dbsave = S->curdb;
		if (deriv) {
			pia = &ia;
			pja = &ja;
			}
		opg1 = opnext;
		i = eread(R, pia, 0);
		ope[0] = opnext;
		if (rv < i)
			rv = i;
		opg = nextop(S, alignnum(5,7));
		af = S->afirst;
		if (ia) {
			if  (ia > af)
				free_a(S, ia);
			if (!(ka = atop))
				ka = new_a(S);
			*deriv = ka;
			db = (derpblock*)mem(sizeof(derpblock));
			*db = dbsave;
			if (!dbsave.nxt)
				note_firstdb(S, db);
			if (ia != ka)
				new_derp(S, ia, ka, S->one);
			db1 = (derpblock*)mem(sizeof(derpblock));
			*db1 = S->curdb;
			S->curdb.de = S->curdb.d0;
			if ((S->curdb.next = db))
				S->curdb.nxt = 1;
			}
		a1 = wlast;
		wlast = a0;
		opg2 = opnext;
		j = eread(R, pja, 0);
		ope[1] = opnext;
		if (ja) {
			if (ja > af)
				free_a(S, ja);
			if (!ka) {
				if (!(ka = atop))
					ka = new_a(S);
				*deriv = ka;
				}
			if (ja != ka)
				new_derp(S, ja, ka, S->one);
			}
		if (rv < j)
			rv = j;
		if (wlast < a1)
			wlast = a1;
		if (rv == wlast)
			wlast = rv + 1;
		if (rv != j) {
			opg3 = nextop(S, 3);
			opg3[0] = OPCOPY;
			opg3[1] = rv;
			opg3[2] = j;
			}
		if (rv != i) {
			opg[0] = OPCOPY;
			opg[1] = rv;
			opg[2] = i;
			opg += 3;
			}
		opalign(opg, 1, OPGOTOalign;)
		opg[0] = OPGOTO;
		*(int**)&opg[1] = opnext;
		if (ia | ja) {
			++S->ncond;
			opalign(op, 6, OPIFnl1align)
			op[0] = nOPIFnl1;
			pop = (int**)&op[6];
			pdb = (derpblock**)&pop[5];
			pdb[0] = db1;
			db1 = db;
			if (ja) {
				db1 = (derpblock*)mem(sizeof(derpblock));
				*db1 = S->curdb;
				}
			pdb[1] = db1;
			S->curdb.de = S->curdb.d0;
			S->curdb.next = 0;
			S->curdb.nxt = op[5] = wlast++;
			}
		else {
			opalign(op, 5, OPIFnl0align)
			op[0] = nOPIFnl0;
			pop = (int**)&op[5];
			}
		op[1] = rv;	/* for fg_write */
		op[2] = k;
		op[3] = i;	/* for fg_write */
		op[4] = j;	/* for fg_write */
		pop[0] = opg1;
		pop[1] = opg2;
		pop[2] = ope[0];	/* for fg_write */
		pop[3] = ope[1];	/* for fg_write */
		pop[4] = opnext;	/* for fg_write */
		return rv;

	  case OPIMPELSE:
		k1 = OPCOPY;
		alignarg(k2 = OPIMPELSE_align;)
		k3 = nOPIMPELSE;
		goto moreifsym;

	  case OPIFSYM:
		k1 = OPCOPYSYM;
		alignarg(k2 = OPIFnl0align;)
		k3 = nOPIFnl0;
 moreifsym:
		k = eread(R, 0, 0);
		op = nextop(S, alignnum(10,16));
		opalign(op, 5, k2)
		op[0] = k3;
		op[2] = k;
		pop = (int**)&op[5];
		rv = a0 = wlast;
		pop[0] = opnext;
		op[3] = i = eread(R, 0, 0);
		pop[2] = opnext;
		if (rv < i)
			rv = i;
		opg = nextop(S, alignnum(5,7));
		pop[1] = opnext;
		a1 = wlast;
		wlast = a0;
		op[4] = j = eread(R, 0, 0);
		pop[3] = opnext;
		if (rv < j)
			rv = j;
		op[1] = rv;
		if (wlast < a1)
			wlast = a1;
		if (rv == wlast)
			wlast = rv + 1;
		if (rv != j) {
			opg2 = nextop(S, 3);
			opg2[0] = k1;
			opg2[1] = rv;
			opg2[2] = j;
			}
		if (rv != i) {
			opg[0] = k1;
			opg[1] = rv;
			opg[2] = i;
			opg += 3;
			}
		opalign(opg, 1, OPGOTOalign)
		opg[0] = OPGOTO;
		*(int**)&opg[1] = pop[4] = opnext;
		return rv;

	  case OPCOUNT:
		k = nOPCOUNT;
		k1 = 1;
		goto more_orlist;
	  case OPNUMBEROF:
		k = nOPNUMBEROF;
		k1 = 1;
		goto more_orlist;
	  case OPNUMBEROFs:
		k = nOPNUMBEROFs;
		k1 = 1;
		goto more_orlist;
	  case OPALLDIFF:
		k = nOPALLDIFF;
		k1 = 1;
		goto more_orlist;
	  case OPSOMESAME:
		k = nOPSOMESAME;
		k1 = 1;
		goto more_orlist;
	  case ANDLIST:
		k = nOPANDLIST;
		k1 = 3;
		goto more_orlist;
	  case ORLIST:
		k = nOPORLIST;
		k1 = 3;
 more_orlist:
		j = 0;
		Xscanf(R, "%d", &j);
		if (j < k1)
			badline(R);
		opg = (int*)new_mblk(htcl(j*sizeof(int)));
		for(i = 0; i < j; ++i)
			opg[i] = eread(R, 0, 0);
 finish_orlist:
		op = nextop(S, j + 3);
		op[0] = k;
		op[1] = rv = wlast++;
		op[2] = j;
		memcpy(op+3, opg, j*sizeof(int));
		del_mblk(opg);
		return rv;

	  case OPSUMLIST:
		k = nOPSUMLIST;
		if (!deriv) {
			k1 = 3;
			goto more_orlist;
			}
		j = 0;
		Xscanf(R, "%d", &j);
		if (j < 3)
			badline(R);
		opg = (int*)new_mblk(htcl(j*sizeof(int)));
		if (atop)
			kaf = 0;
		else
			atop = kaf = new_a(S);
		atop2 = atop;
		af = S->afirst;
		ica = S->one;
		for(i = k2 = 0; i < j; ++i) {
			ia = 0;
			opg[i] = eread(R, &ia, atop2);
			if (ia) {
				++k2;
				if (ia == atop2)
					atop2 = 0;
				else {
					new_derp(S, ia, atop, ica);
					if (ia > af)
						free_a(S, ia);
					}
				}
			}
		if (k2)
			*deriv = atop;
		else if (kaf)
			free_a(S, atop);
		goto finish_orlist;
		}
#endif /* Just_Linear */
	badline(R);
	return 0;
	}

 static int *
eread1(EdRead *R, derpblock **pdb, uint *pia, uint atop)
{
	Static *S;
	int *op, rv;
	uint ia;
#ifndef Just_Linear
	ASL_fg *asl;
	derpblock *db;
	uint ka;
#endif

	S = (Static *)R->S;
	S->opfirst = 0;
	S->opnext0 = opnext;
	opnext = oplast;
	if (pdb) {
		*pdb = 0;
		if (!pia)
			pia = &ia;
		}
	if (pia)
		*pia = 0;
	rv = eread(R, pia, atop);
	op = nextop(S, 2);
	op[0] = nOPRET;
	op[1] = rv;
#ifndef Just_Linear
	if (pia && (ia = *pia)) {
		if (S->curdb.d0 == S->curdb.de && !S->curdb.nxt) {
			ka = new_a(S);
			new_derp(S, ia, ka, S->one);
			}
		if (pdb) {
			asl = S->asl;
			*pdb = db = (derpblock*)mem(sizeof(derpblock));
			*db = S->curdb;
			db_reset(S);
			}
		}
#endif
	return S->opfirst;
	}

#ifndef Just_Linear /*{*/

 static int
compar(const void *a, const void *b, void *v)
{
	Not_Used(v);
	return *(int*)a - *(int*)b;
	}

 static void
derpcopy(Static *S, int kc, int k, int nz, int *ci, derpblock *db)
{
	ASL_fg *asl;
	cexp *ce, *ck, *cx0;
	derp *d0, *de, *first_d;
	derpblock cb, **db0, **db1, *dbo, *dbx;
	int i, na, nn, nv;

	if (!db->nxt)
		note_firstdb(S, db);
	nv = nv0;
	asl = S->asl;
	cx0 = cexps;
	cb = S->curdb;
	first_d = S->first_d;
	if (cb.d0 <= first_d) {
		first_d = (derp*)M1alloc(8190*sizeof(derp));
		cb.d0 = first_d + 8190;
		}
	dbo = (derpblock*)mem(sizeof(derpblock));
	dbo->nxt = 0;
	dbo->next = 0;
	db0 = S->firstdb0;
	db1 = S->firstdb1;
	S->firstdb1 = S->firstdb0;
	while (db1 > db0) {
		db = *--db1;
		db->next = dbo;
		db->nxt = 1;
		}
	dbo->d0 = dbo->de = first_d;
	while(nz > k) {
		i = ci[--nz];
		ce = &cx0[i - nv];
		for(db = ce->dbf; db; db = db->next) {
			d0 = db->d0;
			de = db->de;
			while((nn = de - d0) > 0) {
				if (cb.d0 <= first_d) {
					first_d = (derp*)M1alloc(8190*sizeof(derp));
					cb.d0 = first_d + 8190;
					dbo->nxt = 1;
					dbo->next = dbx = (derpblock*)mem(sizeof(derpblock));
					dbx->nxt = 0;
					dbx->next = 0;
					dbo = dbx;
					dbo->d0 = first_d;
					}
				na = cb.d0 - first_d;
				if (na > nn)
					na = nn;
				memcpy(first_d, d0, na*sizeof(derp));
				d0 += na;
				dbo->de = first_d += na;
				}
			}
		}
	if (!kc)
		goto done;
	if (cb.d0 <= first_d) {
		first_d = (derp*)M1alloc(8190*sizeof(derp));
		cb.d0 = first_d + 8190;
		}
	ck = &cx0[kc - nv];
	ck->dbf = dbo = (derpblock*)mem(sizeof(derpblock));
	dbo->nxt = 0;
	dbo->next = 0;
	dbo->d0 = first_d;
	i = wlast;
	while(k > 0) {
		if (cb.d0 <= first_d) {
			dbo->de = first_d;
			first_d = (derp*)M1alloc(8190*sizeof(derp));
			cb.d0 = first_d + 8190;
			dbo->nxt = 1;
			dbo->next = dbx = (derpblock*)mem(sizeof(derpblock));
			dbx->nxt = 0;
			dbx->next = 0;
			dbo = dbx;
			dbo->d0 = first_d;
			}
		first_d->a = ci[--k];
		first_d->b = kc;
		first_d->c = i++;
		++first_d;
		}
	dbo->de = first_d;
	wlast = i;
 done:
	S->first_d = first_d;
	S->curdb.d0 = cb.d0;
	db_reset(S);
	}
#endif /*} Just_Linear */

 static void
co_read(EdRead *R, cde *d, int k, int wd, int *cexp1st)
{
#ifndef Just_Linear
	ASL_fg *asl;
	Static *S;
	cexp *cx0, *di;
	cexp1 *c1x;
	derpblock *db;
	int *c, *ci, *dv, i, j, j1, j2, ndv0, nv, nz, nz1, *vr, *vre;
	linpart *lp, **plp;
#endif

	d += k;
#ifndef Just_Linear
	S = (Static *)R->S;
	asl = S->asl;
	if (cexp1st) {
		j1 = *S->c1s;
		d->af1 = S->com1 + j1;
		*++S->c1s = j2 = last_cex - nv011;
		d->afn = j2 - j1;
		if ((j2 -= j1) > 0) {
			c1x = cexps1 + j1;
			for(i = j = 0; j < j2; ++j)
				if (c1x[j].lp)
					++i;
			if (i) {
				d->c1lp = plp = (linpart**)mem((i+1)*sizeof(cexp1*));
				do {
					if ((lp = c1x[--j2].lp))
						*plp++ = lp;
					} while(j2 > 0);
				*plp = 0;
				}
			}
		}
	if (!lastj) {
		alast = max_var;
		S->afree = S->afree0;
		}
	lastj = 0;
#endif /* Just_Linear */
	d->o.e = eread1(R, wd ? &d->db : 0, 0, 0);
#ifndef Just_Linear
	if (amax1 < alast)
		amax1 = alast;
	if ((nz = nzc)) {
		nzc = 0;
		nv = nv0;
		c = zc;
		ci = zci;
		for(i = 0; i < nz; ) {
			if ((j = ci[i]) < nv) {
				ci[i] = ci[--nz];
				c[j] = 0;
				}
			else
				++i;
			}
		if (nz) {
			cx0 = cexps;
			ndv0 = asl->i.defvar0;
			for(i = 0; i < nz; ++i) {
				di = &cx0[ci[i] - ndv0];
				if ((vr = di->vref) && vr[2] <= 1) {
					nz1 = vr[0];
					j1 = vr[1];
					vr += 3;
					vre = vr + nz1;
					vr += j1;
					while(vr < vre) {
						if (!c[j2 = *vr++]) {
							ci[nz++] = j2;
							c[j2] = nz;
							}
						}
					}
				}
			if (nz > 1)
				qsortv(ci, nz, sizeof(int), compar, NULL);
			d->dvref = dv = (int*)mem((nz+1)*sizeof(int));
			dv[0] = nz;
			memcpy(dv+1, ci, nz*sizeof(int));
			if ((db = d->db))
				derpcopy(S, 0, 0, nz, ci, db);
			while(nz > 0)
				c[ci[--nz]] = 0;
			}
		}
#endif
	}

#ifndef Just_Linear
 static void
linpt_read(EdRead *R, int nlin, linpart **plp, int a0)
{
	ASL *asl;
	int (*Xscanf)(EdRead*, const char*, ...);
	lincoef *lc, *lce;
	linpart *lp;
	uint u;

	if (nlin <= 0) {
		*plp = 0;
		return;
		}
	asl = R->asl;
	*plp = lp = (linpart*)mem(sizeof(linpart) + (nlin-1)*sizeof(lincoef));
	lc = lp->lc;
	lce = lc + nlin;
	lp->a = a0;
	lp->n = nlin;
	Xscanf = xscanf;
	for(;;) {
		if (Xscanf(R, "%d %lf", &u, &lc->coef) != 2)
			badline(R);
		lc->varno = u;
		if (++lc >= lce)
			break;
		}
	}

 static void
cexp_read(EdRead *R, int kc, int nlin)
{
	ASL_fg *asl;
	Static *S;
	cexp *ce, *cx0;
	derpblock *db;
	int *c, *ci;
	int fmax, i, i1, j, k, nv, nz, nz0, *vr, *vre;
	lincoef *lc, *lce;
	linpart *lp;
	uint ia, ncond0, nd0, nd;

	S = (Static *)R->S;
	asl = S->asl;
	nv = nv0;
	ce = (cx0 = cexps) + (kc - nv);
	linpt_read(R, nlin, &ce->lp, kc);
	alast = S->afirst;
	S->afree = S->afree0;
	ncond0 = S->ncond;
	nd0 = S->nderp;
	ce->o.e = eread1(R, 0, &ia, kc+1);
	if (amax1 < alast)
		amax1 = alast;
	nz = nzc;
	c = zc;
	ci = zci;
	if ((lp = ce->lp)) {
		lc = lp->lc;
		lce = lc + lp->n;
		j = kc + 1;
		do {
			if (!c[i1 = lc->varno]) {
				ci[nz++] = i1;
				c[i1] = 1;
				}
			new_derp(S, i1 + 1, j, numind(S, lc->coef));
			}
			while(++lc < lce);
		}
	db = 0;
	if (S->curdb.d0 != S->curdb.de || S->curdb.nxt) {
		db = (derpblock*)mem(sizeof(derpblock));
		*db = S->curdb;
		db_reset(S);
		}
	ce->db = ce->dbf = db;
	if ((nz0 = nz)) {
		nd = S->nderp - nd0;
		if ((fmax = maxfwd) <= 0)
			goto f1_check;
		for(i = k = 0; i < nz; ++i) {
			if ((j = ci[i]) < nv) {
				if (++k > fmax)
					goto f1_check;
				}
			else if ((vr = cx0[j-nv].vref)) {
				if (vr[2] == 2)
					i1 = vr[1];
				else
					i1 = vr[0];
				vr += 3;
				for(vre = vr + i1; vr < vre; ) {
					if (!c[i1 = *vr++]) {
						c[i1] = 1;
						ci[nz++] = i1;
						}
					}
				}
			}
		if (nd > 3*k && ncond0 == S->ncond) { /* funnelkind 2 */
			if (nz > 1)
				qsortv(ci, nz, sizeof(int), compar, NULL);
			ce->vref = vr = (int*)mem((nz+3)*sizeof(int));
			vr[0] = nz;
			vr[1] = k;
			vr[2] = 2;
			memcpy(vr+3, ci, nz*sizeof(int));
			derpcopy(S, kc, k, nz, ci, db);
			}
		else {
 f1_check:
			if (nz0 > 1)
				qsortv(ci, nz0, sizeof(int), compar, NULL);
			for(k = 0; k < nz0 && ci[k] < nv; ++k);
			ce->vref = vr = (int*)mem((nz0+3)*sizeof(int));
			vr[0] = nz0;
			vr[1] = k;
			vr[2] = 0;
			memcpy(vr+3, ci, nz0*sizeof(int));
			if (nd > 3*nz0 || ncond0 != S->ncond) {
				/* funnelkind 1 */
				vr[2] = 1;
				/* no derpcopy(S, kc, nz0, nz0, ci); */
				}
			/* else no funnel */
			}
		while(nz > 0)
			c[ci[--nz]] = 0;
		}
	nzc = 0;
	S->firstdb1 = S->firstdb0;
	}

 static void
cexp1_read(EdRead *R, int j, int k, int nlin)
{
	ASL_fg *asl;
	Static *S;
	cexp1 *ce;
	linpart *lp;
	uint ia, ka;

	S = (Static *)R->S;
	asl = S->asl;
	ce = cexps1 + (k - nv01);
	linpt_read(R, nlin, &ce->lp, 0);

	if (!lastj) {
		alast = max_var;
		S->afree = S->afree0;
		lastj = j;
		}
	ia = 0;
	ce->o.e = eread1(R, 0, &ia, ka = k + 1);
	if ((lp = ce->lp) || ia) {
		if (ia) {
			if (ia > S->afirst)
				free_a(S, ia);
			else if (ia != ka)
				new_derp(S, ia, ka, S->one);
			}
		if (lp)
			lp->a = ka - 1;
		}
	S->c1a[k] = ka;
	last_cex = k;
	}

#endif /* Just_Linear */

 static cexp **
funnels(Static *S, cexp *cx0, cexp *cx1, int ginit)
{
	ASL_fg *asl;
	cexp *cx, **rv, **t;
	int n, *vr;

	asl = S->asl;
	for(n = 0, cx = cx0; cx < cx1; ++cx) {
		if ((vr = cx->vref) && vr[2] == 2)
			++n;
		}
	if (!n)
		return 0;
	asl->i.x0kindinit |= ginit;
	t = rv = (cexp**)mem((n+1)*sizeof(cde*));
	for(cx = cx0; cx < cx1; ++cx) {
		if ((vr = cx->vref) && vr[2] == 2)
			*t++ = cx;
		}
	*t = 0;
	return rv;
	}

 static void
adjust_compl_rhs(Static *S)
{
	ASL_fg *asl;
	cde *C;
	int *Cvar, i, j, nc, *o, stride;
	real *L, *U, t, t1;

	asl = S->asl;
	L = LUrhs;
	if ((U = Urhsx))
		stride = 1;
	else {
		U = L + 1;
		stride = 2;
		}
	C = con_de;
	Cvar = cvar;
	nc = n_con;
	for(i = nlc; i < nc; i++) {
		if (Cvar[i] && (o = C[i].o.e) && *o == 0
		&& (t = S->htvals_end[o[1]]) != 0.) {
			t1 = t;
			if (L[j = stride*i] > negInfinity) {
				L[j] -= t;
				t1 = 0.;
				}
			if (U[j] < Infinity) {
				U[j] -= t;
				t1 = 0.;
				}
			o[1] = numind(S, t1);
			}
		}
	}

 static void
adjust(Static *S, int flags)
{
	ASL_fg *asl;
	cexp *cx0, *cx1;
	size_t L;
	void **vp, *vp1;

	asl = S->asl;
	if (S->pinvd) {
		asl->i.invd = (void**)mem(L = asl->i.nfinv * sizeof(Invd1*));
		memcpy(asl->i.invd, S->pinvd, L);
		del_mblk(S->pinvd);
		}
	if (S->k_seen) {
		if (!A_vals)
			goff_comp_ASL((ASL*)asl);
		else if (Fortran)
			colstart_inc_ASL((ASL*)asl);
		}
	if (n_cc > nlcc && nlc < n_con
	 && !(flags & ASL_no_linear_cc_rhs_adjust))
		adjust_compl_rhs(S);
	if (ncom0) {
		cx0 = cexps;
		if (comb) {
			asl->I.dvfb = funnels(S, cx0, cx1 = cx0 + comb, ASL_need_comba);
			cx0 = cx1;
			}
		if (comc) {
			asl->I.dvfc = funnels(S, cx0, cx1 = cx0 + comc, ASL_need_comca);
			cx0 = cx1;
			}
		if (como)
			asl->I.dvfo = funnels(S, cx0, cx0 + como, ASL_need_comoa);
		}
#ifndef Just_Linear
	asl->i.derplen = amax1 * sizeof(real);
#endif
	asl->i.wlen = wlast*sizeof(real);
	asl->i.numlen = L = S->numht_n * sizeof(real);
	for(vp = S->numhtf0;; vp = (void**)vp1) {
		vp1 = *vp;
		free(vp);
		if (!vp1)
			break;
		}
	memcpy(asl->i.numvals = (real*)M1alloc(L), S->htvals_end - S->numht_n, L);
	free(S->htvals);
	free(S->numht);
	}

 static void
br_read(EdRead *R, int nc, real *L, real *U, int *Cvar, int nv)
{
	int i, inc, j, k;
	int (*Xscanf)(EdRead*, const char*, ...);
	ASL *asl = R->asl;

	if (U)
		inc = 1;
	else {
		U = L + 1;
		inc = 2;
		}
	Xscanf = xscanf;
	Xscanf(R, ""); /* purge line */
	for(i = 0; i < nc; i++, L += inc, U += inc) {
		switch(edag_peek(R) - '0') {
		  case 0:
			if (Xscanf(R,"%lf %lf",L,U)!= 2)
				badline(R);
			break;
		  case 1:
			if (Xscanf(R, "%lf", U) != 1)
				badline(R);
			*L = negInfinity;
			break;
		  case 2:
			if (Xscanf(R, "%lf", L) != 1)
				badline(R);
			*U = Infinity;
			break;
		  case 3:
			*L = negInfinity;
			*U = Infinity;
			Xscanf(R, ""); /* purge line */
			break;
		  case 4:
			if (Xscanf(R, "%lf", L) != 1)
				badline(R);
			*U = *L;
			break;
		  case 5:
			if (Cvar) {
				if (Xscanf(R, "%d %d", &k, &j) == 2
				 && j > 0 && j <= nv) {
					Cvar[i] = j;
					*L = k & 2 ? negInfinity : 0.;
					*U = k & 1 ?    Infinity : 0.;
					break;
					}
				}
		  default:
			badline(R);
		  }
		}
	}

#ifndef Just_Linear

 static int
aholread(EdRead *R)
{
	FILE *nl = R->nl;
	Static *S = (Static *)R->S;
	char *s1;
	int i, k, *op, rv;

	k = getc(nl);
	if (k < '1' || k > '9')
		badline(R);
	i = k - '0';
	while((k = getc(nl)) != ':') {
		if (k < '0' || k > '9')
			badline(R);
		i = 10*i + k - '0';
		}
	k = i/sizeof(int) + 4;
	op = nextop(S, k);
	op[0] = nOPHOL;
	op[1] = rv = wlast++;
	op[2] = k;
	for(s1= (char*)&op[3];;) {
		if ((k = getc(nl)) < 0) {
			fprintf(Stderr,
				 "Premature end of file in aholread, line %ld of %s\n",
					R->Line, R->asl->i.filename_);
				exit_ASL(R,1);
			}
		if (k == '\n') {
			R->Line++;
			if (!i)
				break;
			}
		if (--i < 0)
			badline(R);
		*s1++ = k;
		}
	*s1 = 0;
	return rv;
	}

 static int
bholread(EdRead *R)
{
	ASL_fg *asl;
	Static *S;
	char *s;
	int i, k, *op, rv;

	S = (Static *)R->S;
	asl = S->asl;
	if (xscanf(R, "%d", &i) != 1)
		badline(R);
	k = i/sizeof(int) + 4;
	op = nextop(S, k);
	op[0] = nOPHOL;
	op[1] = rv = wlast++;
	op[2] = k;
	s = (char*)&op[3];
	if (fread(s, i, 1, R->nl) != 1)
		badline(R);
	s[i] = 0;
	for(;;) switch(*s++) {
			case 0: goto break2; /* so we return at end of fcn */
			case '\n': R->Line++;
			}
 break2:
	return rv;
	}

 static real
Missing_func(arglist *al)
{
	AmplExports *ae = al->AE;
	func_info *fi = (func_info*)al->funcinfo;

	char *s = (char*)(*ae->Tempmem)(al->TMI, strlen(fi->name) + 64);
	(*ae->SprintF)(al->Errmsg = s,
		"Attempt to call unavailable function %s.",
		fi->name);
	return 0.;
	}
#endif /* Just_Linear */

 int
fg_read_ASL(ASL *a, FILE *nl, int flags)
{
	ASL_fg *asl;
	EdRead ER, *R;
	Static SS, *S;
	Jmp_buf JB;
	cgrad *cg, **cgp;
	int i, i1, j, k, *ka, kseen, nc, nc0, nco, nlcon, no;
	int nv1, nvc, nvo, nvr, nxv, *o, readall;
	int (*Xscanf)(EdRead*, const char*, ...);
	ograd *og, **ogp;
	real *oc, t;
	size_t LL[3], *kaz, nz;
	unsigned x;
#ifdef Just_Linear
#define ASL_readtype ASL_read_f
#else /* Just_Linear */
#define ASL_readtype ASL_read_fg
#undef ncom1
	int ncom, ncom1, nv;
	func_info *fi;
	char fname[128];
	int nlin;
#endif /* Just_Linear */

	ASL_CHECK(a, ASL_readtype, who);
	flagsave_ASL(a, flags); /* includes allocation of LUv, LUrhs, A_vals or Cgrad, etc. */
	S = ed_reset(&SS, a);
	asl = (ASL_fg*)a;
	i = comc ? ASL_need_concom : 0;
	if (como)
		i |= ASL_need_objcom;
	asl->i.x0kindinit = i;
	ncom0 = combc + como;
#ifndef Just_Linear
	SS.com1 = a->i.n_var0 + ncom0;
	ncom = comb + comc + como + comc1 + como1;
#endif
	nvr = n_var; /* nv for reading */
	nxv = a->i.nsufext[ASL_Sufkind_var];
	nv0 = nv1 = nvr + nxv;
#ifndef Just_Linear
	SS._wlast = nv1 + ncom;
	nv01 = a->i.n_var0 + ncom0;
#endif
	SS.numht_mask = (1 << 10) - 1;
	nz = SS.numht_mask + 1;
#ifndef Just_Linear
	asl->i.ncom1_ = ncom1 = comc1 + como1;
#endif
	SS.numht = (Numhash**)Malloc(LL[0] = nz*sizeof(Numhash*));
	memset(SS.numht, 0, nz*sizeof(Numhash*));
	SS.htvals = (real*)Malloc(LL[1] = nz*sizeof(real));
	SS.htvals_end = SS.htvals + nz;
	SS.numhtpth = &SS.numhtfirst;
	SS.numhtf0 = (void**)Malloc(LL[2] = nz*sizeof(void*) + ncom1*sizeof(uint));
	asl->i.temp_rd_bytes = LL[0] + LL[1] + LL[2];
#ifndef Just_Linear
	SS.c1a = (uint*)(SS.numhtf0 + nz) - nv01; /* for convenience in eread() */
#endif
	*SS.numhtf0 = 0;
	SS.numhtf = (Numhash*)(SS.numhtf0 + 1);
	SS.numhtfend = SS.numhtf + ((nz-1)*sizeof(void*)/sizeof(Numhash));
	numind(S, 0.); /* for mpec_adj(): w[-1] = 0. */
#ifndef Just_Linear
	SS.one = numind(S, 1.);
#endif
#if 1 /*TEMP DEBUG*/
	SS.R =
#endif
	R = EdReadInit_ASL(&ER, a, nl, S);
	if (flags & ASL_return_read_err) {
		a->i.err_jmp_ = &JB;
		i = setjmp(JB.jb);
		if (i) {
			a->i.err_jmp_ = 0;
			return i;
			}
		}
	nlcon = a->i.n_lcon_;
	if (nlcon && !(flags & ASL_allow_CLP)) {
		if (a->i.err_jmp_)
			return ASL_readerr_CLP;
		sorry_CLP(R, "logical constraints");
		}
	readall = flags & ASL_keep_all_suffixes;
	Xscanf = a->i.xscanf_;
	ER.lineinc = 1;
#ifndef Just_Linear
	SS.firstdb0 = SS.firstdb1 = (derpblock**)new_mblk(SS.kfirstdb = 3);
	SS.firstdbe = SS.firstdb0 + 8;
	if (nfunc)
		func_add(a);
	if (binary_nl)
		holread = bholread;
	else
		holread = aholread;

#endif /* Just_Linear */
	nc0 = n_con;
	nc = nc0 + a->i.nsufext[ASL_Sufkind_con];
	no = n_obj;
	nvc = c_vars;
	nvo = o_vars;
	nco = nc + no + nlcon;
	if (no < 0 || nco <= 0)
		scream(R, ASL_readerr_corrupt,
			"edagread: nc = %d, no = %d, nlcon = %d\n",
			nc0, no, nlcon);
	if (pi0) {
		memset(pi0, 0, nc*sizeof(real));
		if (havepi0)
			memset(havepi0, 0, nc);
		}
#ifdef Just_Linear
	x = nco*sizeof(cde) + no*sizeof(ograd *) + no;
#else
	S->afirst = alast = max_var = asl->i.maxvar = nv = nv1 + ncom;
	nv0b = nv1 + comb;
	nv0c = nv0b + comc;
	last_cex = nv011 = nv01 - 1;
	x = nco*sizeof(cde) + no*sizeof(ograd *)
		+ nv*(2*sizeof(int))
		+ ncom0*sizeof(cexp)
		+ ncom1*sizeof(cexp1)
		+ nfunc*sizeof(func_info *)
		+ no;
	if (ncom1)
		x += (nco + 1)*sizeof(int);
	SS.nvar0 = a->i.n_var0;
	asl->i.defvar0 = SS.nvar0 + nxv;
	if (!(SS.nvinc = a->i.n_var_ - SS.nvar0 + nxv))
		SS.nvar0 += ncom0 + ncom1;
#endif /* Just_Linear */
	if (X0)
		memset(X0, 0, nvr*sizeof(real));
	if (havex0)
		memset(havex0, 0, nvr);
	asl->i.objconst = oc = (real*)M1zapalloc(x + no*sizeof(real));
	con_de = (cde *)(oc + no);
	lcon_de = con_de + nc;
	obj_de = lcon_de + nlcon;
	Ograd = (ograd **)(obj_de + no);
#ifdef Just_Linear
	objtype = (char *)(Ograd + no);
#else
	cexps = (cexp *)(Ograd + no);
	cexps1 = (cexp1 *)(cexps + ncom0);
	funcs = (func_info **)(cexps1 + ncom1);
	zc = (int*)(funcs + nfunc);
	zci = zc + nv;
	ka = zci + nv;
	c_cexp1st = l_cexp1st = o_cexp1st = 0;
	if (ncom1) {
		SS.c1s = ka;
		*ka = 0;
		if (comc1) {
			c_cexp1st = ka;
			ka += nc;
			l_cexp1st = ka;
			ka += nlcon;
			}
		if (como1) {
			*ka = comc1;
			o_cexp1st = ka;
			ka += no;
			}
		++ka;
		}
	objtype = (char *)ka;
#endif
	if (n_cc && !cvar)
		cvar = (int*)M1alloc(nc*sizeof(int));
	if (cvar)
		memset(cvar, 0, nc*sizeof(int));
	ka = 0;
	kaz = 0;
	nz = 0;
	j = kseen = 0;
	for(;;) {
		ER.can_end = 1;
		i = edag_peek(R);
		if (i == EOF) {
			adjust(S, flags);
			nzjac = nz;
			fclose(nl);
			a->i.opunused = (SS._oplast - SS._opnext)*sizeof(int);
#ifndef Just_Linear
			a->i.derpunused = (SS.curdb.d0 - SS.first_d)*sizeof(derp);
			a->p.Objval  = a->p.Objval_nomap  = obj1val_ew_ASL;
			a->p.Objgrd  = a->p.Objgrd_nomap  = obj1grd_ew_ASL;
			a->p.Conval  = con1val_ew_ASL;
			a->p.Jacval  = jac1val_ew_ASL;
			a->p.Conival = a->p.Conival_nomap = con1ival_ew_ASL;
			a->p.Congrd  = a->p.Congrd_nomap  = con1grd_ew_ASL;
			a->p.Lconval = lcon1val_ew_ASL;
			a->p.Xknown  = x1known_ew_ASL;
#endif /* Just_Linear */
			a->p.EWalloc = ewalloc1_ASL;
			return prob_adj_ASL(a);
			}
		ER.can_end = 0;
		k = -1;
		switch(i) {
			case 'C':
				Xscanf(R, "%d", &k);
				if (k < 0 || k >= nc0)
					badline(R);
				co_read(R,con_de,k,want_derivs,c_cexp1st);
				break;
#ifdef Just_Linear
			case 'F':
			case 'V':
				sorry_nonlin(R);
#else
			case 'F':
				if (Xscanf(R, "%d %d %d %127s",
						&i, &j, &k, fname) != 4
				|| i < 0 || i >= nfunc)
					badline(R);
				if ((fi = func_lookup(a, fname,0))) {
					if (fi->nargs != k && fi->nargs >= 0
					 && (k >= 0 || fi->nargs < -(k+1)))
						scream(R, ASL_readerr_argerr,
				"function %s: disagreement of nargs: %d and %d\n",
					 		fname,fi->nargs, k);
					}
				else {
					fi = (func_info *)mem(sizeof(func_info));
					fi->ftype = j;
					fi->nargs = k;
					fi->funcp = 0;
					fi->name = (Const char *)strcpy((char*)
						mem(strlen(fname)+1), fname);
					}
				if (!fi->funcp && !(fi->funcp = dynlink(fname))){
					if (!(flags & ASL_allow_missing_funcs))
					    scream(R, ASL_readerr_unavail,
						"function %s not available\n",
						fname);
					fi->funcp = Missing_func;
					fi->funcinfo = (Char*)fi;
					}
				funcs[i] = fi;
				break;
			case 'L':
				Xscanf(R, "%d", &k);
				if (k < 0 || k >= nlcon)
					badline(R);
				co_read(R, lcon_de, k, 0, l_cexp1st);
				break;
			case 'V':
				if (Xscanf(R, "%d %d %d", &k, &nlin, &j) != 3)
					badline(R);
				if (k >= SS.nvar0)
					k += SS.nvinc;
				if (k < nvr || k >= nv)
					badline(R);
				if (j)
					cexp1_read(R, j, k, nlin);
				else
					cexp_read(R, k, nlin);
				break;
#endif /* Just_Linear */
			case 'G':
				if (Xscanf(R, "%d %d", &j, &k) != 2
				|| j < 0 || j >= no || k <= 0 || k > nvo)
					badline(R);
				ogp = Ograd + j;
				while(k--) {
					*ogp = og = (ograd *)mem(sizeof(ograd));
					ogp = &og->next;
					if (Xscanf(R, "%d %lf", &i, &og->coef) != 2)
						badline(R);
					og->varno = i;
					}
				*ogp = 0;
				break;
			case 'J':
				if (Xscanf(R, "%d %d", &j, &k) != 2
				|| j < 0 || j >= nc0 || k <= 0 || k > nvc)
					badline(R);
				nz += k;
				if (A_vals) {
					j += Fortran;
					if (ka) {
						while(k--) {
							if (Xscanf(R, "%d %lf",
								&i, &t) != 2)
								badline(R);
							i1 = ka[i]++;
							A_vals[i1] = t;
							A_rownos[i1] = j;
							}
						}
					else {
						while(k--) {
							if (Xscanf(R, "%d %lf",
								&i, &t) != 2)
								badline(R);
							i1 = kaz[i]++;
							A_vals[i1] = t;
							A_rownos[i1] = j;
							}
						}
					break;
					}
				cgp = Cgrad + j;
				while(k--) {
					*cgp = cg = (cgrad *)mem(sizeof(cgrad));
					cgp = &cg->next;
					if (kseen) {
						if (Xscanf(R, "%d %lf", &i, &cg->coef) != 2)
							badline(R);
						}
					else
						if (Xscanf(R, "%d %d %lf", &i, &j,
							    &cg->coef) != 3)
							badline(R);
					cg->varno = i;
					cg->goff = j;
					}
				*cgp = 0;
				break;
			case 'O':
				if (Xscanf(R, "%d %d", &k, &j) != 2
				 || k < 0 || k >= no)
					badline(R);
				objtype[k] = j;
				co_read(R,obj_de,k,want_derivs,o_cexp1st);
				o = obj_de[k].o.e;
				if (o[0] == nOPRET && (i = o[1]) < 0)
					oc[k] = S->htvals_end[i];
				break;
			case 'S':
				Suf_read_ASL(R, readall);
				break;
			case 'r':
				br_read(R, asl->i.n_con0, LUrhs, Urhsx, cvar, nvr);
				break;
			case 'b':
				br_read(R, asl->i.n_var0, LUv, Uvx, 0, 0);
				break;
			case 'K':
			case 'k':
				SS.k_seen = ++kseen;
				if (ka_read_ASL(a, R, i, &ka, &kaz))
					badline(R);
				break;
			case 'x':
				if (!Xscanf(R,"%d",&k)
				|| k < 0 || k > nvr)
					badline(R);
				if (!X0 && want_xpi0 & 1) {
					x = nv1*sizeof(real);
					if (want_xpi0 & 4)
						x += nv1;
					X0 = (real *)M1zapalloc(x);
					if (want_xpi0 & 4)
						havex0 = (char*)(X0 + nv1);
					}
				while(k--) {
					if (Xscanf(R, "%d %lf", &j, &t) != 2
					 || j < 0 || j >= nvr)
						badline(R);
					if (X0) {
						X0[j] = t;
						if (havex0)
							havex0[j] = 1;
						}
					}
				break;
			case 'd':
				if (!Xscanf(R,"%d",&k)
				|| k < 0 || k > nc0)
					badline(R);
				if (!pi0 && want_xpi0 & 2) {
					x = nc*sizeof(real);
					if (want_xpi0 & 4)
						x += nc;
					pi0 = (real *)M1zapalloc(x);
					if (want_xpi0 & 4)
						havepi0 = (char*)(pi0 + nc);
					}
				while(k--) {
					if (Xscanf(R, "%d %lf", &j, &t) != 2
					 || j < 0 || j >= nc0)
						badline(R);
					if (pi0) {
						pi0[j] = t;
						if (havepi0)
							havepi0[j] = 1;
						}
					}
				break;
			default:
				badline(R);
			}
		}
	}
