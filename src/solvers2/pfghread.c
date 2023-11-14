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

#define PSHVREAD
#define pfg_read_ASL pfgh_read_ASL

#ifdef DEBUG
#include "assert.h"
#else
#define assert(x) /*nothing*/
#endif

#ifdef PSHVREAD
#include "jacpdim.h"
#include "opnos.hd"
#define PSHV(x) x
#define asltype ASL_read_pfgh
#define who "pfgh_read"
/*#define ASLTYPE ASL_pfgh*/
#else
#define PSHV(x) /*nothing*/
#define ETYPE expr
#include "asl_pfg.h"
#define asltype ASL_read_pfg
#define who "pfg_read"
#define Varval Varval1
/*#define ASLTYPE ASL_pfg*/
#endif

#include "opcode.hd"
#include "opno2.h"

#ifdef X64_bit_pointers
#define needalign(op,n) ((size_t)&op[n]) & (sizeof(void*)-1)
#define alignnum(a,b) (b)
#else
#define needalign(op,n) 0
#define alignnum(a,b) (a)
#endif
#define Rnumberof(x) (((x) + sizeof(real) - 1)/sizeof(real))

 struct
expr {
	int op;
	pei L, R;
	};

 typedef struct
exprc {
	int op;
	pei L, R;
	int c;
	} exprc;

 typedef struct
expr_if {
	int op;
	expr *test, *tval, *fval;
	} expr_if;

 struct
la_ref {
	la_ref *next;
	expr **ep;
	real c;
	real scale;
	};

#ifdef __cplusplus
extern "C" {
#endif

extern real conpival_nomap_ew_ASL(EvalWorkspace*, int, real*, fint*);
extern void conpgrd_nomap_ew_ASL(EvalWorkspace*, int, real *, real*, fint*);

#ifdef PSHVREAD
#define Static Static_pshv
#else
#define Static Static_psed
#endif
typedef struct Static Static;
 static ograd *cotermwalk(Static*, expr **ep, ps_func *f, int wantg, int omitdv);

#define NFHASH 23

#undef f_OPNUM
#include "r_opn0.hd"

#undef nzc

#if defined(IEEE_MC68k) || defined(IBM)
  enum { W0 = 0, W1 = 1 };
#else
  enum { W0 = 1, W1 = 0 };
#endif

/*******/

 typedef struct
expr_vx {
	expr_nv v;
	linarg *la;
	int a0;
	int a1;
	} expr_vx;

 typedef struct
Elemtemp {
	uint esize;
	int nmax;
	int k;
	void **mp;
	} Elemtemp;

 typedef struct
PSfind {
	ps_func *f;
	Elemtemp *b, *g;
	int nb, ng;
	} PSfind;

 typedef struct
Numhash { struct Numhash *nhnext, *nhthread; real *nhval; } Numhash;

 typedef union { real r; uint ui[2]; } UIR;
 typedef struct
Wkinit {
	int k, n, nmax;
	int x[1]; /* really x[nmax] */
	} Wkinit;

 struct
Static {
	ASLTYPE *asl;
	ASL *a;
	Elemtemp *_last_psb_elem;
	Numhash **numht, *numhtf, *numhtfend, *numhtfirst, **numhtpth;
	Wkinit *wk0, *wk2, *wkm1;
	char *digc, *efnext, *eflast, *etofree;
	derpblock curdb, *emptydb, **firstdb0, **firstdb1, **firstdbe;
	derp *first_d;
	expr *(*_holread)(EdRead*);
	expr *expr_free;
	expr_if *exprif_free;
	expr_nv *exprn_free;
	expr **slscratch;	/* used in awalk(sumlist) */
	expr *_last_e;
	exprc *exprc_free;
	expr_nv **_larvlist;
	expr_nv *var_e, **_varp;
	int *_oplast, *_opnext, *_zc, *_zci, *_zl, *c1s, *dvspb;
	int *opfirst, *opnext0, *opprev, *zc1, *zci1;
	real *htvals, *htvals_end;
	int _amax1, _last_cex, _lastj, _max_var;
	int _nv0b, _nv0c, _wlast, com1, combco, dv0, dvsp1, dvsp2, fmax, gscrx;
	int k_afree, k_seen, kfirstdb, maxspvar;
	int ncomo, ndvbcosp, negone, nhop, nopblks, nvar0, nvinc, one, two;
	size_t _lthashmask, _nrange, _rangehashmask, uhlen;
	int klthash, krangehash;
	int _allJ;
	int _cexp_k;
	int _cexp_n;
	int _conno;
	int _groupno;
	int _k_Elemtemp;
	int _kzc;
	int _lasta;
	int _ncom;
	int _ncom_togo;
	int _nocopy;
	int _nndv;
	int _nsce;
	int _nv0x;
	int _nzclim;
	int slmax;
	int slscratchlev;
	int slscratchmax;
	int _termno;		/* current term number */
	int _wantCgroups;
	int _wantOgroups;
	int _zc_lim;
	int *atc, **dop, **dop0;
	la_ref *_laref_free;
	linarg **larep;
	linarg **_lthash;		/* hash table */
	linarg *_ltfree;		/* free list */
	linarg *_tlist;		/* linargs in this term */
	ograd *_freeog;		/* free list */
	range **_rangehash;
	real *_rnz;
	real _lt_scale;
	size_t atlen, diglen;
	tfinfo **ptfi;
	uint afirst, _alast, _nzc, alast0, ncond, nderp, numht_mask, numht_n;
	uint *afree, *afree0, *afree1, *c1a, *oval;
	void **numhtf0;
	EdRead *R;
	};

#define alast		S->_alast
#define allJ		S->_allJ
#define amax1		S->_amax1
#define cexp_k		S->_cexp_k
#define cexp_n		S->_cexp_n
#define Conno		S->_conno
#define freeog		S->_freeog
#define Groupno		S->_groupno
#define holread		S->_holread
#define if2list		S->_if2list
#define if2list_end	S->_if2list_end
#define k_Elemtemp	S->_k_Elemtemp
#define kzc		S->_kzc
#define laref_free	S->_laref_free
#define larvlist	S->_larvlist
#define last_e		S->_last_e
#define last_psb_elem	S->_last_psb_elem
#define lasta		S->_lasta
#define lt_scale	S->_lt_scale
#define ltfree		S->_ltfree
#define lthash		S->_lthash
#define lthashmask	S->_lthashmask
#define max_var		S->_max_var
#define max_var1	asl->P.max_var1_
#define Ncom		S->_ncom
#define ncom_togo	S->_ncom_togo
#define nocopy		S->_nocopy
#define nndv		S->_nndv
#define nrange		S->_nrange
#define nsce		S->_nsce
#define nv0b		S->_nv0b
#define nv0c		S->_nv0c
#define nv0x		S->_nv0x
#define nzc		S->_nzc
#define nzclim		S->_nzclim
#define oplast		S->_oplast
#define opnext		S->_opnext
#define rangehash	S->_rangehash
#define rangehashmask	S->_rangehashmask
#define Termno		S->_termno
#define tlist		S->_tlist
#define varp		S->_varp
#define wantCgroups	S->_wantCgroups
#define wantOgroups	S->_wantOgroups
#define wlast		S->_wlast
#define zc		S->_zc
#define zc_lim		S->_zc_lim
#define zci		S->_zci
#define zl		S->_zl

#ifdef PSHVREAD
 static void hes_setup(Static *);
#else
#define hes_setup(S) /*nothing*/
#endif

 static Static *
S_init(Static *S, ASLTYPE *asl)
{
	memset(S, 0, sizeof(Static));
	S->asl = asl;
	S->a = (ASL*)asl;
	return S;
	}

 static void *
new_mblkzap(ASLTYPE *asl, int k)
{
	void *rv = new_mblk_ASL((ASL*)asl,k);
	memset(rv, 0, sizeof(void*)<<k);
	return rv;
	}

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
	ASLTYPE *asl;
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

 static void
ed_reset(ASLTYPE *asl)
{
	asl->i.memLast = asl->i.memNext = 0;
	memset(asl->i.mblk_free, 0, MBLK_KMAX*sizeof(char*));
	memset(&asl->P.merge + 1, 0, sizeof(ps_info) - sizeof(long));
	}

 static void
wk_add(ASLTYPE *asl, Wkinit **pwk, int i)
{
	Wkinit *wk, *wk1;
	int k, *x;

	if (!(wk = *pwk)) {
		k = htcl(sizeof(Wkinit) + 5*sizeof(int));
		wk = (Wkinit*)new_mblk(k);
		wk->n = 0;
		goto set_nmax;
		}
	else if (wk->n >= wk->nmax) {
		wk1 = (Wkinit*)new_mblk(k = wk->k + 1);
		memcpy(wk1, wk, sizeof(Wkinit) + (wk->n-1)*sizeof(int));
		del_mblk(wk);
		wk = wk1;
 set_nmax:
		wk->k = k;
		wk->nmax = ((sizeof(void*)<<k) - sizeof(Wkinit))/sizeof(int) + 1;
		*pwk = wk;
		}
	x = wk->x;
	x[wk->n++] = i;
	}

#ifdef PSHVREAD
 static int*
wksave(ASLTYPE *asl, Wkinit **pwk)
{
	Wkinit *wk;
	int n, *z;

	if (!(wk = *pwk))
		return 0;
	n = wk->n;
	z = (int*)mem((n + 1)*sizeof(int));
	*z = wk->n;
	memcpy(z+1, wk->x, n*sizeof(int));
	del_mblk(wk);
	*pwk = 0;
	return z;
	}
#endif

#ifdef Double_Align
#define memadj(x) x
#else
#define memadj(x) (((x) + (sizeof(long)-1)) & ~(sizeof(long)-1))
#endif
 static Elemtemp *
new_Elemtemp(Static *S, uint esize, void **mp)
{
	Elemtemp *e;
	int k;
	ASL *asl = S->a;

	e = (Elemtemp *)new_mblk(k_Elemtemp);
	e->esize = esize;
	e->mp = mp;
	e->k = k = htcl(8*esize);
	*mp = new_mblk(k);
	e->nmax = (sizeof(void*) << k) / esize;
	return e;
	}

 static void
del_Elemtemp(Static *S, Elemtemp *e)
{
	ASL *asl = S->a;
	del_mblk(*e->mp);
	del_mblk(e);
	}

 static void
upgrade_Elemtemp(Static *S, Elemtemp *e)
{
	void *m, *m0;
	int k;
	ASL *asl = S->a;

	k = e->k++;
	memcpy(m = new_mblk(e->k), m0 = *e->mp, e->esize * e->nmax);
	del_mblk(m0);
	*e->mp = m;
	e->nmax = (sizeof(void*) << ++k) / e->esize;
	}

 static void
free_laref(Static *S, la_ref **L)
{
	la_ref *L1, *L2;

	if ((L1 = *L)) {
		while((L2 = L1->next))
			L1 = L2;
		L1->next = laref_free;
		laref_free = *L;
		*L = 0;
		}
	}

 static la_ref *
new_laref(Static *S, la_ref *nxt)
{
	la_ref *rv;

	if ((rv = laref_free))
		laref_free = rv->next;
	else
		rv = (la_ref *)mem_ASL(S->a, sizeof(la_ref));
	rv->next = nxt;
	return rv;
	}

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
	d->a = a - 1;
	d->b = b - 1;
	d->c = c;
	}

 static int opblock_gulp = 16384;

 static int*
new_opblock(Static *S, int need)
{
	ASL *asl;
	int n, *op, *opnew;

	asl = S->a;
	++need;
	++S->nopblks;
	if (!S->opfirst) {
		if ((op = S->opnext0)) {
			if (op + need <= oplast) {
				S->opprev = op;
				*op++ = OPRET;
				return S->opfirst = op;
				}
			opnext = 0;
			}
		}
	need += 4;
	if ((n = opblock_gulp) < need)
		n = need;
	opnew = (int*)M1alloc(n*sizeof(int));
	if ((op = opnext)) {
		if (needalign(op, 1))
			*op++ = OP_NEXTBLKalign;
		*op++ = OP_NEXTBLK;
		*(int**)op = opnew;
		}
	oplast = opnew + n - alignnum(3,5);
	if (!S->opfirst) {
		S->opprev = opnew;
		*opnew++ = OPRET;
		S->opfirst = opnew;
		}
	return opnew;
	}

#ifdef DEBUG
	/*DEBUG*/int zorkheswork, zorkheswork1, *zorkopnext;
#endif

 static int*
nextop(Static *S, int n)
{
	int *op, *op1;

	op = opnext;
	op1 = op + n;
	if (op1 >= oplast) {
		op = new_opblock(S, n);
		op1 = op + n;
		}
	opnext = op1;
#ifdef DEBUG
	/*DEBUG*/ if (op <= zorkopnext && op1 > zorkopnext)
	/*DEBUG*/	printf("");
#endif
	return op;
	}

 static int*
nextopp(Static *S, int n)
{
#ifdef X64_bit_pointers
	int d, *op, **x;
	long long dL;

	op = nextop(S,n+5);
	d = (int)(dL = op - S->opprev);
	if (d != dL) {
		if (needalign(op,1)) {
			op[0] = OPGOTOBalign;
			op[1] = OPGOTOB;
			x = (int**)&op[2];
			}
		else {
			op[0] = OPGOTOB;
			x = (int**)&op[1];
			}
		x[0] = S->opprev;
		S->opprev = op;
		op = (int*)&x[1];
		}
	opnext = &op[n];
#else
	int *op = nextop(S,n);
#endif
	op[1] = op - S->opprev;
	S->opprev = op;
	if (S->dop) {
		if (S->dop == S->dop0) {
			if (op == S->opfirst)
				S->opfirst = op + 1;
			*op++ = OPRETB;
			S->opprev = op;
			op[1] = 1;
			opnext = op + n;
			}
		*S->dop = op;
		S->dop = 0;
		}
	return op;
	}

 static void
efree(Static *S, expr *e)
{
	expr **args, **argse, *e1;
	expr_if *eif;
 top:
	switch(optype[e->op]) {

	 case 2: /* binary */
		efree(S, e->R.e);
		/* no break */

	 case 1: /* unary */
	 case 4: /*OPPLTERM*/
		e1 = e->L.e;
		e->L.e = S->expr_free;
		S->expr_free = e;
		e = e1;
		goto top;

	 case 3: /* OP1POW */
		e1 = e->L.e;
		e->L.e = (expr*)S->exprc_free;
		S->exprc_free = (exprc*)e;
		e = e1;
		goto top;

	 case 6: /* sumlist */
	 case 12: /* minlist,maxlist*/
		args = e->L.ep;
		argse = e->R.ep;
		while(args < argse)
			efree(S, *args++);
		e->L.e = S->expr_free;
		S->expr_free = e;
		break;

	 case 5: /*OPIFnl*/
		eif = (expr_if*)e;
		efree(S, eif->test);
		efree(S, eif->tval);
		e = eif->fval;
		eif->test = (expr*)S->exprif_free;
		S->exprif_free = eif;
		goto top;

	 case 9:
		*(expr_nv**)e = S->exprn_free;
		S->exprn_free = (expr_nv*)e;
	 }
	}

 static expr *
exprmem(Static *S, size_t L)
{
	char **b;
	expr *e, *efree;
	size_t Lf;

	L = (L + sizeof(void*) - 1) & ~(sizeof(void*) - 1);
	Lf = S->eflast - S->efnext;
	e = (expr*)S->efnext;
	if (Lf < L) {
		efree = S->expr_free;
		while(Lf >= sizeof(expr)) {
			e->L.e = efree;
			efree = e;
			e = (expr*)(S->efnext += sizeof(expr));
			Lf -= sizeof(expr);
			}
		S->expr_free = efree;
		Lf = 16384*sizeof(expr*);
		if (Lf <= L)
			Lf = L + sizeof(char**);
		S->a->i.temp_rd_bytes += Lf;
		b = (char**)Malloc(Lf);
		S->eflast = (char*)b + Lf;
		*b = S->etofree;
		S->etofree = (char*)b;
		e = (expr*)(S->efnext = (char*)&b[1]);
		}
	S->efnext += L;
	return e;
	}

 static expr *
new_expr(Static *S, int o, expr *L, expr *R)
{
	expr *LR, *rv;
	exprc *rvc;

	switch(o) {
	  case f_OPPOW:
		if (R->op == f_OPNUM) {
			if (((expr_nv *)R)->varno == S->two) {
				o = f_OP2POW;
				R = 0;
				break;
				}
			o = f_OP1POW;
			LR = R;
			}
		else if (L->op == f_OPNUM) {
			o = f_OPCPOW;
			LR = L;
			L = R;
			}
		else
			break;
 rvc_ret:
		if ((rvc = S->exprc_free))
			S->exprc_free = (exprc*)rvc->L.e;
		else
			rvc = (exprc*)exprmem(S, sizeof(exprc));
		rvc->op = o;
		rvc->L.e = L;
		rvc->R.e = 0;
		rvc->c = ((expr_nv*)LR)->varno;
		return (expr*)rvc;
	  case OP_atan2:
		if (R->op == f_OPNUM) {
			o = OPatan210;
			LR = R;
			goto rvc_ret;
			}
		if (L->op == f_OPNUM) {
			o = OPatan201;
			LR = L;
			L = R;
			goto rvc_ret;
			}
	  }
	if ((rv = S->expr_free))
		S->expr_free = rv->L.e;
	else
		rv = exprmem(S, sizeof(expr));
	rv->op = o;
	rv->L.e = L;
	rv->R.e = R;
	return rv;
	}

 static expr*
new_expr_if(Static *S, int o, expr *test, expr *tval, expr *fval)
{
	expr_if *e;

	if ((e = S->exprif_free))
		S->exprif_free = (expr_if*)e->test;
	else
		e = (expr_if*)exprmem(S, sizeof(expr_if));
	e->op = o;
	e->test = test;
	e->tval = tval;
	e->fval = fval;
	return (expr*)e;
	}

 static expr*
new_expr_n(Static *S, real t)
{
	expr_nv *en;

	if ((en = S->exprn_free))
		S->exprn_free = *(expr_nv**)en;
	else
		en = (expr_nv*)exprmem(S, sizeof(expr_nv));
	en->op = f_OPNUM;
	en->varno = numind(S, t);
	return (expr*)en;
	}

 static expr *
eread(EdRead *R)
{
	ASLTYPE *asl;
	Static *S;
	char *wd;
	expr *F, *L, *T, *arg, **args, **argse, **argss;
	fint L1;
	func_info *fi;
	int *at, i, j, k, numargs, op, symargs;
	int (*Xscanf)(EdRead*, const char*, ...);
	plterm *p;
	real *b;
	real r;
	tfinfo *tfi;

	S = (Static *)R->S;
	asl = S->asl;
	Xscanf = xscanf;
	switch(edag_peek(R)) {
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
			tfi = (tfinfo*)mem(sizeof(tfinfo));
			tfi->fi = fi;
			args = (expr**)exprmem(S, j*(sizeof(expr*)+sizeof(int)+1));
			argse = argss = args + j;
			tfi->at = at = (int*)argse;
			wd = (char*)(at + j);
			i = k = numargs = symargs = 0;
			for(; i < j; ++i) {
				arg = eread(R);
				if ((op = arg->op) == f_OPHOL || op == f_OPIFSYM) {
					*--argss = arg;
					at[i] = --symargs;
					}
				else {
					if (op == f_OPNUM)
						wd[numargs] = 1;
					else {
						wd[numargs] = 0;
						k++;
						}
					at[i] = numargs;
					args[numargs++] = arg;
					}
				}
			if ((symargs = argse - argss)) {
				while(--argse > argss) {
					arg = *argss;
					*argss++ = *argse;
					*argse = arg;
					}
				}
			if (symargs && !(fi->ftype & 1))
				fscream(R, fi->name, symargs, "symbolic ");
			if (k && k < numargs)
				S->diglen += numargs;
			if (asl->i.ra_max < numargs)
				asl->i.ra_max = numargs;
			if (asl->i.sa_max < symargs)
				asl->i.sa_max = symargs;
			tfi->ali = asl->i.nfinv++;
			S->atlen += tfi->n = numargs + symargs;
			tfi->nr = numargs;
			tfi->nd = k;
			tfi->wd = wd;
			return new_expr(S, f_OPFUNCALL, (expr*)tfi, (expr*)args);

		case 'h':
			return holread(R);
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
			return new_expr_n(S, r);

		case 'o':
			break;

		case 'v':
			if (Xscanf(R,"%d",&k) != 1 || k < 0)
				badline(R);
			if (k >= S->nvar0)
				k += S->nvinc;
			if (k > max_var)
				badline(R);
			return (expr *)(S->var_e + k);

		default:
			badline(R);
		}

	if (Xscanf(R, asl->i.opfmt, &k) != 1 || k < 0 || k >= N_OPS)
		badline(R);
	switch(optype[k]) {

		case 1:	/* unary */
			return new_expr(S, k, eread(R), 0);

		case 2:	/* binary */
			L = eread(R);
			return new_expr(S, k, L, eread(R));

		case 12: /* vararg (min, max) */
			i = -1;
			Xscanf(R, "%d", &i);
			if (i <= 0)
				badline(R);
 more_3:
			args = argse = (expr**)exprmem(S, i*sizeof(expr*));
			do *argse++ = eread(R); while(--i > 0);
			return new_expr(S, k, (expr*)args, (expr*)argse);

		case 4: /* piece-wise linear */
			i = -1;
			Xscanf(R, "%d", &i);
			if (i <= 1)
				badline(R);
			plterms++;
			j = 2*i - 1;
			p = (plterm *)mem(sizeof(plterm) + (j-1)*sizeof(real));
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
			if (b[-2] <= 0.)
				p->z = 2*i - 2;
			else {
				b = p->bs + 1;
				while(*b <= 0.)
					b += 2;
				p->z = (b - p->bs) - 1;
				}
			return new_expr(S, k, eread(R), (expr*)p);

		case 5: /* if */
			L = eread(R);
			T = eread(R);
			F = eread(R);
			return new_expr_if(S, k, L, T, F);

		case 11: /* OPCOUNT */
		case 6: /* sumlist */
			i = 0;
			Xscanf(R, "%d", &i);
			if (i <= 2 && (optype[k] == 6 || i < 1))
				badline(R);
			if (S->slmax < i)
				S->slmax = i;
			goto more_3;
			}
	badline(R);
	return 0;
	}

 static void
ogfree(Static *S, ograd *og)
{
	ograd *og1;

	if ((og1 = og)) {
		while(og1->next)
			og1 = og1->next;
		og1->next = freeog;
		freeog = og;
		}
	}

 static real
ogfree1(Static *S, ograd **ogp)
{
	ograd *og = *ogp;
	*ogp = og->next;
	og->next = freeog;
	freeog = og;
	return og->coef;
	}

 static ograd *
new_ograd(Static *S, ograd *next, int varno, real coef)
{
	ograd *rv;

	if ((rv = freeog))
		freeog = rv->next;
	else
		rv = (ograd *)mem_ASL(S->a, sizeof(ograd));
	rv->next = next;
	rv->varno = varno;
	rv->coef = coef;
	return rv;
	}

 static linarg *
lahash(Static *S, linarg *la)
{
	ASLTYPE *asl;
	Ulong x;
	int i, k, nnz, *ov, *ov1;
	linarg *la1, *la2, **lap, **lap1, **q, **q0, **qe;
	real *oc, *oc1;
	size_t mask;
	union {real r; uint L[2];} u;

	asl = S->asl;
	x = nnz = la->nnz;
	ov = la->ov;
	oc = la->oc;
	for(i = 0; i < nnz; ++i) {
		u.r = oc[i];
		x = (x << 1 | (x & 0x80000000) >> 31) ^
			  (ov[i]*101 + u.L[W0] + 257*u.L[W1]);
		}
	for(lap = &lthash[x & lthashmask]; (la1 = *lap); lap = &la1->hnext) {
		if (la1->nnz == nnz) {
			ov1 = la1->ov;
			oc1 = la1->oc;
			for(i = 0;; ++i) {
				if (i >= nnz)
					return la1;
				if (ov[i] != ov1[i] || oc[i] != oc1[i])
					break;
				}
			}
		}
	*lap = la;
	la->hnext = 0;
	if (++S->asl->P.nlttot > lthashmask) {
		mask = lthashmask;
		k = S->klthash;
		q = q0 = lthash;
		qe = q + mask;
		lthashmask = mask = (mask << 1) | 1;
		lap = lthash = (linarg**)new_mblkzap(asl, S->klthash = k + 1);
		while(q <= qe) {
			for(la2 = *q++; la2; la2 = la1) {
				la1 = la2->hnext;
				x = nnz = la2->nnz;
				ov = la2->ov;
				oc = la2->oc;
				for(i = 0; i < nnz; ++i) {
					u.r = oc[i];
					x = (x << 1 | (x & 0x80000000) >> 31) ^
						  (ov[i]*101 + u.L[W0] + 257*u.L[W1]);
					}
				lap1 = lap + (x & mask);
				la2->hnext = *lap1;
				*lap1 = la2;
				}
			}
		del_mblk(q0);
		}
	return la;
	}

 static range *
new_range(Static *S, range *r, range **rp)
{
	ASLTYPE *asl;
	int k, len, n, uilen;
	range *r1, *r2, *r3, **rp1, **rq, **rq0, **rqe;
	size_t mask;

	asl = S->asl;
	uilen = r->nv*sizeof(int);
	len = sizeof(range) + uilen;
	r1 = (range*)mem(len);
	r1->nintv = 0;
	n = r1->n = r->n;
	if (n < (r1->nv = r->nv)) {
		r1->uhlen = 0;
		S->nhop += n*n;
		}
#ifdef PSHVREAD
	else {
		S->uhlen += r1->uhlen = sizeof(uHeswork) + (r->n-1)*sizeof(Ogptrs);
		if (asl->P.rnmax < n)
			asl->P.rnmax = n;
		}
#endif
	r1->chksum = r->chksum;
	r1->refs = 0;
	r1->lasttermno = -1;
	r1->hnext = r1->hunext = 0;
	if (uilen)
		memcpy(r1->ui = (int*)(r1+1), r->ui, uilen);
	r1->lap = (linarg**)new_mblk(htcl(len = r->n*sizeof(linarg*)));
	memcpy(r1->lap, r->lap, len);
	r2 = r1->rlist.next = asl->P.rlist.next;
	r1->rlist.prev = r2->rlist.prev;
	r2->rlist.prev = asl->P.rlist.next = r1;
	*rp = r1;
	r1->irange = nrange;
	if (++nrange > rangehashmask) {
		mask = rangehashmask;
		k = S->krangehash;
		rq = rq0 = rangehash;
		rqe = rq + mask;
		rangehashmask = mask = (mask << 1) | 1;
		rp = rangehash = (range**)new_mblkzap(asl, S->krangehash = k + 1);
		while(rq <= rqe) {
			for(r2 = *rq++; r2; r2 = r3) {
				r3 = r2->hunext;
				rp1 = rp + (r2->chksum & mask);
				r2->hunext = *rp1;
				*rp1 = r2;
				}
			}
		del_mblk(rq0);
		}
	return r1;
	}

 static int
lacompar(const void *a, const void *b, void *v)
{
	int i, j, nza, nzb, *ova, *ovb;
	linarg *la, *lb;
	real *oca, *ocb, t;

	if (a == b)
		return 0;
	Not_Used(v);
	la = *(linarg **)a;
	lb = *(linarg **)b;
	ova = la->ov;
	ovb = lb->ov;
	oca = la->oc;
	ocb = lb->oc;
	nza = la->nnz;
	nzb = lb->nnz;
	for(j = 0;; ++j) {
		if (j >= nza) {
			if (j < nzb)
				return -1;
			break;
			}
		if (j >= nzb)
			return 1;
		if ((i = ova[j] - ovb[j]))
			return i;
		if ((t = oca[j] - ocb[j]))
			return t > 0. ? 1 : -1;
		}
	return 0;
	}

 static int
ndiff(range *r, range *r1)
{
	/* Return 0 if r and r1 have the same linargs; else return 1. */

	linarg **la, **la1, **la1e, **lae;

	la = r->lap;
	lae = la + r->n;
	la1 = r1->lap;
	la1e = la1 + r1->n;
	for(;;++la, ++la1) {
		if (la >= lae)
			return la1 < la1e;
		if (la1 >= la1e)
			return 1;
		if (lacompar((char*)la, (char *)la1, NULL))
			return 1;
		}
	return 0;
	}

 static range *
uhash(Static *S, range *r)
{
	ASLTYPE *asl;
	int len, n, nv, *ui, *uie;
	range *r1, **rp;
	size_t L;

	asl = S->asl;
	L = 0;
	nv = r->nv;
	ui = r->ui;
	uie = ui + nv;
	len = sizeof(int)*nv;
	while (ui < uie)
		L = 37*L + *ui++;
	r->chksum = L;
	ui = r->ui;
	rp = rangehash + (L & rangehashmask);
	n = r->n;
	if (asl->P.merge)
	    while((r1 = *rp)) {
		if (r1->nv == nv && r1->n == n && (len == 0 || !memcmp(ui, r1->ui, len))) {
			if (n == 0 || !memcmp(r->lap, r1->lap, n*sizeof(linarg*)))
				return *rp;
			if (!ndiff(r, r1))
				return *rp;
			}
		rp = &r1->hunext;
		}
	return new_range(S, r, rp);
	}

 static range *
rhash(Static *S, range *r, int addnew)
{
	ASLTYPE *asl;
	int len, n, *ov, *ove;
	linarg *la, **lae, **lap;
	range *r1, **rp;
	unsigned long L;

	asl = S->asl;
	L = 0;
	lap = r->lap;
	n = r->n;
	lae = lap + n;
	len = n * sizeof(linarg*);
	while(lap < lae) {
		L *= 37;
		la = *lap++;
		ov = la->ov;
		for(ove = ov + la->nnz; ov < ove; ++ov)
			L = 101*L + *ov;
		}
	r->chksum = L;
	lap = r->lap;
	rp = rangehash + (L & rangehashmask);
	if (asl->P.merge)
	    while((r1 = *rp)) {
		if (r1->n == n && !memcmp(lap, r1->lap, len))
			return r1;
		rp = &r1->hnext;
		}
	if (addnew)
		return new_range(S, r, rp);
	return *rp;
	}

 static int
compar(const void *a0, const void *b0, void *v)
{
	int a, b, c;
	Static *S = (Static *)v;

	if ((a = *(int*)a0) >= max_var) {
		a = zl[a-nv0x];
		if ((b = *(int*)b0) >= max_var) {
			if ((c = a - zl[b-nv0x]))
				return c;
			return *(int*)a0 - b;
			}
		if (a == b)
			return -1;
		}
	else if ((b = *(int*)b0) >= max_var) {
		b = zl[b-nv0x];
		if (a == b)
			return 1;
		}
	return a - b;
	}

 static void
zcsort(Static *S, int *c, int *ci, int i, int n, int p)
{
	int j;

	if (n < nzclim || p < 0)
		qsortv(ci, n, sizeof(int), compar, S);
	else for(j = 0; i < p; i++)
		if (c[i])
			ci[j++] = i;
	}

 static int
hscompar(const void *a0, const void *b0, void *v)
{
	int a, b, c;
	Static *S = (Static *)v;

	if ((a = *(int*)a0) >= Ncom) {
		a = zl[a];
		if ((b = *(int*)b0) >= Ncom) {
			if ((c = a - zl[b]))
				return c;
			return *(int*)a0 - b;
			}
		a -= nv0x;
		if (a == b)
			return -1;
		}
	else if ((b = *(int*)b0) >= Ncom) {
		b = zl[b] - nv0x;
		if (a == b)
			return 1;
		}
	return a - b;
	}

 static int
compar1(const void *a0, const void *b0, void *v)
{
	ASLTYPE *asl;
	int a, a1, b, b1, c, maxv, maxvsp;
	Static *S = (Static *)v;

	maxv = max_var;
	maxvsp = S->maxspvar;
	a = *(int*)a0;
	b = *(int*)b0;
	if (a < maxvsp && b < maxvsp) {
		if (a >= maxv) {
			a1 = S->dvspb[a-nv0x];
			if (b >= maxv) {
				b1 = S->dvspb[b-nv0x];
				if ((c = b1 - a1))
					return c;
				}
			else {
				asl = S->asl;
				if (a1 == b && cexps[b - S->dv0].splitvno)
					return 1;
				a = a1;
				}
			}
		else if (b >= maxv) {
			b1  = S->dvspb[b-nv0x];
			asl = S->asl;
			if (a == b1 && cexps[a - S->dv0].splitvno)
				return -1;
			b = b1;
			}
		}
	return b - a;
	}

 static ograd*
compress(Static *S, ograd *og, real *cp, int *comvar)
{
	ograd *og1;
	int i, ix, j, j1, k, nzc1, *zc1, *zci1;
	real c, *rnz, t;

	rnz = S->_rnz;
	c = og->varno < 0 ? ogfree1(S,&og) : 0.;
	ix = nzc1 = 0;
	zc1 = S->zc1;
	zci1 = S->zci1;
	for(og1 = og; og1; og1 = og1->next) {
		zc1[i = og1->varno] = 1;
		zci1[nzc1++] = i;
		rnz[i] = og1->coef;
		if (ix < i)
			ix = i;
		}
	if (ix < nv0x) {
		*cp = c;
		*comvar = 0;
		for(og1 = og; og1; og1 = og1->next)
			zc1[og1->varno] = 0;
		return og;
		}
	*comvar = 1;
	for(i = 0; i < nzc1; )
		if ((j = zci1[i]) < nv0x)
			i++;
		else {
			if (!zc[j]++)
				zci[nzc++] = j;
			t = rnz[j];
			j1 = j - nv0x;
			if ((og1 = (S->asl->P.dv + j1)->ll)) {
				if (og1->varno < 0) {
					c += t*og1->coef;
					og1 = og1->next;
					}
				for(; og1; og1 = og1->next)
					if (!zc1[k = og1->varno]++) {
						zci1[nzc1++] = k;
						rnz[k] = t*og1->coef;
						}
					else
						rnz[k] += t*og1->coef;
				}
			zc1[j] = 0;
			zci1[i] = zci1[--nzc1];
			}
	*cp = c;
	ogfree(S, og);
	og = 0;
	if (nzc1 > 0) {
		zcsort(S, zc1, zci1, 0, nzc1, max_var);
		while(nzc1 > 0) {
			i = zci1[--nzc1];
			zc1[i] = 0;
			if ((t = rnz[i])) {
				og = new_ograd(S, og, i, t);
				if (!zc[i]++)
					zci[nzc++] = i;
				}
			}
		}
	return og;
	}

 static linarg *
afree(Static *S, ograd *og, expr **ep)
{
	ASLTYPE *asl;
	int comvar, i, nnz, *ov, ov0[32];
	la_ref *refs;
	linarg *la, *la1, *rv;
	ograd *og1, *ogx;
	real c, *oc, oc0[32], s, t;
	size_t L;

	asl = S->asl;
	rv = 0;
	if (!og || !(og = compress(S,og,&c,&comvar)))
		goto done;

	if ((la = ltfree))
		ltfree = la->hnext;
	else {
		la = (linarg *)mem(sizeof(linarg));
		la->refs = 0;
		la->termno = 0;
		}

	s = og->coef;
	if (s < 0)
		s = -s;
	ogx = og1 = og;
	nnz = 1;
	while((og1 = og1->next)) {
		++nnz;
		t = og1->coef;
		if (t < 0)
			t = -t;
		if (s < t) {
			s = t;
			ogx = og1;
			}
		}

	la->nnz = nnz;
	if ((s = ogx->coef) != 1.)
		for(og1 = og; og1; og1 = og1->next)
			og1->coef /= s;
	lt_scale = s;
	if (nnz <= 32) {
		ov = ov0;
		oc = oc0;
		}
	else {
		oc = (real*)new_mblk(htcl(nnz*(sizeof(real) + sizeof(int))));
		ov = (int*)(oc + nnz);
		}
	la->ov = ov;
	la->oc = oc;
	for(og1 = og, i = 0; og1; og1 = og1->next, ++i) {
		ov[i] = og1->varno;
		oc[i] = og1->coef;
		}
	la1 = lahash(S, la);
	if (la1 == la) {
		L = nnz * (sizeof(real) + sizeof(int));
		memcpy(la->oc = (real*)mem(L), oc, nnz*sizeof(real));
		memcpy(la->ov = (int*)(la->oc + nnz), ov, nnz*sizeof(int));
		la->refs = 0;
		la->u.pv = 0;
		la->termno = Termno;
		la->tnext = tlist;
		tlist = la;
		la->lnext = asl->P.lalist;
		asl->P.lalist = la;
		la->hnext = 0;
		}
	else {
		/* !!? Eventually create or adjust la1->v to save cycles. */
		/* This requires possible adjustment for differing */
		/* constants and scale factors. */
		if (la1->termno == Termno)
			asl->P.ndupst++;	/* research statistic */
		else {
			free_laref(S, &la->refs);
			la1->termno = Termno;
			la1->tnext = tlist;
			tlist = la1;
			asl->P.ndupdt++;	/* research statistic */
			}
		la->hnext = ltfree;
		ltfree = la;
		la = la1;
		}
	ogfree(S, og);
	if (oc != oc0)
		del_mblk(oc);
	if (ep && (nnz > 1 || comvar)) {
		la->refs = refs = new_laref(S, la->refs);
		refs->ep = ep;
		refs->c = c;
		refs->scale = s;
		}
	if (nnz > 1)
		rv = la;
 done:
	return rv;
	}

 static ograd *
af_sum(Static *S, ograd *Log, ograd *Rog)
{
	ograd *oL, *oR, *oR1, *og, **ogp;

	oL = Log;
	oR = Rog;
	ogp = &og;
	for(;;) {
		if (!oL) {
			*ogp = oR;
			break;
			}
		if (!oR) {
			*ogp = oL;
			break;
			}
		if (oL->varno > oR->varno) {
			*ogp = oR;
			ogp = &oR->next;
			oR = *ogp;
			}
		else {
			if (oL->varno == oR->varno) {
				oL->coef += oR->coef;
				oR1 = oR->next;
				oR->next = 0;
				ogfree(S, oR);
				oR = oR1;
				if (oL->coef == 0.) {
					oR1 = oL->next;
					oL->next = 0;
					ogfree(S, oL);
					oL = oR1;
					continue;
					}
				}
			*ogp = oL;
			ogp = &oL->next;
			oL = *ogp;
			}
		}
	return og;
	}

 static cgrad *
cf_sum(Static *S, cgrad *Lcg, ograd *Rog)
{
	cgrad *cg, **cgp, *oL;
	ograd *oR, *oR1;

	oL = Lcg;
	oR = Rog;
	cgp = &cg;
	cg = 0;
	for(;;) {
		if (!oL) {
			assert(!oR);
			break;
			}
		if (!oR) {
			*cgp = oL;
			break;
			}
		assert(oL->varno <= oR->varno);
		if (oL->varno == oR->varno) {
			oL->coef += oR->coef;
			oR1 = oR->next;
			oR->next = 0;
			ogfree(S, oR);
			oR = oR1;
			}
		*cgp = oL;
		cgp = &oL->next;
		oL = *cgp;
		}
	return cg;
	}

 static void
sumlist_afree(Static *S, ograd *Laf, expr *e, expr **argso, expr **sls, expr **slscr)
{
	int nlin, nnl;
	expr *e1, *e2, **ep, **ep1;

	nlin = sls - slscr;
	ep = e->L.ep;
	nnl = argso - ep;
	switch(nlin) {
	  case 1:
		e1 = slscr[0];
		break;
	  case 2:
		e1 = new_expr(S, f_OPPLUS, slscr[0], slscr[1]);
		break;
	  default:
		ep1 = (expr**)exprmem(S, nlin*sizeof(expr*));
		e1 = new_expr(S, f_OPSUMLIST, (expr*)ep1, (expr*)(ep1+nlin));
		memcpy(ep1, slscr, nlin*sizeof(expr*));
		}
	switch(nnl) {
	  case 1:
		e2 = *ep;
		break;
	  case 2:
		e2 = new_expr(S, f_OPPLUS, ep[0], ep[1]);
		break;
	  default:
		ep1 = (expr**)exprmem(S, nnl*sizeof(expr*));
		e2 = new_expr(S, f_OPSUMLIST, (expr*)ep1, (expr*)(ep1+nnl));
		memcpy(ep1, ep, nnl*sizeof(expr*));
		}
	e->op = f_OPPLUS;
	e->L.e = e1;
	e->R.e = e2;
	afree(S, Laf, &e->L.e);
	}

 static void
sdvcite(Static *S, int k)
{
	range *r;
	split_ce *cs;
	linarg *la, **lap, **lape;

	cs = S->asl->P.Split_ce + (k - max_var);
	r = cs->r;
	lap = r->lap;
	lape = lap + r->n;
	while(lap < lape) {
		la = *lap++;
		if (la->termno != Termno) {
			free_laref(S, &la->refs);
			la->termno = Termno;
			la->tnext = tlist;
			tlist = la;
			}
		}
	}

 static ograd *
awalk(Static *S, expr *e)		/* return 0 if e is not linear */
{
	ASLTYPE *asl;
	expr **args, **argse, **argso, **sls, **slscr;
	expr_if *rvif;
	int k, k1, k2, kscr;
	linarg *la, **nl;
	ograd *Laf, *Raf, *rv, *taf;
	real t;
	tfinfo *tfi;

	switch(optype[k = e->op]) {

		case 1:	/* unary */
		case 3: /* exprc */
			Laf = awalk(S, e->L.e);
			if (k == f_OPUMINUS) {
				if ((taf = Laf)) {
					do taf->coef = -taf->coef;
						while((taf = taf->next));
					return Laf;
					}
				}
			if (Laf)
				afree(S, Laf, &e->L.e);
			break;

		case 2:	/* binary */
			Laf = awalk(S, e->L.e);
			Raf = awalk(S, e->R.e);
			if (Laf && Raf)
			 switch(k) {
			  case f_OPMINUS:
				taf = Raf;
				do taf->coef = -taf->coef;
					while((taf = taf->next));
				/* no break; */
			  case f_OPPLUS:
				return af_sum(S, Laf, Raf);

			  case f_OPMULT:
				if (Raf->varno < 0 && !Raf->next) {
 swap:
					taf = Laf;
					Laf = Raf;
					Raf = taf;
					goto scale;
					}
				if (Laf->varno < 0 && !Laf->next) {
					taf = Raf;
 scale:					t = Laf->coef;
					if (t == 0.) {
						ogfree(S, Raf);
						return Laf;
						}
					do taf->coef *= t;
						while((taf = taf->next));
					ogfree(S, Laf);
					return Raf;
					}
				break;

			  case f_OPDIV:
				if (Raf->varno < 0 && !Raf->next) {
					Raf->coef = 1. / Raf->coef;
					goto swap;
					}
			  }
			afree(S, Laf, &e->L.e);
			afree(S, Raf, &e->R.e);
			break;

		case 4: /* piece-wise linear */
			afree(S, awalk(S,e->L.e), &e->L.e);
			break;

		case 5: /* if */
			rvif = (expr_if *)e;
			/*afree(S, awalk(S,rvif->test), &rvif->test);*/
			/* The above may introduce "range" rows that */
			/* are not involved in any derivatives... */
			afree(S, awalk(S,rvif->tval), &rvif->tval);
			afree(S, awalk(S,rvif->fval), &rvif->fval);
			break;

		case 6: /* sumlist */
			args = e->L.ep;
			argse = e->R.ep;
			k1 = argse - args;
			while(!(Laf = awalk(S,*args++)))
				if (args >= argse)
					return 0;
			kscr = -1;
			asl = S->asl;
			if (!S->slscratchlev++) {
				if (k1 > S->slscratchmax) {
					if (S->slscratch)
						del_mblk(S->slscratch);
					k2 = htcl(k1*sizeof(expr*));
					S->slscratch = (expr**)new_mblk(k2);
					S->slscratchmax = 1 << k2;
					}
				slscr = S->slscratch;
				}
			else {
				kscr = htcl(k1*sizeof(expr*));
				slscr = (expr**)new_mblk(kscr);
				}
			sls = slscr;
			argso = args - 1;
			*sls++ = *argso;
			while(args < argse)
				if ((Raf = awalk(S,*args))) {
					*sls++ = *args++;
					Laf = af_sum(S, Laf, Raf);
					}
				else
					*argso++ = *args++;
			rv = 0;
			if (argso == e->L.ep) {
				rv = Laf;
				goto delscratch;
				}
			sumlist_afree(S, Laf, e, argso, sls, slscr);
 delscratch:
			--S->slscratchlev;
			if (kscr >= 0)
				del_mblk(slscr);
			return rv;

		case 7: /* function call */
			tfi = (tfinfo *)e->L.e;
			args = e->R.ep;
			for(argse = args + tfi->n; args < argse; args++)
				afree(S, awalk(S,*args), args);

		case 8: /* OPHOL (Hollerith) */
			break;

		case 9: /* OPNUM */
			return new_ograd(S, 0, -1, S->htvals_end[((expr_nv*)e)->varno]);

		case 10:/* OPVARVAL */
			asl = S->asl;
			if ((expr_nv *)e >= S->var_e && (expr_nv *)e < S->var_e + max_var)
				k = (expr_nv *)e - S->var_e;
			else if ((k = ((expr_vx*)e)->a0) < 0) {
				if ((la = ((expr_vx*)e)->la)
				  && la->termno != Termno) {
					la->termno = Termno;
					la->tnext = tlist;
					tlist = la;
					}
				return 0;
				}
			if ((k1 = (k - nv0x)) < 0) {
				if (!zc[k]++)
					zci[nzc++] = k;
				goto ogret;
				}
			cexps[k1].db = (derpblock*)1; /* indicate cexp used */
			if (k1 >= Ncom) {
				if (!zc[k = ((expr_vx*)e)->a1]++)
					zci[nzc++] = k;
				sdvcite(S,k);
				return 0;
				}
			if ((nl = (asl->P.dv + k1)->nl)) {
				if (!zc[k]++)
					zci[nzc++] = k;
				while((la = *nl++))
					if (la->termno != Termno) {
						la->termno = Termno;
						la->tnext = tlist;
						tlist = la;
						}
				return 0;
				}
			if (!zc[k]++)
				zci[nzc++] = k;
		ogret:
			return new_ograd(S, 0, k, 1.);

		case 11: /* OPCOUNT */
			args = e->L.ep;
			for(argse = e->R.ep; args < argse; args++)
				afree(S, awalk(S, *args), args);
			break;

		case 12: /* vararg (min, max) */
			args = e->L.ep;
			for(argse = e->R.ep; args < argse; ++args)
				afree(S, awalk(S, *args), args);
			break;

		default:
			scream(S->R, ASL_readerr_bug,
				"awalk: unexpected optype[%d] = %d\n",
				k, optype[k]);
			}
	return 0;
	}

 static void
cexp_upgrade(Static *S, int t)
{
	int k, n1, n2;
	cexp *ce;
	int *z;
	expr_nv **vp;
	split_ce *cs;
	ASLTYPE *asl = S->asl;

	n2 = t - Ncom;
	k = htcl(t*(sizeof(cexp) + sizeof(int) + sizeof(expr_nv*))
		+ n2*sizeof(split_ce));
	memset(ce = (cexp *)new_mblk(k), 0, n1 = sizeof(char*) << k);
	n1 = (n1 + Ncom*sizeof(split_ce))
		/ (sizeof(cexp) + sizeof(int) + sizeof(split_ce)
			+ sizeof(expr_nv*));
	n2 = n1 - Ncom;
	cs = (split_ce *)(ce + n1);
	vp = (expr_nv **)(cs + n2);
	z = (int *)(vp + n1);
	if (cexps) {
		if (nsce)
			memcpy(cs, asl->P.Split_ce, nsce*sizeof(split_ce));
		memcpy(ce, cexps, cexp_n*sizeof(cexp));
		memcpy(z, zl, cexp_n*sizeof(int));
		memcpy(vp, varp, cexp_n*sizeof(expr_nv*));
		del_mblk(cexps);
		}
	nsce = n2;
	asl->P.Split_ce = cs;
	cexps = ce;
	zl = z;
	cexp_k = k;
	cexp_n = n1;
	varp = vp;
	}

 static void
zc_upgrade(Static *S)
{
	int k, n, n0;
	int *z;
	ASLTYPE *asl = S->asl;

	k = htcl(sizeof(int)*(max_var1 + 1)) + 1;
	z = (int*)new_mblk(k);
	n = (sizeof(char*)/sizeof(int)) << (k-1);
	memset(z + n, 0, n*sizeof(int));
	if (zci) {
		n0 = (sizeof(char*)/sizeof(int)) << (kzc - 1);
		memcpy(z, zci, n0*sizeof(int));
		memcpy(z+n, zci+n0, n0*sizeof(int));
		del_mblk(zci);
		}
	kzc = k;
	zci = z;
	zc = z + n + 1;
	zc_lim = n;
	}

 static void
la_replace(Static *S, linarg *la)
{
	la_ref *r;
	expr_nv *v;
	expr_vx *vx;
	expr *eout;

	if (la->nnz > 1) {
		if ((vx = la->u.pv))
			v = &vx->v;
		else {
			la->u.pv = vx = (expr_vx *)exprmem(S, sizeof(expr_vx));
			vx->la = la;
			vx->a0 = vx->a1 = -1;
			v = &vx->v;
			v->op = f_OPVARVAL;
			}
		}
	else
		v = S->var_e + la->ov[0];
	for(r = la->refs; r; r = r->next) {
		efree(S, *r->ep);
		eout = (expr *)v;
		if (r->scale != 1.) {
			if (r->scale == -1.)
				eout = new_expr(S, f_OPUMINUS, eout, 0);
			else
				eout = new_expr(S, f_OPMULT, eout, new_expr_n(S, r->scale));
			}
		if (r->c)
			eout = new_expr(S, f_OPPLUS, eout, new_expr_n(S, r->c));
		*r->ep = eout;
		}
	free_laref(S, &la->refs);
	}

 static void
tlistgen(Static *S, ps_func *f)
{
	int *ci, *cie, i, *ov, *ove, t;
	linarg *la, **lap, **lape, *tl;
	range *r;
	psb_elem *b, *be;

	b = f->pi.b;
	be = f->pi.be;
	t = ++Termno;
	tl = 0;
	for(; b < be; b++) {
		if ((ci = b->ce)) {
			cie = ci + *ci;
			do {
				if (!zc[i = nv0x + *++ci])
					zci[nzc++] = i;
				}
				while(ci < cie);
			}
		r = b->U;
		lap = r->lap;
		lape = lap + r->n;
		while(lap < lape) {
			la = *lap++;
			if (la->termno != t) {
				la->termno = t;
				la->tnext = tl;
				tl = la;
				ov = la->ov;
				for(ove = ov + la->nnz; ov < ove; ++ov)
					if (!zc[*ov]++)
						zci[nzc++] = *ov;
				}
			}
		}
	tlist = tl;
	}

 static int
might_expand(Static *S, expr *e)
{
 top:
	switch(e->op) {
		case f_OPPLUS:
		case f_OPMINUS:
			if (e->R.e->op == f_OPNUM) {
				e = e->L.e;
				goto top;
				}
			if (e->L.e->op == f_OPNUM) {
				e = e->R.e;
				goto top;
				}
			/* no break */
		case f_OPSUMLIST:
			return 1;
		case f_OPUMINUS:
			e = e->L.e;
			goto top;
		case f_OPMULT:
			if (e->R.e->op == f_OPNUM) {
				e = e->L.e;
				goto top;
				}
			if (e->L.e->op == f_OPNUM) {
				e = e->R.e;
				goto top;
				}
			break;
		case f_OPVARVAL:
			if (((expr_nv*)e)->varno >= nv0x)
				return 1;
		}
	return 0;
	}

 static void
ce_split(Static *S, int i, ps_func *f)
{
	ASLTYPE *asl;
	cexp *c, *c1, *ce;
	expr **ep;
	expr_nv **vp, **vp0;
	expr_vx *vx;
	int j, j1, je, k, n;
	psb_elem *b;
	split_ce *cs;

	asl = S->asl;
	asl->P.ndvsplit++;
	n = f->pi.be - f->pi.b;
	j1 = asl->P.ndvspout;
	j = k = j1 + Ncom;
	asl->P.ndvspout += n;
	je = j + n;
	if (je > cexp_n)
		cexp_upgrade(S, je);
	c = c1 = cexps + j;
	b = f->pi.b;
	cs = asl->P.Split_ce + (j-Ncom);
	for(ce = c + n; c < ce; c++, b++, cs++) {
		c->o.e = b->o.e;
#ifdef PSHVREAD
		c->o.f = b->o.f;
		c->o.b = b->o.b;
#endif
		cs->r = b->U;
		cs->ce = b->ce;
		}
	c = cexps + i;
	c->splitvno = S->dv0 + j1;
	vp0 = vp = varp + j;
	j = max_var + nndv;
	je = j + n;
	nndv += n;
	j1 += max_var;
	while(j < je) {
		vx = (expr_vx *)mem(sizeof(expr_vx));
		vx->la = 0;
		*vp++ = (expr_nv*)vx;
#ifdef PSHVREAD
		(c1++)->varno =
#endif
			vx->v.varno = vx->a0 = j++;
		vx->a1 = j1++;
		vx->v.op = f_OPVARVAL;
		}
	if (n == 2)
		c->o.e = (int*)new_expr(S, f_OPPLUS, (expr*)vp0[0], (expr*)vp0[1]);
	else {
		ep = (expr**)new_mblk(htcl(n*sizeof(expr*)));
		memcpy(ep, vp0, n*sizeof(expr*));
		c->o.e = (int*)new_expr(S, f_OPSUMLIST, (expr*)ep, (expr*)(ep+n));
		}
	if ((max_var1 += n) >= zc_lim)
		zc_upgrade(S);

	/* adjust zl for use in compar() */
	i += nv0x;
	j = k + nv0x;
	n += k;
	do {
		zl[k++] = i;
		zci[nzc++] = j++;	/* for crefs */
		}
		while(k < n);
	}

 static void
linpart_augment(Static *S, cexp *c, ograd *og, ps_func *f)
{
	ASL *asl;
	int n;
	lincoef *lc, *lc0;
	linpart *L;
	ograd *og1;
	real t;

	asl = S->a;
	if (og->varno < 0) {
		if ((t = ogfree1(S, &og)))
			f->pi.b->o.e = (int*)new_expr(S, f_OPPLUS,
				(expr*)f->pi.b->o.e, new_expr_n(S, t));
		if (!og)
			return;
		}
	if ((L = c->lp)) {
		lc0 = L->lc;
		lc = lc0 + L->n;
		og1 = 0;
		do {
			--lc;
			og1 = new_ograd(S, og1, lc->varno, lc->coef);
			}
			while(lc > lc0);
		del_mblk(L);
		og = af_sum(S, og, og1);
		}
	og1 = og;
	for(n = 0; og; og = og->next)
		n++;
	c->lp = L = (linpart*)new_mblk(htcl(sizeof(linpart)+(n-1)*sizeof(lincoef)));
	L->n = n;
	lc = L->lc;
	for(og = og1; og; og = og->next, ++lc) {
		lc->varno = og->varno;
		lc->coef = og->coef;
		}
	ogfree(S, og1);
	}

 static ograd *
linterms(Static *S, linpart *L, real scale)
{
	lincoef *lc, *lc0;
	ograd *og;

	og = 0;
	lc0 = L->lc;
	lc = lc0 + L->n;
	while(lc > lc0) {
		--lc;
		og = new_ograd(S, og, lc->varno, scale*lc->coef);
		}
	return og;
	}

 static void
dvwalk(Static *S, int i)
{
	ASLTYPE *asl;
	cexp *c;
	dv_info *dvi;
	int n, nb, split, *vr, *z, *zi, *zie;
	linarg *tl, **tp;
	linpart *lp;
	ograd *og, *og0;
	ps_func f;

	asl = S->asl;
	Termno++;
	tlist = 0;
	c = cexps + i;
	split = zl[i] & 2;
	zl[i] = 0;
	asl->P.dvsp0[i] = Ncom + asl->P.ndvspout;
	if (split && might_expand(S, (expr*)c->o.e)) {
		asl->P.ndvspcand++;
		if ((og = cotermwalk(S, (expr**)&c->o.e, &f, 0, 0)) && f.pi.be > f.pi.b)
			linpart_augment(S, c, og, &f);
		nb = f.pi.be - f.pi.b;
		if (nb > 1) {
			asl->P.ndvsp[i] = zl[i] = nb;
			ce_split(S, i, &f);
			c = cexps + i;
			og = 0;
			}
		else if (nb == 1) {
			c->o.e = f.pi.b->o.e;
			og = 0;
			}
		tlistgen(S, &f);
		del_Elemtemp(S, last_psb_elem);
		}
	else
		og = awalk(S, (expr*)c->o.e);
	if ((og0 = og) && og->varno < 0 && !og->next)
		og0 = 0;
	if ((lp = c->lp))
		og = af_sum(S, og, linterms(S,lp,1.));
	dvi = asl->P.dv + i;
	if (og0) {
		dvi->ll = og;
		dvi->nl = 0;
		}
	else {
		if ((c->la = dvi->lt = afree(S, og, 0)))
			dvi->scale = lt_scale;
		for(n = 1, tl = tlist; tl; tl = tl->tnext)
			n++;
		dvi->ll = 0;
		dvi->nl = tp = (linarg **)mem(n*sizeof(linarg*));
		for(tl = tlist; tl; tl = tl->tnext)
			la_replace(S, *tp++ = tl);
		*tp = 0;
		}
	if ((n = nzc)) {
		nzc = 0;
		zi = zci;
		c->vref = vr = (int*)exprmem(S, (n+1)*sizeof(int));
		*vr++ = n;
		memcpy(vr, zi, n*sizeof(int));
		zie = zi + n;
		z = zc;
		do z[*zi++] = 0;
		while(zi < zie);
		}
	}

 static int
colindvref(Static *S, expr *e, int ndv)
{
	expr **a, **ae;
	int j, k, rv = 0;

 top:
	switch(e->op) {
		case f_OPPLUS:
		case f_OPMINUS:
			rv |= colindvref(S, e->R.e, ndv);
			/* no break */
		case f_OPUMINUS:
			e = e->L.e;
			goto top;
		case f_OPSUMLIST:
			a = e->L.ep;
			ae = e->R.ep;
			while(a < ae)
				rv |= colindvref(S, *a++, ndv);
			break;
		case f_OPVARVAL:
			k = ((expr_nv *)e)->varno;
			if ((k -= nv0x) < 0)
				break;
			if (zl[k]) {
				rv |= zl[k];
				break;
				}
			zl[k] = 1;
			if ((j = colindvref(S, (expr*)(S->cexps + k)->o.e, k))) {
				rv |= j;
				zl[k] |= j;
				}
			break;
		case f_OPMULT:
			if (e->R.e->op == f_OPNUM) {
				e = e->L.e;
				goto top;
				}
			if (e->L.e->op == f_OPNUM) {
				e = e->R.e;
				goto top;
				}
			/* no break */
		default:
			if (ndv >= 0)
				rv = zl[ndv] |= 2;
		}
	return rv;
	}

 static void
dv_walk(Static *S)
{
	ASLTYPE *asl;
	expr_nv *v;
	int *b, bi, *dvsp0, i, j, k, m, n, nv0, *ndvsp;

	asl = S->asl;
	m = asl->i.n_con0;
	if ((n = Ncom)) {
		nv0 = nv0x;
		for(i = 0; i < n_obj; i++)
			colindvref(S, (expr*)obj_de[i].o.e, -1);
		for(i = 0; i < m; i++)
			colindvref(S, (expr*)con_de[i].o.e, -1);
		larvlist = &v;
		k = S->combco;
		for(i = 0; i < k; ++i)
			dvwalk(S, i);
		S->dvsp1 = n + nndv;
		for(; i < n; ++i)
			dvwalk(S, i);
		S->dvsp2 = m = n + nndv;
		S->ndvbcosp = asl->P.ndvspout;
		if ((i = asl->P.ndvspout))
			nv0b = nv0 + Ncom + i;
		larvlist = 0;
		S->dvspb = b = (int*)exprmem(S, m*sizeof(int));
		dvsp0 = asl->P.dvsp0;
		ndvsp = asl->P.ndvsp;
		for(i = 0; i < n; ++i) {
			b[i] = bi = nv0 + i;
			if ((k = ndvsp[i]))
				for(j = dvsp0[i], k += j; j < k; ++j)
					b[j] = bi;
			}
		}
	}

 static int	/* adjust zci, return k s.t. zci[i] < nv0 for i < k */
nzcperm(Static *S)
{
	int i, j, k;

	k = nzc;
	for(i = 0; i < k; )
		if (zci[i] >= nv0x) {
			j = zci[--k];
			zci[k] = zci[i];
			zci[i] = j;
			}
		else
			i++;
	return k;
	}

 static void
sumlist_adj(ASLTYPE *asl, expr *e, expr *e1)
{
	int k, n;
	expr **ep, **ep0, **ep1;

	ep0 = e->L.ep;
	ep = e->R.ep;
	k = htcl((n = ep - ep0)*sizeof(expr*));
	if (n == (1 << k)) {
		ep1 = (expr**)new_mblk(k+1);
		memcpy(ep1, ep0, n*sizeof(expr*));
		del_mblk(ep0);
		e->L.ep = ep1;
		ep = ep1 + n;
		}
	*ep++ = e1;
	e->R.ep = ep;
	}

 static void
termwalk(Static *S, expr **ep, PSfind *p)
{
	ASLTYPE *asl;
	Elemtemp *bt;
	expr *e, *e1, *e2;
	int *cp, *cp0, *ce, *cee, i, j, j1, k, n, ncp, nzc2, *ov, *ove, *ui, *vr, *vre;
	linarg **lap, **lap1;
	linarg *tl;
	ps_func *f;
	psb_elem *b;
	range *r;
	size_t len;
	split_ce *cs;

	asl = S->asl;
	Termno++;
	tlist = 0;
	afree(S, awalk(S, *ep), ep);

	i = k = nzcperm(S);
	f = p->f;
	if (!larvlist)
		for(; i < nzc; i++) {
			if ((j = zci[i]) < max_var) {
				if ((j -= nv0x) >= 0 && (vr = cexps[j].vref))
					for(vre = vr + *vr; ++vr <= vre; )
						if ((j1 = *vr) >= nv0x && !zc[j1]++)
							zci[nzc++] = j1;
				}
			else {
				cs = asl->P.Split_ce + (j-max_var);
				if ((ce = cs->ce)) {
					cee = ce + *ce;
					while(ce++ < cee)
						if (!zc[j = *ce + nv0x]++)
							zci[nzc++] = j;
					}
				}
			}

	r = (range *)S->_rnz;	/* scratch */
	if ((ncp = nzc - k)) {
		zcsort(S, zc, zci+k, nv0x, ncp, -1 /*max_var*/);
		cp = cp0 = (int*)(r+1);
		i = k;
		*cp = ncp;
		do {
			zc[j = zci[i++]] = 0;
			*++cp = j - nv0x;
			} while(i < nzc);
		}
	else
		cp0 = 0;

	for(n = 0, tl = tlist; tl; tl = tl->tnext) {
		n++;
		ov = tl->ov;
		for(ove = ov + tl->nnz; ov < ove; ++ov)
			if (!zc[*ov]++)
				zci[k++] = *ov;
		}
	nzc2 = k;
	if (zc[-1])
		--nzc2;	/* ignore constant */
	if (n <= 0 && (*ep)->op == f_OPNUM)
		goto done;

	r->n = n;
	r->nv = nzc2;
	r->lap = lap1 = lap = (linarg**)
		new_mblk(htcl(n*sizeof(linarg*)));
	for(tl = tlist; tl; tl = tl->tnext)
		la_replace(S, *lap1++ = tl);
	if (n > 1)
		qsortv(lap, n, sizeof(linarg*), lacompar, NULL);
	r->ui = ui = cp0 ? cp+1 : (int*)(r+1);
	zcsort(S, zc, zci, 0, nzc2, nv0x);
	for(j = 0; j < nzc2; j++)
		*ui++ = zci[j];
	r = n >= nzc2 ? uhash(S,r) : rhash(S,r,1);
	del_mblk(lap);

	e1 = *ep;
	if (!r || (i = r->lasttermno) == -1 || r->lastgroupno != Groupno) {
		bt = p->b;
		if ((i = p->nb++) >= bt->nmax)
			upgrade_Elemtemp(S, bt);
		b = f->pi.be = f->pi.b + i;
		b->conno = Conno;
		b->termno = i;
		b->groupno = Groupno;
		if ((b->U = r)) {
			r->lasttermno = i;
			r->lastgroupno = Groupno;
			}
		b->o.e = (int*)e1;
		if (cp0) {
			cp = (int*)mem(len = (ncp+1)*sizeof(int));
			memcpy(cp, cp0, len);
			cp0 = cp;
			}
		b->ce = cp0;
		while(k > 0)
			zc[zci[--k]]= 0;
		}
	else {
		b = f->pi.b + i;
		while(k > 0)
			zc[zci[--k]] = 0;
		if ((cp = b->ce)
		 && (*cp != ncp || memcmp(cp+1, cp0+1, ncp*sizeof(int)))) {
			/* fix up f->ce */
			n = *cp;
			while(k < n)
				zc[zci[k++] = *++cp] = 1;
			for(n = 0; n < ncp; )
				if (!zc[j = cp0[++n]]++)
					zci[k++] = j;
			qsortv(zci, k, sizeof(int), hscompar, S);
			b->ce = cp = (int*)mem((k+1)*sizeof(int));
			*cp++ = k;
			n = 0;
			while(n < k)
				zc[*cp++ = zci[n++]] = 0;
			}
		e = (expr*)b->o.e;
		switch(e->op) {
		  case f_OPPLUS:
			ep = (expr**)new_mblk(2);
			ep[0] = e->L.e;
			ep[1] = e->R.e;
			ep[2] = e1;
			e->L.ep = ep;
			e->R.ep = ep + 3;
			e->op = f_OPSUMLIST;
			break;

		  case f_OPSUMLIST:
			sumlist_adj(asl, e, e1);
			break;
		  default:
			b->o.e = (int*)(e2 = new_expr(S, 0, e, e1));
			e2->op = f_OPPLUS;
		  }
		}
 done:
	nzc = 0;
	}

 static void
co_finish(ps_func *f)
{
	range *r;
	psb_elem *b, *be;
	psg_elem *g, *ge;

	b = f->pi.b;
	for(be = f->pi.be; b < be; b++)
		if ((r = b->U))
			r->lasttermno = -1;
	g = f->g;
	for(ge = f->ge; g < ge; g++)
		for(b = g->pi.b, be = g->pi.be; b < be; b++)
			if ((r = b->U))
				r->lasttermno = -1;
	}

 static expr *
ecopy(Static *S, expr *e)
{
	expr **a, **a1, **ae;
	int n, op;

	switch(op = e->op) {
		case f_OPPLUS:
		case f_OPMINUS:
			e = new_expr(S, op, ecopy(S, e->L.e), ecopy(S, e->R.e));
			break;
		case f_OPMULT:
			if (e->L.e->op == f_OPNUM) {
				e = new_expr(S, op,
					new_expr_n(S, S->htvals_end[((expr_nv*)e->L.e)->varno]),
					ecopy(S, e->R.e));
				break;
				}
			assert(e->R.e->op == f_OPNUM);
			e = new_expr(S, op,
				new_expr_n(S, S->htvals_end[((expr_nv*)e->R.e)->varno]),
				ecopy(S, e->L.e));
			break;
		case f_OPUMINUS:
			e = new_expr(S, op, ecopy(S, e->L.e), 0);
			break;
		case f_OPSUMLIST:
			a = e->L.ep;
			ae = e->R.ep;
			a1 = (expr**)new_mblk_ASL(S->a, htcl((n = ae-a)*sizeof(expr*)));
			e = new_expr(S, op, (expr*)a1, (expr*)(a1+n));
			while(a < ae)
				*a1++ = ecopy(S, *a++);
			break;
		case f_OPVARVAL:
			break;
#ifdef DEBUG
		default:
			fprintf(Stderr, "Impossible case in ecopy!\n");
			exit(1);
#endif
		}
	return e;
	}

 static int getgroup(Static *S, real scale, expr *e, PSfind *p);

 static ograd *
ltermwalk(Static *S, real scale, expr **ep, PSfind *p)
{
	ASLTYPE *asl;
	cexp *ce, *ce1, *ce1e;
	dv_info *dvi;
	expr **a, **ae, *e, *er;
	int k, k0, n;
	linpart *lp;
	ograd *og, *rv;
	real t;

	asl = S->asl;
	rv = 0;
 top:
	e = *ep;
	switch(e->op) {
		case f_OPPLUS:
			rv = af_sum(S, rv, ltermwalk(S, scale, &e->L.e, p));
			ep = &e->R.e;
			goto top;
		case f_OPMINUS:
			rv = af_sum(S, rv, ltermwalk(S, scale, &e->L.e, p));
			ep = &e->R.e;
			scale = -scale;
			goto top;
		case f_OPUMINUS:
			ep = &e->L.e;
			scale = -scale;
			goto top;
		case f_OPSUMLIST:
			a = e->L.ep;
			ae = e->R.ep;
			while(a < ae)
				rv = af_sum(S, rv, ltermwalk(S, scale, a++, p));
			break;
		case f_OPVARVAL:
			k = ((expr_nv*)e)->varno;
			if (k < nv0x) {
				rv = af_sum(S, rv, new_ograd(S,0,k,scale));
				break;
				}
			k0 = k;
			if ((k -= nv0x) >= Ncom)
				goto dflt;
			if (zl[k]) {
				ce = cexps + k;
				if (!ce->db) {
					ce->db = (derpblock*)1; /* indicate used */
					if ((n = asl->P.ndvsp[k])) {
						ce1 = cexps + asl->P.dvsp0[k];
						for(ce1e = ce1 + n; ce1 < ce1e; ++ ce1)
							ce1->db = (derpblock*)1;
						}
					}
				*ep = ecopy(S, (expr*)ce->o.e);
				if ((lp = ce->lp))
					rv = af_sum(S, rv, linterms(S,lp,scale));
				goto top;
				}
			if (!zc[k0]++)
				zci[nzc++] = k0;
			dvi = asl->P.dv + k;
			if (dvi->nl)
				goto dflt;
			og = new_ograd(S, 0, k0, scale);
			rv = af_sum(S, rv, og);
			break;
		case f_OPNUM:
			rv = af_sum(S, rv,
				new_ograd(S, 0, -1, scale*S->htvals_end[((expr_nv *)e)->varno]));
			break;
		case f_OPMULT:
			er = e->R.e;
			if (er->op == f_OPNUM) {
				*ep = e->L.e;
		case_opnum:
				t = S->htvals_end[((expr_nv*)er)->varno];
				if (t == 0.) {
					efree(S,*ep);
					*ep = er;
					}
				else {
					rv = af_sum(S, rv, ltermwalk(S, scale*t, ep, p));
					*(expr_nv**)er = S->exprn_free;
					S->exprn_free = (expr_nv*)er;
					}
				e->L.e = S->expr_free;
				S->expr_free = e;
				break;
				}
			er = e->L.e;
			if (er->op == f_OPNUM) {
				*ep = e->R.e;
				goto case_opnum;
				}
		default:
		dflt:
			if (p->g && getgroup(S, scale, e, p)) {
				if (p->nb + p->ng)
					++S->nhop;
				break;
				}
			if (scale != 1.)
				*ep = scale == -1.
					? new_expr(S, f_OPUMINUS, e, 0)
					: new_expr(S, f_OPMULT, e, new_expr_n(S, scale));
			termwalk(S, ep, p);
			asl->P.ns0++;
		}
	return rv;
	}

 static void
PSfind_init(Static *S, ps_func *f, PSfind *psf, int wantg)
{
	f->g = f->ge = 0;
	psf->nb = psf->ng = 0;
	psf->f = f;
	psf->b = new_Elemtemp(S, sizeof(psb_elem), (void**)&f->pi.b);
	f->pi.be = f->pi.b;
	if (wantg) {
		psf->g = new_Elemtemp(S, sizeof(psg_elem), (void**)&f->g);
		f->ge = f->g;
		}
	else {
		psf->g = 0;
		last_psb_elem = psf->b;
		}
	}

 static psb_elem *
psb_copy(psb_elem *b, psb_elem *b0, int n)
{
	range *r;
	psb_elem *be;

	memcpy(b, b0, n*sizeof(psb_elem));
	for(be = b + n; b < be; b++) {
		/* *b = *b0++; */ /* DEC Alpha gcc got this wrong */
		if (b->conno != -1 && (r = b->U)) {
			b->next = r->refs;
			r->refs = b;
			}
		}
	return be;
	}

 static int
getgroup(Static *S, real scale, expr *e, PSfind *p)
{
	ps_func *f, f1;
	PSfind p1;
	expr *e0, *e1, *e2;
	ograd *og, *og1;
	psb_elem *b, *be;
	psg_elem *g;
	Elemtemp *gt;
	linarg *la, **lap, **lape;
	lincoef *lc;
	linpart *L;
	range *U;
	int i, nb, nzc1, *ov, *ove, *zc1, *zci1;
	ASL *asl = S->a;

	e0 = 0;
	for(e1 = e; (optype[e1->op] & ~2) == 1; e1 = e1->L.e)
		e0 = e1;
	if (!e0 || !might_expand(S, e1))
		return 0;
	PSfind_init(S, &f1, &p1, 0);
	f = p->f;
	gt = p->g;
	if ((i = p->ng++) >= gt->nmax)
		upgrade_Elemtemp(S, gt);
	Groupno = p->ng;
	memset(g = f->g + i, 0, sizeof(psg_elem));
	g->g.e = e;
	g->ge.e = e0;
	g->scale = scale;
	if ((og = ltermwalk(S, 1., &e0->L.e, &p1)))
		og = compress(S, og, &g->g0, &i);
	for(e1 = e; e1 != e0; e1 = e2) {
		e2 = e1->L.e;
		e2->R.e = e1;	/* back pointer, formerly used in psderprop */
		}
	zc1 = S->zc1;
	zci1 = S->zci1;
	Groupno = nzc1 = 0;
	if ((og1 = og)) {
		for(i = 1; (og = og->next); i++);
		g->L = L = (linpart*)mem(sizeof(linpart) + (i-1)*sizeof(lincoef));
		L->n = i;
		lc = L->lc;
		for(og = og1;; ++lc) {
			zc1[zci1[nzc1++] = lc->varno = og->varno]= 1;
			lc->coef = og->coef;
			if (!(og = og->next))
				break;
			}
		ogfree(S, og1);
		}
	b = be = 0;
	if ((nb = p1.nb)) {
		b = (psb_elem*)mem(nb * sizeof(psb_elem));
		be = psb_copy(b, f1.pi.b, nb);
		}
	g->pi.b = b;
	g->pi.be = be;
	del_Elemtemp(S, p1.b);
	if (!b && !nzc1) {
		--p->ng;
		return 1;
		}
	for(; b < be; b++) {
		if (!(U = b->U))
			continue;
		lap = U->lap;
		lape = lap + U->n;
		while(lap < lape) {
			la = *lap++;
			ov = la->ov;
			for(ove = ov + la->nnz; ov < ove; ++ov)
				if (!zc1[*ov]++)
					zci1[nzc1++] = *ov;
			}
		}
	zcsort(S, zc1, zci1, 0, nzc1, nv0x);
	g->nov = nzc1;
	memcpy(g->ov = (int*)mem(nzc1*sizeof(int)), zci1, nzc1*sizeof(int));
	while(nzc1 > 0) {
		i = zci1[--nzc1];
		zc1[i] = 0;
		}
	return 1;
	}

 static ograd *
cotermwalk(Static *S, expr **ep, ps_func *f, int wantg, int omitdv)
{
	int comvar, nb, ng, x;
	ograd *og;
	real t;
	PSfind psf;
	psb_elem *b;
	psg_elem *g;

	PSfind_init(S, f, &psf, wantg);
	t = 0.;
	og = ltermwalk(S, 1., ep, &psf);
	if (omitdv && og)
		og = compress(S, og, &t, &comvar);
	b = f->pi.b;
	nb = psf.nb;
	g = f->g;
	ng = psf.ng;
	f->ge = g + ng;
	if (nb + ng == 0)
		*ep = new_expr_n(S, t);
	else if (t) {
		if (nb)
			b->o.e = (int*)new_expr(S, f_OPPLUS, (expr*)b->o.e, new_expr_n(S, t));
		else {
			nb = 1;
			memset(b, 0, sizeof(psb_elem));
			b->o.e = (int*)new_expr_n(S, t);
			}
		}
	f->pi.be = b + nb;
	co_finish(f);
	if (!larvlist) {
		if ((x = nb * sizeof(psb_elem) + ng * sizeof(psg_elem))) {
			ASL *asl = S->a;
			g = (psg_elem*)(x >= 256 ? M1alloc(x) : mem(x));
			b = (psb_elem*)(g + ng);
			if (nb)
				psb_copy(b, f->pi.b, nb);
			if (ng)
				memcpy(g, f->g, ng*sizeof(psg_elem));
			}
		del_Elemtemp(S, psf.b);
		if (wantg)
			del_Elemtemp(S, psf.g);
		if (!nb)
			b = 0;
		if (!ng)
			g = 0;
		f->pi.b = b;
		f->pi.be = b + nb;
		f->g = g;
		f->ge = g + ng;
		}
	return og;
	}

 static void
psfind(Static *S)
{
	ASLTYPE *asl;
	expr_vx *vx;
	int i, j, k, m, nla;
	linarg *la, **lap;
	ograd *og;
	ps_func *f, *f1;
	range *r, *r0;
#ifdef PSHVREAD
	int n;
#endif
	asl = S->asl;
	asl->P.krnmax = asl->P.rnmax = 0;
	asl->P.dvsp0 = (int*)M1zapalloc(2*(Ncom)*sizeof(int));
	asl->P.ndvsp = asl->P.dvsp0 + Ncom;
	m = asl->i.n_con0;
	dv_walk(S);
	asl->P.ndvspin = asl->P.ns0;
	asl->P.ns0 = 0;
	PSHV(n = 0);
	f = asl->P.ops;
	f1 = asl->P.cps;
	for(i = 0; i < n_obj; i++, f++) {
		Conno = -2 - i;
		og = cotermwalk(S, (expr**)&obj_de[i].o.e, f, wantOgroups, 1);
		Ograd[i] = af_sum(S, Ograd[i], og);
		PSHV(n += f->ge - f->g);
		}
	PSHV(asl->P.nobjgroups = n);
	PSHV(n = 0);
	for(i = 0; i < m; i++, f1++) {
		Conno = i;
		og = cotermwalk(S, (expr**)&con_de[i].o.e, f1, wantCgroups, 1);
		Cgrad[i] = cf_sum(S, Cgrad[i], og);
		PSHV(n += f1->ge - f1->g);
		}
	PSHV(asl->P.ncongroups = n);

	/* Compute varno values for linargs */

	S->larep = lap = lthash;
	lthash = 0;
	S->maxspvar = j = k = max_var + nndv;
	for(la = asl->P.lalist; la; la = la->lnext) {
		if ((vx = la->u.pv)) {
			vx->v.varno = j++;
			*lap++ = la;
			}
		}
	nndv += nla = j - k;
	if ((max_var1 += nla) >= zc_lim)
		zc_upgrade(S);

	/* Compute nintv value for ranges */
	/* and set lasttermno to first var. */

	r0 = (range*)&asl->P.rlist;
	for(r = asl->P.rlist.next; r != r0; r = r->rlist.next) {
		i = 0;
		if ((j = r->n) > 0) {
			lap = r->lap;
			r->lasttermno = (*lap)->ov[0];
			while(i < j) {
				la = lap[i];
				if (la->u.pv)
					i++;
				else {
					lap[i] = lap[--j];
					lap[j] = la;
					}
				}
			}
		r->nintv = i;
		}
	del_mblk(rangehash);
	rangehash = 0;
	}

/*******/

 static void
db_reset(Static *S)
{
	S->curdb.de = S->curdb.d0;
	S->curdb.next = 0;
	S->curdb.nxt = 0;
	}

 static void
new_dop(Static *S)
{
	int **dop, *op;

	op = nextop(S, 4);
	if (needalign(op,1)) {
		op[0] = OPGOTOFalign;
		op[1] = OPGOTOF;
		dop = (int**)&op[2];
		}
	else {
		op[0] = OPGOTOF;
		dop = (int**)&op[1];
		}
	S->dop = dop;
	opnext = (int*)&dop[1];
	}

 static derpblock*
db_save(Static *S, derpblock *dbsave)
{
	derpblock *db;

	if (dbsave->d0 < dbsave->de || dbsave->nxt) {
		db = (derpblock*)mem_ASL(S->a, sizeof(derpblock));
		*db = *dbsave;
		if (!db->nxt)
			note_firstdb(S, db);
		}
	else if (!(db = S->emptydb)) {
		S->emptydb = db = (derpblock*)mem_ASL(S->a, sizeof(derpblock));
		memset(db, 0, sizeof(derpblock));
		}
	return db;
	}

 static int
ewalk(Static *S, expr *e, uint *deriv, uint atop)
{
	ASL *asl;
	Condptrs *cp, cp0[2];
	Minmaxptrs *mmp, *mmp1;
	char *wd;
	derpblock dbsave, *db, *db1, *db2, **pdb;
	expr **args, **argse;
	expr_if *eif;
	expr_nv *eh;
	int a0, a1, a2, a3, **dop[2], *fh, i, i1, ica, j, j1, jca, k, k1, k2, k3;
	int n, nr, *op, *opg, opg0[8], *opg1, *opg2, *opg3, *oprev[2], opif;
	int **pf, **pf1, ***pg, **pop, rv;
	plterm *p, **pp;
	real *b, t;
	tfinfo **ptfi, *tfi;
	uint af, atop2, ia, ja, ka, kaf, *pia, *pja;

	asl = S->a;
	ica = jca = 0;
	switch(k = e->op) {
	  case f_OPFUNCALL:
		tfi = (tfinfo*)e->L.e;
		if ((n = tfi->n)) {
			memcpy(S->atc, tfi->at, n*sizeof(int));
			tfi->at = S->atc;
			S->atc += n;
			}
		args = (expr**)e->R.e;
		k = 4*n;
		opg = opg0;
		if (k > 8)
			opg = (int*)new_mblk(htcl(k*sizeof(int)));
		opg3 = opg + n;
		pja = (uint*)(opg3 + 2*n);
		wd = tfi->wd;
		ia = 0;
		if (!S->dop)
			new_dop(S);
		pia = 0;
		if (deriv) {
			*deriv = 0;
			pia = &ia;
			}
		nr = tfi->nr;
		for(i = k = 0; i < nr; ++i) {
			opg[i] = ewalk(S, args[i], pia, 0);
			if ((pja[i] = ia)) {
				k2 = 2*k++;
				opg3[k2] = opg[i];
				opg3[k2+1] = i;
				ia = 0;
				}
			else
				wd[i] = 1;
			}
		for(; i < n; ++i)
			opg[i] = ewalk(S, args[i], 0, 0);
		rv = wlast;
		tfi->wd = 0;
		k1 = tfi->nr;
		if ((tfi->nd = k)) {
			op = nextopp(S, n + 4 + sizeof(func_info*)/sizeof(int));
			if (needalign(op, 3)) {
				op[0] = OP_FUNCALL1align;
				ptfi = (tfinfo**)(op + 4);
				}
			else {
				op[0] = OP_FUNCALL1;
				ptfi = (tfinfo**)(op + 3);
				}
			++op;
			tfi->doff = (int*)mem(k*(k+2)*sizeof(int));
			memcpy(tfi->doff, opg3, 2*k*sizeof(int));
			tfi->fh = fh = tfi->doff + 2*k;
			wlast = rv + 4 + k1*(k1 + 3)/2;
			if (!(ja = atop))
				ja = new_a(S);
			*deriv = ja;
			k2 = rv + 4;
			af = S->afirst;
			for(i = 0; i < nr; ++i) {
				if ((ia = pja[i])) {
					new_derp(S, ia, ja, k2 + i);
					if (ia > af)
						free_a(S, ia);
					}
				}
			if (k < k1) {
				memcpy(tfi->wd = S->digc, wd, k1);
				S->digc += k1;
				}
			for(i = 0; i < k; ++i) {
				i1 = opg3[2*i+1];
				for(j = 0; j < k; ++j) {
					j1 = opg3[2*j+1];
					*fh++ = i1 >= j1 ? (i1*(i1+1)>>1) + j1
							 : (j1*(j1+1)>>1) + i1;
					}
				}
			}
		else {
			op = nextop(S, n + 3 + sizeof(func_info*)/sizeof(int));
			wlast = rv + 1;
			if (needalign(op,2)) {
				op[0] = OP_FUNCALL0align;
				ptfi = (tfinfo**)(op + 3);
				}
			else {
				op[0] = OP_FUNCALL0;
				ptfi = (tfinfo**)(op + 2);
				}
			}
		op[1] = rv;
		S->ptfi[tfi->ali] = *ptfi = tfi;
		op = (int*)(ptfi+1);
		memcpy(op, opg, n*sizeof(int));
		if (opg != opg0)
			del_mblk(opg);
		opnext = op + n;
		break;
	  case OPVARVAL:
		rv = 4*(k = ((expr_nv*)e)->varno);
		if (deriv)
			*deriv = k + 1;
		if (!zc[k]) {
			zc[k] = 1;
			zci[nzc++] = k;
			}
		break;
	  case FLOOR:
		k = n_FLOOR;
		goto ceil1;
	  case CEIL:
		k = n_CEIL;
 ceil1:
		i = ewalk(S, e->L.e, 0, 0);
		rv = wlast;
 uret0:
		if (!S->dop)
			new_dop(S);
		op = nextop(S, 3);
 uret0cp:
		op[0] = k;
		op[1] = rv;
		op[2] = i;
		wlast = rv + 1;
		break;
	  case ABS:
		i = ewalk(S, e->L.e, deriv, 0);
		rv = wlast;
		k = n_OPABS0;
		if (deriv && (ia = *deriv)) {
			if (ia > S->afirst)
				free_a(S, ia);
			if (!(ja = atop))
				ja = new_a(S);
			k2 = 5;
			goto uret1;
			}
		goto uret0;
	  case OPUMINUS:
		i = ewalk(S, e->L.e, deriv, 0);
		rv = wlast;
		k = OPUMINUS0;
		if (deriv && (ia = *deriv)) {
			if (ia > S->afirst)
				free_a(S, ia);
			if (!S->negone)
				S->negone = numind(S, -1.);
			if (!(ja = atop))
				ja = new_a(S);
			new_derp(S, ia, *deriv = ja, S->negone);
			k2 = 4;
			goto uret2;
			}
		goto uret0;
	  case OP_tanh:
		k = OPtanh0;
		goto unop;
	  case OP_tan:
		k = OP_tan0;
		goto unop;
	  case OP_sqrt:
		k = OP_sqrt0;
		goto unop;
	  case OP_sinh:
		k = OP_sinh0;
		goto unop;
	  case OP_sin:
		k = OP_sin0;
		goto unop;
	  case OP_log10:
		k = OP_log100;
		goto unop;
	  case OP_log:
		k = OP_log0;
		goto unop;
	  case OP_exp:
		k = OP_exp0;
		i = ewalk(S, e->L.e, deriv, 0);
		rv = wlast;
		if (deriv && (ia = *deriv)) {
			if (ia > S->afirst)
				free_a(S, ia);
			if (!(ja = atop))
				ja = new_a(S);
			new_derp(S, ia, *deriv = ja, rv);
			k2 = 4;
			goto uret2;
			}
		goto uret0;
	  case OP_cosh:
		k = OP_cosh0;
		goto unop;
	  case OP_cos:
		k = OP_cos0;
		goto unop;
	  case OP_atanh:
		k = OP_atanh0;
		goto unop;
	  case OP_atan:
		k = OP_atan0;
		goto unop;
	  case OP_asinh:
		k = OP_asinh0;
		goto unop;
	  case OP_asin:
		k = OP_asin0;
		goto unop;
	  case OP_acosh:
		k = OP_acosh0;
		goto unop;
	  case OP_acos:
		k = OP_acos0;
 unop:
		i = ewalk(S, e->L.e, deriv, 0);
		rv = wlast;
		if (!deriv || !(ia = *deriv))
			goto uret0;
		k2 = 6;
 uret1:
		if (ia > S->afirst)
			free_a(S, ia);
		if (!(ja = atop))
			ja = new_a(S);
		new_derp(S, ia, *deriv = ja, rv + 4);
 uret2:
		op = nextopp(S, 4);
 uret3:
		wlast = rv + k2;
		op[0] = k+1;
		op[2] = rv;
		op[3] = i;
		break;
	  case OPNOT:
		if (deriv) {
			*deriv = 0;
			deriv = 0;
			}
		k = n_OPNOT;
		goto unop;
	  case OPPLUS:
		k = OPPLUS0;
		jca = S->one;
		atop2 = atop;
 bret2:
		ica = S->one;
		if (!deriv) {
 mult0:
			i = ewalk(S, e->L.e, 0, 0);
			j = ewalk(S, e->R.e, 0, 0);
 bret0:
			if (!S->dop)
				new_dop(S);
			op = nextop(S, 4);
			op[0] = k;
			op[2] = i;
			op[3] = j;
			rv = op[1] = wlast;
			wlast = rv + 1;
			break;
			}
		k2 = 4;
		ia = ja = 0;
		*deriv = 0;
		if (atop)
			kaf = 0;
		else
			atop = kaf = new_a(S);
		i = ewalk(S, e->L.e, &ia, atop);
		j = ewalk(S, e->R.e, &ja, ia == atop2 ? 0 : atop2);
		if (!(ia | ja)) {
			if (kaf)
				free_a(S, atop);
			goto bret0;
			}
		af = S->afirst;
		if (ja > af && ja != atop)
			free_a(S, ja);
		if (ia > af && ia != atop)
			free_a(S, ia);
		if (!(ka = atop))
			ka = new_a(S);
		*deriv = ka;
		if (ka == ja && !atop2)
			new_derp(S, ja, ka, jca);
		if (ia && ia != ka)
			new_derp(S, ia, ka, ica);
		if (ja && ja != ka)
			new_derp(S, ja, ka, jca);
		goto bret1;
	  case OPMINUS:
		if (deriv && !(jca = S->negone))
			S->negone = jca = numind(S, -1.);
		k = OPMINUS0;
		atop2 = 0;
		goto bret2;
	  case OPMULT:
		k = OPMULT0;
		if (!deriv)
			goto mult0;
		ia = ja = 0;
		*deriv = 0;
		i = ewalk(S, e->L.e, &ia, 0);
		j = ewalk(S, e->R.e, &ja, 0);
		if (!(ia | ja))
			goto bret0;
		ica = j;
		jca = i;
		k2 = 4;
		rv = wlast;
		goto have_ijca;
	  case OPDIV:
		k = OPDIV0;
		if (!deriv)
			goto mult0;
		ia = ja = 0;
		*deriv = 0;
		i = ewalk(S, e->L.e, &ia, 0);
		j = ewalk(S, e->R.e, &ja, 0);
		if (!(ia | ja))
			goto bret0;
		rv = wlast;
		if (ia) {
			if (ja) {
				/* we use dL2 for dR2 */
				ica = rv + 4;
				jca = rv + 6;
				k2 = 8;
				}
			else {
				ica = rv + 4;
				jca = 0;
				k2 = 5;
				}
			}
		else {
			k2 = 6;
			jca = rv + 4;
			}
 have_ijca:
		af = S->afirst;
		if (ja > af)
			free_a(S, ja);
		if (ia > af)
			free_a(S, ia);
		if (!(ka = atop))
			ka = new_a(S);
		*deriv = ka;
		if (ka == ja)
			new_derp(S, ja, ka, jca);
		if (ia)
			new_derp(S, ia, ka, ica);
		if (ja && ja != ka)
			new_derp(S, ja, ka, jca);
 bret1:
		rv = wlast;
 bret:
		if (ia)
			++k;
		if (ja)
			k += 2;
 bret3:
		wlast = rv + k2;
		op = nextopp(S, 5);
		op[0] = k;
		op[2] = rv;
		op[3] = i;
		op[4] = j;
		break;
	  case OPREM:
		k = n_OPREM0;
		if (!deriv)
			goto mult0;
		ia = ja = 0;
		*deriv = 0;
		i = ewalk(S, e->L.e, &ia, 0);
		j = ewalk(S, e->R.e, &ja, 0);
		if (!(ia | ja))
			goto bret0;
		af = S->afirst;
		if (ja > af)
			free_a(S, ja);
		if (ia > af)
			free_a(S, ia);
		if (!(ka = atop))
			ka = new_a(S);
		*deriv = ka;
		rv = wlast;
		k2 = 4;
		if (ja) {
			k2 = 5;
			new_derp(S, ja, ka, rv+4);
			}
		if (ia && ia != ka)
			new_derp(S, ia, ka, S->one);
		goto bret;
	  case f_OP2POW:
		ia = 0;
		i = ewalk(S, e->L.e, deriv ? &ia : 0, 0);
		rv = wlast;
		if (ia) {
			op = nextopp(S, 4);
			op[0] = OP_2POW1;
			++op;
			wlast = rv + 5;
			if (ia > S->afirst)
				free_a(S, ia);
			if (!(ka = atop))
				ka = new_a(S);
			*deriv = ka;
			new_derp(S, ia, ka, rv+4);
			}
		else {
			if (!S->dop)
				new_dop(S);
			op = nextop(S, 3);
			op[0] = OP_2POW0;
			wlast = rv + 1;
			}
		op[1] = rv;
		op[2] = i;
		break;
	  case f_OP1POW:
	  case f_OPCPOW:
	  case OPPOW:
		if (!deriv) {
			k = n_OPPOW0;
			goto mult0;
			}
		ia = ja = 0;
		*deriv = 0;
		i = ewalk(S, e->L.e, &ia, 0);
		switch(k) {
		  case f_OPCPOW:
			ja = ia;
			j = i;
			ia = 0;
			i = ((exprc*)e)->c;
			k = n_OPPOW0;
			break;
		  case f_OP1POW:
			j = ((exprc*)e)->c;
			goto more_OP1POW;
		  default:
			j = 0;
			k = n_OPPOW0;
			if (e->R.e)
				j = ewalk(S, e->R.e, &ja, 0);
		  }
		if (!ja) {
 more_OP1POW:
			if (j < 0) {
				t = S->htvals_end[j];
				rv = wlast;
				if (!ia) {
					if (!S->dop)
						new_dop(S);
					op = nextop(S, 6);
					if (needalign(op,3)) {
						k = OPCPOW0align;
						b = (real*)&op[4];
						}
					else {
						k = nOPCPOW0;
						b = (real*)&op[3];
						}
					b[0] = t;
					opnext = (int*)&b[1];
					goto uret0cp;
					}
				if (ia > S->afirst)
					free_a(S, ia);
				if (!(ka = atop))
					ka = new_a(S);
				*deriv = ka;
				new_derp(S, ia, ka, rv+4);
				k2 = 6;
				if (t == floor(t)) {
					k = nOPPOW1i;
					goto bret3;
					}
				op = nextopp(S, 9);
				if (needalign(op,4)) {
					k = OPCPOW0align;
					b = (real*)&op[5];
					}
				else {
					k = nOPCPOW0;
					b = (real*)&op[4];
					}
				b[0] = t;
				opnext = (int*)&b[1];
				goto uret3;
				}
			if (!ia)
				goto bret0;
			}
 binop:
		k2 = 9;
		rv = wlast;
		if (ia) {
			ica = rv + 4;
			if (ja)
				jca = rv + 6;
			else {
				k2 = 6;
				jca = 0;
				}
			}
		else {
			k2 = 6;
			jca = rv + 4;
			}
		goto have_ijca;
	  case OPLESS:
		k = n_OPLESS0;
		if (!deriv)
			goto mult0;
		ia = ja = 0;
		*deriv = 0;
		dbsave = S->curdb;
		i = ewalk(S, e->L.e, &ia, 0);
		j = ewalk(S, e->R.e, &ja, 0);
		if (!(ia | ja))
			goto bret0;
		++S->ncond;
		af = S->afirst;
		if (ja > af)
			free_a(S, ja);
		if (ia > af)
			free_a(S, ia);
		if (!(ka = atop))
			ka = new_a(S);
		*deriv = ka;
		op = nextopp(S, 6 + 2*sizeof(derpblock*)/sizeof(int));
		if (needalign(op,5)) {
			pdb = (derpblock**)&op[6];
			k = OPLESS01align - 1;
			}
		else {
			pdb = (derpblock**)&op[5];
			opnext = (int*)&pdb[2];
			}
		rv = wlast;
		wlast = rv + 6;
		jca = 4;
		if (ia) {
			if (ja) {
				wlast = rv + 7;
				k += 3;
				jca = 5;
				}
			else
				++k;
			}
		else
			k += 2;
		op[0] = k;
		op[2] = rv;
		op[3] = i;
		op[4] = j;
		pdb[1] = db_save(S, &dbsave);
		if (ia)
			new_derp(S, ia, ka, rv + 4);
		if (ja)
			new_derp(S, ja, ka, rv + jca);
		db = pdb[0] = (derpblock*)mem(sizeof(derpblock));
		*db = S->curdb;
		S->curdb.next = 0;
		S->curdb.nxt = rv + jca + 1;
		S->curdb.de = S->curdb.d0;
		break;
	  case OPOR:
		k = n_OPOR;
		goto andor;
	  case OPAND:
		k = n_OPAND;
 andor:
		if (!S->dop)
			new_dop(S);
		i = ewalk(S, e->L.e, 0, 0);
		op = nextop(S, 4 + sizeof(void*)/sizeof(int));
		if (needalign(op,3)) {
			pop = (int**)&op[4];
			*op++ = k + OP_ANDalign - n_OPAND;
			}
		else {
			pop = (int**)&op[3];
			opnext = (int*)&pop[1];
			}
		op[0] = k;
		op[1] = rv = wlast++;
		op[2] = i;
		j = ewalk(S, e->R.e, 0, 0);
		if (rv != j) {
			opg = nextop(S, 3);
			opg[0] = OPCOPY0;
			opg[1] = rv;
			opg[2] = j;
			}
		pop[0] = opnext;
		break;
	  case LT:
		k = n_OPLT;
 bnoderiv:
		if (!S->dop)
			new_dop(S);
		if (deriv) {
			*deriv = 0;
			deriv = 0;
			}
		i = ewalk(S, e->L.e, 0, 0);
		j = ewalk(S, e->R.e, 0, 0);
		goto bret0;
	  case LE:
		k = n_OPLE;
		goto bnoderiv;
	  case EQ:
		k = OPEQ;
		goto bnoderiv;
	  case GE:
		k = n_OPGE;
		goto bnoderiv;
	  case GT:
		k = n_OPGT;
		goto bnoderiv;
	  case NE:
		k = n_OPNE;
		goto bnoderiv;
	  case OPatan201:
	  case OPatan210:
	  case OP_atan2:
		if (!deriv) {
			k = OP_atan20;
			goto mult0;
			}
		ia = ja = 0;
		*deriv = 0;
		i = ewalk(S, e->L.e, &ia, 0);
		switch(k) {
		  case OPatan201:
			ja = ia;
			j = i;
			ia = 0;
			i = ((exprc*)e)->c;
			break;
		  case OPatan210:
			j = ((exprc*)e)->c;
			break;
		  default:
			j = ewalk(S, e->R.e, &ja, 0);
		  }
		k = OP_atan20;
		if (!(ia | ja))
			goto bret0;
		goto binop;
	  case OPintDIV:
		k = OPintdiv;
		goto bnoderiv;
	  case OPprecision:
		k = OP_precision;
		goto bnoderiv;
	  case OPround:
		k = OP_round;
		goto bnoderiv;
	  case OPtrunc:
		k = OP_trunc;
		goto bnoderiv;
	  case OPATLEAST:
		k = n_OPATLEAST;
		goto bnoderiv;
	  case OPATMOST:
		k = n_OPATMOST;
		goto bnoderiv;
	  case OPEXACTLY:
		k = n_OPEXACTLY;
		goto bnoderiv;
	  case OPNOTATLEAST:
		k = n_OPNOTATLEAST;
		goto bnoderiv;
	  case OPNOTATMOST:
		k = n_OPNOTATMOST;
		goto bnoderiv;
	  case OPNOTEXACTLY:
		k = n_OPNOTEXACTLY;
		goto bnoderiv;
	  case OP_IFF:
		k = n_OP_IFF;
		goto bnoderiv;
	  case MINLIST:
		k = OPMINLIST0;
		goto maxminlist;
	  case MAXLIST:
		k = OPMAXLIST0;
 maxminlist:
		args  = e->L.ep;
		argse = e->R.ep;
		n = argse - args;
		j = htcl(n*(sizeof(int) + sizeof(Minmaxptrs) + sizeof(int**)));
		mmp = (Minmaxptrs*)new_mblk(j);
		pg = (int***)(mmp + n);
		opg = (int*)(pg + n);
		pf = 0;
		if (deriv) {
			if (!S->opprev) {
				op = nextop(S, 1);
				opnext = op;
				}
			oprev[0] = S->opprev;
			dop[0] = S->dop;
			op = nextop(S, 6);
			if (needalign(op,1)) {
				op[0] = OPGOTOF2align;
				op[1] = OPGOTOF2;
				pf = (int**)&op[2];
				}
			else {
				op[0] = OPGOTOF2;
				pf = (int**)&op[1];
				}
			dop[1] = pf;
			pf1 = pf + 1;
			pf[0] = pf[1] = 0;
			opnext = (int*)&pf[2];
			ica = S->one;
			*deriv = 0;
			dbsave = S->curdb;
			S->curdb.de = dbsave.d0;
			db = db2 = 0;
			if ((ka = atop))
				kaf = 0;
			else
				ka = kaf = new_a(S);
			af = S->afirst;
			for(j = 0; j < n; ++j) {
				ia = 0;
				mmp[j].f = 0;
				S->dop = &mmp[j].f;
				S->opprev = oprev[0];
				opg[j] = ewalk(S, args[j], &ia, ka);
				mmp[j].b = S->opprev;
				if (ia) {
					if (S->opprev == oprev[0]) {/* simple variable */
						mmp[j].b = opg1 = nextopp(S, 3);
						opg1[0] = OPVARREF;
						opg1[2] = opg[j];
						}
					if (pf1) {
						*pf1 = mmp[j].f;
						pf1 = 0;
						}
					if (ia != ka) {
						if (ia > af)
							free_a(S, ia);
						new_derp(S, ia, ka, ica);
						}
					if (!db) {
						if (!(db = db2))
							db = db2 = db_save(S, &dbsave);
						*deriv = ka;
						}
					mmp[j].db = db1 = (derpblock*)mem(sizeof(derpblock));
					*db1 = S->curdb;
					S->curdb.de = S->curdb.d0;
					if ((db1->next = db) != S->emptydb)
						db1->nxt = 1;
					opg1 = nextop(S, 6);
					if (needalign(opg1,1)) {
						opg1[0] = OPGOTOF2align;
						opg1[1] = OPGOTOF2;
						pop = (int**)&opg1[2];
						}
					else {
						opg1[0] = OPGOTOF2;
						pop = (int**)&opg1[1];
						}
					pop[0] = pop[1] = 0;
					pg[j] = pop;
					opnext = (int*)&pop[2];
					mmp[j].d = opg[j];
					}
				else {
					if (!db2)
						db2 = db_save(S, &dbsave);
					mmp[j].d = -1;
					mmp[j].db = db2;
					pg[j] = 0;
					}
				}
			if (db) {
				if (dop[0])
					*dop[0] = op;
				++S->ncond;
				k1 = n*sizeof(Minmaxptrs)/sizeof(int);
				S->dop = 0;
				S->opprev = oprev[0];
				pf[0] = op = nextopp(S, 7 + k1 + n);
				opg1 = op + 4;
				++k;
				if (needalign(op,4+n)) {
					++opg1;
					if (k == OPMAXLIST1)
						k = OPMAXLISTalign;
					else
						k = OPMINLISTalign;
					}
				op[2] = rv = wlast;
				wlast = rv + 6;
				mmp1 = (Minmaxptrs*)&opg1[n];
				opg2 = (int*)&mmp1[n];
				opg2[0] = OPGOTOMM;
				opg2[1] = rv;
				opnext = opg3 = opg2 + 2;
				op[0] = k;
				op[3] = n;
				pf = 0;
				for(j = 0; j < n; ++j) {
					if ((pop = pg[j])) {
						pop[0] = opg2;
						if (pf)
							*pf = mmp[j].f;
						pf = &pop[1];
						}
					else
						mmp[j].f = opg3;
					}
				*pf = op;
				memcpy(opg1, opg, n*sizeof(int));
				memcpy(mmp1, mmp, n*sizeof(Minmaxptrs));
				S->curdb.de = S->curdb.d0;
				S->curdb.next = 0;
				S->curdb.nxt = rv + 5;
				del_mblk(mmp);
				break;
				}
			if (kaf)
				free_a(S, kaf);
			if (needalign(op,1))
				op[0] = OPGOTOF2nalign;
			else
				op[0] = OPGOTOF2n;
			if (!(S->dop = dop[0])) {
				S->dop = pf;
				pf = 0;
				}
			}
		else {
			for(j = 0; j < n; ++j)
				opg[j] = ewalk(S, args[j], 0, 0);
			}
		op = nextop(S, n + 3);
		op[0] = k;
		op[1] = rv = wlast++;
		op[2] = n;
		memcpy(op + 3, opg, n*sizeof(int));
		del_mblk(mmp);
		if (pf)
			pf[0] = pf[1] = opnext;
		break;

	  case OPPLTERM: /* piece-wise linear */
		p = (plterm*)e->R.e;
		if (deriv)
			*deriv = 0;
		i = ewalk(S, e->L.e, deriv, 0);
		++plterms;
		rv = wlast;
		if (deriv && (ia = *deriv)) {
			wlast = rv + 5;
			if (ia > S->afirst)
				free_a(S, ia);
			if (!(ka = atop))
				ka = new_a(S);
			*deriv = ka;
			new_derp(S, ia, ka, rv+4);
			op = nextopp(S, 5 + sizeof(void*)/sizeof(int));
			if (needalign(op,4)) {
				k = OP_PLTERM1align;
				pp = (plterm**)&op[5];
				}
			else {
				k = n_OPPLTERM1;
				pp = (plterm**)&op[4];
				opnext = (int*)&pp[1];
				}
			op[0] = k;
			op[2] = rv;
			op[3] = i;
			}
		else {
			wlast = rv + 1;
			if (!S->dop)
				new_dop(S);
			op = nextop(S, 4 + sizeof(void*)/sizeof(int));
			if (needalign(op,3)) {
				k = OP_PLTERM0align;
				pp = (plterm**)&op[4];
				}
			else {
				k = n_OPPLTERM0;
				pp = (plterm**)&op[3];
				opnext = (int*)&pp[1];
				}
			op[0] = k;
			op[1] = rv;
			op[2] = i;
			}
		*pp = p;
		break;
	  case OPIFnl:
		eif = (expr_if*)e;
		k = ewalk(S, eif->test, 0, 0);
		oprev[0] = S->opprev;
		ia = ja = ka = 0;
		db = 0;
		dbsave = S->curdb;
		dop[0] = dop[1] = S->dop;
		if (deriv) {
			pia = &ia;
			pja = &ja;
			new_dop(S);
			dop[1] = S->dop;
			op = nextopp(S, 6 + (2*sizeof(Condptrs) + sizeof(int*))/sizeof(int));
			if (!(ka = atop))
				ka = new_a(S);
			if (needalign(op,5)) {
				cp = (Condptrs*)&op[6];
				opif = OPIF1align;
				}
			else {
				cp = (Condptrs*)&op[5];
				opif = nOPIF1;
				}
			opnext = (int*)((int**)&cp[2] + 1);
			}
		else {
			pia = pja = 0;
			op = nextop(S, 3 + 2*sizeof(int*)/sizeof(int));
			if (needalign(op,2))
				pop = (int**)&op[3];
			else
				pop = (int**)&op[2];
			opnext = (int*)&pop[2];
			cp = 0;
			opif = nOPIF0;
			}
		cp0[0].e = opnext;
		cp0[0].f = 0;
		S->dop = &cp0[0].f;
		S->opprev = oprev[0];
		a0 = wlast;
		i = ewalk(S, eif->tval, pia, ka);
		cp0[0].b = oprev[1] = S->opprev;
		a1 = wlast;
		S->opprev = oprev[0];
		af = S->afirst;
		pop = 0;
		cp0[0].bder = cp0[1].bder = -1;
		if (ia) {
			cp0[0].bder = i;
			opg = nextop(S, alignnum(9,10));
			if (needalign(opg,5))
				pop = (int**)&opg[6];
			else
				pop = (int**)&opg[5];
			opnext = (int*)&pop[2];
			if  (ia > af && ia != ka)
				free_a(S, ia);
			if (ia != ka)
				new_derp(S, ia, ka, S->one);
			cp0[0].db = (derpblock*)mem(sizeof(derpblock));
			*cp0[0].db = S->curdb;
			S->curdb.de = S->curdb.d0;
			S->curdb.next = db = db_save(S, &dbsave);
			if (db != S->emptydb)
				S->curdb.nxt = 1;
			}
		else {
			opg = nextop(S, alignnum(6,7));
			if (needalign(opg,4))
				pop = (int**)&opg[5];
			else
				pop = (int**)&opg[4];
			opnext = (int*)&pop[1];
			}
		wlast = a0;
		cp0[1].e = opnext;
		cp0[1].f = 0;
		S->dop = &cp0[1].f;
		S->opprev = oprev[0];
		j = ewalk(S, eif->fval, pja, 0);
		cp0[1].b = S->opprev;
		if (ja) {
			cp0[1].bder = j;
			if (ja > af)
				free_a(S, ja);
			if (ja != ka)
				new_derp(S, ja, ka, S->one);
			cp0[1].db = (derpblock*)mem(sizeof(derpblock));
			*cp0[1].db = S->curdb;
			if (!ia)
				cp0[0].db = db_save(S, &dbsave);
			}
		else if (!ia) {
			if (ka > af)
				free_a(S, ka);
			ka = 0;
			}
		else
			cp0[1].db = db;
		a3 = a2 = wlast;
		if (a1 > a2)
			a3 = a1;
		if (ia && ja)
			goto newrv;
		if (i >= af) {
			if (i == j)
				rv = i;
			else if (i < j) {
				if ( a1 <= j)
					rv = j;
				else
					goto newrv;
				}
			else if (j < (int)af || a2 <= i)
				rv = i;
			else
				goto newrv;
			}
		else if (j >= af)
			rv = j;
		else {
 newrv:
			rv = a3;
			a3 += 4;
			}
		if (rv != j) {
			if (ja) {
				cp0[1].b = opg3 = nextopp(S, 4);
				opg3[0] = ia ? OPCOPY1a : OPCOPY1;
				opg3[2] = rv;
				opg3[3] = j;
				}
			else {
				opg3 = nextop(S, 3);
				opg3[0] = OPCOPY0;
				opg3[1] = rv;
				opg3[2] = j;
				}
			}
		if (rv != i) {
			if (ia) {
				opg[0] = OPCOPY1;
				opg[1] = opg - oprev[1];
				opg[2] = rv;
				opg[3] = i;
				cp0[0].b = opg;
				opg += 4;
				}
			else {
				opg[0] = OPCOPY0;
				opg[1] = rv;
				opg[2] = i;
				opg += 3;
				}
			}
		if (ka) {
			++S->ncond;
			if (needalign(opg, 1)) {
				pop = (int**)&opg[2];
				opg[0] = OP_GOTOalign;
				if (ia && ja) {
					opg[0] = OPGOTO2align;
					pop[1] = cp0[1].f;
					}
				}
			else {
				pop = (int**)&opg[1];
				opg[0] = OP_GOTO;
				if (ia && ja) {
					opg[0] = OPGOTO2;
					pop[1] = cp0[1].f;
					}
				}
			pop[0] = opnext;
			op[0] = opif;
			op[2] = rv;
			op[3] = k;
			op[4] = a3;
			*deriv = ka;
			S->opprev = op;
			if (dop[0])
				*dop[0] = op;
			S->dop = 0;
			opg2 = opnext;
			if (!cp0[0].f)
				cp0[0].f = opg2;
			if (!cp0[1].f)
				cp0[1].f = opg2;
			if (ia)
				pop[1] = cp0[1].f;
			memcpy(cp, cp0, 2*sizeof(Condptrs));
			*(int**)&cp[2] = opg2;
			S->curdb.de = S->curdb.d0;
			S->curdb.next = 0;
			S->curdb.nxt = a3;
			a3 += alignnum(1,2);
			}
		else {
			if (needalign(opg, 1)) {
				opg[0] = OP_GOTOalign;
				pop = (int**)&opg[2];
				}
			else {
				opg[0] = OP_GOTO;
				pop = (int**)&opg[1];
				}
			*pop = opnext;
			if (needalign(op, 2)) {
				op[0] = OPIF0align;
				++op;
				}
			else
				op[0] = nOPIF0;
			pop = (int**)&op[2];
			op[1] = k;
			pop[0] = cp0[0].e;
			pop[1] = cp0[1].e;
			if (dop[0])
				*(S->dop = dop[0]) = 0;
			else
				S->dop = dop[1];
			S->opprev = oprev[0];
			}
		wlast = a3;
		break;

	  case OPIMPELSE:
		k1 = OPCOPY0;
		alignarg(k2 = OP_IMPELSE_align;)
		k3 = n_OPIMPELSE;
		goto moreifsym;

	  case OPIFSYM:
		k1 = OP_COPYSYM;
		alignarg(k2 = OPIF0align;)
		k3 = nOPIF0;
 moreifsym:
		eif = (expr_if*)e;
		if (!S->dop)
			new_dop(S);
		k = ewalk(S, eif->test, 0, 0);
		op = nextop(S, alignnum(4,7));
		if (needalign(op, 2)) {
			op[0] = k2;
			++op;
			pop = (int**)&op[2];
			}
		else {
			op[0] = k3;
			pop = (int**)&op[2];
			opnext = (int*)&pop[2];
			}
		op[1] = k;
		rv = a0 = wlast;
		pop[0] = opnext;
		i = ewalk(S, eif->tval, 0, 0);
		if (rv < i)
			rv = i;
		opg = nextop(S, alignnum(5,7));
		pop[1] = opnext;
		a1 = wlast;
		wlast = a0;
		j = ewalk(S, eif->fval, 0, 0);
		if (rv < j)
			rv = j;
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
		if (needalign(opg,1)) {
			opg[0] = OP_GOTOalign;
			pop = (int**)&opg[2];
			}
		else {
			opg[0] = OP_GOTO;
			pop = (int**)&opg[1];
			}
		pop[0] = opnext;
		break;

	  case OPCOUNT:
		k = n_OPCOUNT;
		goto more_orlist;
	  case OPNUMBEROF:
		k = n_OPNUMBEROF;
		goto more_orlist;
	  case OPNUMBEROFs:
		k = n_OPNUMBEROFs;
		goto more_orlist;
	  case OPALLDIFF:
		k = n_OPALLDIFF;
		goto more_orlist;
	  case OPSOMESAME:
		k = n_OPSOMESAME;
		goto more_orlist;
	  case ANDLIST:
		k = n_OPANDLIST;
		goto more_orlist;
	  case ORLIST:
		k = n_OPORLIST;
 more_orlist:
		if (!S->dop)
			new_dop(S);
		args  = e->L.ep;
		argse = e->R.ep;
		n = argse - args;
		opg = opg0;
		if (n > 8)
			opg = (int*)new_mblk(htcl(n*sizeof(int)));
		for(i = 0; i < n; ++i)
			opg[i] = ewalk(S, args[i], 0, 0);
		op = nextop(S, n + 3);
		op[0] = k;
		op[1] = rv = wlast++;
		op[2] = n;
		memcpy(op+3, opg, n*sizeof(int));
		if (opg != opg0)
			del_mblk(opg);
		break;

	  case OPSUMLIST:
		if (!deriv) {
			k = OPSUMLIST0;
			goto more_orlist;
			}
		args  = e->L.ep;
		argse = e->R.ep;
		n = argse - args;
		opg = opg0;
		if (n > 8)
			opg = (int*)new_mblk(htcl(n*sizeof(int)));
		if (atop)
			kaf = 0;
		else
			atop = kaf = new_a(S);
		atop2 = atop;
		af = S->afirst;
		ica = S->one;
		for(i = k2 = 0; i < n; ++i) {
			ia = 0;
			opg[i] = ewalk(S, args[i], &ia, atop2);
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
		rv = wlast;
		if (k2) {
			wlast = rv + 4;
			*deriv = atop;
			op = nextopp(S, n + 4);
			op[0] = OPSUMLIST1;
			++op;
			}
		else {
			wlast = rv + 1;
			if (kaf)
				free_a(S, atop);
			if (!S->dop)
				new_dop(S);
			op = nextop(S, n + 3);
			op[0] = OPSUMLIST0;
			}
		op[1] = rv;
		op[2] = n;
		memcpy(op+3, opg, n*sizeof(int));
		if (opg != opg0)
			del_mblk(opg);
		break;

	  case OPNUM:
		rv = ((expr_nv*)e)->varno;
		break;

	  case f_OPHOL:
		eh = (expr_nv*)e;
		ia = strlen(wd = (char*)(eh+1));
		op = nextop(S, k = 4 + ia/sizeof(int));
		op[0] = n_OPHOL;
		op[1] = rv = wlast++;
		op[2] = k;
		strcpy((char*)(op+3), wd);
		break;
 default:
		/*DEBUG*/fprintf(Stderr, "Bad k = %d in ewalk\n", k);
		/*DEBUG*/exit(1);
		rv = -1;
		}
	return rv;
	}

#ifndef PSHVREAD
#define Ops Ops1
#endif

 static void
ewalk1(Static *S, expr *e, uint *pia, Ops *o)
{
	int *on0, *op, nopblk, rv;
	uint atop, ia;
#ifdef PSHVREAD
	int k;
	o->b = o->f = 0;
	S->dop = S->dop0 = &o->f;
#endif
	S->opfirst = S->opprev = 0;
	S->opnext0 = on0 = opnext;
	nopblk = S->nopblks;
	opnext = oplast;
	atop = 0;
	if (pia) {
		atop = *pia;
		*pia = 0;
		}
	rv = ewalk(S, e, pia, atop);
#ifdef PSHVREAD
	o->b = S->opprev;
	if ((op = o->f)) {
		k = *op;
		if (k == nOPIF1 alignarg(|| k == OPIF1align)) {
			if (k == *o->b)
				*o->b = k += 3;
			else
				k += 1;
			*op = k;
			}
		else {
			op = o->b;
			k = *op;
			if (k == nOPIF1 alignarg(|| k == OPIF1align))
				*op = k + 2;
			}
		}
#endif
	if (atop && (ia = *pia) && ia != atop) {
		new_derp(S, ia, atop, S->one);
		*pia = atop;
		}
	if (rv < 0 && S->htvals_end[rv] == 0.) {
		o->e = 0;
		if (nopblk == S->nopblks)
			opnext = on0;
		}
	else {
		op = nextop(S, 2);
		op[0] = OPRET;
		op[1] = rv;
		o->e = S->opfirst;
		}
	}

 static void
ewalkL(Static *S, cde *cd)
{
	ASLTYPE *asl;
	cexp *ce, *cx0;
	int *ci, *cz, *dvsp0, i, j, m, nc, *ndvsp, nz;

	if ((nc = cd->afn)) {
		asl = S->asl;
		cx0 = cexps;
		dvsp0 = asl->P.dvsp0;
		ndvsp = asl->P.ndvsp;
		for(nc += i = cd->af1 - nv0x; i < nc; ++i) {
			if ((m = ndvsp[i])) {
				m += j = dvsp0[i];
				do {
					ce = cx0 + j;
					ewalk1(S, (expr*)ce->o.e, 0, &ce->o);
					} while(++j < m);
				}
			ce = cx0 + i;
			ewalk1(S, (expr*)ce->o.e, 0, &ce->o);
			}
		}
	ewalk1(S, (expr*)cd->o.e, 0, &cd->o);
	if ((nz = nzc)) {
		nzc = 0;
		cz = zc;
		ci = zci;
		while(nz > 0)
			cz[ci[--nz]] = 0;
		}
	}

 static void
derpcopy(Static *S, int kc, int k, int nz, int *ci, int nva)
{
	ASLTYPE *asl;
	cexp *ce, *ck, *cx0;
	derp *d0, *de, *first_d;
	derpblock cb, *db, **db0, **db1, *dbo, *dbx;
	int *dvsp0, i, j, je, mxv, na, ncom, *ndvsp, nn;

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
	dvsp0 = asl->P.dvsp0;
	ndvsp = asl->P.ndvsp;
	ncom = Ncom;
	mxv = max_var - nva;
	while(nz > k) {
		i = ci[--nz] - nva;
		if (i >= mxv)
			continue;
		ce = &cx0[i];
		for(db = ce->db; db; db = db->next) {
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
		if (!nva || i >= ncom)
			continue;
		j = dvsp0[i];
		je = j + ndvsp[i];
		while(j < je) {
			ce = &cx0[--je];
			for(db = ce->db; db; db = db->next) {
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
		}
	if (cb.d0 <= first_d) {
		first_d = (derp*)M1alloc(8190*sizeof(derp));
		cb.d0 = first_d + 8190;
		}
	ck = &cx0[kc];
	ck->db = dbo = (derpblock*)mem(sizeof(derpblock));
	dbo->nxt = 0;
	dbo->next = 0;
	dbo->d0 = first_d;
	i = wlast;
#ifdef PSHVREAD
	kc = ck->varno;
#else
	kc = asl->i.n_var0 + kc + S->nvinc;
#endif
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
	S->first_d = first_d;
	S->curdb.d0 = cb.d0;
	db_reset(S);
	}

 static int
lpcompar(const void *a, const void *b, void *v)
{
	Not_Used(v);
	return ((lincoef *)a)->varno - ((lincoef *)b)->varno;
	}

 static linpart *
linpt_read(EdRead *R, int nlin)
{
	ASL *asl;
	int i0, j, needsort;
	int (*Xscanf)(EdRead*, const char*, ...);
	lincoef *lc;
	linpart *rv;
	uint u;

	asl = R->asl;
	if ((j = nlin) <= 0)
		return 0;
	Xscanf = xscanf;
	rv = (linpart *)new_mblk(htcl(sizeof(linpart) + (nlin-1)*sizeof(lincoef)));
	rv->n = nlin;
	lc = rv->lc;
	i0 = needsort = 0;
	do {
		if (Xscanf(R, "%d %lf", &u, &lc->coef) != 2)
			badline(R);
		lc->varno = u;
		if (i0 > lc->varno)
			needsort++;
		i0 = lc->varno;
		++lc;
		}
		while(--j > 0);
	if (needsort)
		qsortv(&rv->lc, nlin, sizeof(lincoef), lpcompar, NULL);
	return rv;
	}

 static void
cexp_read(EdRead *R, int k, int nlin, int co)
{
	cde *c;
	cexp *ce;
	expr *e;
	Static *S = (Static *)R->S;
	ASLTYPE *asl = S->asl;

	ce = cexps + k - nv0x;
	if (co) {
		if (--co < n_con)
			c = con_de + co;
		else {
			co -= n_con;
			if (co < n_obj)
				c = obj_de + co;
			else
				c = lcon_de + (co - n_obj);
			}
		if (!c->afn) {
			c->afn = 1;
			c->af1 = k;
			}
		else
			++c->afn;
		}
	if (nlin)
		ce->lp = linpt_read(R, nlin);
	e = eread(R);
	ce->o.e = (int*)e;
	}

 static derpblock*
db_finish(Static *S, uint *pia)
{
	ASLTYPE *asl;
	derpblock *db;
	uint ia, ka;

	if (pia && (ia = *pia)) {
		if (S->curdb.d0 == S->curdb.de && !S->curdb.nxt) {
			/* ensure top derblock is nonempty */
			if (ia > S->afirst)
				ka = ia;
			else {
				ka = new_a(S);
				*pia = ka;
				}
			new_derp(S, ia, ka, S->one);	/* so the derpblock is nonempty */
			}
		}
	if (S->curdb.d0 != S->curdb.de || S->curdb.nxt) {
		asl = S->asl;
		db = (derpblock*)mem(sizeof(derpblock));
		*db = S->curdb;
		if (!db->nxt)
			note_firstdb(S, db); /* undo possible new_derprop in cexp_walk() */
		db_reset(S);
		}
	else {
		S->curdb.d0 = S->curdb.de;
		db = 0;
		}
	return db;
	}

 static int
intcomp(const void *a0, const void *b0, void *v)
{ return *(const int*)a0 - *(const int*)b0; }

 static void
cexp_walk(Static *S, int kc, int wantdb, int wantdc)
{
	ASLTYPE *asl;
	cexp *ce, *cx0;
	dv_info *dvi;
	expr_vx *vx;
	int *ci, *cz, fmax, i, j, k, mxv, n, nv, *o, *vr, *vre;
	linarg *la, **lap;
	lincoef *lc, *lce;
	linpart *lp;
	ograd *og;
	uint i1, ia, ncond0, nd, nd0, nz, nz0, *pia;

	/* for now, assume all common exprs need to be differentiated */

	asl = S->asl;
	cx0 = cexps;
	ce = cx0 + kc;
	nocopy = 0;
	nd0 = S->nderp;
	nz = nzc;
	cz = zc;
	ci = zci;
	if ((vr = ce->vref)) {
		ce->vref = 0;
		n = *vr;
		for(i = 1; i <= n; ++i) {
			if (!cz[j = vr[i]]) {
				cz[j] = 1;
				ci[nz++] = j;
				}
			}
		}
	ncond0 = S->ncond;
	nv = nv0x;
#ifdef PSHVREAD
	ia = ce->varno + 1;
#else
	ia = kc + max_var + 1;
#endif
	if ((lp = ce->lp)) {
		lc = lp->lc;
		if ((n = lp->n) == 1)
			new_derp(S, lc->varno + 1, ia, numind(S, lc->coef));
		else if (kc < Ncom) {
			dvi = asl->P.dv + kc;
			if (dvi->lt)
				new_derp(S, dvi->lt->u.pv->v.varno + 1, ia, numind(S, dvi->scale));
			if ((lap = dvi->nl)) {
				while((la = *lap)) {
					++lap;
					if ((vx = la->u.pv) && !cz[i1 = vx->v.varno]) {
						ci[nz++] = i1;
						cz[i1] = 1;
						}
					}
				}
			}
		lce = lc + n;
		do {
			if (!cz[i1 = lc->varno]) {
				ci[nz++] = i1;
				cz[i1] = 1;
				}
			}
			while(++lc < lce);
		}
	nzc = nz;
	if (ce->db) {
		pia = &ia;
		*pia = ia;
		ce->db = 0;
		}
	else {
		wantdb = 0;
		pia = 0;
		}
	ewalk1(S, (expr*)ce->o.e, pia, &ce->o);
	if (!ce->o.f && (o = ce->o.e) && o[0] == OPRET && (k = o[1]) >= 0) {
		o = nextop(S, 6);
		*o++ = OPRETB;
		o[0] = OPCOPY1;
		o[1] = 1;
		o[2] = 4*(S->dv0 + kc);
		o[3] = k;
		o[4] = OPRET;
		ce->o.b = ce->o.f = o;
		}
	if (wantdb)
		ce->dbf = ce->db = db_finish(S, pia);
	if ((nz0 = nz = nzc)) {
		nd = S->nderp - nd0;
		if ((fmax = maxfwd) <= 0 || !wantdc)
			goto f1_check;
		mxv = S->maxspvar;
		for(i = k = 0; i < nz; ++i) {
			if ((j = ci[i]) < nv) {
				if (++k > fmax)
					goto f1_check;
				}
			else if (j >= mxv)
				continue;
			else if ((vr = cx0[j-nv].vref)) {
				if (vr[2] == 2)
					i1 = vr[1];
				else
					i1 = vr[0];
				vr += 3;
				for(vre = vr + i1; vr < vre; ) {
					if (!cz[i1 = *vr++]) {
						cz[i1] = 1;
						ci[nz++] = i1;
						}
					}
				}
			}
		if (nd > 3*k && ncond0 == S->ncond && kc < Ncom) { /* funnelkind 2 */
			if (nz > 1)
				qsortv(ci, nz, sizeof(int), intcomp, S);
			ce->vref = vr = (int*)mem((nz+3)*sizeof(int));
			vr[0] = nz;
			vr[1] = k;
			vr[2] = 2;
			memcpy(vr+3, ci, nz*sizeof(int));
			derpcopy(S, kc, k, nz, ci, nv);
			}
		else {
 f1_check:
			if (nz0 > 1)
				qsortv(ci, nz0, sizeof(int), intcomp, S);
			for(k = 0; k < nz0 && ci[k] < nv; ++k);
			ce->vref = vr = (int*)mem((nz0+3)*sizeof(int));
			vr[0] = nz0;
			vr[1] = k;
			vr[2] = 0;
			memcpy(vr+3, ci, nz0*sizeof(int));
			if (wantdc && (nd > 3*nz0 || ncond0 != S->ncond) && kc < Ncom) {
				/* funnelkind 1 */
				vr[2] = 1;
				derpcopy(S, kc, nz0, nz0, ci, nv);
				}
			/* else no funnel */
			}
		while(nz > 0)
			cz[ci[--nz]] = 0;
		}
	else if (kc < Ncom && (og = asl->P.dv[kc].ll) && og->varno < 0)
		ce->o.e = 0;
	nzc = 0;
	S->firstdb1 = S->firstdb0;
	}

 static void
cg_zzap(ASLTYPE *asl)
{
	cgrad *cg, **cgp,**cgp1, **cgpe;

	cgp1 = Cgrad;
	cgpe = cgp1 + asl->i.n_con0;
	while(cgp1 < cgpe) {
		cgp = cgp1++;
		while((cg = *cgp))
			if (cg->coef)
				cgp = &cg->next;
			else
				*cgp = cg->next;
		}
	}

 static cexp **
funnels(Static *S, cexp *cx0, cexp *cx1, int ginit)
{
	ASLTYPE *asl;
	cexp *cx, **rv, **t;
	int n, *vr;

	asl = S->asl;
	for(n = 0, cx = cx0; cx < cx1; ++cx) {
		if ((vr = cx->vref) && vr[2])
			++n;
		}
	if (!n)
		return 0;
	asl->i.x0kindinit |= ginit;
	t = rv = (cexp**)mem((n+1)*sizeof(cde*));
	for(cx = cx0; cx < cx1; ++cx) {
		if ((vr = cx->vref) && vr[2])
			*t++ = cx;
		}
	*t = 0;
	return rv;
	}

 static void
adjust_compl_rhs(Static *S)
{
	ASLTYPE *asl;
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
	for(i = nlc; i < nc; ++i)
		if (Cvar[i]
		&& (o = C[i].o.e)
		&& o[0] == 0
		&& (j = o[1]) < 0
		&& (t = S->htvals_end[j]) != 0.) {
			t1 = t;
			if (L[j = stride*i] > negInfinity) {
				L[j] -= t;
				t1 = 0.;
				}
			if (U[j] < Infinity) {
				U[j] -= t;
				t1 = 0.;
				}
			if (t != t1)
				o[1] = numind(S, t1);
			}
	}

 static void
adjust(Static *S, int flags)
{
	ASLTYPE *asl = S->asl;
	cde *c, *ce;
	cexp *cx0, *cx1;
	int *c1, k;
	size_t L;
	void **vp, *vp1;

	if (S->slscratch)
		del_mblk(S->slscratch);
	if (asl->i.n_con0 && !allJ)
		cg_zzap(asl);
	if (S->k_seen) {
		if (!A_vals)
			goff_comp_ASL((ASL*)asl);
		else if (Fortran)
			colstart_inc_ASL((ASL*)asl);
		}
	if (n_cc > nlcc && nlc < n_con
	 && !(flags & ASL_no_linear_cc_rhs_adjust))
		adjust_compl_rhs(S);
	if ((k = ncom0)) {
		cx0 = cexps;
		if (comb) {
			cx1 = cx0 + comb;
			asl->I.dvfb = funnels(S, cx0, cx1, ASL_need_comba);
			cx0 = cx1;
			}
		if (comc) {
			cx1 = cx0 + comc;
			asl->I.dvfc = funnels(S, cx0, cx1, ASL_need_comca);
			cx0 = cx1;
			}
		if (como) {
			cx1 = cx0 + como;
			asl->I.dvfo = funnels(S, cx0, cx1, ASL_need_comoa);
			}
		}
	cx0 = cexps + k;
	if (comc1) {
		c1 = c_cexp1st;
		c = con_de;
		for(ce = c + n_con + n_lcon; c < ce; ++c) {
			*c1++ = k;
			k += c->afn;
			}
		*c1 = k;
		}
	if (como1) {
		c1 = o_cexp1st;
		c = obj_de;
		for(ce = c + n_obj; c < ce; ++c) {
			*c1++ = k;
			k += c->afn;
			}
		*c1 = k;
		}
	del_mblk(S->_rnz);
	asl->i.derplen = amax1 * sizeof(real);
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
	hes_setup(S);
	asl->i.wlen = wlast*sizeof(real);
	}

 static void
br_read(EdRead *R, int nc, real *L, real *U, int *Cvar, int nv)
{
	ASL *asl;
	int i, inc, j, k;
	int (*Xscanf)(EdRead*, const char*, ...);

	asl = R->asl;
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

 static expr *
aholread(EdRead *R)
{
	ASL *asl;
	FILE *nl;
	char *s1;
	expr_nv *rvh;
	int i, k;

	asl = R->asl;
	nl = R->nl;
	k = getc(nl);
	if (k < '1' || k > '9')
		badline(R);
	i = k - '0';
	while((k = getc(nl)) != ':') {
		if (k < '0' || k > '9')
			badline(R);
		i = 10*i + k - '0';
		}
	rvh = (expr_nv *)new_mblk(htcl(sizeof(expr_nv) + i + 1));
	rvh->op = f_OPHOL;
	rvh->varno = i;
	for(s1 = (char*)(rvh+1);;) {
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
	return (expr *)rvh;
	}

 static expr *
bholread(EdRead *R)
{
	ASL *asl;
	FILE *nl;
	char *s;
	expr_nv *rvh;
	int i;

	asl = R->asl;
	nl = R->nl;
	if (xscanf(R, "%d", &i) != 1)
		badline(R);
	rvh = (expr_nv *)new_mblk(htcl(sizeof(expr_nv) + i + 1));
	rvh->op = f_OPHOL;
	rvh->varno = i;
	s = (char*)(rvh+1);
	if (fread(s, i, 1, nl) != 1)
		badline(R);
	s[i] = 0;
	for(;;) switch(*s++) {
			case 0: goto break2; /* so we return at end of fcn */
			case '\n': R->Line++;
			}
 break2:
	return (expr *)rvh;
	}

 static int
qwalk(Static *S, expr *e)
{
	expr **args, **argse;
	int i, j, k;
	tfinfo *tfi;

	if (!e)
		return 0;
 top:
	switch(optype[k = e->op]) {

		case 1:	/* unary */
			switch(k) {
			  case f_OPCPOW:
				if (qwalk(S, e->R.e))
					return 3;
				return 0;
			  case f_OPUMINUS:
				e = e->L.e;
				goto top;
			  case f_OP2POW:
				i = qwalk(S, e->L.e);
				if (i < 2)
					return i << 1;
			  }
			return 3;

		case 2:	/* binary */
			switch(k) {
			 case f_OPPLUS:
			 case f_OPMINUS:
				i = qwalk(S, e->L.e);
				if (i == 3)
					break;
				j = qwalk(S, e->R.e);
				if (i < j)
					i = j;
				return i;
			 case f_OPDIV:
				if (qwalk(S, e->R.e))
					return 3;
				e = e->L.e;
				goto top;
			 case f_OPMULT:
				i = qwalk(S, e->L.e);
				if (i < 3) {
					i += qwalk(S, e->R.e);
					if (i < 3)
						return i;
					}
				}
			return 3;

		case 6: /* sumlist */
			j = 0;
			args = e->L.ep;
			argse = e->R.ep;
			while(args < argse) {
				i = qwalk(S, *args++);
				if (j < i) {
					j = i;
					if (j == 3)
						break;
					}
				}
			return j;

		case 7: /* function call */
			tfi = (tfinfo*)e->L.e;
			args = e->R.ep;
			argse = args + tfi->n;
			while(args < argse) {
				if (qwalk(S, *args++))
					return 3;
				}

		case 9: /* OPNUM */
			return 0;
		case 10:/* OPVARVAL */
			k = (expr_nv *)e - S->var_e;
			if (k < 0)
				goto k_adj;
			if (k < nv0x)
				return 1;
			if (k >= max_var) {
		k_adj:
				k = ((expr_vx*)e)->a1;
				return k < 0 ? 1 : S->asl->I.v_class[k-nv0x];
				}
			return S->asl->I.v_class[k - nv0x];
		}
	return 3;
	}

 static int
ewalkg(Static *S, expr *e, int k1)
{
	int k, *op, rv;

	rv = wlast;
	wlast = rv + 3;
	switch(e->op) {
	  case ABS:
		k = OPABS_g;
		wk_add(S->asl, &S->wk0, rv+2);
		break;
	  case OPUMINUS:
		k = OPUMINUS0;
		wk_add(S->asl, &S->wkm1, rv+1);
		wk_add(S->asl, &S->wk0, rv+2);
		break;
	  case OP_tanh:
		k = OPtanh_g;
		break;
	  case OP_tan:
		k = OPtan_g;
		break;
	  case OP_sqrt:
		k = OPsqrt_g;
		break;
	  case OP_sinh:
		k = OPsinh_g;
		break;
	  case OP_sin:
		k = OPsin_g;
		break;
	  case OP_log10:
		k = OPlog10_g;
		break;
	  case OP_log:
		k = OPlog_g;
		break;
	  case OP_exp:
		k = OPexp_g;
		break;
	  case OP_cosh:
		k = OPcosh_g;
		break;
	  case OP_cos:
		k = OPcos_g;
		break;
	  case OP_atanh:
		k = OPatanh_g;
		break;
	  case OP_atan:
		k = OPatan_g;
		break;
	  case OP_asinh:
		k = OPasinh_g;
		break;
	  case OP_asin:
		k = OPasin_g;
		break;
	  case OP_acosh:
		k = OPacosh_g;
		break;
	  case OP_acos:
		k = OPacos_g;
		break;
	  case OP1POW:
		k = OP1POW_g;
		op = nextop(S, 4);
		op[3] = ((exprc*)e)->c;
		goto have_op;
	  case f_OP2POW:
		k = OP2POW_g;
		wk_add(S->asl, &S->wk2, rv+2);
		break;
	  case f_OPCPOW:
		k = OPCPOW_g;
		op = nextop(S, 4);
		op[3] = k1;
		k1 = ((exprc*)e)->c;
		goto have_op;
	  case OPatan201:
		k = OPatan201_g;
		op = nextop(S, 4);
		op[3] = k1;
		k1 = ((exprc*)e)->c;
		goto have_op;
	  case OPatan210:
		k = OPatan210_g;
		op = nextop(S, 4);
		op[3] = ((exprc*)e)->c;
		goto have_op;
	  default:/*DEBUG*/
		fprintf(Stderr, "bad e->op = %d in ewalkg\n", e->op);
		exit(1);
		k = 0; /* not reached */
	  }
	op = nextop(S, 3);
 have_op:
	op[0] = k;
	op[1] = rv;
	op[2] = k1;
	return rv;
	}

 static void
Psb_walk(Static *S, Psbinfo *pi)
{
	ASLTYPE *asl;
	cexp *cx0, *cx;
	derpblock *db, **pdb;
	int *ce, *cee, *ci, *cz, i, ica, j, k, mxv;
	int ndb, nlp, nv, nz, tno, *vr;
	linarg *la, **lap, **larep, *tl;
	psb_elem *b, *be;
	uint af, ia, ja;

	asl = S->asl;
	ci = zci;
	cz = zc;
	cx0 = cexps;
	nv = nv0x;
	pdb = 0;
	mxv = S->maxspvar;
	b = pi->b;
	be = pi->be;
	ja = 0;
	if (be > b + 1) {
		ja = new_a(S);
		af = S->afirst;
		ica = S->one;
		for(; b < be; b++) {
			ia = 0;
			ewalk1(S, (expr*)b->o.e, &ia, &b->o);
			if (ia) {
				new_derp(S, ia, ja, ica);
				if (ia > af)
					free_a(S, ia);
				}
			}
		}
	else if (b < be)
		ewalk1(S, (expr*)b->o.e, &ja, &b->o);
	nz = nzc;
	nzc = 0;
	if (!(db = db_finish(S,&ja))) {
		while(nz > 0)
			cz[ci[--nz]] = 0;
		return;
		}
	for(i = j = 0; i < nz; ++i) {
		k = ci[i];
		if (k < nv)
			cz[k] = 0;
		else
			ci[j++] = k;
		}
	ndb = 1;
	if (!(nz = j))
		goto dbsave;
	larep = S->larep;
	nlp = 0;
	tl = 0;
	tno = ++Termno;
	for(i = 0; i < nz; ++i) {
		k = ci[i];
		if (k < mxv) {
			cx = &cx0[k-nv];
			if (cx->db)
				++ndb;
			if ((vr = cx->vref)) {
				ce = vr + 3;
				cee = ce + vr[0];
				ce += vr[1];
				while(ce < cee) {
					if (!cz[j = *ce++]) {
						ci[nz++] = j;
						cz[j] = 1;
						}
					if (j >= mxv) {
						la = larep[j-mxv];
						if (la->termno != tno) {
							++nlp;
							la->termno = tno;
							la->tnext = tl;
							tl = la;
							}
						}
					}
				}
			}
		else {
			la = larep[k-mxv];
			if (la->termno != tno) {
				++nlp;
				la->termno = tno;
				la->tnext = tl;
				tl = la;
				}
			}
		}
	if (nlp) {
		pi->lap = lap = (linarg**)mem((nlp+1)*sizeof(linarg*));
		do *lap++ = tl;
		   while((tl = tl->tnext));
		*lap = 0;
		}
 dbsave:
	pdb = pi->pdb = (derpblock**)mem((ndb+1)*sizeof(derpblock*));
	*pdb++ = db;
	mxv -= nv;
	if (nz) {
		pi->ce = ce = (int*)mem((nz+1)*sizeof(int));
		*ce++ = nz;
		if (nz > 1)
			qsortv(ci, nz, sizeof(int), compar1, S);
		for(i = 0; i < nz; ++i) {
			cz[k = ci[i]] = 0;
			ce[i] = k -= nv;
			if (k < mxv) {
				cx = &cx0[k];
				if ((*pdb = cx->db))
					++pdb;
				}
			}
		}
	*pdb = 0;
	}

 static int
dvuse(Static *S, ps_func *f)
{
	int *ce, *cee, k, n1, n2, n3, n4;
	psb_elem *b, *be;
	psg_elem *g, *ge;

	ge = f->ge;
	n1 = S->combco;
	n2 = Ncom;
	n3 = S->dvsp1;
	n4 = S->dvsp2;
	for(g = f->g; g < ge; ++g) {
		if ((b = g->pi.b))
			for(be = g->pi.be; b < be; ++b)
				if ((ce = b->ce))
					for(cee = ce + *ce; ce < cee;) {
						k = *++ce;
						if (k >= n1 && k < n2)
							return 1;
						if (k >= n3 && k < n4)
							return 1;
						}
		}
	return 0;
	}

 static int
co_walkloop(Static *S, ps_func *f, int n, cde *cd, char *c, ograd **o)
{
	ASLTYPE *asl;
	expr *e, *ee;
	int *dvsp0, i, j, k, k1, kx, m, nc, *ndvsp, nv, nz, *op, wantdb, *x;
	ps_func *fe;
	psb_elem *b, *be;
	psg_elem *g, *ge;

	kx = nz = 0;
	asl = S->asl;
	nv = nv0x;
	dvsp0 = asl->P.dvsp0;
	ndvsp = asl->P.ndvsp;
	for(fe = f + n; f < fe; ++cd, ++f) {
		if (amax1 < alast)
			amax1 = alast;
		alast = S->alast0;
		S->afree = S->afree0;
		if ((nc = cd->afn)) {
			wantdb = dvuse(S, f);
			for(nc += i = cd->af1 - nv; i < nc; ++i) {
				if ((m = ndvsp[i])) {
					m += j = dvsp0[i];
					do cexp_walk(S, j, wantdb, 0);
					   while(++j < m);
					}
				cexp_walk(S, i, wantdb, 0);
				}
			}
		if (c) {
			k = *o++ != 0;
			for(g = f->g, ge = f->ge; g < ge; g++) {
				if (g->g.e->op != f_OP2POW) {
					k = 3;
					goto have_c;
					}
				if (g->L)
					k = 2;
				for(b = g->pi.b, be = g->pi.be; b < be; b++) {
					if ((k1 = qwalk(S, (expr*)b->o.e)) > 1) {
						k = 3;
						goto have_c;
						}
					k = 2;
					}
				}
			for(b = f->pi.b, be = f->pi.be; b < be; b++) {
				if ((k1 = qwalk(S, (expr*)b->o.e)) > k) {
					k = k1;
					if (k == 3)
						goto have_c;
					}
				}
 have_c:		*c++ = (char)k;
			if (kx < k)
				kx = k;
			}
		if (!f->pi.b && !f->g) { /* linear */
			op = nextop(S, 2);
			op[0] = OPRET;
			op[1] = ((expr_nv*)cd->o.e)->varno;
			cd->o.e = op;
			continue;
			}
		if (f->pi.b)
			Psb_walk(S, &f->pi);
		for(g = f->g, ge = f->ge; g < ge; g++) {
			g->gm = wlast;
			wlast += g->nov + 3;
			for(e = g->g.e, ee = g->ge.e, k = 1; e != ee; ++k, e = e->L.e);
			if (S->gscrx < k)
				S->gscrx = k;
			g->nu = i = k;
			g->g.i = x = (int*)mem(k*sizeof(int));
			g->ge.i = x + k - 1;
			k1 = g->gm;
			S->dop = 0;
			S->opfirst = 0;
			S->opnext0 = opnext;
			opnext = oplast;
			opnext = g->o = nextop(S, 3);
			for(;;) {
				x[--i] = k1 = ewalkg(S, e, k1);
				if (!i)
					break;
				e = e->R.e;
				}
			op = nextop(S, 2);
			op[0] = OPRET;
			op[1] = k1;
			Psb_walk(S, &g->pi);
			}
		}
	return kx;
	}

 static void
do_ewalk(Static *S)
{
	ASLTYPE *asl;
	cde *Lde;
	cexp *c;
	char *v, *ve;
	expr_vx *vx;
	int *ci, cmb, *cz, *dvsp0, i, j, k, m, n, *ndvsp, nlogc, nz, wdc;
	linarg *la, **lap;
	size_t L;
#ifdef PSHVREAD
	Ihinfo *ihi, *ihi0, *ihi1;
	int h, ihdlim, ihdlim1, ihdmax, ihdmin;
	range *r, *r0;
#endif
/*
// After psfind() detects partially separable structure,
//
// qwalk() returns nonlinearity: 0 = const, 1 = linear, 2 = quadratic, 3 = general nonlinear
//
// to populate vclass().  Then
//
// cexp_walk() generates derps; calls ewalk1(), which turns pei.e into pei.i
// and computes derps.
//
// co_walkloop() uses qwalk() to populate c_class and o_class, then calls
// ewalk1() for psg_elem unary ops, then calls
// co_walk(), which calls ewalk1().

// Plan for later: do qwalk and ewalk1 of cexp1 common expressions after (or before) the
// read of the associated constraint or objective.  Do the constraint or objective's
// qwalk and ewalk1 after reading the constraint or objective.
*/

	asl = S->asl;
	psfind(S);
	amax1 = lasta = (S->alast0 = asl->i.maxvar = S->afirst = max_var + nndv) + 1;
	wlast = amax1 * sizeof(Varval)/sizeof(real);

	if (Ncom && asl->I.o_class) {
		i = max_var1 - nv0x;
		v = asl->I.v_class = (char*)M1zapalloc(i);
		ve = v + (S->maxspvar - nv0x);
		c = cexps;
		while(v < ve)
			*v++ = (char)((i = qwalk(S, (expr*)(c++)->o.e)) ? i : 1);
		if ((i = max_var1 - S->maxspvar)) {
			ve += i;
			do *v++ = 1; while(v < ve);
			}
		}

	if (S->asl->i.nfinv) {
		L = asl->i.nfinv*sizeof(tfinfo*) + S->atlen*sizeof(int) + S->diglen;
		S->ptfi = (tfinfo**)mem(L);
		asl->i.invd = (void**)S->ptfi;
		S->atc = (int*)(S->ptfi + asl->i.nfinv);
		S->digc = (char*)(S->atc + S->atlen);
		}
	ci = zci;
	cz = zc;
	nz = nzc;
	while(nz > 0)
		cz[ci[--nz]] = 0;
	nzc = 0;
	n = S->combco;
	dvsp0 = asl->P.dvsp0;
	ndvsp = asl->P.ndvsp;
	cmb = comb;
	wdc = 1;
	for(i = 0; i < n; i++) {
		if (i >= cmb)
			wdc = 0;
		if ((m = ndvsp[i])) {
			m += j = dvsp0[i];
			do {
				if (amax1 < alast)
					amax1 = alast;
				alast = S->alast0;
				S->afree = S->afree0;
				cexp_walk(S, j, 1, wdc);
				} while(++j < m);
			}
		if (amax1 < alast)
			amax1 = alast;
		alast = S->alast0;
		S->afree = S->afree0;
		cexp_walk(S, i, 1, wdc);
		}

	asl->I.c_class_max = co_walkloop(S, asl->P.cps, asl->i.n_con0, con_de,
					asl->I.c_class, (ograd**)Cgrad);
	asl->I.o_class_max = co_walkloop(S, asl->P.ops, n_obj, obj_de,
					asl->I.o_class, Ograd);
	if ((nlogc = n_lcon)) {
		Lde = lcon_de;
		for(i = 0; i < nlogc; ++i)
			ewalkL(S, Lde + i);
		}
	if (amax1 < alast)
		amax1 = alast;
#ifdef PSHVREAD
	r0 = (range*)&asl->P.rlist;
	if ((ihdlim = ihd_limit) > 0) {
		ihdmin = ihdlim1 = ihdlim + 1;
		n = ihdlim1*sizeof(Ihinfo);
		asl->P.ihi = ihi0 = (Ihinfo*)M1zapalloc(n);
		ihdmax = 0;
		h = wlast;
		for(r = asl->P.rlist.next; r != r0; r = r->rlist.next)
			if ((n = r->n) > 0) {
				if (n > r->nv)
					n = r->nv;
				if (n > ihdlim)
					n = ihdlim1;
				else {
					if (ihdmax < n)
						ihdmax = n;
					if (ihdmin > n)
						ihdmin = n;
					}
				ihi = ihi0 + n - 1;
				r->rlist.prev = ihi->r;
				ihi->r = r;
				ihi->nr++;
				}
		asl->P.ihdmax = ihdmax;
		asl->P.ihdmin = asl->P.ndhmax = ihdmin;
		ihi1 = ihi = ihi0 + ihdlim;
		ihi->ihd = ihdlim1;	/* sentinel */
		for(i = ihdlim; i > 0; --i) {
			if ((n = (--ihi)->nr)) {
				ihi->next = ihi1;
				ihi1 = ihi;
				ihi->ihd = i;
				k = (i*(i+1)) >> 1;
				for(r = ihi->r; r; r = r->rlist.prev) {
					r->hest = h;
					h += k;
					}
				}
			}
		asl->P.ihi1 = ihi1;
		wlast = h;
		}
#endif
	del_mblk(zci);
	zci = zc = 0;	/* for debugging */
	lap = &asl->P.lalist;
	k = 0;
	while((la = *lap)) {
		if ((vx = la->u.pv)) {
			lap = &la->lnext;
			la->u.v = vx->v.varno;
			la->termno = ++k;	/* formerly for psgcomp */
			}
		else {
			assert(la->nnz == 1);
			assert(la->ov[0] < nv0x);
			la->u.v = la->ov[0];
			*lap = la->lnext;
			la->termno = 0;
			}
		}
	del_mblk(S->larep);
	S->larep = 0;
#ifdef PSHVREAD
	asl->P.wkinit0 = wksave(asl, &S->wk0);
	asl->P.wkinit2 = wksave(asl, &S->wk2);
	asl->P.wkinitm1 = wksave(asl, &S->wkm1);
	asl->I.gscrx = S->gscrx;
	asl->I.nhop = S->nhop;
	asl->I.uhlen = S->uhlen;
	asl->P.nran = nrange;
	asl->P.zlsave = zl;
#endif
	}

 int
pfg_read_ASL(ASL *a, FILE *nl, int flags)
{
	ASLTYPE *asl;
	EdRead ER, *R;
	Jmp_buf JB;
	Static SS, *S;
	cde *Cde, *Lde, *Ode;
#ifdef PSHVREAD
	cexp *ce;
#endif
	cgrad *cg, **cgp;
	char *etofree, *etofree1, fname[128];
	expr_nv *var_ex;
	func_info *fi;
	int Ncom1, allG, i, i1, j, k, *ka, kseen, nc, nc0, nco, nlcon, nlin, nlogc, no;
	int nv, nv01, nvc, nvo, nvr, nvextra, nvx, *o, readall;
	int (*Xscanf)(EdRead*, const char*, ...);
	ograd *og, **ogp;
	real *oc, t;
	size_t L, LL[3], *kaz, nz;
	unsigned x;

	ASL_CHECK(a, asltype, who);
	flagsave_ASL(a, flags); /* includes allocation of LUv, LUrhs, A_vals or Cgrad, etc. */
	asl = (ASLTYPE*)a;
	S = S_init(&SS, asl);
	ed_reset(asl);
	i = comc ? ASL_need_concom: 0;
	if (como)
		i |= ASL_need_objcom;
	asl->i.x0kindinit = i;
	ncom0 = combc + como;
	SS.com1 = a->i.n_var0 + ncom0;
	asl->P.ncom = Ncom = (SS.combco = comb + comc + como) + comc1 + como1;
	nvr = n_var; /* nv for reading */
	nvextra = a->i.nsufext[ASL_Sufkind_var];
	nv01 = a->i.n_var0 + ncom0;
	SS.numht_mask = (1 << 10) - 1;
	nz = SS.numht_mask + 1;
	Ncom1 = ncom1 = comc1 + como1;
	SS.numht = (Numhash**)Malloc(LL[0] = nz*sizeof(Numhash*));
	memset(SS.numht, 0, nz*sizeof(Numhash*));
	SS.htvals = (real*)Malloc(LL[1] = nz*sizeof(real));
	SS.htvals_end = SS.htvals + nz;
	SS.numhtpth = &SS.numhtfirst;
	SS.numhtf0 = (void**)Malloc(LL[2] = nz*sizeof(void*) + Ncom1*sizeof(uint));
	asl->i.temp_rd_bytes = LL[0] + LL[1] + LL[2];
	SS.c1a = (uint*)(SS.numhtf0 + nz) - nv01; /* for convenience in eread() */
	*SS.numhtf0 = 0;
	SS.numhtf = (Numhash*)(SS.numhtf0 + 1);
	SS.numhtfend = SS.numhtf + ((nz-1)*sizeof(void*)/sizeof(Numhash));
	numind(S, 0.); /* for mpec_adj(): w[-1] = 0. */
	SS.one = numind(S, 1.);
	SS.two = numind(S, 2.);
	SS.negone = numind(S, -1.);
	SS.fmax = maxfwd;
	SS.R = R = EdReadInit_ASL(&ER, a, nl, S);
	SS.firstdb0 = SS.firstdb1 = (derpblock**)new_mblk(SS.kfirstdb = 3);
	SS.firstdbe = SS.firstdb0 + 8;
	if (flags & ASL_return_read_err) {
		a->i.err_jmp_ = &JB;
		i = setjmp(JB.jb);
		if (i) {
			a->i.err_jmp_ = 0;
			return i;
			}
		}
	if ((nlogc = a->i.n_lcon_) && !(flags & ASL_allow_CLP)) {
		if (a->i.err_jmp_)
			return ASL_readerr_CLP;
		sorry_CLP(R, "logical constraints");
		}
	if (!(flags & ASL_find_default_no_groups))
		flags |= ASL_findgroups;
	Xscanf = xscanf;
	readall = flags & ASL_keep_all_suffixes;
	PSHV(asl->P.pshv_g1 = 1;)
	k_Elemtemp = htcl(sizeof(Elemtemp));
	allJ = (flags & ASL_J_zerodrop) == 0;
	allG = (flags & ASL_G_zerodrop) == 0;
	wantCgroups = flags & ASL_findCgroups;
	wantOgroups = flags & ASL_findOgroups;
	asl->P.rlist.next = asl->P.rlist.prev = (range*)&asl->P.rlist;
	if (nfunc)
		func_add(a);
	if (binary_nl)
		holread = bholread;
	else
		holread = aholread;

	nc0 = n_con;
	nc = nc0 + a->i.nsufext[ASL_Sufkind_con];
	no = n_obj;
	SS.ncomo = Ncom + no;
	nvc = c_vars;
	nvo = o_vars;
	nco = nc + no + nlogc;
	if (no < 0 || nco <= 0)
		scream(R, ASL_readerr_corrupt,
			"pshvread: nc = %d, no = %d, nlogc = %d\n",
			nc0, no, nlogc);
	if (pi0)
		memset(pi0, 0, nc*sizeof(real));
	if (havepi0)
		memset(havepi0, 0, nc);
	asl->i.defvar0 = nv0x = nvx = nvr + nvextra;
	max_var1 = max_var = nv = nvr + Ncom + nvextra;
	combc = comb + comc;
	ncom0 = ncom_togo = combc + como;
	nzclim = nv >> 3;
	nv0b = nvr + comb;
	nv0c = nv0b + comc;
	x = (nc + no)*sizeof(ps_func)
		+ nco*sizeof(cde)
		+ no*sizeof(ograd *)
		+ nfunc*sizeof(func_info *)
		+ no;
	if (Ncom1)
		x += (nco + 1)*sizeof(int);
	SS.nvar0 = a->i.n_var0;
	if (!(SS.nvinc = (SS.dv0 = a->i.n_var_ + nvextra) - SS.nvar0)) {
		SS.dv0 = SS.nvar0;
		SS.nvar0 += ncom0 + Ncom1;
		}
	if (flags & ASL_find_co_class)
		x += nco;
	if (X0)
		memset(X0, 0, nvr*sizeof(real));
	if (havex0)
		memset(havex0, 0, nvr);
	asl->i.objconst = oc = (real*)M1zapalloc(x + no*sizeof(real));
	con_de = Cde = (cde *)(oc + no);
	lcon_de = Lde = Cde + nc;
	obj_de = Ode = Lde + nlogc;
	Ograd = (ograd **)(obj_de + no);
	funcs = (func_info **)(Ograd + no);
	asl->P.ops = (ps_func*)(funcs + nfunc);
	asl->P.cps = asl->P.ops + no;
	ka = (int *)(asl->P.cps + nc);
	c_cexp1st = l_cexp1st = o_cexp1st = 0;
	nlcon = n_lcon;
	if (Ncom1) {
		SS.c1s = ka;
		*ka = 0;
		if (comc1) {
			c_cexp1st = ka;
			ka += nc;
			if (nlcon) {
				l_cexp1st = ka;
				ka += nlcon;
				}
			}
		if (como1) {
			*ka = comc1;
			o_cexp1st = ka;
			ka += no;
			}
		++ka;
		}
	objtype = (char *)ka;
	k = 3;
	if ((j = n_var) > 7) {
		j >>= 3;
		do ++k; while(j >>= 1);
		}
	SS.klthash = SS.krangehash = k;
	SS._lthash = (linarg**)new_mblkzap(asl, k);
	SS._rangehash = (range**)new_mblkzap(asl, k);
	L = 1;
	SS._lthashmask = SS._rangehashmask = (L << k) - 1;
	SS._nrange = 0;
	L = nv*sizeof(real);
	if (L < sizeof(range))
		L = sizeof(range);
	k = htcl(L + (2*nv+1)*sizeof(int));
	SS._rnz = (real *)new_mblk(k);
	SS.zci1 = (int *)((char*)SS._rnz + L);
	SS.zc1 = S->zci1 + nv + 1; /* +1 for cp0 in termwalk() */
	memset(SS.zc1, 0, nv*sizeof(int));

	SS._conno = -1;	/* do not record refs for split defined vars */

	SS.var_e = (expr_nv*)exprmem(S, nv*sizeof(expr_nv));
	for(i = 0; i < nv; ++i) {
		SS.var_e[i].op = f_OPVARVAL;
		SS.var_e[i].varno = i;
		}
	zc_upgrade(&SS); /* initially allow nv+1 entries in zci and zc[-1] */
	if (Ncom) {
		var_ex = SS.var_e + a->i.n_var0;
		cexps = 0;
		asl->P.ndvspout = 0;
		memset(asl->P.dv = (dv_info*)mem(Ncom*sizeof(dv_info)),
			0, Ncom*sizeof(dv_info));
		cexp_upgrade(S, Ncom);
		j = asl->i.defvar0;
#ifdef PSHVREAD
		ce = cexps;
#endif
		for(k = 0; k < Ncom; k++) {
			varp[k] = &var_ex[k];
#ifdef PSHVREAD
			ce[k].varno = j + k;
#endif
			}
		}
	if (flags & ASL_find_co_class) {
		asl->I.o_class = objtype + no;
		asl->I.c_class = asl->I.o_class + no;
		}
	if (n_cc && !cvar)
		cvar = (int*)M1alloc(nc*sizeof(int));
	if (cvar)
		memset(cvar, 0, nc*sizeof(int));
	ka = 0;
	kaz = 0;
	kseen = 0;
	nz = 0;
	for(;;) {
		ER.can_end = 1;
		i = edag_peek(R);
		if (i == EOF) {
			fclose(nl);
			do_ewalk(S);
			for(etofree1 = SS.etofree; (etofree = etofree1); ) {
				etofree1 = *(char**)etofree;
				free(etofree);
				}
			/* Make amax long enough for nlc to handle */
			/* var_e[i].a for common variables i. */
			adjust(S, flags);
			nzjac = nz;
#ifdef PSHVREAD
			a->p.Conival = conpival_ew_ASL;
			a->p.Conival_nomap = conpival_nomap_ew_ASL;
			a->p.Congrd  = conpgrd_ew_ASL;
			a->p.Congrd_nomap  = conpgrd_nomap_ew_ASL;
			a->p.Objval = a->p.Objval_nomap = objpval_ew_ASL;
			a->p.Objgrd = a->p.Objgrd_nomap = objpgrd_ew_ASL;
			a->p.Conval = conpval_ew_ASL;
			a->p.EWalloc = ewalloc2_ASL;
			a->p.Jacval = jacpval_ew_ASL;
			a->p.Lconval= lconpval_ew_ASL;
			a->p.Hvcomp = a->p.Hvcomp_nomap = hvpcomp_ew_ASL;
			a->p.Hvcompd = hvpcompd_ew_ASL;
			a->p.Hvcompde = hvpcompde_ew_ASL;
			a->p.Hvcompe = hvpcompe_ew_ASL;
			a->p.Hvcomps = hvpcomps_ew_ASL;
			a->p.Hvcompse = hvpcompse_ew_ASL;
			a->p.Hvinit = a->p.Hvinit_nomap = hvpinit_ew_ASL;
			a->p.Hvinite = a->p.Hvinite_nomap = hvpinite_ew_ASL;
			a->p.Xknown = xp2known_ew_ASL;
			a->p.Duthes = a->p.Duthes_nomap = duthes_ew_ASL;
			a->p.Duthese = a->p.Duthese_nomap = duthese_ew_ASL;
			a->p.Fulhes = a->p.Fulhes_nomap = fullhes_ew_ASL;
			a->p.Fulhese = a->p.Fulhese_nomap = fullhese_ew_ASL;
			a->p.Sphes  = a->p.Sphes_nomap  = sphes_ew_ASL;
			a->p.Sphese = a->p.Sphese_nomap	= sphese_ew_ASL;
			a->p.Sphset = a->p.Sphset_nomap = sphes_setup_ew_ASL;
#else
			a->p.Xknown = xp1known_ew_ASL;
#endif
			return prob_adj_ASL(a);
			}
		ER.can_end = 0;
		k = -1;
		switch(i) {
			case 'C':
				Xscanf(R, "%d", &k);
				if (k < 0 || k >= nc0)
					badline(R);
				Cde[k].o.e = (int*)eread(R);
				break;
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
						mem(memadj(strlen(fname)+1)),
						fname);
					}
				if (!fi->funcp && !(fi->funcp = dynlink(fname)))
					scream(R, ASL_readerr_unavail,
						"function %s not available\n",
						fname);
				funcs[i] = fi;
				break;
			case 'L':
				Xscanf(R, "%d", &k);
				if (k < 0 || k >= nlogc)
					badline(R);
				Lde[k].o.e = (int*)eread(R);
				break;
			case 'V':
				if (Xscanf(R, "%d %d %d", &k, &nlin, &j) != 3)
					badline(R);
				if (k >= SS.nvar0)
					k += SS.nvinc;
				if (k < nvr || k >= nv)
					badline(R);
				cexp_read(R, k, nlin, j);
				break;
			case 'G':
				if (Xscanf(R, "%d %d", &j, &k) != 2
				|| j < 0 || j >= no || k <= 0 || k > nvo)
					badline(R);
				ogp = Ograd + j;
				while(k--) {
					if (Xscanf(R, "%d %lf", &i, &t) != 2)
						badline(R);
					if (allG || t) {
						*ogp = og = (ograd *)
							mem(sizeof(ograd));
						ogp = &og->next;
						og->varno = i;
						og->coef = t;
						}
					}
				*ogp = 0;
				break;
			case 'J':
				if (Xscanf(R, "%d %d", &j, &k) != 2
				|| j < 0 || j >= nc || k <= 0 || k > nvc)
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
				j = 0;
				while(k--) {
					if (kseen) {
						if (Xscanf(R, "%d %lf", &i, &t) != 2)
							badline(R);
						}
					else
						if (Xscanf(R, "%d %d %lf", &i, &j, &t) != 3)
							badline(R);
					*cgp = cg = (cgrad *)
						mem(sizeof(cgrad));
					cgp = &cg->next;
					cg->varno = i;
					cg->goff = j;
					cg->coef = t;
					}
				*cgp = 0;
				break;
			case 'O':
				if (Xscanf(R, "%d %d", &k, &j) != 2
				 || k < 0 || k >= no)
					badline(R);
				objtype[k] = j;
				o = Ode[k].o.e = (int*)eread(R);
				if (*o == f_OPNUM)
					oc[k] = SS.htvals_end[o[1]];
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
					x = nvx*sizeof(real);
					if (want_xpi0 & 4)
						x += nvx;
					X0 = (real *)M1zapalloc(x);
					if (want_xpi0 & 4)
						havex0 = (char*)(X0 + nvx);
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
#undef asl

#ifdef PSHVREAD

 static int
heswork(int *o, int *oe)
{
	Condptrs *cp;
	Minmaxptrs *mmp;
	int i, j, k, n, *z;
	plterm **pp;
	real *bs;
	tfinfo **ptfi, *tfi;

	n = 0;
	while(o != oe) {
#ifdef DEBUG
		/*DEBUG*/ if (++zorkheswork == zorkheswork1)
		/*DEGUG*/	printf("");
#endif
	    switch(*o) {

		case OPCOPY0:
			o += 3;
			break;

		case OPCOPY1:
		case OPCOPY1a:
			o += 4;
			break;

/*		case Hv_timesR: */
/*		case Hv_timesL: */
		case OPMULT01:
		case OPMULT10:
			n += 4;
			o += 5;
			break;

/*		case Hv_negate: */
/*		case Hv_plterm: */
		case OPUMINUS1:
			n += 4;
			o += 4;
			break;
#ifdef X64_bit_pointers
		case OP_PLTERM1align:
			pp = (plterm**)(o+5);
			goto more_plterm;
#endif
		case n_OPPLTERM1:
			pp = (plterm**)(o+5);
alignarg(more_plterm:)
			n += 4;
			o = (int*)&pp[1];
			break;

/*		case Hv_binaryR:*/
#ifdef X64_bit_pointers
		case OPCPOW1align:
		bs = (real*)&o[5];
		goto more_CPOW1;
#endif
		case nOPCPOW1:
		bs = (real*)&o[4];
alignarg(more_CPOW1:)
			n += 6;
			o = (int*)&bs[1];
			break;
		case OPDIV01:
		case OPDIV10:
		case nOPPOW01:
		case nOPPOW10:
		case nOPREM10:
			n += 6;
			o += 5;
			break;

		case nOPLESS01:
		case nOPLESS10:
			n += 6;
			o += 5 + 2*sizeof(derpblock*)/sizeof(int);
			break;

#ifdef X64_bit_pointers
		case OPLESS01align:
		case OPLESS10align:
			n += 6;
			o += 6 + 2*sizeof(derpblock*)/sizeof(int);
			break;
#endif
/*		case Hv_unary:	*/
		case n_OPABS1:
		case OP_acos1:
		case OP_acosh1:
		case OP_asin1:
		case OP_asinh1:
		case OP_atan1:
		case OP_atan201:
		case OP_atan210:
		case OP_atanh1:
		case OP_cos1:
		case OP_cosh1:
		case OP_exp1:
		case OP_log101:
		case OP_log1:
		case OP_sin1:
		case OP_sinh1:
		case OP_sqrt1:
		case OP_tan1:
		case OPtanh1:
		case OP_2POW1:
			n += 6;
			o += 4;
			break;

		case nOPPOW1i:
			n += 6;
			o += 5;
			break;

/*		case Hv_timesLR: */
		case OPMULT2:
			n += 10;
			o += 5;
			break;

/*		case Hv_binaryLR: */
		case OPDIV2:
		case n_OPPOW2:
		case OP_atan22:
			n += 14;
			o += 5;
			break;

/*		case Hv_vararg: */
#ifdef X64_bit_pointers
		case OPMINLISTalign:
		case OPMAXLISTalign:
			z = o + 5;
			goto more_minmax;
#endif
		case OPMINLIST1:
		case OPMAXLIST1:
			z = o + 4;
alignarg(more_minmax:)
			n = o[3];
			mmp = (Minmaxptrs*)&z[n];
			o = (int*)&mmp[n];
			for(i = k = 0; i < n; ++i) {
				if (mmp[i].f != o) {
					j = heswork(mmp[i].f, mmp[i].b);
					if (k < j)
						k = j;
					}
				}
			n = k + 2;
			break;

/*		case Hv_if: */
#ifdef X64_bit_pointers
		case OPIF1align:
		case OPIF11align:
		case OPIF12align:
		case OPIF13align:
			cp = (Condptrs*)&o[6];
			goto more_OPIF;
#endif
		case nOPIF1:
		case nOPIF11:
		case nOPIF12:
		case nOPIF13:
			cp = (Condptrs*)&o[5];
alignarg(more_OPIF:)
			o = *(int**)&cp[2];
			i = 0;
			if (cp[0].bder >= 0)
				i = heswork(cp[0].f, cp[0].b);
			if (cp[1].bder >= 0) {
				j = heswork(cp[1].f, cp[1].b);
				if (i < j)
					i = j;
				}
			n += i + 2;
			break;

/*		case Hv_sumlist: */
		case OPSUMLIST1:
			n += o[3];
			o += o[3] + 4;
			break;

/*		case Hv_func: */
#ifdef X64_bit_pointers
		case OP_FUNCALL1align:
			ptfi = (tfinfo**)(o+4);
			goto more_func;
#endif
		case OP_FUNCALL1:
			ptfi = (tfinfo**)(o+3);
alignarg(more_func:)
			tfi = *ptfi;
			i = tfi->nd;
			n += i*(i+4);
			o = (int*)(ptfi+1) + tfi->n;
			break;

/*		case Hv_plusR:	*/
/*		case Hv_plusL:	*/
/*		case Hv_minusR:	*/
		case OPPLUS01:
		case OPPLUS10:
		case OPMINUS01:
		case OPMINUS10:
			n += 2;
			o += 5;
			break;

/*		case Hv_plusLR:	*/
/*		case Hv_minusLR:*/
		case OPPLUS2:
		case OPMINUS2:
			n += 3;
			o += 5;
			break;

#ifdef X64_bit_pointers
		case OP_GOTOalign:
		case OPGOTOFalign:
			o = *(int**)&o[2];
			break;
#endif
		case OP_GOTO:
		case OPGOTOF:
			o = *(int**)&o[1];
			break;

#ifdef X64_bit_pointers
		case OPGOTO2align:
		case OPGOTOF2align:
		case OPGOTOF2nalign:
#endif
		case OPRET:
		case OPGOTO2:
		case OPGOTOF2:
		case OPGOTOF2n:
			goto done;

		case OPVARREF:
			++n;
			o += 3;
			break;

		case OP_NEXTBLKalign:
			++o;
		case OP_NEXTBLK:
			o = *(int**)(o+1);
			break;

		default:/*DEBUG*/
			fprintf(Stderr, "bad *o = %d in heswork\n", *o);
			exit(1);
		}
	    }
 done:
	return n;
	}

 static cexp *
hesfunnel(Static *S, int *ip, int nrefs, ograd **ogp, linarg ***lapp, linarg ***lapep)
{
	ASLTYPE *asl;
	cexp *c;
	int hw, i, k, m, n, nldv, *ov, *ove, *vr, *vre;
	linarg *la, **lap, **lape;
	ograd *og;
	dv_info *dvi;
	range *r;

	asl = S->asl;
	i = *ip;
	c = cexps + i;
	vr = c->vref;
#if 1 /*20170615*/
	if (vr && vr[1] < vr[0])
		return 0;
#endif
	n = 0;
	if (i >= Ncom) {
		vr = 0;
		r = asl->P.Split_ce[i-Ncom].r;
		lap = r->lap;
		lape = lap + (n = r->n);
		}
	else if (!(lap = (dvi = asl->P.dv + i)->nl)) {
		asl->P.linmultr++;
		og = dvi->ll;
		if (og->varno < 0)
			og = og->next;
		for(*ogp = og; og; og = og->next)
			n++;
		if (n > 1 && asl->p.hffactor > 0) {
			asl->P.linhesfun++;
			*ip = n;
			return c;
			}
		return 0;
		}
	else {
		lape = lap;
		while(*lape)
			lape++;
		n = lape - lap;
		}
	if (!c->o.f)
		return 0;
	*lapp = lap;
	*lapep = lape;
	*ogp = 0;
	asl->P.nlmultr++;
	while(lap < lape) {
		la = *lap++;
		ov = la->ov;
		for(ove = ov + la->nnz; ov < ove; ++ov)
			if (!zc[*ov]++)
				zci[nzc++] = *ov;
		}
	m = nzc;
	while(nzc > 0)
		zc[zci[--nzc]] = 0;
	if (vr) {
		vre = vr + vr[0] + 3;
		nldv = asl->i.defvar0 + asl->P.ncom + asl->P.ndvspout;
		for(vr += vr[1] + 3; vr < vre && *vr < nldv; ++vr)
			++n;
		}
	if (m > n)
		m = n;
	if ((k = m*nrefs - n) <= 0)
		return 0;
	hw = heswork(c->o.f, c->o.b);
	if (k*hw <= nrefs*n*(n+3)*asl->p.hffactor)
		return 0;
	*ip = n;
	asl->P.nlhesfun++;
	return c;
	}

 ASL_pfgh *
pscheck_ASL(ASL *a, const char *who1)
{
	ASL_CHECK(a, ASL_read_pfgh, who1);
	return (ASL_pfgh*)a;
	}

 static void
hes_setup(Static *S)
{
	ASL_pfgh *asl;
	cexp *C;
	hes_fun *hf, *hfth;
	int *c, *ce, *cei, *dvsp0, h, h0, i, i0, k, k1, k2, n, n0, ncom;
	int *ndvsp, nmax, ns, nvx, *vp, *vr, *vre, *z;
	linarg *la, **lap, **lape;
	ograd *og;
	psb_elem *b;
	range *r, *r0;

	asl = S->asl;
	ncom = asl->P.ncom;
	nvx = asl->i.defvar0;
	nzclim = max_var1 >> 3;
	zl = asl->P.zlsave;

	zc_upgrade(S);	/* make zc and zci available */
	n = nmax = 0;
	r0 = (range*)&asl->P.rlist;
	h = wlast;
	dvsp0 = asl->P.dvsp0;
	ndvsp = asl->P.ndvsp;
	for(r = asl->P.rlist.next; r != r0; r = r->rlist.next) {
		if (r->n <= 0) {
			r->cei = 0;
			continue;
			}
		if (nmax < r->n)
			nmax = r->n;
		n0 = 0;
		cei = 0;
		for(b = r->refs; b; b = b->next) {
			if ((c = b->ce)) {
				if (!n0) {
					cei = c;
					n0 = *c;
					}
				ce = c + *c;
				while(c < ce) {
					if (!zc[i = *++c]++) {
						zci[n++] = i;
						if (i < ncom && (ns = ndvsp[i])) {
							k1 = dvsp0[i];
							for(k2 = k1 + ns; k1 < k2; ++k1) {
								if (!zc[k1]++)
									zci[n++] = k1;
								}
							}
						}
					}
				}
			}
		if (n0 < n)
			cei = 0;
		else if ((r->cei = cei)) {
			while(n > 0)
				zc[zci[--n]] = 0;
			continue;
			}
		if (!n)
			continue;
		r->cei = cei = (int*)mem((n+1)*sizeof(int));
		*cei++ = n;
		if (n > 1)
			qsortv(zci, n, sizeof(int), hscompar, S);
		cei += n;
		do zc[*--cei = zci[--n]] = 0;
			while(n > 0);
		}
	z = zc + nvx;
	k = -1;
	for(r = asl->P.rlist.next; r != r0; r = r->rlist.next) {
		if (!(cei = r->cei))
			continue;
		c = cei + 1;
		ce = c + *cei;
		do {
			if (z[i = *c++]++ && k < i)
				k = i;
			} while(c < ce);
		}
	hfth = 0;
	asl->P.linmultr = asl->P.linhesfun = 0;
	asl->P.nlmultr = asl->P.nlhesfun = 0;
	lap = lape = 0; /* silence false "used uninitialized" warning */
	for(; k >= 0; --k)
		if (z[i = k] > 1
		 && (C = hesfunnel(S,&i,z[i],&og,&lap,&lape))) {
			hf = (hes_fun*)mem(sizeof(hes_fun));
			hf->hfthread = hfth;
			C->hfun = hfth = hf;
			hf->c = C;
			hf->n = hf->nd = i;
			if ((hf->og = og)) {
				hf->vp = 0;
				hf->grdhes = 0;
				}
			else {
				i0 = i;
				if ((vr = C->vref))
					hf->nd = i += vr[0] - vr[1];
				hf->vp = vp = (int *)mem(hf->nd*sizeof(int));
				while(lap < lape) {
					la = *lap++;
					*vp++ = la->u.v;
					}
				if (vr) {
					for(vre = vr + *vr + 3, vr += vr[1] + 3; vr < vre; ++vr)
						*vp++ = *vr;
					}
				hf->grdhes = h;
				h += i + i0*i0;
				}
			}
	if ((asl->I.hesthread = hfth))
		asl->i.x0kindinit |= ASL_need_funnel;
	del_mblk(zci);
	asl->P.khesoprod = 5;
	asl->P.nmax = nmax;
	asl->P.rtodo = h0 = h;
	h += n = (nvx*sizeof(range*) + sizeof(real) - 1) / sizeof(real);
	asl->P.utodo = h;
	h += n;
	asl->P.otodo = h;
	h += n;
	asl->P.dOscratch = h;
	h += nmax;
	asl->P.iOscratch = h;
	h += (nmax*sizeof(int) + sizeof(real) - 1) / sizeof(real);
	wlast = h;
	asl->P.zaplen = (h - h0)*sizeof(real);
	}
#endif	/* PSHVREAD */

#ifdef __cplusplus
}
#endif
