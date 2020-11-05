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

#include "asl.h"
#include "errchk.h"

 void
repwhere_ASL(EvalWorkspace *ew, int jv)
{
	ASL *asl;
	FILE *f;
	char *b, buf[512];
	int i, j, k, k1;
	static const char *what[2] = { "constraint", "objective" };
	static const char *which[4] = { "function: ", "gradient: ", "Hessian: ", "???: " };

	asl = ew->asl;
	fflush(stdout);
	need_nl = 0;
	fprintf(Stderr, "Error evaluating ");

#define next_line fgets(buf,sizeof(buf),f)

	if ((i = ew->cv_index)) {
		strcpy(stub_end, ".fix");
		j = 0;
		if ((f = fopen(filename, "r"))) {
			for(;;) {
				if (!next_line)
					goto eof;
				for(b = buf; *b; b++)
					if (*b == '=') {
						while(++j < i)
							if (!next_line)
								goto eof;
						b = buf;
						while(*b && *b != '=')
							b++;
						if (*b != '=' || b < buf + 2)
							j = 0;
						else
							b[-1] = 0;
						goto eof;
						}
				}
 eof:
			fclose(f);
			}
		if (j == i)
			fprintf(Stderr, "var %s: ", buf);
		else
			fprintf(Stderr, "\"var =\" definition %d: ", i);
		goto ret;
		}

	k = k1 = 0;
	if ((i = ew->co_index) < 0) {
		k = 1;
		i = asl->i.n_con0 - i;
		if (n_obj <= 1)
			k1 = 1;
		}
	else
		++i;
	fprintf(Stderr, "%s ", what[k]);
	if (maxrownamelen) {
		strcpy(stub_end, ".row");
		if ((f = fopen(filename, "r"))) {
			for(j = 0; j <= i; j++)
				if (!next_line)
					break;
			fclose(f);
			if (j >= i) {
				for(b = buf; *b; b++)
					if (*b == '\n') {
						*b = 0;
						break;
						}
				fprintf(Stderr, "%s ", buf);
				}
			}
		}
	else
		if (k1 == 0)
			fprintf(Stderr, "%d ", i);
	fprintf(Stderr, "%s", which[jv-1]);
 ret:
	errno = 0;	/* in case it was set by fopen */
	fflush(Stderr);
	}

 void
report_where_ASL(ASL *asl)
{
	EvalWorkspace *ew;

	if (!(ew = asl->i.Ew0)) {
		switch(asl->i.ASLtype) {
		  case ASL_read_f:
		  case ASL_read_fg:
			ew = ewalloc1_ASL(asl);
			break;
		  default:
			ew = ewalloc2_ASL(asl);
		  }
		asl->i.Ew0 = ew;
		}
	repwhere_ASL(ew, 1);
	}

 static void
jmp_check(Jmp_buf *J, int jv)
{
	if (J)
		longjmp(J->jb, jv);
	}

 static void
Errprint(const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
#ifndef NO_PERROR
	if (errno)
		fprintf(Stderr, "\n%s: ", strerror(errno));
#endif
	vfprintf(Stderr, fmt, ap);
	va_end(ap);
	fflush(Stderr);
	}

 typedef struct DerrRecord DerrRecord;
 typedef void (*DerrPrint)(EvalWorkspace*, DerrRecord*);

 struct
DerrRecord {
	DerrPrint errprint;
	const char *fmt, *who;
	real a;
	union { const char *s; real b; } u;
	int jv;
	int dv;
	};

 static void
derrprint1(EvalWorkspace *ew, DerrRecord *R)
{
	fprintf(Stderr, R->fmt, R->who, R->a);
	}

 static void
derrprint2(EvalWorkspace *ew, DerrRecord *R)
{
	fprintf(Stderr, R->fmt, R->who, R->a, R->u.b);
	}

 static void
derrprintf(EvalWorkspace *ew, DerrRecord *R)
{
	fprintf(Stderr, R->fmt, R->who, R->u.s);
	}

 typedef struct DerrMblock DerrMblock;
 struct
DerrMblock {
	DerrMblock *next;
	size_t len;
	real align[1]; /* would prefer align[0], but some older compilers would complain */
	};

 struct
DerivErrInfo {
	DerrMblock *curmb, *freemb;
	char *mbnext, *mblast;
	DerrRecord **R;
	int *busy;
	int nbusy;
	};

 void
deriv_errchk_ASL(EvalWorkspace *ew, int coi, int n, int jv)
{
	ASL *asl;
	DerivErrInfo *D;
	DerrRecord *R, **Rp, **Rpe;
	int k;

	asl = ew->asl;
	D = ew->Derrs;
	if ((k = coi) < 0) {
		k = -(k + 1);
		if (k >= nlo)
			return;
		k += nlc;
		}
	else if (k >= nlc)
		return;
	Rp = D->R + k;
	for(Rp = D->R + k, Rpe = Rp + n; Rp < Rpe; ++Rp, ++coi)
		if ((R = *Rp) && R->jv <= jv) {
			jmp_check(ew->err_jmpw, R->jv);
			if (ew == asl->i.Ew0)
				jmp_check(asl->i.err_jmp_, R->jv);
			ew->co_index = coi;
			ew->cv_index = R->dv;
			repwhere_ASL(ew, R->jv);
			R->errprint(ew,R);
			fflush(Stderr);
			jmp_check(ew->err_jmpw1, R->jv);
			if (ew == asl->i.Ew0)
				jmp_check(asl->i.err_jmp1_, R->jv);
			exit(1);
			}
	}

 void
deriv2_errchk_ASL(EvalWorkspace *ew, int jv)
{
	ASL *asl;
	DerivErrInfo *D;
	DerrRecord *R, **Rp;
	int coi, k, ke;

	asl = ew->asl;
	D = ew->Derrs;
	for(ke = nlc + nlo, k = 0; k < ke; ++k) {
		Rp = D->R + k;
		if ((R = *Rp) && R->jv <= jv) {
			jmp_check(ew->err_jmpw, R->jv);
			if (ew == asl->i.Ew0)
				jmp_check(asl->i.err_jmp_, R->jv);
			if ((coi = k) >= nlc)
				coi = nlc - k - 1;
			ew->co_index = coi;
			ew->cv_index = R->dv;
			repwhere_ASL(ew, R->jv);
			R->errprint(ew,R);
			fflush(Stderr);
			jmp_check(ew->err_jmpw1, R->jv);
			if (ew == asl->i.Ew0)
				jmp_check(asl->i.err_jmp1_, R->jv);
			exit(1);
			}
		}
	}

 static DerivErrInfo *
new_DerrMblock(EvalWorkspace *ew, size_t len)
{
	ASL *asl;
	DerivErrInfo *D;
	DerrMblock *M, **Mp;
	char *s;
	int nlco;
	size_t L, L1;

	asl = ew->asl;
	len = len < 4096
		? 4096
		: (len + sizeof(real) - 1) & ~(sizeof(real) - 1);
	if (!(D = ew->Derrs)) {
		if ((D = ew->Derrs0)) {
			ew->Derrs = D;
			M = D->curmb;
			if (M->len >= len)
				return D;
			}
		else {
			nlco = nlc + nlo;
			L = sizeof(DerivErrInfo)
				+ nlco*(sizeof(int) + sizeof(DerrRecord*));
			L = (L + sizeof(real) - 1) & ~(sizeof(real) - 1);
			L1 = L + (sizeof(DerrMblock) - sizeof(real)) + len;
			D = (DerivErrInfo*)M1alloc(L1);
			memset(D, 0, L);
			ew->Derrs = ew->Derrs0 = D;
			D->R = (DerrRecord**)(D+1);
			D->busy = (int*)(D->R + nlco);
			M = (DerrMblock*)((char*)D + L);
			M->len = len;
			goto have_M;
			}
		}
	for(Mp = &D->freemb;; Mp = &M->next) {
		if (!(M = *Mp)) {
			M = (DerrMblock*)M1alloc((sizeof(DerrMblock) - sizeof(real)) + len);
			M->len = len;
			break;
			}
		if (M->len >= len) {
			*Mp = M->next;
			break;
			}
		}
 have_M:
	M->next = D->curmb;
	D->curmb = M;
	D->mbnext = s = (char*)M->align;
	D->mblast = s + M->len;
	return D;
	}

 static DerrRecord *
getDR(EvalWorkspace *ew, int jv)
{
	ASL *asl;
	DerivErrInfo *D;
	DerrRecord *R;
	int i, j, je, k;
	size_t L;

	asl = ew->asl;
	if ((k = ew->co_index) < 0) {
		k = -(k + 1);
		if (k >= nlo)
			return 0;
		k += nlc;
		}
	else if (k >= nlc)
		return 0;
	L = (sizeof(DerrRecord) + sizeof(real)-1) & ~(sizeof(real)-1);
	R = 0;
	if ((D = ew->Derrs)) {
		if ((R = D->R[k]) && R->jv <= jv)
			return 0;
		if (L <= D->mblast - D->mbnext)
			goto have_D;
		}
	D = new_DerrMblock(ew, L);
 have_D:
	if (!R) {
		R = (DerrRecord*)(D->mblast - L);
		D->mblast = (char*)R;
		}
	D->R[k] = R;
	D->busy[D->nbusy++] = k;
	if ((R->dv = i = ew->cv_index)) {
		j = 0;
		je = nlc + nlo;
		if (i > comb) {
			if (i <= combc)
				je = nlc;
			else if (i <= ncom0)
				j = combc;
			}
		for(; j < je; ++j) {
			if (!D->R[j]) {
				D->R[j] = R;
				D->busy[D->nbusy++] = j;
				}
			}
		}
	return R;
	}

 void
deriv_errclear_ASL(EvalWorkspace *ew)
{
	DerivErrInfo *D;
	DerrMblock *M, *M0, *M1;
	DerrRecord **R;
	char *s;
	int *b, *be;

	D = ew->Derrs;
	ew->Derrs = 0;
	R = D->R;
	for(b = D->busy, be = b + D->nbusy; b < be; ++b)
		R[*b] = 0;
	D->nbusy = 0;
	M0 = D->freemb;
	for(M = D->curmb; M; M0 = M, M = M1) {
		M1 = M->next;
		M->next = M0;
		}
	D->freemb = M0->next;
	M0->next = 0;
	D->curmb = M0;
	D->mbnext = s = (char*)M0->align;
	D->mblast = s + M0->len;
	}

 void
introuble_ASL(EvalWorkspace *ew, const char *who, real a, int jv)
{
	ASL *asl;
	static const char fmt[] = "can't evaluate %s(%g).\n";
#ifndef ASL_OLD_DERIV_CHECK /*{*/
	DerrRecord *R;

	if (jv > 1 && !(ew->wantderiv & 2)) {
		if ((R = getDR(ew,jv))) {
			R->errprint = derrprint1;
			R->a = a;
			R->jv = jv;
			R->fmt = fmt;
			R->who = who;
			}
		return;
		}
#endif /*}*/
	jmp_check(ew->err_jmpw, jv);
	asl = ew->asl;
	if (ew == asl->i.Ew0)
		jmp_check(asl->i.err_jmp_, jv);
	repwhere_ASL(ew, jv);
	Errprint(fmt, who, a);
	jmp_check(ew->err_jmpw1, jv);
	if (ew == asl->i.Ew0)
		jmp_check(asl->i.err_jmp1_, jv);
	exit(1);
	}

 void
introuble2_ASL(EvalWorkspace *ew, const char *who, real a, real b, int jv)
{
	ASL *asl;
	static const char fmt[] = "can't evaluate %s(%g,%g).\n";
#ifndef ASL_OLD_DERIV_CHECK /*{*/
	DerrRecord *R;

	if (jv > 1 && !(ew->wantderiv & 2)) {
		if ((R = getDR(ew,jv))) {
			R->errprint = derrprint2;
			R->a = a;
			R->u.b = b;
			R->jv = jv;
			R->fmt = fmt;
			R->who = who;
			}
		return;
		}
#endif /*}*/
	jmp_check(ew->err_jmpw, jv);
	asl = ew->asl;
	if (ew == asl->i.Ew0)
		jmp_check(asl->i.err_jmp_, jv);
	repwhere_ASL(ew, jv);
	Errprint(fmt, who, a, b);
	jmp_check(ew->err_jmpw1, jv);
	if (ew == asl->i.Ew0)
		jmp_check(asl->i.err_jmp1_, jv);
	exit(1);
	}

 void
zero_div_ASL(EvalWorkspace *ew, real L, const char *op)
{
	ASL *asl;

	errno_set(EDOM);
	jmp_check(ew->err_jmpw, 1);
	asl = ew->asl;
	if (ew == asl->i.Ew0)
		jmp_check(asl->i.err_jmp_, 1);
	repwhere_ASL(ew, 4);
	fprintf(Stderr, "can't compute %g%s0.\n", L, op);
	fflush(Stderr);
	jmp_check(ew->err_jmpw1, 1);
	if (ew == asl->i.Ew0)
		jmp_check(asl->i.err_jmp1_, 1);
	exit(1);
	}

 void
fintrouble_ASL(EvalWorkspace *ew, func_info *fi, const char *s, TMInfo *T)
{
	ASL *asl;
	TMInfo *T1, *T1prev;
	int jv;
	static const char fmt[] = "Error in function %s:\n\t%s\n";

	jv = 1;
	switch(*s) {
	 case '\'':
		jv = 2;
		goto inc_s;
	 case '"':
		jv = 3;
 inc_s:
		++s;
	 }
#ifndef ASL_OLD_DERIV_CHECK /*{*/
	if (jv > 1 && !(ew->wantderiv & 2)) {
		DerivErrInfo *D;
		DerrRecord *R;
		size_t L;

		if ((R = getDR(ew,jv))) {
			D = ew->Derrs;
			L = strlen(s) + 1;
			if (L > D->mblast - D->mbnext)
				D = new_DerrMblock(ew, L);
			memcpy(D->mbnext, s, L);
			R->u.s = D->mbnext;
			D->mbnext += L;
			R->errprint = derrprintf;
			R->jv = jv;
			R->fmt = fmt;
			R->who = fi->name;
			}
		return;
		}
#endif /*}*/
	jmp_check(ew->err_jmpw, jv);
	asl = ew->asl;
	if (ew == asl->i.Ew0)
		jmp_check(asl->i.err_jmp_, jv);
	repwhere_ASL(ew, jv);
	fprintf(Stderr, fmt, fi->name, s);
	fflush(Stderr);
	for(T1 = T->u.prev; T1; T1 = T1prev) {
		T1prev = T1->u.prev;
		free(T1);
		}
	jmp_check(ew->err_jmpw1, jv);
	if (ew == asl->i.Ew0)
		jmp_check(asl->i.err_jmp1_, jv);
	exit(1);
	}
