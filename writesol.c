/****************************************************************
Copyright (C) 1997, 1999, 2000 Lucent Technologies
All Rights Reserved

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the name of Lucent or any of its entities
not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.
****************************************************************/

#include "getstub.h"

 typedef struct
SufHead {
	char sufid[8];
	fint kind;
	fint n;
	fint namelen;
	fint tablen;
	} SufHead;

#ifdef __cplusplus
extern "C" {
#endif

typedef char *(*Name) ANSI((ASL*,int));

#ifdef __cplusplus
}
#endif

 static void
#ifdef KR_headers
getsufhead(asl, d, sh, np, zp) ASL *asl; SufDesc *d; SufHead *sh; int *np, **zp;
#else
getsufhead(ASL *asl, SufDesc *d, SufHead *sh, int *np, int **zp)
#endif
{
	int i, *ip, *ipe, n, nz;
	real *rp, *rpe;

	memcpy(sh->sufid, "\nSuffix\n", 8);
	sh->kind = d->kind &
		(ASL_Sufkind_mask | ASL_Sufkind_real | ASL_Sufkind_iodcl);
	*np = n = (&asl->i.n_var_)[i = d->kind & ASL_Sufkind_mask];
	*zp = i < 2 ? asl->i.z[i] : 0;
	nz = 0;
	if (d->kind & ASL_Sufkind_real) {
		rp = d->u.r;
		rpe = rp + n;
		while(rp < rpe)
			if (*rp++)
				nz++;
		}
	else {
		ip = d->u.i;
		ipe = ip + n;
		while(ip < ipe)
			if (*ip++)
				nz++;
		}
	sh->n = nz;
	sh->namelen = strlen(d->sufname) + 1;
	sh->tablen = 0;
	if (d->table)
		sh->tablen = strlen(d->table) + 1;
	}

 static long
#ifdef KR_headers
tablines(s) char *s;
#else
tablines(char *s)
#endif
{
	long n;
	if (!s)
		return 0;
	n = 1;
	while(*s)
		if (*s++ == '\n')
			n++;
	return n;
	}

 static void
#ifdef KR_headers
showsol(asl, x, n, n0, z, name, what, pfix)
	ASL *asl; real *x; int n, n0, *z; Name name; char *what, *pfix;
#else
showsol(ASL *asl, real *x, int n, int n0, int *z, Name name, char *what, char *pfix)
#endif
{
	int i, j, k, k0;

	if (!x || n <= 0)
		return;
	k0 = k = strlen(what);
	for(i = 0; i < n; i++)
		if ((j = strlen((*name)(asl,i))) > k)
			k = j;
	k += 2;
	printf("\n%s%*s%svalue\n", what, k-k0, "", pfix);
	for(i = 0; i < n0; i++)
		if ((j = z ? z[i] : i) >= 0)
			printf("%-*s%.g\n", k, (*name)(asl,j), x[i]);
	}

 static real *
#ifdef KR_headers
scale(x, s, y, n) real *x, *s, *y; int n;
#else
scale(real *x, real *s, real *y, int n)
#endif
{
	real *y0, *xe;

	y0 = y;
	xe = x + n;
	while(x < xe)
		*y++ = *s++ * *x++;
	return y0;
	}

 static void
#ifdef KR_headers
copyup(n, z, x) int n; int *z; real *x;
#else
copyup(int n, int *z, real *x)
#endif
{
	int j;

	while(--n >= 0)
		x[n] = (j = z[n]) >= 0 ? x[j] : 0.;
	}

 static void
#ifdef KR_headers
equ_adjust1(ip, LU, U, n) int *ip; real *LU; real *U; int n;
#else
equ_adjust1(int *ip, real *LU, real *U, int n)
#endif
{
	int i = 0;
	if (U) {
		for(; i < n; i++)
			if (LU[i] == U[i] && (ip[i] == 3 || ip[i] == 4))
				ip[i] = 5;
		}
	else if (LU)
		for(; i < n; i++, LU += 2)
			if (LU[0] == LU[1] && (ip[i] == 3 || ip[i] == 4))
				ip[i] = 5;
	}

 void
#ifdef KR_headers
equ_adjust_ASL(asl, cstat, rstat) ASL *asl; int *cstat; int *rstat;
#else
equ_adjust_ASL(ASL *asl, int *cstat, int *rstat)
#endif
{
	if (cstat)
		equ_adjust1(cstat, LUv, Uvx, n_var);
	if (rstat)
		equ_adjust1(rstat, LUrhs, Urhsx, n_con);
	}

 void
#ifdef KR_headers
write_sol_ASL(asl, msg, x, y, oi)
	ASL *asl; char *msg; double *x, *y; Option_Info *oi;
#else
write_sol_ASL(ASL *asl, char *msg, double *x, double *y, Option_Info *oi)
#endif
{
	FILE *f;
	int N, binary, i, i1, *ip, j, k, n, tail, wantsol, *zz;
	char buf[80], *s, *s1, *s2;
	static char *wkind[] = {"w", "wb"};
	ftnlen L[6];
	fint J[2], m, z[4];
	size_t nn;
	real *rp, *y1, *xycopy;
	SufDesc *d;
	SufHead sh;

	if (!asl || asl->i.ASLtype < 1 || asl->i.ASLtype > 5)
		badasl_ASL(asl,0,"write_sol");

	xycopy = 0;
	if ((wantsol = oi ? oi->wantsol : 1) || amplflag) {
		k = 0;
		if (x && (asl->i.vscale || asl->i.z[0]))
			k = asl->i.n_var0;
		if (y && (asl->i.cscale || asl->i.z[1]))
			k += asl->i.n_con0;
		if (k)
			y1 = xycopy = (real*)Malloc(k*sizeof(real));
		if (x) {
			if (asl->i.vscale) {
				x = scale(x, asl->i.vscale, y1, n_var);
				y1 += asl->i.n_var0;
				}
			else if (asl->i.z[0]) {
				memcpy(y1, x, n_var*sizeof(real));
				x = y1;
				y1 += asl->i.n_var0;
				}
			if (asl->i.z[0])
				copyup(asl->i.n_var0, asl->i.z[0], x);
			}
		z[0] = m = asl->i.n_con0;
		if (!y)
			m = 0;
		else {
			if (asl->i.cscale)
				y = scale(y, asl->i.cscale, y1, n_con);
			else if (asl->i.z[1]) {
				memcpy(y1, y, n_con*sizeof(real));
				y = y1;
				}
			if (asl->i.z[1])
				copyup(asl->i.n_con0, asl->i.z[1], y);
			}
		}
	if (!amplflag && !(wantsol & 1))
		goto write_done;
	tail = 0;
	if (obj_no || solve_code != -1)
		tail = 1;
	else  {
		for(i1 = 0; i1 < 4; i1++)
		    for(d = asl->i.suffixes[i1]; d; d = d->next)
			if (d->kind & ASL_Sufkind_output
			 && (d->kind & ASL_Sufkind_real
					? (int*)d->u.r : d->u.i)) {
				tail = 1;
				goto break2;
				}
		}
 break2:
	binary = binary_nl & 1;
	strcpy(stub_end, ".sol");
	f = fopen(filename, wkind[binary]);
	if (!f) {
		fprintf(Stderr, "can't open %s\n", filename);
		exit(2);
		}
	z[1] = m;
	z[2] = n = asl->i.n_var0;
	if (!x)
		n = 0;
	z[3] = n;
	k = (int)ampl_options[0];
	if (binary) {
		L[0] = 6;
		L[1] = strlen(msg);
		L[2] = 0;
		L[3] = (ampl_options[0] + 5)*sizeof(fint) + 7;
		L[4] = m*sizeof(double);
		L[5] = n*sizeof(double);
		fwrite(L, sizeof(ftnlen), 1, f);
		fwrite("binary", 6, 1, f);
		fwrite(L, sizeof(ftnlen), 2, f);
		if (L[1]) {
			fwrite(msg, L[1], 1, f);
			fwrite(L+1, sizeof(ftnlen), 2, f);
			}
		if (k) {
			fwrite(L+2, sizeof(ftnlen), 2, f);
			fwrite("Options",7,1,f);
			nn = (size_t)ampl_options[0]+1;
			if (ampl_options[2] == 3)
				ampl_options[0] += 2;
			fwrite(ampl_options, sizeof(fint), nn, f);
			fwrite(z, sizeof(fint), 4, f);
			if (ampl_options[2] == 3)
				fwrite(&ampl_vbtol, sizeof(real), 1, f);
			fwrite(L+3, sizeof(ftnlen), 2, f);
			}
		else {
			fwrite(L+2, sizeof(ftnlen), 1, f);
			fwrite(L+4, sizeof(ftnlen), 1, f);
			}
		if (y)
			fwrite(y, sizeof(double), m, f);
		fwrite(L+4, sizeof(ftnlen), 2, f);
		if (x)
			fwrite(x, sizeof(double), n, f);
		fwrite(L+5, sizeof(ftnlen), 1, f);
		if (tail)
		  switch(asl->i.flags & 1) {
		    case 0:
			if (obj_no) {
				L[0] = L[2] = sizeof(fint);
				L[1] = obj_no;
				fwrite(L, sizeof(fint), 3, f);
				}
			break;
		    case 1:
			L[0] = L[3] = 2*sizeof(fint);
			L[1] = obj_no;
			L[2] = solve_code;
			fwrite(L, sizeof(fint), 4, f);
			for(i1 = 0; i1 < 4; i1++)
			  for(d = asl->i.suffixes[i1]; d; d = d->next)
			    if (d->kind & ASL_Sufkind_output
			     && (d->kind & ASL_Sufkind_real
					? (int*)d->u.r : d->u.i)) {
				getsufhead(asl, d, &sh, &N, &zz);
				L[0] = sizeof(sh) + sh.namelen + sh.tablen
					+ sh.n*(sizeof(int) +
						(d->kind & ASL_Sufkind_real
						? sizeof(real) : sizeof(int)));
				fwrite(L, sizeof(fint), 1, f);
				fwrite(&sh, sizeof(sh), 1, f);
				fwrite(d->sufname, sh.namelen, 1, f);
				if (sh.tablen)
					fwrite(d->table, sh.tablen, 1, f);
				i = j = 0;
				if (d->kind & ASL_Sufkind_real)
				    for(rp = d->u.r; i < N; i++) {
					if (rp[i]) {
						if (zz) {
							while(zz[j] < i)
								j++;
							J[0] = j;
							}
						else
							J[0] = i;
						fwrite(J, sizeof(fint), 1, f);
						fwrite(rp+i,sizeof(real),1,f);
						}
					}
				else
				    for(ip = d->u.i; i < N; i++) {
					if (J[1] = ip[i]) {
						if (zz) {
							while(zz[j] < i)
								j++;
							J[0] = j;
							}
						else
							J[0] = i;
						fwrite(J, sizeof(fint), 2, f);
						}
					}
				fwrite(L, sizeof(fint), 1, f);
				}
			}
		}
	else {
		if (*(s = msg)) {
			for(s2 = s + strlen(s); s2 > s && s2[-1] == '\n'; --s2);
			while (s < s2) {
				for(s1 = s; *s1 != '\n' && ++s1 < s2;);
				fprintf(f, s == s1 ? " \n" : "%.*s\n",s1-s,s);
				s = s1 + 1;
				}
			}
		fprintf(f, "\n");
		if (k = (int)ampl_options[0]) {
			if (ampl_options[2] == 3)
				ampl_options[0] += 2;
			fprintf(f, "Options\n");
			for(i = 0; i <= k; i++)
				fprintf(f,"%ld\n",(long)ampl_options[i]);
			fprintf(f,"%ld\n%ld\n%ld\n%ld\n",
				(long)z[0],(long)z[1],(long)z[2],(long)z[3]);
			if (ampl_options[2] == 3) {
				g_fmtp(buf, ampl_vbtol, 0);
				fprintf(f, "%s\n", buf);
				}
			}
		y1 = y;
		while(--m >= 0) {
			g_fmtp(buf, *y1++, 0);
			fprintf(f,"%s\n", buf);
			}
		y1 = x;
		while(--n >= 0) {
			g_fmtp(buf, *y1++, 0);
			fprintf(f, "%s\n", buf);
			}
		if (tail)
		  switch(asl->i.flags & 1) {
		    case 0:
			if (obj_no)
				fprintf(f, "objno %d\n", obj_no);
			break;
		    case 1:
			fprintf(f, "objno %d %d\n", obj_no, solve_code);
			for(i1 = 0; i1 < 4; i1++)
			  for(d = asl->i.suffixes[i1]; d; d = d->next)
			    if (d->kind & ASL_Sufkind_output
			     && (d->kind & ASL_Sufkind_real
					? (int*)d->u.r : d->u.i)) {
				getsufhead(asl, d, &sh, &N, &zz);
				fprintf(f, "suffix %ld %ld %ld %ld %ld\n%s\n",
					(long)sh.kind, (long)sh.n,
					(long)sh.namelen, (long)sh.tablen,
					tablines(d->table), d->sufname);
				if (sh.tablen)
					fprintf(f, "%s\n", d->table);
				i = j = 0;
				if (d->kind & ASL_Sufkind_real)
				    for(rp = d->u.r; i < N; i++) {
					if (rp[i]) {
						if (zz)
							while(zz[j] < i)
								j++;
						else
							j = i;
						fprintf(f, "%d %.g\n",
							j, rp[i]);
						}
					}
				else
				    for(ip = d->u.i; i < N; i++) {
					if (ip[i]) {
						if (zz)
							while(zz[j] < i)
								j++;
						else
							j = i;
						fprintf(f, "%d %d\n",
							j, ip[i]);
						}
					}
				}
			}
		}
	fclose(f);
 write_done:
	if (i = need_nl)
		if (i > sizeof(buf)-1 || i < 0)
			printf("\n");
		else {
			buf[i] = 0;
			do buf[--i] = '\b';
				while(i > 0);
			printf(buf);
			}
	if (!amplflag) {
		if (!(wantsol & 8))
			printf("%s\n", msg);
		if (wantsol & 2)
			showsol(asl, x, n_var, asl->i.n_var0, asl->i.z[0],
				var_name_ASL, "variable", "");
		if (wantsol & 4)
			showsol(asl, y, n_con, asl->i.n_con0, asl->i.z[1],
				con_name_ASL, "constraint", "dual ");
		}
	if (xycopy)
		free(xycopy);
	}
