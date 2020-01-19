/****************************************************************
Copyright (C) 1997, 1999 Lucent Technologies
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

#include "asl.h"

 void
#ifdef KR_headers
name_map_ASL(n, no, z, nam) int n; int no; int *z; char **nam;
#else
name_map_ASL(int n, int no, int *z, char **nam)
#endif
{
	int i, j, k;

	for(i = k = 0; i < n; i++) {
		if ((j = z[i]) >= 0)
			nam[k = j] = nam[i];
		}
	n += no;
	while(no--)
		nam[++k] = nam[i++];
	while(++k < n)
		nam[k] = 0;
	}

 static char **
#ifdef KR_headers
get_names(asl, suf, no, n, z) ASL *asl; char *suf; int no, n, *z;
#else
get_names(ASL *asl, char *suf, int no, int n, int *z)
#endif
{
	FILE *f;
	char buf[512], **nam, **ne, **rv, *s;
	int nt = n + no;

	nam = rv = (char**)mem(nt*sizeof(char*));
	ne = nam + nt;
	strcpy(stub_end, suf);
	if (f = fopen(filename, "r")) {
		while(nam < ne && fgets(buf,sizeof(buf),f)) {
			for(s = buf; *s && *s != '\n'; s++);
			*s = 0;
			strcpy(*nam++ = (char*)mem(s-buf+1), buf);
			}
		fclose(f);
		}
	while(nam < ne)
		*nam++ = 0;
	if (z)
		name_map_ASL(n, no, z, rv);
	return rv;
	}

 char *
#ifdef KR_headers
con_name_ASL(asl, n) ASL *asl; int n;
#else
con_name_ASL(ASL *asl, int n)
#endif
{
	char buf[32], **np, *rv;
	if (n < 0 || n >= n_con)
		return "**con_name(bad n)**";
	if (!asl->i.conames)
		asl->i.conames = get_names(asl, ".row", n_obj,
					asl->i.n_con0, asl->i.z[1]);
	np = asl->i.conames + n;
	if (!(rv = *np)) {
		*np = rv = (char*)mem(Sprintf(buf,"_scon[%d]",n+1)+1);
		strcpy(rv, buf);
		}
	return rv;
	}

 char *
#ifdef KR_headers
obj_name_ASL(asl, n) ASL *asl; int n;
#else
obj_name_ASL(ASL *asl, int n)
#endif
{
	char buf[32], **np, *rv;
	if (n < 0 || n >= n_obj)
		return "**obj_name(bad n)**";
	if (!asl->i.conames)
		asl->i.conames = get_names(asl, ".row", n_obj,
					asl->i.n_con0, asl->i.z[1]);
	np = asl->i.conames + (n + n_con);
	if (!(rv = *np)) {
		*np = rv = (char*)mem(Sprintf(buf,"_sobj[%d]",n+1)+1);
		strcpy(rv, buf);
		}
	return rv;
	}

 char *
#ifdef KR_headers
var_name_ASL(asl, n) ASL *asl; int n;
#else
var_name_ASL(ASL *asl, int n)
#endif
{
	char buf[32], **np, *rv;
	if (n < 0 || n >= n_var)
		return "**var_name(bad n)**";
	if (!asl->i.varnames)
		asl->i.varnames = get_names(asl, ".col", 0,
					asl->i.n_var0, asl->i.z[0]);
	np = asl->i.varnames + n;
	if (!(rv = *np)) {
		*np = rv = (char*)mem(Sprintf(buf,"_svar[%d]",n+1)+1);
		strcpy(rv, buf);
		}
	return rv;
	}
