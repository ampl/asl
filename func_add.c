/****************************************************************
Copyright (C) 1997-2000 Lucent Technologies
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

#ifdef __cplusplus
extern "C" {
#endif
extern ASLhead ASLhead_ASL;

char *i_option_ASL;
static int n_added;

 func_info *
#ifdef KR_headers
func_lookup(asl, s, add) ASL *asl; register char *s; int add;
#else
func_lookup(ASL *asl, register const char *s, int add)
#endif
{
	register unsigned x = 0;
	func_info *fi, **finext;
	Const char *s0 = s;

	while(*s)
		x = 31*x + *s++;
	finext = &fhash[x % NFHASH];
	for(fi = *finext; fi; fi = fi->next)
		if (!strcmp(s0, fi->name)) {
			if (add) {
				fprintf(Stderr,
				"addfunc: duplicate function %s\n", s0);
				fi = 0;
				}
			return fi;
			}
	if (add) {
		fi = (func_info *)mem_ASL(asl, sizeof(func_info));
		fi->next = *finext;
		*finext = fi;
		fi->name = s0;
		}
	return fi;
	}

 void
#ifdef KR_headers
addfunc_ASL(fname, f, ftype, nargs, funcinfo, ae)
	char *fname, *funcinfo; ufunc *f; AmplExports *ae;
#else
addfunc_ASL(const char *fname, ufunc *f, int ftype, int nargs, void *funcinfo, AmplExports *ae)
#endif
{
	register func_info *fi;
	ASL *asl = (ASL*)ae->asl;
	if (ftype && ftype != 1) {
#ifndef COMPLAIN_AT_BAD_FTYPE
		if (ftype < 0 || ftype > 6)
#endif
		{
		fprintf(Stderr, "function %s: ftype = %d; expected 0 or 1\n",
			fname, ftype);
		exit(1);
		}
#ifndef COMPLAIN_AT_BAD_FTYPE
		return;
#endif
		}
	if (fi = func_lookup(asl, fname, 1)) {
		n_added++;
		fi->funcp = f;
		fi->ftype = ftype;
		fi->nargs = nargs;
		fi->funcinfo = funcinfo;
		if (!funcsfirst)
			funcsfirst = fi;
		else
			funcslast->fnext = fi;
		funcslast = fi;
		fi->fnext = 0;
		}
	}

enum { NEFB = 5, NEFB0 = 2 };

 static Exitcall a_e_info[NEFB0];
 static Exitcall *a_e_next = a_e_info;
 static Exitcall *a_e_last = a_e_info + NEFB0;
 static Exitcall *a_e_prev;

 static void
#ifdef KR_headers
AtReset(ae, ef, v) AmplExports *ae; Exitfunc *ef; char *v;
#else
AtReset(AmplExports *ae, Exitfunc *ef, void *v)
#endif
{
	Exitcall *ec;
	ASL *asl = (ASL*)ae->asl;
	if (asl->i.arnext >= asl->i.arlast) {
		asl->i.arnext = (Exitcall*)M1alloc(NEFB*sizeof(Exitcall));
		asl->i.arlast = asl->i.arnext + NEFB;
		}
	asl->i.arnext->prev = asl->i.arprev;
	asl->i.arprev = ec = asl->i.arnext++;
	ec->ef = ef;
	ec->v = v;
	}

 void
#ifdef KR_headers
at_end_ASL(ec) Exitcall *ec;
#else
at_end_ASL(Exitcall *ec)
#endif
{
	while(ec) {
		(*ec->ef)(ec->v);
		ec = ec->prev;
		}
	}

 void
at_exit_ASL(VOID)
{
	Exitcall *ec;
	ASLhead *h, *h0;

	h0 = &ASLhead_ASL;
	h = ASLhead_ASL.next;
	h0->next = h0->prev = h0;
	for(; h != h0; h = h->next)
		if (ec = ((ASL*)h)->i.arprev)
			at_end_ASL(ec);
	if (ec = a_e_prev) {
		a_e_prev = 0;
		at_end_ASL(ec);
		}
	}

 static void
#ifdef KR_headers
AtExit(ae, ef, v) AmplExports *ae; Exitfunc *ef; char *v;
#else
AtExit(AmplExports *ae, Exitfunc *ef, void *v)
#endif
{
	Exitcall *ec;
	Not_Used(ae);
#ifndef NO_ONEXIT
	if (!a_e_prev)
		atexit(at_exit_ASL); /* in case mainexit() is bypassed */
#endif
	if (a_e_next >= a_e_last) {
		a_e_next = (Exitcall*)mymalloc(NEFB*sizeof(Exitcall));
		a_e_last = a_e_next + NEFB;
		}
	a_e_next->prev = a_e_prev;
	a_e_prev = ec = a_e_next++;
	ec->ef = ef;
	ec->v = v;
	}

 struct
TMInfo {
	union {
		TMInfo *prev;
		double align;
		} u;
	};

 static Char *
#ifdef KR_headers
Tempmem(T, L) TMInfo *T; unsigned long L;
#else
Tempmem(TMInfo *T, unsigned long L)
#endif
{
	TMInfo *T1 = (TMInfo *)mymalloc(L + sizeof(TMInfo));
	T1->u.prev = T->u.prev;
	T->u.prev = T1;
	return (Char*)(T1+1);
	}

#ifdef KR_headers
 static void
No_table_handler(Dbread, Dbwrite, hname, flags, vinfo)
	int (*Dbread)(), (*Dbwrite)(), flags;
	char *hname; Char *vinfo;
{}

 static cryptblock*
No_crypto(key, scrbytes) char *key; Uint scrbytes;
#else
 static void
No_table_handler(
	int (*Dbread)(AmplExports*, TableInfo*),
	int (*Dbwrite)(AmplExports*, TableInfo*),
	char *hname,
	int flags,
	void *vinfo)
{}

 static cryptblock*
No_crypto(char *key, Uint scrbytes)
#endif
{ return 0; }

typedef void Funcadd ANSI((AmplExports*));

#ifndef CLOSE_AT_RESET
static Funcadd *Fa0[4], **Fa = Fa0;
static int nFa = 0, nFamax = 4;
#endif

#ifdef SYMANTEC
#define No_popen_or_pclose
typedef char *(*Tempnamtype)(const char*, const char*);
#define Tempnam_cast (Tempnamtype)
#endif

#ifdef WATCOM
#define tempnam _tempnam
#endif

#ifdef NO_tempnam

/* If the system does not provide a true tempnam function */
/* the AMPL/solver interface library will not do so either. */

 static char *
tempnam(const char *dir, const char *pfix)
{ return 0; }
#endif /* NO_tempnam */

#ifdef _WIN32
#define popen _popen
#define pclose _pclose
#else
#ifdef MSDOS
#undef  No_popen_or_pclose
#define No_popen_or_pclose
#endif
#endif

#ifdef No_popen_or_pclose
#undef popen
#define popen no_popen
#undef pclose
#define pclose no_pclose

 static int
no_pclose(FILE*f) { return 1; }

 static FILE*
no_popen(const char*cmd, const char*type) { return 0; }
#endif

 static AmplExports AE;

#ifdef clearerr
 static void
#ifdef KR_headers
myclearerr(f) FILE *f;
#else
myclearerr(FILE *f)
#endif
{ clearerr(f); }
#undef clearerr
#define clearerr myclearerr
#endif /*clearerr*/

#ifdef feof
 static int
#ifdef KR_headers
myfeof(f) FILE *f;
#else
myfeof(FILE *f)
#endif
{ return feof(f); }
#undef feof
#define feof myfeof
#endif /*feof*/

#ifdef ferror
 static int
#ifdef KR_headers
myferror(f) FILE *f;
#else
myferror(FILE *f)
#endif
{ return ferror(f); }
#undef ferror
#define ferror myferror
#endif /*ferror*/

#ifdef _fileno
#undef fileno
#define fileno _fileno
#endif

#ifdef fileno
 static int
#ifdef KR_headers
myfileno(f) FILE *f;
#else
myfileno(FILE *f)
#endif
{ return fileno(f); }
#undef fileno
#define fileno myfileno
#endif /* fileno */

#ifndef Tempnam_cast
#define Tempnam_cast /*nothing*/
#endif

 void
#ifdef KR_headers
func_add(asl) ASL *asl;
#else
func_add(ASL *asl)
#endif
{
	AmplExports *ae;

	if (need_funcadd) {
		if (!i_option_ASL
		 && !(i_option_ASL = getenv("ampl_funclibs")))
		      i_option_ASL = getenv("AMPLFUNC");
		if (!AE.PrintF) {
			AE.StdIn = stdin;
			AE.StdOut = stdout;
			AE.StdErr = Stderr;
			AE.ASLdate = ASLdate_ASL;
			AE.Addfunc = addfunc_ASL;
			AE.PrintF = printf;
			AE.FprintF = fprintf;
			AE.SprintF = sprintf;
			AE.VfprintF = vfprintf;
			AE.VsprintF = vsprintf;
			AE.Strtod = strtod;
			AE.AtExit = AtExit;
			AE.AtReset = AtReset;
			AE.Tempmem = Tempmem;
			AE.Add_table_handler = No_table_handler;
			AE.Crypto = No_crypto;
			AE.Qsortv = qsortv;
			AE.Clearerr = clearerr;
			AE.Fclose = fclose;
			AE.Fdopen = fdopen;
			AE.Feof = feof;
			AE.Ferror = ferror;
			AE.Fflush = fflush;
			AE.Fgetc = fgetc;
			AE.Fgets = fgets;
			AE.Fileno = fileno;
			AE.Fopen = fopen;
			AE.Fputc = fputc;
			AE.Fputs = fputs;
			AE.Fread = fread;
			AE.Freopen = freopen;
			AE.Fscanf = fscanf;
			AE.Fseek = fseek;
			AE.Ftell = ftell;
			AE.Fwrite = fwrite;
			AE.Pclose = pclose;
			AE.Perror = perror;
			AE.Popen = popen;
			AE.Puts = puts;
			AE.Rewind = rewind;
			AE.Scanf = scanf;
			AE.Setbuf = setbuf;
			AE.Setvbuf = setvbuf;
			AE.Sscanf = sscanf;
			AE.Tempnam = Tempnam_cast tempnam;
			AE.Tmpfile = tmpfile;
			AE.Tmpnam = tmpnam;
			AE.Ungetc = ungetc;
			AE.Getenv = getenv_ASL;
			}
		if (AE.asl)
			memcpy(ae = (AmplExports*)M1alloc(sizeof(AmplExports)),
				&AE, sizeof(AmplExports));
		else
			ae = &AE;
		asl->i.ae = ae;
		ae->asl = (Char*)asl;
		auxinfo_ASL(ae);
#ifndef CLOSE_AT_RESET
		if (nFa > 0) {
			/* not the first nl_reader call */
			int i = 0;
			while(i < nFa)
				(*Fa[i++])(ae);
			}
		else
#endif
			funcadd(ae);
		need_funcadd = 0;
		}
	}

 void
#ifdef KR_headers
show_funcs_ASL(asl) ASL *asl;
#else
show_funcs_ASL(ASL *asl)
#endif
{
	func_info *fi;
	int nargs;
	char *atleast;

	func_add(asl);
	fprintf(Stderr, "Available nonstandard functions:%s\n",
		(fi = funcsfirst) ? "" : " none");
	for(; fi; fi = fi->fnext) {
		if ((nargs = fi->nargs) >= 0)
			atleast = "";
		else {
			nargs = -(1 + nargs);
			atleast = "at least ";
			}
		fprintf(Stderr, "\t%s(%s%d %sarg%s)\n", fi->name,
			atleast, nargs, fi->ftype ? "" : "real ",
			nargs == 1 ? "" : "s");
		}
	fflush(Stderr);
	}

 int
#ifdef KR_headers
aflibname_ASL(ae, fullname, name, nlen, fa, save_fa)
	AmplExports *ae; char *fullname; char *name; int nlen; Funcadd *fa;
	int save_fa;
#else
aflibname_ASL(AmplExports *ae, char *fullname, char *name, int nlen,
	Funcadd *fa, int save_fa)
#endif
{
	Not_Used(fullname);
	Not_Used(name);
	Not_Used(nlen);
	n_added = 0;
	(*fa)(ae);
#ifndef CLOSE_AT_RESET
	if (n_added && save_fa) {
		if (++nFa >= nFamax) {
			Funcadd **Fa1;
			nFamax <<= 1;
			Fa1 = (Funcadd**)Malloc(nFamax * sizeof(Funcadd*));
			memcpy(Fa1, Fa, nFa*sizeof(Funcadd*));
			if (Fa != Fa0)
				free(Fa);
			Fa = Fa1;
			}
		Fa[nFa-1] = fa;
		}
#endif /* !CLOSE_AT_RESET */
	return n_added;
	}
#ifdef __cplusplus
}
#endif
/* Last relevant change to asl.h: 19991013. */
