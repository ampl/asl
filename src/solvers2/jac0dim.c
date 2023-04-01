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

#include "asl.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void fpinit_ASL(void);

 extern int Sscanf(char*, const char*, ...);

#undef Want_bswap
#ifdef IEEE_MC68k
#define Want_bswap
#endif
#ifdef IEEE_8087
#define Want_bswap
#endif
jmp_buf Jb;
#ifdef Want_bswap

 void
bswap_ASL(void *x, size_t L)
{
	char *s = (char*)x;
	int t;
	switch(L) {
	 case 2:
		t = s[0]; s[0] = s[1]; s[1] = t;
		break;
	 case 4:
		t = s[0]; s[0] = s[3]; s[3] = t;
		t = s[1]; s[1] = s[2]; s[2] = t;
		break;
	 case 8:
		t = s[0]; s[0] = s[7]; s[7] = t;
		t = s[1]; s[1] = s[6]; s[6] = t;
		t = s[2]; s[2] = s[5]; s[5] = t;
		t = s[3]; s[3] = s[4]; s[4] = t;
	 }
	}

#endif /* Want_bswap */

 static void
badfmt(EdRead *R)
{
	if (R->return_nofile & 2)
		longjmp(Jb, 1);
	badread(R);
	fprintf(Stderr, "Unrecognized binary format.\n");
	exit(1);
	}

 static void
badints(EdRead *R, int got, int wanted)
{
	if (R->return_nofile & 2)
		longjmp(Jb, 1);
	badread(R);
	fprintf(Stderr, "got only %d integers; wanted %d\n", got, wanted);
	exit(1);
	}

 static void
read2(EdRead *R, int *x, int *y)
{
	char *s;
	int k;

	s = read_line(R);
	k = Sscanf(s, " %d %d", x, y);
	if (k != 2)
		badints(R, k, 2);
	}

#define ndcc asl->i.ndcc_
#define nzlb asl->i.nzlb_

 FILE *
jac0dim_ASL(ASL *asl, const char *stub, ftnlen stub_len)
{
	FILE *nl;
	int i, k, nlv;
	char *s, *se;
	const char *opfmt;
	EdRead ER, *R;

	if (!asl)
		badasl_ASL(asl,0,"jac0dim");
	fpinit_ASL();	/* get IEEE arithmetic, if possible */

	if (stub_len <= 0)
		for(i = 0; stub[i]; i++);
	else
		for(i = stub_len; stub[i-1] == ' ' && i > 0; --i);
	filename = (char *)M1alloc(i + 5);
	s = stub_end = filename + i;
	strncpy(filename, stub, i);
	strcpy(s, ".nl");
	nl = fopen(filename, "rb");
	if (!nl && i > 3 && !strncmp(s-3, ".nl", 3)) {
		*s = 0;
		stub_end = s - 3;
		nl = fopen(filename, "rb");
		}
	if (!nl) {
		if (return_nofile & 1)
			return 0;
		fflush(stdout);
		what_prog();
		fprintf(Stderr, "can't open %s\n", filename);
		exit(1);
		}
	if (return_nofile & 2) {
		if (setjmp(Jb))
			return 0;
		}
	R = EdReadInit_ASL(&ER, asl, nl, 0);
	R->Line = 0;
	s = read_line(R);
	binary_nl = 0;
	opfmt = "%d";
	switch(*s) {
#ifdef DEPRECATED
		case 'E':	/* deprecated "-oe" format */
			{int ncsi = 0;
			k = Sscanf(s, "E%d %d %d %d %d %d", &n_var, &n_con,
				&n_obj, &maxrownamelen, &maxcolnamelen, &ncsi);
			if (k < 5)
				badints(R, k, 5);
			if (ncsi) {
				if (ncsi != 6) {
					if (return_file & 2)
						return 0;
					badread(R);
					fprintf(Stderr,
					 "expected 6th integer to be 0 or 6, not %d\n",
						ncsi);
					exit(1);
					}
				s = read_line(R);
				k = Sscanf(s, " %d %d %d %d %d %d",
					&comb, &comc, &como, &comc1, &como1, &nfunc);
				if (k != 6)
					badints(R, k, 6);
				combc = comb + comc;
				}
			}
			break;
#endif
		case 'z':
		case 'Z':
			opfmt = "%hd";
		case 'B':
		case 'b':
			binary_nl = 1;
			xscanf = bscanf;
			goto have_xscanf;
		case 'h':
		case 'H':
			opfmt = "%hd";
			binary_nl = 1;
			xscanf = hscanf;
			goto have_xscanf;
		case 'G':
		case 'g':
			xscanf = ascanf;
 have_xscanf:
			if ((k = ampl_options[0] = strtol(++s, &se, 10))) {
				if (k > 9) {
					fprintf(Stderr,
					"ampl_options = %d is too large\n", k);
					exit(1);
					}
				for(i = 1; i <= k && se > s; i++)
					ampl_options[i] = strtol(s = se,&se,10);
				if (ampl_options[2] == 3)
					ampl_vbtol = strtod(s = se, &se);
				}
			s = read_line(R);
			n_eqn = -1;
			k = Sscanf(s, " %d %d %d %d %d %d", &n_var, &n_con,
				&n_obj, &nranges, &n_eqn, &n_lcon);
			if (k < 3)
				badints(R,k,3);
			nclcon = n_con + n_lcon;

			/* formerly read2(R, &nlc, &nlo); */
			s = read_line(R);
			n_cc = nlcc = ndcc = nzlb = 0;
			k = Sscanf(s, " %d %d %d %d %d %d", &nlc, &nlo, &n_cc, &nlcc,
					&ndcc, &nzlb);
			if (k < 2)
				badints(R,k,2);
			asl->i.nlc0 = nlc;
			asl->i.nlo0 = nlo;
			if ((n_cc += nlcc) > 0 && k < 6)
				ndcc = -1; /* indicate unknown */

			read2(R, &nlnc, &lnc);
			nlvb = -1;
			s = read_line(R);
			k = Sscanf(s, " %d %d %d", &nlvc, &nlvo, &nlvb);
			if (k < 2)
				badints(R,k,2);

			/* read2(R, &nwv, &nfunc); */
			s = read_line(R);
			asl->i.flags = 0;
			k = Sscanf(s, " %d %d %d %d", &nwv, &nfunc, &i,
				&asl->i.flags);
			if (k < 2)
				badints(R,k,2);
			else if (k >= 3 && i != Arith_Kind_ASL && i) {
#ifdef Want_bswap
				if (i > 0 && i + Arith_Kind_ASL == 3) {
					asl->i.iadjfcn = asl->i.dadjfcn = bswap_ASL;
					binary_nl = i << 1;
					}
				else
#endif
					badfmt(R);
				}

			if (nlvb < 0)	/* ampl versions < 19930630 */
				read2(R, &nbv, &niv);
			else {
				s = read_line(R);
				k = Sscanf(s, " %d %d %d %d %d", &nbv, &niv,
					&nlvbi, &nlvci, &nlvoi);
				if (k != 5)
					badints(R,k,5);
				}
			/* formerly read2(R, &nzc, &nzo); */
			s = read_line(R);
			k = Sscanf(s, " %zd %zd", &nZc, &nZo);
			if (k != 2)
				badints(R, k, 2);
			nzc = nZc;
			nzo = nZo;
			read2(R, &maxrownamelen, &maxcolnamelen);
			s = read_line(R);
			k = Sscanf(s, " %d %d %d %d %d", &comb, &comc, &como,
					&comc1, &como1);
			if (k != 5)
				badints(R,k,5);
			combc = comb + comc;
		}
	student_check_ASL(asl);
	if (n_con < 0 || n_var <= 0 || n_obj < 0) {
		what_prog();
		fprintf(Stderr,
		"jacdim: got M = %d, N = %d, NO = %d\n", n_con, n_var, n_obj);
		exit(1);
		}
	asl->i.opfmt = opfmt;
	asl->i.n_var0 = asl->i.n_var1 = n_var;
	asl->i.n_con0 = asl->i.n_con1 = n_con;
	if ((nlv = nlvc) < nlvo)
		nlv = nlvo;
	x0len = nlv * sizeof(real);
	n_conjac[0] = 0;
	n_conjac[1] = n_con;
	c_vars = o_vars = n_var;	/* confusion arises otherwise */
	return nl;
	}

 FILE *
jac_dim_ASL(ASL *asl, const char *stub, fint *M, fint *N, fint *NO, fint *NZ,
		fint *MXROW, fint *MXCOL, fint stub_len)
{
	FILE *nl;

	if ((nl = jac0dim_ASL(asl, stub, stub_len))) {
		*M = n_con;
		*N = n_var;
		*NO = n_obj;
		*NZ = nzc;
		*MXROW = maxrownamelen;
		*MXCOL = maxcolnamelen;
		}
	return nl;
	}

#ifdef __cplusplus
	}
#endif
