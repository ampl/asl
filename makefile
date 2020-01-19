# /****************************************************************
# Copyright (C) 1997 Lucent Technologies
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that both that the copyright notice and this
# permission notice and warranty disclaimer appear in supporting
# documentation, and that the name of Lucent or any of its entities
# not be used in advertising or publicity pertaining to
# distribution of the software without specific, written prior
# permission.
#
# LUCENT DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.
# IN NO EVENT SHALL LUCENT OR ANY OF ITS ENTITIES BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
# IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,
# ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
# THIS SOFTWARE.
# ****************************************************************/

.SUFFIXES: .c .o
CC = cc
CFLAGS = -O
SHELL=/bin/sh

# Add -DKR_headers to CFLAGS if your C compiler does not
# understand ANSI C function headers, e.g.
#	CFLAGS = -O -DKR_headers
# You may also need to add
#	strerror.c \
# to the "a =" assignment below (if you get an error message about
# strerror not being found when you try to link a solver).
# If things don't run right, you may need to change -O to -g
# so you can poke around with a debugger.
# You may also need to add -D_SIZE_T to CFLAGS, or to
# comment out the definition of size_t in nlp.h .

# For SunOS systems, try
# CFLAGS = -O -DKR_headers -D_SIZE_T

# For the DEC Alpha, try
# CFLAGS = -g -ieee_with_no_inexact

# For HP, try
# CFLAGS = -O -Aa
# or (if the compiler does not recognize -Aa for ANSI syntax)
# CFLAGS = -O -DKR_headers

# For IBM RS6000 machines, add
#	-qnomaf
# to CFLAGS (to avoid surprises from fused mutiply-add instructions).

.c.o:
	$(CC) -c $(CFLAGS) $*.c

all: arith.h stdio1.h amplsolver.a funcadd0.o

a = \
	asldate.c \
	atof.c \
	b_search.c \
	basename.c \
	bscanf.c \
	com2eval.c \
	comeval.c \
	con1ival.c \
	con2ival.c \
	con2val.c \
	conadj.c \
	conpval.c \
	conscale.c \
	conval.c \
	der0prop.c \
	derprop.c \
	dtoa1.c \
	duthes.c \
	dynlink.c \
	f_read.c \
	fg_read.c \
	fgh_read.c \
	fpecatch.c \
	fullhes.c \
	func_add.c \
	g_fmt.c \
	getenv.c \
	getstub.c \
	htcl.c \
	jac0dim.c \
	jacdim.c \
	jacinc.c \
	mach.c \
	mainexit.c \
	mip_pri.c \
	misc.c \
	mypow.c \
	names.c \
	nl_obj.c \
	nqpcheck.c \
	obj2val.c \
	obj_prec.c \
	objconst.c \
	objval.c \
	objval_.c \
	op_type.c \
	pfg_read.c \
	pfghread.c \
	printf.c \
	pshvprod.c \
	punknown.c \
	qp_read.c \
	qpcheck.c \
	qsortv.c \
	readsol.c \
	repwhere.c \
	rops.c \
	rops2.c \
	sphes.c \
	sscanf.c \
	stderr.c \
	value.c \
	writesol.c \
	wrtsol_.c \
	ws_desc.c \
	wsu_desc.c \
	x2check.c \
	xectim.c \
	xp1known.c \
	xp2known.c

amplsolver.a: $a
	$(CC) -c $(CFLAGS) $?
	x=`echo $? | sed 's/\.c/.o/g'` && ar ruv amplsolver.a $$x && rm $$x
	ranlib amplsolver.a || true
# If your system lacks ranlib, add a dummy ranlib to your
# search path, e.g.
#	exec true
# or just comment out the ranlib invocation above.

Aslh = asl.h stdio1.h
dtoa1.o ed2read.o edagread.o g_fmt.o mach.o pfg_read.o pfghread.o: arith.h
dtoa1.o: dtoa.c
conadj.o conscale.o conval.o con1ival.o jacdim.o jacinc.o names.o nqpcheck.o\
 objval.o qpcheck.o wrtsol_.o: nlp.h $(Aslh)
com2eval.o rops2.o: nlp2.h $(Aslh)
comeval.o ed0read.o edagread.o qp_read.o mip_pri.o rops.o: nlp.h $(Aslh)
con2ival.o con2val.o ed2read.o obj2val.o x2check.o: jac2dim.h nlp2.h $(Aslh)
conpval.o xp2known.o: jacpdim.h asl_pfg.h nlp2.h $(Aslh)
der0prop.o derprop.o objval_.o repwhere.o: $(Aslh)
duthes.o pfg_read.o pfghread.o: asl_pfg.h nlp2.h $(Aslh)
dynlink.o: funcadd.h
ed2read.o pfg_read.o: dvalue.hd op_type.hd opnos.hd
ed0read.o edagread.o: nlp.h $(Aslh)
ed0read.o: edagread.c
edagread.o qp_read.o: op_type.hd dvalue.hd r_opn.hd
funcadd.o: funcadd.h
func_add.o getstub.o misc.o value.o: $(Aslh)
objconst.o: $(Aslh) nlp.h nlp2.h asl_pfg.h
op_type.o: op_type.hd
pfg_read.o pfghread.o: r_opn0.hd
pfghread.o: pfg_read.c
readsol.o writesol.o: nlp.h $(Aslh)
rops.o rops2.o: r_op.hd
nqpcheck.o qpcheck.o: r_qp.hd
printf.o: stdio1.h
xp1known.o: asl_pfg.h nlp.h $(Aslh)

# Use CFLAGS in compiling arithchk.c in case something in CFLAGS affects
# the number of bits in integral data types.  (It's probably best not to
# add such options to CFLAGS.)

arith.h: arithchk.c
	cc $(CFLAGS) arithchk.c
	a.out >arith.h
	rm -f a.out arithchk.o

### Alternative to arithchk.c: copy arith.h0 to arith.h, then edit
### arith.h to activate the appropriate #define line, as explained
### in the comments at the top.  For systems with IBM-mainframe
### arithmetic, see README.  You'll need to #define Arith_Kind_ASL
### suitably.  If "make arith.h" works, use the #define it gives
### (and don't fool with arith.h0).  Otherwise use
###	#define Arith_Kind_ASL 0

# If compiling dtoa1.c reveals that your system lacks float.h, malloc.h
# or memory.h, you could try
#
#	  make float.h
#
#         make malloc.h
# and/or
#         make memory.h
#
# as appropriate.

float.h: float.h0
	ln float.h0 float.h

malloc.h:
	echo 'extern char *malloc();' >malloc.h

memory.h:
	echo 'extern char *memcpy();' >memory.h

stdio1.h: stdio1.h0
	cat stdio1.h0 >stdio1.h

### The rule above arranges for amplsolver.a to use printf, fprintf,
### and sprintf (renamed Printf, Fprintf, and Sprintf) as described
### in the comments at the start of printf.c, rather than the
### system-supplied routines (whose sprintf has the wrong return
### type and values on some systems).  In your solver, say
### #include "stdio1.h" rather than <stdio.h> for consistency
### with amplsolver.a .  To use the system-supplied printf (etc.),
### say "make systemprintf", and change printf.c to sprintf.c
### in the "a =" assignment above.

systemprintf:
	echo '#define NO_STDIO1' >stdio1.h
	cat stdio1.h0 >>stdio1.h

# make xsum.out to check for transmission errors.
# This assumes you have the xsum program, whose source
# you can get by asking research!netlib to
#	send xsum.c from f2c/src

xs0 = \
	README \
	README.f77 \
	amplsolv.lbc \
	amplsolv.sy \
	arith.ibm \
	arith.h0 \
	arithchk.c \
	asl.h \
	asl_pfg.h \
	asl_pfgh.h \
	asldate.c \
	atof.c \
	b_search.c \
	basename.c \
	bscanf.c \
	com2eval.c \
	comeval.c \
	con1ival.c \
	con2ival.c \
	con2val.c \
	conadj.c \
	conpval.c \
	conscale.c \
	conval.c \
	der0prop.c \
	derprop.c \
	dtoa.c \
	dtoa1.c \
	duthes.c \
	dvalue.hd \
	dynlink.c \
	f_read.c \
	fg_read.c \
	fgh_read.c \
	float.h0 \
	fpecatch.c \
	fullhes.c \
	func_add.c \
	funcadd.c \
	funcadd.h \
	funcadd0.c \
	funcaddk.c \
	g_fmt.c \
	getenv.c \
	getstub.c \
	getstub.h \
	htcl.c \
	jac0dim.c \
	jac2dim.h \
	jacdim.c \
	jacinc.c \
	jacpdim.h \
	mach.c \
	mainexit.c \
	makefile \
	makefile.sy \
	makefile.vc \
	makefile.wat \
	mip_pri.c \
	misc.c \
	mypow.c \
	names.c \
	nl_obj.c \
	nlp.h \
	nlp2.h \
	nqpcheck.c \
	obj2val.c \
	obj_prec.c \
	objconst.c \
	objval.c \
	objval_.c \
	op_type.c \
	op_type.hd \
	opnos.hd \
	pfg_read.c \
	pfghread.c \
	printf.c \
	pshvprod.c \
	psinfo.h \
	punknown.c \
	qp_read.c \
	qpcheck.c \
	qsortv.c \
	r_op.hd \
	r_opn.hd \
	r_opn0.hd \
	r_qp.hd \
	readsol.c \
	repwhere.c \
	rnd_prod.s \
	rops.c \
	rops2.c \
	sjac0dim.c \
	sphes.c \
	sprintf.c \
	sscanf.c \
	stderr.c \
	stdio1.h0 \
	strerror.c \
	value.c \
	writesol.c \
	wrtsol_.c \
	ws_desc.c \
	wsu_desc.c \
	x2check.c \
	xectim.c \
	xp1known.c \
	xp2known.c

xsum.out: xsum0.out $(xs0)
	xsum $(xs0) >xsum1.out
	cmp xsum0.out xsum1.out && mv xsum1.out xsum.out || diff xsum[01].out
