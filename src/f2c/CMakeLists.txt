cmake_minimum_required(VERSION 3.3)
project(f2c)


set(MISC 
f77vers.c	i77vers.c	s_rnge.c	abort_.c
exit_.c		getarg_.c	iargc_.c	getenv_.c
signal_.c	s_paus.c	s_stop.c	system_.c
cabs.c		derf_.c		derfc_.c	erf_.c
erfc_.c		sig_die.c	uninit.c  main.c)

if(NOT MSVC)
	set(MISC ${MISC} ctype.c)
endif()

set(POW 
pow_ci.c	pow_dd.c	pow_di.c	pow_hh.c
pow_ii.c	pow_ri.c	pow_zi.c	pow_zz.c)

set(CX 
c_abs.c 	c_cos.c		c_div.c 	c_exp.c 	
c_log.c		c_sin.c		c_sqrt.c)

set(DCX 
z_abs.c		z_cos.c		z_div.c		z_exp.c
z_log.c		z_sin.c		z_sqrt.c)

set(REAL 
r_abs.c		r_acos.c	r_asin.c	r_atan.c	
r_atn2.c	r_cnjg.c	r_cos.c		r_cosh.c	
r_dim.c		r_exp.c		r_imag.c	r_int.c		
r_lg10.c	r_log.c		r_mod.c		r_nint.c	
r_sign.c	r_sin.c		r_sinh.c	r_sqrt.c	
r_tan.c		r_tanh.c)

set(DBL 
d_abs.c 	d_acos.c 	d_asin.c 	d_atan.c
d_atn2.c	d_cnjg.c	d_cos.c		d_cosh.c	
d_dim.c		d_exp.c		d_imag.c	d_int.c		
d_lg10.c	d_log.c		d_mod.c		d_nint.c	
d_prod.c	d_sign.c	d_sin.c		d_sinh.c	
d_sqrt.c	d_tan.c		d_tanh.c)

set(INT 
i_abs.c		i_dim.c		i_dnnt.c	i_indx.c
i_len.c		i_mod.c		i_nint.c	i_sign.c
lbitbits.c	lbitshft.c)

set(HALF
h_abs.c		h_dim.c		h_dnnt.c	h_indx.c	
h_len.c		h_mod.c		h_nint.c	h_sign.c)

set(CMP 
l_ge.c		l_gt.c		l_le.c		l_lt.c
hl_ge.c		hl_gt.c		hl_le.c		hl_lt.c)

set(EFL
ef1asc_.c	ef1cmc_.c)

set(CHAR
f77_aloc.c	s_cat.c		s_cmp.c	s_copy.c)

set(QINT
backspac.c	close.c		dfe.c	dolio.c)

set(TIME
dtime_.c etime_.c)

set(I77
close.c		due.c		endfile.c	err.c		
fmt.c		fmtlib.c	ftell_.c	iio.c
ilnw.c		inquire.c	lread.c		lwrite.c
open.c		rdfmt.c		rewind.c	rsfe.c
rsli.c		rsne.c		sfe.c		sue.c
typesize.c	uio.c		util.c		wref.c
wrtfmt.c	wsfe.c		wsle.c		wsne.c
xwsne.c   dolio.c)

set(QINT
pow_qq.c qbitbits.c qbitshft.c ftell64_.c)

# For INTEGER*8 support (which requires system-dependent adjustments to
# f2c.h), add $(QINT) to the F2C_SOURCES assignment below...
set(F2C_SOURCES ${MISC} ${POW} ${CX} ${DCX} ${REAL} 
${DBL} ${INT}	${HALF} ${CMP} ${EFL} ${CHAR} ${I77} ${TIME})

  # NOTE main.c defines the main function, that then gives 
  # multiple defined symbol error when linking to solvers
  # I have modified it and placed it outside the src directory
  # waiting for a workaround


# Generation of f2c.h
# The following implements what is written in the README file
if(NOT MSVC)
  if(${ARCH} MATCHES 64)
    set(CMAKE_REQUIRED_FLAGS "-m64")
  else()
    set(CMAKE_REQUIRED_FLAGS "-m32")
  endif()
endif()
include(CheckTypeSize)
check_type_size("int" INT_SIZE)
check_type_size("long" LONG_SIZE)
check_type_size("double" DOUBLE_SIZE)
if(INT_SIZE EQUAL 4 AND LONG_SIZE EQUAL 8 AND DOUBLE_SIZE EQUAL 8) 
  file(READ f2c.h0 f2c_header)
  string(REGEX REPLACE "long " "" f2c_header "${f2c_header}")
  file(WRITE ${GENERATED_INCLUDE_DIR}/f2c.h "${f2c_header}")
else()
  configure_file(f2c.h0 ${GENERATED_INCLUDE_DIR}/f2c.h)
endif()

# Note TODO: generating math.h in the include dir causes infinite
# recursion in includes
configure_file(math.hvc math.h)
configure_file(signal1.h0 ${GENERATED_INCLUDE_DIR}/signal1.h)
configure_file(sysdep1.h0 ${GENERATED_INCLUDE_DIR}/sysdep1.h)

add_library(f2c ${F2C_SOURCES})
add_dependencies(f2c arith-h) # To make sure that we generate arith.h
setarchitectureflags(${ARCH} f2c)
target_include_directories(f2c PUBLIC "${GENERATED_INCLUDE_DIR}")

if(MSVC)
  target_compile_definitions(f2c PRIVATE USE_CLOCK MSDOS NO_ONEXIT NO_My_ctype NO_ISATTY)
  target_compile_options(f2c PRIVATE /MT$<$<CONFIG:Debug>:d> /wd4244 /wd4996 /wd4554 /wd4267 /wd4244)
elseif (MINGW)
  add_definitions(f2c PRIVATE USE_CLOCK MSDOS NO_My_ctype)
elseif (LINUX)

else()
endif()


