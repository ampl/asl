# This cmake file builds the libraries asl and asl2.
# It also builds asl-mt and asl2-mt, that are the same library with appropriate
# compiler flags to support multiple threads.
# On windows it also builds asl-dynrt, asl2-dynrt, aslcpp-dynrt, built against the
# dynamic runtime (they are used to compile mp, for example).
# It exports the static libraries as targets with the corresponding name.
# Define GENERATE_ARITH to force old-style generation of the arith.h file, useful for
# non-x86 CPUs. No cross-compiling when this is enabled, as it sets the CPU arithmetics 
# constants to the machine in which is run.
# The following options (default off) control what is built:
# BUILD_CPP to build the cpp wrapper to the ASL library
# BUILD_SHARED_LIBS to build dynamic libraries instead of static
# BUILD_DYNRT_LIBS to build the libraries linked against the dynamic runtime
# BUILD_MT_LIBS to build the multithreaded libraries 
# BUILD_F2C to build the fortran to c converter (https://www.netlib.org/f2c/f2c.pdf)
# BUILD_EXAMPLES to build the examples

if (${CMAKE_VERSION} VERSION_GREATER "3.13.0")
  cmake_policy(SET CMP0077 NEW)
endif()
if (${CMAKE_VERSION} VERSION_GREATER "3.15.0")
  cmake_policy(SET CMP0091 NEW)
endif()
cmake_policy(SET CMP0069 NEW) 
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
project(ASL)

cmake_minimum_required(VERSION 3.0)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
# Set the path to CMake modules.
set(AMPL_CMAKE_MODULE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/support/cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${AMPL_CMAKE_MODULE_DIR})

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(BUILD_MCMODELLARGE "Build the library linked with mcmodel=large flag, used for compiling lgo" OFF)
option(BUILD_LICCHECK_PRINT "Build libraries with licence check flag" OFF)
option(BUILD_DYNRT_LIBS "Build the libraries linked against the dynamic runtime (windows only)" OFF)
option(BUILD_MT_LIBS "Build the multithreaded library" OFF)
option(BUILD_CPP "Build the cpp interface" OFF)
option(BUILD_F2C "Build the f2c library" OFF)
option(USE_LTO "Use link time optimization" ON)

if(NOT BUILD_ASL_EXAMPLES)
  option(BUILD_ASL_EXAMPLES "Build the examples" OFF)
else()
  option(BUILD_F2C "Build the f2c library" ON)
  option(BUILD_ASL_EXAMPLES "Build the examples" ON)
  set(BUILD_F2C ON)
endif()

# Base directory
set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(GENERATED_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include)

if(USE_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT isSupported)
  if(isSupported)
    message("INTERPROCEDURAL OPTIMIZATION is supported, enabled by default")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  endif()
endif()

# Test if I have a parent project, in which case put generated targets 
# in the asl folder on supported generators
get_directory_property(hasParent PARENT_DIRECTORY)

function(addToASLFolder folder)
    if(hasParent)
        set(prefix asl/)
    else()
        set(prefix "")
    endif()
    add_to_folder("${prefix}${folder}" ${ARGN})
endfunction()

# Set output directories.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Get rid of useless warnings
if(MSVC)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

# For addPrefix and add_to_folder
include(addPrefix)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
# Define the ARCH variable (32 or 64)
include(initArchitecture)
# Set generic architecture compiler flags
include(setArchitecture)
getArchitectureFlags(${ARCH} CCFLAGS LLFLAGS)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CCFLAGS}")

# ################# Generate dynamic files ##################
file(READ ${SRCDIR}/solvers/details.c0 DETAILS)
string(REPLACE "System_details"
               "${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}" DETAILS
               "${DETAILS}")
file(WRITE ${GENERATED_INCLUDE_DIR}/details.c "${DETAILS}")

configure_file(${SRCDIR}/solvers/stdio1.h0 ${GENERATED_INCLUDE_DIR}/stdio1.h
               COPYONLY)

if(CMAKE_CROSSCOMPILING AND CMAKE_SYSTEM_PROCESSOR MATCHES "^x86")
  include(CheckTypeSize)
  check_type_size(double DOUBLE_SIZE)
  check_type_size(long LONG_SIZE)
  set(ARITH_H "#define IEEE_8087\n#define Arith_Kind_ASL 1\n")
  set(ARITH_INT) # Integer type used by arith.h checks.
  math(EXPR LONG_SIZE_X2 "${LONG_SIZE} * 2")
  if(DOUBLE_SIZE EQUAL LONG_SIZE_X2)
    set(ARITH_INT long)
    set(ARITH_INT_SIZE ${LONG_SIZE})
  else()
    check_type_size(int INT_SIZE)
    math(EXPR INT_SIZE_X2 "${INT_SIZE} * 2")
    if(DOUBLE_SIZE EQUAL INT_SIZE_X2)
      set(ARITH_H "${ARITH_H}#define Long int\n#define Intcast (int)(long)\n")
      set(ARITH_INT int)
      set(ARITH_INT_SIZE ${INT_SIZE})
    endif()
  endif()
  if(ARITH_INT)
    check_type_size("struct { double d\; ${ARITH_INT} L\; }[2]" STRUCT_SIZE)
    math(EXPR DOUBLE_PLUS_INT_SIZE_X2
         "2 * (${DOUBLE_SIZE} + ${ARITH_INT_SIZE})")
    if(STRUCT_SIZE GREATER DOUBLE_PLUS_INT_SIZE_X2)
      set(ARITH_H "${ARITH_H}#define Double_Align\n")
    endif()
  endif()
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(ARITH_H "${ARITH_H}#define X64_bit_pointers\n")
  endif()
  # Check long long.
  check_type_size("long long" LONG_LONG_SIZE)
  if(LONG_LONG_SIZE GREATER LONG_SIZE AND LONG_LONG_SIZE EQUAL
                                          CMAKE_SIZEOF_VOID_P)
    set(ARITH_H "${ARITH_H}#define LONG_LONG_POINTERS\n")
  endif()
  if(LONG_LONG_SIZE LESS 8)
    set(ARITH_H "${ARITH_H}#define NO_LONG_LONG\n")
  endif()
  # Check ssize_t.
  check_type_size(size_t SIZE_T_SIZE)
  check_type_size(ssize_t SSIZE_T_SIZE)
  if(NOT SSIZE_T_SIZE)
    if(SIZE_T_SIZE EQUAL LONG_SIZE)
      set(ARITH_SSIZE_T long)
    elseif(SIZE_T_SIZE EQUAL INT_SIZE)
      set(ARITH_SSIZE_T int)
    elseif(SIZE_T_SIZE EQUAL LONG_LONG_SIZE)
      set(ARITH_SSIZE_T "long long")
    else()
      set(ARITH_SSIZE_T "signed size_t")
    endif()
    set(ARITH_H "${ARITH_H}#define ssize_t ${ARITH_SSIZE_T}\n")
  elseif(NOT SIZE_T_SIZE EQUAL SSIZE_T_SIZE)
    set(ARITH_H "${ARITH_H}/* sizeof(size_t) = ${SIZE_T_SIZE}")
    set(ARITH_H "${ARITH_H} but sizeof(ssize_t) = ${SSIZE_T_SIZE} */\n")
  endif()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/arith.h
       "${ARITH_H}#define QNaN0 0x0\n#define QNaN1 0xfff80000\n")
else()
  if(NOT WIN32)
    set(FPINIT ${SRCDIR}/solvers/fpinit.c)
  endif()
  if(GENERATE_ARITH)
    add_executable(arithchk ${SRCDIR}/solvers/arithchk.c ${FPINIT})
    if(WIN32)
      target_compile_definitions(arithchk PRIVATE NO_FPINIT NO_SSIZE_T)
    else()
      target_compile_definitions(arithchk PRIVATE ASL_NO_FPINITMT)
      target_link_libraries(arithchk m)
    endif()
    add_custom_command(
      OUTPUT ${GENERATED_INCLUDE_DIR}/arith.h
      DEPENDS arithchk
      COMMAND ${WINE} $<TARGET_FILE:arithchk> > ${GENERATED_INCLUDE_DIR}/arith.h
      COMMENT "Writing arith.h to ${GENERATED_INCLUDE_DIR}/arith.h")
    if(NOT MSVC)
      string(REPLACE "-Wall -Wextra -pedantic" "" CMAKE_C_FLAGS
                     "${CMAKE_C_FLAGS}")
      target_compile_options(arithchk PUBLIC -Wno-format-security)
    endif()
  else()
    add_custom_command(
      OUTPUT ${GENERATED_INCLUDE_DIR}/arith.h
      COMMAND ${CMAKE_COMMAND} -E copy ${SRCDIR}/solvers/arith.h1
              ${GENERATED_INCLUDE_DIR}/arith.h
      COMMENT "Copying arith.h1 to ${GENERATED_INCLUDE_DIR}/arith.h")
  endif(GENERATE_ARITH)
endif()

# Use a custom target for arith.h, because including a file generated by a
# custom command in more than one target may result in it being generated
# multiple times in a parallel build which doesn't work in msbuild. See:
# http://www.cmake.org/pipermail/cmake/2008-October/024492.html
add_custom_target(arith-h DEPENDS ${GENERATED_INCLUDE_DIR}/arith.h)
addToASLFolder(support arith-h arithchk)

# Extract version information from the content of asldate.c
file(READ ${SRCDIR}/solvers/asldate.c ASL_VER)
string(REGEX MATCH "= ([0-9]*)" _ ${ASL_VER})
# Check if the version information was found
if(CMAKE_MATCH_1)
    set(ASL_VERSION ${CMAKE_MATCH_1})
    if(hasParent)
    set(ASL_VERSION ${CMAKE_MATCH_1} PARENT_SCOPE)
    endif()
    message(STATUS "Found ASL version: ${ASL_VERSION}")
else()
    message(WARNING "Failed to extract ASL version from asldate.c")
endif()

# ASL 1 sources
set(ASL_SOURCE_DIR ${SRCDIR}/solvers)
set(ASL_SOURCES ${GENERATED_INCLUDE_DIR}/details.c)

set(ASL_HEADERS ${GENERATED_INCLUDE_DIR}/arith.h)
add_prefix(ASL_HEADERS ${ASL_SOURCE_DIR}/
  asl.h asl_pfg.h asl_pfgh.h avltree.h errchk.h funcadd.h getstub.h 
  jac2dim.h jacpdim.h nlp.h nlp2.h obj_adj.h psinfo.h)
add_prefix(ASL_SOURCES ${ASL_SOURCE_DIR}/
  asldate.c       atof.c          auxinfo.c       avltree.c
  b_search.c      basename.c      bscanf.c        com2eval.c
  comeval.c       con1ival.c      con2ival.c      con2val.c
  conadj.c        conpval.c       conscale.c      conval.c
  degree.c 
  derprop.c       dtoa1.c         duthes.c        dynlink.c
  f_read.c        fg_read.c       fg_write.c      fgh_read.c
  fpecatch.c      fpinit.c        fullhes.c       func_add.c
  funcadd1.c      g_fmt.c         genrowno.c      getenv.c
  getstub.c       htcl.c          indic_cons.c    jac0dim.c
  jac2dim.c       jacdim.c        jacinc.c        jacinc1.c
  libnamsave.c    mach.c          mainexit.c      mip_pri.c
  misc.c          mpec_adj.c      mqpcheckv.c     mypow.c
  names.c         nl_obj.c        nqpcheck.c      obj2val.c
  obj_adj.c       obj_prec.c      objconst.c      objval.c
  objval_.c       op_type.c       pfg_read.c      pfghread.c
  printf.c        pshvprod.c      punknown.c      qp_read.c
  qpcheck.c       qsortv.c        readsol.c       repwhere.c
  rops.c          rops2.c         sigcatch.c      sos_add.c
  sphes.c         sscanf.c        stderr.c        studchk0.c
  suf_sos.c       value.c         writesol.c      wrtsol_.c
  ws_desc.c       wsu_desc.c      x2check.c       xectim.c
  xp1known.c      xp2known.c)

# ASL 2 sources
set(ASL2_SOURCE_DIR ${SRCDIR}/solvers2)
set(ASL2_HEADERS ${GENERATED_INCLUDE_DIR}/arith.h)
add_prefix(ASL2_HEADERS ${ASL2_SOURCE_DIR}/
  asl.h asl_pfg.h asl_pfgh.h avltree.h errchk.h funcadd.h getstub.h 
  jac2dim.h jacpdim.h nlp.h nlp2.h  obj_adj.h opno2.h psinfo.h)
set(ASL2_SOURCES ${GENERATED_INCLUDE_DIR}/details.c)
add_prefix(ASL2_SOURCES  ${ASL2_SOURCE_DIR}/
  asldate.c       atof.c          auxinfo.c       avltree.c
  b_search.c      basename.c      bscanf.c        conscale.c
  degree.c        derprop.c       dtoa1.c         duthes.c
  dynlink.c       eval1.c         eval2.c         ewalloc1.c
  ewalloc2.c      f_read.c        fg_read.c       fg_write.c
  fpecatch.c      fpinit.c        fullhes.c       func_add.c
  funcadd1.c      g_fmt.c         genrowno.c      getenv.c
  getstub.c       htcl.c          indic_cons.c    jac0dim.c
  jac2dim.c       jacdim.c        jacinc.c        jacinc1.c
  libnamsave.c    mach.c          mainexit.c      mip_pri.c
  misc.c          mpec_adj.c      mqpcheckv.c     mypow.c
  names.c         nl_obj.c        nqpcheck.c      nqpcheckZ.c
  obj_adj.c       obj_prec.c      objconst.c      objval_.c
  op_type.c       pfghread.c      printf.c        pshvprod.c
  punknown.c      qpcheck.c       qpcheckZ.c      qsortv.c
  readsol.c       repwhere.c      sigcatch.c      sos_add.c
  sphes.c         sscanf.c        stderr.c        studchk0.c
  suf_sos.c       value.c         writesol.c      wrtsol_.c
  ws_desc.c       wsu_desc.c      xectim.c        xp2known.c)

# Set the properties in the parent project, if present
if(hasParent)
  set(ASL1_SOURCE_DIR ${ASL_SOURCE_DIR} PARENT_SCOPE)
  set(ASL1_INCLUDE_DIRS ${GENERATED_INCLUDE_DIR}
    ${ASL_SOURCE_DIR} PARENT_SCOPE)
  set(ASL1_SOURCES ${ASL_SOURCES} PARENT_SCOPE)
    
  set(ASL2_SOURCE_DIR ${ASL_SOURCE_DIR} PARENT_SCOPE)
  set(ASL2_INCLUDE_DIRS ${GENERATED_INCLUDE_DIR}
    ${ASL2_SOURCE_DIR} PARENT_SCOPE)  
  set(ASL2_SOURCES ${ASL2_SOURCES} PARENT_SCOPE)
endif()

# Set preprocessor flags
include(CheckSymbolExists)
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  # Changing the floating point precision is not supported on x64.
  set(ASL_COMPILE_DEFINITIONS No_Control87)
endif()
check_symbol_exists(strtoull stdlib.h HAVE_STRTOULL)
if(NOT HAVE_STRTOULL)
  check_symbol_exists(_strtoui64 stdlib.h HAVE_STRTOUI64)
  if(HAVE_STRTOUI64)
    set(ASL_COMPILE_DEFINITIONS ${ASL_COMPILE_DEFINITIONS} strtoull=_strtoui64)
  endif()
endif()
if(NOT WIN32)
  check_symbol_exists(getrusage "sys/resource.h" HAVE_RUSAGE)
  if(NOT HAVE_RUSAGE)
    set(ASL_COMPILE_DEFINITIONS ${ASL_COMPILE_DEFINITIONS} NO_RUSAGE=1)
  endif()
endif()

macro(createSingleASL name sourcedir sources)
  cmake_parse_arguments(PAR "" "" "DEFINITIONS;INCLUDE;OPTIONS" ${ARGN})
  if (BUILD_SHARED_LIBS)
    add_library(${name} SHARED ${${sources}} ${GENERATED_INCLUDE_DIR}/arith.h)
    message("-- Building shared ${name} library")
    set_target_properties(${name} PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})
  else()
    add_library(${name} STATIC ${${sources}} ${GENERATED_INCLUDE_DIR}/arith.h)
    message("-- Building static ${name} library")
  endif()

  
  if(PAR_DEFINITIONS)
    set(MY_COMPILE_DEFINITIONS ${PAR_DEFINITIONS} ${ASL_COMPILE_DEFINITIONS})
  else()
    set(MY_COMPILE_DEFINITIONS ${ASL_COMPILE_DEFINITIONS})
  endif()
 
  target_include_directories(${name} PUBLIC $<BUILD_INTERFACE:${GENERATED_INCLUDE_DIR}>
                                            $<BUILD_INTERFACE:${sourcedir}>
                                            $<INSTALL_INTERFACE:include>)
  if(PAR_INCLUDE)
    target_include_directories(${name} PUBLIC ${PAR_INCLUDE})
  endif()
  target_compile_definitions(${name} PRIVATE ${MY_COMPILE_DEFINITIONS})
  
  add_dependencies(${name} arith-h)
  target_link_libraries(${name} PUBLIC ${CMAKE_DL_LIBS})
  if(NOT WIN32)
    target_link_libraries(${name} PRIVATE m)
  endif()
  if(MSVC)
    if(${name} MATCHES "dynrt")
	set_property(TARGET ${name} PROPERTY
            MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
    else()
	set_property(TARGET ${name} PROPERTY
	    MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
    target_compile_options(${name} PRIVATE 
      /wd4013 /wd4018 /wd4101 /wd4244 /wd4273 /wd4267 /wd4996)
  else()
    target_compile_options(${name} PRIVATE -Wno-unused-result -Wno-parentheses)
  endif()
  set_property(TARGET ${name} PROPERTY POSITION_INDEPENDENT_CODE ON)

   if(PAR_OPTIONS)
    target_compile_options(${name} PUBLIC ${PAR_OPTIONS})
  endif()
  addToASLFolder("" ${name})
endmacro()

set(ASL_SOURCES ${ASL_SOURCES} ${ASL_HEADERS})
set(ASL2_SOURCES ${ASL2_SOURCES} ${ASL2_HEADERS})

# Create ASL 1 
createSingleASL(asl ${ASL_SOURCE_DIR} ASL_SOURCES)
# Create ASL 2
createSingleASL(asl2 ${ASL2_SOURCE_DIR} ASL2_SOURCES)

if(BUILD_LICCHECK_PRINT)
  createSingleASL(asl-lc ${ASL_SOURCE_DIR} ASL_SOURCES 
                        DEFINITIONS LICCHECKVERSION=${LICCHECKVERSION})
  createSingleASL(asl2-lc ${ASL2_SOURCE_DIR} ASL2_SOURCES
                        DEFINITIONS LICCHECKVERSION=${LICCHECKVERSION})
endif()

# Create libs with support for multiple threads
if(BUILD_MT_LIBS)
  set(aslnames asl-mt asl2-mt)
  createSingleASL(asl-mt ${ASL_SOURCE_DIR} ASL_SOURCES
                              DEFINITIONS ALLOW_OPENMP)
  createSingleASL(asl2-mt ${ASL2_SOURCE_DIR} ASL2_SOURCES
                              DEFINITIONS ALLOW_OPENMP)

  if(BUILD_LICCHECK_PRINT)
    createSingleASL(asl-mt-lc ${ASL_SOURCE_DIR} ASL_SOURCES
                            DEFINITIONS LICCHECKVERSION=${LICCHECKVERSION} ALLOW_OPENMP)
    createSingleASL(asl2-mt-lc ${ASL2_SOURCE_DIR} ASL2_SOURCES
                            DEFINITIONS LICCHECKVERSION=${LICCHECKVERSION} ALLOW_OPENMP)
    set(aslnames ${aslnames} asl-mt-lc asl2-mt-lc)
  endif()

  include(FindOpenMP)
  foreach(name ${aslnames})
    if(APPLE) 
      target_include_directories(${name} PUBLIC ${OpenMP_C_INCLUDE_DIRS})
    else()
      target_compile_options(${name}  PUBLIC ${OpenMP_C_FLAGS}) 
    endif()
  endforeach()
endif()


if(MSVC AND BUILD_DYNRT_LIBS)
  # Create ASLs linked to dynamic runtime on windows
  # createSingleASL adds the appropriate public flags if the library name
  # matches "dynrt"
  createSingleASL(asl-dynrt ${ASL_SOURCE_DIR} ASL_SOURCES)
  createSingleASL(asl2-dynrt ${ASL2_SOURCE_DIR} ASL2_SOURCES)
  if(BUILD_LICCHECK_PRINT)
    createSingleASL(asl-dynrt-lc ${ASL_SOURCE_DIR} ASL_SOURCES
                                DEFINITIONS LICCHECKVERSION=${LICCHECKVERSION})
    createSingleASL(asl2-dynrt-lc ${ASL2_SOURCE_DIR} ASL2_SOURCES
                                DEFINITIONS LICCHECKVERSION=${LICCHECKVERSION})
  endif()
endif()
if(BUILD_MCMODELLARGE)
 createSingleASL(asl2-large ${ASL2_SOURCE_DIR} ASL2_SOURCES)
 if(${ARCH} EQUAL 64)
 target_compile_options(asl2-large PUBLIC "-mcmodel=large")
 endif()
 if(BUILD_LICCHECK_PRINT)
 createSingleASL(asl2-large-lc ${ASL2_SOURCE_DIR} ASL2_SOURCES
                               DEFINITIONS LICCHECKVERSION=${LICCHECKVERSION})
 if(${ARCH} EQUAL 64)
  target_compile_options(asl2-large-lc PUBLIC "-mcmodel=large")
 endif()
 endif()
endif()

# Cpp wrapper, only if explicitly enabled
if(BUILD_CPP)
    set(ASLCPP_SRCDIR ${SRCDIR}/cpp)
    set(ASL_CPP_HEADERS ${ASLCPP_SRCDIR}/aslinterface.h)
    set(ASL_CPP_SOURCES ${ASLCPP_SRCDIR}/aslinterface.cc)
    add_library(aslcpp STATIC ${ASL_CPP_SOURCES} ${ASL_CPP_HEADERS} ${ASL_HEADERS})
    target_link_libraries(aslcpp PUBLIC asl ${CMAKE_DL_LIBS})
    if(NOT WIN32)
      target_link_libraries(aslcpp PRIVATE m)
    endif()
    addToASLFolder("" aslcpp)
endif()

if(BUILD_F2C)
  add_subdirectory(${SRCDIR}/f2c)
  addToASLFolder(support f2c)
endif()
 if(BUILD_ASL_EXAMPLES)
  add_subdirectory(${SRCDIR}/examples)
endif()

if(NOT ASL_SKIP_INSTALL)
  install(FILES ${ASL_HEADERS} ${ASL_SOURCE_DIR}/opcode.hd ${ASL_SOURCE_DIR}/r_opn.hd
                ${GENERATED_INCLUDE_DIR}/stdio1.h ${GENERATED_INCLUDE_DIR}/arith.h
                DESTINATION include/asl  COMPONENT asl) 
  install(FILES ${ASL2_HEADERS} ${ASL2_SOURCE_DIR}/opcode.hd ${ASL2_SOURCE_DIR}/r_opn.hd
                ${GENERATED_INCLUDE_DIR}/stdio1.h ${GENERATED_INCLUDE_DIR}/arith.h
                DESTINATION include/asl2  COMPONENT asl)
  install(TARGETS asl asl2 EXPORT ampl-asl-config DESTINATION lib COMPONENT asl)
  if(BUILD_MT_LIBS)
    install(TARGETS asl-mt asl2-mt EXPORT ampl-asl-config DESTINATION lib COMPONENT asl)
  endif()
  if(MSVC AND BUILD_DYNRT_LIBS)
    install(TARGETS asl-dynrt asl2-dynrt EXPORT ampl-asl-config DESTINATION lib COMPONENT asl)
  endif()
  if(TARGET aslcpp)
      install(FILES ${ASL_CPP_HEADERS} DESTINATION include/aslcpp  COMPONENT asl) 
      install(TARGETS aslcpp EXPORT ampl-asl-config DESTINATION lib  COMPONENT asl)
  endif()
  install(EXPORT ampl-asl-config DESTINATION share/ampl-asl COMPONENT asl)
endif()
