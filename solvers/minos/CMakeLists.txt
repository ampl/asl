cmake_minimum_required(VERSION 3.3)

# Minos library
set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src/csrc)
add_prefix(MINOS_LIB_SOURCES ${SRCDIR}/
mi02lib.c mi10unix.c mi11sys.c 
mi12solv.c mi15blas.c mi17util.c 
mi20amat.c mi21amat.c mi25bfac.c mi26bfac.c
mi27lu.c mi30spec.c mi35mps.c mi36mps.c 
mi40bfil.c mi50lp.c mi60srch.c mi65rmod.c 
mi70nobj.c mi80ncon.c mi81ncon.c)
add_library(minos551_lib ${MINOS_LIB_SOURCES})
setArchitectureFlags(${ARCH} minos551_lib)
if(MSVC)
  target_compile_options(minos551_lib PRIVATE /MT$<$<CONFIG:Debug>:d>)
else()
  target_link_libraries(minos551_lib m f2c)
endif()
target_include_directories(minos551_lib PUBLIC $<TARGET_PROPERTY:f2c,INCLUDE_DIRECTORIES>) # for f2c.h


