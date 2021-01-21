# LOQO solver library
set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_prefix(LOQOSLVR_SOURCES ${SRCDIR}/
choles.c linalg.c options.c optlp.c iolp.c nonlin.c nlmps.c 
hash.c strdup.c cputime.c prep.c sortsearch.c version.c)
# checks.c  was included in the mkfile but not in the sources
add_prefix(LOQOSLVR_INCLUDES ${SRCDIR}/
loqo.h myalloc.h choles.h hash.h)
set_source_files_properties(${LOQOSLVR_INCLUDES} HEADER_FILE_ONLY)
add_library(libloqo OBJECT ${LOQOSLVR_SOURCES} ${LOQOSLVR_INCLUDES})
target_include_directories(libloqo PUBLIC ${SRCDIR})
setArchitectureFlags(${ARCH} libloqo)
if(MSVC)
  target_compile_definitions(libloqo PRIVATE PC)
  target_compile_options(libloqo PRIVATE /MT$<$<CONFIG:Debug>:d>)
else()
	target_compile_options(libloqo PRIVATE -Wno-implicit-function-declaration)
	if(APPLE)
		target_compile_definitions(libloqo PUBLIC GNU)
	endif()
endif()

