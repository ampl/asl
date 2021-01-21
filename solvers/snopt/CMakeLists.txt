# SNOPT solver interface,
# target must be defined in its own directory because of the -no-pie below, 
# that does not seem to work reliably on cmake yet

set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR})
set(LIBS ${snopt_LIBRARY} f2c)
if (LINUX)
  set(LIBS ${LIBS} gfortran)
  if(ARCH EQUAL 64)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-no-pie" HasNopie)
    if(HasNopie)
      # Idiomatic way to specify this does not seem to work
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -no-pie")
    endif()
  endif()
endif()
set(snopt_DRIVER_EXTERNALBUILD FALSE)
createSolverTarget(snopt)

if(WIN32 AND ARCH EQUAL 32)
  target_link_options(snopt PRIVATE "/SAFESEH:NO" "/NODEFAULTLIB:\"libc.lib\"") 
endif()