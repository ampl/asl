# This initializes the ARCH environment variable to 32 or 64,
# that is used by other scripts (notable setArchitecture) 
# to appropriately set compiler options, definitions and sources
# In case the build is Windows/MinGW, it sets the variable
# MINGW to 1.

include(setArchitecture)
targetArchitecture(TARGETARCH)

if(MSVC)
  # Bittage in VS depends on the generator used, so it does not need to be specified
  if(${TARGETARCH} MATCHES "i386")
    set(ARCH 32)
  else()
    set(ARCH 64)
  endif()
endif()
if(NOT ARCH) # not set on command line, get the default for the current compiler
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(ARCH 64)
  else ()
    set(ARCH 32)
  endif ()
  message("Defaulting to ${ARCH} bits build. To build 32 bits binaries on a 64 bits compiler, define the architecture of the build, calling CMake with the parameter -DARCH=32")
else()
  message("Building for ${ARCH} bits")
endif()

if(WIN32)
  if(${CMAKE_GENERATOR} MATCHES "MinGW")
    set(MINGW 1)
  endif()
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(LINUX 1)
endif()

message(STATUS "System: ${CMAKE_SYSTEM}")
message(STATUS "Build type: ${CMAKE_C_COMPILER_ID}")
message(STATUS "Compiler version: ${CMAKE_C_COMPILER_VERSION}")
message(STATUS "Generator version: ${CMAKE_GENERATOR}")
message(STATUS "Platform toolset: ${CMAKE_VS_PLATFORM_TOOLSET}" )
message(STATUS "Architecture: ${TARGETARCH}")