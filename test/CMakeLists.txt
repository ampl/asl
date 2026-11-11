# Numerical regression tests: every reader of both libraries must give
# derivatives consistent with finite differences and values that agree
# with the other readers, on small .nl fixtures (see README.md).

set(ASL_TEST_NL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/nl)
set(ASL_TEST_OUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/out)
file(MAKE_DIRECTORY ${ASL_TEST_OUT_DIR})

# Checker executables and the readers each library provides.
set(ASL_TEST_CONFIGS)
function(aslTestChecker name lib)
  add_executable(${name} derivcheck.c)
  target_link_libraries(${name} PRIVATE ${lib})
  target_compile_definitions(${name} PRIVATE ${ARGN})
  if(UNIX)
    target_link_libraries(${name} PRIVATE m dl)
  endif()
  addToASLFolder(tests ${name})
endfunction()

if(BUILD_ASL_LIB)
  aslTestChecker(derivcheck-asl asl)
  # pfg_read is left out: the solvers library does not install evaluation
  # routines for ASL_read_pfg ("objval called before fg_read").
  list(APPEND ASL_TEST_CONFIGS asl:fg asl:fgh asl:pfgh)
endif()
if(BUILD_ASL2_LIB)
  aslTestChecker(derivcheck-asl2 asl2 DERIVCHECK_ASL2)
  list(APPEND ASL_TEST_CONFIGS asl2:fg asl2:pfgh)
endif()

# One test per fixture, reader and point.  The first configuration writes
# the reference values the others are compared with.
function(aslDerivTests fixture)
  foreach(point 0v 1v 1g)
    string(SUBSTRING ${point} 0 1 seed)
    string(SUBSTRING ${point} 1 1 order)
    set(refname)
    foreach(config ${ASL_TEST_CONFIGS})
      string(REPLACE ":" ";" parts ${config})
      list(GET parts 0 lib)
      list(GET parts 1 reader)
      set(name deriv.${fixture}.${lib}-${reader}.${point})
      set(dumpfile ${ASL_TEST_OUT_DIR}/${fixture}.${lib}-${reader}.${point})
      if(refname)
        add_test(NAME ${name} COMMAND derivcheck-${lib}
          ${ASL_TEST_NL_DIR}/${fixture} ${reader} ${seed} ${order}
          ${dumpfile} ${reffile})
        set_tests_properties(${name} PROPERTIES FIXTURES_REQUIRED ${refname})
      else()
        add_test(NAME ${name} COMMAND derivcheck-${lib}
          ${ASL_TEST_NL_DIR}/${fixture} ${reader} ${seed} ${order} ${dumpfile})
        set(refname ref.${fixture}.${point})
        set(reffile ${dumpfile})
        set_tests_properties(${name} PROPERTIES FIXTURES_SETUP ${refname})
      endif()
    endforeach()
  endforeach()
endfunction()

# A reader must reject (or accept) a fixture when asked only to read it.
function(aslReadTest fixture lib reader expectation)
  set(name read.${fixture}.${lib}-${reader})
  add_test(NAME ${name} COMMAND derivcheck-${lib}
    ${ASL_TEST_NL_DIR}/${fixture} ${reader} 0 r)
  set_tests_properties(${name} PROPERTIES
    PASS_REGULAR_EXPRESSION "${expectation}")
endfunction()

file(GLOB ASL_TEST_FIXTURES RELATIVE ${ASL_TEST_NL_DIR} ${ASL_TEST_NL_DIR}/*.nl)
foreach(file ${ASL_TEST_FIXTURES})
  string(REGEX REPLACE "\\.nl$" "" fixture ${file})
  if(fixture MATCHES "^owner_")
    # Defined variables whose owner field contradicts the header: rejected,
    # except by the solvers pfg/pfgh reader, which does not use the field.
    set(rejected "contrary to the header")
    if(BUILD_ASL_LIB)
      aslReadTest(${fixture} asl fg "${rejected}")
      aslReadTest(${fixture} asl fgh "${rejected}")
      aslReadTest(${fixture} asl pfgh "read n=")
    endif()
    if(BUILD_ASL2_LIB)
      aslReadTest(${fixture} asl2 fg "${rejected}")
      aslReadTest(${fixture} asl2 pfgh "${rejected}")
    endif()
  else()
    aslDerivTests(${fixture})
  endif()
endforeach()
