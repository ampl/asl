# Numerical regression tests: every reader of both libraries must give
# derivatives consistent with finite differences and values that agree
# with the other readers, on small .nl fixtures (see README.md).

set(ASL_TEST_NL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/nl)
set(ASL_TEST_OUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/out)
file(MAKE_DIRECTORY ${ASL_TEST_OUT_DIR})

# Checker executables and the readers each library provides.
set(ASL_TEST_CONFIGS)
function(aslTestChecker name lib source)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE ${lib})
  target_compile_definitions(${name} PRIVATE ${ARGN})
  if(UNIX)
    target_link_libraries(${name} PRIVATE m dl)
  endif()
  addToASLFolder(tests ${name})
endfunction()

if(BUILD_ASL_LIB)
  aslTestChecker(derivcheck-asl asl derivcheck.c)
  aslTestChecker(hvseq-asl asl hvseq.c)
  # pfg_read is left out: the solvers library does not install evaluation
  # routines for ASL_read_pfg ("objval called before fg_read").
  list(APPEND ASL_TEST_CONFIGS asl:fg asl:fgh asl:pfgh)
endif()
if(BUILD_ASL2_LIB)
  aslTestChecker(derivcheck-asl2 asl2 derivcheck.c DERIVCHECK_ASL2)
  aslTestChecker(hvseq-asl2 asl2 hvseq.c)
  list(APPEND ASL_TEST_CONFIGS asl2:fg asl2:pfgh)
endif()

# One test per fixture, reader and point.  The first configuration writes
# the reference values the others are compared with.  A fixture with a
# "<fixture>.err" file is also evaluated exactly at its start point
# ("0e"), where the routines listed there must report an error.
function(aslDerivTests fixture)
  set(points 0v 1v 1g)
  if(EXISTS ${ASL_TEST_NL_DIR}/${fixture}.err)
    list(APPEND points 0e)
  endif()
  foreach(point ${points})
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

# One Hessian-vector product at each point of a sequence (pfgh_read): a
# product computed with the Hessian of an earlier point is a mismatch.
# Modes: see hvseq.c.  hvcomp, and modes m and h, for both libraries; the
# other modes for the solvers2 library only: in the solvers library they
# still give wrong products with defined variables.  Modes q Q r R y Y
# (first calls, no pass that finds the points with evaluation errors) are
# not run for fixtures named skip*, which have such points.
function(aslHvSeqTests fixture)
  foreach(lib asl asl2)
    if(TARGET hvseq-${lib})
      add_test(NAME hvseq.${fixture}.${lib}-pfgh COMMAND hvseq-${lib}
        ${ASL_TEST_NL_DIR}/${fixture})
      if(lib STREQUAL "asl2")
        set(modes z d s D S m o h F G)
        if(NOT fixture MATCHES "^skip")
          list(APPEND modes q Q r R y Y)
        endif()
      else()
        set(modes m h)
      endif()
      foreach(mode ${modes})
        add_test(NAME hvseq.${fixture}.${lib}-pfgh-${mode} COMMAND hvseq-${lib}
          ${ASL_TEST_NL_DIR}/${fixture} 12 ${mode})
      endforeach()
    endif()
  endforeach()
endfunction()

# hvseq with some constraint weights 0 (zeromask: bit i for constraint i),
# solvers2: the defined variables of functions with weight 0 need not have
# values and must not be visited.
function(aslHvSeqZeroTest fixture mode zeromask)
  if(TARGET hvseq-asl2)
    add_test(NAME hvseq.${fixture}.asl2-pfgh-${mode}-zero${zeromask}
      COMMAND hvseq-asl2 ${ASL_TEST_NL_DIR}/${fixture} 12 ${mode} ${zeromask})
  endif()
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
  elseif(fixture MATCHES "^skipall")
    # Every hvseq point has an evaluation error: hvseq must fail (exit 4).
    if(TARGET hvseq-asl2)
      add_test(NAME hvseq.${fixture}.asl2-pfgh COMMAND hvseq-asl2
        ${ASL_TEST_NL_DIR}/${fixture})
      set_tests_properties(hvseq.${fixture}.asl2-pfgh PROPERTIES
        PASS_REGULAR_EXPRESSION "too few points compared: 0")
    endif()
  else()
    aslDerivTests(${fixture})
    aslHvSeqTests(${fixture})
  endif()
endforeach()
aslHvSeqZeroTest(dvfirst1 y 7)
aslHvSeqZeroTest(dvfirst1 Y 7)
aslHvSeqZeroTest(dvfirst2 y 2)
aslHvSeqZeroTest(zw2 Y 7)
