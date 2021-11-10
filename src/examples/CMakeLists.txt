if(MSVC)
    add_compile_options(/MT$<$<CONFIG:Debug>:d>)
endif()

function(aslExample name)
  add_executable(${name} ${name}.c ${ARGN})
  target_link_libraries(${name} PRIVATE asl)
  if(UNIX)
    target_link_libraries(${name} PRIVATE m dl)
  endif()
  addToASLFolder(examples ${name})
  install(TARGETS ${name} DESTINATION examples COMPONENT examples)
endfunction()

function(portExample name)
  aslExample(${name} ${ARGN})
  target_link_libraries(${name} PRIVATE cport)
endfunction()

add_library(funcadd0 OBJECT ../solvers/funcadd0.c)
add_dependencies(funcadd0 arith-h)
get_target_property(asl_INCLUDES asl INCLUDE_DIRECTORIES)
target_include_directories(funcadd0 PUBLIC ${asl_INCLUDES})
addToASLFolder(examples funcadd0)

function(funcadd0Example name)
    aslExample(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE funcadd0)
endfunction()

add_subdirectory(PORT/cport)
addToASLFolder(examples cport)

set(ASLEXAMPLES evalchk gjh nlcopy et)
foreach(example ${ASLEXAMPLES})
  aslExample(${example})
endforeach()
if(LINUX)
  target_link_libraries(et PRIVATE rt)
endif()

set(FUNCADDEXAMPLES lin1 lin2 lin3 linrc dualconv qtest)
foreach(example ${FUNCADDEXAMPLES})
  funcadd0Example(${example})
endforeach()


portExample(mng1)
portExample(nl21)
portExample(mng rvmsg.c keywds.c)
portExample(nl2 rvmsg.c keywds.c)
portExample(mnh rvmsg.c keywds.c)
target_compile_definitions(mnh PRIVATE Std_dev)

