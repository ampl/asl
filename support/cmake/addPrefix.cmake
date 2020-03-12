# Adds a prefix to arguments.
function (add_prefix var prefix)
  set(result ${${var}})
  foreach (arg ${ARGN})
    set(result ${result} "${prefix}${arg}")
  endforeach ()
  set(${var} ${result} PARENT_SCOPE)
endfunction ()


# Add files specified as additional arguments to the specified folder.
function (add_to_folder folder)
  foreach (target ${ARGN})
	if(TARGET ${target})
		set_property(TARGET ${target} PROPERTY FOLDER ${folder})
	endif()
  endforeach ()
endfunction ()

# Add the specified prefix to each argument, add the resulting list to a source
# group and append it to ${var}.
# The optional flag [HEADER_ONLY] can be specified to
# declare the whole source group as inert
function (add_to_source_group var group prefix)
  set(options HEADER_ONLY)
	cmake_parse_arguments(add_to_source_group "${options}" "" "" ${ARGN} )
  set(files )
  foreach (f ${ARGN})
  if(${f} MATCHES "HEADER_ONLY")
    continue()
  endif()
  set(files ${files} ${prefix}${f})
  endforeach ()
  set(${var} ${${var}} ${files} PARENT_SCOPE)
  source_group(${group} FILES ${files})
  if(add_to_source_group_HEADERONLY)
  set_source_files_properties(${group} HEADER_FILE_ONLY TRUE)
  endif()
endfunction ()