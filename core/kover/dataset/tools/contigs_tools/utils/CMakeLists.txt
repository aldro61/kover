
################################################################################
#  EXAMPLES 
################################################################################

# We add the compilation options for the library
add_definitions (${gatb-core-flags} " -Wno-ambiguous-member-template")

include_directories (${gatb-core-includes})

list (APPEND PROGRAMS multidsk)

FOREACH (program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${gatb-core-libraries} progressbar)
ENDFOREACH (program)

################################################################################
#  INSTALLATION 
################################################################################
INSTALL (TARGETS multidsk DESTINATION bin)
