set(TESTS mult_test mult_test_small)
set(ROTESTS polymul_test polymul_test_small)

foreach(X IN LISTS TESTS)
  add_executable("${X}" "${X}.cpp")
  target_link_libraries("${X}" neon)
  target_include_directories("${X}" PUBLIC ../include)
endforeach()

foreach(X IN LISTS ROTESTS)
  add_executable("${X}" "${X}.c")
  target_link_libraries("${X}" rader_outer)
  target_include_directories("${X}" PUBLIC ../include)
endforeach()
