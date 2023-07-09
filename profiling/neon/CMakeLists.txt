set(PROFILES ntt_10_profile ntt_9_profile intt_10_x10_profile intt_9_x9_profile crt_profile center_poly_profile div90_profile base_mul_profile mult_low_profile mult_profile)

foreach(X IN LISTS PROFILES)
  add_executable("${X}" "${X}.cpp")
  target_link_libraries("${X}" neon)
  target_include_directories("${X}" PUBLIC ../../include)
endforeach()
