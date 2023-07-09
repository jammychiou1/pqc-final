set(PROFILES ntt_17_profile intt_17_profile ntt_3_2_profile intt_3_2_profile base_mul_profile polymul_profile)

foreach(X IN LISTS PROFILES)
  add_executable("ro_${X}" "${X}.c")
  target_link_libraries("ro_${X}" rader_outer)
  target_include_directories("ro_${X}" PUBLIC ../../include)
endforeach()
