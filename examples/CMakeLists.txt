macro(ADDEXAMPLE name)
    add_executable(example_${name} example_${name}.f90)
    target_link_libraries(example_${name} "${PROJECT_NAME}")
    add_test(NAME ${name}
             COMMAND $<TARGET_FILE:example_${name}> ${CMAKE_CURRENT_BINARY_DIR}
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(ADDEXAMPLE)

set(examples
    "hybrd"
    "hybrd1"
    "lmder1"
    "lmdif1"
    "primes"
)

foreach(example IN LISTS examples)
  ADDEXAMPLE(${example})
endforeach()
