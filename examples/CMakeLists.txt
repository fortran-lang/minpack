include_directories(${PROJECT_BINARY_DIR}/src)

add_executable(example_lmder1 example_lmder1.f90)
target_link_libraries(example_lmder1 minpack)
