cmake_minimum_required(VERSION 3.14.0)

# Include overwrites before setting up the project
set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_CURRENT_SOURCE_DIR}/config/DefaultFlags.cmake)

project(minpack
        LANGUAGES Fortran
	VERSION "2.0.0"
        DESCRIPTION "Minpack includes software for solving nonlinear equations and nonlinear least squares problems"
)

include(CTest)

# Follow GNU conventions for installation directories
include(GNUInstallDirs)

# --- CMake specific configuration and package data export
add_subdirectory(config)

# --- compiler feature checks
include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)

add_subdirectory(src)

if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(test)
  add_subdirectory(examples)
endif()

install(EXPORT ${PROJECT_NAME}-targets
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
