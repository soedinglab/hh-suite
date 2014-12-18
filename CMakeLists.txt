cmake_minimum_required(VERSION 2.4)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

project("hhsuite")

set(CMAKE_BUILD_TYPE Release)

# find compiler
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  message("-- Compilter is clang(++)")
  set(CMAKE_COMPILER_IS_CLANG 1)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  message("-- Compilter is GNU ")
  set(CMAKE_COMPILER_IS_GNUCXX 1)
  set(CMAKE_COMPILER_IS_GNUCCC 1)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  set(CMAKE_COMPILER_IS_ICC 1)
endif()

set (HHSUITE_VERSION_MAJOR 2)
set (HHSUITE_VERSION_MINOR 0)
set (HHSUITE_VERSION_PATCH 16)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
	"src/hhsuite_config.h.in"
	"${PROJECT_BINARY_DIR}/hhsuite_config.h"
)

include_directories(src)
include_directories(lib/ffindex/src)

add_subdirectory(src)
add_subdirectory(data)
add_subdirectory(scripts)
add_subdirectory(doc)
add_subdirectory(lib/ffindex)
