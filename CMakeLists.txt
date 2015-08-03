cmake_minimum_required(VERSION 2.8)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)
SET(HAVE_SSSE3 0 CACHE BOOL "Have SSSE3")
SET(HAVE_AVX2 0 CACHE BOOL "Have AVX2")
project("hhsuite")

#set(CMAKE_BUILD_TYPE Release)

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

set (HHSUITE_VERSION_MAJOR 3)
set (HHSUITE_VERSION_MINOR 0)
set (HHSUITE_VERSION_PATCH 0)
set (HHSUITE_DATE "15-03-2015")

include_directories(src)
include_directories(lib/ffindex/src)

add_subdirectory(src)
add_subdirectory(data)
add_subdirectory(scripts)
#add_subdirectory(doc)
add_subdirectory(lib/ffindex)
