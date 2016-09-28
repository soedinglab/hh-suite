cmake_minimum_required(VERSION 2.8.12)
project(hhsuite CXX)

SET(HAVE_SSSE3 0 CACHE BOOL "Have SSSE3")
SET(HAVE_AVX2 0 CACHE BOOL "Have AVX2")

#set(CMAKE_BUILD_TYPE Release)

# find compiler
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  message("-- Compiler is clang(++)")
  set(CMAKE_COMPILER_IS_CLANG 1)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  message("-- Compiler is GNU ")
  set(CMAKE_COMPILER_IS_GNUCXX 1)
  set(CMAKE_COMPILER_IS_GNUCCC 1)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  message("-- Compiler is ICC ")
  set(CMAKE_COMPILER_IS_ICC 1)
endif()

set (HHSUITE_VERSION_MAJOR 3)
set (HHSUITE_VERSION_MINOR 0)
set (HHSUITE_VERSION_PATCH 0)
set (HHSUITE_DATE "15-03-2015")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

INSTALL(FILES
  hhsuite-userguide.pdf
  LICENSE
  README
  PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
  DESTINATION ./
)

include_directories(src)
include_directories(lib/ffindex/src)

add_subdirectory(src)
add_subdirectory(data)
add_subdirectory(scripts)
add_subdirectory(test)
#add_subdirectory(doc)
add_subdirectory(lib/ffindex)


SET (CPACK_GENERATOR "DEB;TGZ;ZIP;RPM")
SET (CPACK_DEBIAN_PACKAGE_MAINTAINER "Markus Meier")

set(CPACK_SOURCE_IGNORE_FILES
  ".synctex.gz$"
  ".toc$"
  ".blg$"
  ".bbl$"
  ".aux$"
  ".log$"
  ".dvi$"
  ".out.ps$"
  "/rpm/"
  "^${PROJECT_SOURCE_DIR}/build/"
  "^${PROJECT_SOURCE_DIR}/.git/"
  "^${PROJECT_SOURCE_DIR}/.gitignore"
  "^${PROJECT_SOURCE_DIR}/.gitattributes"
  "^${PROJECT_SOURCE_DIR}/.gitmodules"
  "^${PROJECT_SOURCE_DIR}/.pydevproject"
  "^${PROJECT_SOURCE_DIR}/.project"
  "^${PROJECT_SOURCE_DIR}/.cproject"
  "^${PROJECT_SOURCE_DIR}/lib/ffindex/.git"
  "^${PROJECT_SOURCE_DIR}/lib/ffindex/.git"
  "^${PROJECT_SOURCE_DIR}/.autotools"
  "^${PROJECT_SOURCE_DIR}/README.md"
)

set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set (CPACK_PACKAGE_VERSION_MAJOR "${HHSUITE_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${HHSUITE_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${HHSUITE_VERSION_PATCH}")
INCLUDE(CPack)

enable_testing()
