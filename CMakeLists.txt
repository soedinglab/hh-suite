cmake_minimum_required(VERSION 2.8.12)
project(hhsuite CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(HAVE_SSE2 0 CACHE BOOL "Have SSE2")
set(HAVE_AVX2 0 CACHE BOOL "Have AVX2")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)

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
endif ()

message("-- Processor is ${CMAKE_SYSTEM_PROCESSOR}")
if (CMAKE_SYSTEM_PROCESSOR MATCHES "amd64.*|x86_64.*|AMD64.*")
    set(X86_64 1)
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "i686.*|i386.*|x86.*|amd64.*|AMD64.*")
    set(X86 1)
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm.*|ARM.*)")
    set(ARM 1)
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*)")
    set(AARCH64 1)
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "PPC64*|ppc64*|powerpc64*")
    set(PPC64 1)
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "PPC*|ppc*|powerpc*")
    set(PPC 1)
endif ()

if (NOT (X86 OR X86_64))
    message("!! NOT X86 OR X86_64 !!")
endif ()

if (${HAVE_SANITIZER})
    include(FindUBSan)
    include(FindASan)
    include(FindMSan)
    include(FindTSan)
endif ()

set(HHSUITE_VERSION_MAJOR 3)
set(HHSUITE_VERSION_MINOR 0)
set(HHSUITE_VERSION_PATCH 3)
set(HHSUITE_DATE "14-07-2017")

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

add_subdirectory(lib/ffindex)
add_subdirectory(src)
add_subdirectory(data)
add_subdirectory(scripts)
add_subdirectory(test)
#add_subdirectory(doc)

set(CPACK_GENERATOR "DEB;TGZ;ZIP;RPM")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Markus Meier")

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

set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${HHSUITE_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${HHSUITE_VERSION_MINOR}-beta")
set(CPACK_PACKAGE_VERSION_MINOR "${HHSUITE_VERSION_PATCH}")
include(CPack)
