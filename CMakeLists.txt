cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

project("hhsuite")

set (HHSUITE_VERSION_MAJOR 2)
set (HHSUITE_VERSION_MINOR 0)
set (HHSUITE_VERSION_PATCH 16)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
	"src/hhsuite_config.h.in"
	"${PROJECT_BINARY_DIR}/hhsuite_config.h"
)

add_subdirectory(src)
add_subdirectory(data)
add_subdirectory(scripts)
add_subdirectory(doc)
add_subdirectory(lib/ffindex)
