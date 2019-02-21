include(${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
check_function_exists(fmemopen HAVE_FMEMOPEN)

if (NOT HAVE_FMEMOPEN)
    set(EXT_SRC fmemopen.c fmemopen.h)
    add_library(ext ${EXT_SRC})
endif ()
