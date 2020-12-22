set(CHECK_MPI 1 CACHE BOOL "Check MPI availability")
set(HAVE_AVX2 0 CACHE BOOL "Have CPU with AVX2")
set(HAVE_SSE4_1 0 CACHE BOOL "Have CPU with SSE4.1")
set(HAVE_SSE2 0 CACHE BOOL "Have CPU with SSE2")
set(HAVE_POWER9 0 CACHE BOOL "Have POWER9 CPU")
set(HAVE_POWER8 0 CACHE BOOL "Have POWER8 CPU")
set(HAVE_ARM8 0 CACHE BOOL "Have ARMv8 CPU")
set(NATIVE_ARCH 1 CACHE BOOL "Assume native architecture for SIMD. Use one of the HAVE_* options or set CMAKE_CXX_FLAGS to the appropriate flags if you disable this.")
set(EXE_SUFFIX "" CACHE STRING "Optional executable suffix, can be used to to indicate SIMD level corresponding to manually specified build flags.")

# see https://wiki.debian.org/ArchitectureSpecificsMemo for char signedness
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fsigned-char")

# SIMD instruction sets support
set(ARCH_FLAGS "")
if (HAVE_AVX2)
    if (CMAKE_COMPILER_IS_CLANG)
        set(ARCH_FLAGS "${ARCH_FLAGS} -mavx2")
    else ()
        set(ARCH_FLAGS "${ARCH_FLAGS} -mavx2 -Wa,-q")
    endif ()
    set(X64 1)
elseif (HAVE_SSE4_1)
    set(ARCH_FLAGS "${ARCH_FLAGS} -msse4.1")
    set(X64 1)
elseif (HAVE_SSE2)
    set(ARCH_FLAGS "${ARCH_FLAGS} -msse2")
    set(X64 1)
elseif (HAVE_POWER9)
    set(ARCH_FLAGS "${ARCH_FLAGS} -mcpu=power9 -mvsx")
    set(PPC64 1)
elseif (HAVE_POWER8)
    set(ARCH_FLAGS "${ARCH_FLAGS} -mcpu=power8 -mvsx")
    set(PPC64 1)
elseif (HAVE_ARM8)
    set(ARCH_FLAGS "${ARCH_FLAGS} -march=armv8-a+simd")
    set(ARM 1)
endif ()

if (NATIVE_ARCH AND (ARCH_FLAGS STREQUAL ""))
    if (CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm.*|ARM.*|aarch64.*|AARCH64.*)")
        set(ARM 1)
    elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "PPC64*|ppc64*|powerpc64*")
        set(PPC64 1)
    elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "amd64|AMD64")
        set(X64 1)
    elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "x86|X86")
        set(X86 1)
    else ()
        message(WARNING "CPU without native SIMD instructions. Performance will be bad.")
    endif ()
    if (EMSCRIPTEN)
        set(ARCH_FLAGS "-msimd128 -s WASM=1 -s ASSERTIONS=1")
    elseif (X86 OR X64)
        include(CheckSSEFeatures)
        if (NOT HAVE_SSE4_1_EXTENSIONS)
            if (NOT HAVE_SSE2_EXTENSIONS)
                message(FATAL_ERROR "At least SSE2 is needed to compile")
            endif ()
            message(WARNING "At least SSE4.1 is needed for best performance")
        endif ()
        # clang has a problem with march=native on travis
        if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.0.0")
            set(ARCH_FLAGS "${SSE_FLAGS}")
        else()
            set(ARCH_FLAGS "-march=native")
        endif()
    else ()
        if (PPC64)
            set(ARCH_FLAGS "-mcpu=native")
        else ()
            set(ARCH_FLAGS "-march=native")
        endif ()
    endif ()
endif ()

if (NOT (ARCH_FLAGS STREQUAL ""))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ARCH_FLAGS}")
endif ()

if (CMAKE_COMPILER_IS_ICC)
    # default -fp-model results in inconsistent results in profile search
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model precise")
endif ()

if (X86)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse")
endif ()

find_package(OpenMP)
if (OPENMP_FOUND)
    message("-- Found OpenMP")
    add_definitions(-DOPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing")

# pass some of the CMake settings to the source code
configure_file("hhsuite_config.h.in" "hhsuite_config.h")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(cs)
include_directories(cs)

set(HH_SOURCE
        hash.h
        hhblits.h
        hhblits.cpp
        hhdecl.h
        hhdecl.cpp
        hhhit.h
        hhhit-inl.h
        hhhit.cpp
        hhmatrices.h
        hhmatrices.cpp
        hhsearch.h
        hhsearch.cpp
        hhalign.h
        hhalign.cpp
        hhhitlist.h
        hhhitlist-inl.h
        hhhitlist.cpp
        hhposteriordecoder.h
        hhposteriordecoder.cpp
        hhutil.h
        hhutil-inl.h
        hhutil.cpp
        util.h
        util-inl.h
        util.cpp
        hhalignment.h
        hhalignment.cpp
        hhforwardalgorithm.cpp
        hhhmm.h
        hhhmm.cpp
        hhposteriordecoderrunner.h
        hhposteriordecoderrunner.cpp
        hhviterbialgorithm.cpp
        hhfullalignment.h
        hhfullalignment.cpp
        hhhmmsimd.h
        hhhmmsimd.cpp
        hhposteriormatrix.h
        hhposteriormatrix.cpp
        hhviterbi.h
        hhviterbi.cpp
        hhbacktracemac.cpp
        hhmacalgorithm.cpp
        hhprefilter.h
        hhprefilter.cpp
        hhviterbimatrix.h
        hhviterbimatrix-inl.h
        hhviterbimatrix.cpp
        hhbackwardalgorithm.cpp
        ffindexdatabase.h
        ffindexdatabase.cpp
        hhdatabase.h
        hhdatabase.cpp
        hhhalfalignment.h
        hhhalfalignment.cpp
        hhviterbirunner.h
        hhviterbirunner.cpp
        hhfunc.h
        hhfunc.cpp
        list.h
        log.h
        )

add_library(hhviterbialgorithm_with_celloff hhviterbialgorithm.cpp)
set_property(TARGET hhviterbialgorithm_with_celloff PROPERTY COMPILE_FLAGS "-DVITERBI_CELLOFF=1")

add_library(hhviterbialgorithm_with_celloff_and_ss hhviterbialgorithm.cpp)
set_property(TARGET hhviterbialgorithm_with_celloff_and_ss PROPERTY COMPILE_FLAGS "-DVITERBI_CELLOFF=1 -DVITERBI_SS_SCORE=1")

add_library(hhviterbialgorithm_and_ss hhviterbialgorithm.cpp)
set_property(TARGET hhviterbialgorithm_and_ss PROPERTY COMPILE_FLAGS "-DVITERBI_SS_SCORE=1")

add_library(HH_OBJECTS ${HH_SOURCE})
add_dependencies(HH_OBJECTS generated)
target_link_libraries(HH_OBJECTS
        ffindex
        CS_OBJECTS
        hhviterbialgorithm_with_celloff
        hhviterbialgorithm_and_ss
        hhviterbialgorithm_with_celloff_and_ss)
if (OPENMP_FOUND AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(HH_OBJECTS ${OpenMP_CXX_LIBRARIES})
endif()

add_executable(hhblits${EXE_SUFFIX} hhblits_app.cpp)
target_link_libraries(hhblits${EXE_SUFFIX} HH_OBJECTS)

add_executable(hhmake${EXE_SUFFIX} hhmake.cpp)
target_link_libraries(hhmake${EXE_SUFFIX} HH_OBJECTS)

add_executable(hhfilter${EXE_SUFFIX} hhfilter.cpp)
target_link_libraries(hhfilter${EXE_SUFFIX} HH_OBJECTS)

add_executable(hhsearch${EXE_SUFFIX} hhblits_app.cpp)
target_link_libraries(hhsearch${EXE_SUFFIX} HH_OBJECTS)
set_property(TARGET hhsearch${EXE_SUFFIX} PROPERTY COMPILE_FLAGS "-DHHSEARCH=1")

add_executable(hhalign${EXE_SUFFIX} hhblits_app.cpp)
target_link_libraries(hhalign${EXE_SUFFIX} HH_OBJECTS)
set_property(TARGET hhalign${EXE_SUFFIX} PROPERTY COMPILE_FLAGS "-DHHALIGN=1")

add_executable(hhconsensus${EXE_SUFFIX} hhconsensus.cpp)
target_link_libraries(hhconsensus${EXE_SUFFIX} HH_OBJECTS)

add_library(A3M_COMPRESS a3m_compress.cpp)
if (OPENMP_FOUND AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(A3M_COMPRESS ${OpenMP_CXX_LIBRARIES})
endif()

add_executable(a3m_extract${EXE_SUFFIX} a3m_extract.cpp)
target_link_libraries(a3m_extract${EXE_SUFFIX} ffindex A3M_COMPRESS)

add_executable(a3m_reduce${EXE_SUFFIX} a3m_reduce.cpp)
target_link_libraries(a3m_reduce${EXE_SUFFIX} ffindex A3M_COMPRESS)

add_executable(a3m_database_reduce${EXE_SUFFIX} a3m_database_reduce.cpp )
target_link_libraries(a3m_database_reduce${EXE_SUFFIX} ffindex A3M_COMPRESS)

add_executable(a3m_database_extract${EXE_SUFFIX} a3m_database_extract.cpp)
target_link_libraries(a3m_database_extract${EXE_SUFFIX} ffindex A3M_COMPRESS)

add_executable(a3m_database_filter${EXE_SUFFIX} a3m_database_filter.cpp)
target_link_libraries(a3m_database_filter${EXE_SUFFIX} ffindex A3M_COMPRESS)

add_executable(cstranslate${EXE_SUFFIX} cs/cstranslate_app.cc)
target_link_libraries(cstranslate${EXE_SUFFIX} HH_OBJECTS A3M_COMPRESS)

INSTALL(TARGETS
        hhblits${EXE_SUFFIX}
        hhmake${EXE_SUFFIX}
        hhfilter${EXE_SUFFIX}
        hhsearch${EXE_SUFFIX}
        hhalign${EXE_SUFFIX}
        hhconsensus${EXE_SUFFIX}
        a3m_extract${EXE_SUFFIX}
        a3m_reduce${EXE_SUFFIX}
        a3m_database_reduce${EXE_SUFFIX}
        a3m_database_extract${EXE_SUFFIX}
        a3m_database_filter${EXE_SUFFIX}
        cstranslate${EXE_SUFFIX}
        DESTINATION bin
        )

if (OPENMP_FOUND)
    add_executable(hhblits_omp${EXE_SUFFIX} hhblits_omp.cpp)
    target_link_libraries(hhblits_omp${EXE_SUFFIX} HH_OBJECTS)

    add_executable(hhsearch_omp${EXE_SUFFIX} hhblits_omp.cpp)
    target_link_libraries(hhsearch_omp${EXE_SUFFIX} HH_OBJECTS)
    set_property(TARGET hhsearch_omp${EXE_SUFFIX} PROPERTY COMPILE_FLAGS "-DHHSEARCH=1")

    add_executable(hhalign_omp${EXE_SUFFIX} hhblits_omp.cpp)
    target_link_libraries(hhalign_omp${EXE_SUFFIX} HH_OBJECTS)
    set_property(TARGET hhalign_omp${EXE_SUFFIX} PROPERTY COMPILE_FLAGS "-DHHALIGN=1")

    add_executable(hhblits_ca3m${EXE_SUFFIX} hhblits_ca3m.cpp)
    target_link_libraries (hhblits_ca3m${EXE_SUFFIX} HH_OBJECTS)

    INSTALL(TARGETS hhblits_omp${EXE_SUFFIX} hhsearch_omp${EXE_SUFFIX} hhalign_omp${EXE_SUFFIX} hhblits_ca3m${EXE_SUFFIX} DESTINATION bin)
endif ()

if (CHECK_MPI)
find_package(MPI QUIET)
if (MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})

    add_executable(hhblits_mpi${EXE_SUFFIX} hhblits_mpi.cpp)
    target_link_libraries(hhblits_mpi${EXE_SUFFIX} HH_OBJECTS mpq ${MPI_CXX_LIBRARIES})
    set_target_properties(hhblits_mpi${EXE_SUFFIX} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    set_target_properties(hhblits_mpi${EXE_SUFFIX} PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    add_executable(hhsearch_mpi${EXE_SUFFIX} hhblits_mpi.cpp)
    target_link_libraries(hhsearch_mpi${EXE_SUFFIX} HH_OBJECTS mpq ${MPI_CXX_LIBRARIES})
    set_target_properties(hhsearch_mpi${EXE_SUFFIX} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} -DHHSEARCH=1")
    set_target_properties(hhsearch_mpi${EXE_SUFFIX} PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    add_executable(hhalign_mpi${EXE_SUFFIX} hhblits_mpi.cpp)
    target_link_libraries(hhalign_mpi${EXE_SUFFIX} HH_OBJECTS mpq ${MPI_CXX_LIBRARIES})
    set_target_properties(hhalign_mpi${EXE_SUFFIX} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} -DHHALIGN=1")
    set_target_properties(hhalign_mpi${EXE_SUFFIX} PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    add_executable(cstranslate_mpi${EXE_SUFFIX} cs/cstranslate_mpi_app.cc)
    target_link_libraries(cstranslate_mpi${EXE_SUFFIX} HH_OBJECTS A3M_COMPRESS mpq ${MPI_CXX_LIBRARIES})
    set_target_properties(cstranslate_mpi${EXE_SUFFIX} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    set_target_properties(cstranslate_mpi${EXE_SUFFIX} PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    install(TARGETS hhblits_mpi${EXE_SUFFIX} hhsearch_mpi${EXE_SUFFIX} hhalign_mpi${EXE_SUFFIX} cstranslate_mpi${EXE_SUFFIX} DESTINATION bin)
endif ()
endif ()

