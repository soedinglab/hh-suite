set(CHECK_MPI 1 CACHE BOOL "Check MPI availability")
set(HAVE_AVX2 0 CACHE BOOL "Have CPU with AVX2")
set(HAVE_SSE4_1 0 CACHE BOOL "Have CPU with SSE4.1")
set(HAVE_SSE2 0 CACHE BOOL "Have CPU with SSE2")
set(HAVE_POWER9 0 CACHE BOOL "Have POWER9 CPU")
set(HAVE_POWER8 0 CACHE BOOL "Have POWER8 CPU")
set(HAVE_ARM8 0 CACHE BOOL "Have ARMv8 CPU")
set(NATIVE_ARCH 1 CACHE BOOL "Assume native architecture for SIMD. Use one of the HAVE_* options or set CMAKE_CXX_FLAGS to the appropriate flags if you disable this.")

# see https://wiki.debian.org/ArchitectureSpecificsMemo for char signedness
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fsigned-char")
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

# SIMD instruction sets support
set(ARCH_FLAGS "")
if (HAVE_AVX2)
    if (CMAKE_COMPILER_IS_CLANG)
        set(ARCH_FLAGS "${ARCH_FLAGS} -mavx2")
    else ()
        set(ARCH_FLAGS "${ARCH_FLAGS} -mavx2 -Wa,-q")
    endif ()
elseif (HAVE_SSE4_1)
    set(ARCH_FLAGS "${ARCH_FLAGS} -msse4.1")
elseif (HAVE_SSE2)
    set(ARCH_FLAGS "${ARCH_FLAGS} -msse2")
elseif (HAVE_POWER9)
    set(ARCH_FLAGS "${ARCH_FLAGS} -mcpu=power9 -mvsx")
elseif (HAVE_POWER8)
    set(ARCH_FLAGS "${ARCH_FLAGS} -mcpu=power8 -mvsx")
elseif (HAVE_ARM8)
    set(ARCH_FLAGS "${ARCH_FLAGS} -march=armv8-a+simd")
endif ()

if (NATIVE_ARCH AND (ARCH_FLAGS STREQUAL ""))
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

add_executable(hhblits hhblits_app.cpp)
target_link_libraries(hhblits HH_OBJECTS)

add_executable(hhmake hhmake.cpp)
target_link_libraries(hhmake HH_OBJECTS)

add_executable(hhfilter hhfilter.cpp)
target_link_libraries(hhfilter HH_OBJECTS)

add_executable(hhsearch hhblits_app.cpp)
target_link_libraries(hhsearch HH_OBJECTS)
set_property(TARGET hhsearch PROPERTY COMPILE_FLAGS "-DHHSEARCH=1")

add_executable(hhalign hhblits_app.cpp)
target_link_libraries(hhalign HH_OBJECTS)
set_property(TARGET hhalign PROPERTY COMPILE_FLAGS "-DHHALIGN=1")

add_executable(hhconsensus hhconsensus.cpp)
target_link_libraries(hhconsensus HH_OBJECTS)

add_library(A3M_COMPRESS a3m_compress.cpp)
if (OPENMP_FOUND AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(A3M_COMPRESS ${OpenMP_CXX_LIBRARIES})
endif()

add_executable(a3m_extract a3m_extract.cpp)
target_link_libraries(a3m_extract ffindex A3M_COMPRESS)

add_executable(a3m_reduce a3m_reduce.cpp)
target_link_libraries(a3m_reduce ffindex A3M_COMPRESS)

add_executable(a3m_database_reduce a3m_database_reduce.cpp )
target_link_libraries(a3m_database_reduce ffindex A3M_COMPRESS)

add_executable(a3m_database_extract a3m_database_extract.cpp)
target_link_libraries(a3m_database_extract ffindex A3M_COMPRESS)

add_executable(a3m_database_filter a3m_database_filter.cpp)
target_link_libraries(a3m_database_filter ffindex A3M_COMPRESS)

add_executable(cstranslate cs/cstranslate_app.cc)
target_link_libraries(cstranslate HH_OBJECTS A3M_COMPRESS)

INSTALL(TARGETS
        hhblits
        hhmake
        hhfilter
        hhsearch
        hhalign
        hhconsensus
        a3m_extract
        a3m_reduce
        a3m_database_reduce
        a3m_database_extract
        a3m_database_filter
        cstranslate
        DESTINATION bin
        )

if (OPENMP_FOUND)
    add_executable(hhblits_omp hhblits_omp.cpp)
    target_link_libraries(hhblits_omp HH_OBJECTS)

    add_executable(hhsearch_omp hhblits_omp.cpp)
    target_link_libraries(hhsearch_omp HH_OBJECTS)
    set_property(TARGET hhsearch_omp PROPERTY COMPILE_FLAGS "-DHHSEARCH=1")

    add_executable(hhalign_omp hhblits_omp.cpp)
    target_link_libraries(hhalign_omp HH_OBJECTS)
    set_property(TARGET hhalign_omp PROPERTY COMPILE_FLAGS "-DHHALIGN=1")

    add_executable(hhblits_ca3m hhblits_ca3m.cpp)
    target_link_libraries (hhblits_ca3m HH_OBJECTS)

    INSTALL(TARGETS hhblits_omp hhsearch_omp hhalign_omp hhblits_ca3m DESTINATION bin)
endif ()

if (${CHECK_MPI})
find_package(MPI QUIET)
if (MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})

    add_executable(hhblits_mpi hhblits_mpi.cpp)
    target_link_libraries(hhblits_mpi HH_OBJECTS mpq ${MPI_CXX_LIBRARIES})
    set_target_properties(hhblits_mpi PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    set_target_properties(hhblits_mpi PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    add_executable(hhsearch_mpi hhblits_mpi.cpp)
    target_link_libraries(hhsearch_mpi HH_OBJECTS mpq ${MPI_CXX_LIBRARIES})
    set_target_properties(hhsearch_mpi PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} -DHHSEARCH=1")
    set_target_properties(hhsearch_mpi PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    add_executable(hhalign_mpi hhblits_mpi.cpp)
    target_link_libraries(hhalign_mpi HH_OBJECTS mpq ${MPI_CXX_LIBRARIES})
    set_target_properties(hhalign_mpi PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} -DHHALIGN=1")
    set_target_properties(hhalign_mpi PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    add_executable(cstranslate_mpi cs/cstranslate_mpi_app.cc)
    target_link_libraries(cstranslate_mpi HH_OBJECTS A3M_COMPRESS mpq ${MPI_CXX_LIBRARIES})
    set_target_properties(cstranslate_mpi PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    set_target_properties(cstranslate_mpi PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    install(TARGETS hhblits_mpi hhsearch_mpi hhalign_mpi cstranslate_mpi DESTINATION bin)
endif ()
endif ()

