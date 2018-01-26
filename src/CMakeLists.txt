if (${HAVE_AVX2})
    ADD_DEFINITIONS("-DAVX2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -Wa,-q")
elseif (${HAVE_SSE2})
    ADD_DEFINITIONS("-DSSE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -msse2")
else ()
    include(CheckSSEFeatures)
    ADD_DEFINITIONS(${CXX_DFLAGS})
    if (${HAVE_AVX2_EXTENSIONS})
        ADD_DEFINITIONS("-DAVX2")
    elseif (${HAVE_SSE2_EXTENSIONS})
        ADD_DEFINITIONS("-DSSE")
    else ()
        message(FATAL_ERROR "At least SSE2 instruction set support is required! CMake will exit.")
    endif (${HAVE_AVX2_EXTENSIONS})
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
        simd.h
        )

add_library(HH_OBJECTS ${HH_SOURCE})

add_library(hhviterbialgorithm_with_celloff hhviterbialgorithm.cpp)
set_property(TARGET hhviterbialgorithm_with_celloff PROPERTY COMPILE_FLAGS "-DVITERBI_CELLOFF=1")

add_library(hhviterbialgorithm_with_celloff_and_ss hhviterbialgorithm.cpp)
set_property(TARGET hhviterbialgorithm_with_celloff_and_ss PROPERTY COMPILE_FLAGS "-DVITERBI_CELLOFF=1 -DVITERBI_SS_SCORE=1")

add_library(hhviterbialgorithm_and_ss hhviterbialgorithm.cpp)
set_property(TARGET hhviterbialgorithm_and_ss PROPERTY COMPILE_FLAGS "-DVITERBI_SS_SCORE=1")

add_executable(hhblits hhblits_app.cpp)
target_link_libraries(hhblits
        ffindex CS_OBJECTS HH_OBJECTS
        hhviterbialgorithm_with_celloff
        hhviterbialgorithm_and_ss
        hhviterbialgorithm_with_celloff_and_ss
        )

add_executable(hhmake hhmake.cpp)
target_link_libraries(hhmake
        ffindex CS_OBJECTS HH_OBJECTS
        hhviterbialgorithm_with_celloff
        hhviterbialgorithm_and_ss
        hhviterbialgorithm_with_celloff_and_ss
        )
add_executable(hhfilter hhfilter.cpp)
target_link_libraries(hhfilter
        ffindex CS_OBJECTS HH_OBJECTS
        hhviterbialgorithm_with_celloff
        hhviterbialgorithm_and_ss
        hhviterbialgorithm_with_celloff_and_ss
        )

add_executable(hhsearch hhsearch_app.cpp)
target_link_libraries(hhsearch
        ffindex CS_OBJECTS HH_OBJECTS
        hhviterbialgorithm_with_celloff
        hhviterbialgorithm_and_ss
        hhviterbialgorithm_with_celloff_and_ss
        )

add_executable(hhalign hhalign_app.cpp)
target_link_libraries(hhalign
        ffindex CS_OBJECTS HH_OBJECTS
        hhviterbialgorithm_with_celloff
        hhviterbialgorithm_and_ss
        hhviterbialgorithm_with_celloff_and_ss
        )

add_executable(hhconsensus hhconsensus.cpp)
target_link_libraries(hhconsensus
        ffindex CS_OBJECTS HH_OBJECTS
        hhviterbialgorithm_with_celloff
        hhviterbialgorithm_and_ss
        hhviterbialgorithm_with_celloff_and_ss
        )

add_library(A3M_COMPRESS a3m_compress.cpp)

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
target_link_libraries(cstranslate
        ffindex CS_OBJECTS HH_OBJECTS A3M_COMPRESS
        hhviterbialgorithm_with_celloff
        hhviterbialgorithm_and_ss
        hhviterbialgorithm_with_celloff_and_ss
        )

INSTALL(PROGRAMS
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
    target_link_libraries(hhblits_omp
            ffindex CS_OBJECTS HH_OBJECTS
            hhviterbialgorithm_with_celloff
            hhviterbialgorithm_and_ss
            hhviterbialgorithm_with_celloff_and_ss
            )
    INSTALL(PROGRAMS hhblits_omp DESTINATION bin)

    add_executable(hhsearch_omp hhsearch_omp.cpp)
    target_link_libraries(hhsearch_omp
            ffindex CS_OBJECTS HH_OBJECTS
            hhviterbialgorithm_with_celloff
            hhviterbialgorithm_and_ss
            hhviterbialgorithm_with_celloff_and_ss
            )
    INSTALL(PROGRAMS hhsearch_omp DESTINATION bin)
endif ()

find_package(MPI QUIET)
if (MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})

    add_executable(hhblits_mpi hhblits_mpi.cpp)
    target_link_libraries(hhblits_mpi
            ffindex mpq CS_OBJECTS HH_OBJECTS
            hhviterbialgorithm_with_celloff
            hhviterbialgorithm_and_ss
            hhviterbialgorithm_with_celloff_and_ss
            ${MPI_CXX_LIBRARIES}
            )
    set_target_properties(hhblits_mpi PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    set_target_properties(hhblits_mpi PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    install(PROGRAMS hhblits_mpi DESTINATION bin)

    add_executable(cstranslate_mpi cs/cstranslate_mpi_app.cc)
    target_link_libraries(cstranslate_mpi
            ffindex mpq CS_OBJECTS HH_OBJECTS A3M_COMPRESS
            hhviterbialgorithm_with_celloff
            hhviterbialgorithm_and_ss
            hhviterbialgorithm_with_celloff_and_ss
            ${MPI_CXX_LIBRARIES}
            )
    set_target_properties(cstranslate_mpi PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    set_target_properties(cstranslate_mpi PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    install(PROGRAMS cstranslate_mpi DESTINATION bin)
endif ()

