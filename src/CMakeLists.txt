list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

#Sanitizers
if(${HAVE_SANITIZER})
	include(FindUBSan)
	include(FindASan)
	include(FindMSan)
	include(FindTSan)
endif()


if(X86 OR X86_64)

#SSE
if(${HAVE_AVX2})
    ADD_DEFINITIONS("-DAVX2")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -Wa,-q")
elseif(${HAVE_SSE2})
    ADD_DEFINITIONS("-DSSE")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -msse2")
else()
    include(CheckSSEFeatures)
    ADD_DEFINITIONS(${CXX_DFLAGS})
    if(${HAVE_AVX2_EXTENSIONS})
        ADD_DEFINITIONS("-DAVX2")
    elseif(${HAVE_SSE2_EXTENSIONS})
        ADD_DEFINITIONS("-DSSE")
    else()
        message( FATAL_ERROR "SSE2 is needed to run compile! CMake will exit." )
    endif(${HAVE_AVX2_EXTENSIONS})
endif()

else() #

message("!! NOT X86 OR X86_64 !!")

endif() # NOT X86 OR X86_64



find_package(OpenMP)
if (OPENMP_FOUND)
	message("-- Found OpenMP")
	add_definitions(-DOPENMP)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing")

add_subdirectory(cs)
include_directories(cs)
include_directories("${PROJECT_BINARY_DIR}")

set(HH_SOURCE
	${HH_SOURCE}
	${CMAKE_CURRENT_SOURCE_DIR}/hhblits.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhblits.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhdecl.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhdecl.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhhit.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhhit-inl.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhhit.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhmatrices.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhmatrices.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhsearch.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhsearch.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhalign.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhalign.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhhitlist.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhhitlist-inl.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhhitlist.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhposteriordecoder.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhposteriordecoder.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/hhutil.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhutil-inl.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhutil.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/util.h
	${CMAKE_CURRENT_SOURCE_DIR}/util.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhalignment.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhalignment.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhforwardalgorithm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhhmm.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhhmm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhposteriordecoderrunner.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhposteriordecoderrunner.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhviterbialgorithm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhfullalignment.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhfullalignment.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhhmmsimd.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhhmmsimd.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhposteriormatrix.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhposteriormatrix.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhviterbi.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhviterbi.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhbacktracemac.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/hhmacalgorithm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhprefilter.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhprefilter.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhviterbimatrix.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhviterbimatrix-inl.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhviterbimatrix.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/hhbackwardalgorithm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhdatabase.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhdatabase.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhhalfalignment.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhhalfalignment.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhviterbirunner.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhviterbirunner.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/hhfunc.h
	${CMAKE_CURRENT_SOURCE_DIR}/hhfunc.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/simd.h
)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/hhsuite_config.h.in"
  "${PROJECT_BINARY_DIR}/hhsuite_config.h"
)

add_library(HH_OBJECTS OBJECT ${HH_SOURCE})

add_library(hhviterbialgorithm_with_celloff OBJECT hhviterbialgorithm.cpp)
set_property(TARGET hhviterbialgorithm_with_celloff PROPERTY COMPILE_FLAGS "-DVITERBI_CELLOFF=1")

add_library(hhviterbialgorithm_with_celloff_and_ss OBJECT hhviterbialgorithm.cpp)
set_property(TARGET hhviterbialgorithm_with_celloff_and_ss PROPERTY COMPILE_FLAGS "-DVITERBI_CELLOFF=1 -DVITERBI_SS_SCORE=1")

add_library(hhviterbialgorithm_and_ss OBJECT hhviterbialgorithm.cpp)
set_property(TARGET hhviterbialgorithm_and_ss PROPERTY COMPILE_FLAGS "-DVITERBI_SS_SCORE=1")


add_library(CS_OBJECTS OBJECT ${CS_SOURCE})

# add the executable
add_executable(hhblits 
    hhblits_app.cpp	$<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS>
    $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
    $<TARGET_OBJECTS:hhviterbialgorithm_and_ss>
    $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff_and_ss>
)
target_link_libraries (hhblits ffindex)

if (OPENMP_FOUND)
    add_executable(hhblits_omp
        hhblits_omp.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS>
        $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
        $<TARGET_OBJECTS:hhviterbialgorithm_and_ss>
        $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff_and_ss>
    )
    target_link_libraries (hhblits_omp ffindex)

    add_executable(hhsearch_omp
        hhsearch_omp.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS>
        $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
        $<TARGET_OBJECTS:hhviterbialgorithm_and_ss>
        $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff_and_ss>
    )
    target_link_libraries (hhsearch_omp ffindex)

    INSTALL(PROGRAMS 
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hhblits_omp
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hhsearch_omp
        DESTINATION bin
    )
endif()

add_executable(hhmake
  hhmake.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS>
  $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
  $<TARGET_OBJECTS:hhviterbialgorithm_and_ss>
  $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff_and_ss>
)
target_link_libraries (hhmake ffindex)

add_executable(hhfilter
  hhfilter.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS>
  $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
  $<TARGET_OBJECTS:hhviterbialgorithm_and_ss>
  $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff_and_ss>
)
target_link_libraries (hhfilter ffindex)

add_executable(hhsearch
  hhsearch_app.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS>
  $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
  $<TARGET_OBJECTS:hhviterbialgorithm_and_ss>
  $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff_and_ss>
)
target_link_libraries (hhsearch ffindex)

add_executable(hhalign
  hhalign_app.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS>
  $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
  $<TARGET_OBJECTS:hhviterbialgorithm_and_ss>
  $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff_and_ss>
)
target_link_libraries (hhalign ffindex)

add_executable(hhconsensus
  hhconsensus.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS>
  $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
  $<TARGET_OBJECTS:hhviterbialgorithm_and_ss>
  $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff_and_ss>
)
target_link_libraries (hhconsensus ffindex)


add_library(A3M_COMPRESS_OBJECT OBJECT a3m_compress.cpp)

add_executable(a3m_extract
	a3m_extract.cpp $<TARGET_OBJECTS:A3M_COMPRESS_OBJECT>
)
target_link_libraries (a3m_extract ffindex)

add_executable(a3m_reduce
  a3m_reduce.cpp $<TARGET_OBJECTS:A3M_COMPRESS_OBJECT>
)
target_link_libraries (a3m_reduce ffindex)

add_executable(a3m_database_reduce
  a3m_database_reduce.cpp $<TARGET_OBJECTS:A3M_COMPRESS_OBJECT>
)
target_link_libraries (a3m_database_reduce ffindex)

add_executable(a3m_database_extract
  a3m_database_extract.cpp $<TARGET_OBJECTS:A3M_COMPRESS_OBJECT>
)
target_link_libraries (a3m_database_extract ffindex)

add_executable(a3m_database_filter
  a3m_database_filter.cpp $<TARGET_OBJECTS:A3M_COMPRESS_OBJECT>
)
target_link_libraries (a3m_database_filter ffindex)

add_executable(cstranslate
  cs/cstranslate_app.cc
        $<TARGET_OBJECTS:CS_OBJECTS>
        $<TARGET_OBJECTS:HH_OBJECTS>
        $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
        $<TARGET_OBJECTS:hhviterbialgorithm_and_ss>
        $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff_and_ss>
        $<TARGET_OBJECTS:A3M_COMPRESS_OBJECT>
)
target_link_libraries (cstranslate ffindex)

find_package(MPI QUIET)

if(MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})

    #hhblits_mpi
    add_executable(hhblits_mpi
	         hhblits_mpi.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS>
	         $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
	         $<TARGET_OBJECTS:hhviterbialgorithm_and_ss>
	         $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff_and_ss>
	       )
    target_link_libraries(hhblits_mpi ffindex mpq)
    target_link_libraries(hhblits_mpi ${MPI_CXX_LIBRARIES})
    set_target_properties(hhblits_mpi PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    set_target_properties(hhblits_mpi PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    install(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hhblits_mpi DESTINATION bin)

    #cstranslate_mpi
    add_executable(cstranslate_mpi
            cs/cstranslate_mpi_app.cc
            $<TARGET_OBJECTS:CS_OBJECTS>
            $<TARGET_OBJECTS:HH_OBJECTS>
            $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
            $<TARGET_OBJECTS:hhviterbialgorithm_and_ss>
            $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff_and_ss>
            $<TARGET_OBJECTS:A3M_COMPRESS_OBJECT>
            )
    target_link_libraries(cstranslate_mpi ffindex mpq)
    target_link_libraries(cstranslate_mpi ${MPI_CXX_LIBRARIES})
    set_target_properties(cstranslate_mpi PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    set_target_properties(cstranslate_mpi PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")

    install(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cstranslate_mpi DESTINATION bin)
endif()

INSTALL(PROGRAMS 
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hhblits 
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hhmake 
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hhfilter 
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hhsearch 
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hhalign 
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hhconsensus 
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/a3m_extract
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/a3m_database_reduce
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/a3m_database_extract
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/a3m_database_filter 
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cstranslate
	DESTINATION bin
)

