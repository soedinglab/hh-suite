#SSE
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
include(CheckSSEFeatures)
ADD_DEFINITIONS(${CXX_DFLAGS})
if(${HAVE_AVX2_EXTENSIONS})
    ADD_DEFINITIONS("-DAVX2")
else()
if(${HAVE_SSE2_EXTENSIONS})
    ADD_DEFINITIONS("-DSSE")
endif(${HAVE_SSE2_EXTENSIONS})
endif(${HAVE_AVX2_EXTENSIONS})


include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-fopenmp SUPPORT_OMP)
if (SUPPORT_OMP)
    message("-- Found OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DOPENMP -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOPENMP -fopenmp")
endif()


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

set (CMAKE_CXX_FLAGS "-fno-strict-aliasing")

add_library(HH_OBJECTS OBJECT ${HH_SOURCE})

add_library(hhviterbialgorithm_with_celloff OBJECT hhviterbialgorithm.cpp)
set_property(TARGET hhviterbialgorithm_with_celloff PROPERTY COMPILE_FLAGS "-DVITERBI_CELLOFF=1")

add_library(CS_OBJECTS OBJECT ${CS_SOURCE})

# add the executable
add_executable(hhblits 
	hhblits_app.cpp	$<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS> $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff> 
)
target_link_libraries (hhblits ffindex)
if (SUPPORT_OMP)
add_executable(hhblits_omp
  hhblits_omp.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS> $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
)
target_link_libraries (hhblits_omp ffindex)
endif()
find_package(MPI QUIET)
IF(MPI_FOUND)
	add_executable(hhblits_mpi
  	hhblits_mpi.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS> $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
	)
	set_property(TARGET hhblits_mpi PROPERTY COMPILE_FLAGS ${MPI_LINK_FLAGS})
	target_include_directories(hhblits_mpi PUBLIC ${MPI_INCLUDE_PATH})
	target_link_libraries (hhblits_mpi ${MPI_LIBRARIES} ffindex) 

	INSTALL(PROGRAMS 
  	${CMAKE_CURRENT_BINARY_DIR}/hhblits_mpi 
  	DESTINATION bin
	)
ENDIF()

add_executable(hhmake
  hhmake.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS> $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
)
target_link_libraries (hhmake ffindex)

add_executable(hhfilter
  hhfilter.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS> $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
)
target_link_libraries (hhfilter ffindex)

add_executable(hhsearch
  hhsearch_app.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS> $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
)
target_link_libraries (hhsearch ffindex)

add_executable(hhalign
  hhalign_app.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS> $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
)
target_link_libraries (hhalign ffindex)

add_executable(hhconsensus
  hhconsensus.cpp $<TARGET_OBJECTS:CS_OBJECTS> $<TARGET_OBJECTS:HH_OBJECTS> $<TARGET_OBJECTS:hhviterbialgorithm_with_celloff>
)
target_link_libraries (hhconsensus ffindex)


add_library(A3M_COMPRESS_OBJECT OBJECT a3m_compress.cpp)

add_executable(a3m_extractor
	a3m_extractor.cpp $<TARGET_OBJECTS:A3M_COMPRESS_OBJECT>
)
target_link_libraries (a3m_extractor ffindex)

add_executable(a3m_reducer
  a3m_reducer.cpp $<TARGET_OBJECTS:A3M_COMPRESS_OBJECT>
)
target_link_libraries (a3m_reducer ffindex)

add_executable(a3m_database_reducer
  a3m_database_reducer.cpp $<TARGET_OBJECTS:A3M_COMPRESS_OBJECT>
)
target_link_libraries (a3m_database_reducer ffindex)

add_executable(a3m_database_filter
  a3m_database_filter.cpp $<TARGET_OBJECTS:A3M_COMPRESS_OBJECT>
)
target_link_libraries (a3m_database_filter ffindex)

add_executable(cstranslate
  cs/cstranslate_app.cc $<TARGET_OBJECTS:CS_OBJECTS>
)

INSTALL(PROGRAMS 
	${CMAKE_CURRENT_BINARY_DIR}/hhblits 
	${CMAKE_CURRENT_BINARY_DIR}/hhblits_omp 
	${CMAKE_CURRENT_BINARY_DIR}/hhmake 
	${CMAKE_CURRENT_BINARY_DIR}/hhfilter 
	${CMAKE_CURRENT_BINARY_DIR}/hhsearch 
	${CMAKE_CURRENT_BINARY_DIR}/hhalign 
	${CMAKE_CURRENT_BINARY_DIR}/hhconsensus 
	${CMAKE_CURRENT_BINARY_DIR}/a3m_extractor 
	${CMAKE_CURRENT_BINARY_DIR}/a3m_database_reducer
	${CMAKE_CURRENT_BINARY_DIR}/a3m_database_filter 
	${CMAKE_CURRENT_BINARY_DIR}/cstranslate 
	DESTINATION bin
)
