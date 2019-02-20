set(CMAKE_C_FLAGS "-std=c99 ${CMAKE_C_FLAGS}")

# sets HAVE_FMEMOPEN
add_subdirectory(ext)

add_library (ffindex ffindex.c ffutil.c)

target_include_directories (ffindex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library (ffindex_shared SHARED ffindex.c ffutil.c)

if(NOT HAVE_FMEMOPEN)
        target_link_libraries(ffindex ext)
        target_link_libraries(ffindex_shared ext)
endif()



add_executable(ffindex_apply
        ffindex_apply_mpi.c
)
target_link_libraries (ffindex_apply ffindex)
set_property(TARGET ffindex_apply PROPERTY COMPILE_FLAGS "-UHAVE_MPI")

find_package(MPI)
if(MPI_C_FOUND AND HAVE_MPI)
    add_executable(ffindex_apply_mpi
      ffindex_apply_mpi.c
    )
    target_link_libraries (ffindex_apply_mpi ffindex)
    add_subdirectory(mpq)

    set_property(TARGET ffindex_apply_mpi PROPERTY COMPILE_FLAGS "-DHAVE_MPI=1 ${MPI_C_COMPILE_FLAGS}")
    set_property(TARGET ffindex_apply_mpi PROPERTY LINK_FLAGS "${MPI_C_LINK_FLAGS}")

    target_include_directories(ffindex_apply_mpi PUBLIC ${MPI_C_INCLUDE_PATH})
    target_link_libraries (ffindex_apply_mpi mpq ${MPI_C_LIBRARIES})
    install(TARGETS
        ffindex_apply_mpi
        DESTINATION bin
    )
endif()


add_executable(ffindex_reduce
        ffindex_reduce.c
)
target_link_libraries (ffindex_reduce ffindex)

add_executable(ffindex_build
  ffindex_build.c
)
target_link_libraries (ffindex_build ffindex)


add_executable(ffindex_from_fasta
  ffindex_from_fasta.c
)
target_link_libraries (ffindex_from_fasta ffindex)


add_executable(ffindex_get
  ffindex_get.c
)
target_link_libraries (ffindex_get ffindex)


add_executable(ffindex_modify
  ffindex_modify.c
)
target_link_libraries (ffindex_modify ffindex)


add_executable(ffindex_unpack
  ffindex_unpack.c
)
target_link_libraries (ffindex_unpack ffindex)


add_executable(ffindex_order
  ffindex_order.c
)
target_link_libraries (ffindex_order ffindex)


add_executable(ffindex_from_fasta_with_split
    ffindex_from_fasta_with_split.c
)
target_link_libraries (ffindex_from_fasta_with_split ffindex)


install(PROGRAMS 
        ffindex.h 
        ffutil.h
        DESTINATION include
)

install(TARGETS ffindex
  LIBRARY DESTINATION lib64 COMPONENT libraries
        ARCHIVE DESTINATION lib64 COMPONENT libraries
)


install(TARGETS ffindex_shared
  LIBRARY DESTINATION lib64 COMPONENT libraries
        ARCHIVE DESTINATION lib64 COMPONENT libraries
)


INSTALL( 
  TARGETS ffindex_reduce
  ffindex_apply
  ffindex_build
  ffindex_from_fasta
  ffindex_get
  ffindex_modify
  ffindex_unpack
  ffindex_order
  ffindex_from_fasta_with_split
  DESTINATION bin
)
