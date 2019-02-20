find_package(MPI)
IF(MPI_FOUND)
	set(MPQ_SRC
	  mpq.c
	  mpq.h
	)
	add_library(mpq ${MPQ_SRC})
	set_property(TARGET mpq PROPERTY COMPILE_FLAGS ${MPI_LINK_FLAGS})
	target_include_directories(mpq PUBLIC ${MPI_INCLUDE_PATH})
ENDIF()