apply_include_paths("${GRALE_INTERNAL_INCLUDES}")
apply_include_paths("${GRALE_EXTERNAL_INCLUDES}")

add_subdirectory(single)
add_subdirectory(clientserver)

find_package(MPI)
if (MPI_FOUND)
	add_subdirectory(mpi)
endif()

