apply_include_paths("${GRALE_INTERNAL_INCLUDES}")
apply_include_paths("${GRALE_EXTERNAL_INCLUDES}")

add_subdirectory(single)
add_subdirectory(clientserver)
add_subdirectory(newga)

if (GRALE_SUPPORT_MPI)
	add_subdirectory(mpi)
endif()

