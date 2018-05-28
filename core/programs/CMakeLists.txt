apply_include_paths("${GRALE_INTERNAL_INCLUDES}")
apply_include_paths("${GRALE_EXTERNAL_INCLUDES}")

find_package(Threads)
if (THREADS_FOUND)
	add_subdirectory(threads)
endif()

if (GRALE_SUPPORT_MPI)
	add_subdirectory(mpi)
endif()

add_subdirectory(opencl)
