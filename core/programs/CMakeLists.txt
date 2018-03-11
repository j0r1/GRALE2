apply_include_paths("${GRALE_INTERNAL_INCLUDES}")
apply_include_paths("${GRALE_EXTERNAL_INCLUDES}")

find_package(Threads)
if (THREADS_FOUND)
	add_subdirectory(threads)
endif()

find_package(MPI)
if (MPI_FOUND)
	add_subdirectory(mpi)
endif()

find_package(OpenCL)
if (OPENCL_FOUND)
	add_subdirectory(opencl)
endif()
