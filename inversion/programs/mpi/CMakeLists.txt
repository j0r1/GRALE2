include_directories("..")
include_directories("${PROJECT_SOURCE_DIR}/core/programs/")

find_package(MPI)
if (MPI_FOUND)

	add_executable(grale_invert_mpi grale_invert_mpi.cpp ../inversioncommunicator.cpp
		"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
		)
	target_link_libraries(grale_invert_mpi graleinv-static gralecore-static ${GRALE_LINK_LIBS} ${MPI_LIBRARIES})
	add_mpi_flags(grale_invert_mpi)

	install(TARGETS grale_invert_mpi DESTINATION bin)
endif ()
