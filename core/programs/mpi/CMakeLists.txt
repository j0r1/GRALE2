include_directories(${MPI_INCLUDE_PATH})
include_directories("..")

add_executable(grale_lensplane_mpi grale_lensplane_mpi.cpp ../inputoutput.cpp ../communicator.cpp)
target_link_libraries(grale_lensplane_mpi ${GRALE_PROGRAMS_LIBRARY_CORE} ${GRALE_LINK_LIBS} ${MPI_LIBRARIES})
add_mpi_flags(grale_lensplane_mpi)

add_executable(grale_massdens_mpi grale_massdens_mpi.cpp ../inputoutput.cpp ../communicator.cpp)
target_link_libraries(grale_massdens_mpi ${GRALE_PROGRAMS_LIBRARY_CORE} ${GRALE_LINK_LIBS} ${MPI_LIBRARIES})
add_mpi_flags(grale_massdens_mpi)

install(TARGETS grale_lensplane_mpi DESTINATION bin)
install(TARGETS grale_massdens_mpi DESTINATION bin)
