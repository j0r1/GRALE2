include_directories("..")

add_executable(grale_lensplane_openmp grale_lensplane_openmp.cpp ../inputoutput.cpp ../communicator.cpp)
target_link_libraries(grale_lensplane_openmp gralecore-shared ${GRALE_LINK_LIBS})
add_openmp_flags(grale_lensplane_openmp)

add_executable(grale_massdens_openmp grale_massdens_openmp.cpp ../inputoutput.cpp ../communicator.cpp)
target_link_libraries(grale_massdens_openmp gralecore-shared ${GRALE_LINK_LIBS})
add_openmp_flags(grale_massdens_openmp)

install(TARGETS grale_lensplane_openmp DESTINATION bin)
install(TARGETS grale_massdens_openmp DESTINATION bin)

