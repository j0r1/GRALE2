include_directories("..")

add_executable(grale_lensplane_threads grale_lensplane_threads ../inputoutput.cpp ../communicator.cpp)
target_link_libraries(grale_lensplane_threads gralecore-static ${GRALE_LINK_LIBS} ${CMAKE_THREAD_LIBS_INIT})

add_executable(grale_massdens_threads grale_massdens_threads.cpp ../inputoutput.cpp ../communicator.cpp)
target_link_libraries(grale_massdens_threads gralecore-static ${GRALE_LINK_LIBS} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS grale_lensplane_threads DESTINATION bin)
install(TARGETS grale_massdens_threads DESTINATION bin)

