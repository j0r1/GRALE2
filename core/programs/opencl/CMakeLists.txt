include_directories("..")

add_executable(grale_lensplane_opencl grale_lensplane_opencl.cpp ../inputoutput.cpp ../communicator.cpp)
target_link_libraries(grale_lensplane_opencl ${GRALE_PROGRAMS_LIBRARY_CORE} ${GRALE_LINK_LIBS} ${GRALE_DLLIBS})

add_executable(ocltst ocltst.cpp)
target_link_libraries(ocltst ${GRALE_PROGRAMS_LIBRARY_CORE} ${GRALE_LINK_LIBS} ${GRALE_DLLIBS})

install(TARGETS grale_lensplane_opencl DESTINATION bin)
