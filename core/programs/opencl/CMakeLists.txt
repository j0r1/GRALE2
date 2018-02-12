include_directories(${OPENCL_INCLUDE_DIRS})
check_include_file(OpenCL/cl.h GRALE_FOUND_OPENCLCL_H)
if (NOT GRALE_FOUND_OPENCLCL_H)
	check_include_file(CL/cl.h GRALE_FOUND_CLCL_H)
	if (NOT GRALE_FOUND_CLCL_H)
		message(WARNING "OpenCL detected but neither OpenCL/cl.h nor CL/cl.h were found, assuming CL/cl.h")
	endif (NOT GRALE_FOUND_CLCL_H)
else()
	add_definitions("-DGRALE_OPENCL_OPENCL_CL_H")
endif (NOT GRALE_FOUND_OPENCLCL_H)

include_directories("..")

add_executable(grale_lensplane_opencl grale_lensplane_opencl.cpp openclkernel.cpp ../inputoutput.cpp ../communicator.cpp)
target_link_libraries(grale_lensplane_opencl gralecore-shared ${GRALE_LINK_LIBS} ${OPENCL_LIBRARIES})

install(TARGETS grale_lensplane_opencl DESTINATION bin)

