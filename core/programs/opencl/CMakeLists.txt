include_directories("..")

check_cxx_source_compiles("#include <windows.h>\nint main(void)\n{\nHMODULE x = LoadLibrary(\"x.dll\");\nreturn 0;\n}\n" TMP_LOADLIBRARY_CHECK)
if (TMP_LOADLIBRARY_CHECK)
	add_definitions(-DGRALE_LOADLIBRARY)
	set(GRALE_DLLIBS "")	
else (TMP_LOADLIBRARY_CHECK)
	set(GRALE_DLLIBS "-ldl")	
endif (TMP_LOADLIBRARY_CHECK)

add_executable(grale_lensplane_opencl grale_lensplane_opencl.cpp openclkernel.cpp ../inputoutput.cpp ../communicator.cpp)
target_link_libraries(grale_lensplane_opencl gralecore-static ${GRALE_LINK_LIBS} ${GRALE_DLLIBS})

install(TARGETS grale_lensplane_opencl DESTINATION bin)
