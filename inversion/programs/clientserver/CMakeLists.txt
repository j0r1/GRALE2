include_directories("..")
include_directories("${PROJECT_SOURCE_DIR}/core/programs/")

add_executable(grale_invert_clientserver grale_invert_clientserver.cpp ../inversioncommunicator.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_clientserver graleinv-shared gralecore-static ${GRALE_LINK_LIBS})

install(TARGETS grale_invert_clientserver DESTINATION bin)

