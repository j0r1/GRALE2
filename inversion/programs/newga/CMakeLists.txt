include_directories("..")
include_directories("${PROJECT_SOURCE_DIR}/core/programs/")

add_executable(grale_invert_newga grale_invert_newga.cpp ../inversioncommunicator.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_newga graleinv-shared gralecore-shared ${GRALE_LINK_LIBS})

install(TARGETS grale_invert_newga DESTINATION bin)

