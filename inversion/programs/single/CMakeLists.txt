include_directories("..")
include_directories("${PROJECT_SOURCE_DIR}/core/programs/")

add_executable(grale_invert_single grale_invert_single.cpp ../inversioncommunicator.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_single graleinv-shared gralecore-shared ${GRALE_LINK_LIBS})

add_executable(grale_invert_usage grale_invert_usage.cpp ../inversioncommunicator.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_usage graleinv-shared gralecore-shared ${GRALE_LINK_LIBS})

add_executable(grale_invert_confparamdefaults grale_invert_confparamdefaults.cpp ../inversioncommunicator.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_confparamdefaults graleinv-shared gralecore-shared ${GRALE_LINK_LIBS})

add_executable(grale_invert_calcfitness grale_invert_calcfitness.cpp ../inversioncommunicator.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_calcfitness graleinv-shared gralecore-shared ${GRALE_LINK_LIBS})

install(TARGETS grale_invert_single DESTINATION bin)
install(TARGETS grale_invert_usage DESTINATION bin)
install(TARGETS grale_invert_confparamdefaults DESTINATION bin)
install(TARGETS grale_invert_calcfitness DESTINATION bin)


