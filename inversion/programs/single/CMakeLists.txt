include_directories("..")
include_directories("${PROJECT_SOURCE_DIR}/core/programs/")

add_executable(grale_invert_single grale_invert_single.cpp ../inversioncommunicator.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_single graleinv-shared gralecore-shared ${GRALE_LINK_LIBS})

add_executable(grale_invert_usage_new grale_invert_usage_new.cpp ../inversioncommunicatornewga.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_usage_new graleinv-shared gralecore-shared ${GRALE_LINK_LIBS})

add_executable(grale_invert_confparamdefaults_new grale_invert_confparamdefaults_new.cpp ../inversioncommunicatornewga.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_confparamdefaults_new graleinv-shared gralecore-shared ${GRALE_LINK_LIBS})


add_executable(grale_invert_calcfitness grale_invert_calcfitness.cpp ../inversioncommunicator.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_calcfitness graleinv-shared gralecore-shared ${GRALE_LINK_LIBS})

install(TARGETS grale_invert_single DESTINATION bin)
install(TARGETS grale_invert_usage_new DESTINATION bin)
install(TARGETS grale_invert_confparamdefaults_new DESTINATION bin)
install(TARGETS grale_invert_calcfitness DESTINATION bin)


