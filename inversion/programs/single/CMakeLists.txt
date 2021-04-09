include_directories("..")
include_directories("${PROJECT_SOURCE_DIR}/core/programs/")

add_executable(grale_invert_usage_new grale_invert_usage_new.cpp ../inversioncommunicatornewga.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_usage_new graleinv-static gralecore-static ${GRALE_LINK_LIBS})

add_executable(grale_invert_confparamdefaults_new grale_invert_confparamdefaults_new.cpp ../inversioncommunicatornewga.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_confparamdefaults_new graleinv-static gralecore-static ${GRALE_LINK_LIBS})

add_executable(grale_invert_calcfitness_new grale_invert_calcfitness_new.cpp ../inversioncommunicatornewga.cpp
	"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
target_link_libraries(grale_invert_calcfitness_new graleinv-static gralecore-static ${GRALE_LINK_LIBS})

install(TARGETS grale_invert_usage_new DESTINATION bin)
install(TARGETS grale_invert_confparamdefaults_new DESTINATION bin)
install(TARGETS grale_invert_calcfitness_new DESTINATION bin)


