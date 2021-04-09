include_directories("..")
include_directories("${PROJECT_SOURCE_DIR}/core/programs/")

add_executable(grale_invert_newga grale_invert_newga.cpp ../inversioncommunicatornewga.cpp
		"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
)
target_link_libraries(grale_invert_newga graleinv-static gralecore-static ${GRALE_LINK_LIBS})

install(TARGETS grale_invert_newga DESTINATION bin)

if (GRALE_SUPPORT_MPI AND EATK_SUPPORT_MPI)
	add_executable(grale_invert_newgampi grale_invert_newgampi.cpp ../inversioncommunicatornewga.cpp
			"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
	target_link_libraries(grale_invert_newgampi graleinv-static gralecore-static ${GRALE_LINK_LIBS})

	install(TARGETS grale_invert_newgampi DESTINATION bin)
endif()

