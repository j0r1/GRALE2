include_directories("..")
include_directories("${PROJECT_SOURCE_DIR}/core/programs/")

add_executable(grale_invert_newga grale_invert_newga.cpp ../inversioncommunicatornewga.cpp
		"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
)
target_link_libraries(grale_invert_newga ${GRALE_PROGRAMS_LIBRARY_INV} ${GRALE_PROGRAMS_LIBRARY_CORE} ${GRALE_LINK_LIBS})

add_executable(grale_invert_parametric grale_invert_parametric.cpp ../inversioncommunicatornewga.cpp
		"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
)
target_link_libraries(grale_invert_parametric ${GRALE_PROGRAMS_LIBRARY_INV} ${GRALE_PROGRAMS_LIBRARY_CORE} ${GRALE_LINK_LIBS})

install(TARGETS grale_invert_newga DESTINATION bin)
install(TARGETS grale_invert_parametric DESTINATION bin)

if (GRALE_SUPPORT_MPI AND EATK_SUPPORT_MPI)
	add_executable(grale_invert_newgampi grale_invert_newgampi.cpp ../inversioncommunicatornewga.cpp
			"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
	target_link_libraries(grale_invert_newgampi ${GRALE_PROGRAMS_LIBRARY_INV} ${GRALE_PROGRAMS_LIBRARY_CORE} ${GRALE_LINK_LIBS})

	add_executable(grale_invert_parametricmpi grale_invert_parametricmpi.cpp ../inversioncommunicatornewga.cpp
			"${PROJECT_SOURCE_DIR}/core/programs/inputoutput.cpp"
	)
	target_link_libraries(grale_invert_parametricmpi ${GRALE_PROGRAMS_LIBRARY_INV} ${GRALE_PROGRAMS_LIBRARY_CORE} ${GRALE_LINK_LIBS})

	install(TARGETS grale_invert_newgampi DESTINATION bin)
	install(TARGETS grale_invert_parametricmpi DESTINATION bin)
endif()

