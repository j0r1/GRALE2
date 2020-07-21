add_definitions(-DGRALE_COMPILING)

set(HEADERS
	galens/backprojectmatrixnew.h
	galens/deflectionmatrix.h
	galens/fitnesshistory.h
	galens/galensmodule.h
	galens/gridlensinversiongafactorybase.h
	galens/gridlensinversiongafactoryparams.h
	galens/gridlensinversiongenomebase.h
	galens/imagesbackprojector.h
	galens/lensfitnessobject.h
	galens/multifitnesshistory.h
	galens/projectedimagesinterface.h
	galens/precalculatedbackprojector.h
	prob/randomnumbergenerator.h
	prob/probabilitydistribution.h
	prob/simpleuniformdistribution.h
	)

set(SOURCES
	galens/backprojectmatrixnew.cpp
	galens/deflectionmatrix.cpp
	galens/fitnesshistory.cpp
	galens/galensmodule.cpp
	galens/gridlensinversiongafactorybase.cpp
	galens/gridlensinversiongafactoryparams.cpp
	galens/gridlensinversiongenomebase.cpp
	galens/imagesbackprojector.cpp
	galens/multifitnesshistory.cpp
	galens/projectedimagesinterface.cpp
	galens/precalculatedbackprojector.cpp
	prob/randomnumbergenerator.cpp
	prob/probabilitydistribution.cpp
	)

if (MSVC)
	set(CMAKE_DEBUG_POSTFIX _d)
endif (MSVC)

if (NOT MSVC OR GRALE_COMPILE_STATIC)
	add_library(graleinv-static STATIC ${SOURCES} ${HEADERS})
	set_target_properties(graleinv-static PROPERTIES OUTPUT_NAME grale2inv)
	set_target_properties(graleinv-static PROPERTIES CLEAN_DIRECT_OUTPUT 1)
	set(GRALE_INSTALLTARGETS graleinv-static)
	target_link_libraries(graleinv-static ${GRALE_LINK_LIBS})
endif (NOT MSVC OR GRALE_COMPILE_STATIC)

if (NOT MSVC OR NOT GRALE_COMPILE_STATIC)
	add_library(graleinv-shared SHARED ${SOURCES} ${HEADERS})
	set_target_properties(graleinv-shared PROPERTIES OUTPUT_NAME grale2inv)
	set_target_properties(graleinv-shared PROPERTIES VERSION ${VERSION})
	set_target_properties(graleinv-shared PROPERTIES CLEAN_DIRECT_OUTPUT 1)
	set(GRALE_INSTALLTARGETS ${GRALE_INSTALLTARGETS} graleinv-shared)
	if (GRALE_LINK)
		target_link_libraries(graleinv-shared gralecore-shared ${GRALE_LINK_LIBS})
	endif (GRALE_LINK)
endif (NOT MSVC OR NOT GRALE_COMPILE_STATIC)

apply_include_paths("${GRALE_INTERNAL_INCLUDES}")
apply_include_paths("${GRALE_EXTERNAL_INCLUDES}")

install(FILES ${HEADERS} DESTINATION include/grale)
install(TARGETS ${GRALE_INSTALLTARGETS} DESTINATION ${LIBRARY_INSTALL_DIR})


