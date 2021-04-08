add_definitions(-DGRALE_COMPILING)

set(HEADERS
	galens/backprojectmatrix.h
	galens/deflectionmatrix.h
	galens/fitnesshistory.h
	galens/lensinversiongafactorysingleplanecpu.h
	galens/lensinversiongafactoryparamssingleplanecpu.h
	galens/lensinversiongafactorymultiplanegpu.h
	galens/imagesbackprojector.h
	galens/lensfitnessobject.h
	galens/multifitnesshistory.h
	galens/projectedimagesinterface.h
	galens/precalculatedbackprojector.h
	galens/mpcudabackprojector.h
	galens/lensinversiongafactorycommon.h
	prob/randomnumbergenerator.h
	prob/probabilitydistribution.h
	prob/simpleuniformdistribution.h

	galens/pointgroupstorage.h
	galens/fitnesscomponent.h
	galens/fitnesscomponent_caustic.h
	galens/fitnesscomponent_null.h
	galens/fitnesscomponent_overlap.h
	galens/fitnesscomponent_time.h
	galens/lensfitnessgeneral.h
	galens/numericgradientcalculator.h
	galens/fitnesscomponent_defl.h
	galens/fitnesscomponent_dens.h
	galens/fitnesscomponent_weak.h
	galens/fitnessutil.h
	galens/inversionregistry.h
	galens/lensgaindividual.h
	galens/lensgagenomemutation.h
	galens/lensgagenomecrossover.h
	galens/lensgafitnesscomparison.h
	galens/lensgasingleobjectivecrossover.h
	galens/lensgastopcriterion.h
	galens/lensgacrossoverbase.h
	galens/lensgamultiobjectivecrossover.h
	galens/lensgacalculatorregistry.h
	galens/lensgagenomecalculator.h
	)

set(SOURCES
	galens/backprojectmatrix.cpp
	galens/deflectionmatrix.cpp
	galens/fitnesshistory.cpp
	galens/lensinversiongafactorysingleplanecpu.cpp
	galens/lensinversiongafactoryparamssingleplanecpu.cpp
	galens/lensinversiongafactorymultiplanegpu.cpp
	galens/imagesbackprojector.cpp
	galens/multifitnesshistory.cpp
	galens/projectedimagesinterface.cpp
	galens/precalculatedbackprojector.cpp
	galens/mpcudabackprojector.cpp
	galens/lensinversiongafactorycommon.cpp
	prob/randomnumbergenerator.cpp
	prob/probabilitydistribution.cpp

	galens/lensfitnessobject.cpp
	galens/lensfitnessgeneral.cpp
	galens/fitnesscomponent.cpp
	galens/fitnessutil_overlap.cpp
	galens/fitnessutil_null.cpp
	galens/fitnessutil_weak.cpp
	galens/fitnessutil_time.cpp
	galens/fitnessutil_caustic.cpp
	galens/pointgroupstorage.cpp
	galens/numericgradientcalculator.cpp
	galens/fitnesscomponent_overlap.cpp
	galens/fitnesscomponent_null.cpp
	galens/fitnesscomponent_weak.cpp
	galens/fitnesscomponent_time.cpp
	galens/fitnesscomponent_dens.cpp
	galens/fitnesscomponent_caustic.cpp
	galens/fitnesscomponent_defl.cpp	

	galens/inversionregistry.cpp
	galens/lensfitnessobject_defaults.cpp
	galens/lensgaindividual.cpp
	galens/lensgagenomemutation.cpp
	galens/lensgagenomecrossover.cpp
	galens/lensgafitnesscomparison.cpp
	galens/lensgasingleobjectivecrossover.cpp
	galens/lensgastopcriterion.cpp
	galens/lensgacrossoverbase.cpp
	galens/lensgamultiobjectivecrossover.cpp
	galens/lensgacalculatorregistry.cpp
	galens/lensgaoldfactorywrapper.cpp
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


