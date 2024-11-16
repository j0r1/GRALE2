add_definitions(-DGRALE_COMPILING)

set(HEADERS
	galens/backprojectmatrix.h
	galens/deflectionmatrix.h
	galens/fitnesshistory.h
	galens/lensinversiongafactorysingleplanecpu.h
	galens/lensinversiongafactorymultiplanegpu.h
	galens/imagesbackprojector.h
	galens/multifitnesshistory.h
	galens/projectedimagesinterface.h
	galens/precalculatedbackprojector.h
	galens/mpcudabackprojector.h
	galens/lensinversiongafactorycommon.h
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
	galens/openclcalculator.h
	galens/oclcalculatedbackprojector.h
	galens/lensdemutationcrossover.h
	galens/lensdeevolver.h
	galens/populationdump.h
	galens/lensnsga2evolver.h
	galens/eaevolverfunctions.h
	galens/openclsingleplanedeflection.h
	galens/lensgaparametricsingleplanecalculator.h
	prob/randomnumbergenerator.h
	prob/xoshiro128plus.h
	fitness/lensfitnessobject.h
	fitness/numericgradientcalculator.h
	fitness/pointgroupstorage.h
	fitness/fitnesscomponent.h
	fitness/fitnesscomponent_caustic.h
	fitness/fitnesscomponent_null.h
	fitness/fitnesscomponent_overlap.h
	fitness/fitnesscomponent_time.h
	fitness/lensfitnessgeneral.h
	fitness/fitnesscomponent_defl.h
	fitness/fitnesscomponent_dens.h
	fitness/fitnesscomponent_weak.h
	fitness/fitnesscomponent_parity.h
	fitness/fitnessutil.h
	)

set(SOURCES
	galens/backprojectmatrix.cpp
	galens/deflectionmatrix.cpp
	galens/fitnesshistory.cpp
	galens/lensinversiongafactorysingleplanecpu.cpp
	galens/lensinversiongafactorymultiplanegpu.cpp
	galens/imagesbackprojector.cpp
	galens/multifitnesshistory.cpp
	galens/projectedimagesinterface.cpp
	galens/precalculatedbackprojector.cpp
	galens/mpcudabackprojector.cpp
	galens/lensinversiongafactorycommon.cpp
	galens/lensgaindividual.cpp
	galens/lensgagenomemutation.cpp
	galens/lensgagenomecrossover.cpp
	galens/lensgafitnesscomparison.cpp
	galens/lensgasingleobjectivecrossover.cpp
	galens/lensgastopcriterion.cpp
	galens/lensgacrossoverbase.cpp
	galens/lensgamultiobjectivecrossover.cpp
	galens/lensgacalculatorregistry.cpp
	galens/openclcalculator.cpp
	galens/oclcalculatedbackprojector.cpp
	galens/lensdemutationcrossover.cpp
	galens/lensdeevolver.cpp
	galens/populationdump.cpp
	galens/lensnsga2evolver.cpp
	galens/openclsingleplanedeflection.cpp
	galens/lensgaparametricsingleplanecalculator.cpp
	prob/randomnumbergenerator.cpp
	fitness/numericgradientcalculator.cpp
	fitness/lensfitnessobject.cpp
	fitness/lensfitnessgeneral.cpp
	fitness/fitnesscomponent.cpp
	fitness/fitnessutil_overlap.cpp
	fitness/fitnessutil_null.cpp
	fitness/fitnessutil_weak.cpp
	fitness/fitnessutil_time.cpp
	fitness/fitnessutil_caustic.cpp
	fitness/pointgroupstorage.cpp
	fitness/fitnesscomponent_overlap.cpp
	fitness/fitnesscomponent_null.cpp
	fitness/fitnesscomponent_weak.cpp
	fitness/fitnesscomponent_time.cpp
	fitness/fitnesscomponent_dens.cpp
	fitness/fitnesscomponent_caustic.cpp
	fitness/fitnesscomponent_defl.cpp	
	fitness/fitnesscomponent_parity.cpp
	fitness/lensfitnessobject_defaults.cpp
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


