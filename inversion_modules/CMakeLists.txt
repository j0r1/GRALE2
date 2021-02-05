cmake_minimum_required(VERSION 3.1)
if (POLICY CMP0042)
	cmake_policy(SET CMP0042 NEW)
endif()

project(gamods)

include(${PROJECT_SOURCE_DIR}/cmake/Macros.cmake)
get_install_directory(LIBRARY_INSTALL_DIR)

include(CheckCXXCompilerFlag)
include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckCXXSourceCompiles)

CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
if(COMPILER_SUPPORTS_CXX14)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else()
	# Test if the compiler already supports C++14
	check_cxx_source_compiles("int main(void)\n{\n\tauto f = [](auto x)\n\t{\n\t\treturn x*x;\n\t};\n\n\treturn 0;\n}\n" C14AVAIL)
	if (NOT C14AVAIL)
		message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support.")
	endif()
endif()

find_package(GRALECORE REQUIRED)
find_package(GRALEINV REQUIRED)
set(GRALE_LIBRARIES "${GRALECORE_LIBRARIES}" "${GRALEINV_LIBRARIES}")

add_additional_stuff(EXTRA_INCLUDES EXTRA_LIBS)

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${GRALECORE_INCLUDE_DIRS})
include_directories(${GRALEINV_INCLUDE_DIRS})

if (NOT UNIX)
	add_definitions("-DGAMODS_EXPORT=__declspec(dllexport)")
else()
	add_definitions("-DGAMODS_EXPORT=")
endif()

#add_library(galens_pointoverlap SHARED src/galensmodule1/pointoverlap.cpp src/lensfitnesspointoverlap.cpp src/fitnessutil.cpp 
# 	                               src/fitnesscomponent.cpp src/pointgroupstorage.cpp)
#target_link_libraries(galens_pointoverlap ${EXTRA_LIBS} ${GRALE_LIBRARIES})
#install(TARGETS galens_pointoverlap DESTINATION ${LIBRARY_INSTALL_DIR})

#add_library(galens_pointoverlapnull SHARED src/galensmodule1/pointoverlapnull.cpp src/lensfitnesspointoverlap.cpp src/fitnessutil.cpp 
#	                                   src/fitnesscomponent.cpp src/pointgroupstorage.cpp)
#target_link_libraries(galens_pointoverlapnull ${EXTRA_LIBS} ${GRALE_LIBRARIES})
#install(TARGETS galens_pointoverlapnull DESTINATION ${LIBRARY_INSTALL_DIR})

#add_library(galens_simplerect SHARED src/galensmodule1/galens1.cpp src/lensfitnesssimplerectangles.cpp
#	                             src/pointgroupstorage.cpp src/fitnessutil.cpp src/fitnesscomponent.cpp)
#target_link_libraries(galens_simplerect ${EXTRA_LIBS} ${GRALE_LIBRARIES})
#install(TARGETS galens_simplerect DESTINATION ${LIBRARY_INSTALL_DIR})

#add_library(galens_j1004 SHARED src/j10044112/j1004test.cpp src/lensfitnessj1004.cpp src/pointgroupstorage.cpp src/fitnesscomponent.cpp src/fitnessutil.cpp)
#target_link_libraries(galens_j1004 ${EXTRA_LIBS} ${GRALE_LIBRARIES})
#install(TARGETS galens_j1004 DESTINATION ${LIBRARY_INSTALL_DIR})

add_library(galens_general SHARED src/general/general.cpp src/lensfitnessgeneral.cpp src/fitnesscomponent.cpp
				  src/fitnessutil_overlap.cpp src/fitnessutil_null.cpp src/fitnessutil_weak.cpp
				  src/fitnessutil_time.cpp src/fitnessutil_caustic.cpp
	                          src/pointgroupstorage.cpp src/numericgradientcalculator.cpp
				  src/fitnesscomponent_overlap.cpp src/fitnesscomponent_null.cpp src/fitnesscomponent_weak.cpp
				  src/fitnesscomponent_time.cpp src/fitnesscomponent_dens.cpp src/fitnesscomponent_caustic.cpp
				  src/fitnesscomponent_defl.cpp
				  )
target_link_libraries(galens_general ${EXTRA_LIBS} ${GRALE_LIBRARIES})
install(TARGETS galens_general DESTINATION ${LIBRARY_INSTALL_DIR})

add_library(galens_general_gpu SHARED src/general/general_gpu.cpp src/lensfitnessgeneral.cpp src/fitnesscomponent.cpp
				  src/fitnessutil_overlap.cpp src/fitnessutil_null.cpp src/fitnessutil_weak.cpp
				  src/fitnessutil_time.cpp src/fitnessutil_caustic.cpp
	                          src/pointgroupstorage.cpp src/numericgradientcalculator.cpp
				  src/fitnesscomponent_overlap.cpp src/fitnesscomponent_null.cpp src/fitnesscomponent_weak.cpp
				  src/fitnesscomponent_time.cpp src/fitnesscomponent_dens.cpp src/fitnesscomponent_caustic.cpp
				  src/fitnesscomponent_defl.cpp
				  )
target_link_libraries(galens_general_gpu ${EXTRA_LIBS} ${GRALE_LIBRARIES})
install(TARGETS galens_general_gpu DESTINATION ${LIBRARY_INSTALL_DIR})

add_subdirectory(tests)
