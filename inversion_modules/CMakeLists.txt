cmake_minimum_required(VERSION 2.8)

project(gamods)

include(${PROJECT_SOURCE_DIR}/cmake/Macros.cmake)
get_install_directory(LIBRARY_INSTALL_DIR)

include(CheckCXXCompilerFlag)
include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckCXXSourceCompiles)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
	# Test if the compiler already supports C++11
	check_cxx_source_compiles("#include <map>\nusing namespace std;\nint main(void) \n {\nmap<int, double> m;\n auto it = m.begin();\n return 0; \n }" C11AVAIL)
	if (NOT C11AVAIL)
		message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
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

add_library(galens_pointoverlap SHARED src/galensmodule1/pointoverlap.cpp src/lensfitnesspointoverlap.cpp src/fitnessutil.cpp 
	                               src/fitnesscomponent.cpp src/pointgroupstorage.cpp)
target_link_libraries(galens_pointoverlap ${EXTRA_LIBS} ${GRALE_LIBRARIES})
install(TARGETS galens_pointoverlap DESTINATION ${LIBRARY_INSTALL_DIR})

add_library(galens_pointoverlapnull SHARED src/galensmodule1/pointoverlapnull.cpp src/lensfitnesspointoverlap.cpp src/fitnessutil.cpp 
	                                   src/fitnesscomponent.cpp src/pointgroupstorage.cpp)
target_link_libraries(galens_pointoverlapnull ${EXTRA_LIBS} ${GRALE_LIBRARIES})
install(TARGETS galens_pointoverlapnull DESTINATION ${LIBRARY_INSTALL_DIR})

add_library(galens_simplerect SHARED src/galensmodule1/galens1.cpp src/lensfitnesssimplerectangles.cpp
	                             src/pointgroupstorage.cpp src/fitnessutil.cpp src/fitnesscomponent.cpp)
target_link_libraries(galens_simplerect ${EXTRA_LIBS} ${GRALE_LIBRARIES})
install(TARGETS galens_simplerect DESTINATION ${LIBRARY_INSTALL_DIR})

add_library(galens_j1004 SHARED src/j10044112/j1004test.cpp src/lensfitnessj1004.cpp src/pointgroupstorage.cpp src/fitnesscomponent.cpp src/fitnessutil.cpp)
target_link_libraries(galens_j1004 ${EXTRA_LIBS} ${GRALE_LIBRARIES})
install(TARGETS galens_j1004 DESTINATION ${LIBRARY_INSTALL_DIR})

add_library(galens_weak SHARED src/weak/weak.cpp src/lensfitnessweak.cpp src/fitnessutil.cpp src/fitnesscomponent.cpp src/pointgroupstorage.cpp)
target_link_libraries(galens_weak ${EXTRA_LIBS} ${GRALE_LIBRARIES})
install(TARGETS galens_weak DESTINATION ${LIBRARY_INSTALL_DIR})

add_library(galens_pointoverlapnullweak SHARED src/weak/pointoverlapnullweak.cpp src/lensfitnesspointoverlapweak.cpp 
	                                       src/fitnessutil.cpp src/fitnesscomponent.cpp src/pointgroupstorage.cpp)
target_link_libraries(galens_pointoverlapnullweak ${EXTRA_LIBS} ${GRALE_LIBRARIES})
install(TARGETS galens_pointoverlapnullweak DESTINATION ${LIBRARY_INSTALL_DIR})

add_library(galens_extendednullweak SHARED src/weak/extendednullweak.cpp src/lensfitnessextendedweak.cpp src/pointgroupstorage.cpp src/fitnessutil.cpp src/fitnesscomponent.cpp)
target_link_libraries(galens_extendednullweak ${EXTRA_LIBS} ${GRALE_LIBRARIES})
install(TARGETS galens_extendednullweak DESTINATION ${LIBRARY_INSTALL_DIR})

add_library(galens_general SHARED src/general/general.cpp src/lensfitnessgeneral.cpp src/fitnesscomponent.cpp src/fitnessutil.cpp
	                          src/pointgroupstorage.cpp)
target_link_libraries(galens_general ${EXTRA_LIBS} ${GRALE_LIBRARIES})
install(TARGETS galens_general DESTINATION ${LIBRARY_INSTALL_DIR})
