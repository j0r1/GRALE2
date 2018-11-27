add_executable(experimentaltdcheck experimentaltdcheck.cpp ../src/fitnessutil.cpp ../src/fitnesscomponent.cpp ../src/pointgroupstorage.cpp)
add_executable(reltdcheck reltdcheck.cpp ../src/fitnessutil.cpp ../src/fitnesscomponent.cpp ../src/pointgroupstorage.cpp)
target_link_libraries(experimentaltdcheck ${GRALE_LIBRARIES})
target_link_libraries(reltdcheck ${GRALE_LIBRARIES})
