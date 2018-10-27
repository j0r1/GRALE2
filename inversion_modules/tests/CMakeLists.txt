add_executable(experimentaltdcheck experimentaltdcheck.cpp ../src/fitnessutil.cpp ../src/fitnesscomponent.cpp ../src/pointgroupstorage.cpp)
target_link_libraries(experimentaltdcheck ${GRALE_LIBRARIES})
