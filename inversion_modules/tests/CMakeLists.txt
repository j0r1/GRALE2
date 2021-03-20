#foreach(i "experimentaltdcheck" "experimentaltdcheck2" "reltdcheck")
#	add_executable(${i} ${i}.cpp ../src/fitnessutil.cpp ../src/fitnesscomponent.cpp ../src/pointgroupstorage.cpp )
#	target_link_libraries(${i} ${GRALE_LIBRARIES})
#endforeach(i)

include_directories(../src/old)
add_executable(genometest genometest.cpp ../src/old/lensfitnesssimplerectangles.cpp)
target_link_libraries(genometest ${GRALE_LIBRARIES})

add_executable(bayeswlcheck bayeswlcheck.cpp)
target_link_libraries(bayeswlcheck ${GRALE_LIBRARIES})

