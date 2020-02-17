#foreach(i "experimentaltdcheck" "experimentaltdcheck2" "reltdcheck")
#	add_executable(${i} ${i}.cpp ../src/fitnessutil.cpp ../src/fitnesscomponent.cpp ../src/pointgroupstorage.cpp )
#	target_link_libraries(${i} ${GRALE_LIBRARIES})
#endforeach(i)

add_executable(genometest genometest.cpp ../src/fitnessutil.cpp ../src/fitnesscomponent.cpp ../src/pointgroupstorage.cpp 
	                  ../src/lensfitnesssimplerectangles.cpp )
target_link_libraries(genometest ${GRALE_LIBRARIES})

