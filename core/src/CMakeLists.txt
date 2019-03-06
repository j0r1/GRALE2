add_definitions(-DGRALE_COMPILING)

set(HEADERS
	${PROJECT_BINARY_DIR}/core/src/general/graleconfig.h 
	general/constants.h 
	general/configurationparameters.h
	general/log.h
	lens/gravitationallens.h 
	lens/multipleplummerlens.h 
	lens/gausslens.h 
	lens/nsislens.h 
	lens/nsielens.h 
	lens/plummerlens.h 
	lens/pointmasslens.h 
	lens/sielens.h 
	lens/sislens.h 
	lens/symmetriclens.h 
	lens/real2dlensfunction.h 
	lens/squarelens.h 
	lens/multiplesquarelens.h 
	lens/multiplegausslens.h 
	lens/squarelensinfo.h 
	lens/gausslensinfo.h 
	lens/masssheetlens.h 
	lens/compositelens.h 
	lens/massdisklens.h 
	lens/profilelens.h 
	lens/polynomialmassprofilelens.h 
	lens/plummerlensinfo.h 
	lens/multiplewendlandlens.h 
	lens/deflectiongridlens.h 
	lens/circularlensprofile.h 
	lens/ellipticlens.h 
	lens/nfwlens.h 
	lens/ellipticnfwlens.h 
	lens/sersiclens.h 
	lens/ellipticsersiclens.h 
	lens/piemdlens.h
	lens/pimdlens.h
	lens/alphapotlens.h
	math/mathfunctions.h
	math/real1dfunction.h
	math/real1dfunctionintegrator.h
	math/real2dderivablefunction.h
	math/real2dfunction.h
	math/vector2d.h
	math/gridfunction.h
	math/triangle2d.h
	math/triangleindices.h
	math/line2d.h
	math/polygon2d.h
	math/rectangle2d.h
	math/contourfinder.h
	image/imagesdata.h
	image/imagesdataextended.h
	image/lensplane.h
	image/imageplane.h
	source/sourceimage.h
	source/circularsource.h
	source/ellipticalsource.h
	source/pointsource.h
	source/polygonsource.h
	source/discretesource.h
	galens/grid.h
	galens/gridlensinversionparameters.h
	galens/gaparameters.h
	thirdparty/WildMagic5/Wm5Assert.h
	thirdparty/WildMagic5/Wm5Box2.h
	thirdparty/WildMagic5/Wm5Box2.inl
	thirdparty/WildMagic5/Wm5ContMinBox2.h
	thirdparty/WildMagic5/Wm5ConvexHull.h
	thirdparty/WildMagic5/Wm5ConvexHull1.h
	thirdparty/WildMagic5/Wm5ConvexHull2.h
	thirdparty/WildMagic5/Wm5CoreLIB.h
	thirdparty/WildMagic5/Wm5FileIO.h
	thirdparty/WildMagic5/Wm5FileIO.inl
	thirdparty/WildMagic5/Wm5Math.h
	thirdparty/WildMagic5/Wm5Math.inl
	thirdparty/WildMagic5/Wm5MathematicsLIB.h
	thirdparty/WildMagic5/Wm5Query.h
	thirdparty/WildMagic5/Wm5Query2.h
	thirdparty/WildMagic5/Wm5Query2.inl
	thirdparty/WildMagic5/Wm5Tuple.h
	thirdparty/WildMagic5/Wm5Tuple.inl
	thirdparty/WildMagic5/Wm5Vector2.h
	thirdparty/WildMagic5/Wm5Vector2.inl
	thirdparty/WildMagic5/Wm5Intersector.h
	thirdparty/WildMagic5/Wm5IntrTriangle2Triangle2.h
	thirdparty/WildMagic5/Wm5Triangle2.h
	thirdparty/WildMagic5/Wm5Vector3.h
	thirdparty/WildMagic5/Wm5Vector3.inl
	thirdparty/WildMagic5/Wm5Intersector1.h
	thirdparty/WildMagic5/Wm5Triangle2.inl
	)

set(SOURCES
	general/debugnew.cpp 
	general/configurationparameters.cpp
	general/gslcblasdummy.cpp
	general/log.cpp
	math/real1dfunctionintegrator.cpp
	math/gridfunction.cpp
	math/contourfinder.cpp
	lens/gravitationallens.cpp 
	lens/gausslens.cpp 
	lens/multipleplummerlens.cpp 
	lens/nsielens.cpp 
	lens/nsislens.cpp 
	lens/plummerlens.cpp 
	lens/pointmasslens.cpp 
	lens/sielens.cpp 
	lens/sislens.cpp 
	lens/symmetriclens.cpp 
	lens/squarelens.cpp 
	lens/multiplesquarelens.cpp 
	lens/multiplegausslens.cpp 
	lens/masssheetlens.cpp 
	lens/compositelens.cpp 
	lens/massdisklens.cpp 
	lens/profilelens.cpp 
	lens/polynomialmassprofilelens.cpp 
	lens/multiplewendlandlens.cpp 
	lens/deflectiongridlens.cpp 
	lens/ellipticlens.cpp 
	lens/nfwlens.cpp 
	lens/ellipticnfwlens.cpp 
	lens/sersiclens.cpp 
	lens/ellipticsersiclens.cpp 
	lens/piemdlens.cpp
	lens/pimdlens.cpp
	lens/alphapotlens.cpp
	image/imagesdata.cpp
	image/imagesdataextended.cpp
	image/lensplane.cpp
	image/imageplane.cpp
	source/sourceimage.cpp
	source/circularsource.cpp
	source/ellipticalsource.cpp
	source/pointsource.cpp
	source/polygonsource.cpp
	source/discretesource.cpp
	galens/grid.cpp
	galens/gridlensinversionparameters.cpp
	galens/gaparameters.cpp
	thirdparty/WildMagic5/Wm5Assert.cpp
	thirdparty/WildMagic5/Wm5ConvexHull.cpp
	thirdparty/WildMagic5/Wm5ConvexHull2.cpp
	thirdparty/WildMagic5/Wm5Intersector1.cpp
	thirdparty/WildMagic5/Wm5IntrTriangle2Triangle2.cpp
	thirdparty/WildMagic5/Wm5Vector3.cpp
	thirdparty/WildMagic5/Wm5Intersector.cpp
	thirdparty/WildMagic5/Wm5ContMinBox2.cpp
	thirdparty/WildMagic5/Wm5ConvexHull1.cpp
	thirdparty/WildMagic5/Wm5Delaunay.cpp
	thirdparty/WildMagic5/Wm5Delaunay1.cpp
	thirdparty/WildMagic5/Wm5Delaunay2.cpp
	thirdparty/WildMagic5/Wm5Endian.cpp
	thirdparty/WildMagic5/Wm5FileIO.cpp
	thirdparty/WildMagic5/Wm5Math.cpp
	thirdparty/WildMagic5/Wm5Query.cpp
	thirdparty/WildMagic5/Wm5VEManifoldMesh.cpp
	thirdparty/WildMagic5/Wm5Vector2.cpp
	#thirdparty/ttl/src/halfedge/HeTriang.cpp
	#thirdparty/ttl/src/utils/HandleId.cpp
	)

if (MSVC)
	set(CMAKE_DEBUG_POSTFIX _d)
endif (MSVC)

if (NOT MSVC OR GRALE_COMPILE_STATIC)
	add_library(gralecore-static STATIC ${SOURCES} ${HEADERS})
	set_target_properties(gralecore-static PROPERTIES OUTPUT_NAME grale2core)
	set_target_properties(gralecore-static PROPERTIES CLEAN_DIRECT_OUTPUT 1)
	set(GRALE_INSTALLTARGETS gralecore-static)
	target_link_libraries(gralecore-static ${GRALE_LINK_LIBS})
endif (NOT MSVC OR GRALE_COMPILE_STATIC)

if (NOT MSVC OR NOT GRALE_COMPILE_STATIC)
	add_library(gralecore-shared SHARED ${SOURCES} ${HEADERS})
	set_target_properties(gralecore-shared PROPERTIES OUTPUT_NAME grale2core)
	set_target_properties(gralecore-shared PROPERTIES VERSION ${VERSION})
	set_target_properties(gralecore-shared PROPERTIES CLEAN_DIRECT_OUTPUT 1)
	set(GRALE_INSTALLTARGETS ${GRALE_INSTALLTARGETS} gralecore-shared)
	if (GRALE_LINK)
		target_link_libraries(gralecore-shared ${GRALE_LINK_LIBS})
	endif (GRALE_LINK)
endif (NOT MSVC OR NOT GRALE_COMPILE_STATIC)

apply_include_paths("${GRALE_INTERNAL_INCLUDES}")
apply_include_paths("${GRALE_EXTERNAL_INCLUDES}")

install(FILES ${HEADERS} DESTINATION include/grale)
install(TARGETS ${GRALE_INSTALLTARGETS} DESTINATION ${LIBRARY_INSTALL_DIR})


