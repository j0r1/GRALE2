set(GRALE_EXTERNAL_INCLUDES "")
set(GRALE_INTERNAL_INCLUDES "")
set(GRALE_LINK_LIBS "")

add_additional_stuff(INV GRALE_EXTERNAL_INCLUDES GRALE_LINK_LIBS)

find_package(ErrUt REQUIRED)
list(APPEND GRALE_EXTERNAL_INCLUDES "${ERRUT_INCLUDE_DIRS}")
list(APPEND GRALE_LINK_LIBS "${ERRUT_LIBRARIES}")

find_package(SerUt REQUIRED)
list(APPEND GRALE_EXTERNAL_INCLUDES "${SERUT_INCLUDE_DIRS}")
list(APPEND GRALE_LINK_LIBS "${SERUT_LIBRARIES}")

find_package(GSL REQUIRED)
list(APPEND GRALE_EXTERNAL_INCLUDES "${GSL_INCLUDE_DIR}")
list(APPEND GRALE_LINK_LIBS ${GSL_LIBRARIES})
list(APPEND GRALE_LINK_LIBS ${GSLCBLAS_LIBRARIES})

find_package(MOGAL REQUIRED)
list(APPEND GRALE_EXTERNAL_INCLUDES "${MOGAL_INCLUDE_DIRS}")
list(APPEND GRALE_LINK_LIBS "${MOGAL_LIBRARIES}")

list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_BINARY_DIR}/core/src/general")
list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/core/src/general")
list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/core/src/cosmology")
list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/core/src/image")
list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/core/src/lens")
list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/core/src/math")
list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/core/src/source")
list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/core/src/galens")
list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/core/src/thirdparty/WildMagic5")
list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/core/src/")
list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/inversion/src/galens")
list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/inversion/src/prob")
list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/inversion/src/")
#list(APPEND GRALE_INTERNAL_INCLUDES "${PROJECT_SOURCE_DIR}/inversion_modules/src/")

apply_include_paths("${GRALE_INTERNAL_INCLUDES}")
apply_include_paths("${GRALE_EXTERNAL_INCLUDES}")

add_executable(backproject backproject.cpp)
target_link_libraries(backproject graleinv-static gralecore-static ${GRALE_LINK_LIBS})

add_executable(precalcbackproj precalcbackproj.cpp)
target_link_libraries(precalcbackproj graleinv-static gralecore-static ${GRALE_LINK_LIBS})

add_executable(mpcu mpcu.cpp)
target_link_libraries(mpcu gralecore-static ${GRALE_LINK_LIBS})

add_executable(cosmcheck cosmcheck.cpp)
target_link_libraries(cosmcheck gralecore-static ${GRALE_LINK_LIBS})

add_executable(mpcubp mpcubp.cpp)
target_link_libraries(mpcubp graleinv-static gralecore-static ${GRALE_LINK_LIBS})
