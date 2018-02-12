find_path(GSL_INCLUDE_DIR "gsl/gsl_rng.h")
set(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR})

find_library(GSL_LIBRARY gsl)
find_library(GSL_CBLAS_LIBRARY gslcblas)
if (NOT GSL_CBLAS_LIBRARY)
	find_library(GSL_CBLAS_LIBRARY cblas)
endif()

if (GSL_LIBRARY AND GSL_CBLAS_LIBRARY)
	set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GSL DEFAULT_MSG GSL_INCLUDE_DIRS GSL_LIBRARIES)
