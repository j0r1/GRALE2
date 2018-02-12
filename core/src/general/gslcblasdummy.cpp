#include <gsl/gsl_cblas.h>

namespace grale
{

// A call to a gslcblas function to avoid -dead_strip_dylibs removing
// the reference to this library
double _dummyGSLCBLASCall()
{
	const int N = 16;
	double X[N], Y[N];
	return cblas_ddot(N, X, 1, Y, 1);
}

} // end namespace
