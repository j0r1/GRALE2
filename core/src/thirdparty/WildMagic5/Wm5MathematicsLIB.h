// Geometric Tools, LLC
// Copyright (c) 1998-2011
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
//
// File Version: 5.0.4 (2011/07/09)

#ifndef WM5MATHEMATICSLIB_H
#define WM5MATHEMATICSLIB_H

#include "graleconfig.h"
#include "Wm5CoreLIB.h"

#define WM5_MATHEMATICS_ITEM GRALE_IMPORTEXPORT

// Enable this define if you want the Rational class to assert when the
// constructor is passed a floating-point infinity or NaN.
//#define WM5_ASSERT_ON_RATIONAL_CONVERT_NAN

// Enable this define if you want Vector2<Real>::GetBarycentrics to assert
// when the input triangle is degenerate.
//#define WM5_ASSERT_ON_BARYCENTRIC2_DEGENERATE

// Enable this define if you want Vector3<Real>::GetBarycentrics to assert
// when the input tetrahedron is degenerate.
//#define WM5_ASSERT_ON_BARYCENTRIC3_DEGENERATE

// Enable this define if you want index range checking in GMatrix<Real>.
#define WM5_ASSERT_GMATRIX_OUT_OF_RANGE

#endif
