/*

  This file is a part of GRALE, a library to facilitate the simulation
  and inversion of gravitational lenses.

  Copyright (C) 2008-2012 Jori Liesenborgs

  Contact: jori.liesenborgs@gmail.com
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
  
*/

#ifndef GRALE_MATHFUNCTIONS_H

#define GRALE_MATHFUNCTIONS_H

#include <string.h>
#include <cmath>
#include <algorithm>

namespace grale
{
	template<class T>
	inline T CalculateDotProduct(const T *a, const T *b, int num)
	{
		T v = 0;
		for (int i = 0 ; i < num ; i++)
			v += a[i]*b[i];
		return v;
	}

	template<class T>
	inline void CalculateAddProductC(T *pDest, const T *pSrc, T C, int num, T *tmpbuf)
	{
		for (int i = 0 ; i < num ; i++)
			pDest[i] = pDest[i] + pSrc[i]*C;
	}

	template<class T>
	inline void CopyVector(T *pDest, const T *pSrc, int num)
	{
		memcpy(pDest, pSrc, num*sizeof(T));
	}

	template<class T>
	inline void SquareVector(T *pVect, int num)
	{
		for (int i = 0 ; i < num ; i++)
			pVect[i] *= pVect[i];
	}

	template<class T>
	inline void SquareVector(T *pDest, const T *pSrc, int num)
	{
		for (int i = 0 ; i < num ; i++)
			pDest[i] = pSrc[i]*pSrc[i];
	}
	template<class T>
	inline void MultiplyVector(T *pSrcDst, const T *pSrc, int num)
	{
		for (int i = 0 ; i < num ; i++)
			pSrcDst[i] *= pSrc[i];
	}

	template<class T>
	inline void MultiplyVector(T *pDst, const T *pSrc, const T factor, int num)
	{
		for (int i = 0 ; i < num ; i++)
			pDst[i] = pSrc[i]*factor;
	}

	template<class T>
	inline void MultiplyVector(T *pSrcDst, const T factor, int num)
	{
		for (int i = 0 ; i < num ; i++)
			pSrcDst[i] *= factor;
	}
	
	template<class T>
	inline void MultiplyVector(T *pDst, const T *pSrc1, const T *pSrc2, int num)
	{
		for (int i = 0 ; i < num ; i++)
			pDst[i] = pSrc1[i]*pSrc2[i];
	}

	template<class T>
	inline void AbsVector(T *pVect, int num)
	{
		for (int i = 0 ; i < num ; i++)
			pVect[i] = std::abs(pVect[i]);
	}

	template<class T>
	inline void AddVector(T *pSrcDst, const T *pSrc, int num)
	{
		for (int i = 0 ; i < num ; i++)
			pSrcDst[i] += pSrc[i];
	}

	template<class T>
	inline void AddVector(T *pSrcDst, T value, int num)
	{
		for (int i = 0 ; i < num ; i++)
			pSrcDst[i] += value;
	}

	template<class T>
	inline void AddVector(T *pDst, const T *pSrc1, const T *pSrc2, int num)
	{
		for (int i = 0 ; i < num ; i++)
			pDst[i] = pSrc1[i] + pSrc2[i];
	}

	template<class T>
	inline void SubVector(T *pSrcDst, const T *pSrc, int num)
	{
		for (int i = 0 ; i < num ; i++)
			pSrcDst[i] -= pSrc[i];
	}

	template<class T>
	inline void SubVector(T *pDst, const T *pSrc1, const T *pSrc2, int num)
	{
		for (int i = 0 ; i < num ; i++)
			pDst[i] = pSrc1[i] - pSrc2[i];
	}
}

template<class T> inline T ABS(T x) { return std::abs(x); }
#define SQRT(x)		std::sqrt(x)
#define COS(x)		std::cos(x)
#define SIN(x)		std::sin(x)
#define TAN(x)		std::tan(x)
#define ACOS(x)		std::acos(x)
#define ASIN(x)		std::asin(x)
#define ATAN(x)		std::atan(x)
#define LN(x)		std::log(x)
#define ASINH(x)	std::asinh(x)
#define ACOSH(x)	std::acosh(x)
#define ATANH(x)	std::atanh(x)
#define POW(x,y)	std::pow(x,y)
#define EXP(x)		std::exp(x)
#define SINH(x)		std::sinh(x)
#define COSH(x)		std::cosh(x)
#define TANH(x)		std::tanh(x)
#define MIN(x,y)	std::min(x,y)
#define MAX(x,y)	std::max(x,y)

#endif // GRALE_MATHFUNCTIONS_H

