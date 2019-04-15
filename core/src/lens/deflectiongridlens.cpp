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

#include "graleconfig.h"
#include "deflectiongridlens.h"
#include "gridfunction.h"
#include "constants.h"

#include "debugnew.h"

namespace grale
{

// Lens parameters

DeflectionGridLensParams::DeflectionGridLensParams()
{
}

DeflectionGridLensParams::DeflectionGridLensParams(const std::vector<double> &alphaX, const std::vector<double> &alphaY, 
			         int width, int height, Vector2D<double> bottomLeft, Vector2D<double> topRight)
{
	m_alphaX = alphaX;
	m_alphaY = alphaY;
	m_width = width;
	m_height = height;
	m_bottomLeft = bottomLeft;
	m_topRight = topRight;
}

GravitationalLensParams *DeflectionGridLensParams::createCopy() const
{
	int totalSize = m_width*m_height;
	
	if (!(m_alphaX.size() == totalSize && m_alphaY.size() == totalSize))
	{
		setErrorString("Data length doesn't match the specified dimensions");
		return 0;
	}

	return new DeflectionGridLensParams(m_alphaX, m_alphaY, m_width, m_height, m_bottomLeft, m_topRight);
}

bool DeflectionGridLensParams::write(serut::SerializationInterface &si) const
{
	// Perform a check on the validity of the data
	
	int totalSize = m_width*m_height;
	
	if (!(m_alphaX.size() == totalSize && m_alphaY.size() == totalSize))
	{
		setErrorString("Data length doesn't match the specified dimensions");
		return false;
	}

	int32_t dimensions[2];

	dimensions[0] = m_width;
	dimensions[1] = m_height;

	if (!si.writeInt32s(dimensions, 2))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeDoubles(m_alphaX))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeDoubles(m_alphaY))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	double region[4];
	
	region[0] = m_bottomLeft.getX();
	region[1] = m_bottomLeft.getY();
	region[2] = m_topRight.getX();
	region[3] = m_topRight.getY();

	if (!si.writeDoubles(region, 4))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	return true;
}

bool DeflectionGridLensParams::read(serut::SerializationInterface &si)
{
	int32_t dimensions[2];

	if (!si.readInt32s(dimensions, 2))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	int totalSize = dimensions[0]*dimensions[1];

	std::vector<double> alphaX(totalSize);
	std::vector<double> alphaY(totalSize);

	if (!si.readDoubles(alphaX))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.readDoubles(alphaY))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	double region[4];

	if (!si.readDoubles(region, 4))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	m_width = dimensions[0];
	m_height = dimensions[1];
	m_bottomLeft = Vector2D<double>(region[0], region[1]);
	m_topRight = Vector2D<double>(region[2], region[3]);
	m_alphaX = alphaX;
	m_alphaY = alphaY;

	return true;
}

// Lens itself

DeflectionGridLens::DeflectionGridLens() : GravitationalLens(DeflectionGrid)
{
	m_pAxFunction = 0;
	m_pAyFunction = 0;
}

DeflectionGridLens::~DeflectionGridLens()
{
	if (m_pAxFunction)
		delete m_pAxFunction;
	if (m_pAyFunction)
		delete m_pAyFunction;
}

bool DeflectionGridLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const DeflectionGridLensParams *pParams = dynamic_cast<const DeflectionGridLensParams *>(pLensParams);

	if (pParams == 0)
	{
		setErrorString("Parameters are not of type 'DeflectionGridLensParams'");
		return false;
	}

	Vector2D<double> bottomLeft = pParams->getBottomLeft();
	Vector2D<double> topRight = pParams->getTopRight();
	m_pixelWidth = ABS((topRight.getX() - bottomLeft.getX())/(double)(pParams->getWidth()-1));
	m_pixelHeight = ABS((topRight.getY() - bottomLeft.getY())/(double)(pParams->getHeight()-1));

	m_alphaX = pParams->getAlphaX();
	m_alphaY = pParams->getAlphaY();

	m_pAxFunction = new GridFunction(&(m_alphaX[0]), bottomLeft, topRight, pParams->getWidth(), 
			                 pParams->getHeight());
	m_pAyFunction = new GridFunction(&(m_alphaY[0]), bottomLeft, topRight, pParams->getWidth(), 
			                 pParams->getHeight());
	m_densFactor = (SPEED_C*SPEED_C)/(getLensDistance()*8.0*CONST_PI*CONST_G);

	m_x0 = bottomLeft.getX();
	m_x1 = topRight.getX();
	m_y0 = bottomLeft.getY();
	m_y1 = topRight.getY();

	if (m_x0 > m_x1) std::swap(m_x0, m_x1);
	if (m_y0 > m_y1) std::swap(m_y0, m_y1);

	return true;
}

bool DeflectionGridLens::getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const
{
	if (theta.getX() < m_x0 || theta.getX() > m_x1 ||
		theta.getY() < m_y0 || theta.getY() > m_y1)
	{
		setErrorString("Theta position doesn't lie inside area of deflection grid");
		return false;
	}

	double alphaX = (*m_pAxFunction)(theta);
	double alphaY = (*m_pAyFunction)(theta);

	*pAlpha = Vector2D<double>(alphaX, alphaY);
	return true;
}

double DeflectionGridLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	if (theta.getX() < m_x0 || theta.getX() > m_x1 ||
		theta.getY() < m_y0 || theta.getY() > m_y1)
		return std::numeric_limits<double>::quiet_NaN();

	Vector2D<double> x2 = theta + Vector2D<double>(m_pixelWidth/2.0, 0);
	Vector2D<double> x1 = theta - Vector2D<double>(m_pixelWidth/2.0, 0);
	Vector2D<double> y2 = theta + Vector2D<double>(0, m_pixelHeight/2.0);
	Vector2D<double> y1 = theta - Vector2D<double>(0, m_pixelHeight/2.0);

	double daxx = (*m_pAxFunction)(x2) - (*m_pAxFunction)(x1);
	double dayy = (*m_pAyFunction)(y2) - (*m_pAyFunction)(y1);

	return m_densFactor*(daxx/m_pixelWidth + dayy/m_pixelHeight);
}

bool DeflectionGridLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	if (theta.getX() < m_x0 || theta.getX() > m_x1 ||
		theta.getY() < m_y0 || theta.getY() > m_y1)
	{
		setErrorString("Theta position doesn't lie inside area of deflection grid");
		return false;
	}

	Vector2D<double> x2 = theta + Vector2D<double>(m_pixelWidth/2.0, 0);
	Vector2D<double> x1 = theta - Vector2D<double>(m_pixelWidth/2.0, 0);
	Vector2D<double> y2 = theta + Vector2D<double>(0, m_pixelHeight/2.0);
	Vector2D<double> y1 = theta - Vector2D<double>(0, m_pixelHeight/2.0);

	double daxx = (*m_pAxFunction)(x2) - (*m_pAxFunction)(x1);
	double dayy = (*m_pAyFunction)(y2) - (*m_pAyFunction)(y1);
	double daxy = (*m_pAxFunction)(y2) - (*m_pAxFunction)(y1);
	double dayx = (*m_pAyFunction)(x2) - (*m_pAyFunction)(x1);

	axx = daxx/m_pixelWidth;
	ayy = dayy/m_pixelHeight;
	axy = 0.5*(daxy/m_pixelHeight + dayx/m_pixelWidth);

	return true;
}

} // end namespace

