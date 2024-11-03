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
#include <cmath>
#include <iostream>
#include <limits>

using namespace std;

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

std::unique_ptr<GravitationalLensParams> DeflectionGridLensParams::createCopy() const
{
	int totalSize = m_width*m_height;
	
	if (!((int)m_alphaX.size() == totalSize && (int)m_alphaY.size() == totalSize))
	{
		setErrorString("Data length doesn't match the specified dimensions");
		return 0;
	}

	return std::make_unique<DeflectionGridLensParams>(m_alphaX, m_alphaY, m_width, m_height, m_bottomLeft, m_topRight);
}

bool DeflectionGridLensParams::write(serut::SerializationInterface &si) const
{
	// Perform a check on the validity of the data
	
	int totalSize = m_width*m_height;
	
	if (!((int)m_alphaX.size() == totalSize && (int)m_alphaY.size() == totalSize))
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
}

DeflectionGridLens::~DeflectionGridLens()
{
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

	m_pAxFunction = make_unique<GridFunction>(&(m_alphaX[0]), bottomLeft, topRight, pParams->getWidth(), 
			                 pParams->getHeight());
	m_pAyFunction = make_unique<GridFunction>(&(m_alphaY[0]), bottomLeft, topRight, pParams->getWidth(), 
			                 pParams->getHeight());
	m_densFactor = (SPEED_C*SPEED_C)/(getLensDistance()*8.0*CONST_PI*CONST_G);

	m_x0 = bottomLeft.getX();
	m_x1 = topRight.getX();
	m_y0 = bottomLeft.getY();
	m_y1 = topRight.getY();

	if (m_x0 > m_x1) std::swap(m_x0, m_x1);
	if (m_y0 > m_y1) std::swap(m_y0, m_y1);

	// Build phi values
	// ax = (phi_next-phi_current)/dx => phi_next = ax*dx+phi_current

	int W = pParams->getWidth();
	int H = pParams->getHeight();
	double gridW = topRight.getX()-bottomLeft.getX();
	double gridH = topRight.getY()-bottomLeft.getY();
	double dx = gridW/(W-1);
	double dy = gridH/(H-1);

	{
		// First process each row separately
		m_phiFromX.resize(H*(W+1));
		for (int y = 0 ; y < H ; y++)
		{
			m_phiFromX[y*(W+1)+0] = 0;
			for (int x = 0 ; x < W ; x++)
				m_phiFromX[y*(W+1)+x+1] = m_phiFromX[y*(W+1)+x] + m_pixelWidth*m_alphaX[y*W+x];

			// Then try to introduce an offset in the rows so that the Y deflection matches
			if (y == 0)
				continue;

			double requiredDiff = 0;
			double requiredDiff2 = 0;
			for (int x = 1 ; x < W ; x++)
			{
				// ay = (phi_cur+reqDiff-phi_prev)/Dy => ay*Dy+phi_prev-phi_cur = reqDiff
				double ay = (m_alphaY[y*W+x]+m_alphaY[y*W+x-1])*0.5;
				double reqDiff = ay*m_pixelHeight + m_phiFromX[(y-1)*(W+1)+x] - m_phiFromX[y*(W+1)+x];
				requiredDiff += reqDiff;
				requiredDiff2 += reqDiff*reqDiff;
			}
			requiredDiff /= (double)(W-1);
			requiredDiff2 /= (double)(W-1);

			//double stdDev = SQRT(requiredDiff2-requiredDiff*requiredDiff);
			//cerr << requiredDiff << " " << requiredDiff2 << endl;

			for (int x = 1 ; x < W+1 ; x++)
				m_phiFromX[y*(W+1)+x] += requiredDiff;
		}

		Vector2Dd bottomLeftExtraX { bottomLeft.getX()-dx*0.5, bottomLeft.getY() };
		Vector2Dd topRightExtraX { topRight.getX()+dx*0.5, topRight.getY() };

		m_pPhiFromXFunction = make_unique<GridFunction>(&(m_phiFromX[0]), bottomLeftExtraX, topRightExtraX, W+1, H);
	}
	
	{
		// Do the same for the columns
		m_phiFromY.resize((H+1)*W);
		for (int x = 0 ; x < W ; x++)
		{
			m_phiFromY[0*W+x] = 0;
			for (int y = 0 ; y < H ; y++)
				m_phiFromY[(y+1)*W+x] = m_phiFromY[y*W+x] + m_pixelHeight*m_alphaY[y*W+x];

			// Then try to introduce an offset in the columns so that the X deflection matches
			if (x == 0)
				continue;

			double requiredDiff = 0;
			double requiredDiff2 = 0;
			for (int y = 1 ; y < H ; y++)
			{
				// ax = (phi_cur+reqDiff-phi_prev)/Dx => ax*Dx+phi_prev-phi_cur = reqDiff
				double ax = (m_alphaX[y*W+x]+m_alphaX[(y-1)*W+x])*0.5;
				double reqDiff = ax*m_pixelWidth + m_phiFromY[y*W+x-1] - m_phiFromY[y*W+x];
				requiredDiff += reqDiff;
				requiredDiff2 += reqDiff*reqDiff;
			}
			requiredDiff /= (double)(H-1);
			requiredDiff2 /= (double)(H-1);

			//double stdDev = SQRT(requiredDiff2-requiredDiff*requiredDiff);
			//cerr << requiredDiff << " " << requiredDiff2 << endl;

			for (int y = 1 ; y < H+1 ; y++)
				m_phiFromY[y*W+x] += requiredDiff;
		}

		Vector2Dd bottomLeftExtraY { bottomLeft.getX(), bottomLeft.getY()-dy*0.5 };
		Vector2Dd topRightExtraY { topRight.getX(), topRight.getY()+dy*0.5 };

		m_pPhiFromYFunction = make_unique<GridFunction>(&(m_phiFromY[0]), bottomLeftExtraY, topRightExtraY, W, H+1);
	}

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

bool DeflectionGridLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	if (theta.getX() < m_x0 || theta.getX() > m_x1 ||
		theta.getY() < m_y0 || theta.getY() > m_y1)
	{
		setErrorString("Theta position doesn't lie inside area of deflection grid");
		return false;
	}

	double phiFromX = (*m_pPhiFromXFunction)(theta);
	double phiFromY = (*m_pPhiFromYFunction)(theta);
	double phi = 0.5*(phiFromX+phiFromY);
	*pPotentialValue = (D_ds/D_s)*phi;
	return true;
}

bool DeflectionGridLens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	// TODO: look at grid contents, grid size?
	double angle = ANGLE_ARCSEC;

	*pDeflectionScale = angle;
	*pPotentialScale = angle*angle;

	return true;
}

bool DeflectionGridLens::getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const
{
	// TODO? Anything else?
	*pNumIntParams = 2; // numX, numY for deflection points
	*pNumFloatParams = 4 // corners 
	                 + m_alphaX.size()*2; // same number for alphaY
	return true;
}

bool DeflectionGridLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	const int w = m_pAxFunction->getWidth();
	const int h = m_pAxFunction->getHeight();
	
	pIntParams[0] = w;
	pIntParams[1] = h;

	pFloatParams[0] = (float)(m_pAxFunction->getBottomLeft().getX()/deflectionScale);
	pFloatParams[1] = (float)(m_pAxFunction->getBottomLeft().getY()/deflectionScale);
	pFloatParams[2] = (float)(m_pAxFunction->getTopRight().getX()/deflectionScale);
	pFloatParams[3] = (float)(m_pAxFunction->getTopRight().getY()/deflectionScale);
	
	pFloatParams += 4;

	for (int y = 0 ; y < h ; y++)
	{
		for (int x = 0 ; x < w ; x++)
		{
			const int pos = (x+y*w)*2;
			pFloatParams[pos + 0] = (float)(m_pAxFunction->getValue(x, y)/deflectionScale);
			pFloatParams[pos + 1] = (float)(m_pAyFunction->getValue(x, y)/deflectionScale);
		}
	}
	return true;
}

std::unique_ptr<GravitationalLensParams> DeflectionGridLens::createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const
{
	vector<double> alphaX, alphaY;
	const int w = m_pAxFunction->getWidth();
	const int h = m_pAxFunction->getHeight();

	double blX = ((double)pFloatParams[0])*deflectionScale;
	double blY = ((double)pFloatParams[1])*deflectionScale;
	double trX = ((double)pFloatParams[2])*deflectionScale;
	double trY = ((double)pFloatParams[3])*deflectionScale;

	pFloatParams += 4;

	for (int y = 0 ; y < h ; y++)
	{
		for (int x = 0 ; x < w ; x++)
		{
			const int pos = (x+y*w)*2;
			double ax = (double)pFloatParams[pos+0] * deflectionScale;
			double ay = (double)pFloatParams[pos+1] * deflectionScale;
			alphaX.push_back(ax);
			alphaY.push_back(ay);
		}
	}

	return make_unique<DeflectionGridLensParams>(alphaX, alphaY, w, h, Vector2Dd(blX, blY), Vector2Dd(trX, trY));
}

string DeflectionGridLens::getCLProgram(double deflectionScale, double potentialScale, string &subRoutineName, bool derivatives, bool potential) const
{
	subRoutineName = "clDeflectionGridLensProgram";
	string program = R"XYZ(

void clDeflectionGridLensProgram_getPosAndFrac(float x, float x0, float x1, int num, int *pPix, float *pFrac)
{
	float f = (x-x0)/(x1-x0)*(num-1);
	int fInt = (int)f;
	if (fInt > num-2)
		fInt = num-2;

	*pPix = fInt;
	*pFrac = f - (float)fInt;
}

LensQuantities clDeflectionGridLensProgram(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)
{
	const int numX = pIntParams[0];
	const int numY = pIntParams[1];
	const float left = pFloatParams[0];
	const float bottom = pFloatParams[1];
	const float right = pFloatParams[2];
	const float top = pFloatParams[3];
	
	pFloatParams += 4;

	LensQuantities r;

	const float x = coord.x;
	const float y = coord.y;
	if (x < left || x > right || y < bottom || y > top)
	{
		r.alphaX = nan((uint)0);
		r.alphaY = nan((uint)0);
	}
	else
	{
		int X, Y;
		float t, u;

		clDeflectionGridLensProgram_getPosAndFrac(x, left, right, numX, &X, &t);
		clDeflectionGridLensProgram_getPosAndFrac(y, bottom, top, numY, &Y, &u);
		//printf("x = %g, left = %g, right = %g, numX = %d, X = %d, t = %g\n",x,left,right,numX,X,t);

		int pos = (X+Y*numX)*2;
		float2 corner0 = (float2)(pFloatParams[pos + 0], pFloatParams[pos + 1]);
		float2 corner1 = (float2)(pFloatParams[pos + 2 + 0], pFloatParams[pos + 2 + 1]);
		pos += numX*2;
		float2 corner2 = (float2)(pFloatParams[pos + 0], pFloatParams[pos + 1]);
		float2 corner3 = (float2)(pFloatParams[pos + 2 + 0], pFloatParams[pos + 2 + 1]);

		r.alphaX = (1.0f-t)*(1.0-u)*corner0.x + t*(1.0f-u)*corner1.x + t*u*corner3.x + (1.0f-t)*u*corner2.x;
		r.alphaY = (1.0f-t)*(1.0-u)*corner0.y + t*(1.0f-u)*corner1.y + t*u*corner3.y + (1.0f-t)*u*corner2.y;
	}
)XYZ";

	if (potential)
		program += R"XYZ(
	r.potential = nan((uint)0); // TODO?
)XYZ";

	if (derivatives)
		program += R"XYZ(
	r.axx = nan((uint)0); // TODO?
	r.ayy = nan((uint)0); // TODO?
	r.axy = nan((uint)0); // TODO?
)XYZ";

	program += R"XYZ(
	return r;
}
)XYZ";

	return program;
}

} // end namespace

