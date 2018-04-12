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

#ifndef GRALE_IMAGEPLANE_H

#define GRALE_IMAGEPLANE_H

#include "graleconfig.h"
#include "vector2d.h"
#include "triangle2d.h"
#include "gridfunction.h"
#include <serut/serializationinterface.h>
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <string>

namespace grale
{

class GravitationalLens;
class LensPlane;
class SourceImage;

class GRALE_IMPORTEXPORT ImagePlane : public errut::ErrorBase
{
public:
	ImagePlane();
	virtual ~ImagePlane();

	bool init(const LensPlane *lensplane, double Ds, double Dds);
	bool isInit() const										{ if (m_pRealLens == 0) return false; return true; }
	const GravitationalLens *getLens() const				{ return m_pRealLens; }

    double getDs() const                                    { return m_Ds; }
    double getDds() const                                   { return m_Dds; }
	
	int getNumXPoints() const								{ return numx; }
	int getNumYPoints() const								{ return numy; }
	int getNumXPixels() const								{ return numx-1; }
	int getNumYPixels() const								{ return numy-1; }
	Vector2Dd getBottomLeft() const							{ return bottomleft; }
	Vector2Dd getTopRight() const							{ return topright; }
	double getXStep() const									{ return xstep; }
	double getYStep() const									{ return ystep; }
	
	double getImageIntensityAccurate(const std::vector<SourceImage *> &sources, 
			                         int xPixel, int yPixel, int subSamples = 10) const;
	double getSourceIntensityAccurate(const std::vector<SourceImage *> &sources, 
			                          int xPixel, int yPixel, int subSamples = 10) const;
		
	double getInverseMagnification(int xpos,int ypos) const					{ double axx = m_alphaxx[xpos+ypos*numx]; double ayy = m_alphayy[xpos+ypos*numx]; double axy = m_alphaxy[xpos+ypos*numx]; return (1.0-axx)*(1.0-ayy)-axy*axy; }
	void getInverseMagnificationComponents(int xpos,int ypos, double &kappa, double &gamma1,  double &gamma2) const				{ double axx = m_alphaxx[xpos+ypos*numx]; double ayy = m_alphayy[xpos+ypos*numx]; double axy = m_alphaxy[xpos+ypos*numx]; kappa = 0.5*(axx+ayy); gamma1 = 0.5*(axx-ayy); gamma2 = axy; }
	
	const std::vector<std::vector<Vector2Dd> > &getCriticalLineSegments() const	{ return m_criticalLines; }
	const std::vector<std::vector<Vector2Dd> > &getCausticSegments() const		{ return m_caustics; }
	
	Vector2Dd getIndexCoordinate(int xpos, int ypos) const;
	Vector2Dd getPixelCoordinate(int xPixel, int yPixel) const;

	Vector2Dd traceIndexCoordinate(int xpos, int ypos) const;
	bool getIndices(Vector2Dd v,int *xpos,int *ypos) const;

	static bool staticTraceBetaApproximately(Vector2Dd beta, std::vector<Vector2Dd> &thetaPoints,
											 const std::vector<Vector2Dd> &betaMap,
	                                         Vector2Dd bottomLeft, Vector2Dd topRight, int numX, int numY,
											 std::string &errorString);

	bool traceBeta(Vector2Dd beta, std::vector<Vector2Dd> &thetaPoints, int numIterations = 4, bool check = true) const;

	bool traceThetaApproximately(Vector2Dd theta, Vector2Dd *pBeta) const;
private:
	void searchCritCaust(double *inversemagnifications);
	bool refinePosition(Vector2Dd beta, Vector2Dd startTheta, Vector2Dd &theta, int numIterations) const;

	static bool isSourceInRange(const std::vector<SourceImage *> &sources, Vector2Dd beta, double radius);
	static void getPointSourceIntensities(const std::vector<SourceImage *> &sources, Triangle2Dd area, 
                                          std::vector<Vector2DdPlus> &pointSourceInfo);
	static double getSurfaceBrightness(const std::vector<SourceImage *> &sources, Vector2Dd beta);

	GravitationalLens *m_pRealLens;
	GravitationalLens *m_pGridLens;
	double m_Dds, m_Ds;
	
	std::vector<Vector2Dd > m_betas;
	std::vector<double> m_alphaxx;
	std::vector<double> m_alphayy;
	std::vector<double> m_alphaxy;
	GridFunction *m_pAlphaxxFunction;
	GridFunction *m_pAlphayyFunction;
	GridFunction *m_pAlphaxyFunction;
	std::vector<std::vector<Vector2Dd> > m_criticalLines;
	std::vector<std::vector<Vector2Dd> > m_caustics;
	
	Vector2Dd bottomleft;
	Vector2Dd topright;
	double xstep,ystep;
	int numx,numy;
};

inline Vector2Dd ImagePlane::getIndexCoordinate(int xpos,int ypos) const
{ 
	return Vector2Dd((((double)xpos))*xstep+bottomleft.getX(),(((double)ypos))*ystep+bottomleft.getY());
}

inline Vector2Dd ImagePlane::getPixelCoordinate(int xPixel, int yPixel) const
{
	return Vector2Dd((((double)xPixel))*xstep+bottomleft.getX(),(((double)yPixel))*ystep+bottomleft.getY())+Vector2Dd(xstep/2.0, ystep/2.0);
}

inline Vector2Dd ImagePlane::traceIndexCoordinate(int xpos,int ypos) const
{ 
	return m_betas[xpos+ypos*numx];
}

inline bool ImagePlane::getIndices(Vector2Dd v,int *xpos,int *ypos) const
{
	int x,y;

	x = (int)(((v.getX()-bottomleft.getX())/xstep)+0.5);
	y = (int)(((v.getY()-bottomleft.getY())/ystep)+0.5);

	if (x < 0 || y < 0 || x >= numx || y >= numy)
		return false;
	*xpos = x;
	*ypos = y;
	return true;
}

/*
inline double ImagePlane::getInverseMagnificationApprox(Vector2Dd theta) const
{
	double axx = (*m_pAlphaxxFunction)(theta);
	double ayy = (*m_pAlphayyFunction)(theta);
	double axy = (*m_pAlphaxyFunction)(theta);
	double inverseMagnification = ABS((1.0-axx)*(1.0-ayy)-axy*axy);

	return inverseMagnification;
}
*/

} // end namespace

#endif // GRALE_IMAGEPLANE_H

