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
#include "imageplane.h"
#include "lensplane.h"
#include "sourceimage.h"
#include "gravitationallens.h"
#include "contourfinder.h"
#include <string.h>
#include <iostream>

#include "debugnew.h"

using namespace std;

namespace grale
{

ImagePlane::ImagePlane()
{
	m_pRealLens = 0;
	m_pGridLens = 0;

	numx = 0;
	numy = 0;
	xstep = 0;
	ystep = 0;
	m_pAlphaxxFunction = 0;
	m_pAlphayyFunction = 0;
	m_pAlphaxyFunction = 0;

	m_Dds = 0;
	m_Ds = 0;
}

ImagePlane::~ImagePlane()
{
	delete m_pRealLens;
	delete m_pGridLens;
	delete m_pAlphaxxFunction;
	delete m_pAlphayyFunction;
	delete m_pAlphaxyFunction;
}

bool ImagePlane::init(const LensPlane *lensplane, double Ds, double Dds)
{
	if (m_pRealLens)
	{
		setErrorString("Already initialized");
		return false;
	}
	if (Ds < 0 || Dds < 0)
	{
		setErrorString("Distances must be positive");
		return false;
	}

	if (!lensplane || !lensplane->isInit())
	{
		setErrorString("Can't use an non-existant or uninitialized lens plane");
		return false;
	}
	
	m_pRealLens = lensplane->getLens()->createCopy();
	if (m_pRealLens == 0)
	{
		setErrorString("Unable to create copy of the lens: " + lensplane->getLens()->getErrorString());
		return false;
	}
	m_pGridLens = lensplane->createDeflectionGridLens();
	if (m_pGridLens == 0)
	{
		setErrorString("Unable to create deflection grid based lens: " + m_pGridLens->getErrorString());
		return false;
	}

	bottomleft = lensplane->getBottomLeft();
	topright = lensplane->getTopRight();
	xstep = lensplane->getXStep();
	ystep = lensplane->getYStep();
	numx = lensplane->getNumXPoints();
	numy = lensplane->getNumYPoints();

	m_betas.resize(numx*numy);
	m_alphaxx.resize(numx*numy);
	m_alphayy.resize(numx*numy);
	m_alphaxy.resize(numx*numy);
	vector<double> inversemagnifications(numx*numy);

	m_Ds = Ds;
	m_Dds = Dds;
	double factor = Dds/Ds;

	//setFeedbackStatus("Creating map from lens plane");
	
	for (int y = 0 ; y < numy ; y++)
	{
		for (int x = 0 ; x < numx ; x++)
		{
			Vector2Dd theta = getIndexCoordinate(x,y);
			m_betas[x+y*numx] = theta - factor*lensplane->getAlpha(x,y);
			
			double axx,ayy,axy;
			
			lensplane->getAlphaDerivatives(x,y,axx,ayy,axy);
			axx *= factor;
			ayy *= factor;
			axy *= factor;
			
			m_alphaxx[x+y*numx] = axx;
			m_alphayy[x+y*numx] = ayy;
			m_alphaxy[x+y*numx] = axy;

			// This is to calculate the critical lines and caustics, will not be
			// stored permanently
			inversemagnifications[x+y*numx] = (1.0-axx)*(1.0-ayy)-axy*axy;
		}
	}

	ContourFinder contourFinder(inversemagnifications, bottomleft, topright, numx, numy);
	m_criticalLines = contourFinder.findContour(0);
	m_caustics.clear();
	for (auto &critPart : m_criticalLines)
	{
		vector<Vector2Dd> caustPart;

		for (auto theta : critPart)
		{
			Vector2Dd beta;
			bool inplane = false;

			if (m_pGridLens->traceTheta(m_Ds, m_Dds, theta, &beta))
				caustPart.push_back(beta);
			else
			{
				if (!caustPart.empty())
				{
					m_caustics.push_back(caustPart);
					caustPart.clear();
				}
			}
		}
		if (!caustPart.empty())
			m_caustics.push_back(caustPart);
	}

	m_pAlphaxxFunction = new GridFunction(&(m_alphaxx[0]), bottomleft, topright, numx, numy);
	m_pAlphayyFunction = new GridFunction(&(m_alphayy[0]), bottomleft, topright, numx, numy);
	m_pAlphaxyFunction = new GridFunction(&(m_alphaxy[0]), bottomleft, topright, numx, numy);
	
	//setFeedbackStatus("Finished");
	return true;	
}

double ImagePlane::getImageIntensityAccurate(const vector<SourceImage *> &sources, 
		                                     int xPixel, int yPixel, int subSamples) const
{
	int index[4];
	double axx[4];
	double ayy[4];
	double axy[4];
	Vector2Dd beta[4];
	Vector2Dd theta[4];
	
	theta[0] = getIndexCoordinate(xPixel, yPixel);
	theta[1] = getIndexCoordinate(xPixel + 1, yPixel);
	theta[2] = getIndexCoordinate(xPixel, yPixel + 1);
	theta[3] = getIndexCoordinate(xPixel + 1, yPixel + 1);

	index[0] = xPixel + yPixel*numx;
	index[1] = index[0] + 1;
	index[2] = index[0] + numx;
	index[3] = index[2] + 1;

	for (int i = 0 ; i < 4 ; i++)
	{
		beta[i] = m_betas[index[i]];
		axx[i] = m_alphaxx[index[i]];
		axy[i] = m_alphaxy[index[i]];
		ayy[i] = m_alphayy[index[i]];
	}

	Vector2Dd v1 = (beta[0]-beta[3])/2.0;
	Vector2Dd v2 = (beta[1]-beta[2])/2.0;
	Vector2Dd b = (beta[0]+beta[1]+beta[2]+beta[3])/4.0;

	if (!isSourceInRange(sources, b, MAX(v1.getLength(), v2.getLength())*1.1))
		return 0;

	// First we'll calculate the surface brightness, using subsampling if requested
	
	int numSub = (int)(SQRT((double)subSamples)+0.5);

	if (numSub < 1)
		numSub = 1;

	double sum = 0;
	double ySubStep = ystep/(double)numSub;
	double xSubStep = xstep/(double)numSub;

	for (int i = 0 ; i < numSub ; i++)
	{
		double dy = ((double)i)*ySubStep + ySubStep/2.0;
		double thetay = theta[0].getY() + dy;

		for (int j = 0 ; j < numSub ; j++)
		{
			double dx = ((double)j)*xSubStep + xSubStep/2.0;
			double thetax = theta[0].getX() + dx;

			Vector2Dd b[4];

			for (int k = 0 ; k < 4 ; k++)
			{
				double dthetax = thetax - theta[k].getX();
				double dthetay = thetay - theta[k].getY();

				double bx = beta[k].getX() + (1.0 - axx[k])*dthetax - axy[k]*dthetay;
				double by = beta[k].getY() - axy[k]*dthetax + (1.0 - ayy[k])*dthetay;
				
				b[k] = Vector2Dd(bx, by);
			}
			
			double u = dx/xstep;
			double v = dy/ystep;

			Vector2Dd B = b[0]*(1.0-u)*(1.0-v) + b[1]*u*(1.0-v) + b[2]*(1.0-u)*v + b[3]*u*v;

			sum += getSurfaceBrightness(sources, B); // surface brightness is conserved
		}
	}

	sum /= (double)(numSub*numSub);
	
	// Then we'll process the point sources
	
	int refIdx[2] = {0, 3};

	for (int i = 0 ; i < 2 ; i++)
	{
		int k = refIdx[i];
		vector<Vector2DdPlus> pointSourceInfo;
		
		getPointSourceIntensities(sources, Triangle2D<double>(beta[k], beta[1], beta[2]), pointSourceInfo);
		
		Vector2Dd p = beta[1]-beta[k];
		Vector2Dd q = beta[2]-beta[k];
		double daxx1 = axx[1]-axx[k];
		double daxx2 = axx[2]-axx[k];
		double dayy1 = ayy[1]-ayy[k];
		double dayy2 = ayy[2]-ayy[k];
		double daxy1 = axy[1]-axy[k];
		double daxy2 = axy[2]-axy[k];

		double det = p.getX()*q.getY()-p.getY()*q.getX();
		double Ax = q.getY()/det;
		double Ay = -q.getX()/det;
		double Bx = -p.getY()/det;
		double By = p.getX()/det;
		
		for (auto it = pointSourceInfo.begin() ; it != pointSourceInfo.end() ; it++)
		{
			double dx = (*it).getX() - beta[k].getX();
			double dy = (*it).getY() - beta[k].getY();

			double u = Ax*dx+Ay*dy;
			double w = Bx*dx+By*dy;
			
			double axx0 = axx[k] + daxx1*u + daxx2*w;
			double ayy0 = ayy[k] + dayy1*u + dayy2*w;
			double axy0 = axy[k] + daxy1*u + daxy2*w;

			double inverseMagnification = ABS((1.0-axx0)*(1.0-ayy0)-axy0*axy0);

			sum += (*it).getValue()/inverseMagnification;
		}
	}
	
	return sum;
}

double ImagePlane::getSourceIntensityAccurate(const vector<SourceImage *> &sources, 
		                                      int xPixel, int yPixel, int subSamples) const
{
	int index[4];
	Vector2Dd theta[4];
	
	theta[0] = getIndexCoordinate(xPixel, yPixel);
	theta[1] = getIndexCoordinate(xPixel + 1, yPixel);
	theta[2] = getIndexCoordinate(xPixel, yPixel + 1);
	theta[3] = getIndexCoordinate(xPixel + 1, yPixel + 1);

	Vector2Dd v1 = (theta[0]-theta[3])/2.0;
	Vector2Dd v2 = (theta[1]-theta[2])/2.0;
	Vector2Dd t = (theta[0]+theta[1]+theta[2]+theta[3])/4.0;

	if (!isSourceInRange(sources, t, MAX(v1.getLength(), v2.getLength())*1.1))
		return 0;

	// First we'll calculate the surface brightness, using subsampling if requested
	
	int numSub = (int)(SQRT((double)subSamples)+0.5);

	if (numSub < 1)
		numSub = 1;

	double sum = 0;
	double ySubStep = ystep/(double)numSub;
	double xSubStep = xstep/(double)numSub;

	for (int i = 0 ; i < numSub ; i++)
	{
		double dy = ((double)i)*ySubStep + ySubStep/2.0;
		double thetay = theta[0].getY() + dy;

		for (int j = 0 ; j < numSub ; j++)
		{
			double dx = ((double)j)*xSubStep + xSubStep/2.0;
			double thetax = theta[0].getX() + dx;

			sum += getSurfaceBrightness(sources, Vector2Dd(thetax, thetay));
		}
	}

	sum /= (double)(numSub*numSub);
	
	// Then we'll process the point sources
	
	vector<Vector2DdPlus> pointSourceInfo;
	
	int refIdx[2] = {0, 3};

	for (int i = 0 ; i < 2 ; i++)
	{
		int k = refIdx[i];
		
		getPointSourceIntensities(sources,Triangle2D<double>(theta[k], theta[1], theta[2]), pointSourceInfo);
		
		for (auto it = pointSourceInfo.begin() ; it != pointSourceInfo.end() ; it++)
			sum += (*it).getValue();
	}
	
	return sum;
}

bool ImagePlane::traceBeta(Vector2Dd beta, vector<Vector2Dd> &thetaPoints, int numIterations, bool check) const
{
	thetaPoints.clear();

	for (int y = 0 ; y < numy-1 ; y++)
	{
		for (int x = 0 ; x < numx-1 ; x++)
		{
			Vector2Dd b1 = traceIndexCoordinate(x, y);
			Vector2Dd b2 = traceIndexCoordinate(x+1, y);
			Vector2Dd b3 = traceIndexCoordinate(x, y+1);
			Vector2Dd b4 = traceIndexCoordinate(x+1, y+1);
			Vector2Dd t1 = getIndexCoordinate(x, y);
			Vector2Dd t2 = getIndexCoordinate(x+1, y);
			Vector2Dd t3 = getIndexCoordinate(x, y+1);
			Vector2Dd t4 = getIndexCoordinate(x+1, y+1);

			Triangle2D<double> triangle1(b1, b2, b3);
			Triangle2D<double> triangle2(b4, b2, b3);

			if (triangle1.isInside(beta) || triangle2.isInside(beta))
			{
				Vector2Dd theta;
				
				if (!refinePosition(beta, (t1+t4)/2.0, theta, numIterations))
					return false;

				thetaPoints.push_back(theta);
			}
		}
	}

	if (check)
	{
		double pixelScale = SQRT(xstep*xstep+ystep*ystep);

		auto it = thetaPoints.begin();
		
		while (it != thetaPoints.end())
		{
			Vector2Dd theta = *it;
			Vector2Dd newBeta;
			bool ok = true;
		
			if (!m_pRealLens->traceTheta(m_Ds, m_Dds, theta, &newBeta))
				ok = false;
			else
			{
				//double invMag = m_pLens->getInverseMagnification(sourceplane->getD_s(), sourceplane->getD_ds(), theta);
				//double lengthResize = SQRT(invMag); // magnification relates areas, we'll take the square root to compare distances
				
				Vector2Dd diff = newBeta-beta;

				//if (diff.getLength() > lengthResize*pixelScale*2.0) // to be on the safe side
				if (diff.getLength() > pixelScale*2.0) // to be on the safe side
					ok = false;
			}

			if (ok)
				it++;
			else
				it = thetaPoints.erase(it);
		}
	}
	
	return true;
}


bool ImagePlane::refinePosition(Vector2Dd beta, 
                                Vector2Dd startTheta,
	 		        Vector2Dd &theta,
				int numIterations) const
{	
	grale::Vector2Dd newTheta = startTheta;
	grale::Vector2Dd newBeta;
	
	double prevDiff = 1e30;
	bool done = false;
		
	if (!m_pRealLens->traceTheta(m_Ds, m_Dds, newTheta, &newBeta))
	{
		setErrorString(string("Unable to trace a point to the source plane: ") + m_pRealLens->getErrorString());
		return false;
	}

	double fraction = m_Dds/m_Ds;
	grale::Vector2Dd diffBeta = beta - newBeta;
	grale::Vector2Dd closestPoint = newTheta;
	prevDiff = diffBeta.getLength();

	for (int i = 0 ; i < numIterations ; i++)
	{
		double mxx, myy, mxy;
		
		if (!m_pRealLens->getAlphaVectorDerivatives(newTheta, mxx, myy, mxy))
		{
			setErrorString(m_pRealLens->getErrorString());
			return false;
		}

		mxx = (1.0-fraction*mxx);
		myy = (1.0-fraction*myy);
		mxy = -fraction*mxy;

		double denom = mxx*myy-mxy*mxy;
		double dx = (myy*diffBeta.getX() - mxy*diffBeta.getY())/denom;
		double dy = (-mxy*diffBeta.getX() + mxx*diffBeta.getY())/denom;

		newTheta += grale::Vector2Dd(dx, dy);
		if (!m_pRealLens->traceTheta(m_Ds, m_Dds, newTheta, &newBeta))
		{	
			setErrorString(string("Unable to trace a point to the source plane: ") + m_pRealLens->getErrorString());
			return false;
		}

		diffBeta = beta - newBeta;

		if (diffBeta.getLength() < prevDiff)
		{
			prevDiff = diffBeta.getLength();
			closestPoint = newTheta;
		}
	}

	theta = closestPoint;
	return true;
}

double ImagePlane::getSurfaceBrightness(const vector<SourceImage *> &sources, Vector2Dd beta)
{
	double sum = 0;

	for (auto it = sources.begin() ; it != sources.end() ; it++)
	{
		const SourceImage *pSrc = *it;
		
		if (pSrc->getSourceType() != SourceImage::Point) // Don't consider point sources here
			sum += pSrc->getIntensity(beta);
	}
	return sum;
}

void ImagePlane::getPointSourceIntensities(const vector<SourceImage *> &sources, Triangle2Dd area, 
                                           vector<Vector2DdPlus> &pointSourceInfo)
{
	pointSourceInfo.clear();
	
	for (auto it = sources.begin() ; it != sources.end() ; it++)
	{
		const SourceImage *pSrc = (*it);
		
		if (pSrc->getSourceType() == SourceImage::Point) // Only consider point sources here
		{
			Vector2Dd position = pSrc->getAngularPosition();

			if (area.isInside(position))
				pointSourceInfo.push_back(Vector2DdPlus(position.getX(), position.getY(), pSrc->getIntensity(position)));
		}
	}
}

bool ImagePlane::isSourceInRange(const vector<SourceImage *> &sources, Vector2Dd beta, double radius)
{
	for (auto it = sources.begin() ; it != sources.end() ; it++)
		if ((*it)->isSourceInRange(beta, radius))
			return true;

	return false;
}

inline void calculateTrianglePointCoefficients(Vector2Dd V, Vector2Dd v1, Vector2Dd v2, Vector2Dd v3, double &a, double &b)
{
	a = 0;
	b = 0;

	Vector2Dd D1 = v2-v1;
	Vector2Dd D2 = v3-v1;
	Vector2Dd D = V-v1;

	double dx1 = D1.getX();
	double dy1 = D1.getY();
	double dx2 = D2.getX();
	double dy2 = D2.getY();
	double X = D.getX();
	double Y = D.getY();

	double det = dx1*dy2-dx2*dy1;
	if (det == 0)
		return; // we've set a=0, b=0, so hopefully the effect won't be too dramatic

	a = (dy2*X - dx2*Y)/det;
	b = (-dy1*X + dx1*Y)/det;
}

inline void log(Vector2Dd beta, Vector2Dd b1, Vector2Dd b2, Vector2Dd b3,
		        Vector2Dd theta, Vector2Dd t1, Vector2Dd t2, Vector2Dd t3)
{
	cout << "# Beta" << endl;
	cout << b1.getX() << " " << b1.getY() << endl;
	cout << beta.getX() << " " << beta.getY() << endl;
	cout << endl << endl;
	cout << "# Beta triang" << endl;
	cout << b1.getX() << " " << b1.getY() << endl;
	cout << b2.getX() << " " << b2.getY() << endl;
	cout << b3.getX() << " " << b3.getY() << endl;
	cout << b1.getX() << " " << b1.getY() << endl;
	cout << endl << endl;

	cout << "# Theta" << endl;
	cout << t1.getX() << " " << t1.getY() << endl;
	cout << theta.getX() << " " << theta.getY() << endl;
	cout << endl << endl;
	cout << "# Theta triang" << endl;
	cout << t1.getX() << " " << t1.getY() << endl;
	cout << t2.getX() << " " << t2.getY() << endl;
	cout << t3.getX() << " " << t3.getY() << endl;
	cout << t1.getX() << " " << t1.getY() << endl;
	cout << endl << endl;
}

bool ImagePlane::staticTraceBetaApproximately(Vector2Dd beta, vector<Vector2Dd> &thetaPoints,
										     const vector<Vector2Dd> &betaMap,
	                                         Vector2Dd bottomLeft, Vector2Dd topRight, int numX, int numY,
											 std::string &errorString)
{
	if (numX < 2 || numY < 2)
	{
		errorString = "Both X and Y dimensions must be at least 2";
		return false;
	}
	if (numX*numY != (int)betaMap.size())
	{
		errorString = "Specified dimensions are not compatible with the beta map size";
		return false;
	}

	thetaPoints.clear();

	double W = (topRight.getX() - bottomLeft.getX())/(numX-1);
	double H = (topRight.getY() - bottomLeft.getY())/(numY-1);

	for (int y = 0 ; y < numY-1 ; y++)
	{
		for (int x = 0 ; x < numX-1 ; x++)
		{
			Vector2Dd b1 = betaMap[x + y*numX];
			Vector2Dd b2 = betaMap[x+1 + y*numX];
			Vector2Dd b3 = betaMap[x + (y+1)*numX];
			Vector2Dd b4 = betaMap[(x+1) + (y+1)*numX];
			Vector2Dd t1(x, y);
			Vector2Dd t2(x+1, y);
			Vector2Dd t3(x, y+1);
			Vector2Dd t4(x+1, y+1);
			double a, b;

			Triangle2D<double> triangle1(b1, b2, b3);
			Triangle2D<double> triangle2(b4, b2, b3);

			if (triangle1.isInside(beta))
			{
				calculateTrianglePointCoefficients(beta, b1, b2, b3, a, b);
				Vector2Dd theta = t1 + (t2-t1)*a + (t3-t1)*b;
				//log(beta, b1, b2, b3, theta, t1, t2, t3);
				thetaPoints.push_back(Vector2Dd(theta.getX()*W + bottomLeft.getX(), theta.getY()*H + bottomLeft.getY()));
			}

			if (triangle2.isInside(beta))
			{
				calculateTrianglePointCoefficients(beta, b4, b2, b3, a, b);
				Vector2Dd theta = t4 + (t2-t4)*a + (t3-t4)*b;
				//log(beta, b4, b2, b3, theta, t4, t2, t3);
				thetaPoints.push_back(Vector2Dd(theta.getX()*W + bottomLeft.getX(), theta.getY()*H + bottomLeft.getY()));
			}
		}
	}

	return true;
}

} // end namespace

