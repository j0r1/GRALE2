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
#include "lensplane.h"
#include "gravitationallens.h"
#include "compositelens.h"
#include "deflectiongridlens.h"
#include "constants.h"
#include <serut/fileserializer.h>
#include <serut/dummyserializer.h>
#include <serut/memoryserializer.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <iostream>

#include "debugnew.h"

using namespace serut;
using namespace std;

namespace grale
{

LensPlane::LensPlane()
{
	m_pLens = 0;
	
	m_numx = 0;
	m_numy = 0;
	m_xstep = 0;
	m_ystep = 0;
	m_init = false;
}

LensPlane::~LensPlane()
{
	delete m_pLens;
}

GravitationalLens *LensPlane::commonInitChecks(const GravitationalLens *pLens, Vector2Dd bl, Vector2Dd tr, 
                                               int xp, int yp, bool copyLens)
{
	if (m_init)
	{
		setErrorString("Already initialized");
		return 0;
	}
	
	if (xp < 2 || yp < 2)
	{
		setErrorString("Number of X and Y points should be at least 2");
		return 0;
	}
	if (xp > 100000 && yp > 100000)
	{
		setErrorString("Number of X and Y points should be at most 100000");
		return 0;
	}
	if (pLens == 0)
	{
		setErrorString("No lens was specified");
		return 0;
	}
	
	double width,height;

	width = tr.getX()-bl.getX();
	height = tr.getY()-bl.getY();

	if (width < 0.0 || height < 0.0)
	{
		setErrorString("Specified coordinates imply a negative width or height");
		return 0;
	}

	GravitationalLens *pLensCopy = 0;
	if (copyLens)
	{
		pLensCopy = pLens->createCopy();
		if (pLensCopy == 0)
		{
			setErrorString("Unable to create a copy of the lens");
			return 0;
		}
	}
	else
		pLensCopy = (GravitationalLens *)pLens; // TODO: make this cleaner, just to signal ok

	m_numx = xp;
	m_numy = yp;

	m_alphas.resize(m_numx*m_numy);
	m_alphaxx.resize(m_numx*m_numy);
	m_alphayy.resize(m_numx*m_numy);
	m_alphaxy.resize(m_numx*m_numy);

	//cerr << m_numx << " " << m_numy << " " << m_alphas.size() << endl;

	m_xstep = width/((double)(m_numx-1));
	m_ystep = height/((double)(m_numy-1));

	m_bottomleft = bl;
	m_topright = tr;

	return pLensCopy;
}

bool LensPlane::init(const GravitationalLens *pLens, Vector2Dd bl, Vector2Dd tr, int xp, int yp)
{
	GravitationalLens *pLensCopy = commonInitChecks(pLens, bl, tr, xp, yp, true);
	if (!pLensCopy)
		return false;

	return initInternal(pLensCopy, bl, tr, xp, yp);
}

bool LensPlane::init(SerializationInterface &lensData, Vector2Dd bl, Vector2Dd tr, int xp, int yp)
{
	GravitationalLens *pLensCopy = 0;
	string errStr;

	if (!GravitationalLens::read(lensData, &pLensCopy, errStr))
	{
		setErrorString("Unable to create lens from data: " + errStr);
		return false;
	}
	if (!commonInitChecks(pLensCopy, bl, tr, xp, yp, false))
	{
		delete pLensCopy;
		return false;
	}

	return initInternal(pLensCopy, bl, tr, xp, yp);
}

bool LensPlane::initInternal(GravitationalLens *pLensCopy, Vector2Dd bl, Vector2Dd tr, int xp, int yp)
{
	double width,height;

	width = tr.getX()-bl.getX();
	height = tr.getY()-bl.getY();

	double minstep = m_xstep;

	if (m_ystep < m_xstep)
		minstep = m_ystep;

	pLensCopy->setDerivativeAngularDistanceScale(minstep/100.0);

	int numtotal = m_numx*m_numy;

	// fill in the arrays
	
	int index = 0;
	int prevPct = -1;
	
	setFeedbackStatus("Building map");
	
	for (int y = 0 ; y < m_numy ; y++)
	{
		for (int x = 0 ; x < m_numx ; x++, index++)
		{
			Vector2Dd theta = getIndexCoordinate(x, y);
			Vector2Dd a;
			double axx,ayy,axy;
			
			if (!pLensCopy->getAlphaVector(theta,&a))
			{
				setErrorString(std::string("Error raytracing theta vector: ") + m_pLens->getErrorString());
				delete pLensCopy;
				return false;
			}
			
			//cerr << a.getX() << " " << a.getY() << endl;

			m_alphas[index] = a;
			pLensCopy->getAlphaVectorDerivatives(theta,axx,ayy,axy);
			m_alphaxx[index] = axx;
			m_alphayy[index] = ayy;
			m_alphaxy[index] = axy;
		}
		
		int pct = (int)(((((double)index)/((double)(numtotal-1))))*100.0+0.5);
		if (pct > 100)
			pct = 100;

		if (pct != prevPct)
		{
			setFeedbackPercentage(pct);
			prevPct = pct;
		}
	}

	setFeedbackStatus("Finished");

	m_init = true;
	m_pLens = pLensCopy;
	
	return true;
}

bool LensPlane::init(const GravitationalLens *pLens, Vector2Dd bl, Vector2Dd tr, int xp, int yp,
					 SerializationInterface &renderedData)
{
	GravitationalLens *pLensCopy = commonInitChecks(pLens, bl, tr, xp, yp, true);
	if (!pLensCopy)
		return false;

	return initInternal(pLensCopy, bl, tr, xp, yp, renderedData);
}
	
bool LensPlane::init(SerializationInterface &lensData, Vector2Dd bl, Vector2Dd tr, int xp, int yp,
                     SerializationInterface &renderedData)
{
	GravitationalLens *pLensCopy = 0;
	string errStr;

	if (!GravitationalLens::read(lensData, &pLensCopy, errStr))
	{
		setErrorString("Unable to create lens from data: " + errStr);
		return false;
	}

	if (!commonInitChecks(pLensCopy, bl, tr, xp, yp, false))
	{
		delete pLensCopy;
		return false;
	}

	return initInternal(pLensCopy, bl, tr, xp, yp, renderedData);
}


bool LensPlane::initInternal(GravitationalLens *pLensCopy, Vector2Dd bl, Vector2Dd tr, int xp, int yp,
		                     SerializationInterface &renderedData)
{
	double width,height;

	width = tr.getX()-bl.getX();
	height = tr.getY()-bl.getY();

	double minstep = m_xstep;

	if (m_ystep < m_xstep)
		minstep = m_ystep;

	pLensCopy->setDerivativeAngularDistanceScale(minstep/100.0);

	int numtotal = m_numx*m_numy;

	// fill in the arrays
	
	int index = 0;
	
	setFeedbackStatus("Importing rendered data");
	
	for (int y = 0 ; y < m_numy ; y++)
	{
		for (int x = 0 ; x < m_numx ; x++, index++)
		{
			double d[5];

			if (!renderedData.readDoubles(d, 5))
			{
				setErrorString("Unable to read deflection data or derivatives from rendered data: " + renderedData.getErrorString());
				delete pLensCopy;
				return false;
			}

			m_alphas[index] = Vector2Dd(d[0], d[1]);
			m_alphaxx[index] = d[2];
			m_alphayy[index] = d[3];
			m_alphaxy[index] = d[4];
		}
	}

	setFeedbackStatus("Finished");

	m_init = true;
	m_pLens = pLensCopy;
	
	return true;
}



#define LENSPLANEID			0x414C504C

bool LensPlane::write(SerializationInterface &si) const
{
	if (!isInit())
	{
		setErrorString("Lens plane is not initialized");
		return false;
	}

	if (!si.writeInt32(LENSPLANEID))
	{
		setErrorString(std::string("Error writing lens plane ID: ") + si.getErrorString());
		return false;
	}
	if (!si.writeDoubles(m_bottomleft.getComponents(), 2))
	{
		setErrorString(std::string("Error writing lens plane bottom-left coordinates: ") + si.getErrorString());
		return false;
	}
	if (!si.writeDoubles(m_topright.getComponents(), 2))
	{
		setErrorString(std::string("Error writing lens plane top-right coordinates: ") + si.getErrorString());
		return false;
	}
	if (!si.writeInt32(m_numx))
	{
		setErrorString(std::string("Error writing lens plane number of X points: ") + si.getErrorString());
		return false;
	}
	if (!si.writeInt32(m_numy))
	{
		setErrorString(std::string("Error writing lens plane number of Y points: ") + si.getErrorString());
		return false;
	}

	int numtotal = m_numx*m_numy;
	int i;
	for (i = 0 ; i < numtotal ; i++)
	{
		double array[5];

		array[0] = m_alphas[i].getX();
		array[1] = m_alphas[i].getY();
		array[2] = m_alphaxx[i];
		array[3] = m_alphayy[i];
		array[4] = m_alphaxy[i];

		if (!si.writeDoubles(array,5))
		{
			setErrorString(std::string("Error writing lens plane map: ") + si.getErrorString());
			return false;
		}
	}

	if (!m_pLens->write(si))
	{
		setErrorString(std::string("Error writing the lens plane's associated lens: ") + m_pLens->getErrorString());
		return false;
	}

	return true;
}

bool LensPlane::read(SerializationInterface &si,LensPlane **ip,std::string &errstr)
{
	int32_t id;
	Vector2Dd bottomleft,topright;
	int32_t numx,numy;
	std::vector<Vector2Dd > alphas;
	std::vector<double> alphaxx, alphayy, alphaxy;
	
	if (!si.readInt32(&id))
	{
		errstr = std::string(std::string("Error reading lens plane ID: ") + si.getErrorString());
		return false;
	}
	if (id != LENSPLANEID)
	{
		errstr = std::string("Read invalid lens plane ID");
		return false;
	}
	if (!si.readDoubles(bottomleft.getComponents(), 2))
	{
		errstr = std::string(std::string("Error reading lens plane bottom-left coordinates: ") + si.getErrorString());
		return false;
	}
	if (!si.readDoubles(topright.getComponents(), 2))
	{
		errstr = std::string(std::string("Error reading lens plane top-right coordinates: ") + si.getErrorString());
		return false;
	}
	if (!si.readInt32(&numx))
	{
		errstr = std::string(std::string("Error reading lens plane number of X points: ") + si.getErrorString());
		return false;
	}
	if (!si.readInt32(&numy))
	{
		errstr = std::string(std::string("Error reading lens plane number of Y points: ") + si.getErrorString());
		return false;
	}

	int32_t numtotal = numx*numy;
	int32_t i;

	alphas.resize(numtotal);
	alphaxx.resize(numtotal);
	alphayy.resize(numtotal);
	alphaxy.resize(numtotal);
	
	for (i = 0 ; i < numtotal ; i++)
	{
		double array[5];

		if (!si.readDoubles(array,5))
		{
			errstr = std::string(std::string("Error reading lens plane map: ") + si.getErrorString());
			return false;
		}
		alphas[i] = Vector2Dd(array[0],array[1]);
		alphaxx[i] = array[2];
		alphayy[i] = array[3];
		alphaxy[i] = array[4];
	}
	
	GravitationalLens *lens;
	std::string errstr2;

	if (!GravitationalLens::read(si,&lens,errstr2))
	{
		errstr = std::string(std::string("Error reading the lens plane's associated lens: ") + errstr2);
		return false;
	}

	double width,height;

	LensPlane *plane = new LensPlane();
	
	width = topright.getX()-bottomleft.getX();
	height = topright.getY()-bottomleft.getY();
	
	plane->m_xstep = width/((double)(numx-1));
	plane->m_ystep = height/((double)(numy-1));
	plane->m_numx = numx;
	plane->m_numy = numy;
	plane->m_bottomleft = bottomleft;
	plane->m_topright = topright;
	plane->m_alphas = alphas;
	plane->m_alphaxx = alphaxx;
	plane->m_alphayy = alphayy;
	plane->m_alphaxy = alphaxy;
	plane->m_pLens = lens;
	plane->m_init = true;

	*ip = plane;
	
	return true;
}

bool LensPlane::load(const std::string &fname,LensPlane **ip,std::string &errstr)
{
	FileSerializer fs;
	bool status;

	if (!fs.open(fname, FileSerializer::ReadOnly))
	{
		errstr = std::string("Error opening lens plane file for reading: ") + fs.getErrorString();
		return false;
	}
	return LensPlane::read(fs,ip,errstr);
}

bool LensPlane::save(const std::string &fname) const
{
	FileSerializer fs;
	bool status;

	if (!fs.open(fname, FileSerializer::WriteOnly))
	{
		setErrorString(std::string("Error opening lens plane file for writing: ") + fs.getErrorString());
		return false;
	}
	return write(fs);
}

bool LensPlane::scaleDeflections(double factor)
{
	CompositeLensParams newParams;

	if (!newParams.addLens(factor, Vector2Dd(0, 0), 0, *m_pLens))
	{
		setErrorString(std::string("Couldn't create scaled version of the lens: ") + newParams.getErrorString());
		return false;
	}

	CompositeLens *pNewLens = new CompositeLens();

	if (!pNewLens->init(m_pLens->getLensDistance(), &newParams))
	{
		setErrorString(std::string("Couldn't create scaled version of the lens: ") + newParams.getErrorString());
		delete pNewLens;
		return false;
	}

	// Ok, got the new lens. Now store it and scale the deflections

	delete m_pLens;
	m_pLens = pNewLens;

	int index = 0;
	for (int y = 0 ; y < m_numy ; y++)
	{
		for (int x = 0 ; x < m_numx ; x++)
		{
			m_alphas[index] *= factor;
			m_alphaxx[index] *= factor;
			m_alphayy[index] *= factor;
			m_alphaxy[index] *= factor;
			
			index++;
		}
	}

	return true;
}

GravitationalLens *LensPlane::createDeflectionGridLens() const
{
	int totalNum = m_numx*m_numy;
	std::vector<double> alphaX(totalNum), alphaY(totalNum);

	for (int i = 0 ; i < totalNum ; i++)
	{
		alphaX[i] = m_alphas[i].getX();
		alphaY[i] = m_alphas[i].getY();
	}

	DeflectionGridLensParams lensParams(alphaX, alphaY, m_numx, m_numy, m_bottomleft, m_topright);
	DeflectionGridLens *pNewLens = new DeflectionGridLens();

	if (!pNewLens->init(m_pLens->getLensDistance(), &lensParams))
	{
		setErrorString(std::string("Couldn't initialize deflection grid lens: ") + pNewLens->getErrorString());
		delete pNewLens;
		return 0;
	}

	return pNewLens;
}

} // end namespace

