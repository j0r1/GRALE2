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

#ifndef GRALE_LENSPLANE_H

#define GRALE_LENSPLANE_H

#include "graleconfig.h"
#include "vector2d.h"
#include <serut/serializationinterface.h>
#include <stdint.h>
#include <stdio.h>
#include <list>
#include <string>
#include <vector>
#include <memory>

namespace grale
{

class GravitationalLens;
class OpenCLKernel;

class GRALE_IMPORTEXPORT LensPlane : public errut::ErrorBase
{
public:
	LensPlane();
	virtual ~LensPlane();

	bool init(serut::SerializationInterface &lensData, Vector2Dd bottomLeft, Vector2Dd topRight,
	          int xpoints,int ypoints);

	bool init(serut::SerializationInterface &lensData, Vector2Dd bottomLeft, Vector2Dd topRight,
			  int xpoints, int ypoints, serut::SerializationInterface &renderedData);

	bool isInit() const									{ return m_init; }
	const GravitationalLens *getLens() const						{ return m_pLens.get(); }

	GravitationalLens *createDeflectionGridLens() const;
	
	bool scaleDeflections(double factor);
	int getNumXPoints() const								{ return m_numx; }
	int getNumYPoints() const								{ return m_numy; }
	Vector2Dd getBottomLeft() const							{ return m_bottomleft; }
	Vector2Dd getTopRight() const							{ return m_topright; }
	double getXStep() const									{ return m_xstep; }
	double getYStep() const									{ return m_ystep; }
	
	Vector2Dd getAlpha(int xpos, int ypos) const					{ return m_alphas[xpos+ypos*m_numx]; }
	void getAlphaDerivatives(int xpos, int ypos, double &axx, double &ayy, double &axy) const { axx = m_alphaxx[xpos+ypos*m_numx]; ayy = m_alphayy[xpos+ypos*m_numx]; axy = m_alphaxy[xpos+ypos*m_numx]; }
	
	Vector2Dd getIndexCoordinate(int xpos,int ypos) const;
	bool getIndices(Vector2Dd v,int *xpos,int *ypos) const;

	bool write(serut::SerializationInterface &si) const;
	static bool read(serut::SerializationInterface &si,LensPlane **ip,std::string &errstr);
	static bool load(const std::string &fname,LensPlane **ip,std::string &errstr);
	bool save(const std::string &fname) const;

	const std::vector<Vector2Dd> &getAlphas() const						{ return m_alphas; }
	const std::vector<double> &getAlphaXXs() const						{ return m_alphaxx; }
	const std::vector<double> &getAlphaXYs() const						{ return m_alphaxy; }
	const std::vector<double> &getAlphaYYs() const						{ return m_alphayy; }
protected:
	virtual void setFeedbackStatus(const std::string &msg)					{ }
	virtual void setFeedbackPercentage(int pct)						{ }
private:
	bool commonInitChecks(const GravitationalLens *pLens, Vector2Dd bl, Vector2Dd tr, 
                                        int xp, int yp);

	bool initInternal(std::unique_ptr<GravitationalLens> pLensCopy, Vector2Dd bl, Vector2Dd tr, int xp, int yp);
	bool initInternal(std::unique_ptr<GravitationalLens> pLensCopy, Vector2Dd bl, Vector2Dd tr, int xp, int yp,
	                  serut::SerializationInterface &renderedData);

	std::unique_ptr<GravitationalLens> m_pLens;
	
	std::vector<Vector2Dd > m_alphas;
	std::vector<double> m_alphaxx, m_alphayy, m_alphaxy;
	
	Vector2Dd m_bottomleft;
	Vector2Dd m_topright;
	double m_xstep, m_ystep;
	int m_numx, m_numy;
	bool m_init;
};

inline Vector2Dd LensPlane::getIndexCoordinate(int xpos,int ypos) const
{ 
	return Vector2Dd((((double)xpos))*m_xstep + m_bottomleft.getX(),
			        (((double)ypos))*m_ystep + m_bottomleft.getY());
}

inline bool LensPlane::getIndices(Vector2Dd v,int *xpos,int *ypos) const
{
	int x, y;

	x = (int)(((v.getX()-m_bottomleft.getX())/m_xstep)+0.5);
	y = (int)(((v.getY()-m_bottomleft.getY())/m_ystep)+0.5);

	if (x < 0 || y < 0 || x >= m_numx || y >= m_numy)
		return false;

	*xpos = x;
	*ypos = y;
	return true;
}

} // end namespace

#endif // GRALE_LENSPLANE_H

