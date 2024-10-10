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

/**
 * \file compositelens.h
 */

#ifndef GRALE_COMPOSITELENS_H

#define GRALE_COMPOSITELENS_H

#include "graleconfig.h"
#include "gravitationallens.h"
#include <vector>
#include <map>

namespace grale
{

class GRALE_IMPORTEXPORT CompositeLensParams : public GravitationalLensParams
{
public:
	class LensInfo;
public:
	CompositeLensParams();
	~CompositeLensParams();
	bool addLens(double factor, Vector2D<double> position, double angle, const GravitationalLens &lens);
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	const std::vector<std::shared_ptr<LensInfo>> &getLensInfo() const						{ return m_lensInfo; }
	std::unique_ptr<GravitationalLensParams> createCopy() const;
private:
	std::vector<std::shared_ptr<LensInfo>> m_lensInfo;
};

// TODO
// WARNING!!! No check on the different angular diameter distances of the sublenses is performed!!!
// 
/** This class allows you to combine several other lenses.
 *  \f[ 
 *  	\left(\begin{array}{cc}
 *  	  \hat{\alpha}_{x,x} & \hat{\alpha}_{x,y} \\
 *  	  \hat{\alpha}_{x,y} & \hat{\alpha}_{y,y}
 *  	\end{array}\right)
 *      =
 *  	\left(\begin{array}{cc}
 *  	  \cos\theta & -\sin\theta \\
 *  	  \sin\theta & \cos\theta 
 *  	\end{array}\right)
 *  	\left(\begin{array}{cc}
 *  	  \hat{\alpha}'_{x,x} & \hat{\alpha}'_{x,y} \\
 *  	  \hat{\alpha}'_{x,y} & \hat{\alpha}'_{y,y}
 *  	\end{array}\right)
 *  	\left(\begin{array}{cc}
 *  	  \cos\theta & \sin\theta \\
 *  	  -\sin\theta & \cos\theta 
 *  	\end{array}\right)
 *  \f]
 *
 *  \f{eqnarray*}
 *  	\hat{\alpha}_{x,x} & = & \hat{\alpha}'_{x,x}\cos^2\theta + \hat{\alpha}'_{y,y}\sin^2\theta - 2 \hat{\alpha}'_{x,y} \sin\theta\cos\theta \\
 *  	\hat{\alpha}_{y,y} & = & \hat{\alpha}'_{x,x}\sin^2\theta + \hat{\alpha}'_{y,y}\cos^2\theta + 2 \hat{\alpha}'_{x,y} \sin\theta\cos\theta \\
 *  	\hat{\alpha}_{x,y} & = & (\hat{\alpha}'_{x,x}-\hat{\alpha}'_{y,y})\sin\theta\cos\theta + \hat{\alpha}'_{x,y}(2\cos^2\theta-1)
 *  \f}
 */
class GRALE_IMPORTEXPORT CompositeLens : public GravitationalLens
{
public:
	CompositeLens();
	~CompositeLens();

	static CompositeLens *cast(GravitationalLens *pLens);
	static const CompositeLens *cast(const GravitationalLens *pLens);

	bool getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	bool getAlphaVectorSecondDerivatives(Vector2D<double> theta, double &axxx, double &ayyy, double &axxy, double &ayyx) const;
	void setDerivativeAngularDistanceScale(double distanceScale);
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;

	int getNumberOfSubLenses() const								{ return m_lenses.size(); }
	const GravitationalLens *getSubLens(int i) const				{ return m_lenses[i].get(); }
	Vector2D<double> getSubLensPosition(int i) const				{ return m_positions[i]; }
	double getSubLensAngle(int i) const								{ return m_angles[i]; }
	double getSubLensAngleOriginal(int i) const						{ return m_origAngles[i]; }
	double getSubLensFactor(int i) const							{ return m_factors[i]; }

	bool getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const;
	bool getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const;
	bool getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const;
	std::string getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives = true, bool potential = true) const override;

	static std::string getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, const std::vector<std::string> &otherRoutineNames, int maxRecursionCount, 
		                     bool derivatives, bool potential);

	// Same function with different name - newer cython seemed to get confused...
	static std::string getCLProgram_static(double deflectionScale, double potentialScale, std::string &subRoutineName, const std::vector<std::string> &otherRoutineNames, int maxRecursionCount,
		                     bool derivatives, bool potential) { return getCLProgram(deflectionScale, potentialScale, subRoutineName, otherRoutineNames, maxRecursionCount, derivatives, potential); }

	// Returns maxRecursionCount
	int findCLSubroutines(double deflectionScale, double potentialScale, std::map<std::string,std::string> &subRoutineCodes, std::vector<std::string> &otherRoutineNames, bool derivatives, bool potential) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	static std::string getCLProgram(double deflectionScale, double potentialScale, const std::vector<std::string> &subRoutineNames, int recursionLevel, int maxRecursion, bool derivatives, bool potential);
	void findCLSubroutines(double deflectionScale, double potentialScale, std::map<std::string,std::string> &subRoutineCodes, std::vector<std::string> &otherRoutineNames, int recursionLevel, int &maxRecursionLevel, bool derivatives, bool potential) const;

	std::vector<std::shared_ptr<GravitationalLens>> m_lenses;
	std::vector<Vector2D<double> > m_positions;
	std::vector<double> m_angles;
	std::vector<double> m_origAngles;
	std::vector<double> m_factors;
};
	
} // end namespace

#endif // GRALE_COMPOSITELENS_H
