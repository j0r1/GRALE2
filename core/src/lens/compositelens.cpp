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
#include "compositelens.h"
#include "multiplanecontainer.h"
#include "constants.h"
#include <serut/dummyserializer.h>
#include <serut/memoryserializer.h>
#include <string.h>
#include <stdio.h>
#include <stdexcept>
#include <iostream>

namespace grale
{

class CompositeLensParams::LensInfo
{
public:
	LensInfo(double factor, Vector2D<double> position, double angle, std::vector<uint8_t> &pData);
	~LensInfo();
	double getFactor() const;
	double getAngle() const;
	Vector2D<double> getPosition() const;
	const uint8_t *getData() const;
	int getDataLength() const;
private:
	double m_factor;
	double m_angle;
	Vector2D<double> m_position;
	std::vector<uint8_t> m_pData;
};

CompositeLensParams::LensInfo::LensInfo(double factor, Vector2D<double> position, double angle, std::vector<uint8_t> &pData)
{
	m_factor = factor;
	m_position = position;
	m_angle = angle;
	swap(m_pData, pData);
}
		
CompositeLensParams::LensInfo::~LensInfo()
{
}

double CompositeLensParams::LensInfo::getFactor() const
{ 
	return m_factor; 
}

double CompositeLensParams::LensInfo::getAngle() const
{ 
	return m_angle; 
}

Vector2D<double> CompositeLensParams::LensInfo::getPosition() const
{ 
	return m_position; 
}

const uint8_t *CompositeLensParams::LensInfo::getData() const
{ 
	return m_pData.data(); 
}

int CompositeLensParams::LensInfo::getDataLength() const
{ 
	return (int)m_pData.size(); 
}

CompositeLensParams::CompositeLensParams()
{
}

CompositeLensParams::~CompositeLensParams()
{
}

bool CompositeLensParams::addLens(double factor, Vector2D<double> position, double angle, const GravitationalLens &lens)
{
	if (dynamic_cast<const MultiPlaneContainer*>(&lens) != 0)
	{
		setErrorString("A MultiPlaneContainer lens cannot be used inside a CompositeLens");
		return false;
	}

	serut::DummySerializer dumSer;

	if (!lens.write(dumSer))
	{
		setErrorString(std::string("Couldn't store lens info: ") + dumSer.getErrorString());
		return false;
	}

	std::vector<uint8_t> pData(dumSer.getBytesWritten());
	serut::MemorySerializer memSer(0, 0, pData.data(), (int)pData.size());
	
	lens.write(memSer);
	
	m_lensInfo.push_back(std::make_shared<LensInfo>(factor, position, angle, pData));
	
	return true;
}

bool CompositeLensParams::write(serut::SerializationInterface &si) const
{
	if (!si.writeInt32(m_lensInfo.size()))
	{
		setErrorString(std::string("Couldn't write number of lenses: ") + si.getErrorString());
		return false;
	}

	for (auto it = m_lensInfo.begin() ; it != m_lensInfo.end() ; ++it)
	{
		const LensInfo *pLensInfo = it->get();

		if (!si.writeDouble(pLensInfo->getFactor()))
		{
			setErrorString(std::string("Couldn't write scale factor: ") + si.getErrorString());
			return false;
		}

		if (!si.writeDouble(pLensInfo->getAngle()))
		{
			setErrorString(std::string("Couldn't write rotation angle: ") + si.getErrorString());
			return false;
		}

		if (!si.writeDoubles(pLensInfo->getPosition().getComponents(), 2))
		{
			setErrorString(std::string("Couldn't write position: ") + si.getErrorString());
			return false;
		}

		if (!si.writeInt32(pLensInfo->getDataLength()))
		{
			setErrorString(std::string("Couldn't write data length: ") + si.getErrorString());
			return false;
		}

		if (!si.writeBytes(pLensInfo->getData(), pLensInfo->getDataLength()))
		{
			setErrorString(std::string("Couldn't write lens data: ") + si.getErrorString());
			return false;
		}
	}

	return true;
}

bool CompositeLensParams::read(serut::SerializationInterface &si)
{
	m_lensInfo.clear();

	int32_t numLenses = 0;
			
	if (!si.readInt32(&numLenses))
	{
		setErrorString(std::string("Couldn't read number of lenses: ") + si.getErrorString());
		return false;
	}	

	for (int32_t i = 0 ; i < numLenses ; i++)
	{
		double factor, angle;
		Vector2D<double> position;
		int32_t dataLength;

		if (!si.readDouble(&factor))
		{
			setErrorString(std::string("Couldn't read scale factor: ") + si.getErrorString());
			return false;
		}

		if (!si.readDouble(&angle))
		{
			setErrorString(std::string("Couldn't read rotation angle: ") + si.getErrorString());
			return false;
		}
		
		if (!si.readDoubles(position.getComponents(), 2))
		{
			setErrorString(std::string("Couldn't read scale position: ") + si.getErrorString());
			return false;
		}

		if (!si.readInt32(&dataLength))
		{
			setErrorString(std::string("Couldn't read data length: ") + si.getErrorString());
			return false;
		}

		std::vector<uint8_t> pData(dataLength);

		if (!si.readBytes(pData.data(), (int)pData.size()))
		{
			setErrorString(std::string("Couldn't read lens data: ") + si.getErrorString());
			return false;
		}

		m_lensInfo.push_back(std::make_shared<LensInfo>(factor, position, angle, pData));
	}

	return true;
}

std::unique_ptr<GravitationalLensParams> CompositeLensParams::createCopy() const
{
	std::unique_ptr<CompositeLensParams> pParams = std::make_unique<CompositeLensParams>();

	for (auto it = m_lensInfo.begin() ; it != m_lensInfo.end() ; ++it)
	{
		const LensInfo *pLensInfo = it->get();
		std::vector<uint8_t> pData(pLensInfo->getDataLength());

		memcpy(pData.data(), pLensInfo->getData(), pLensInfo->getDataLength());
		pParams->m_lensInfo.push_back(std::make_shared<LensInfo>(pLensInfo->getFactor(), pLensInfo->getPosition(), pLensInfo->getAngle(),
		                                           pData));
	}

	return pParams;
}

CompositeLens *CompositeLens::cast(GravitationalLens *pLens)
{
	return dynamic_cast<CompositeLens*>(pLens);
}

const CompositeLens *CompositeLens::cast(const GravitationalLens *pLens)
{
	return dynamic_cast<const CompositeLens*>(pLens);
}

CompositeLens::CompositeLens() : GravitationalLens(GravitationalLens::Composite)
{
}

CompositeLens::~CompositeLens()
{
}

bool CompositeLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const CompositeLensParams *pParams = dynamic_cast<const CompositeLensParams *>(pLensParams);

	if (!pParams)
	{
		setErrorString("Parameters are not of type 'CompositeLensParams'");
		return false;
	}

	auto lensInfo = pParams->getLensInfo();
	
	for (auto it = lensInfo.begin() ; it != lensInfo.end() ; ++it)
	{
		const CompositeLensParams::LensInfo *pLensInfo = it->get();
		serut::MemorySerializer mSer(pLensInfo->getData(), pLensInfo->getDataLength(), 0, 0);
		std::unique_ptr<GravitationalLens> pLens;
		std::string errStr;
		
		if (!GravitationalLens::read(mSer, pLens, errStr))
		{
			setErrorString(std::string("Couldn't initialize a lens component: ") + errStr);
			return false;
		}

		m_lenses.push_back(move(pLens));
		m_positions.push_back(pLensInfo->getPosition());
		m_factors.push_back(pLensInfo->getFactor());
		m_origAngles.push_back(pLensInfo->getAngle());
		m_angles.push_back(pLensInfo->getAngle()*CONST_PI/180.0);
	}

	return true;
}

bool CompositeLens::getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const
{
	Vector2D<double> alphaSum;
	
	for (int i = 0 ; i < (int)m_lenses.size() ; i++)
	{
		double rotationangle = m_angles[i];
		Vector2D<double> theta0 = theta - m_positions[i];
		
		theta0 = Vector2D<double>(theta0.getX()*COS(rotationangle)+theta0.getY()*SIN(rotationangle),
			                 -theta0.getX()*SIN(rotationangle)+theta0.getY()*COS(rotationangle));

		Vector2D<double> alpha0;
		
		if (!m_lenses[i]->getAlphaVector(theta0, &alpha0))
		{
			setErrorString(m_lenses[i]->getErrorString());
			return false;
		}

		Vector2D<double> alpha = alpha0 * m_factors[i];
		
		alpha = Vector2D<double>(alpha.getX()*COS(rotationangle)-alpha.getY()*SIN(rotationangle),
	                                 alpha.getX()*SIN(rotationangle)+alpha.getY()*COS(rotationangle));

		alphaSum += alpha;
	}
	
	*pAlpha = alphaSum;

	return true;
}

double CompositeLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	double sum = 0;

	for (int i = 0 ; i < (int)m_lenses.size() ; i++)
	{
		double rotationangle = m_angles[i];
		Vector2D<double> theta0 = theta - m_positions[i];

		theta0 = Vector2D<double>(theta0.getX()*COS(rotationangle)+theta0.getY()*SIN(rotationangle),
			                 -theta0.getX()*SIN(rotationangle)+theta0.getY()*COS(rotationangle));

		sum += m_factors[i] * m_lenses[i]->getSurfaceMassDensity(theta0);
	}
		
	return sum;
}

bool CompositeLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *potentialval) const
{
	double sum = 0;

	for (int i = 0 ; i < (int)m_lenses.size() ; i++)
	{
		double rotationangle = m_angles[i];
		Vector2D<double> theta0 = theta - m_positions[i];
		double phi = 0;

		theta0 = Vector2D<double>(theta0.getX()*COS(rotationangle)+theta0.getY()*SIN(rotationangle),
			                 -theta0.getX()*SIN(rotationangle)+theta0.getY()*COS(rotationangle));

		if (!m_lenses[i]->getProjectedPotential(D_s, D_ds, theta0, &phi))
		{
			setErrorString(std::string("Sublens error: ") + m_lenses[i]->getErrorString());
			return false;
		}
		
		sum += m_factors[i] * phi;
	}
		
	*potentialval = sum;
	
	return true;
}

bool CompositeLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	double axxSum = 0;
	double ayySum = 0;
	double axySum = 0;

	for (int i = 0 ; i < (int)m_lenses.size() ; i++)
	{
		Vector2D<double> theta0 = theta - m_positions[i];
		double rotationangle = m_angles[i];
		double axx0, ayy0, axy0;
		double ca = COS(rotationangle);
		double sa = SIN(rotationangle);
		double ca2 = ca*ca;
		double sa2 = sa*sa;
		double csa = ca*sa;
	
		theta0 = Vector2D<double>(theta0.getX()*ca+theta0.getY()*sa,
			                 -theta0.getX()*sa+theta0.getY()*ca);

		if (!m_lenses[i]->getAlphaVectorDerivatives(theta0, axx0, ayy0, axy0))
		{
			setErrorString(m_lenses[i]->getErrorString());
			return false;
		}

		axx0 *= m_factors[i];
		axy0 *= m_factors[i];
		ayy0 *= m_factors[i];
		
		double axx1 = axx0*ca2-2.0*axy0*csa+ayy0*sa2;
		double ayy1 = axx0*sa2+2.0*axy0*csa+ayy0*ca2;
		double axy1 = (axx0-ayy0)*csa+axy0*(2.0*ca2-1.0);

		axxSum += axx1;
		ayySum += ayy1;
		axySum += axy1;
	}
	
	axx = axxSum;
	ayy = ayySum;
	axy = axySum;
	
	return true;
}

bool CompositeLens::getAlphaVectorSecondDerivatives(Vector2D<double> theta, double &axxx, double &ayyy, double &axxy, double &ayyx) const
{
	double axxxSum = 0;
	double ayyySum = 0;
	double axxySum = 0;
	double ayyxSum = 0;

	for (size_t i = 0 ; i < m_lenses.size() ; i++)
	{
		Vector2D<double> theta0 = theta - m_positions[i];
		double rotationangle = m_angles[i];
		double axxx0, ayyy0, axxy0, ayyx0;
		double ca = COS(rotationangle);
		double sa = SIN(rotationangle);
		double factor = m_factors[i];

		theta0 = Vector2D<double>(theta0.getX()*ca+theta0.getY()*sa,
			                      -theta0.getX()*sa+theta0.getY()*ca);

		if (!m_lenses[i]->getAlphaVectorSecondDerivatives(theta0, axxx0, ayyy0, axxy0, ayyx0))
		{
			setErrorString(m_lenses[i]->getErrorString());
			return false;
		}

		axxx0 *= factor;
		axxy0 *= factor;
		ayyx0 *= factor;
		ayyy0 *= factor;

		double R[2][2] = { { ca, -sa}, { sa, ca } };
		auto A = [axxx0, ayyy0, axxy0, ayyx0](size_t a, size_t b, size_t c)
		{
			assert(a == 0 || a == 1);
			assert(b == 0 || b == 1);
			assert(c == 0 || c == 1);
			size_t num = a | (b<<1) | (c<<2);
			switch (num)
			{
			case 0: // xxx
				return axxx0;
			case 1: // yxx
			case 2: // xyx
			case 4: // xxy
				return axxy0;
			case 3: // yyx
			case 5: // yxy
			case 6: // yyx
				return ayyx0;
			case 7: // yyy
				return ayyy0;
			}
			throw std::runtime_error("Internal error, num = " + std::to_string(num));
		};
		
		auto sumIt = [&R,&A](size_t u, size_t v, size_t w)
		{
			double s = 0;
			for (size_t a = 0 ; a < 2 ; a++)
				for (size_t b = 0 ; b < 2 ; b++)
					for(size_t c = 0 ; c < 2 ; c++)
						s += R[u][a]*R[v][b]*R[w][c]*A(a,b,c);
			return s;
		};

		double axxx1 = sumIt(0, 0, 0);
		double ayyy1 = sumIt(1, 1, 1);
		double axxy1 = sumIt(0, 0, 1);
		double ayyx1 = sumIt(1, 1, 0);
		
		axxxSum += axxx1;
		ayyySum += ayyy1;
		axxySum += axxy1;
		ayyxSum += ayyx1;
	}

	axxx = axxxSum;
	ayyy = ayyySum;
	axxy = axxySum;
	ayyx = ayyxSum;

	return true;
}

void CompositeLens::setDerivativeAngularDistanceScale(double dist)
{
	GravitationalLens::setDerivativeAngularDistanceScale(dist);

	for (int i = 0 ; i < (int)m_lenses.size() ; i++)
		m_lenses[i]->setDerivativeAngularDistanceScale(dist);
}

bool CompositeLens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	double deflectionScale = 0;
	double potentialScale = 0;

	for (int i = 0 ; i < (int)m_lenses.size() ; i++)
	{
		double d, p;

		if (!m_lenses[i]->getSuggestedScales(&d, &p))
		{
			setErrorString(std::string("Can't get suggested scale for a sublens: ") + m_lenses[i]->getErrorString());
			return false;
		}

		deflectionScale += d*d*m_factors[i];
		potentialScale += p;
	}

	deflectionScale /= (double)m_lenses.size();
	deflectionScale = SQRT(deflectionScale);

	potentialScale /= (double)m_lenses.size();

	*pDeflectionScale = deflectionScale;
	*pPotentialScale = potentialScale;

	return true;
}

bool CompositeLens::getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const
{
	int intCount = 1; // al least one for the number of sublenses
	int floatCount = 0;

	for (int i = 0 ; i < (int)m_lenses.size() ; i++)
	{
		int numSubInt = 0;
		int numSubFloat = 0;

		if (!m_lenses[i]->getCLParameterCounts(&numSubInt, &numSubFloat))
		{
			setErrorString(m_lenses[i]->getErrorString());
			return false;
		}

		intCount += 3; // lens type, number of int parameters and float parameters for this lens
		intCount += numSubInt;

		floatCount += 4; // position, angle and scale factor
		floatCount += numSubFloat;
	}

	*pNumIntParams = intCount;
	*pNumFloatParams = floatCount;

	return true;
}

bool CompositeLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	pIntParams[0] = m_lenses.size();
	
	int intOffset = 1;
	int floatOffset = 0;

	for (int i = 0 ; i < (int)m_lenses.size() ; i++)
	{
		int numSubInt = 0;
		int numSubFloat = 0;

		if (!m_lenses[i]->getCLParameterCounts(&numSubInt, &numSubFloat))
		{
			setErrorString(m_lenses[i]->getErrorString());
			return false;
		}

		int lensNumber = m_lenses[i]->getLensType();

		pIntParams[intOffset++] = lensNumber;
		pIntParams[intOffset++] = numSubInt;
		pIntParams[intOffset++] = numSubFloat;

		pFloatParams[floatOffset++] = m_positions[i].getX()/deflectionScale;
		pFloatParams[floatOffset++] = m_positions[i].getY()/deflectionScale;
		pFloatParams[floatOffset++] = (float)(m_angles[i]);
		pFloatParams[floatOffset++] = (float)(m_factors[i]);

		if (!m_lenses[i]->getCLParameters(deflectionScale, potentialScale, pIntParams+intOffset, pFloatParams+floatOffset))
		{
			setErrorString(m_lenses[i]->getErrorString());
			return false;
		}

		intOffset += numSubInt;
		floatOffset += numSubFloat;
	}

	return true;
}

std::unique_ptr<GravitationalLensParams> CompositeLens::createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const
{
	auto params = std::make_unique<CompositeLensParams>();
	size_t floatOffset = 0;

	for (int i = 0 ; i < (int)m_lenses.size() ; i++)
	{
		int numSubInt = 0;
		int numSubFloat = 0;

		if (!m_lenses[i]->getCLParameterCounts(&numSubInt, &numSubFloat))
		{
			setErrorString(m_lenses[i]->getErrorString());
			return nullptr;
		}

		double x = ((double)pFloatParams[floatOffset++])*deflectionScale;
		double y = ((double)pFloatParams[floatOffset++])*deflectionScale;
		double angle = ((double)pFloatParams[floatOffset++])*180.0/CONST_PI;
		double factor = (double)pFloatParams[floatOffset++];

		auto newLens = m_lenses[i]->createLensFromCLFloatParams(deflectionScale, potentialScale, pFloatParams+floatOffset);
		if (!newLens)
		{
			setErrorString("Unable to create a sublens from OpenCL float params: " + m_lenses[i]->getErrorString());
			return nullptr;
		}

		floatOffset += numSubFloat;

		if (!params->addLens(factor, {x, y}, angle, *newLens))
		{
			setErrorString("Unable to add a lens to the composite lens parameters: " + params->getErrorString());
			return nullptr;
		}
	}
	return params;
}

std::vector<CLFloatParamInfo> CompositeLens::getCLAdjustableFloatingPointParameterInfo() const
{
	std::vector<CLFloatParamInfo> allParamInfo;
	size_t floatOffset = 0;

	for (int i = 0 ; i < (int)m_lenses.size() ; i++)
	{
		int numSubInt = 0;
		int numSubFloat = 0;

		if (!m_lenses[i]->getCLParameterCounts(&numSubInt, &numSubFloat))
		{
			std::cerr << "Warning: couldn't get CL parameter counts" << m_lenses[i]->getErrorString() << std::endl;
			return {};
		}

		//std::string name_prefix = "comppart_" + std::to_string(i) + ",";
		allParamInfo.push_back({.name = "x_" + std::to_string(i) + "_scaled", .offset = floatOffset++});
		allParamInfo.push_back({.name = "y_" + std::to_string(i) + "_scaled", .offset = floatOffset++});
		allParamInfo.push_back({.name = "angle_" + std::to_string(i), .offset = floatOffset++});
		allParamInfo.push_back({.name = "factor_" + std::to_string(i), .offset = floatOffset++, .hardMin = 0 }); // Allow negatives?

		auto paramInfo = m_lenses[i]->getCLAdjustableFloatingPointParameterInfo();
		for (auto s : paramInfo)
		{
			s.name = "lens_" + std::to_string(i) + "," + s.name;
			s.offset += floatOffset;
			allParamInfo.push_back(s);
		}

		floatOffset += numSubFloat;
	}
	return allParamInfo;
}

std::string CompositeLens::getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives, bool potential) const
{
	std::vector<std::string> otherRoutineNames;
	std::map<std::string, std::string> subCodes;
	int maxRecursionCount = findCLSubroutines(deflectionScale, potentialScale, subCodes, otherRoutineNames, derivatives, potential);

	std::string prog;
	for (auto &kv : subCodes)
		prog += kv.second;

	prog += getCLProgram(deflectionScale, potentialScale, subRoutineName, otherRoutineNames, maxRecursionCount, derivatives, potential);

	return prog;
}

std::string CompositeLens::getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, const std::vector<std::string> &otherRoutineNames, int maxRecursionCount, 
		                                bool derivatives, bool potential)
{
	if (otherRoutineNames.size() != MaxLensType)
		return "ERROR: otherRoutineNames must be of length " + std::to_string(MaxLensType);

	subRoutineName = "clCompositeLensProgram";

	std::string prog;
	for (int i = maxRecursionCount ; i >= 0 ; i--)
		prog += getCLProgram(deflectionScale, potentialScale, otherRoutineNames, i, maxRecursionCount, derivatives, potential);
	return prog;
}

int CompositeLens::findCLSubroutines(double deflectionScale, double potentialScale, std::map<std::string,std::string> &subRoutineCodes, std::vector<std::string> &otherRoutineNames, bool derivatives, bool potential) const
{
	int maxRecursionCount = 0;
	otherRoutineNames.resize(MaxLensType);
	otherRoutineNames[getLensType()] = "clCompositeLensProgram";
	findCLSubroutines(deflectionScale, potentialScale, subRoutineCodes, otherRoutineNames, 0, maxRecursionCount, derivatives, potential);
	return maxRecursionCount;
}

void CompositeLens::findCLSubroutines(double deflectionScale, double potentialScale, std::map<std::string,std::string> &subCodes, std::vector<std::string> &otherRoutineNames, int recursionLevel, int &maxRecursionLevel, bool derivatives, bool potential) const
{
	if (recursionLevel > maxRecursionLevel)
		maxRecursionLevel = recursionLevel;

	for (int i = 0 ; i < (int)m_lenses.size() ; i++)
	{
		int lensNumber = m_lenses[i]->getLensType();

		if (lensNumber == getLensType()) // another composite lens
		{
			const CompositeLens *pNextCompLens = (const CompositeLens *)(m_lenses[i].get());

			pNextCompLens->findCLSubroutines(deflectionScale, potentialScale, subCodes, otherRoutineNames, recursionLevel+1, maxRecursionLevel, derivatives, potential);
		}
		else
		{
			if (otherRoutineNames[lensNumber].length() == 0)
			{
				std::string subName;
				std::string prog = m_lenses[i]->getCLProgram(deflectionScale, potentialScale, subName, derivatives, potential);

				otherRoutineNames[lensNumber] = subName;
				subCodes[subName] = prog;
			}
		}
	}
}

std::string CompositeLens::getCLProgram(double deflectionScale, double potentialScale, const std::vector<std::string> &subRoutineNames, 
                                        int recursionLevel, int maxRecursion,
										bool derivatives, bool potential)
{
	std::string prog;

	prog += "LensQuantities clCompositeLensProgram";
	if (recursionLevel != 0)
	{
		char str[256];

		sprintf(str, "%d", recursionLevel);

		prog += std::string(str);
	}
	prog += "(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)\n";
	prog += "{\n";
	prog += "	LensQuantities zeroQuant = { 0 };\n";
	prog += "	LensQuantities r = zeroQuant ;\n";
	prog += "	int numSubLenses = pIntParams[0];\n";
	prog += "	int intOffset = 1;\n";
	prog += "	int floatOffset = 0;\n";
	prog += "\n";
	prog += "	for (int i = 0 ; i < numSubLenses ; i++)\n";
	prog += "	{\n";
	prog += "		int lensNumber = pIntParams[intOffset];\n";
	prog += "		int numIntParams = pIntParams[intOffset+1];\n";
	prog += "		int numFloatParams = pIntParams[intOffset+2];\n";
	prog += "\n";
	prog += "		intOffset += 3;\n";
	prog += "\n";
	prog += "		float xPos = pFloatParams[floatOffset];\n";
	prog += "		float yPos = pFloatParams[floatOffset+1];\n";
	prog += "		float angle = pFloatParams[floatOffset+2];\n";
	prog += "		float scaleFactor = pFloatParams[floatOffset+3];\n";
	prog += "\n";
	prog += "		floatOffset += 4;\n";
	prog += "\n";
	prog += "		float dx0 = coord.x-xPos;\n";
	prog += "		float dy0 = coord.y-yPos;\n";
	prog += "		float ca = cos(angle);\n";
	prog += "		float sa = sin(angle);\n";
	prog += "		float ca2 = ca*ca;\n";
	prog += "		float sa2 = sa*sa;\n";
	prog += "		float csa = ca*sa;\n";
	prog += "\n";
	prog += "		float2 newCoord = { dx0*ca+dy0*sa, -dx0*sa+dy0*ca };\n";
	prog += "\n";
	prog += "		LensQuantities s;\n";
	prog += "\n";
	prog += "		if (lensNumber == -1)\n";
	prog += "				return zeroQuant;\n";
	for (int i = 0 ; i < (int)subRoutineNames.size() ; i++)
	{
		if (subRoutineNames[i].length() > 0)
		{
			char str[256];

			sprintf(str, "%d", i);

			if (subRoutineNames[i] != std::string("clCompositeLensProgram"))
			{
				prog += "		else if (lensNumber == " + std::string(str) +")\n";
				prog += "			s = " + subRoutineNames[i] + "(newCoord, &(pIntParams[intOffset]), &(pFloatParams[floatOffset]));\n";
			}
			else
			{
				if (recursionLevel < maxRecursion)
				{
					prog += "		else if (lensNumber == " + std::string(str) +")\n";

					sprintf(str, "%d", recursionLevel+1);

					prog += "			s = " + subRoutineNames[i] + std::string(str) + "(newCoord, &(pIntParams[intOffset]), &(pFloatParams[floatOffset]));\n";
				}
			}
		}
	}
	prog += "\n";
	prog += "		intOffset += numIntParams;\n";
	prog += "		floatOffset += numFloatParams;\n";
	prog += "\n";
	prog += "		float ax = s.alphaX * scaleFactor;\n";
	prog += "		float ay = s.alphaY * scaleFactor;\n";
	prog += "\n";
	if (derivatives)
	{
		prog += "		float axx0 = s.axx * scaleFactor;\n";
		prog += "		float axy0 = s.axy * scaleFactor;\n";
		prog += "		float ayy0 = s.ayy * scaleFactor;\n";
	}
	prog += "\n";
	prog += "		r.alphaX += (ax*ca-ay*sa);\n";
	prog += "		r.alphaY += (ax*sa+ay*ca);\n";
	if (potential)
		prog += "		r.potential += s.potential * scaleFactor;\n";
	if (derivatives)
	{
		prog += "		r.axx += axx0*ca2-2.0*axy0*csa+ayy0*sa2;\n";
		prog += "		r.ayy += axx0*sa2+2.0*axy0*csa+ayy0*ca2;\n";
		prog += "		r.axy += (axx0-ayy0)*csa+axy0*(2.0*ca2-1.0);\n";
	}
	prog += "	}\n";
	prog += "\n";
	prog += "	return r;\n";
	prog += "}\n";

	return prog;
}

} // end namespace

