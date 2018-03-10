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
#include "constants.h"
#include <serut/dummyserializer.h>
#include <serut/memoryserializer.h>
#include <string.h>
#include <stdio.h>

#include "debugnew.h"

namespace grale
{

class CompositeLensParams::LensInfo
{
public:
	LensInfo(double factor, Vector2D<double> position, double angle, uint8_t *pData, int dataLength);
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
	int m_dataLength;
	uint8_t *m_pData;
};

CompositeLensParams::LensInfo::LensInfo(double factor, Vector2D<double> position, double angle, uint8_t *pData, int dataLength)
{
	m_factor = factor;
	m_position = position;
	m_angle = angle;
	m_pData = pData;
	m_dataLength = dataLength;
}
		
CompositeLensParams::LensInfo::~LensInfo()
{
	delete [] m_pData;
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
	return m_pData; 
}

int CompositeLensParams::LensInfo::getDataLength() const
{ 
	return m_dataLength; 
}

CompositeLensParams::CompositeLensParams()
{
}

CompositeLensParams::~CompositeLensParams()
{
	for (size_t i = 0 ; i < m_lensInfo.size() ; i++)
		delete m_lensInfo[i];
	m_lensInfo.clear();
}

bool CompositeLensParams::addLens(double factor, Vector2D<double> position, double angle, const GravitationalLens &lens)
{
	serut::DummySerializer dumSer;

	if (!lens.write(dumSer))
	{
		setErrorString(std::string("Couldn't store lens info: ") + dumSer.getErrorString());
		return false;
	}

	uint8_t *pData = new uint8_t[dumSer.getBytesWritten()];
	serut::MemorySerializer memSer(0, 0, pData, dumSer.getBytesWritten());
	
	lens.write(memSer);
	
	m_lensInfo.push_back(new LensInfo(factor, position, angle, pData, dumSer.getBytesWritten()));
	
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
		const LensInfo *pLensInfo = (*it);

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
	for (size_t i = 0 ; i < m_lensInfo.size() ; i++)
		delete m_lensInfo[i];
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

		uint8_t *pData = new uint8_t[dataLength];

		if (!si.readBytes(pData, dataLength))
		{
			delete [] pData;
			setErrorString(std::string("Couldn't read lens data: ") + si.getErrorString());
			return false;
		}

		m_lensInfo.push_back(new LensInfo(factor, position, angle, pData, dataLength));
	}

	return true;
}

GravitationalLensParams *CompositeLensParams::createCopy() const
{
	CompositeLensParams *pParams = new CompositeLensParams();

	for (auto it = m_lensInfo.begin() ; it != m_lensInfo.end() ; ++it)
	{
		const LensInfo *pLensInfo = (*it);
		uint8_t *pData = new uint8_t[pLensInfo->getDataLength()];

		memcpy(pData, pLensInfo->getData(), pLensInfo->getDataLength());
		pParams->m_lensInfo.push_back(new LensInfo(pLensInfo->getFactor(), pLensInfo->getPosition(), pLensInfo->getAngle(),
		                                           pData, pLensInfo->getDataLength()));
	}

	return pParams;
}

CompositeLens::CompositeLens() : GravitationalLens(GravitationalLens::Composite)
{
}

CompositeLens::~CompositeLens()
{
	for (int i = 0 ; i < m_lenses.size() ; i++)
		delete m_lenses[i];
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
		const CompositeLensParams::LensInfo *pLensInfo = (*it);
		serut::MemorySerializer mSer(pLensInfo->getData(), pLensInfo->getDataLength(), 0, 0);
		GravitationalLens *pLens = 0;
		std::string errStr;
		
		if (!GravitationalLens::read(mSer, &pLens, errStr))
		{
			setErrorString(std::string("Couldn't initialize a lens component: ") + errStr);
			return false;
		}

		m_lenses.push_back(pLens);
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
	
	for (int i = 0 ; i < m_lenses.size() ; i++)
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

	for (int i = 0 ; i < m_lenses.size() ; i++)
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

	for (int i = 0 ; i < m_lenses.size() ; i++)
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

	for (int i = 0 ; i < m_lenses.size() ; i++)
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

void CompositeLens::setDerivativeAngularDistanceScale(double dist)
{
	GravitationalLens::setDerivativeAngularDistanceScale(dist);

	for (int i = 0 ; i < m_lenses.size() ; i++)
		m_lenses[i]->setDerivativeAngularDistanceScale(dist);
}

bool CompositeLens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	double deflectionScale = 0;
	double potentialScale = 0;

	for (int i = 0 ; i < m_lenses.size() ; i++)
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

	for (int i = 0 ; i < m_lenses.size() ; i++)
	{
		int numSubInt = 0;
		int numSubFloat = 0;

		if (!m_lenses[i]->getCLParameterCounts(&numSubInt, &numSubFloat))
		{
			setErrorString(m_lenses[i]->getErrorString());
			return false;
		}

		intCount += 4; // sublenses, lens type, number of int parameters and float parameters for this lens
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

	for (int i = 0 ; i < m_lenses.size() ; i++)
	{
		int numSubInt = 0;
		int numSubFloat = 0;

		if (!m_lenses[i]->getCLParameterCounts(&numSubInt, &numSubFloat))
		{
			setErrorString(m_lenses[i]->getErrorString());
			return false;
		}

		int lensNumber = m_lenses[i]->getLensType();

		pIntParams[intOffset++] = m_lenses[i]->getCLSubLenses();
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

std::string CompositeLens::getCLProgram(std::string &subRoutineName) const
{
	int numOtherSubRoutines = GravitationalLens::MaxLensType;
	int maxRecursionCount = 0;
	std::string prog;
	std::vector<std::string> otherRoutineNames(numOtherSubRoutines);

	otherRoutineNames[getLensType()] = "clCompositeLensProgram";
	subRoutineName = otherRoutineNames[getLensType()];

	findCLSubroutines(prog, otherRoutineNames, 0, maxRecursionCount);

	for (int i = maxRecursionCount ; i >= 0 ; i--)
		prog += getCLProgram(otherRoutineNames, i, maxRecursionCount);

	return prog;
}

void CompositeLens::findCLSubroutines(std::string &prog, std::vector<std::string> &otherRoutineNames, int recursionLevel, int &maxRecursionLevel) const
{
	if (recursionLevel > maxRecursionLevel)
		maxRecursionLevel = recursionLevel;

	for (int i = 0 ; i < m_lenses.size() ; i++)
	{
		int lensNumber = m_lenses[i]->getLensType();

		if (lensNumber == getLensType()) // another composite lens
		{
			const CompositeLens *pNextCompLens = (const CompositeLens *)(m_lenses[i]);

			pNextCompLens->findCLSubroutines(prog, otherRoutineNames, recursionLevel+1, maxRecursionLevel);
		}
		else
		{
			if (otherRoutineNames[lensNumber].length() == 0)
			{
				std::string subName;

				prog += m_lenses[i]->getCLProgram(subName);

				otherRoutineNames[lensNumber] = subName;
			}
		}
	}
}

std::string CompositeLens::getCLProgram(const std::vector<std::string> &subRoutineNames, 
                                        int recursionLevel, int maxRecursion)
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
	prog += "	LensQuantities zeroQuant = { 0, 0, 0, 0, 0, 0};\n";
	prog += "	LensQuantities r = zeroQuant ;\n";
	prog += "	int numSubLenses = pIntParams[0];\n";
	prog += "	int intOffset = 1;\n";
	prog += "	int floatOffset = 0;\n";
	prog += "\n";
	prog += "	for (int i = 0 ; i < numSubLenses ; i++)\n";
	prog += "	{\n";
	prog += "		int numCLSublenses = pIntParams[intOffset];\n";
	prog += "		int lensNumber = pIntParams[intOffset+1];\n";
	prog += "		int numIntParams = pIntParams[intOffset+2];\n";
	prog += "		int numFloatParams = pIntParams[intOffset+3];\n";
	prog += "\n";
	prog += "		intOffset += 4;\n";
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
	for (int i = 0 ; i < subRoutineNames.size() ; i++)
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
	prog += "		float axx0 = s.axx * scaleFactor;\n";
	prog += "		float axy0 = s.axy * scaleFactor;\n";
	prog += "		float ayy0 = s.ayy * scaleFactor;\n";
	prog += "\n";
	prog += "		r.alphaX += (ax*ca-ay*sa);\n";
	prog += "		r.alphaY += (ax*sa+ay*ca);\n";
	prog += "		r.potential += s.potential * scaleFactor;\n";
	prog += "		r.axx += axx0*ca2-2.0*axy0*csa+ayy0*sa2;\n";
	prog += "		r.ayy += axx0*sa2+2.0*axy0*csa+ayy0*ca2;\n";
	prog += "		r.axy += (axx0-ayy0)*csa+axy0*(2.0*ca2-1.0);\n";
	prog += "	}\n";
	prog += "\n";
	prog += "	return r;\n";
	prog += "}\n";

	return prog;
}

int CompositeLens::getCLSubLenses() const
{
	int sum = 0;

	for (int i = 0 ; i < m_lenses.size() ; i++)
	{
		int num = m_lenses[i]->getCLSubLenses();

		if (num < 1) // can't handle the num == 0 case yet
			return -1;

		sum += num;
	}
	return sum;
}

} // end namespace

