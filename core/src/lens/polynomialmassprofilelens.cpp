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
#include "polynomialmassprofilelens.h"
#include "constants.h"
#include <stdio.h>

namespace grale
{

PolynomialPart::PolynomialPart()
{
	m_xOffset = 0;
	m_yOffset = 0;
	m_xScale = 0;
	m_yScale = 0;
	m_xEnd = 0;
}

PolynomialPart::PolynomialPart(double xOffset, double yOffset, double xScale, double yScale, double xEnd, const std::vector<double> &coefficients)
{
	m_xOffset = xOffset;
	m_yOffset = yOffset;
	m_xScale = xScale;
	m_yScale = yScale;
	m_xEnd = xEnd;
	m_coefficients = coefficients;

	calculateOffsetCoefficients();
}

PolynomialPart::PolynomialPart(const PolynomialPart &polynomialPart)
{
	m_xOffset = polynomialPart.m_xOffset;
	m_yOffset = polynomialPart.m_yOffset;
	m_xScale = polynomialPart.m_xScale;
	m_yScale = polynomialPart.m_yScale;
	m_xEnd = polynomialPart.m_xEnd;
	m_coefficients = polynomialPart.m_coefficients;
	m_offsetCoefficients = polynomialPart.m_offsetCoefficients;
}

PolynomialPart &PolynomialPart::operator=(const PolynomialPart &polynomialPart)
{
	m_xOffset = polynomialPart.m_xOffset;
	m_yOffset = polynomialPart.m_yOffset;
	m_xScale = polynomialPart.m_xScale;
	m_yScale = polynomialPart.m_yScale;
	m_xEnd = polynomialPart.m_xEnd;
	m_coefficients = polynomialPart.m_coefficients;
	m_offsetCoefficients = polynomialPart.m_offsetCoefficients;

	return *this;
}

PolynomialPart::~PolynomialPart()
{
}

void PolynomialPart::calculateOffsetCoefficients()
{
	if (m_coefficients.size() == 0)
		return;
	
	int M = m_coefficients.size()-1; // degree of the polynomial
	double b = m_xOffset/m_xScale;
	
	if (b == 0)
		m_offsetCoefficients = m_coefficients;
	
	// we'll calculate and store the necessary binomial coefficients here
	
	std::vector<std::vector<int> > binomialCoefficients; 
	
	binomialCoefficients.resize(M+1);
	binomialCoefficients[0].resize(1);
	binomialCoefficients[0][0] = 1;
	
	for (int i = 1 ; i <= M ; i++)
	{
		binomialCoefficients[i].resize(i+1);
		
		for (int k = 0 ; k <= i ; k++)
		{
			int value1 = 0;
			int value2 = 0;

			if (k > 0)
				value1 = binomialCoefficients[i-1][k-1];
			if (k < i)
				value2 = binomialCoefficients[i-1][k];
		
			binomialCoefficients[i][k] = value1 + value2;
		}
	}

	// Now we can proceed with the calculation of the new coefficients
		
	m_offsetCoefficients.resize(M+1);
	
	for (int k = 0 ; k <= M ; k++)
	{
		double fk = 0;
		double blk = 1;
		
		for (int l = k ; l <= M ; l++)
		{
			fk += blk*m_coefficients[l]*binomialCoefficients[l][k];
			blk *= -b;
		}
		m_offsetCoefficients[k] = fk;
	}
}
	
PolynomialMassProfileLensParams::PolynomialMassProfileLensParams()
{
}

PolynomialMassProfileLensParams::~PolynomialMassProfileLensParams()
{
}

void PolynomialMassProfileLensParams::addPolynomialPart(const PolynomialPart &p)
{
	m_polynomialParts.push_back(p);
}

bool PolynomialMassProfileLensParams::write(serut::SerializationInterface &si) const
{
	if (!si.writeInt32(m_polynomialParts.size()))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	for (auto it = m_polynomialParts.begin() ; it != m_polynomialParts.end() ; it++)
	{
		double info[5];
		const PolynomialPart &polyPart = (*it);

		info[0] = polyPart.getXOffset();
		info[1] = polyPart.getYOffset();
		info[2] = polyPart.getXScale();
		info[3] = polyPart.getYScale();
		info[4] = polyPart.getEndPosition();

		if (!si.writeDoubles(info, 5))
		{
			setErrorString(si.getErrorString());
			return false;
		}

		if (!si.writeInt32(polyPart.getCoefficients().size()))
		{
			setErrorString(si.getErrorString());
			return false;
		}

		if (!si.writeDoubles(polyPart.getCoefficients()))
		{
			setErrorString(si.getErrorString());
			return false;
		}
	}

	return true;
}

bool PolynomialMassProfileLensParams::read(serut::SerializationInterface &si)
{
	m_polynomialParts.clear();

	int32_t numParts;

	if (!si.readInt32(&numParts))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	for (int32_t i = 0 ; i < numParts ; i++)
	{
		double info[5];
		int32_t numCoefficients;
		std::vector<double> coefficients;

		if (!si.readDoubles(info, 5))
		{
			setErrorString(si.getErrorString());
			return false;
		}

		if (!si.readInt32(&numCoefficients))
		{
			setErrorString(si.getErrorString());
			return false;
		}

		coefficients.resize(numCoefficients);

		if (!si.readDoubles(coefficients))
		{
			setErrorString(si.getErrorString());
			return false;
		}

		m_polynomialParts.push_back(PolynomialPart(info[0], info[1], info[2], info[3], info[4], coefficients));
	}

	return true;
}

PolynomialZeroMassLensParams::PolynomialZeroMassLensParams(double D_d, double densityScale, double angularRadius, double zeroPoint)
{
	std::vector<double> a(5);
	std::vector<double> b(4);
	double m = zeroPoint;
	double norm = m*m/4.0;
	
	a[0] = 0;
	a[1] = 0;
	a[2] = 2.0/(m*m)*norm;
	a[3] = 0;
	a[4] = -1.0/(m*m*m*m)*norm;
	
	double d = (m-1.0)*(m-1.0)*(m-1.0);
	
	b[0] = (3.0*m-1)/d*norm;
	b[1] = -6.0*m/d*norm;
	b[2] = 3.0*(1.0+m)/d*norm;
	b[3] = -2.0/d*norm;
	
	PolynomialPart part1(0, 0, angularRadius, densityScale*D_d*D_d*angularRadius*angularRadius*2.0*CONST_PI, zeroPoint*angularRadius, a);
	PolynomialPart part2(0, 0, angularRadius, densityScale*D_d*D_d*angularRadius*angularRadius*2.0*CONST_PI, angularRadius, b);

	addPolynomialPart(part1);
	addPolynomialPart(part2);
}
	
PolynomialZeroMassLensParams::~PolynomialZeroMassLensParams()
{
}

PolynomialTimeDelayLensParams::PolynomialTimeDelayLensParams(double theta1, double theta2, double timeDiff, double zLens)
{
	std::vector<double> a(1);
	std::vector<double> b(6);

	a[0] = 0;

	double Yscale = (SPEED_C*SPEED_C*SPEED_C/((1.0+zLens)*4.0*CONST_G))*timeDiff;
	double H = theta1/(theta2-theta1);

	b[0] = 0;
	b[1] = 0;
	b[2] = 30.0*H;
	b[3] = 30.0-60.0*H;
	b[4] = 30.0*H-60.0;
	b[5] = 30.0;

	PolynomialPart part1(0, 0, 1, 1, theta1, a);
	PolynomialPart part2(theta1, 0, theta2-theta1, Yscale, theta2, b);

	addPolynomialPart(part1);
	addPolynomialPart(part2);
}

PolynomialTimeDelayLensParams::~PolynomialTimeDelayLensParams()
{
}

std::unique_ptr<GravitationalLensParams> PolynomialMassProfileLensParams::createCopy() const
{
	std::unique_ptr<PolynomialMassProfileLensParams> pNewParams = std::make_unique<PolynomialMassProfileLensParams>();
	pNewParams->m_polynomialParts = m_polynomialParts;
	return pNewParams;
}

PolynomialMassProfileLens::PolynomialMassProfileLens() : SymmetricLens(PolynomialMassProfile)
{
}

PolynomialMassProfileLens::~PolynomialMassProfileLens()
{
}

bool PolynomialMassProfileLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const PolynomialMassProfileLensParams *pParams = dynamic_cast<const PolynomialMassProfileLensParams *>(pLensParams);

	if (!pParams)
	{
		setErrorString("Parameters are not of type 'PolynomialMassProfileLensParams'");
		return false;
	}

	double prevEndRadius = 0;
	const auto &polynomialParts = pParams->getPolynomialParts();

	for (auto it = polynomialParts.begin() ; it != polynomialParts.end() ; it++)
	{
		double curEndRadius = (*it).getEndPosition();

		if (curEndRadius <= prevEndRadius)
		{
			setErrorString("The specified end positions don't increase strictly monotonically");
			return false;
		}

		prevEndRadius = curEndRadius;
	}

	// Okay, now store the polynomial data
	
	m_polynomialParts.resize(polynomialParts.size());

	auto it = polynomialParts.begin();
	int pos = 0;
	for (  ; it != polynomialParts.end() ; it++, pos++)
		m_polynomialParts[pos] = (*it);

	// just need to calculate the total mass, prevEndRadius contains the furthest end position;
	
	double thetaLength = prevEndRadius;

	m_endPosition = 0;

	if (m_polynomialParts.empty())
		m_totalMass = 0;
	else
	{
		int polyIndex = m_polynomialParts.size()-1;
		double z = (thetaLength - m_polynomialParts[polyIndex].getXOffset())/m_polynomialParts[polyIndex].getXScale();
		const std::vector<double> &coefficients = m_polynomialParts[polyIndex].getCoefficients();
		double sum = 0;
		double zPower = 1;

		for (int i = 0 ; i < (int)coefficients.size() ; i++)
		{
			sum += zPower*coefficients[i];
			zPower *= z;
		}

		m_totalMass = sum*m_polynomialParts[polyIndex].getYScale() + m_polynomialParts[polyIndex].getYOffset();

		m_endPosition = m_polynomialParts[polyIndex].getEndPosition();
	}

	if (m_polynomialParts.empty())
		m_totalPotential = 0;
	else
	{
		double total = 0;
		
		m_potentialOffsets.resize(m_polynomialParts.size()+1);
	
		m_potentialOffsets[0] = 0;
		for (int i = 0 ; i < (int)m_polynomialParts.size() ; i++)
		{
			double partIntegral = 0;
				
			const std::vector<double> &offsetCoefficients = m_polynomialParts[i].getOffsetCoefficients();
		 
			if (offsetCoefficients.size() != 0)
			{
				double z1 = 0;
				double z2 = (m_polynomialParts[i].getEndPosition() - m_polynomialParts[i].getXOffset())/m_polynomialParts[i].getXScale();
				double b = m_polynomialParts[i].getXOffset()/m_polynomialParts[i].getXScale();
				
				if (i > 0)
				{
					double theta1 = m_polynomialParts[i-1].getEndPosition();

					z1 = (theta1 - m_polynomialParts[i].getXOffset())/m_polynomialParts[i].getXScale();
				}

				double part1 = z1+b;
				double part2 = z2+b;
				double logTerm = 0;
				double logCoeff = m_polynomialParts[i].getYScale()*offsetCoefficients[0] + m_polynomialParts[i].getYOffset();
				
				if (logCoeff != 0)
					logTerm = logCoeff*LN(part2/part1);
					
				double part1k = part1;
				double part2k = part2;
				double sum = 0;

				for (int k = 1 ; k < (int)offsetCoefficients.size() ; k++)
				{
					sum += (offsetCoefficients[k]/((double)k))*(part2k-part1k);
					
					part1k *= part1;
					part2k *= part2;
				}

				sum *= m_polynomialParts[i].getYScale();
				sum += logTerm;

				partIntegral = sum;
			}

			total += partIntegral;
			m_potentialOffsets[i+1] = total;
			m_totalPotential = total;
		}
	}
	
	return true;
}

double PolynomialMassProfileLens::getMassInside(double thetaLength) const
{
	int polyIndex = -1;
	int num = (int)m_polynomialParts.size();
	
	if (num == 0)
		return 0;

	for (int i = 0 ; polyIndex < 0 && i < num ; i++)
	{
		if (thetaLength < m_polynomialParts[i].getEndPosition())
			polyIndex = i;
	}

	if (polyIndex < 0)
		return m_totalMass;

	double z = (thetaLength - m_polynomialParts[polyIndex].getXOffset())/m_polynomialParts[polyIndex].getXScale();
	const std::vector<double> &coefficients = m_polynomialParts[polyIndex].getCoefficients();
	double sum = 0;
	double zPower = 1;

	for (int i = 0 ; i < (int)coefficients.size() ; i++)
	{
		sum += zPower*coefficients[i];
		zPower *= z;
	}

	double mass = sum*m_polynomialParts[polyIndex].getYScale() + m_polynomialParts[polyIndex].getYOffset();

	return mass;
}

double PolynomialMassProfileLens::getProfileSurfaceMassDensity(double thetaLength) const
{
	int polyIndex = -1;
	int num = m_polynomialParts.size();
	
	if (num == 0)
		return 0;

	for (int i = 0 ; polyIndex < 0 && i < num ; i++)
	{
		if (thetaLength < m_polynomialParts[i].getEndPosition())
			polyIndex = i;
	}

	if (polyIndex < 0)
		return 0;

	double z = (thetaLength - m_polynomialParts[polyIndex].getXOffset())/m_polynomialParts[polyIndex].getXScale();

#if 0
	const std::vector<double> &coefficients = m_polynomialParts[polyIndex].getCoefficients();
	double sum = 0;
	double zPower = 1;

	for (int i = 1 ; i < coefficients.size() ; i++)
	{
		sum += zPower*coefficients[i]*(double)i;
		zPower *= z;
	}

	double derivative = sum*m_polynomialParts[polyIndex].getYScale()/m_polynomialParts[polyIndex].getXScale();
	double density = derivative*(1.0/(thetaLength*2.0*CONST_PI*getLensDistance()*getLensDistance()));
#else
	const std::vector<double> &offsetCoefficients = m_polynomialParts[polyIndex].getOffsetCoefficients();
	double sum = 0;
	double zbPower = 1;
	double b = m_polynomialParts[polyIndex].getXOffset()/m_polynomialParts[polyIndex].getXScale();

	for (int k = 2 ; k < (int)offsetCoefficients.size() ; k++)
	{
		sum += ((double)k)*offsetCoefficients[k]*zbPower;
		zbPower *= (z + b);
	}

	if (offsetCoefficients.size() > 1 && offsetCoefficients[1] != 0)
		sum += offsetCoefficients[1]/(z + b);

	double density = m_polynomialParts[polyIndex].getYScale()/(m_polynomialParts[polyIndex].getXScale()*m_polynomialParts[polyIndex].getXScale()) 
	                 * sum * 1.0/(2.0*CONST_PI*getLensDistance()*getLensDistance());
#endif
	return density;
}

bool PolynomialMassProfileLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	double thetaLength = theta.getLength();
	int polyIndex = -1;
	int num = (int)m_polynomialParts.size();
	
	if (num == 0)
		return 0;

	for (int i = 0 ; polyIndex < 0 && i < num ; i++)
	{
		if (thetaLength < m_polynomialParts[i].getEndPosition())
			polyIndex = i;
	}

	double factor = ((4.0*CONST_G/(SPEED_C*SPEED_C))/getLensDistance())*(D_ds/D_s);

	if (polyIndex < 0)
	{
		double logPart = 0;

		if (m_totalMass != 0)
			logPart = m_totalMass*LN(thetaLength/m_endPosition);

		*pPotentialValue = factor*(m_totalPotential + logPart);
		return true;
	}

	const std::vector<double> &offsetCoefficients = m_polynomialParts[polyIndex].getOffsetCoefficients();
	double potentialOffset = m_potentialOffsets[polyIndex];
 
	if (offsetCoefficients.size() == 0)
	{
		*pPotentialValue = factor*potentialOffset;
		return true;
	}
	
	double z1 = 0;
	double z2 = (thetaLength - m_polynomialParts[polyIndex].getXOffset())/m_polynomialParts[polyIndex].getXScale();
	double b = m_polynomialParts[polyIndex].getXOffset()/m_polynomialParts[polyIndex].getXScale();
	
	if (polyIndex > 0)
	{
		double theta1 = m_polynomialParts[polyIndex-1].getEndPosition();

		z1 = (theta1 - m_polynomialParts[polyIndex].getXOffset())/m_polynomialParts[polyIndex].getXScale();
	}

	double part1 = z1+b;
	double part2 = z2+b;
	double logTerm = 0;
	double logCoeff = m_polynomialParts[polyIndex].getYScale()*offsetCoefficients[0] + m_polynomialParts[polyIndex].getYOffset();
	
	if (logCoeff != 0)
		logTerm = logCoeff*LN(part2/part1);
		
	double part1k = part1;
	double part2k = part2;
	double sum = 0;

	for (int k = 1 ; k < (int)offsetCoefficients.size() ; k++)
	{
		sum += (offsetCoefficients[k]/((double)k))*(part2k-part1k);
		
		part1k *= part1;
		part2k *= part2;
	}

	sum *= m_polynomialParts[polyIndex].getYScale();
	sum += logTerm;

	*pPotentialValue = factor*(potentialOffset + sum);
	
	return true;
}

} // end namespace


