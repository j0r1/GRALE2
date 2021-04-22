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
#include "multiplewendlandlens.h"
#include "constants.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

namespace grale
{

// Parameters

std::unique_ptr<GravitationalLensParams> MultipleWendlandLensParams::createCopy() const
{
	return std::make_unique<MultipleWendlandLensParams>(m_phiXInfo, m_phiYInfo);
}

bool MultipleWendlandLensParams::write(serut::SerializationInterface &si) const
{
	if (!si.writeInt32(m_phiXInfo.size()))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	for (auto it = m_phiXInfo.begin() ; it != m_phiXInfo.end() ; it++)
	{
		WendlandLensInfo inf = *it;
		double params[4];

		params[0] = inf.getHeightFactor();
		params[1] = inf.getAngularScale();
		params[2] = inf.getAngularPosition().getX();
		params[3] = inf.getAngularPosition().getY();

		if (!si.writeDoubles(params, 4))
		{
			setErrorString(si.getErrorString());
			return false;
		}
	}

	if (!si.writeInt32(m_phiYInfo.size()))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	for (auto it = m_phiYInfo.begin() ; it != m_phiYInfo.end() ; it++)
	{
		WendlandLensInfo inf = *it;
		double params[4];

		params[0] = inf.getHeightFactor();
		params[1] = inf.getAngularScale();
		params[2] = inf.getAngularPosition().getX();
		params[3] = inf.getAngularPosition().getY();

		if (!si.writeDoubles(params, 4))
		{
			setErrorString(si.getErrorString());
			return false;
		}
	}

	return true;
}

bool MultipleWendlandLensParams::read(serut::SerializationInterface &si)
{
	int32_t num;

	m_phiXInfo.clear();
	m_phiYInfo.clear();

	if (!si.readInt32(&num))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	for (int32_t i = 0 ; i < num ; i++)
	{
		double params[4];

		if (!si.readDoubles(params, 4))
		{
			setErrorString(si.getErrorString());
			return false;
		}

		m_phiXInfo.push_back(WendlandLensInfo(params[0], params[1], Vector2D<double>(params[2], params[3])));
	}

	if (!si.readInt32(&num))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	for (int32_t i = 0 ; i < num ; i++)
	{
		double params[4];

		if (!si.readDoubles(params, 4))
		{
			setErrorString(si.getErrorString());
			return false;
		}

		m_phiYInfo.push_back(WendlandLensInfo(params[0], params[1], Vector2D<double>(params[2], params[3])));
	}

	return true;
}

bool MultipleWendlandLensParams::matchDeflections(const std::vector<Vector2D<double> > &deflectionPoints,
			      const std::vector<Vector2D<double> > &deflectionAngles,
			      double angularScale)
{
	if (deflectionAngles.size() != deflectionPoints.size())
	{
		setErrorString("Number of deflection points must match the number of deflection angles");
		return false;
	}
	if (deflectionAngles.size() < 1)
	{
		setErrorString("No deflection information was given");
		return false;
	}

	int numDataPoints = deflectionPoints.size()*2; // data for both x and y

	gsl_matrix *pPhiMatrix = gsl_matrix_alloc(numDataPoints, numDataPoints);
	gsl_vector *pAnglesVector = gsl_vector_alloc(numDataPoints);
	gsl_vector *pCoeffVector = gsl_vector_alloc(numDataPoints);

	for (int i = 0 ; i < numDataPoints ; i += 2)
	{
		Vector2D<double> point = deflectionPoints[i/2];

		for (int j = 0 ; j < numDataPoints ; j += 2)
		{
			Vector2D<double> centerPoint = deflectionPoints[j/2];
			Vector2D<double> scaledDiff = (point-centerPoint)/angularScale;
			double r2 = scaledDiff.getLengthSquared();
			double r = scaledDiff.getLength();

			if (r2 < 1.0)
			{
				gsl_matrix_set(pPhiMatrix, i, j, MultipleWendlandLens::phiXX(r, r2, scaledDiff.getX(), scaledDiff.getY()));
				gsl_matrix_set(pPhiMatrix, i, j+1, MultipleWendlandLens::phiXY(r, r2, scaledDiff.getX(), scaledDiff.getY()));
				gsl_matrix_set(pPhiMatrix, i+1, j, MultipleWendlandLens::phiXY(r, r2, scaledDiff.getX(), scaledDiff.getY()));
				gsl_matrix_set(pPhiMatrix, i+1, j+1, MultipleWendlandLens::phiYY(r, r2, scaledDiff.getX(), scaledDiff.getY()));
			}
			else
			{
				gsl_matrix_set(pPhiMatrix, i, j, 0);
				gsl_matrix_set(pPhiMatrix, i, j+1, 0);
				gsl_matrix_set(pPhiMatrix, i+1, j, 0);
				gsl_matrix_set(pPhiMatrix, i+1, j+1, 0);
			}
		}

		gsl_vector_set(pAnglesVector, i, deflectionAngles[i/2].getX()/angularScale);
		gsl_vector_set(pAnglesVector, i+1, deflectionAngles[i/2].getY()/angularScale);
	}

	gsl_permutation *pPermut = gsl_permutation_alloc(numDataPoints);
	int sigNum;
	int status;
	char str[1024];

	status = gsl_linalg_LU_decomp(pPhiMatrix, pPermut, &sigNum);
	if (status != GSL_SUCCESS)
	{
		sprintf(str, "Can't calculate LU decomposition (error code %d)", (int)status);
		setErrorString(str);
		return false;
	}
	status = gsl_linalg_LU_solve(pPhiMatrix, pPermut, pAnglesVector, pCoeffVector);
	if (status != GSL_SUCCESS)
	{
		sprintf(str, "Can't solve linear system using LU decomposition (error code %d)", (int)status);
		setErrorString(str);
		return false;
	}

	m_phiXInfo.clear();
	m_phiYInfo.clear();

	for (int i = 0 ; i < deflectionPoints.size() ; i++)
	{
		double ca = gsl_vector_get(pCoeffVector, i*2);
		double cb = gsl_vector_get(pCoeffVector, i*2+1);

		m_phiXInfo.push_back(WendlandLensInfo(ca, angularScale, deflectionPoints[i]));
		m_phiYInfo.push_back(WendlandLensInfo(cb, angularScale, deflectionPoints[i]));
	}

	gsl_matrix_free(pPhiMatrix);
	gsl_vector_free(pAnglesVector);
	gsl_vector_free(pCoeffVector);
	gsl_permutation_free(pPermut);
	return true;
}

// Lens

MultipleWendlandLens::MultipleWendlandLens() : GravitationalLens(GravitationalLens::MultipleWendland)
{
}

MultipleWendlandLens::~MultipleWendlandLens()
{
}

bool MultipleWendlandLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const MultipleWendlandLensParams *pParams = dynamic_cast<const MultipleWendlandLensParams *>(pLensParams);
	
	if (!pParams)
	{
		setErrorString("Parameters not of type MultipleWendlandLensParams");
		return false;
	}

	std::vector<WendlandLensInfo>::const_iterator it;
	int i;

	m_phiXInfo.resize(pParams->getPhiXInfo().size());
	m_phiYInfo.resize(pParams->getPhiYInfo().size());

	for (i = 0, it = pParams->getPhiXInfo().begin() ; i < m_phiXInfo.size() ; i++, it++)
		m_phiXInfo[i] = *it;

	for (i = 0, it = pParams->getPhiYInfo().begin() ; i < m_phiYInfo.size() ; i++, it++)
		m_phiYInfo[i] = *it;

	return true;
}

bool MultipleWendlandLens::getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const
{
	Vector2D<double> sum(0, 0);

	for (int i = 0 ; i < m_phiXInfo.size() ; i++)
	{
		Vector2D<double> scaledDiffTheta = (theta-m_phiXInfo[i].getAngularPosition())/m_phiXInfo[i].getAngularScale();
		double r2 = scaledDiffTheta.getLengthSquared();

		if (r2 < 1.0)
		{
			double r = SQRT(r2);

			sum += m_phiXInfo[i].getAngularScale()*m_phiXInfo[i].getHeightFactor()*Vector2D<double>(phiXX(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY()), phiXY(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY()));
		}
	}

	for (int i = 0 ; i < m_phiYInfo.size() ; i++)
	{
		Vector2D<double> scaledDiffTheta = (theta-m_phiYInfo[i].getAngularPosition())/m_phiYInfo[i].getAngularScale();
		double r2 = scaledDiffTheta.getLengthSquared();

		if (r2 < 1.0)
		{
			double r = SQRT(r2);

			sum += m_phiYInfo[i].getAngularScale()*m_phiYInfo[i].getHeightFactor()*Vector2D<double>(phiXY(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY()), phiYY(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY()));
		}
	}

	*pAlpha = sum;
	return true;
}

double MultipleWendlandLens::getSurfaceMassDensity(Vector2D<double> theta) const 
{
	double densFactor = 0.5*SPEED_C*SPEED_C/((4.0*CONST_PI*CONST_G)*getLensDistance());
	double sum = 0;	

	for (int i = 0 ; i < m_phiXInfo.size() ; i++)
	{
		Vector2D<double> scaledDiffTheta = (theta-m_phiXInfo[i].getAngularPosition())/m_phiXInfo[i].getAngularScale();
		double r2 = scaledDiffTheta.getLengthSquared();

		if (r2 < 1.0)
		{
			double r = SQRT(r2);

			sum += m_phiXInfo[i].getHeightFactor()*(phiXXX(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY()) + phiXYY(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY()));
		}
	}

	for (int i = 0 ; i < m_phiYInfo.size() ; i++)
	{
		Vector2D<double> scaledDiffTheta = (theta-m_phiYInfo[i].getAngularPosition())/m_phiYInfo[i].getAngularScale();
		double r2 = scaledDiffTheta.getLengthSquared();

		if (r2 < 1.0)
		{
			double r = SQRT(r2);

			sum += m_phiYInfo[i].getHeightFactor()*(phiXXY(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY()) + phiYYY(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY()));
		}
	}

	return densFactor * sum;
}

bool MultipleWendlandLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	double factor = D_ds/D_s;
	double sum = 0;	

	for (int i = 0 ; i < m_phiXInfo.size() ; i++)
	{
		Vector2D<double> scaledDiffTheta = (theta-m_phiXInfo[i].getAngularPosition())/m_phiXInfo[i].getAngularScale();
		double r2 = scaledDiffTheta.getLengthSquared();

		if (r2 < 1.0)
		{
			double r = SQRT(r2);

			sum += m_phiXInfo[i].getAngularScale()*m_phiXInfo[i].getAngularScale()*m_phiXInfo[i].getHeightFactor()*phiX(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY());
		}
	}

	for (int i = 0 ; i < m_phiYInfo.size() ; i++)
	{
		Vector2D<double> scaledDiffTheta = (theta-m_phiYInfo[i].getAngularPosition())/m_phiYInfo[i].getAngularScale();
		double r2 = scaledDiffTheta.getLengthSquared();

		if (r2 < 1.0)
		{
			double r = SQRT(r2);

			sum += m_phiYInfo[i].getAngularScale()*m_phiYInfo[i].getAngularScale()*m_phiYInfo[i].getHeightFactor()*phiY(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY());
		}
	}

	*pPotentialValue = factor * sum;
	return true;
}

bool MultipleWendlandLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	double axxSum = 0;
	double ayySum = 0;
	double axySum = 0;

	for (int i = 0 ; i < m_phiXInfo.size() ; i++)
	{
		Vector2D<double> scaledDiffTheta = (theta-m_phiXInfo[i].getAngularPosition())/m_phiXInfo[i].getAngularScale();
		double r2 = scaledDiffTheta.getLengthSquared();

		if (r2 < 1.0)
		{
			double r = SQRT(r2);

			axxSum += m_phiXInfo[i].getHeightFactor()*phiXXX(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY());
			ayySum += m_phiXInfo[i].getHeightFactor()*phiXYY(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY());
			axySum += m_phiXInfo[i].getHeightFactor()*phiXXY(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY());
		}
	}

	for (int i = 0 ; i < m_phiYInfo.size() ; i++)
	{
		Vector2D<double> scaledDiffTheta = (theta-m_phiYInfo[i].getAngularPosition())/m_phiYInfo[i].getAngularScale();
		double r2 = scaledDiffTheta.getLengthSquared();

		if (r2 < 1.0)
		{
			double r = SQRT(r2);

			axxSum += m_phiYInfo[i].getHeightFactor()*phiXXY(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY());
			ayySum += m_phiYInfo[i].getHeightFactor()*phiYYY(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY());
			axySum += m_phiYInfo[i].getHeightFactor()*phiXYY(r, r2, scaledDiffTheta.getX(), scaledDiffTheta.getY());
		}
	}

	axx = axxSum;
	ayy = ayySum;
	axy = axySum;

	return true;
}

double MultipleWendlandLens::phiX(double r, double r2, double x, double y)
{
	double A = r-1.0;
	double A2 = A*A;
	double A4 = A2*A2;
	double A7 = A4*A2*A;
	return 22.0*A7*x*(16.0*r2+7.0*r+1.0);
}

double MultipleWendlandLens::phiY(double r, double r2, double x, double y)
{
	double A = r-1.0;
	double A2 = A*A;
	double A4 = A2*A2;
	double A7 = A4*A2*A;
	return 22.0*A7*y*(16.0*r2+7.0*r+1.0);
}

double MultipleWendlandLens::phiXX(double r, double r2, double x, double y)
{
	double A = r-1.0;
	double A2 = A*A;
	double A4 = A2*A2;
	double A6 = A4*A2;
	double x2 = x*x;
	return 22.0*A6*(16.0*r*r2-9.0*r2+144.0*x2*r-6.0*r+24.0*x2-1.0);
}

double MultipleWendlandLens::phiXY(double r, double r2, double x, double y)
{
	double A = r-1.0;
	double A2 = A*A;
	double A4 = A2*A2;
	double A6 = A4*A2;
	return 528.0*A6*(6.0*r+1.0)*x*y;
}

double MultipleWendlandLens::phiYY(double r, double r2, double x, double y)
{
	double A = r-1.0;
	double A2 = A*A;
	double A4 = A2*A2;
	double A6 = A4*A2;
	double y2 = y*y;
	return 22.0*A6*(16.0*r*r2-9.0*r2+144.0*y2*r-6.0*r+24.0*y2-1.0);
}
 
double MultipleWendlandLens::phiXXX(double r, double r2, double x, double y)
{
	double A = r-1.0;
	double A2 = A*A;
	double A5 = A2*A2*A;
	double x2 = x*x;
	double y2 = y*y;
	return 1584.0*A5*x*(20.0*x2-5.0*r-1.0+6.0*y2);
}

double MultipleWendlandLens::phiXXY(double r, double r2, double x, double y)
{
	double A = r-1.0;
	double A2 = A*A;
	double A5 = A2*A2*A;
	double x2 = x*x;
	double y2 = y*y;
	return 528.0*A5*y*(48.0*x2-5.0*r-1.0+6.0*y2);
}

double MultipleWendlandLens::phiXYY(double r, double r2, double x, double y)
{
	double A = r-1.0;
	double A2 = A*A;
	double A5 = A2*A2*A;
	double x2 = x*x;
	double y2 = y*y;
	return 528.0*A5*x*(6.0*x2-5.0*r-1.0+48.0*y2);
}

double MultipleWendlandLens::phiYYY(double r, double r2, double x, double y)
{
	double A = r-1.0;
	double A2 = A*A;
	double A5 = A2*A2*A;
	double x2 = x*x;
	double y2 = y*y;
	return 1584.0*A5*y*(20.0*y2-5.0*r-1.0+6.0*x2);
}

} // end namespace
