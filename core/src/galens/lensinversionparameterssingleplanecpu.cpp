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
#include "lensinversionparameterssingleplanecpu.h"
#include "imagesdataextended.h"
#include "gravitationallens.h"
#include "configurationparameters.h"
#include "plummerlens.h"
#include "squarelens.h"
#include "gausslens.h"
#include "masssheetlens.h"
#include <errut/booltype.h>
#include <assert.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <memory>

using namespace std;
using namespace errut;

namespace grale
{

void LensInversionParametersSinglePlaneCPU::printBasisLenses()
{
	if (!getenv("GRALE_DEBUG_PRINTSINGLEPLANEBASISLENSES"))
		return;

	auto oldPrec = cerr.precision();
	cerr.precision(10);
	cerr << "DEBUG: " << m_basisLenses.size() << " basis lenses, mass scale is " << m_massScale/MASS_SOLAR << " M_sun" << endl;

	for (size_t i = 0 ; i < m_basisLenses.size() ; i++)
	{
		cerr << "  Lens " << i << " is ";
		auto &bl = m_basisLenses[i];

		const PlummerLensParams *pPlummerParams = dynamic_cast<const PlummerLensParams *>(bl.m_pLens->getLensParameters());
		if (pPlummerParams)
		{
			cerr << "Plummer lens with parameters:" << endl;
			cerr << "     mass:  " << pPlummerParams->getLensMass()/MASS_SOLAR << " M_sun" << endl;
			cerr << "     width: " << pPlummerParams->getAngularWidth()/ANGLE_ARCSEC << " arcsec" << endl;
		}
		else
			cerr << "not a plummer lens" << endl;

		cerr << "   Center:        " << bl.m_center.getX()/ANGLE_ARCSEC << ", " << bl.m_center.getY()/ANGLE_ARCSEC << " arcsec" << endl;
		cerr << "   Relevant mass: " << bl.m_relevantLensingMass/MASS_SOLAR << " M_sun" << endl;
	}
	cerr.precision(oldPrec);
}

LensInversionParametersSinglePlaneCPU::LensInversionParametersSinglePlaneCPU()
{
	zero();
}

void LensInversionParametersSinglePlaneCPU::commonConstructor(
			const vector<shared_ptr<ImagesDataExtended>> &images,
			double D_d,
			double z_d,
			double massScale,
			bool allowNegativeValues,
			const GravitationalLens *pBaseLens,
			const GravitationalLens *pSheetLens,
			const ConfigurationParameters *pFitnessObjectParams,
			const ScaleSearchParameters &massScaleSearchParams,
			bool randomizeImagePositions,
			uint64_t initialUncertSeed
			)
{
	zero();

	m_allowNegative = allowNegativeValues;

	m_Dd = D_d;
	m_zd = z_d;
	m_massScale = massScale; 
	m_images = images;
	
	m_pBaseLens.reset();
	if (pBaseLens)
		m_pBaseLens = shared_ptr<GravitationalLens>(pBaseLens->createCopy());

	m_pSheetLens.reset();
	if (pSheetLens)
		m_pSheetLens = shared_ptr<GravitationalLens>(pSheetLens->createCopy());

	m_pParams.reset();
	if (pFitnessObjectParams)
		m_pParams = make_shared<ConfigurationParameters>(*pFitnessObjectParams);

	m_scaleSearchParams = massScaleSearchParams;
	m_randomizeImagePositions = randomizeImagePositions;
	m_initialUncertSeed = initialUncertSeed;
}

LensInversionParametersSinglePlaneCPU::LensInversionParametersSinglePlaneCPU(
			const vector<shared_ptr<ImagesDataExtended>> &images,
			const vector<LensInversionBasisLensInfo> &basisLenses,
			double D_d,
			double z_d,
			double massScale,
			bool allowNegativeValues,
			const GravitationalLens *pBaseLens,
			const GravitationalLens *pSheetLens,
			const ConfigurationParameters *pFitnessObjectParams,
			const ScaleSearchParameters &massScaleSearchParams,
			bool randomizeImagePositions,
			uint64_t initialUncertSeed
			)
{
	commonConstructor(images, D_d, z_d, massScale, allowNegativeValues, pBaseLens,
			          pSheetLens, pFitnessObjectParams, massScaleSearchParams,
					  randomizeImagePositions, initialUncertSeed);

	m_basisLenses = basisLenses;
	printBasisLenses();
}

LensInversionParametersSinglePlaneCPU::LensInversionParametersSinglePlaneCPU(
		const vector<shared_ptr<ImagesDataExtended>> &images,
		const vector<GridSquare> &gridsquares,
		double D_d, double z_d, double massScale,
		bool useweights,
		BasisFunctionType b, bool allowNegativeValues,
		const GravitationalLens *pBaseLens,
		const GravitationalLens *pSheetLens,
		const ConfigurationParameters *pFitnessObjectParams,
		const ScaleSearchParameters &massScaleSearchParams,
		bool randomizeImagePositions,
		uint64_t initialUncertSeed
		)
{
	commonConstructor(images, D_d, z_d, massScale, allowNegativeValues, pBaseLens,
			          pSheetLens, pFitnessObjectParams, massScaleSearchParams,
					  randomizeImagePositions, initialUncertSeed);

	buildBasisLenses(gridsquares, b, useweights);
	printBasisLenses();
}

LensInversionParametersSinglePlaneCPU::~LensInversionParametersSinglePlaneCPU()
{
	clear();
}

bool LensInversionParametersSinglePlaneCPU::write(serut::SerializationInterface &si) const
{
	if (!si.writeDouble(m_Dd))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeDouble(m_zd))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeDouble(m_massScale))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	int32_t sizes[2];

	sizes[0] = m_images.size();
	sizes[1] = m_basisLenses.size();

	if (!si.writeInt32s(sizes,2))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	for (auto pImg : m_images)
	{
		assert(pImg.get());

		if (!pImg->write(si))
		{
			setErrorString(pImg->getErrorString());
			return false;
		}
	}

	for (const auto &bl : m_basisLenses)
	{
		assert(bl.m_pLens.get());
		if (!bl.m_pLens->write(si))
		{
			setErrorString("Unable to write a basis lens model: " + bl.m_pLens->getErrorString());
			return false;
		}

		const double centerAndMass[3] = { bl.m_center.getX(), bl.m_center.getY(), bl.m_relevantLensingMass };
		if (!si.writeDoubles(centerAndMass, 3))
		{
			setErrorString("Unable to write basis lens model position or mass: " + si.getErrorString());
			return false;
		}
	}

	int32_t n = (m_allowNegative)?1:0;

	if (!si.writeInt32(n))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	auto writeLens = [](const shared_ptr<GravitationalLens> pLens, serut::SerializationInterface &si) -> bool_t
	{
		int32_t val = (pLens.get() != nullptr)?1:0;
		if (!si.writeInt32(val))
			return "Couldn't write lens flag: " + si.getErrorString();

		if (pLens.get())
		{
			if (!pLens->write(si))
				return "Couldn't write lens: " + pLens->getErrorString();
		}
		return true;
	};

	bool_t r;
	if (!(r = writeLens(m_pBaseLens, si)))
	{
		setErrorString("Couldn't write base lens: " + r.getErrorString());
		return false;
	}
	if (!(r = writeLens(m_pSheetLens, si)))
	{
		setErrorString("Couldn't write sheet lens: " + r.getErrorString());
		return false;
	}

	int32_t val = (m_pParams.get() == nullptr)?0:1;
	if (!si.writeInt32(val)) // write marker for additional parameters
	{
		setErrorString(si.getErrorString());
		return false;
	}

	if (m_pParams.get())
	{
		if (!m_pParams->write(si))
		{
			setErrorString(m_pParams->getErrorString());
			return false;
		}
	}

	if (!m_scaleSearchParams.write(si))
	{
		setErrorString(m_scaleSearchParams.getErrorString());
		return false;
	}

	vector<int32_t> vals = { (m_randomizeImagePositions)?1:0,
		                     (int32_t)((uint32_t)((m_initialUncertSeed >> 0) & 0xffffffff)),
							 (int32_t)((uint32_t)((m_initialUncertSeed >> 32) & 0xffffffff)) };

	if (!si.writeInt32s(vals))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	cerr << "JORI: wrote m_initialUncertSeed = " << m_initialUncertSeed << endl;

	return true;
}

void LensInversionParametersSinglePlaneCPU::clear()
{
	m_images.clear();
	m_basisLenses.clear();

	m_pBaseLens.reset();
	m_pSheetLens.reset();
	m_pParams.reset();

	zero();
}

void LensInversionParametersSinglePlaneCPU::zero()
{
	m_Dd = 0;
	m_massScale = 0;
	m_zd = 0;
	m_allowNegative = false;
	m_pBaseLens.reset();
	m_pSheetLens.reset();
	m_pParams.reset();

	m_randomizeImagePositions = false;
	m_initialUncertSeed = 0;
}

bool LensInversionParametersSinglePlaneCPU::read(serut::SerializationInterface &si)
{
	clear();
	
	if (!si.readDouble(&m_Dd))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.readDouble(&m_zd))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.readDouble(&m_massScale))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	int32_t sizes[2];

	if (!si.readInt32s(sizes,2))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	int32_t numimages = sizes[0];
	int32_t numBasisFunctions = sizes[1];
	
	for (int32_t i = 0 ; i < numimages ; i++)
	{
		shared_ptr<ImagesDataExtended> img = make_shared<ImagesDataExtended>();
		if (!img->read(si))
		{
			setErrorString(img->getErrorString());
			return false;
		}
		m_images.push_back(img);
	}

	for (int32_t i = 0 ; i < numBasisFunctions ; i++)
	{
		string errStr;
		unique_ptr<GravitationalLens> pBasisFunction;
		if (!GravitationalLens::read(si, pBasisFunction, errStr))
		{
			setErrorString("Unable to read basis lens model: " + errStr);
			return false;
		}
		shared_ptr<GravitationalLens> basisLens(move(pBasisFunction));

		double centerAndMass[3];
		if (!si.readDoubles(centerAndMass, 3))
		{
			setErrorString("Unable to read center or mass of a basis lens model:" + si.getErrorString());
			return false;
		}

		m_basisLenses.push_back( { basisLens, { centerAndMass[0], centerAndMass[1] }, centerAndMass[2] });
	}
	
	int32_t n;

	if (!si.readInt32(&n))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	m_allowNegative = (n == 0)?false:true;

	auto readLens = [](shared_ptr<GravitationalLens> &pLens, serut::SerializationInterface &si) -> bool_t
	{
		int32_t val;

		if (!si.readInt32(&val))
			return "Couldn't read lens flag: " + si.getErrorString();

		if (val != 0)
		{
			string errStr;
			unique_ptr<GravitationalLens> pTmpLens;

			if (!GravitationalLens::read(si, pTmpLens, errStr))
				return "Couldn't read lens: " + errStr; 
			pLens = shared_ptr<GravitationalLens>(move(pTmpLens));
		}
		return true;
	};

	bool_t r;
	if (!(r = readLens(m_pBaseLens, si)))
	{
		setErrorString("Couldn't read base lens: " + r.getErrorString());
		return false;
	}
	if (!(r = readLens(m_pSheetLens, si)))
	{
		setErrorString("Couldn't read sheet lens: " + r.getErrorString());
		return false;
	}

	int32_t val = 0;
	if (!si.readInt32(&val))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	if (val != 0)
	{
		m_pParams = make_shared<ConfigurationParameters>();
		if (!m_pParams->read(si))
		{
			setErrorString(m_pParams->getErrorString());
			return false;
		}
	}

	if (!m_scaleSearchParams.read(si))
	{
		setErrorString(m_scaleSearchParams.getErrorString());
		return false;
	}

	vector<int32_t> vals(3);
	if (!si.readInt32s(vals))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	m_randomizeImagePositions = (vals[0] == 0)?false:true;
	m_initialUncertSeed = (((uint64_t)((uint32_t)vals[2])) << 32 ) | (((uint64_t)((uint32_t)vals[1])) << 0 );
	cerr << "JORI: read m_initialUncertSeed = " << m_initialUncertSeed << endl;
	return true;
}

unique_ptr<LensInversionParametersSinglePlaneCPU> LensInversionParametersSinglePlaneCPU::createCopy() const
{
	vector<LensInversionBasisLensInfo> copiedBasisLenses;

	for (const auto &bl : m_basisLenses)
	{
		assert(bl.m_pLens);
		shared_ptr<GravitationalLens> newLens(bl.m_pLens->createCopy());
		if (!newLens.get())
		{
			setErrorString("Unable to create a copy from one of the basis lenses: " + bl.m_pLens->getErrorString());
			return nullptr;
		}

		copiedBasisLenses.push_back( { newLens, bl.m_center, bl.m_relevantLensingMass } );
	}

	vector<shared_ptr<ImagesDataExtended>> imagesCopy;

	for (const auto &img : m_images)
	{
		assert(img.get());
		shared_ptr<ImagesDataExtended> imgCopy = make_shared<ImagesDataExtended>(*img);
		if (!imgCopy.get())
		{
			setErrorString("Unable to create copy of an images data set: " + imgCopy->getErrorString());
			return nullptr;
		}
		imagesCopy.push_back(imgCopy);
	}

	return make_unique<LensInversionParametersSinglePlaneCPU>(imagesCopy, copiedBasisLenses,
			                               m_Dd, m_zd, m_massScale, m_allowNegative,
										   m_pBaseLens.get(), m_pSheetLens.get(), m_pParams.get(),
										   m_scaleSearchParams, m_randomizeImagePositions,
										   m_initialUncertSeed);
}

void LensInversionParametersSinglePlaneCPU::buildBasisLenses(const vector<GridSquare> &squares,
                                                   BasisFunctionType basisFunctionType,
												   bool useMassWeights)
{
	const int numMasses = squares.size();
	vector<double> massWeights(numMasses);

	if (numMasses == 0)
		return;

	if (useMassWeights)
	{
		// cerr << "Using mass weights" << endl;
		
		double sum = 0;

		auto it = squares.begin();
		for (it = squares.begin(); it != squares.end() ; it++)
		{
			double w = (*it).getSize();
			double w2 = (w*w);

			sum += w2;
		}
		
		// on a uniform grid, we'd like each weight to be one
		it = squares.begin();
		for (int i = 0 ; it != squares.end() ; i++, it++)
		{
			double w = (*it).getSize();
			double w2 = (w*w);

			massWeights[i] = ((((double)numMasses)*w2)/sum);
			massWeights[i] /= squares.size();
		}
	}
	else
	{
		int i = 0;

		for (auto it = squares.begin(); it != squares.end() ; it++, i++)
			massWeights[i] = 1.0/squares.size();
	}

	// Create the basis functions

	double massScale = getMassScale();
	double D_d = getD_d();
	int squareNumber = 0;
	vector<LensInversionBasisLensInfo> basisLenses;

	for (auto squareIt = squares.begin() ; squareIt != squares.end() ; squareIt++, squareNumber++)
	{
		double gridSize = squareIt->getSize();
		Vector2D<double> gridCenter = squareIt->getCenter();
		bool returnValue;
		shared_ptr<GravitationalLens> pLens;
	
		if (basisFunctionType == LensInversionParametersSinglePlaneCPU::PlummerBasis)
		{
			PlummerLensParams lensParams(massScale * massWeights[squareNumber], gridSize);

			pLens = make_shared<PlummerLens>();
			returnValue = pLens->init(D_d, &lensParams);
		}
		else if (basisFunctionType == LensInversionParametersSinglePlaneCPU::GaussBasis)
		{
			GaussLensParams lensParams(massScale * massWeights[squareNumber], gridSize);

			pLens = make_shared<GaussLens>();
			returnValue = pLens->init(D_d, &lensParams);
		}
		else // Squares
		{
			assert(basisFunctionType == LensInversionParametersSinglePlaneCPU::SquareBasis);

			SquareLensParams lensParams(massScale * massWeights[squareNumber], gridSize);

			pLens = make_shared<SquareLens>();
			returnValue = pLens->init(D_d, &lensParams);
		}

		if (!returnValue)
		{
			cerr << "Couldn't initialize basis lens: " << pLens->getErrorString() << endl;
			return;
		}

		basisLenses.push_back({ pLens, gridCenter, massScale*massWeights[squareNumber] });
	}

	m_basisLenses = basisLenses;
}

shared_ptr<GravitationalLens> LensInversionParametersSinglePlaneCPU::createDefaultSheetLens(MassSheetSearchType t, double Dd)
{
	shared_ptr<GravitationalLens> lens(nullptr);
	if (t == Genome)
	{
		MassSheetLensParams params(Dd, 1.1, 1); // should give the same behaviour als the old mass sheet implementation
		lens = make_shared<MassSheetLens>();
		lens->init(Dd, &params);
	}

	return lens;
}

void LensInversionParametersSinglePlaneCPU::copyFrom(const LensInversionParametersSinglePlaneCPU &src)
{
	m_Dd = src.m_Dd;
	m_massScale = src.m_massScale;
	m_zd = src.m_zd;
	m_images = src.m_images;
	m_allowNegative = src.m_allowNegative;
	m_pBaseLens = src.m_pBaseLens;
	m_pSheetLens = src.m_pSheetLens;
	m_pParams = src.m_pParams;
	m_scaleSearchParams = src.m_scaleSearchParams;
	m_basisLenses = src.m_basisLenses;
	m_initialUncertSeed = src.m_initialUncertSeed;
	m_randomizeImagePositions = src.m_randomizeImagePositions;
}

} // end namespace

