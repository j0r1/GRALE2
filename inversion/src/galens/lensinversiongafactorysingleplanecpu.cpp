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

#include "lensinversiongafactorysingleplanecpu.h"
#include "lensinversionparameterssingleplanecpu.h"
#include "compositelens.h"
#include "masssheetlens.h"
#include "lensfitnessobject.h"
#include "imagesdata.h"
#include "imagesdataextended.h"
#include "configurationparameters.h"
#include "positionrandomizationbpwrapper.h"
#include "xoshiro128plus.h"
#include <assert.h>
#include <iostream>
#include <sstream>
#include <vector>

#ifdef SHOWEVOLUTION
#include "gravitationallens.h"
#include <sys/stat.h>
#include <sys/types.h>
#endif // SHOWEVOLUTION

using namespace std;
using namespace errut;

namespace grale
{

class GADeflectionMatrix : public DeflectionMatrix
{
public:
	GADeflectionMatrix(LensInversionGAFactorySinglePlaneCPU *pFactory) : m_pFactory(pFactory)
	{
		assert(pFactory != nullptr);
	}

	~GADeflectionMatrix()
	{
	}
private:
	void log(const string &s) override
	{
		m_pFactory->log(s);
	}

	LensInversionGAFactorySinglePlaneCPU *m_pFactory;
};

LensInversionGAFactorySinglePlaneCPU::LensInversionGAFactorySinglePlaneCPU(unique_ptr<LensFitnessObject> fitObj)
	: LensInversionGAFactoryCommon(move(fitObj))
{
}

void LensInversionGAFactorySinglePlaneCPU::clear()
{
	m_pCurrentParams.reset();
	m_pDeflectionMatrix.reset();
	m_pShortBPMatrix.reset();
	m_pTotalBPMatrix.reset();

	m_basisLenses.clear();
	m_sheetLens.reset();
}

LensInversionGAFactorySinglePlaneCPU::~LensInversionGAFactorySinglePlaneCPU()
{
	clear();
}

bool_t LensInversionGAFactorySinglePlaneCPU::init(const LensInversionParametersBase &p)
{
	if (m_pCurrentParams.get() != 0)
		return "Already initialized";

	const LensInversionParametersSinglePlaneCPU *p2 = dynamic_cast<const LensInversionParametersSinglePlaneCPU *>(&p);
	if (!p2)
		return "Invalid type of GA factory parameters";

	m_pCurrentParams = p2->createCopy();
	assert(m_pCurrentParams.get());

	auto basisLenses = m_pCurrentParams->getBasisLenses();
	if (basisLenses.size() == 0)
	{
		clear();
		return "Unable to get basis lenses: " + m_pCurrentParams->getErrorString();
	}
	
	for (const auto &bl : basisLenses)
	{
		m_basisLenses.push_back( { bl.m_pLens, bl.m_center });
		
		double m = bl.m_relevantLensingMass;
		if (m < 0) // TODO: allow this to be overridden
		{
			clear();
			return "A basis lens has a negative strong lensing mass";
		}
	}

	const GravitationalLens *pSheetLens = p2->getSheetLens();
	int numSheetValues = 0;
	if (pSheetLens)
	{
		m_sheetLens = shared_ptr<GravitationalLens>(pSheetLens->createCopy());
		numSheetValues = 1;
	}
	else // no sheet
		m_sheetLens.reset();

	bool_t r;
	// perform sub-initialization
	if (!(r = localSubInit(p2->getZ_d(), p2->getImages(), m_basisLenses, p2->getBaseLens(), pSheetLens, 
				m_pCurrentParams->getFitnessObjectParameters())))
	{
		clear();
		return r;
	}
	
	bool allowNegativeValues = m_pCurrentParams->allowNegativeValues();
	auto massScaleSearchParams = m_pCurrentParams->getMassScaleSearchParameters();
	double massScale = m_pCurrentParams->getMassScale();

	vector<double> basisFunctionMasses;
	for (const auto &bl : basisLenses)
		basisFunctionMasses.push_back(bl.m_relevantLensingMass);

	if (!(r = setCommonParameters(numSheetValues, allowNegativeValues,
	                         basisFunctionMasses,
							 massScale, massScale,
							 massScaleSearchParams)))
	{
		clear();
		return r;
	}

	m_prevIteration = numeric_limits<size_t>::max();

	return true;
}

unique_ptr<GravitationalLens> LensInversionGAFactorySinglePlaneCPU::createLens(const std::vector<float> &basisFunctionWeights,
	                                      const std::vector<float> &sheetValues,
										  float scaleFactor,
										  std::string &errStr) const
{
	CompositeLensParams lensParams;

	float sheetValue = 0;
	if (sheetValues.size() > 0)
	{
		sheetValue = sheetValues[0];
		if (sheetValues.size() > 1)
		{
			errStr = "Unexpected error: more than one (" + to_string(sheetValues.size()) +") sheet value";
			return nullptr;
		}
	}

	for (int i = 0 ; i < m_basisLenses.size() ; i++)
		lensParams.addLens((double)basisFunctionWeights[i]*(double)scaleFactor, m_basisLenses[i].second, 0, *m_basisLenses[i].first.get());

	if (m_sheetLens.get())
		lensParams.addLens(sheetValue, Vector2D<double>(0,0), 0, *m_sheetLens.get());

	auto pLens = make_unique<CompositeLens>();
	if (!pLens->init(m_basisLenses[0].first->getLensDistance(), &lensParams))
	{
		errStr = "Can't init composite lens: " + pLens->getErrorString();
		return nullptr;
	}

	// Note that this lens does not include any base lens that might be set
	return pLens;
}

// This is some very old code, haven't used this for a long time. I'm just leaving this
// in here in case I'd like to do something similar in the future.
#ifdef SHOWEVOLUTION

void GridLensInversionGAFactoryBase::onSortedPopulation(const std::vector<mogal::GenomeWrapper> &population)
{
#if 0
	char str[1024];
	int generation = getCurrentAlgorithm()->getCurrentGeneration();

	sprintf(str,"g_%d",generation);
	if(mkdir(str,S_IRWXU) != 0)
		return;
	
	for (int i = 0 ; i < (int)population.size() ; i++)
	{
		char str[1024];
		sprintf(str,"g_%d/p_%d_%d_%d.lensdata",generation,i,population[i].getParent1(),population[i].getParent2());

		const LensInversionGenomeBase *g = (const LensInversionGenomeBase *)population[i].getGenome();

		double m = 0;
		std::string errstr;
		GravitationalLens *lens = g->createLens(&m,errstr);
		lens->save(str);
		delete lens;
	}
#else
	int generation = getCurrentAlgorithm()->getCurrentGeneration();
	char str[1024];

	if (generation == 0)
		mkdir("bestgenomes",S_IRWXU);

	sprintf(str,"bestgenomes/gen_%d.lensdata",generation);

	const LensInversionGenomeBase *g = (const LensInversionGenomeBase *)population[0].getGenome();

	double m = 0;
	std::string errstr;
	GravitationalLens *lens = g->createLens(&m,errstr);
	lens->save(str);
	delete lens;
#endif
}

#endif // SHOWEVOLUTION

void getUncertaintyInfo(const vector<ImagesDataExtended*> reducedImagesVector,
					    vector<bool> &srcHaveUncert, vector<double> &allUncerts)
// For randomization tests
{
	for (size_t s = 0 ; s < reducedImagesVector.size() ; s++)
	{
		const ImagesDataExtended &img = *reducedImagesVector[s];

		if (!img.hasProperty(ImagesData::PositionUncertainty))
		{
			srcHaveUncert.push_back(false);
			continue;
		}

		int numImg = img.getNumberOfImages();
		for (int i = 0 ; i < numImg ; i++)
		{
			int numPoints = img.getNumberOfImagePoints(i);
			for (int p = 0 ; p < numPoints ; p++)
				allUncerts.push_back(img.getImagePointProperty(ImagesData::PositionUncertainty, i, p));
		}

		srcHaveUncert.push_back(true);
	}
}

bool_t LensInversionGAFactorySinglePlaneCPU::localSubInit(double z_d, const vector<shared_ptr<ImagesDataExtended>> &images, 
	                  const vector<pair<shared_ptr<GravitationalLens>, Vector2D<double> > > &basisLenses,
                      const GravitationalLens *pBaseLens, const GravitationalLens *pSheetLens, 
					  const ConfigurationParameters *pFitnessObjectParams)
{
	vector<ImagesDataExtended*> reducedImagesVector;
	vector<ImagesDataExtended*> shortImagesVector;
	bool_t r;

	if (!(r = initializeLensFitnessObject(z_d, images, pFitnessObjectParams, reducedImagesVector, shortImagesVector)))
		return r;

	m_pDeflectionMatrix = make_unique<GADeflectionMatrix>(this);
	if (!m_pDeflectionMatrix->startInit())
		return "Can't start deflection matrix initialization: " + m_pDeflectionMatrix->getErrorString();

	std::vector<bool> totalDeflectionFlags, totalDerivativeFlags, totalPotentialFlags, totalSecondDerivFlags;

	LensFitnessObject &fitnessObject = getFitnessObject();
	fitnessObject.getTotalCalcFlags(totalDeflectionFlags, totalDerivativeFlags, totalPotentialFlags, totalSecondDerivFlags);

	vector<bool> srcHaveUncert;
	vector<double> allUncerts;

	if (m_pCurrentParams->getRandomizeInputPositions())
	{
		getUncertaintyInfo(reducedImagesVector, srcHaveUncert, allUncerts);
		bool anyHaveUncert = false;
		for (auto x : srcHaveUncert)
			anyHaveUncert |= x;
		
		if (!anyHaveUncert)
			return "Randomization of input positions requested, but no position uncertainties were specified";

		// Make sure we have the information for the necessary calculations
		assert(totalDerivativeFlags.size() == reducedImagesVector.size());
		assert(totalDeflectionFlags.size() == reducedImagesVector.size());
		for (size_t i = 0 ; i < srcHaveUncert.size() ; i++)
		{
			if (srcHaveUncert[i])
			{
				assert(i < totalDerivativeFlags.size() && i < totalDeflectionFlags.size());
				totalDeflectionFlags[i] = true; // We may need this for adjustments to the lens potential
				totalDerivativeFlags[i] = true; // We need this to calculate small shifts in position, to map them to source plane
			}
		}

		uint64_t posUncertSeed = m_pCurrentParams->getInitialPositionUncertaintySeed();
		if (posUncertSeed == 0)
			return "If random offsets of input positions is requested, the initial seed may not be 0";

		cerr << "INFO: enabling EXPERIMENTAL positional uncertainty with seed " << posUncertSeed << endl;

		// Init the rng for the offsets
		for (size_t i = 0 ; i < 4 ; i++)
			m_rndRngState[i] = (uint32_t)(posUncertSeed >> (i*16));
		xoshiro128plus::jump(m_rndRngState);
	}
	else
	{
		if (m_pCurrentParams->getInitialPositionUncertaintySeed() != 0)
			return "If random offsets of input positions is disabled, the initial seed must be 0";

		cerr << "INFO: NOT enabling EXPERIMENTAL positional uncertainties" << endl;
	}
	
	m_pTotalBPMatrix = make_shared<BackProjectMatrix>();
	if (!(r = m_pTotalBPMatrix->startInit(z_d, basisLenses[0].first->getLensDistance(), m_pDeflectionMatrix.get(),
				         reducedImagesVector, totalDeflectionFlags, totalDerivativeFlags, 
					     totalPotentialFlags, totalSecondDerivFlags, pBaseLens, pSheetLens)))
		return "Can't start back projection matrix initialization: " + m_pTotalBPMatrix->getErrorString();

	if (shortImagesVector.size() > 0) // Ok, can improve speed in mass-scale calculation
	{
		std::vector<bool> shortDeflectionFlags, shortDerivativeFlags, shortPotentialFlags, shortSecondDerivFlags;

		fitnessObject.getShortCalcFlags(shortDeflectionFlags, shortDerivativeFlags, shortPotentialFlags, shortSecondDerivFlags);

		m_pShortBPMatrix = make_shared<BackProjectMatrix>();
		if (!(r = m_pShortBPMatrix->startInit(z_d, basisLenses[0].first->getLensDistance(), m_pDeflectionMatrix.get(),
						 shortImagesVector, shortDeflectionFlags, shortDerivativeFlags, 
						 shortPotentialFlags, shortSecondDerivFlags, pBaseLens, pSheetLens)))
			return "Can't start short back projection matrix initialization: " + r.getErrorString();
	}
	else
		m_pShortBPMatrix = m_pTotalBPMatrix;

	bool err = false;

	if (!(r = m_pDeflectionMatrix->endInit(basisLenses)))
		return "Can't end deflection matrix initialization: " + m_pDeflectionMatrix->getErrorString();

	if (!m_pTotalBPMatrix->endInit())
		return "Can't end total backprojection matrix initialization: " + m_pTotalBPMatrix->getErrorString();

	if (m_pShortBPMatrix.get() != m_pTotalBPMatrix.get())
	{
		if (!m_pShortBPMatrix->endInit())
			return "Can't end short backprojection matrix initialization: " + m_pShortBPMatrix->getErrorString();
	}

	// TODO: should we ignore the small randomizations for the shortbpmatrix?
	//       since that's just use to get a scale, this might make sense

	if (m_pCurrentParams->getRandomizeInputPositions())
	{
		m_rndBpWrapper = make_shared<PositionRandomizationBackprojectWrapper>(m_pTotalBPMatrix);
		size_t numRndPoints = 0;
		if (!(r = m_rndBpWrapper->initRandomization(srcHaveUncert, numRndPoints)))
			return "Can't init randomization wrapper: " + r.getErrorString();

		m_rndSigmas.clear();
		for (auto x : allUncerts) // convert to floats
			m_rndSigmas.push_back((float)(x/m_pTotalBPMatrix->getAngularScale()));
	}

	return true;
}

bool_t LensInversionGAFactorySinglePlaneCPU::onNewCalculationStart(size_t iteration, size_t genomesForThisCalculator, size_t genomesForPopulationCalculator)
{
	if (iteration == m_prevIteration) // Don't think this actually happens
		return true;

	m_prevIteration = iteration;
	if (m_rndBpWrapper.get()) // We have a wrapper for randomized input positions, adjust this
	{
		m_rndOffsets.clear(); // calculate these based on the sigmas
		for (auto s : m_rndSigmas)
		{
			Vector2Df dTheta = xoshiro128plus::next_gaussians(m_rndRngState);
			dTheta *= s;
			m_rndOffsets.push_back(dTheta);
		}

		bool_t r;
		if (!(r = m_rndBpWrapper->setRandomOffsets(m_rndOffsets)))
			return "Can't set input position offsets: " + r.getErrorString();

		//cerr << "Set " << m_rndOffsets.size() << " new offsets" << endl;
	}
	return true;
}

bool_t LensInversionGAFactorySinglePlaneCPU::initializeNewCalculation(const vector<float> &basisFunctionWeights, const vector<float> &sheetValues)
{
	if (sheetValues.size() == 0)
		m_sheetScale = 0;
	else if (sheetValues.size() == 1)
		m_sheetScale = sheetValues[0];
	else
		return "Unexpected: got more (" + to_string(sheetValues.size()) + ") mass sheet contributions than expected";

	m_pDeflectionMatrix->calculateBasisMatrixProducts(basisFunctionWeights, true, true, true, true);
	m_pTotalBPMatrix->storeDeflectionMatrixResults();
	m_pShortBPMatrix->storeDeflectionMatrixResults();
	return true;
}

bool_t LensInversionGAFactorySinglePlaneCPU::calculateMassScaleFitness(float scaleFactor, float &fitness)
{
	LensFitnessObject &fitnessFunction = getFitnessObject();

	m_pShortBPMatrix->calculate(scaleFactor, m_sheetScale);
	if (fitnessFunction.shortNeedInverseMagnifications())
		m_pShortBPMatrix->calculateInverseMagnifications(*(fitnessFunction.getShortInverseMagnificationFlags()));
	if (fitnessFunction.shortNeedShearComponents())
		m_pShortBPMatrix->calculateShearComponents(*(fitnessFunction.getShortShearComponentFlags()));
	if (fitnessFunction.shortNeedConvergence())
		m_pShortBPMatrix->calculateConvergence(*(fitnessFunction.getShortConvergenceFlags()));

	if (!fitnessFunction.calculateMassScaleFitness(*m_pShortBPMatrix, fitness))
		return "Unable to calculate mass scale fitness: " + fitnessFunction.getErrorString();

	return true;
}

bool_t LensInversionGAFactorySinglePlaneCPU::calculateTotalFitness(float scaleFactor, float *pFitnessValues)
{
	LensFitnessObject &fitnessFunction = getFitnessObject();

	m_pTotalBPMatrix->calculate(scaleFactor, m_sheetScale);
	if (fitnessFunction.totalNeedInverseMagnifications())
		m_pTotalBPMatrix->calculateInverseMagnifications(*(fitnessFunction.getTotalInverseMagnificationFlags()));
	if (fitnessFunction.totalNeedShearComponents())
		m_pTotalBPMatrix->calculateShearComponents(*(fitnessFunction.getTotalShearComponentFlags()));
	if (fitnessFunction.totalNeedConvergence())
		m_pTotalBPMatrix->calculateConvergence(*(fitnessFunction.getTotalConvergenceFlags()));

	/*
	auto dumpBp = [](const ProjectedImagesInterface &iface)
	{
		double factor = iface.getAngularScale()/ANGLE_ARCSEC;
		for (int s = 0 ; s < iface.getNumberOfSources() ; s++)
		{
			int numPoints = iface.getNumberOfImagePoints(s);
			for (int p = 0 ; p < numPoints ; p++)
			{
				Vector2Df pt = iface.getBetas(s)[p];
				cout << "  " << pt.getX()*factor << " " << pt.getY()*factor << endl;
			}
		}
		cout << endl;
	};
	dumpBp(*m_pTotalBPMatrix);
	*/

	ProjectedImagesInterface *pIface = nullptr;
	if (m_rndBpWrapper.get())
	{
		pIface = m_rndBpWrapper.get();
		m_rndBpWrapper->clearCachedValues();
		// cerr << "Cleared cached values for new calculation" << endl;
	}
	else
		pIface = m_pTotalBPMatrix.get();

	if (!fitnessFunction.calculateOverallFitness(*pIface, pFitnessValues))
		return "Unable to calculate full fitness: " + fitnessFunction.getErrorString();

	return true;
}

} // end namespace
