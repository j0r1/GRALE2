/*
#if 0

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
#include "configurationparameters.h"
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
		m_pFactory->sendMessage(s);
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

	std::vector<bool> totalDeflectionFlags, totalDerivativeFlags, totalPotentialFlags;

	LensFitnessObject &fitnessObject = getFitnessObject();
	fitnessObject.getTotalCalcFlags(totalDeflectionFlags, totalDerivativeFlags, totalPotentialFlags);
	
	m_pTotalBPMatrix = make_shared<BackProjectMatrix>();
	if (!(r = m_pTotalBPMatrix->startInit(z_d, basisLenses[0].first->getLensDistance(), m_pDeflectionMatrix.get(),
				         reducedImagesVector, totalDeflectionFlags, totalDerivativeFlags, 
					 totalPotentialFlags, pBaseLens, pSheetLens)))
		return "Can't start back projection matrix initialization: " + m_pTotalBPMatrix->getErrorString();

	if (shortImagesVector.size() > 0) // Ok, can improve speed in mass-scale calculation
	{
		std::vector<bool> shortDeflectionFlags, shortDerivativeFlags, shortPotentialFlags;

		fitnessObject.getShortCalcFlags(shortDeflectionFlags, shortDerivativeFlags, shortPotentialFlags);

		m_pShortBPMatrix = make_shared<BackProjectMatrix>();
		if (!(r = m_pShortBPMatrix->startInit(z_d, basisLenses[0].first->getLensDistance(), m_pDeflectionMatrix.get(),
						 shortImagesVector, shortDeflectionFlags, shortDerivativeFlags, 
						 shortPotentialFlags, pBaseLens, pSheetLens)))
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
	
	m_pDeflectionMatrix->calculateBasisMatrixProducts(basisFunctionWeights, true, true, true);
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

	if (!fitnessFunction.calculateOverallFitness(*m_pTotalBPMatrix, pFitnessValues))
		return "Unable to calculate full fitness: " + fitnessFunction.getErrorString();

	return true;
}

} // end namespace
