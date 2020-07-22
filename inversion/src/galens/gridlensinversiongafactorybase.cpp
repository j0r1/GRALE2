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

#include "gridlensinversiongafactorybase.h"
#include "gridlensinversiongafactoryparams.h"
#include "gridlensinversiongenomebase.h"
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

#include "debugnew.h"

using namespace std;

namespace grale
{

class GADeflectionMatrix : public DeflectionMatrix
{
public:
	GADeflectionMatrix(GridLensInversionGAFactoryBase *pFactory) : m_pFactory(pFactory)
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

	GridLensInversionGAFactoryBase *m_pFactory;
};

GridLensInversionGAFactoryBase::GridLensInversionGAFactoryBase()
{
	zero();
}

void GridLensInversionGAFactoryBase::zero()
{
	m_pCurrentParams = nullptr;
	m_pDeflectionMatrix = nullptr;
	m_pShortBPMatrix = nullptr;
	m_pTotalBPMatrix = nullptr;
	m_pFitnessObject = nullptr;
}

void GridLensInversionGAFactoryBase::clear()
{
	delete m_pCurrentParams;
	delete m_pDeflectionMatrix;
	if (m_pShortBPMatrix == m_pTotalBPMatrix)
		m_pShortBPMatrix = nullptr; // Avoid double deletion
	delete m_pShortBPMatrix;
	delete m_pTotalBPMatrix;
	delete m_pFitnessObject;

	m_basisLenses.clear();
	m_sheetLens = nullptr;
	m_queuedMessages.clear();

	zero();
}

GridLensInversionGAFactoryBase::~GridLensInversionGAFactoryBase()
{
	clear();
}

mogal::GAFactoryParams *GridLensInversionGAFactoryBase::createParamsInstance() const
{
	return new GridLensInversionGAFactoryParams();
}

const mogal::GAFactoryParams *GridLensInversionGAFactoryBase::getCurrentParameters() const
{
	return m_pCurrentParams;
}

bool GridLensInversionGAFactoryBase::init(const mogal::GAFactoryParams *p)
{
	if (m_pCurrentParams != 0)
	{
		setErrorString("Already initialized");
		return false;
	}

	if (p == 0)
	{
		setErrorString("Specified parameters can't be null");
		return false;
	}
	
	const GridLensInversionGAFactoryParams *p2 = dynamic_cast<const GridLensInversionGAFactoryParams *>(p);
	if (!p2)
	{
		setErrorString("Invalid type of GA factory parameters");
		return false;
	}

	m_pCurrentParams = p2->createCopy();
	assert(m_pCurrentParams);

	// Determine mass weights

	auto basisLenses = m_pCurrentParams->getBasisLenses();
	if (basisLenses.size() == 0)
	{
		setErrorString("Unable to get basis lenses: " + m_pCurrentParams->getErrorString());
		return false;
	}

	int numMasses = basisLenses.size();

	double massScale = m_pCurrentParams->getMassScale();
	double totalRelevantLensingMass = 0;
	for (const auto &bl : basisLenses)
	{
		m_basisLenses.push_back( { bl.m_pLens, bl.m_center });
		
		double m = bl.m_relevantLensingMass;
		if (m <= 0) // TODO: allow this to be overridden
		{
			setErrorString("A basis lens has a negative lensing mass");
			clear();
			return false;
		}
	}

	vector<float> massWeights; // masses expressed in units of massScale
	for (const auto &bl : basisLenses)
		massWeights.push_back((float)(bl.m_relevantLensingMass/massScale));

	const GravitationalLens *pSheetLens = p2->getSheetLens();
	int numSheetValues = 0;
	if (pSheetLens)
	{
		m_sheetLens = shared_ptr<GravitationalLens>(pSheetLens->createCopy());
		numSheetValues = 1;
	}
	else // no sheet
		m_sheetLens = nullptr;

	// perform sub-initialization
	if (!localSubInit(p2->getZ_d(), p2->getImages(), m_basisLenses, p2->getBaseLens(), pSheetLens, 
				m_pCurrentParams->getFitnessObjectParameters()))
	{
		// Error string is set in subInit
		clear();
		return false;
	}
	
	int maxGenerations = m_pCurrentParams->getMaximumNumberOfGenerations();
	bool allowNegativeValues = p2->allowNegativeValues();
	auto massScaleSearchParams = p2->getMassScaleSearchParameters();

	if (!setCommonParameters(numMasses, numSheetValues, maxGenerations, allowNegativeValues,
	                         m_pDeflectionMatrix->getAngularScale(), massWeights,
							 massScaleSearchParams))
	{
		clear();
		// Error string was already set
		return false;
	}

	if (massScaleSearchParams == ScaleSearchParameters(true))
		sendMessage("Using wide scale factor search");
	else if (massScaleSearchParams == ScaleSearchParameters(false))
		sendMessage("Using normal scale factor search");
	else if (massScaleSearchParams.getNumberOfIterations() <= 0)
		sendMessage("Not using any extra mass scale search");
	else
		sendMessage("Using custom scaling parameters: " + massScaleSearchParams.toString());

	return true;
}

GravitationalLens *GridLensInversionGAFactoryBase::createLens(const GridLensInversionGenomeBase &genome,
                                                              std::string &errStr) const
{
	CompositeLensParams lensParams;
	CompositeLens *pLens;

	const vector<float> &masses = genome.getMasses();
	const float scaleFactor = genome.getScaleFactor();
	const vector<float> &sheetValues = genome.getSheetValues();

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
		lensParams.addLens((double)masses[i]*(double)scaleFactor, m_basisLenses[i].second, 0, *m_basisLenses[i].first.get());

	if (m_sheetLens.get())
		lensParams.addLens(sheetValue, Vector2D<double>(0,0), 0, *m_sheetLens.get());

	pLens = new CompositeLens();
	pLens->init(m_basisLenses[0].first->getLensDistance(), &lensParams);

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

bool GridLensInversionGAFactoryBase::localSubInit(double z_d, const vector<shared_ptr<ImagesDataExtended>> &images, 
	                  const vector<pair<shared_ptr<GravitationalLens>, Vector2D<double> > > &basisLenses,
                      const GravitationalLens *pBaseLens, const GravitationalLens *pSheetLens, 
					  const ConfigurationParameters *pFitnessObjectParams)
{
	std::list<ImagesDataExtended *> reducedImages;
	
	for (auto i :images)
		reducedImages.push_back(i.get());

	std::list<ImagesDataExtended *> shortImages;

	m_pFitnessObject = createFitnessObject();
	if (m_pFitnessObject == 0)
	{
		clear();
		return false; // error string should be set in the createFitnessObject function
	}

	ConfigurationParameters fitnesObjectParams;
	if (pFitnessObjectParams)
		fitnesObjectParams = *pFitnessObjectParams;

	fitnesObjectParams.clearRetrievalMarkers();

	if (!m_pFitnessObject->init(z_d, reducedImages, shortImages, &fitnesObjectParams))
	{
		setErrorString(m_pFitnessObject->getErrorString());
		clear();
		return false;
	}

	// Check that all keys in the fitness object parameters are actually used
	vector<string> unusedKeys;
	fitnesObjectParams.getUnretrievedKeys(unusedKeys);

	if (unusedKeys.size() > 0)
	{
		stringstream ss;

		ss << "Some parameters that were specified for the lens fitness object were not used:";
		for (auto &k : unusedKeys)
			ss << " " << k;
		setErrorString(ss.str());

		clear();
		return false;
	}

	// Transform from list to vector
	std::vector<ImagesDataExtended *> reducedImagesVector { reducedImages.begin(), reducedImages.end() };
	std::vector<ImagesDataExtended *> shortImagesVector { shortImages.begin(), shortImages.end() };
	
	// TODO: this has to be done in a subclass!
	//setNumberOfFitnessComponents(m_pFitnessObject->getNumberOfFitnessComponents());

	// This should at least perform the setNumberOfFitnessComponents call
	if (!subInit(m_pFitnessObject))
	{
		clear();
		return false;
	}

	m_pDeflectionMatrix = new GADeflectionMatrix(this);
	if (!m_pDeflectionMatrix->startInit())
	{
		setErrorString(m_pDeflectionMatrix->getErrorString());
		clear();
		return false;
	}

	bool totalStoreIntens, totalStoreTimeDelay, totalStoreShearInfo;
	std::vector<bool> totalDeflectionFlags, totalDerivativeFlags, totalPotentialFlags;

	m_pFitnessObject->getTotalCalcFlags(totalDeflectionFlags, totalDerivativeFlags, totalPotentialFlags);
	m_pFitnessObject->getTotalStoreFlags(&totalStoreIntens, &totalStoreTimeDelay, &totalStoreShearInfo);

	m_pTotalBPMatrix = new BackProjectMatrixNew();
	if (!m_pTotalBPMatrix->startInit(z_d, basisLenses[0].first->getLensDistance(), m_pDeflectionMatrix,
				         reducedImagesVector, totalDeflectionFlags, totalDerivativeFlags, 
					 totalPotentialFlags, pBaseLens, totalStoreIntens, totalStoreTimeDelay,
					 totalStoreShearInfo,
					 pSheetLens))
	{
		setErrorString(m_pTotalBPMatrix->getErrorString());
		clear();
		return false;
	}

	if (shortImagesVector.size() > 0) // Ok, can improve speed in mass-scale calculation
	{
		bool shortStoreIntens, shortStoreTimeDelay, shortStoreShearInfo;
		std::vector<bool> shortDeflectionFlags, shortDerivativeFlags, shortPotentialFlags;

		m_pFitnessObject->getShortCalcFlags(shortDeflectionFlags, shortDerivativeFlags, shortPotentialFlags);
		m_pFitnessObject->getShortStoreFlags(&shortStoreIntens, &shortStoreTimeDelay, &shortStoreShearInfo);

		m_pShortBPMatrix = new BackProjectMatrixNew();
		if (!m_pShortBPMatrix->startInit(z_d, basisLenses[0].first->getLensDistance(), m_pDeflectionMatrix,
						 shortImagesVector, shortDeflectionFlags, shortDerivativeFlags, 
						 shortPotentialFlags, pBaseLens, shortStoreIntens, shortStoreTimeDelay,
						 shortStoreShearInfo,
						 pSheetLens))
		{
			setErrorString(m_pShortBPMatrix->getErrorString());
			clear();
			return false;
		}
	}
	else
		m_pShortBPMatrix = m_pTotalBPMatrix;

	bool err = false;

	if (!m_pDeflectionMatrix->endInit(basisLenses))
	{
		setErrorString(m_pDeflectionMatrix->getErrorString());
		clear();
		return false;
	}

	if (!m_pTotalBPMatrix->endInit())
	{
		setErrorString(m_pTotalBPMatrix->getErrorString());
		clear();
		return false;
	}

	if (m_pShortBPMatrix != m_pTotalBPMatrix)
	{
		if (!m_pShortBPMatrix->endInit())
		{
			setErrorString(m_pShortBPMatrix->getErrorString());
			clear();
			return false;
		}
	}

	m_pFitnessObject->postInit(reducedImages, shortImages, m_pDeflectionMatrix->getAngularScale());
	return true;
}

bool GridLensInversionGAFactoryBase::initializeNewCalculation(const vector<float> &masses, const vector<float> &sheetValues)
{
	if (sheetValues.size() == 0)
		m_sheetScale = 0;
	else if (sheetValues.size() == 1)
		m_sheetScale = sheetValues[0];
	else
	{
		cerr << "Unexpected: got more (" << sheetValues.size() << ") mass sheet contributions than expected" << endl;
		return false;
	}
	
	m_pDeflectionMatrix->calculateBasisMatrixProducts(masses, true, true, true);
	m_pTotalBPMatrix->storeDeflectionMatrixResults();
	m_pShortBPMatrix->storeDeflectionMatrixResults();
	return true;
}

bool GridLensInversionGAFactoryBase::calculateMassScaleFitness(float scaleFactor, float &fitness)
{
	LensFitnessObject &fitnessFunction = *m_pFitnessObject;

	m_pShortBPMatrix->calculate(scaleFactor, m_sheetScale);
	if (fitnessFunction.shortNeedInverseMagnifications())
		m_pShortBPMatrix->calculateInverseMagnifications(*(fitnessFunction.getShortInverseMagnificationFlags()));
	if (fitnessFunction.shortNeedShearComponents())
		m_pShortBPMatrix->calculateShearComponents(*(fitnessFunction.getShortShearComponentFlags()));
	if (fitnessFunction.shortNeedConvergence())
		m_pShortBPMatrix->calculateConvergence(*(fitnessFunction.getShortConvergenceFlags()));

	if (!fitnessFunction.calculateMassScaleFitness(*m_pShortBPMatrix, fitness))
	{
		cerr << "ERROR: Unable to calculate mass scale fitness: " << fitnessFunction.getErrorString() << endl;
		return false;
	}
	return true;
}

bool GridLensInversionGAFactoryBase::calculateTotalFitness(float scaleFactor, float *pFitnessValues)
{
	LensFitnessObject &fitnessFunction = *m_pFitnessObject;

	m_pTotalBPMatrix->calculate(scaleFactor, m_sheetScale);
	if (fitnessFunction.totalNeedInverseMagnifications())
		m_pTotalBPMatrix->calculateInverseMagnifications(*(fitnessFunction.getTotalInverseMagnificationFlags()));
	if (fitnessFunction.totalNeedShearComponents())
		m_pTotalBPMatrix->calculateShearComponents(*(fitnessFunction.getTotalShearComponentFlags()));
	if (fitnessFunction.totalNeedConvergence())
		m_pTotalBPMatrix->calculateConvergence(*(fitnessFunction.getTotalConvergenceFlags()));

	if (!fitnessFunction.calculateOverallFitness(*m_pTotalBPMatrix, pFitnessValues))
	{
		cerr << "ERROR: Unable to calculate full fitness: " << fitnessFunction.getErrorString() << endl;
		return false;
	}

	return true;
}

} // end namespace
