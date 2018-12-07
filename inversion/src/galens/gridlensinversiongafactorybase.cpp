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
#include "plummerlens.h"
#include "squarelens.h"
#include "gausslens.h"
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
	m_pCurrentParams = 0;
	m_maxGenerations = 0;
	m_pDeflectionMatrix = 0;
	m_pShortBPMatrix = 0;
	m_pTotalBPMatrix = 0;
	m_pFitnessObject = 0;
	m_wideSearch = false;
}

GridLensInversionGAFactoryBase::~GridLensInversionGAFactoryBase()
{
	if (m_pCurrentParams)
		delete m_pCurrentParams;
	for (int i = 0 ; i < m_basisLenses.size() ; i++)
		delete m_basisLenses[i].first;
	if (m_pDeflectionMatrix)
			delete m_pDeflectionMatrix;
	if (m_pShortBPMatrix)
	{
		if (m_pShortBPMatrix != m_pTotalBPMatrix)
			delete m_pShortBPMatrix;
	}
	if (m_pTotalBPMatrix)
		delete m_pTotalBPMatrix;
	if (m_pFitnessObject)
		delete m_pFitnessObject;
}

mogal::GAFactoryParams *GridLensInversionGAFactoryBase::createParamsInstance() const
{
	return new GridLensInversionGAFactoryParams();
}

const mogal::GAFactoryParams *GridLensInversionGAFactoryBase::getCurrentParameters() const
{
	return m_pCurrentParams;
}

void GridLensInversionGAFactoryBase::sendMessage(const std::string &s)
{
	if (getCurrentAlgorithm())
		GAFactory::sendMessage(s);
	else // queue for later
	{
		cerr << "Queue: " << s << endl;
		m_queuedMessages.push_back(s);
	}
}

void GridLensInversionGAFactoryBase::onGeneticAlgorithmStart()
{
	// Send the messages that were previously queued
	for (auto &s : m_queuedMessages)
		GAFactory::sendMessage(s);

	m_queuedMessages.clear();
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
	
	const LensInversionGAFactoryParams *p0 = (const LensInversionGAFactoryParams *)p;
	if (p0->getInversionType() != LensInversionGAFactoryParams::GridInversion)
	{
		setErrorString("This algorithm is not a grid based one");
		return false;
	}

	sendMessage("RNG SEED: " + std::to_string(m_rndGen.getSeed()));

	const GridLensInversionGAFactoryParams *p2 = (const GridLensInversionGAFactoryParams *)p;
	
	m_pCurrentParams = p2->createCopy();
	assert(m_pCurrentParams);

	// Determine mass weights

	const std::vector<GridSquare> &squares = m_pCurrentParams->getGridSquares();
	m_massWeights.resize(squares.size());
	m_numMasses = m_pCurrentParams->getGridSquares().size();

	if (p2->useMassWeights())
	{
		sendMessage("Using mass weights");
		
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

			m_massWeights[i] = (float)((((double)m_numMasses)*w2)/sum);
//			std::cout << " " << m_massWeights[i];
		}
		//std::cout << std::endl;
	}
	else
	{
		int i = 0;

		for (auto it = squares.begin(); it != squares.end() ; it++, i++)
			m_massWeights[i] = 1;
	}

	// Create the basis functions

	std::vector<std::pair<GravitationalLens *, Vector2D<double> > > basisLenses;
	double massScale = p2->getMassScale();
	double D_d = p2->getD_d();
	int squareNumber = 0;

	for (auto squareIt = squares.begin() ; squareIt != squares.end() ; squareIt++, squareNumber++)
	{
		double gridSize = (*squareIt).getSize();
		Vector2D<double> gridCenter = (*squareIt).getCenter();
		bool returnValue;
		GravitationalLens *pLens;
	
		if (m_pCurrentParams->getBasisFunctionType() == GridLensInversionParameters::PlummerBasis)
		{
			PlummerLensParams lensParams(massScale * m_massWeights[squareNumber], gridSize);

			pLens = new PlummerLens();
			returnValue = pLens->init(D_d, &lensParams);
		}
		else if (m_pCurrentParams->getBasisFunctionType() == GridLensInversionParameters::GaussBasis)
		{
			GaussLensParams lensParams(massScale * m_massWeights[squareNumber], gridSize);

			pLens = new GaussLens();
			returnValue = pLens->init(D_d, &lensParams);
		}
		else // Squares
		{
			SquareLensParams lensParams(massScale * m_massWeights[squareNumber], gridSize);

			pLens = new SquareLens();
			returnValue = pLens->init(D_d, &lensParams);
		}

		if (!returnValue)
		{
			setErrorString(std::string("Couldn't initialize basis lens: ") + pLens->getErrorString());
			for (int i = 0 ; i < m_basisLenses.size() ; i++)
				delete m_basisLenses[i].first;
			m_basisLenses.clear();
			delete pLens;
			delete m_pCurrentParams;
			m_pCurrentParams = 0;
			return false;
		}

		m_basisLenses.push_back(std::pair<GravitationalLens *, Vector2D<double> >(pLens, gridCenter));
	}

	bool useSheet = false;

	if (p2->getMassSheetSearchType() == GridLensInversionParameters::Genome)
	{
		useSheet = true;
		// TODO: what is a good value here?
		m_sheetScale = 1.1;
		m_useGenomeSheet = true;
	}
	else // no sheet
	{
		m_sheetScale = 0;
		m_useGenomeSheet = false;
	}

	// perform sub-initialization
	if (!localSubInit(p2->getZ_d(), p2->getImages(), m_basisLenses, p2->getBaseLens(), useSheet, 
				m_pCurrentParams->getFitnessObjectParameters()))
	{
		// Error string is set in subInit
		for (int i = 0 ; i < m_basisLenses.size() ; i++)
			delete m_basisLenses[i].first;
		m_basisLenses.clear();
		delete m_pCurrentParams;
		m_pCurrentParams = 0;
		return false;
	}
	
	m_maxGenerations = m_pCurrentParams->getMaximumNumberOfGenerations();
	m_allowNegativeValues = p2->allowNegativeValues();
	m_wideSearch = p2->useWideSearch();

	if (m_wideSearch)
		sendMessage("Using wide scale factor search");
	else
		sendMessage("Using normal scale factor search");

	return true;
}

mogal::Genome *GridLensInversionGAFactoryBase::createNewGenome() const
{
	return new GridLensInversionGenomeBase((GridLensInversionGAFactoryBase *)this, m_numMasses, m_useGenomeSheet);
}

bool GridLensInversionGAFactoryBase::writeGenome(serut::SerializationInterface &si, const mogal::Genome *g) const
{
	const GridLensInversionGenomeBase *g2 = (const GridLensInversionGenomeBase *)g;

	if (!si.writeFloats(g2->getMasses()))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeFloat(g2->getSheetValue()))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool GridLensInversionGAFactoryBase::readGenome(serut::SerializationInterface &si, mogal::Genome **g) const
{
	std::vector<float> masses(m_numMasses);
	float sheetValue;

	if (!si.readFloats(masses))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.readFloat(&sheetValue))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	*g = new GridLensInversionGenomeBase((GridLensInversionGAFactoryBase *)this, masses, sheetValue);
	return true;
}

bool GridLensInversionGAFactoryBase::writeGenomeFitness(serut::SerializationInterface &si, const mogal::Genome *g) const
{
	const GridLensInversionGenomeBase *g2 = (const GridLensInversionGenomeBase *)g;
	int num = getNumberOfFitnessComponents();
	const float *f = g2->getFitnessValues();

	if (!si.writeFloats(f, num))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeFloat(g2->getScaleFactor()))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeFloat(g2->getSheetFactor()))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool GridLensInversionGAFactoryBase::readGenomeFitness(serut::SerializationInterface &si, mogal::Genome *g) const
{
	GridLensInversionGenomeBase *g2 = (GridLensInversionGenomeBase *)g;
	int num = getNumberOfFitnessComponents();
	float x[GRIDLENSINVERSIONGENOMEBASE_MAXFITNESSCOMP+2];

	if (!si.readFloats(x, num+2))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	g2->setFitnessValues(x);
	g2->setScaleFactor(x[num]);
	g2->setSheetFactor(x[num+1]);
	return true;
}

bool GridLensInversionGAFactoryBase::writeCommonGenerationInfo(serut::SerializationInterface &si) const
{
	return true;
}

bool GridLensInversionGAFactoryBase::readCommonGenerationInfo(serut::SerializationInterface &si)
{
	return true;
}

GravitationalLens *GridLensInversionGAFactoryBase::createLens(const std::vector<float> &masses, float sheetValue, float scaleFactor, 
                                                              double *pTotalMass, std::string &errStr) const
{
	CompositeLensParams lensParams;
	CompositeLens *pLens;

	for (int i = 0 ; i < m_basisLenses.size() ; i++)
		lensParams.addLens((double)masses[i]*(double)scaleFactor, m_basisLenses[i].second, 0, *(m_basisLenses[i].first));

	if (m_sheetScale != 0)
	{
		double D_d = m_basisLenses[0].first->getLensDistance();
		MassSheetLensParams sheetParams(((SPEED_C*SPEED_C)/(4.0*CONST_PI*CONST_G*D_d))*sheetValue*m_sheetScale);
		MassSheetLens sheetLens;

		if (!sheetLens.init(D_d, &sheetParams))
		{
			setErrorString(std::string("Can't initialize mass sheet component: ") + sheetLens.getErrorString());
			return 0;
		}

		lensParams.addLens(1.0, Vector2D<double>(0,0), 0, sheetLens);
	}

	pLens = new CompositeLens();
	pLens->init(m_basisLenses[0].first->getLensDistance(), &lensParams);

	// TODO
	*pTotalMass = 0;

	// TODO: base lens!?
	return pLens;
}

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

bool GridLensInversionGAFactoryBase::localSubInit(double z_d, const std::vector<ImagesDataExtended *> &images, 
	                  const std::vector<std::pair<GravitationalLens *, Vector2D<double> > > &basisLenses,
                      const GravitationalLens *pBaseLens, bool useSheet, 
					  const ConfigurationParameters *pFitnessObjectParams)
{
	std::list<ImagesDataExtended *> reducedImages;
	
	for (auto i :images)
		reducedImages.push_back(i);

	std::list<ImagesDataExtended *> shortImages;

	m_pFitnessObject = createFitnessObject();
	if (m_pFitnessObject == 0)
		return false; // error string should be set in the createFitnessObject function

	ConfigurationParameters fitnesObjectParams;
	if (pFitnessObjectParams)
		fitnesObjectParams = *pFitnessObjectParams;

	fitnesObjectParams.clearRetrievalMarkers();

	if (!m_pFitnessObject->init(z_d, reducedImages, shortImages, &fitnesObjectParams))
	{
		setErrorString(m_pFitnessObject->getErrorString());
		delete m_pFitnessObject;
		m_pFitnessObject = 0;
		return false;
	}

	vector<string> unusedKeys;
	fitnesObjectParams.getUnretrievedKeys(unusedKeys);

	if (unusedKeys.size() > 0)
	{
		stringstream ss;

		ss << "Some parameters that were specified for the lens fitness object were not used:";
		for (auto &k : unusedKeys)
			ss << " " << k;
		setErrorString(ss.str());

		delete m_pFitnessObject;
		m_pFitnessObject = 0;
		return false;
	}

	std::vector<ImagesDataExtended *> reducedImagesVector, shortImagesVector;
	std::list<ImagesDataExtended *>::const_iterator it;

	for (it = reducedImages.begin() ; it != reducedImages.end() ; it++)
		reducedImagesVector.push_back(*it);

	for (it = shortImages.begin() ; it != shortImages.end() ; it++)
		shortImagesVector.push_back(*it);
	
	// TODO: this has to be done in a subclass!
	//setNumberOfFitnessComponents(m_pFitnessObject->getNumberOfFitnessComponents());

	// This should at least perform the setNumberOfFitnessComponents call
	if (!subInit(m_pFitnessObject))
	{
		delete m_pFitnessObject;
		m_pFitnessObject = 0;
		return false;
	}

	m_pDeflectionMatrix = new GADeflectionMatrix(this);
	if (!m_pDeflectionMatrix->startInit())
	{
		setErrorString(m_pDeflectionMatrix->getErrorString());
		delete m_pDeflectionMatrix;
		delete m_pFitnessObject;
		m_pDeflectionMatrix = 0;
		m_pFitnessObject = 0;
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
					 useSheet))
	{
		setErrorString(m_pTotalBPMatrix->getErrorString());
		delete m_pTotalBPMatrix;
		delete m_pDeflectionMatrix;
		delete m_pFitnessObject;
		m_pDeflectionMatrix = 0;
		m_pFitnessObject = 0;
		m_pTotalBPMatrix = 0;
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
						 useSheet))
		{
			setErrorString(m_pShortBPMatrix->getErrorString());
			delete m_pShortBPMatrix;
			delete m_pTotalBPMatrix;
			delete m_pDeflectionMatrix;
			delete m_pFitnessObject;
			m_pDeflectionMatrix = 0;
			m_pFitnessObject = 0;
			m_pTotalBPMatrix = 0;
			m_pShortBPMatrix = 0;
			return false;
		}

	}
	else
		m_pShortBPMatrix = m_pTotalBPMatrix;

	bool err = false;

	if (!m_pDeflectionMatrix->endInit(basisLenses))
	{
		setErrorString(m_pDeflectionMatrix->getErrorString());
		err = true;
	}
	if (!err)
	{
		if (!m_pTotalBPMatrix->endInit())
		{
			setErrorString(m_pTotalBPMatrix->getErrorString());
			err = true;
		}
	}
	if (!err)
	{
		if (m_pShortBPMatrix != m_pTotalBPMatrix)
		{
			if (!m_pShortBPMatrix->endInit())
			{
				setErrorString(m_pShortBPMatrix->getErrorString());
				err = true;
			}
		}
	}

	if (err)
	{
		if (m_pShortBPMatrix != m_pTotalBPMatrix)
		{
			delete m_pShortBPMatrix;
			m_pShortBPMatrix = 0;
		}
		delete m_pTotalBPMatrix;
		delete m_pDeflectionMatrix;
		delete m_pFitnessObject;
		m_pDeflectionMatrix = 0;
		m_pFitnessObject = 0;
		m_pTotalBPMatrix = 0;
		return false;
	}

	m_pFitnessObject->postInit(reducedImages, shortImages, m_pDeflectionMatrix->getAngularScale());

	return true;
}

void GridLensInversionGAFactoryBase::getGenomeCalculationParameters(float &startfactor, float &stopfactor, int &numiterationsteps, int &numiterations, int &numiterationsteps2) const
{ 
	if (!m_wideSearch)
	{
		startfactor = 0.2;
		stopfactor = 5.0;
		numiterationsteps = 20;
		numiterationsteps2 = 20;
		numiterations = 5;
	}
	else
	{
		startfactor = 0.01;
		stopfactor = 100.0;
		numiterationsteps = 50;
		numiterationsteps2 = 5;
		numiterations = 11;
	}
}

void GridLensInversionGAFactoryBase::getGenomeSheetCalculationParameters(float &startfactor, float &stopfactor) const
{ 
	startfactor = 0;
	stopfactor = 2.0;
}

void GridLensInversionGAFactoryBase::onCurrentBest(const list<mogal::Genome *> &bestGenomes)
{
	stringstream ss;
	ss << "Current best:";
	for (auto g : bestGenomes)
		ss << "( " << g->getFitnessDescription() << ")";

	sendMessage(ss.str());
}

} // end namespace
