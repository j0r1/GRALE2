#include "lensinversiongafactorycommon.h"
#include "lensinversiongenome.h"
#include "configurationparameters.h"
#include "utils.h"
#include <limits>
#include <iostream>
#include <sstream>

#include "debugnew.h"

using namespace std;

namespace grale
{

LensInversionGAFactoryCommon::LensInversionGAFactoryCommon()
{
  	m_numBasisFunctions = 0;
	m_numSheetValues = 0;
	m_maxGenerations = 0;
	m_allowNegativeValues = false;
}

LensInversionGAFactoryCommon::~LensInversionGAFactoryCommon()
{

}

bool LensInversionGAFactoryCommon::setCommonParameters(int numSheetValues, 
							int maxGenerations,
							bool allowNeg,
							const std::vector<double> &basisFunctionMasses,
							double massUnit, double targetMass,
							const ScaleSearchParameters &searchParams)
{
	sendMessage("RNG SEED: " + std::to_string(m_rndGen.getSeed()));

	m_numBasisFunctions = (int)basisFunctionMasses.size();
	if (m_numBasisFunctions == 0)
	{
		setErrorString("No basis function masses specified");
		return false;
	}

	m_numSheetValues = numSheetValues;
	m_maxGenerations = maxGenerations;
	m_allowNegativeValues = allowNeg;

	m_basisFunctionMasses.clear();
	for (auto m : basisFunctionMasses)
		m_basisFunctionMasses.push_back((float)(m/massUnit));

	m_targetMass = (float)(targetMass/massUnit);

	m_massScaleSearchParams = searchParams;
	if (m_massScaleSearchParams == ScaleSearchParameters(true))
		sendMessage("Using wide scale factor search");
	else if (m_massScaleSearchParams == ScaleSearchParameters(false))
		sendMessage("Using normal scale factor search");
	else if (m_massScaleSearchParams.getNumberOfIterations() <= 0)
		sendMessage("Not using any extra mass scale search");
	else
		sendMessage("Using custom scaling parameters: " + m_massScaleSearchParams.toString());

	return true;
}

bool LensInversionGAFactoryCommon::initializeLensFitnessObject(double z_d,
	const std::vector<std::shared_ptr<ImagesDataExtended>> &images,
	const ConfigurationParameters *pFitnessObjectParams,
	vector<ImagesDataExtended*> &reducedImagesVector,
	vector<ImagesDataExtended*> &shortImagesVector)
{
	list<ImagesDataExtended *> reducedImages;
	list<ImagesDataExtended *> shortImages;

	for (auto i : images)
		reducedImages.push_back(i.get());

	auto pFitnessObject = createFitnessObject();
	if (pFitnessObject == 0)
		return false; // error string should be set in the createFitnessObject function
	m_fitnessObject.reset(pFitnessObject);

	ConfigurationParameters fitnesObjectParams;
	if (pFitnessObjectParams)
		fitnesObjectParams = *pFitnessObjectParams;

	fitnesObjectParams.clearRetrievalMarkers();
	if (!m_fitnessObject->init(z_d, reducedImages, shortImages, &fitnesObjectParams))
	{
		setErrorString(m_fitnessObject->getErrorString());
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
		return false;
	}

	shortImagesVector.clear();
	reducedImagesVector.clear();

	// Transform from list to vector
	reducedImagesVector.insert(reducedImagesVector.end(), reducedImages.begin(), reducedImages.end());
	shortImagesVector.insert(shortImagesVector.end(), shortImages.begin(), shortImages.end());

	// TODO: for now this has to be done in a subclass, as the actual factory
	//       may have either a single objective GA as parent, or a multi-objective
	//       one. The single objecive factory doesn't have the member function to
	//       set the number of fitness components
	//setNumberOfFitnessComponents(m_pFitnessObject->getNumberOfFitnessComponents());

	// This should at least perform the setNumberOfFitnessComponents call
	if (!subInit(m_fitnessObject.get()))
		return false;

	return true;
}


void LensInversionGAFactoryCommon::sendMessage(const std::string &s)
{
	if (getCurrentAlgorithm())
		GAFactory::sendMessage(s);
	else // queue for later
	{
		cerr << "Queue: " << s << endl;
		m_queuedMessages.push_back(s);
	}
}

void LensInversionGAFactoryCommon::onGeneticAlgorithmStart()
{
	// Send the messages that were previously queued
	for (auto &s : m_queuedMessages)
		GAFactory::sendMessage(s);

	m_queuedMessages.clear();

	// To debug the mass scale search
	string fileName;
	if (getenv("GRALE_DEBUG_MASSCALESEARCH", fileName))
	{
		m_scaleSearchFileStream.open(fileName, ios_base::out);
		if (!m_scaleSearchFileStream.is_open())
			sendMessage("WARNING: couldn't open 'GRALE_DEBUG_MASSCALESEARCH' debug file '" + fileName + "'");
		else
		{
			m_scaleSearchFileStream.precision(10);
			m_scaleSearchFileStream << "# " << m_massScaleSearchParams.toString() << endl;
		}
	}
}

mogal::Genome *LensInversionGAFactoryCommon::createNewGenome() const
{
	assert(m_numBasisFunctions > 0);
	return new LensInversionGenome(const_cast<LensInversionGAFactoryCommon*>(this), m_numBasisFunctions, m_numSheetValues);
}

bool LensInversionGAFactoryCommon::writeGenome(serut::SerializationInterface &si, const mogal::Genome *g) const
{
	const LensInversionGenome *g2 = (const LensInversionGenome *)g;

	if (!si.writeFloats(g2->getBasisFunctionWeights()))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeFloats(g2->getSheetValues()))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool LensInversionGAFactoryCommon::readGenome(serut::SerializationInterface &si, mogal::Genome **g) const
{
	assert(m_numBasisFunctions > 0);
	std::vector<float> basisFunctionWeights(m_numBasisFunctions);
	std::vector<float> sheetValues(m_numSheetValues);

	if (!si.readFloats(basisFunctionWeights))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.readFloats(sheetValues))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	// This moves weights and sheetValues
	*g = new LensInversionGenome(const_cast<LensInversionGAFactoryCommon*>(this), basisFunctionWeights, sheetValues);
	return true;
}

bool LensInversionGAFactoryCommon::writeGenomeFitness(serut::SerializationInterface &si, const mogal::Genome *g) const
{
	const LensInversionGenome *g2 = (const LensInversionGenome *)g;
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
	return true;
}

bool LensInversionGAFactoryCommon::readGenomeFitness(serut::SerializationInterface &si, mogal::Genome *g) const
{
	LensInversionGenome *g2 = (LensInversionGenome *)g;
	int num = getNumberOfFitnessComponents();
	float x[GRIDLENSINVERSIONGENOMEBASE_MAXFITNESSCOMP+1];

	if (!si.readFloats(x, num+1))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	g2->setFitnessValues(x);
	g2->setScaleFactor(x[num]);
	return true;
}

bool LensInversionGAFactoryCommon::writeCommonGenerationInfo(serut::SerializationInterface &si) const
{
	return true;
}

bool LensInversionGAFactoryCommon::readCommonGenerationInfo(serut::SerializationInterface &si)
{
	return true;
}

void LensInversionGAFactoryCommon::onCurrentBest(const list<mogal::Genome *> &bestGenomes)
{
	stringstream ss;
	ss << "Current best:";
	for (auto g : bestGenomes)
		ss << "( " << g->getFitnessDescription() << ")";

	sendMessage(ss.str());
}

static float LogTrans(float x) { return LN(x); }
static float ExpTrans(float x) { return EXP(x); }
static float IdentityTrans(float x) { return x; }

bool LensInversionGAFactoryCommon::calculateFitness(const vector<float> &basisFunctionWeights,
													const vector<float> &sheetValues,
													float &scaleFactor,
													float *pFitnessValues
													)
{
	if (m_scaleSearchFileStream.is_open()) // For debugging
	{
		m_scaleSearchStringStream.clear();
		m_scaleSearchStringStream.precision(10);
		m_searchedPoints.clear();
	}

	int numBasisFunctions = basisFunctionWeights.size();

	if (!initializeNewCalculation(basisFunctionWeights, sheetValues))
		return false;

	// TODO: adjust mass scale depending on sheet value?

	float massum = 0;
	if (allowNegativeValues())
	{
		// TODO: is this enough? or do we need some additional constraints?
		// TODO: should find something better!!
		for (int i = 0 ; i < numBasisFunctions ; i++)
			massum += ABS(basisFunctionWeights[i]*m_basisFunctionMasses[i]);
	}
	else
	{
		for (int i = 0 ; i < numBasisFunctions ; i++)
			massum += basisFunctionWeights[i]*m_basisFunctionMasses[i];
	}

	float startValue = m_massScaleSearchParams.getStartFactor();
	float stopValue = m_massScaleSearchParams.getStopFactor();
	int numiterations = m_massScaleSearchParams.getNumberOfIterations();
	int numiterationsteps = m_massScaleSearchParams.getStepsOnFirstIteration();
	int numiterationsteps2 = m_massScaleSearchParams.getStepsOnSubsequentIterations();

	// we want a scale factor of 1 to correspond to the total mass specified by the
	// massscale used in the BackProjectMatrix

	startValue /= (massum/m_targetMass);
	stopValue /= (massum/m_targetMass);

	// Forward transform
	float (*FT)(float x) = IdentityTrans;
	// Inverse transform
	float (*IT)(float x) = IdentityTrans;

	if (useLogarithmicScaleSearch())
	{
		FT = LogTrans;
		IT = ExpTrans;
	}

	startValue = FT(startValue);
	stopValue = FT(stopValue);

	// Store these start values so that we don't do out of bounds
	float startValue0 = startValue;
	float stopValue0 = stopValue;

	if (numiterations > 0)
	{
		float currentBestFitness = std::numeric_limits<float>::max();
		float currentBestScaleFactor = 1.0f;

 		for (int i = 0 ; i < numiterations ; i++)
		{
			float stepsize = (stopValue-startValue)/((float)(numiterationsteps-1));
			float s = startValue;

			for (int j = 0 ; j < numiterationsteps ; j++, s += stepsize)
			{
				float realScale = IT(s);
				float f;
				
				if (!calculateMassScaleFitness(realScale, f))
					return false;

				if (f < currentBestFitness)
				{
					currentBestFitness = f;
					currentBestScaleFactor = s;
				}

				if (m_scaleSearchFileStream.is_open()) // To debug/illustrate the scale search
					m_searchedPoints.push_back({ realScale, f});

			}

			startValue = currentBestScaleFactor-stepsize;
			stopValue = currentBestScaleFactor+stepsize;

			// Make sure we stay within bounds
			if (startValue < startValue0)
				startValue = startValue0;
			if (stopValue > stopValue0)
				stopValue = stopValue0;

			// After the first loop, we're going to just zoom in on the first located value
			// This can use a different number of steps
			numiterationsteps = numiterationsteps2;
		}

		scaleFactor = IT(currentBestScaleFactor);
	}
	else // numiterations <= 0, no further scale factor determination requested
	{
		scaleFactor = 1.0f;
	}

	// Do the final (possibly multi-component) fitness evaluation
	if (!calculateTotalFitness(scaleFactor, pFitnessValues))
		return false;

	// To debug the mass scale search
	if (m_scaleSearchFileStream.is_open())
	{
		std::sort(m_searchedPoints.begin(), m_searchedPoints.end(), [](auto x, auto y) { return x.first < y.first; });
		m_scaleSearchStringStream << "{" << endl;
		m_scaleSearchStringStream << " 'scaleFactor': " << scaleFactor << "," << endl;
		m_scaleSearchStringStream << " 'startValue0': " << IT(startValue0) << "," << endl;
		m_scaleSearchStringStream << " 'stopValue0': " << IT(stopValue0) << "," << endl;
		m_scaleSearchStringStream << " 'startFactor': " << m_massScaleSearchParams.getStartFactor() << "," << endl;
		m_scaleSearchStringStream << " 'stopFactor': " << m_massScaleSearchParams.getStopFactor() << "," << endl;
		m_scaleSearchStringStream << " 'scaleSearch': [ " << endl;
		for (auto x : m_searchedPoints)
			m_scaleSearchStringStream << "[ " << x.first << ", " << x.second << " ]," << endl;
		m_scaleSearchStringStream << " ]" << endl;
		m_scaleSearchStringStream << "}," << endl;
	}
	return true;
}

void LensInversionGAFactoryCommon::onSortedPopulation(const std::vector<mogal::GenomeWrapper> &population)
{
	if (m_scaleSearchFileStream.is_open())
	{
		m_scaleSearchFileStream << "[" << endl;
		m_scaleSearchFileStream << m_scaleSearchStringStream.str() << endl;
		m_scaleSearchFileStream << "]," << endl;
	}
}

} // end namespace