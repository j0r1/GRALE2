#include "lensinversiongafactorycommon.h"
#include "configurationparameters.h"
#include "utils.h"
#include <limits>
#include <iostream>
#include <sstream>

using namespace std;
using namespace errut;

namespace grale
{

LensInversionGAFactoryCommon::LensInversionGAFactoryCommon(unique_ptr<LensFitnessObject> fitObj)
	: m_fitnessObject(move(fitObj))
{
  	m_numBasisFunctions = 0;
	m_numSheetValues = 0;
	m_allowNegativeValues = false;
}

LensInversionGAFactoryCommon::~LensInversionGAFactoryCommon()
{

}

bool_t LensInversionGAFactoryCommon::setCommonParameters(int numSheetValues, 
							bool allowNeg,
							const std::vector<double> &basisFunctionMasses,
							double massUnit, double targetMass,
							const ScaleSearchParameters &searchParams)
{
	m_numBasisFunctions = (int)basisFunctionMasses.size();
	if (m_numBasisFunctions == 0)
		return "No basis function masses specified";

	m_numSheetValues = numSheetValues;
	m_allowNegativeValues = allowNeg;

	m_basisFunctionMasses.clear();
	for (auto m : basisFunctionMasses)
		m_basisFunctionMasses.push_back((float)(m/massUnit));

	m_targetMass = (float)(targetMass/massUnit);

	m_massScaleSearchParams = searchParams;
	if (m_massScaleSearchParams == ScaleSearchParameters(true))
		log("Using wide scale factor search");
	else if (m_massScaleSearchParams == ScaleSearchParameters(false))
		log("Using normal scale factor search");
	else if (m_massScaleSearchParams.getNumberOfIterations() <= 0)
		log("Not using any extra mass scale search");
	else
		log("Using custom scaling parameters: " + m_massScaleSearchParams.toString());

	return true;
}

bool_t LensInversionGAFactoryCommon::initializeLensFitnessObject(double z_d,
	const std::vector<std::shared_ptr<ImagesDataExtended>> &images,
	const ConfigurationParameters *pFitnessObjectParams,
	vector<ImagesDataExtended*> &reducedImagesVector,
	vector<ImagesDataExtended*> &shortImagesVector)
{
	list<ImagesDataExtended *> reducedImages;
	list<ImagesDataExtended *> shortImages;

	for (auto i : images)
		reducedImages.push_back(i.get());

	ConfigurationParameters fitnesObjectParams;
	if (pFitnessObjectParams)
		fitnesObjectParams = *pFitnessObjectParams;

	fitnesObjectParams.clearRetrievalMarkers();
	if (!m_fitnessObject->init(z_d, reducedImages, shortImages, &fitnesObjectParams))
		return m_fitnessObject->getErrorString();

	// Check that all keys in the fitness object parameters are actually used
	vector<string> unusedKeys;
	fitnesObjectParams.getUnretrievedKeys(unusedKeys);

	if (unusedKeys.size() > 0)
	{
		stringstream ss;

		ss << "Some parameters that were specified for the lens fitness object were not used:";
		for (auto &k : unusedKeys)
			ss << " " << k;
		return ss.str();
	}

	shortImagesVector.clear();
	reducedImagesVector.clear();

	// Transform from list to vector
	reducedImagesVector.insert(reducedImagesVector.end(), reducedImages.begin(), reducedImages.end());
	shortImagesVector.insert(shortImagesVector.end(), shortImages.begin(), shortImages.end());

	return true;
}

// TODO: re-enable this
//       Note that this is now a GenomeCalculator and that multiple instances can exist
#if 0 
void LensInversionGAFactoryCommon::onGeneticAlgorithmStart()
{
	// TODO: move this to constructor

	// To debug the mass scale search
	string fileName;
	if (getenv("GRALE_DEBUG_MASSCALESEARCH", fileName))
	{
		m_scaleSearchFileStream.open(fileName, ios_base::out);
		if (!m_scaleSearchFileStream.is_open())
			cerr << "WARNING: couldn't open 'GRALE_DEBUG_MASSCALESEARCH' debug file '" + fileName + "'" << endl;
		else
		{
			m_scaleSearchFileStream.precision(10);
			m_scaleSearchFileStream << "# " << m_massScaleSearchParams.toString() << endl;
		}
	}
}
#endif

static float LogTrans(float x) { return LN(x); }
static float ExpTrans(float x) { return EXP(x); }
static float IdentityTrans(float x) { return x; }

bool_t LensInversionGAFactoryCommon::createLens(const LensGAGenome &genome, unique_ptr<GravitationalLens> &lens) const
{
	string errStr = "unknown error";

	lens = move(createLens(genome.m_weights, genome.m_sheets, genome.m_scaleFactor, errStr));
	if (!lens.get())
		return errStr;
	return true;
}

bool_t LensInversionGAFactoryCommon::calculate(const eatk::Genome &genome, eatk::Fitness &fitness)
{
	const LensGAGenome &g = static_cast<const LensGAGenome&>(genome);
	LensGAFitness &f = static_cast<LensGAFitness&>(fitness);

	auto r = calculateFitness(g.m_weights, g.m_sheets, f.m_scaleFactor, f.m_fitnesses.data());
	if (!r)
		return "Unable to calculate fitness: " + r.getErrorString();

	// f.setCalculated(true); // Set the flag already so that toString shows something
	// cout << "Genome: " << g.toString() << " / " << f.toString() << endl;
	return true;		
}

float LensInversionGAFactoryCommon::getScalingMassSum(const vector<float> &basisFunctionWeights) const
{
	// TODO: adjust mass scale depending on sheet value?

	float massum = 0;
	int numBasisFunctions = basisFunctionWeights.size();
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
	return massum;
}

float LensInversionGAFactoryCommon::getStepsAndStepSize(pair<float,float> startStopValue, int iteration,
                                                        vector<pair<float,float>> &steps) const
{
	assert(startStopValue.first <= startStopValue.second);
	
	int numiterationsteps = (iteration == 0)?m_massScaleSearchParams.getStepsOnFirstIteration():m_massScaleSearchParams.getStepsOnSubsequentIterations();
	steps.resize(numiterationsteps);

	// Inverse transform
	float (*IT)(float x) = (useLogarithmicScaleSearch())?ExpTrans:IdentityTrans;

	float s = startStopValue.first;
	float stepsize = (startStopValue.second-startStopValue.first)/((float)(numiterationsteps-1));
	for (int j = 0 ; j < numiterationsteps ; j++, s += stepsize)
		steps[j] = { s, IT(s) };

	return stepsize;
}

pair<float,float> LensInversionGAFactoryCommon::getInitialStartStopValues(const vector<float> &basisFunctionWeights) const
{
	float massum = getScalingMassSum(basisFunctionWeights);
	float startValue = m_massScaleSearchParams.getStartFactor();
	float stopValue = m_massScaleSearchParams.getStopFactor();
	// we want a scale factor of 1 to correspond to the total mass specified by the
	// massscale used in the BackProjectMatrix

	startValue /= (massum/m_targetMass);
	stopValue /= (massum/m_targetMass);

	// Forward transform
	float (*FT)(float x) = (useLogarithmicScaleSearch())?LogTrans:IdentityTrans;
	startValue = FT(startValue);
	stopValue = FT(stopValue);

	assert(startValue < stopValue);

	return { startValue, stopValue };
}

void LensInversionGAFactoryCommon::updateStartStopValues(pair<float,float> &startStopValue, pair<float,float> startStopValue0,
										float currentBestScaleFactor, float stepsize) const
{
	assert(stepsize >= 0);
	assert(startStopValue.first <= startStopValue.second);

	startStopValue.first = currentBestScaleFactor-stepsize;
	startStopValue.second = currentBestScaleFactor+stepsize;

	// Make sure we stay within bounds
	if (startStopValue.first < startStopValue0.first)
		startStopValue.first = startStopValue0.first;
	if (startStopValue.second > startStopValue0.second)
		startStopValue.second = startStopValue0.second;

	// After these adjustments, it's possible that first became larger than second, swap them in this case
	if (startStopValue.first > startStopValue.second)
		swap(startStopValue.first, startStopValue.second);
}

bool_t LensInversionGAFactoryCommon::calculateFitness(const vector<float> &basisFunctionWeights,
													const vector<float> &sheetValues,
													float &scaleFactor,
													float *pFitnessValues
													)
{
	// TODO: re-enable something similar
	// if (m_scaleSearchFileStream.is_open()) // For debugging
	// {
	// 	m_scaleSearchStringStream.clear();
	// 	m_scaleSearchStringStream.precision(10);
	// 	m_searchedPoints.clear();
	// }

	bool_t r;

	if (!(r = initializeNewCalculation(basisFunctionWeights, sheetValues)))
		return "Can't initialize new calculation: " + r.getErrorString();

	int numiterations = getNumberOfCalculationIterations();
	if (numiterations > 0) // Scale search requested
	{
		const auto startStopValueInitial = getInitialStartStopValues(basisFunctionWeights);
		auto startStopValue = startStopValueInitial;

		float currentBestFitness = std::numeric_limits<float>::max();
		pair<float,float> currentBestStep;
		vector<pair<float,float>> &steps = m_tmpSteps;

 		for (int i = 0 ; i < numiterations ; i++)
		{
			float stepsize = getStepsAndStepSize(startStopValue, i, steps);

			for (auto step : steps) // .first is scaled step, .second is real value
			{
				float f;
				if (!(r = calculateMassScaleFitness(step.second, f)))
					return "Can't calculate mass scale fitness: " + r.getErrorString();

				if (f < currentBestFitness)
				{
					currentBestFitness = f;
					currentBestStep = step;
				}

				// TODO: re-enable something similar
				// if (m_scaleSearchFileStream.is_open()) // To debug/illustrate the scale search
				// 	m_searchedPoints.push_back({ step.second, f});
			}

			updateStartStopValues(startStopValue, startStopValueInitial, currentBestStep.first, stepsize);
		}

		scaleFactor = currentBestStep.second;

		// TODO: re-enable something similar
		// // To debug the mass scale search
		// if (m_scaleSearchFileStream.is_open())
		// {
		// 	std::sort(m_searchedPoints.begin(), m_searchedPoints.end(), [](auto x, auto y) { return x.first < y.first; });
		// 	m_scaleSearchStringStream << "{" << endl;
		// 	m_scaleSearchStringStream << " 'scaleFactor': " << scaleFactor << "," << endl;
		// 	float (*IT)(float x) = (useLogarithmicScaleSearch())?ExpTrans:IdentityTrans;
		// 	m_scaleSearchStringStream << " 'startValue0': " << IT(startStopValueInitial.first) << "," << endl;
		// 	m_scaleSearchStringStream << " 'stopValue0': " << IT(startStopValueInitial.second) << "," << endl;
		// 	m_scaleSearchStringStream << " 'startFactor': " << m_massScaleSearchParams.getStartFactor() << "," << endl;
		// 	m_scaleSearchStringStream << " 'stopFactor': " << m_massScaleSearchParams.getStopFactor() << "," << endl;
		// 	m_scaleSearchStringStream << " 'scaleSearch': [ " << endl;
		// 	for (auto x : m_searchedPoints)
		// 		m_scaleSearchStringStream << "[ " << x.first << ", " << x.second << " ]," << endl;
		// 	m_scaleSearchStringStream << " ]" << endl;
		// 	m_scaleSearchStringStream << "}," << endl;
		// }
	}
	else // numiterations <= 0, no further scale factor determination requested
	{
		scaleFactor = 1.0f;
	}

	// Do the final (possibly multi-component) fitness evaluation
	if (!(r = calculateTotalFitness(scaleFactor, pFitnessValues)))
		return "Can't calculate total fitness: " + r.getErrorString();

	return true;
}

bool_t LensInversionGAFactoryCommon::initializeNewCalculation(const std::vector<float> &basisFunctionWeights, const std::vector<float> &sheetValues)
{
	return "ERROR: initializeNewCalculation does not have an implementation";
}

bool_t LensInversionGAFactoryCommon::calculateMassScaleFitness(float scaleFactor, float &fitness)
{
	return "ERROR: calculateMassScaleFitness does not have a default implementation";
}

bool_t LensInversionGAFactoryCommon::calculateTotalFitness(float scaleFactor, float *pFitnessValues)
{
	return "ERROR: calculateTotalFitness does not have a default implementation";
}

// TODO: make this available again
/*
void LensInversionGAFactoryCommon::onSortedPopulation(const std::vector<mogal::GenomeWrapper> &population)
{
	if (m_scaleSearchFileStream.is_open())
	{
		m_scaleSearchFileStream << "[" << endl;
		m_scaleSearchFileStream << m_scaleSearchStringStream.str() << endl;
		m_scaleSearchFileStream << "]," << endl;
	}
}
*/

} // end namespace
