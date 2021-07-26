#include "lensinversiongafactorymultiplanegpu.h"
#include "multiplanecontainer.h"
#include "compositelens.h"
#include "masssheetlens.h"
#include "plummerlens.h"
#include "utils.h"
#include <assert.h>
#include <mutex>
#include <map>

using namespace std;
using namespace errut;

namespace grale
{

class OpenCLCalculator
{
public:
    static bool_t initInstance(int devIdx)
	{
		lock_guard<mutex> guard(s_instanceMutex);

		if (s_pInstance.get()) // Already initialized
			return true;

		if (s_initTried)
			return "GPU initialization failed before";
		s_initTried = true;

		unique_ptr<OpenCLCalculator> oclCalc = make_unique<OpenCLCalculator>();
		bool_t r = oclCalc->init(devIdx);
		if (!r)
			return "Can't init GPU: " + r.getErrorString();

		s_pInstance = move(oclCalc);
		return true;
	}

	static OpenCLCalculator &instance()
	{
		OpenCLCalculator *pInst = s_pInstance.get();
		assert(pInst);
		return *pInst;
	}

	bool_t startNewBackprojection(const LensGAGenome &g, size_t &genomeIndex)
	{
		lock_guard<mutex> guard(m_mutex);
		
		if (m_hasCalculated)
		{
			m_hasCalculated = false;
			m_states.clear();
			m_nextGenomeIndex = 0;
			m_uploadInfoCount = 0;
			m_devStatus = Ready;
		}

		size_t idx = m_nextGenomeIndex;
		m_states[&g] = make_unique<State>(m_nextGenomeIndex);
		m_nextGenomeIndex++;

		size_t numBfWeights = g.m_weights.size();
		size_t numSheetWeights = g.m_sheets.size();
		size_t startBfIdx = idx * numBfWeights;
		size_t endBfIdx = startBfIdx + numBfWeights;
		size_t startSheetIdx = idx * numSheetWeights;
		size_t endSheetIdx = startSheetIdx + numSheetWeights;

		if (endBfIdx > m_allBasisFunctionWeights.size())
			m_allBasisFunctionWeights.resize(endBfIdx);
		if (endSheetIdx > m_allSheetWeights.size())
			m_allSheetWeights.resize(endSheetIdx);
		
		memcpy(m_allBasisFunctionWeights.data() + startBfIdx, g.m_weights.data(), sizeof(float)*numBfWeights);
		if (numSheetWeights > 0)
			memcpy(m_allSheetWeights.data() + startSheetIdx, g.m_sheets.data(), sizeof(float)*numSheetWeights);
		
		genomeIndex = idx;
		return true;
	}

	enum DeviceStatus { Ready, UploadingWeights, Unknown };

	bool_t scheduleUploadAndCalculation(size_t genomeIndex, const vector<pair<float,float>> &scaleFactors)
	{
		lock_guard<mutex> guard(m_mutex);

		if (m_devStatus == Ready) // First time this is called after setting up the weights
		{
			m_devStatus = UploadingWeights;
			// TODO: start upload of weights to GPU
		}

		// We're assuming that all genome fitness calculations will run in sync, so we can use scaleFactors.size()
		// to calculate an offset
		size_t startIdx = genomeIndex * scaleFactors.size();
		size_t endIdx = startIdx + scaleFactors.size();

		if (endIdx > m_allFactors.size())
			m_allFactors.resize(endIdx);

		for (size_t i = 0, j = startIdx ; i < scaleFactors.size() ; i++, j++)
		{
		 	assert(j < endIdx && j < m_allFactors.size());
			m_allFactors[j] = scaleFactors[i].second; // .second is the real value to be used
		}

		// TODO: make sure we've got an array somewhere to store the backprojection results
		// TODO: this is different for the scale search points and the final points

		m_uploadInfoCount++;
		if (m_uploadInfoCount == m_states.size()) // Information for all genomes is present!
		{
			m_hasCalculated = true; // make sure we clear some state on next iteration (don't do this too soon, other threads may otherwise still be in the startNewBackprojection stage)

			// TODO: Upload the info, and start the OpenCL calculation
			cout << "Ready to upload and backproject" << endl;
		}
		return true;
	}

	bool_t getGenomeIndex(const LensGAGenome &g, size_t &genomeIndex) const
	{
		// Should we use a mutex here?
		auto it = m_states.find(&g);
		if (it == m_states.end())
			return "Can't find state for genome";
		genomeIndex = it->second->m_genomeIndex;
		return true;
	}

	OpenCLCalculator()
	{
	}

	~OpenCLCalculator()
	{
	}
private:
	bool_t init(int devIdx)
	{
		m_hasCalculated = true; // Causes the internal state to be reset
		return true;
	}

	class State
	{
	public:
		State(size_t genomeIndex) : m_genomeIndex(genomeIndex) { }
		const size_t m_genomeIndex;
	};

	map<const LensGAGenome *, unique_ptr<State>> m_states;
	size_t m_nextGenomeIndex;
	bool m_hasCalculated;
	mutex m_mutex;
	vector<float> m_allBasisFunctionWeights;
	vector<float> m_allSheetWeights;
	vector<float> m_allFactors;
	size_t m_uploadInfoCount;
	DeviceStatus m_devStatus;

	static unique_ptr<OpenCLCalculator> s_pInstance;
	static mutex s_instanceMutex;
	static bool s_initTried;
};

unique_ptr<OpenCLCalculator> OpenCLCalculator::s_pInstance;
mutex OpenCLCalculator::s_instanceMutex;
bool OpenCLCalculator::s_initTried = false;




LensInversionGAFactoryMultiPlaneGPU::LensInversionGAFactoryMultiPlaneGPU(unique_ptr<LensFitnessObject> fitObj)
	: LensInversionGAFactoryCommon(move(fitObj))
{

}

LensInversionGAFactoryMultiPlaneGPU::~LensInversionGAFactoryMultiPlaneGPU()
{

}

bool_t LensInversionGAFactoryMultiPlaneGPU::init(const LensInversionParametersBase &p)
{
	auto pParams = dynamic_cast<const LensInversionParametersMultiPlaneGPU *>(&p);
	if (!pParams)
		return "Specified parameters are not of the correct type";

	bool_t r;

	if (!(r = analyzeLensBasisFunctions(pParams->getLensRedshifts(), pParams->getBasisLenses())))
		return r;

	m_images.clear();
	if (!(r = analyzeSourceImages(pParams->getSourceImages(), pParams->getCosmology(), m_images)))
		return r;

	// These are the ones to actually use
	m_reducedImages.clear();
	m_shortImages.clear();
	if (!(r = initializeLensFitnessObject(numeric_limits<double>::quiet_NaN(), m_images, pParams->getFitnessObjectParameters(),
									 m_reducedImages, m_shortImages)))
		return r;

	// sheet densities, sheet multipliers
	m_sheetDensities.clear();
	m_sheetMultipliers.clear();
	if (pParams->useMassSheetBasisFunctions())
	{
		const auto &redshifts = pParams->getLensRedshifts();
		const Cosmology &cosmology = pParams->getCosmology();
		m_sheetDensities.resize(redshifts.size(), 0); // Just allocate some space

		for (double z : redshifts)
		{
			double Dd = cosmology.getAngularDiameterDistance(z);
			MassSheetLensParams params(Dd, 1.1, 1); // Same as in CPU case
			
			// A value of 1 in the genome will then correspond to this density
			m_sheetMultipliers.push_back((float)params.getDensity());
		}
	}

	// Call setCommonParameters
	if (!(r = setCommonParameters(m_sheetMultipliers.size(),
						pParams->getAllowNegativeWeights(), m_basisFunctionMasses,
						pParams->getMassEstimate(), pParams->getMassEstimate(),
						pParams->getMassScaleSearchParameters())))
		return r;

	// We'll init the GPU part here, so we can use it immediately later. This is
	// a single instance for this process, to be used by the available threads.
	// TODO: what if deviceIndex is different in different threads?
	if (!(r = OpenCLCalculator::initInstance(pParams->getDeviceIndex())))
		return r;

	// TODO: we need to setup the backprojection for this particular layout
	// TODO: we don't really want to do this multiple times (for each thread), doesn't
	//       really matter that much, but doesn't feel very clean

	m_currentParams = make_unique<LensInversionParametersMultiPlaneGPU>(*pParams);
	return true;
}

bool_t LensInversionGAFactoryMultiPlaneGPU::analyzeLensBasisFunctions(const vector<double> redshifts,
						  const vector<vector<shared_ptr<LensInversionBasisLensInfo>>> &basisLenses)
{
	if (redshifts.size() == 0)
		return "No lens planes present";

	if (redshifts.size() != basisLenses.size())
		return "Incompatible number of redshifts (" + to_string(redshifts.size()) + ") and lens planes (" + to_string(basisLenses.size()) + ")";

	// Store the doubles in a float vector
	m_lensRedshifts.clear();
	m_lensRedshifts.insert(m_lensRedshifts.end(), redshifts.begin(), redshifts.end());
	
	m_basisFunctionMasses.clear();
	m_basisLenses.clear();

	m_scaledPlaneWeights.clear();
	m_basePlaneWeights.clear();
	for (auto &plane : basisLenses)
	{
		if (plane.size() == 0)
			return "A lens plane without basis functions is present";

		// Just create the right sizes for these vectors
		m_scaledPlaneWeights.push_back(vector<float>(plane.size(), 0));
		m_basePlaneWeights.push_back(vector<float>(plane.size(), 0));

		m_basisLenses.push_back(vector<PlummerLensInfo>());
		auto &plummers = m_basisLenses.back();

		for (auto &bl : plane)
		{
			if (bl->m_relevantLensingMass < 0)
				return "A basis lens was found to have a negative strong lensing mass";

			const PlummerLens *pPlummerLens = dynamic_cast<const PlummerLens *>(bl->m_pLens.get());
			if (!pPlummerLens)
				return "Not all basis functions appear to be Plummer lens models";

			const PlummerLensParams *pParams = dynamic_cast<const PlummerLensParams *>(pPlummerLens->getLensParameters());
			if (!pParams)
				return "Unexpected: couldn't get Plummer lens parameters";

			plummers.push_back(PlummerLensInfo { pParams->getLensMass(), pParams->getAngularWidth(), bl->m_center });
			m_basisFunctionMasses.push_back(bl->m_relevantLensingMass);
		}
	}

	return true;
}

bool_t LensInversionGAFactoryMultiPlaneGPU::analyzeSourceImages(const vector<shared_ptr<ImagesDataExtended>> &sourceImages,
						 const Cosmology &cosmology,
						 vector<shared_ptr<ImagesDataExtended>> &imagesVector)
{
	if (sourceImages.size() == 0)
		return "No images present";

	imagesVector.clear();
	
	int count = 0;
	for (auto img : sourceImages)
	{
		count++;

		double z;
		if (!img->hasExtraParameter("z") || !img->getExtraParameter("z", z))
			return "Images data set " + to_string(count) + " does not contain the redshift 'z' parameter";

		if (img->getDds() != 0 || img->getDs() != 0)
			return "Images data set " + to_string(count) + " does not have Dds and Ds set to zero, required for the multi-plane inversion";

		// Calculate Ds from redshift, and set it in a copy
		double Ds = cosmology.getAngularDiameterDistance(z);
		shared_ptr<ImagesDataExtended> newImg = make_shared<ImagesDataExtended>(*img);
		newImg->setDs(Ds);

		imagesVector.push_back(newImg);
	}
	return true;
}

unique_ptr<GravitationalLens> LensInversionGAFactoryMultiPlaneGPU::createLens(const std::vector<float> &basisFunctionWeights,
	                                      const std::vector<float> &sheetValues,
										  float scaleFactor,
										  std::string &errStr) const
{
	bool useSheets = m_currentParams->useMassSheetBasisFunctions();
	auto &cosm = m_currentParams->getCosmology();
	auto &basisFunctions = m_currentParams->getBasisLenses();
	auto &redshifts = m_currentParams->getLensRedshifts();
	assert(redshifts.size() == basisFunctions.size());

	vector<float> sheetDensities(sheetValues.size());
	assert((useSheets && sheetValues.size() == redshifts.size()) || (!useSheets && sheetValues.size() == 0));

	convertGenomeSheetValuesToDensities(sheetValues, sheetDensities);
	double scale = scaleFactor;
	int basisFunctionWeightIdx = 0;

	MultiPlaneContainerParams containerParams;
	
	for (int planeIdx = 0 ; planeIdx < basisFunctions.size() ; planeIdx++)
	{
		double z = redshifts[planeIdx];
		double Dd = cosm.getAngularDiameterDistance(z);
		auto planeBasisFunctions = basisFunctions[planeIdx];
		
		CompositeLensParams planeLensParams;

		for (auto &bf : planeBasisFunctions)
		{
			assert(basisFunctionWeightIdx < basisFunctionWeights.size());

			if (!planeLensParams.addLens(basisFunctionWeights[basisFunctionWeightIdx++]*scale,
									bf->m_center, 0, *(bf->m_pLens.get())))
			{
				errStr = "Unable to add basis function to composite lens: " + planeLensParams.getErrorString();
				return nullptr;
			}
		}

		if (useSheets)
		{
			assert(planeIdx < sheetDensities.size());
			MassSheetLensParams sheetParams(sheetDensities[planeIdx]);
			MassSheetLens sheetLens;

			if (!sheetLens.init(Dd, &sheetParams))
			{
				errStr = "Could not initialize a mass sheet lens: " + sheetLens.getErrorString();
				return nullptr;
			}
			if (!planeLensParams.addLens(1.0, Vector2Dd(0, 0), 0, sheetLens))
			{
				errStr = "Unable to add sheet lens to composite lens: " + planeLensParams.getErrorString();
				return nullptr;
			}
		}

		auto compLens = make_shared<CompositeLens>();
		if (!compLens->init(Dd, &planeLensParams))
		{
			errStr = "Unable to create a composite lens for a lens plane: " + compLens->getErrorString();
			return nullptr;
		}
		containerParams.add(compLens, z);
	}

	assert(basisFunctionWeightIdx == basisFunctionWeights.size());
	
	unique_ptr<MultiPlaneContainer> containerLens = make_unique<MultiPlaneContainer>();
	if (!containerLens->init(0, &containerParams)) // Dd must be set to 0!
	{
		errStr = "Unable to initialize the multi-plane container: " + containerLens->getErrorString();
		return nullptr;
	}

	return containerLens;
}

void LensInversionGAFactoryMultiPlaneGPU::convertGenomeSheetValuesToDensities(const vector<float> &sheetValues,
																			  vector<float> &sheetDensities) const
{
	assert(sheetValues.size() == m_sheetMultipliers.size());
	assert(sheetValues.size() == sheetDensities.size());

	for (size_t i = 0 ; i < sheetValues.size() ; i++)
		sheetDensities[i] = sheetValues[i]*m_sheetMultipliers[i];
}


bool_t LensInversionGAFactoryMultiPlaneGPU::startNewCalculation(const eatk::Genome &genome)
{
	const LensGAGenome &g = static_cast<const LensGAGenome&>(genome);
	OpenCLCalculator &oclCalc = OpenCLCalculator::instance();
	bool_t r;
	size_t index;
	
	if (!(r = oclCalc.startNewBackprojection(g, index)))
		return "Error starting new calculation for genome: " + r.getErrorString();

	// Make sure that the state vector is large enough
	lock_guard<mutex> stateVectorMutex(m_calcStateMutex);
	if (index >= m_calcStates.size())
		m_calcStates.resize(index+1);

	m_calcStates[index].reset(getInitialStartStopValues(g.m_weights));
	return true;
}

bool_t LensInversionGAFactoryMultiPlaneGPU::pollCalculate(const eatk::Genome &genome, eatk::Fitness &fitness)
{
	assert(!fitness.isCalculated());

	const LensGAGenome &g = static_cast<const LensGAGenome&>(genome);
	LensGAFitness &f = static_cast<LensGAFitness&>(fitness);
	OpenCLCalculator &oclCalc = OpenCLCalculator::instance();
	bool_t r;
	size_t index;

	if (!(r = oclCalc.getGenomeIndex(g, index)))
		return "Can't get genome index: " + r.getErrorString();

	assert(index < m_calcStates.size());
	State &state = m_calcStates[index];

	state.m_stepSize = getStepsAndStepSize(state.m_startStopValues, state.m_nextIteration, state.m_steps);

	if (!(r = oclCalc.scheduleUploadAndCalculation(index, state.m_steps)))
		return "Unable to schedule upload of steps to be calculated: " + r.getErrorString();

	return "TODO: implement pollCalculate";
}

} // end namespace
