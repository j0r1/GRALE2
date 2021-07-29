#include "lensinversiongafactorymultiplanegpu.h"
#include "multiplanecontainer.h"
#include "compositelens.h"
#include "masssheetlens.h"
#include "plummerlens.h"
#include "openclcalculator.h"
#include "utils.h"
#include <assert.h>

using namespace std;
using namespace errut;

namespace grale
{


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

	// Currently for the mass sheet basis functions, may expand on this later
	m_unscaledBasisLenses.clear();
	if (pParams->useMassSheetBasisFunctions())
	{
		const auto &redshifts = pParams->getLensRedshifts();
		const Cosmology &cosmology = pParams->getCosmology();

		for (double z : redshifts)
		{
			double Dd = cosmology.getAngularDiameterDistance(z);
			MassSheetLensParams params(Dd, 1.1, 1); // Same as in CPU case
			shared_ptr<MassSheetLens> pSheetLens = make_shared<MassSheetLens>();
			if (pSheetLens->init(Dd, &params))
				return "Error creating one of the mass sheet basis functions: " + pSheetLens->getErrorString();
			m_unscaledBasisLenses.push_back(pSheetLens);
		}
	}

	// Call setCommonParameters
	if (!(r = setCommonParameters(m_unscaledBasisLenses.size(),
						pParams->getAllowNegativeWeights(), m_basisFunctionMasses,
						pParams->getMassEstimate(), pParams->getMassEstimate(),
						pParams->getMassScaleSearchParameters())))
		return r;

	// We'll init the GPU part here, so we can use it immediately later. This is
	// a single instance for this process, to be used by the available threads.
	// TODO: what if deviceIndex is different in different threads?
	if (!(r = OpenCLCalculator::initInstance(pParams->getDeviceIndex(),
	                                         m_reducedImages, m_shortImages,
											 m_lensRedshifts,
											 pParams->getCosmology(),
											 pParams->getBasisLenses(),
											 m_unscaledBasisLenses
											 )))
		return r;

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

	for (size_t i = 1 ;  i < redshifts.size() ; i++)
	{
		if (redshifts[i] <= redshifts[i-1])
			return "Lens redshifts must be in strictly increasing order (" + to_string(redshifts[i]) + " <= " + to_string(redshifts[i-1]) + ")";
	}

	// Store the doubles in a float vector
	m_lensRedshifts.clear();
	m_lensRedshifts.insert(m_lensRedshifts.end(), redshifts.begin(), redshifts.end());
	
	m_basisFunctionMasses.clear();

	for (auto &plane : basisLenses)
	{
		if (plane.size() == 0)
			return "A lens plane without basis functions is present";

		for (auto &bl : plane)
		{
			if (bl->m_relevantLensingMass < 0)
				return "A basis lens was found to have a negative strong lensing mass";

			m_basisFunctionMasses.push_back(bl->m_relevantLensingMass);
		}
	}

	// We'll need to re-analyze these later as well, to build up the OpenCL kernel and
	// Upload the parameters

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

	assert((useSheets && sheetValues.size() == redshifts.size()) || (!useSheets && sheetValues.size() == 0));

	double scale = scaleFactor;
	int basisFunctionWeightIdx = 0;

	MultiPlaneContainerParams containerParams;
	
	for (int planeIdx = 0 ; planeIdx < basisFunctions.size() ; planeIdx++)
	{
		double z = redshifts[planeIdx];
		double Dd = cosm.getAngularDiameterDistance(z);
		auto planeBasisFunctions = basisFunctions[planeIdx];
		const GravitationalLens *pUnscaledBf = (useSheets)?m_unscaledBasisLenses[planeIdx].get():nullptr;
		
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
			assert(pUnscaledBf);
			if (!planeLensParams.addLens(1.0, Vector2Dd(0, 0), 0, *pUnscaledBf))
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

bool_t LensInversionGAFactoryMultiPlaneGPU::onNewCalculationStart(size_t genomesInThisThread, size_t genomesInAllThreads)
{
	OpenCLCalculator &oclCalc = OpenCLCalculator::instance();
	oclCalc.setGenomesToCalculate(genomesInAllThreads);
	return true;
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

	if (m_currentParams->getMassScaleSearchParameters().getNumberOfIterations() == 0) // No search requested, just use weights
		m_calcStates[index].reset({ numeric_limits<float>::quiet_NaN(), numeric_limits<float>::quiet_NaN() }, true);
	else
		m_calcStates[index].reset(getInitialStartStopValues(g.m_weights), false);
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

	if (!state.m_calculationScheduled)
	{
		state.m_calculationScheduled = true;
		if (state.m_isFinalCalculation)
			r = oclCalc.scheduleUploadAndCalculation(index, { { numeric_limits<float>::quiet_NaN(), state.m_bestScaleFactor } }, false);
		else
		{
			state.m_stepSize = getStepsAndStepSize(state.m_startStopValues, state.m_nextIteration, state.m_steps);
			r = oclCalc.scheduleUploadAndCalculation(index, state.m_steps, true);
		}
		if (!r)
			return "Unable to schedule upload of steps to be calculated: " + r.getErrorString();
	}
	else
	{
		if (oclCalc.isCalculationDone())
		{
			state.m_calculationScheduled = false;

			if (state.m_isFinalCalculation)
			{
				LensGAFitness &f = static_cast<LensGAFitness&>(fitness);

				// Calculate full fitness

				// TODO: just dummy
				for (auto &x : f.m_fitnesses)
					x = (float)index;

				// Tell opencl module we're done for this genome; if all genomes were calculated, we can
				// prepare for the next iteration (new genomes)
				return "TODO: tell OpenCLCalculator this genome is done";
			}
			else
			{
				// TODO: Calculate fitness values for all the steps for which the backprojection was calculated
				
				// TODO: just a test
				state.m_bestScaleFactor = state.m_steps[state.m_steps.size()/2].second;

				state.m_nextIteration++;
				if (state.m_nextIteration >= m_currentParams->getMassScaleSearchParameters().getNumberOfIterations())
					state.m_isFinalCalculation = true; // we'll calculate the best scale factor results in full next
				else
				{
					//cout << "Updating start/stop values for genome " << index << ", iteration is now " << state.m_nextIteration << endl;
					updateStartStopValues(state.m_startStopValues, state.m_initialStartStopValues, state.m_bestScaleFactor, state.m_stepSize);
				}
			}
		}
		else
		{
			// Nothing to do but wait till GPU is ready
		}
	}
	
	return true;
}

} // end namespace
