#include "lensinversiongafactorymultiplanegpu.h"
#include "multiplanecontainer.h"
#include "compositelens.h"
#include "masssheetlens.h"
#include "plummerlens.h"
#include "openclcalculator.h"
#include "utils.h"
#include <thread>
#include <fstream>
#include <assert.h>
#include <sys/time.h>

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
	if (m_currentParams.get()) // Initialization worked, we can release GPU part now
		OpenCLCalculator::releaseInstance((uint64_t)this);
	m_perNodeCounter.reset();
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
			if (!pSheetLens->init(Dd, &params))
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

	int devIdx = pParams->getDeviceIndex();
	if (devIdx < 0) // Means rotate over the available devices
	{
		string fileName = "/dev/shm/grale_mpopencl_nextdevice.dat";
		getenv("GRALE_OPENCL_AUTODEVICEFILE", fileName); // Doesn't change file name if envvar not set

		m_perNodeCounter = make_unique<PerNodeCounter>(fileName);

		int idx = m_perNodeCounter->getCount();
		if (idx < 0)
		{
			m_perNodeCounter.reset();
			return "Couldn't read per-node device index from file '" + fileName + "': " + m_perNodeCounter->getErrorString();
		}

		OpenCLLibrary clLib;
		if (!clLib.loadLibrary(OpenCLLibrary::getLibraryName()))
			return "Can't open OpenCL library: " + clLib.getErrorString();

		int numDevices = clLib.getDeviceCount();
		if (numDevices < 0)
			return "Error getting device count: " + clLib.getErrorString();
		if (numDevices == 0)
			return "Unexpectedly got zero GPU devices";

		devIdx = idx%numDevices;

		auto GetTimeStamp = []() {
			struct timeval tv;
			gettimeofday(&tv,NULL);
			return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
		};

		cerr << "DEBUG (" << GetTimeStamp() << "): Automatically using new device index " << devIdx << endl;
	}

	// We'll init the GPU part here, so we can use it immediately later. This is
	// a single instance for this process, to be used by the available threads.
	// TODO: what if deviceIndex is different in different threads?
	if (!(r = OpenCLCalculator::initInstance(devIdx,
	                                         m_reducedImages, m_shortImages,
											 m_lensRedshifts,
											 pParams->getCosmology(),
											 pParams->getBasisLenses(),
											 m_unscaledBasisLenses,
											 (uint64_t)this // Use this pointer as an identifier
											 )))
		return r;

	const double angScale = OpenCLCalculator::instance().getAngularScale();
	m_bpAll = make_shared<OclCalculatedBackProjector>();
	if (!(r = m_bpAll->init(m_reducedImages, angScale)))
		return "Error initializing backprojection framework: " + r.getErrorString();

	if (m_shortImages.empty()) // use full images
		m_bpShort = m_bpAll;
	else
	{
		m_bpShort = make_shared<OclCalculatedBackProjector>();
		if (!(r = m_bpShort->init(m_shortImages, angScale)))
			return "Error initializing backprojection framework for short images: " + r.getErrorString();
	}

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
			assert(planeIdx < sheetValues.size());
			if (!planeLensParams.addLens(sheetValues[planeIdx], Vector2Dd(0, 0), 0, *pUnscaledBf))
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
	m_calcStates.clear();
	return true;
}

bool_t LensInversionGAFactoryMultiPlaneGPU::startNewCalculation(const eatk::Genome &genome)
{
	const LensGAGenome &g = static_cast<const LensGAGenome&>(genome);
	OpenCLCalculator &oclCalc = OpenCLCalculator::instance();
	bool_t r;
	
	if (!(r = oclCalc.startNewBackprojection(g)))
		return "Error starting new calculation for genome: " + r.getErrorString();

	State &s = m_calcStates[&g];
	if (m_currentParams->getMassScaleSearchParameters().getNumberOfIterations() == 0) // No search requested, just use weights
		s.reset({ numeric_limits<float>::quiet_NaN(), numeric_limits<float>::quiet_NaN() }, true);
	else
		s.reset(getInitialStartStopValues(g.m_weights), false);
	return true;
}

bool_t LensInversionGAFactoryMultiPlaneGPU::pollCalculate(const eatk::Genome &genome, eatk::Fitness &fitness)
{
	assert(!fitness.isCalculated());

	const LensGAGenome &g = static_cast<const LensGAGenome&>(genome);
	LensGAFitness &f = static_cast<LensGAFitness&>(fitness);
	OpenCLCalculator &oclCalc = OpenCLCalculator::instance();
	bool_t r;

	State &state = m_calcStates[&g];

	if (state.m_calculationIdentifier < 0) // Not scheduled yet
	{
		// Try to schedule calculation
		if (state.m_isFinalCalculation)
		{
			state.m_stepSize = numeric_limits<float>::quiet_NaN();
			state.m_steps = { { numeric_limits<float>::quiet_NaN(), state.m_bestScaleFactor } };
		}
		else
		{
			state.m_stepSize = getStepsAndStepSize(state.m_startStopValues, state.m_nextIteration, state.m_steps);
		}

		if (!(r = oclCalc.scheduleUploadAndCalculation(g, state.m_calculationIdentifier, state.m_steps, !state.m_isFinalCalculation)))
			return "Unable to schedule upload of steps to be calculated: " + r.getErrorString();

		if (state.m_calculationIdentifier >= 0)
		{
			//cout << "Successfully scheduled new calculation for genome " << (void*)&g << ", calc identifier is " << state.m_calculationIdentifier << ", thread is " << std::this_thread::get_id() << endl;
			//cout.flush();
		}
	}
	else
	{
		size_t numPoints = 0, numGenomes = 0, numSteps = 0, genomeIndex = 0;
		bool calcDone = false;
		const float *pAllBetas = nullptr;

		if (!(r = oclCalc.isCalculationDone(g, state.m_calculationIdentifier, &genomeIndex, &pAllBetas, &numGenomes, &numSteps, &numPoints, &calcDone)))
			return "Error checking calculation done: " + r.getErrorString();

		if (calcDone)
		{
			// ofstream dbg("betas.dat");
			// dbg.write((char*)pAllBetas, numSteps*(numPoints*2)*numGenomes*sizeof(float));

			//cout << "calcDone, numPoints = " << numPoints << " numGenomes = " << numGenomes << " numSteps = " << numSteps << " genomeIndex = " << genomeIndex << endl;
			//cout << "          state.m_steps.size() = " << state.m_steps.size() << endl;
			assert(numSteps == state.m_steps.size());
			assert(genomeIndex < numGenomes);

			// cout << "Got results for genome " << (void*)&g << endl;
			// cout << "   m_initialStartStopValues = " << state.m_initialStartStopValues.first << "," << state.m_initialStartStopValues.second << endl;
			// cout << "   m_startStopValues = " << state.m_startStopValues.first << "," << state.m_startStopValues.second << endl;
			// cout << "   m_bestScaleFactor = " << state.m_bestScaleFactor << endl;
			// cout << "   m_isFinalCalculation = " << (int)state.m_isFinalCalculation << endl;
			// cout << "   m_stepSize = " << state.m_stepSize << endl;
			// cout << "   m_nextIteration = " << state.m_nextIteration << endl;
			// cout << "   m_steps = " << endl;
			// for (auto s : state.m_steps)
			// 	cout << "      " << s.first << "," << s.second << endl;

			const float *pBetasForGenome = pAllBetas + numSteps*(numPoints*2)*genomeIndex;
			LensFitnessObject &fitnessFunction = getFitnessObject();

			if (state.m_isFinalCalculation)
			{
				assert(numSteps == 1);
				LensGAFitness &f = static_cast<LensGAFitness&>(fitness);

				m_bpAll->setBetaBuffer(pBetasForGenome, numPoints*2);

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
				dumpBp(*m_bpAll);
				*/

				fitnessFunction.calculateOverallFitness(*m_bpAll, f.m_fitnesses.data());
				f.m_scaleFactor = state.m_bestScaleFactor;
				f.setCalculated(true);

				// cout << "Genome: " << g.toString() << " / " << f.toString() << endl;
			}
			else
			{
				float bestFitness = numeric_limits<float>::infinity();
				int bestStep = -1;
				for (size_t s = 0 ; s < numSteps ; s++)
				{
					const float *pBetasForStep = pBetasForGenome + (numPoints*2)*s;
					m_bpShort->setBetaBuffer(pBetasForStep, numPoints*2);

					float fitness = numeric_limits<float>::quiet_NaN();
					if (!fitnessFunction.calculateMassScaleFitness(*m_bpShort, fitness))
						return "Error calculating mass scale (short) fitness: " + fitnessFunction.getErrorString();
					// cout << "Fitness for step " << s << " " << state.m_steps[s].first << " <=> " << state.m_steps[s].second << " = " << fitness << endl;

					if (fitness < bestFitness)
					{
						bestFitness = fitness;
						bestStep = s;
					}
				}

				state.m_bestScaleFactor = state.m_steps[bestStep].second;

				state.m_nextIteration++;
				if (state.m_nextIteration >= m_currentParams->getMassScaleSearchParameters().getNumberOfIterations())
					state.m_isFinalCalculation = true; // we'll calculate the best scale factor results in full next
				else
				{
					updateStartStopValues(state.m_startStopValues, state.m_initialStartStopValues, state.m_steps[bestStep].first, state.m_stepSize);
				}
			}

			// Tell the OpenCL part that we're done with these results; if that's the case
			// for all genomes, we don't need this anymore
			if (!(r = oclCalc.setCalculationProcessed(g, state.m_calculationIdentifier)))
				return "Error telling OpenCL part that our calculation is processed: " + r.getErrorString();
			state.m_calculationIdentifier = -1;
		}
		else
		{
			// Nothing to do but wait till GPU is ready
		}
	}
	
	return true;
}

} // end namespace
