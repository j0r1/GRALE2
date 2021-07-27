#include "lensinversiongafactorymultiplanegpu.h"
#include "multiplanecontainer.h"
#include "compositelens.h"
#include "masssheetlens.h"
#include "plummerlens.h"
#include "openclkernel.h"
#include "utils.h"
#include <assert.h>
#include <mutex>
#include <map>

using namespace std;
using namespace errut;

namespace grale
{

class OpenCLCalculator : public OpenCLKernel
{
public:
    static bool_t initInstance(int devIdx, const std::vector<ImagesDataExtended *> &allImages,
	                                       const std::vector<ImagesDataExtended *> &shortImages,
										   const std::vector<float> &zds,
										   const Cosmology &cosm)
	{
		lock_guard<mutex> guard(s_instanceMutex);

		if (s_pInstance.get()) // Already initialized
			return true;

		if (s_initTried)
			return "GPU initialization failed before";
		s_initTried = true;

		unique_ptr<OpenCLCalculator> oclCalc = make_unique<OpenCLCalculator>();
		bool_t r = oclCalc->initAll(devIdx, allImages, shortImages, zds, cosm);
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
		// TODO: cleanup!
	}
private:
	bool_t initAll(int devIdx, const std::vector<ImagesDataExtended *> &allImages,
	                        const std::vector<ImagesDataExtended *> &shortImages,
							const std::vector<float> &zds,
							const Cosmology &cosm)
	{
		bool_t r;
		if (!(r = initGPU(devIdx)))
			return "Error initializing GPU: " + r.getErrorString();

		m_angularScale = ANGLE_ARCSEC; // TODO: calculate something better?
		m_potentialScale = ANGLE_ARCSEC*ANGLE_ARCSEC; // TODO

		if (!(r = setupImages(allImages, shortImages)))
			return "Can't setup images for GPU backprojection: " + r.getErrorString();

		if (!(r = setupMultiPlaneDistanceMatrix(cosm, zds)))
			return "Can't setup multi-plane distance matrix: " + r.getErrorString();

		if (!(r = setupAngularDiameterDistances(cosm, zds, allImages, m_pDevDpointsAll, m_pDevUsedPlanesAll)))
			return "Can't setup angular diameter distances for all points: " + r.getErrorString();
		if (!m_shortImagesAreAllImages)
		{
			if (!(r = setupAngularDiameterDistances(cosm, zds, shortImages, m_pDevDpointsShort, m_pDevUsedPlanesShort)))
				return "Can't setup angular diameter distances for short points: " + r.getErrorString();
		}
		else
		{
			m_pDevUsedPlanesShort = m_pDevUsedPlanesAll;
			m_pDevUsedPlanesShort = m_pDevUsedPlanesAll;
		}
		


		m_hasCalculated = true; // Causes the internal state to be reset

		return true;
	}

	bool_t initGPU(int devIdx)
	{
		string library = getLibraryName();
		if (!loadLibrary(library))
			return "Can't load OpenCL library: " + getErrorString();
		if (!OpenCLKernel::init(devIdx))
			return "Can't init specified GPU device: " + getErrorString();
		return true;
	}

	bool_t setupMultiPlaneDistanceMatrix(const Cosmology &cosm, const vector<float> &zds)
	{
		size_t numPlanes = zds.size();
		size_t cols = numPlanes;
		size_t rows = numPlanes+1;
		
		vector<float> Dij(rows*cols, numeric_limits<float>::quiet_NaN());
		for (size_t j = 0 ; j < numPlanes ; j++)
			Dij[0 + j] = (float)(cosm.getAngularDiameterDistance(zds[j])/DIST_MPC);

		for (size_t i = 1 ; i <= numPlanes ; i++)
		{
			for (size_t j = i ; j < numPlanes ; j++)
				Dij[i*cols + j] = (float)(cosm.getAngularDiameterDistance(zds[i-1], zds[j])/DIST_MPC);
		}
		
		cl_int err;
		m_pDevDmatrix = clCreateBuffer(getContext(), CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float)*Dij.size(), Dij.data(), &err);
		if (err != CL_SUCCESS)
			return "Error uploading lens plane distance matrix to GPU";

		cout << "Distance matrix uploaded: " << endl;
		for (size_t i = 0 ; i <= numPlanes ; i++)
		{
			for (size_t j = 0 ; j < numPlanes ; j++)
				cout << "\t" << Dij[cols*i+j];
			cout << endl;
		}
		return true;
	}

	bool_t setupAngularDiameterDistances(const Cosmology &cosm, const vector<float> &zds,
	                                     const vector<ImagesDataExtended*> &images,
										 cl_mem &pDevDpoints, cl_mem &pDevUsedPlanes)
	{
		// Build the image points redshift vector
		// TODO: For now, we're associating a zs to every single point, not just to the entire source
		vector<float> zss;
		for (auto img : images)
		{
			double zs = numeric_limits<double>::quiet_NaN();
			img->getExtraParameter("z", zs);

			int numImages = img->getNumberOfImages();
			for (int i = 0 ; i < numImages ; i++)
			{
				int numPoints = img->getNumberOfImagePoints(i);
				for (int p = 0 ; p < numPoints ; p++)
					zss.push_back((float)zs);
			}
		}

		size_t numPlanes = zds.size();
		size_t numPoints = zss.size();

		size_t numCols = numPlanes+1;
		size_t numRows = numPoints;
		vector<float> Dsources(numRows*numCols, numeric_limits<float>::quiet_NaN());
		vector<cl_int> usedPlanes(numPoints);

		for (size_t p = 0 ; p < numPoints ; p++)
		{
			float zs = zss[p];

			Dsources[p*numCols + 0] = (float)(cosm.getAngularDiameterDistance(zs)/DIST_MPC);
			size_t useCount = 0;
			for (size_t i = 0 ; i < numPlanes ; i++)
			{
				if (zs > zds[i])
				{
					useCount++;
					Dsources[p*numCols + (i+1)] = (float)(cosm.getAngularDiameterDistance(zds[i], zs)/DIST_MPC);
				}
			}
			usedPlanes[p] = useCount;
		}

		cout << "Point distance info:" << endl;
		for (size_t p = 0 ; p < numPoints ; p++)
		{
			cout << p << ") " << usedPlanes[p] << "|\t";
			for (size_t i = 0 ; i < numCols ; i++)
				cout << "\t" << Dsources[p*numCols + i];
			cout << endl;
		}
		
		// Upload info to GPU
		cl_int err;
		pDevDpoints = clCreateBuffer(getContext(), CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float)*Dsources.size(), Dsources.data(), &err);
		if (err != CL_SUCCESS)
			return "Error uploading distances from lens planes to GPU";

		pDevUsedPlanes = clCreateBuffer(getContext(), CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float)*usedPlanes.size(), usedPlanes.data(), &err);
		if (err != CL_SUCCESS)
			return "Error uploading number of used lens planes to GPU";
		
		return true;
	}

	bool_t allocateAndUploadImages(const std::vector<ImagesDataExtended *> &images, cl_mem &pDevImages, size_t &numPoints)
	{
		vector<float> allCoordinates;
		for (auto img : images)
		{
			int numImages = img->getNumberOfImages();
			for (int i = 0 ; i < numImages ; i++)
			{
				int numPoints = img->getNumberOfImagePoints(i);
				for (int p = 0 ; p < numPoints ; p++)
				{
					Vector2Dd pos = img->getImagePointPosition(i, p);
					pos /= m_angularScale;
					allCoordinates.push_back((float)pos.getX());
					allCoordinates.push_back((float)pos.getY());
				}
			}
		}

		cl_int err;
		pDevImages = clCreateBuffer(getContext(), CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float)*allCoordinates.size(), allCoordinates.data(), &err);
		if (err != CL_SUCCESS)
			return "Error uploading image data to GPU";
		
		numPoints = allCoordinates.size()/2;
		cout << "Uploaded coordinates for " << numPoints << " points to GPU" << endl;

		return true;
	}

	bool_t setupImages(const std::vector<ImagesDataExtended *> &allImages,
	                   const std::vector<ImagesDataExtended *> &shortImages)
	{
		bool_t r;
		if (!(r = allocateAndUploadImages(allImages, m_pDevAllImages, m_numAllImagePoints)))
			return r;

		if (shortImages.size() == 0) // indicates that all images should be used
		{
			m_shortImagesAreAllImages = true;
			m_pDevShortImages = m_pDevAllImages;
			m_numShortImagePoints = m_numAllImagePoints;
		}
		else
		{
			m_shortImagesAreAllImages = false;
			if (!(r = allocateAndUploadImages(shortImages, m_pDevShortImages, m_numShortImagePoints)))
				return r;
		}

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

	double m_angularScale, m_potentialScale;

	bool m_shortImagesAreAllImages;
	cl_mem m_pDevAllImages, m_pDevShortImages;
	size_t m_numAllImagePoints, m_numShortImagePoints;

	cl_mem m_pDevDmatrix;
	cl_mem m_pDevDpointsAll, m_pDevDpointsShort;
	cl_mem m_pDevUsedPlanesAll, m_pDevUsedPlanesShort;

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
											 pParams->getCosmology())))
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
