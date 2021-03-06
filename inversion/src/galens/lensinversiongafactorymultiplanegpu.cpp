#include "lensinversiongafactorymultiplanegpu.h"
#include "multiplanecontainer.h"
#include "compositelens.h"
#include "masssheetlens.h"
#include "plummerlens.h"
#include "utils.h"
#include <assert.h>

using namespace std;

namespace grale
{

LensInversionGAFactoryMultiPlaneGPU::LensInversionGAFactoryMultiPlaneGPU()
{

}

LensInversionGAFactoryMultiPlaneGPU::~LensInversionGAFactoryMultiPlaneGPU()
{

}

bool LensInversionGAFactoryMultiPlaneGPU::init(const LensInversionParametersBase *p)
{
	auto pParams = dynamic_cast<const LensInversionParametersMultiPlaneGPU *>(p);
	if (!pParams)
	{
		setErrorString("Specified parameters are not of the correct type");
		return false;
	}

	if (!getenv("GRALE_MPCUDA_LIBRARY", m_libraryPath))
	{
		setErrorString("Environment variable GRALE_MPCUDA_LIBRARY for the helper library is not set");
		return false;
	}

	if (!analyzeLensBasisFunctions(pParams->getLensRedshifts(), pParams->getBasisLenses()))
		return false;

	m_images.clear();
	if (!analyzeSourceImages(pParams->getSourceImages(), pParams->getCosmology(), m_images))
		return false;

	// These are the ones to actually use
	m_reducedImages.clear();
	m_shortImages.clear();
	if (!initializeLensFitnessObject(numeric_limits<double>::quiet_NaN(), m_images, pParams->getFitnessObjectParameters(),
									 m_reducedImages, m_shortImages))
		return false;

	// We're going to init CUDA later, only when we actually need it. The way the code
	// is structured now, extra factory instances could be created (well at least one),
	// causing more CUDA inits to happen than needed
	m_cudaInitAttempted = false;
	m_cudaInitialized = false;

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
	if (!setCommonParameters(m_sheetMultipliers.size(), pParams->getMaximumNumberOfGenerations(),
						pParams->getAllowNegativeWeights(), m_basisFunctionMasses,
						pParams->getMassEstimate(), pParams->getMassEstimate(),
						pParams->getMassScaleSearchParameters()))
		return false;

	m_currentParams = make_unique<LensInversionParametersMultiPlaneGPU>(*pParams);
	return true;
}

bool LensInversionGAFactoryMultiPlaneGPU::analyzeLensBasisFunctions(const vector<double> redshifts,
						  const vector<vector<shared_ptr<LensInversionBasisLensInfo>>> &basisLenses)
{
	if (redshifts.size() == 0)
	{
		setErrorString("No lens planes present");
		return false;
	}

	if (redshifts.size() != basisLenses.size())
	{
		setErrorString("Incompatible number of redshifts (" + to_string(redshifts.size()) + ") and lens planes (" + to_string(basisLenses.size()) + ")");
		return false;
	}

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
		{
			setErrorString("A lens plane without basis functions is present");
			return false;
		}

		// Just create the right sizes for these vectors
		m_scaledPlaneWeights.push_back(vector<float>(plane.size(), 0));
		m_basePlaneWeights.push_back(vector<float>(plane.size(), 0));

		m_basisLenses.push_back(vector<PlummerLensInfo>());
		auto &plummers = m_basisLenses.back();

		for (auto &bl : plane)
		{
			if (bl->m_relevantLensingMass < 0)
			{
				setErrorString("A basis lens was found to have a negative strong lensing mass");
				return false;
			}

			const PlummerLens *pPlummerLens = dynamic_cast<const PlummerLens *>(bl->m_pLens.get());
			if (!pPlummerLens)
			{
				setErrorString("Not all basis functions appear to be Plummer lens models");
				return false;
			}

			const PlummerLensParams *pParams = dynamic_cast<const PlummerLensParams *>(pPlummerLens->getLensParameters());
			if (!pParams)
			{
				setErrorString("Unexpected: couldn't get Plummer lens parameters");
				return false;
			}

			plummers.push_back(PlummerLensInfo { pParams->getLensMass(), pParams->getAngularWidth(), bl->m_center });
			m_basisFunctionMasses.push_back(bl->m_relevantLensingMass);
		}
	}

	return true;
}

bool LensInversionGAFactoryMultiPlaneGPU::analyzeSourceImages(const vector<shared_ptr<ImagesDataExtended>> &sourceImages,
						 const Cosmology &cosmology,
						 vector<shared_ptr<ImagesDataExtended>> &imagesVector)
{
	if (sourceImages.size() == 0)
	{
		setErrorString("No images present");
		return false;
	}

	imagesVector.clear();
	
	int count = 0;
	for (auto img : sourceImages)
	{
		count++;

		double z;
		if (!img->hasExtraParameter("z") || !img->getExtraParameter("z", z))
		{
			setErrorString("Images data set " + to_string(count) + " does not contain the redshift 'z' parameter");
			return false;
		}

		if (img->getDds() != 0 || img->getDs() != 0)
		{
			setErrorString("Images data set " + to_string(count) + " does not have Dds and Ds set to zero, required for the multi-plane inversion");
			return false;
		}

		// Calculate Ds from redshift, and set it in a copy
		double Ds = cosmology.getAngularDiameterDistance(z);
		shared_ptr<ImagesDataExtended> newImg = make_shared<ImagesDataExtended>(*img);
		newImg->setDs(Ds);

		imagesVector.push_back(newImg);
	}
	return true;
}

GravitationalLens *LensInversionGAFactoryMultiPlaneGPU::createLens(const std::vector<float> &basisFunctionWeights,
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
				setErrorString("Unable to add basis function to composite lens: " + planeLensParams.getErrorString());
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
				setErrorString("Could not initialize a mass sheet lens: " + sheetLens.getErrorString());
				return nullptr;
			}
			if (!planeLensParams.addLens(1.0, Vector2Dd(0, 0), 0, sheetLens))
			{
				setErrorString("Unable to add sheet lens to composite lens: " + planeLensParams.getErrorString());
				return nullptr;
			}
		}

		auto compLens = make_shared<CompositeLens>();
		if (!compLens->init(Dd, &planeLensParams))
		{
			setErrorString("Unable to create a composite lens for a lens plane: " + compLens->getErrorString());
			return nullptr;
		}
		containerParams.add(compLens, z);
	}

	assert(basisFunctionWeightIdx == basisFunctionWeights.size());
	
	unique_ptr<MultiPlaneContainer> containerLens = make_unique<MultiPlaneContainer>();
	if (!containerLens->init(0, &containerParams)) // Dd must be set to 0!
	{
		setErrorString("Unable to initialize the multi-plane container: " + containerLens->getErrorString());
		return nullptr;
	}

	return containerLens.release();
}

void LensInversionGAFactoryMultiPlaneGPU::convertGenomeSheetValuesToDensities(const vector<float> &sheetValues,
																			  vector<float> &sheetDensities) const
{
	assert(sheetValues.size() == m_sheetMultipliers.size());
	assert(sheetValues.size() == sheetDensities.size());

	for (size_t i = 0 ; i < sheetValues.size() ; i++)
		sheetDensities[i] = sheetValues[i]*m_sheetMultipliers[i];
}

bool LensInversionGAFactoryMultiPlaneGPU::initializeNewCalculation(const std::vector<float> &basisFunctionWeights, const std::vector<float> &sheetValues)
{
	convertGenomeSheetValuesToDensities(sheetValues, m_sheetDensities);

	// Store the basis function weights with the same multi-plane structure
	int idx = 0;
	for (auto &planeWeights : m_basePlaneWeights)
	{
		for (auto &w : planeWeights)
		{
			assert(idx < basisFunctionWeights.size());
			w = basisFunctionWeights[idx++];
		}
	}
	assert(idx == basisFunctionWeights.size());
	return true;
}

bool LensInversionGAFactoryMultiPlaneGPU::scaleWeights(float scaleFactor)
{
	assert(m_basePlaneWeights.size() == m_scaledPlaneWeights.size());
	for (size_t i = 0 ; i  < m_scaledPlaneWeights.size() ; i++)
	{
		auto &basePlane = m_basePlaneWeights[i];
		auto &scaledPlane = m_scaledPlaneWeights[i];

		assert(basePlane.size() == scaledPlane.size());
		for (size_t j = 0 ; j < scaledPlane.size() ; j++)
			scaledPlane[j] = scaleFactor*basePlane[j];
	}

	return true;
}

bool LensInversionGAFactoryMultiPlaneGPU::calculateMassScaleFitness(float scaleFactor, float &fitness)
{
	if (!checkCUDAInit())
		return false;

	if (!scaleWeights(scaleFactor))
		return false;

	assert(m_cudaBpShort.get());
	if (!m_cudaBpShort->calculateSourcePositions(m_scaledPlaneWeights, m_sheetDensities))
	{
		setErrorString("Error back-projecting images: " + m_cudaBpShort->getErrorString());
		return false;
	}

	LensFitnessObject &fitnessObject = getFitnessObject();
	if (!fitnessObject.calculateMassScaleFitness(*m_cudaBpShort, fitness))
	{
		setErrorString("Unable to calculate fitness for mass scale: " + fitnessObject.getErrorString());
		return false;
	}
	return true;
}

bool LensInversionGAFactoryMultiPlaneGPU::calculateTotalFitness(float scaleFactor, float *pFitnessValues)
{
	if (!checkCUDAInit())
		return false;

	if (!scaleWeights(scaleFactor))
		return false;

	assert(m_cudaBpFull.get());
	if (!m_cudaBpFull->calculateSourcePositions(m_scaledPlaneWeights, m_sheetDensities))
	{
		setErrorString("Error back-projecting images: " + m_cudaBpFull->getErrorString());
		return false;
	}

	LensFitnessObject &fitnessObject = getFitnessObject();
	if (!fitnessObject.calculateOverallFitness(*m_cudaBpFull, pFitnessValues))
	{
		setErrorString("Unable to calculate full fitness: " + fitnessObject.getErrorString());
		return false;
	}

	return true;
}

bool LensInversionGAFactoryMultiPlaneGPU::checkCUDAInit()
{
	if (m_cudaInitialized)
		return true;

	if (m_cudaInitAttempted)
	{
		setErrorString("CUDA initialization has failed previously");
		return false;
	}
	m_cudaInitAttempted = true;

	const LensInversionParametersMultiPlaneGPU *pParams = m_currentParams.get();
	if (!pParams)
	{
		setErrorString("Unexpected: parameters is null, we should already have checked this");
		return false;
	}

	auto createBackProjector = [this, pParams](auto imgs, int devIdx=-1) -> shared_ptr<MPCUDABackProjector>
	{
		vector<float> sourceRedshifts;
		for (auto i : imgs)
		{
			double z;
			if (!i->hasExtraParameter("z") || !i->getExtraParameter("z", z))
			{
				setErrorString("Unexpected: no 'z' parameter in images data instance");
				return nullptr;
			}
			sourceRedshifts.push_back((float)z);
		}

		if (devIdx < 0)
			devIdx = pParams->getDeviceIndex();

		auto bp = make_shared<MPCUDABackProjector>();
		if (!bp->init(m_libraryPath, devIdx, pParams->getCosmology(), 
					  m_lensRedshifts, m_basisLenses, sourceRedshifts, imgs))
		{
			setErrorString("Unable to initialize CUDA based backprojector: " + bp->getErrorString());
			return nullptr;
		}
		return bp;
	};

	// Allocate backprojector for reducedImages
	m_cudaBpFull = createBackProjector(m_reducedImages);
	if (!m_cudaBpFull.get())
		return false;

	if (m_shortImages.size() > 0) // We can speed up the scale search with a smaller set of images
	{
		// Make sure that we're using the same device, we won't be calculating
		// short and full versions at the same time
		m_cudaBpShort = createBackProjector(m_shortImages, m_cudaBpFull->getDeviceIndex());
		if (!m_cudaBpShort.get())
			return false;
	}
	else // use the same
		m_cudaBpShort = m_cudaBpFull;

	m_cudaInitialized = true;
	return true;
}

} // end namespace