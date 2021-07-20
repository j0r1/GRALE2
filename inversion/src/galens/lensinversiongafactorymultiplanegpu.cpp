#include "lensinversiongafactorymultiplanegpu.h"
#include "multiplanecontainer.h"
#include "compositelens.h"
#include "masssheetlens.h"
#include "plummerlens.h"
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

bool_t LensInversionGAFactoryMultiPlaneGPU::initializeNewCalculation(const std::vector<float> &basisFunctionWeights, const std::vector<float> &sheetValues)
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

void LensInversionGAFactoryMultiPlaneGPU::scaleWeights(float scaleFactor)
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
}

bool_t LensInversionGAFactoryMultiPlaneGPU::calculateMassScaleFitness(float scaleFactor, float &fitness)
{
	return "ERROR: calculateMassScaleFitness is no longer used";
}

bool_t LensInversionGAFactoryMultiPlaneGPU::calculateTotalFitness(float scaleFactor, float *pFitnessValues)
{
	return "ERROR: calculateTotalFitness is no longer used";
}

bool_t LensInversionGAFactoryMultiPlaneGPU::startNewCalculation(const eatk::Genome &genome)
{
	return "TODO: implement startNewCalculation";
}

bool_t LensInversionGAFactoryMultiPlaneGPU::pollCalculate(const eatk::Genome &genome, eatk::Fitness &fitness)
{
	return "TODO: implement pollCalculate";
}

} // end namespace
