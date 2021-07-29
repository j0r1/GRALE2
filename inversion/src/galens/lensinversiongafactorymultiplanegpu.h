#pragma once

#include "graleconfig.h"
#include "lensinversiongafactorycommon.h"
#include "randomnumbergenerator.h"
#include "vector2d.h"
#include "lensinversionparametersmultiplanegpu.h"
#include "lensinversionbasislensinfo.h"
#include "plummerlensinfo.h"
#include <vector>
#include <memory>
#include <mutex>

namespace grale
{

class LensInversionParametersMultiPlaneGPU;
class LensInversionGenome;
class ImagesDataExtended;
class GravitationalLens;
class LensFitnessObject;
class ConfigurationParameters;

// NOTE: the virtual inheritance is again very important!
class GRALE_IMPORTEXPORT LensInversionGAFactoryMultiPlaneGPU : public LensInversionGAFactoryCommon
{
public:
	LensInversionGAFactoryMultiPlaneGPU(std::unique_ptr<LensFitnessObject> fitObj);
	~LensInversionGAFactoryMultiPlaneGPU();

	errut::bool_t init(const LensInversionParametersBase &p) override;

	std::unique_ptr<GravitationalLens> createLens(const std::vector<float> &basisFunctionWeights,
	                                      const std::vector<float> &sheetValues,
										  float scaleFactor,
										  std::string &errStr) const override;

private:
	errut::bool_t onNewCalculationStart(size_t genomesInThisThread, size_t genomesInAllThreads) override;
	// These a the EATK functions that we'll override to allow an async calculation
	errut::bool_t startNewCalculation(const eatk::Genome &genome) override;
	errut::bool_t pollCalculate(const eatk::Genome &genome, eatk::Fitness &fitness) override;


	errut::bool_t analyzeLensBasisFunctions(const std::vector<double> redshifts,
								   const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &basisLenses);
	errut::bool_t analyzeSourceImages(const std::vector<std::shared_ptr<ImagesDataExtended>> &sourceImages,
							 const Cosmology &cosmology,
							 std::vector<std::shared_ptr<ImagesDataExtended>> &imagesVector);

	std::unique_ptr<LensInversionParametersMultiPlaneGPU> m_currentParams;

	std::vector<float> m_lensRedshifts;
	std::vector<double> m_basisFunctionMasses;
	std::vector<std::shared_ptr<GravitationalLens>> m_unscaledBasisLenses; // One per lens plane

	std::vector<std::shared_ptr<ImagesDataExtended>> m_images;
	std::vector<ImagesDataExtended *> m_reducedImages, m_shortImages;

	class State
	{
	public:
		State() : m_calculationScheduled(false) { }
		void reset(std::pair<float,float> initialStartStopValues, bool finalCalculation)
		{
			m_nextIteration = 0;
			m_startStopValues = initialStartStopValues;
			m_initialStartStopValues = initialStartStopValues;
			m_isFinalCalculation = finalCalculation;
			m_bestScaleFactor = 1.0f;
		}

		bool m_calculationScheduled;

		std::pair<float,float> m_initialStartStopValues;
		std::pair<float,float> m_startStopValues;
		float m_bestScaleFactor;
		bool m_isFinalCalculation;

		std::vector<std::pair<float,float>> m_steps;
		float m_stepSize;
		int m_nextIteration;
	};

	std::vector<State> m_calcStates;
	std::mutex m_calcStateMutex;
};

} // end namespace
