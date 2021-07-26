#pragma once

#include "graleconfig.h"
#include "lensinversionparameterssingleplanecpu.h"
#include "lensfitnessobject.h"
#include "lensgagenomecalculator.h"
#include "vector2d.h"
#include <vector>
#include <list>
#include <memory>
#include <fstream>
#include <sstream>

namespace grale
{

class LensInversionParametersBase;
class ImagesDataExtended;
class GravitationalLens;
class LensFitnessObject;
class ConfigurationParameters;

class GRALE_IMPORTEXPORT LensInversionGAFactoryCommon : public LensGAGenomeCalculator
{
public:
	LensInversionGAFactoryCommon(std::unique_ptr<LensFitnessObject> fitObj);
	~LensInversionGAFactoryCommon();

	bool allowNegativeValues() const override								{ return m_allowNegativeValues; }
	bool useLogarithmicScaleSearch() const { return true; }

	errut::bool_t createLens(const LensGAGenome &genome, std::unique_ptr<GravitationalLens> &lens) const;
	errut::bool_t calculate(const eatk::Genome &genome, eatk::Fitness &fitness);

	virtual errut::bool_t initializeNewCalculation(const std::vector<float> &basisFunctionWeights, const std::vector<float> &sheetValues);
	virtual errut::bool_t calculateMassScaleFitness(float scaleFactor, float &fitness);
	virtual errut::bool_t calculateTotalFitness(float scaleFactor, float *pFitnessValues);

	size_t getNumberOfBasisFunctions() const override { return (size_t)m_numBasisFunctions; }
	size_t getNumberOfSheets() const override { return (size_t)m_numSheetValues; }

	size_t getNumberOfObjectives() const override { return m_fitnessObject->getNumberOfFitnessComponents(); }
protected:
	virtual std::unique_ptr<GravitationalLens> createLens(const std::vector<float> &basisFunctionWeights,
	                                      const std::vector<float> &sheetValues,
										  float scaleFactor,
										  std::string &errStr) const = 0;

	errut::bool_t setCommonParameters(int numSheetValues,
							 bool allowNeg,
							 const std::vector<double> &basisFunctionMasses,
							 double massUnit, double targetMass,
							 const ScaleSearchParameters &searchParams);

	LensFitnessObject &getFitnessObject() { return *(m_fitnessObject.get()); }

	errut::bool_t initializeLensFitnessObject(double z_d,
		const std::vector<std::shared_ptr<ImagesDataExtended>> &images,
		const ConfigurationParameters *pFitnessObjectParameters,
		std::vector<ImagesDataExtended*> &reducedImages,
		std::vector<ImagesDataExtended*> &shortImages);

	int getNumberOfCalculationIterations() const { return m_massScaleSearchParams.getNumberOfIterations(); }
	std::pair<float,float> getInitialStartStopValues(const std::vector<float> &basisFunctionWeights) const;
	float getStepsAndStepSize(std::pair<float,float> startStopValue, int iteration, std::vector<std::pair<float,float>> &steps) const;
	void updateStartStopValues(std::pair<float,float> &startStopValue, std::pair<float,float> startStopValue0, 
										float currentBestScaleFactor, float stepsize) const;
private:
	errut::bool_t calculateFitness(const std::vector<float> &basisFunctionWeights,
						  const std::vector<float> &sheetValues,
						  float &scaleFactor,
						  float *pFitnessValues);

	float getScalingMassSum(const std::vector<float> &basisFunctionWeights) const;

	// void onGeneticAlgorithmStart(); // TODO: re-enable this

	int m_numBasisFunctions, m_numSheetValues;
	bool m_allowNegativeValues;
	ScaleSearchParameters m_massScaleSearchParams;

	// In units of the specified mass unit
	std::vector<float> m_basisFunctionMasses;
	float m_targetMass;
	std::vector<std::pair<float,float>> m_tmpSteps; // To avoid reallocation

	std::unique_ptr<LensFitnessObject> m_fitnessObject;

	// To be able to debug the scale factor search
	std::fstream m_scaleSearchFileStream;
	std::stringstream m_scaleSearchStringStream;
	std::vector<std::pair<float,float>> m_searchedPoints;
};

} // end namespace

