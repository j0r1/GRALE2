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

class GRALE_IMPORTEXPORT LensInversionGAFactoryCommon : public LensGAGenomeCalculator, public errut::ErrorBase
{
public:
	LensInversionGAFactoryCommon(std::unique_ptr<LensFitnessObject> fitObj);
	~LensInversionGAFactoryCommon();

	bool allowNegativeValues() const override								{ return m_allowNegativeValues; }
	bool useLogarithmicScaleSearch() const { return true; }

	errut::bool_t createLens(const LensGAGenome &genome, std::unique_ptr<GravitationalLens> &lens) const;
	errut::bool_t calculate(const eatk::Genome &genome, eatk::Fitness &fitness);

	virtual bool initializeNewCalculation(const std::vector<float> &basisFunctionWeights, const std::vector<float> &sheetValues) = 0;
	virtual bool calculateMassScaleFitness(float scaleFactor, float &fitness) = 0;
	virtual bool calculateTotalFitness(float scaleFactor, float *pFitnessValues) = 0;

	size_t getNumberOfBasisFunctions() const override { return (size_t)m_numBasisFunctions; }
	size_t getNumberOfSheets() const override { return (size_t)m_numSheetValues; }

	size_t getNumberOfObjectives() const override { return m_fitnessObject->getNumberOfFitnessComponents(); }

	// TODO: use general log method
	void sendMessage(const std::string &s) const;
protected:
	virtual std::unique_ptr<GravitationalLens> createLens(const std::vector<float> &basisFunctionWeights,
	                                      const std::vector<float> &sheetValues,
										  float scaleFactor,
										  std::string &errStr) const = 0;

	bool setCommonParameters(int numSheetValues,
							 bool allowNeg,
							 const std::vector<double> &basisFunctionMasses,
							 double massUnit, double targetMass,
							 const ScaleSearchParameters &searchParams);

	LensFitnessObject &getFitnessObject() { return *(m_fitnessObject.get()); }

	bool initializeLensFitnessObject(double z_d,
		const std::vector<std::shared_ptr<ImagesDataExtended>> &images,
		const ConfigurationParameters *pFitnessObjectParameters,
		std::vector<ImagesDataExtended*> &reducedImages,
		std::vector<ImagesDataExtended*> &shortImages);
private:
	bool calculateFitness(const std::vector<float> &basisFunctionWeights,
						  const std::vector<float> &sheetValues,
						  float &scaleFactor,
						  float *pFitnessValues);

	// void onGeneticAlgorithmStart(); // TODO: re-enable this

	int m_numBasisFunctions, m_numSheetValues;
	bool m_allowNegativeValues;
	ScaleSearchParameters m_massScaleSearchParams;

	// In units of the specified mass unit
	std::vector<float> m_basisFunctionMasses;
	float m_targetMass;

	std::unique_ptr<LensFitnessObject> m_fitnessObject;

	// To be able to debug the scale factor search
	std::fstream m_scaleSearchFileStream;
	std::stringstream m_scaleSearchStringStream;
	std::vector<std::pair<float,float>> m_searchedPoints;
};

} // end namespace

