#pragma once

#include "graleconfig.h"
#include "randomnumbergenerator.h"
#include "lensinversionparameterssingleplanecpu.h"
#include "lensfitnessobject.h"
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

// NOTE: the virtual inheritance is again very important!
class GRALE_IMPORTEXPORT LensInversionGAFactoryCommon : public errut::ErrorBase
{
public:
	LensInversionGAFactoryCommon();
	~LensInversionGAFactoryCommon();

	virtual bool init(const LensInversionParametersBase *p) = 0;

	bool allowNegativeValues() const								{ return m_allowNegativeValues; }
	bool useLogarithmicScaleSearch() const { return true; }

	virtual float getChanceMultiplier() = 0;
	virtual bool useAbsoluteMutation() = 0;
	virtual float getMutationAmplitude() = 0;

	// TODO: vector<shared_prt<GravitationalLens>>, createLenses
	virtual GravitationalLens *createLens(const std::vector<float> &basisFunctionWeights,
	                                      const std::vector<float> &sheetValues,
										  float scaleFactor,
										  std::string &errStr) const = 0;

	void sendMessage(const std::string &s);

	bool calculateFitness(const std::vector<float> &basisFunctionWeights,
						  const std::vector<float> &sheetValues,
						  float &scaleFactor,
						  float *pFitnessValues);

	virtual bool initializeNewCalculation(const std::vector<float> &basisFunctionWeights, const std::vector<float> &sheetValues) = 0;
	virtual bool calculateMassScaleFitness(float scaleFactor, float &fitness) = 0;
	virtual bool calculateTotalFitness(float scaleFactor, float *pFitnessValues) = 0;

	size_t getNumberOfBasisFunctions() const { return (size_t)m_numBasisFunctions; }
	size_t getNumberOfSheets() const { return (size_t)m_numSheetValues; }

	LensFitnessObject &getFitnessObject() { return *(m_fitnessObject.get()); }
	int getMaximumNumberOfGenerations() const						{ return m_maxGenerations; }

	size_t getNumberOfFitnessComponents() { return getFitnessObject().getNumberOfFitnessComponents(); }
protected:
	virtual LensFitnessObject *createFitnessObject() = 0; // implemented in module

	bool initializeLensFitnessObject(double z_d,
		const std::vector<std::shared_ptr<ImagesDataExtended>> &images,
		const ConfigurationParameters *pFitnessObjectParameters,
		std::vector<ImagesDataExtended*> &reducedImages,
		std::vector<ImagesDataExtended*> &shortImages);

	bool setCommonParameters(int numSheetValues, int maxGenerations,
							 bool allowNeg,
							 const std::vector<double> &basisFunctionMasses,
							 double massUnit, double targetMass,
							 const ScaleSearchParameters &searchParams);

private:
	void onGeneticAlgorithmStart();

	RandomNumberGenerator m_rndGen;
	int m_numBasisFunctions, m_numSheetValues, m_maxGenerations;
	bool m_allowNegativeValues;
	ScaleSearchParameters m_massScaleSearchParams;

	// In units of the specified mass unit
	std::vector<float> m_basisFunctionMasses;
	float m_targetMass;

	std::vector<std::string> m_queuedMessages;

	std::unique_ptr<LensFitnessObject> m_fitnessObject;

	// To be able to debug the scale factor search
	std::fstream m_scaleSearchFileStream;
	std::stringstream m_scaleSearchStringStream;
	std::vector<std::pair<float,float>> m_searchedPoints;
};

} // end namespace

