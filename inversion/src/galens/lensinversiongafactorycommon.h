#pragma once

#include "graleconfig.h"
#include "randomnumbergenerator.h"
#include "lensinversionparameterssingleplanecpu.h"
#include "lensinversiongenome.h"
#include "lensfitnessobject.h"
#include "vector2d.h"
#include <mogal/gafactorydefaults.h>
#include <vector>
#include <list>
#include <memory>
#include <fstream>
#include <sstream>

namespace grale
{

class LensInversionGAFactoryParamsSinglePlaneCPU;
class ImagesDataExtended;
class GravitationalLens;
class LensFitnessObject;
class ConfigurationParameters;

// NOTE: the virtual inheritance is again very important!
class GRALE_IMPORTEXPORT LensInversionGAFactoryCommon : public virtual mogal::GAFactory
{
public:
	LensInversionGAFactoryCommon();
	~LensInversionGAFactoryCommon();

	mogal::Genome *createNewGenome() const;

	size_t getMaximalFitnessSize() const							{ return sizeof(float)*(1+GRIDLENSINVERSIONGENOMEBASE_MAXFITNESSCOMP); }
	size_t getMaximalGenomeSize() const								{ return (m_numBasisFunctions+m_numSheetValues)*sizeof(float); }

	bool writeGenome(serut::SerializationInterface &si, const mogal::Genome *g) const;
	bool readGenome(serut::SerializationInterface &si, mogal::Genome **g) const;
	bool writeGenomeFitness(serut::SerializationInterface &si, const mogal::Genome *g) const;
	bool readGenomeFitness(serut::SerializationInterface &si, mogal::Genome *g) const;
	bool writeCommonGenerationInfo(serut::SerializationInterface &si) const;
	bool readCommonGenerationInfo(serut::SerializationInterface &si);

	bool hasFloatingPointFitnessValues() const 						{ return true; }

	const mogal::RandomNumberGenerator *getRandomNumberGenerator() const { return &m_rndGen; }
	bool allowNegativeValues() const								{ return m_allowNegativeValues; }
	bool useLogarithmicScaleSearch() const { return true; }

	virtual float getChanceMultiplier() = 0;
	virtual bool useAbsoluteMutation() = 0;
	virtual float getMutationAmplitude() = 0;

	// TODO: vector<shared_prt<GravitationalLens>>, createLenses
	virtual GravitationalLens *createLens(const LensInversionGenome &genome, std::string &errStr) const = 0;
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

protected:
	virtual LensFitnessObject *createFitnessObject() = 0; // implemented in module
	// This should at least set the number of fitness components
	virtual bool subInit(LensFitnessObject *pFitnessObject) = 0;

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
	void onCurrentBest(const std::list<mogal::Genome *> &bestGenomes) override;
	void onGeneticAlgorithmStart() override;

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
	void onSortedPopulation(const std::vector<mogal::GenomeWrapper> &population) override;
	std::fstream m_scaleSearchFileStream;
	std::stringstream m_scaleSearchStringStream;
	std::vector<std::pair<float,float>> m_searchedPoints;
};

} // end namespace

