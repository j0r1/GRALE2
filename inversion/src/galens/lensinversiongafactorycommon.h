#pragma once

#include "graleconfig.h"
#include "randomnumbergenerator.h"
#include "gridlensinversionparameters.h"
#include "gridlensinversiongenomebase.h"
#include "vector2d.h"
#include <mogal/gafactorydefaults.h>
#include <vector>
#include <memory>

namespace grale
{

class GridLensInversionGAFactoryParams;
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
	size_t getMaximalGenomeSize() const							    { return (m_numMasses+m_numSheetValues)*sizeof(float); }

	bool writeGenome(serut::SerializationInterface &si, const mogal::Genome *g) const;
	bool readGenome(serut::SerializationInterface &si, mogal::Genome **g) const;
	bool writeGenomeFitness(serut::SerializationInterface &si, const mogal::Genome *g) const;
	bool readGenomeFitness(serut::SerializationInterface &si, mogal::Genome *g) const;
	bool writeCommonGenerationInfo(serut::SerializationInterface &si) const;
	bool readCommonGenerationInfo(serut::SerializationInterface &si);

	bool hasFloatingPointFitnessValues() const 						{ return true; }

	const mogal::RandomNumberGenerator *getRandomNumberGenerator() const { return &m_rndGen; }
	bool allowNegativeValues() const							    { return m_allowNegativeValues; }
    const float *getMassWeights() const							    { return &(m_massWeights[0]); }
	
	void getGenomeCalculationParameters(float &startfactor, float &stopfactor, int &numiterationsteps, int &numiterations, int &numiterationsteps2) const;
	bool useLogarithmicScaleSearch() const { return true; }

	virtual float getChanceMultiplier() = 0;
	virtual bool useAbsoluteMutation() = 0;
	virtual float getMutationAmplitude() = 0;

    // TODO: vector<shared_prt<GravitationalLens>>, createLenses
	virtual GravitationalLens *createLens(const GridLensInversionGenomeBase &genome, std::string &errStr) const = 0;

	void sendMessage(const std::string &s);

	virtual bool initializeNewCalculation(const std::vector<float> &masses, const std::vector<float> &sheetValues) = 0;
	virtual bool calculateMassScaleFitness(float scaleFactor, float &fitness) = 0;
	virtual bool calculateTotalFitness(float scaleFactor, float *pFitnessValues) = 0;
protected:
	virtual LensFitnessObject *createFitnessObject() = 0;
	// This should at least set the number of fitness components
	virtual bool subInit(LensFitnessObject *pFitnessObject) = 0;

	// The massWeights are the masses of the basis lenses in units of the total mass
    bool setCommonParameters(int numMasses, int numSheetValues, int maxGenerations,
                             bool allowNeg, double angularScale,
                             const std::vector<float> &massWeights,
                             const ScaleSearchParameters &searchParams);

	int getMaximumNumberOfGenerations() const						{ return m_maxGenerations; }
	int getNumberOfMasses() const                                   { return m_numMasses; }
	int getNumberOfSheetValues() const                              { return m_numSheetValues; }
	double getAngularScale() const                                  { return m_angularScale; }
private:
	void onCurrentBest(const std::list<mogal::Genome *> &bestGenomes) override;
	void onGeneticAlgorithmStart() override;

	RandomNumberGenerator m_rndGen;
	int m_numMasses, m_numSheetValues, m_maxGenerations;
	bool m_allowNegativeValues;
    double m_angularScale;
	ScaleSearchParameters m_massScaleSearchParams;
	std::vector<float> m_massWeights;

	std::vector<std::string> m_queuedMessages;
};

} // end namespace

