#include "lensinversiongafactorycommon.h"
#include "gridlensinversiongenomebase.h"
#include <limits>
#include <iostream>
#include <sstream>

#include "debugnew.h"

using namespace std;

namespace grale
{

LensInversionGAFactoryCommon::LensInversionGAFactoryCommon()
{
  	m_numMasses = 0;
    m_numSheetValues = 0;
    m_maxGenerations = 0;
	m_allowNegativeValues = false;
    m_angularScale = numeric_limits<double>::quiet_NaN();
}

LensInversionGAFactoryCommon::~LensInversionGAFactoryCommon()
{

}

bool LensInversionGAFactoryCommon::setCommonParameters(int numMasses, int numSheetValues, 
                            int maxGenerations,
                            bool allowNeg, double angularScale,
                            const std::vector<float> &massWeights,
                            const ScaleSearchParameters &searchParams)
{
	sendMessage("RNG SEED: " + std::to_string(m_rndGen.getSeed()));

    m_numMasses = numMasses;
    m_numSheetValues = numSheetValues;
    m_maxGenerations = maxGenerations;
    m_allowNegativeValues = allowNeg;
    m_angularScale = angularScale;

    if (numMasses != (int)massWeights.size())
    {
        setErrorString("Specified number of masses is incompatible with number of mass weight");
        return false;
    }
    m_massWeights = massWeights;

    m_massScaleSearchParams = searchParams;
    return true;
}

void LensInversionGAFactoryCommon::sendMessage(const std::string &s)
{
	if (getCurrentAlgorithm())
		GAFactory::sendMessage(s);
	else // queue for later
	{
		cerr << "Queue: " << s << endl;
		m_queuedMessages.push_back(s);
	}
}

void LensInversionGAFactoryCommon::onGeneticAlgorithmStart()
{
	// Send the messages that were previously queued
	for (auto &s : m_queuedMessages)
		GAFactory::sendMessage(s);

	m_queuedMessages.clear();
}

mogal::Genome *LensInversionGAFactoryCommon::createNewGenome() const
{
    assert(m_numMasses > 0);
	return new GridLensInversionGenomeBase(const_cast<LensInversionGAFactoryCommon*>(this), m_numMasses, m_numSheetValues);
}

bool LensInversionGAFactoryCommon::writeGenome(serut::SerializationInterface &si, const mogal::Genome *g) const
{
	const GridLensInversionGenomeBase *g2 = (const GridLensInversionGenomeBase *)g;

	if (!si.writeFloats(g2->getMasses()))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeFloats(g2->getSheetValues()))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool LensInversionGAFactoryCommon::readGenome(serut::SerializationInterface &si, mogal::Genome **g) const
{
    assert(m_numMasses > 0);
	std::vector<float> masses(m_numMasses);
	std::vector<float> sheetValues(m_numSheetValues);

	if (!si.readFloats(masses))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.readFloats(sheetValues))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	// This moves masses and sheetValues
	*g = new GridLensInversionGenomeBase(const_cast<LensInversionGAFactoryCommon*>(this), masses, sheetValues);
	return true;
}

bool LensInversionGAFactoryCommon::writeGenomeFitness(serut::SerializationInterface &si, const mogal::Genome *g) const
{
	const GridLensInversionGenomeBase *g2 = (const GridLensInversionGenomeBase *)g;
	int num = getNumberOfFitnessComponents();
	const float *f = g2->getFitnessValues();

	if (!si.writeFloats(f, num))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeFloat(g2->getScaleFactor()))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool LensInversionGAFactoryCommon::readGenomeFitness(serut::SerializationInterface &si, mogal::Genome *g) const
{
	GridLensInversionGenomeBase *g2 = (GridLensInversionGenomeBase *)g;
	int num = getNumberOfFitnessComponents();
	float x[GRIDLENSINVERSIONGENOMEBASE_MAXFITNESSCOMP+1];

	if (!si.readFloats(x, num+1))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	g2->setFitnessValues(x);
	g2->setScaleFactor(x[num]);
	return true;
}

bool LensInversionGAFactoryCommon::writeCommonGenerationInfo(serut::SerializationInterface &si) const
{
	return true;
}

bool LensInversionGAFactoryCommon::readCommonGenerationInfo(serut::SerializationInterface &si)
{
	return true;
}

void LensInversionGAFactoryCommon::getGenomeCalculationParameters(float &startfactor, float &stopfactor, int &numiterationsteps, int &numiterations, int &numiterationsteps2) const
{
	startfactor = m_massScaleSearchParams.getStartFactor();
	stopfactor = m_massScaleSearchParams.getStopFactor();
	numiterationsteps = m_massScaleSearchParams.getStepsOnFirstIteration();
	numiterations = m_massScaleSearchParams.getNumberOfIterations();
	numiterationsteps2 = m_massScaleSearchParams.getStepsOnSubsequentIterations();
}

void LensInversionGAFactoryCommon::onCurrentBest(const list<mogal::Genome *> &bestGenomes)
{
	stringstream ss;
	ss << "Current best:";
	for (auto g : bestGenomes)
		ss << "( " << g->getFitnessDescription() << ")";

	sendMessage(ss.str());
}

} // end namespace