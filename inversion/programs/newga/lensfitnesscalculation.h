#pragma once

#include "lensinversiongafactorycommon.h"
#include <mogal2/genomefitness.h>

class LensFitnessCalculation : public mogal2::GenomeFitnessCalculation
{
public:
	LensFitnessCalculation(grale::LensInversionGAFactoryCommon &factory) : m_pFactory(&factory) { }
	LensFitnessCalculation(const std::shared_ptr<grale::LensInversionGAFactoryCommon> &factory)
	{
		m_spFactory = factory; // keep a reference to the lib
		m_pFactory = factory.get();
	}

	~LensFitnessCalculation() { }

	errut::bool_t calculate(const mogal2::Genome &genome, mogal2::Fitness &fitness)
	{
		const grale::LensGAGenome &g = dynamic_cast<const grale::LensGAGenome &>(genome);
		grale::LensGAFitness &f = dynamic_cast<grale::LensGAFitness &>(fitness);
		
		// Anything that needs to get communicated back needs to be in the fitness, so
		// we'll store the scalefactor in the fitness and transfer it to the genome in a
		// later stage
		if (!m_pFactory->calculateFitness(g.m_weights, g.m_sheets, f.m_scaleFactor, f.m_fitnesses.data()))
			return "Error calculating fitness: " + m_pFactory->getErrorString();

		return true;
	}
private:
	grale::LensInversionGAFactoryCommon *m_pFactory;
	std::shared_ptr<grale::LensInversionGAFactoryCommon> m_spFactory;
};
