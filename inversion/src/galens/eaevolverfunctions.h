#pragma once

#include "graleconfig.h"
#include "lensgaindividual.h"
#include <eatk/population.h>
#include <eatk/vectorgenomefitness.h>

namespace grale
{

inline errut::bool_t evolverCheck(const eatk::Population &population, bool freeform)
{
	if (population.size() == 0)
		return "Empty population";

	if (freeform)
	{
		for (auto &i : population.individuals())
		{
			auto *ind = dynamic_cast<const grale::LensGAIndividual *>(i.get());
			if (!ind)
				return "Each individual should be of type LensGAIndividual";

			auto *g = dynamic_cast<const grale::LensGAGenome *>(i->genomePtr());
			if (!g)
				return "Each individual should have a LensGAGenome";

			auto *f = dynamic_cast<const grale::LensGAFitness *>(i->fitnessPtr());
			if (!f)
				return "Each individual should have a LensGAFitness";
		}
	}
	else
	{
		for (auto &i : population.individuals())
		{
			auto *g = dynamic_cast<const eatk::FloatVectorGenome *>(i->genomePtr());
			if (!g)
				return "Each individual should have a FloatVectorGenome";

			auto *f = dynamic_cast<const eatk::FloatVectorFitness *>(i->fitnessPtr());
			if (!f)
				return "Each individual should have a FloatVectorFitness";
		}
	}

	return true;
}

inline void copyScaleFactorFromFitnessToGenome(const eatk::Population &population)
{
	for (auto &i : population.individuals())
	{
		assert(dynamic_cast<grale::LensGAGenome *>(i->genomePtr()));
		assert(dynamic_cast<grale::LensGAFitness *>(i->fitnessPtr()));

		auto &g = static_cast<grale::LensGAGenome &>(i->genomeRef());
		auto &f = static_cast<grale::LensGAFitness &>(i->fitnessRef());
		g.m_scaleFactor = f.m_scaleFactor;
	}
}

}
