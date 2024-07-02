#include "lensnsga2evolver.h"
#include "lensgaindividual.h"

using namespace errut;
using namespace std;
using namespace eatk;

namespace grale
{

LensNSGA2Evolver::LensNSGA2Evolver(const shared_ptr<RandomNumberGenerator> &rng,
		const shared_ptr<GenomeCrossover> &genomeCrossover,
		const shared_ptr<GenomeMutation> &genomeMutation,
		const shared_ptr<FitnessComparison> &fitComp, size_t numObjectives)
	: eatk::NSGA2Evolver(rng, genomeCrossover, genomeMutation, fitComp, numObjectives)
{
	// TODO
}

// TODO: is copy-paste from DE evolver, move to common code somewhere
inline bool_t evolverCheck(const eatk::Population &population)
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
	return true;
}

// TODO: same
inline void copyScaleFactorFromFitnessToGenome(const eatk::Population &population)
{
	for (auto &i : population.individuals())
	{
		auto &g = static_cast<grale::LensGAGenome &>(i->genomeRef());
		auto &f = static_cast<grale::LensGAFitness &>(i->fitnessRef());
		g.m_scaleFactor = f.m_scaleFactor;
	}
}

bool_t LensNSGA2Evolver::check(const std::shared_ptr<eatk::Population> &population)
{
	bool_t r = evolverCheck(*population);
	if (!r)
		return r;
	return eatk::NSGA2Evolver::check(population);
}

bool_t LensNSGA2Evolver::createNewPopulation(size_t generation, std::shared_ptr<eatk::Population> &population, size_t targetPopulationSize)
{
	copyScaleFactorFromFitnessToGenome(*population);

//	m_popDump.checkDumpLoad(generation, *population);

	return eatk::NSGA2Evolver::createNewPopulation(generation, population, targetPopulationSize);
}

}