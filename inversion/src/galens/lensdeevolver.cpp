#include "lensdeevolver.h"
#include "lensgaindividual.h"

using namespace errut;
using namespace std;

namespace grale
{

LensDEEvolver::LensDEEvolver(
		const shared_ptr<eatk::RandomNumberGenerator> &rng,
		const shared_ptr<eatk::DifferentialEvolutionMutation> &mut,
		const shared_ptr<eatk::DifferentialEvolutionCrossover> &cross,
		const shared_ptr<eatk::FitnessComparison> &fitComp) 
	: eatk::JADEEvolver(rng, mut, cross, fitComp) // TODO: other parameters
{
}

bool_t LensDEEvolver::check(const std::shared_ptr<eatk::Population> &population)
{
	for (auto &i : population->individuals())
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

	return eatk::JADEEvolver::check(population);
}

// We need to override this function to copy the calculated scale factors to the genomes
bool_t LensDEEvolver::createNewPopulation(size_t generation, shared_ptr<eatk::Population> &population, size_t targetPopulationSize)
{
	for (auto &i : population->individuals())
	{
		auto &g = static_cast<grale::LensGAGenome &>(i->genomeRef());
		auto &f = static_cast<grale::LensGAFitness &>(i->fitnessRef());
		g.m_scaleFactor = f.m_scaleFactor;
	}
	return eatk::JADEEvolver::createNewPopulation(generation, population, targetPopulationSize);
}

} // end namespace
