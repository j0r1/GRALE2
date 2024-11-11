#include "lensnsga2evolver.h"
#include "lensgaindividual.h"
#include "eaevolverfunctions.h"

using namespace errut;
using namespace std;
using namespace eatk;

namespace grale
{

LensNSGA2Evolver::LensNSGA2Evolver(bool freeform, const shared_ptr<RandomNumberGenerator> &rng,
		const shared_ptr<GenomeCrossover> &genomeCrossover,
		const shared_ptr<GenomeMutation> &genomeMutation,
		const shared_ptr<FitnessComparison> &fitComp, size_t numObjectives)
	: eatk::NSGA2Evolver(rng, genomeCrossover, genomeMutation, fitComp, numObjectives),
	  m_popDump(PopulationDump(freeform)),
	  m_freeForm(freeform)
{
}

bool_t LensNSGA2Evolver::check(const std::shared_ptr<eatk::Population> &population)
{
	bool_t r = evolverCheck(*population, m_freeForm);
	if (!r)
		return r;
	return eatk::NSGA2Evolver::check(population);
}

bool_t LensNSGA2Evolver::createNewPopulation(size_t generation, std::shared_ptr<eatk::Population> &population, size_t targetPopulationSize)
{
	if (m_freeForm)
		copyScaleFactorFromFitnessToGenome(*population);

	m_popDump.checkDumpLoad(generation, *population);

	return eatk::NSGA2Evolver::createNewPopulation(generation, population, targetPopulationSize);
}

}
