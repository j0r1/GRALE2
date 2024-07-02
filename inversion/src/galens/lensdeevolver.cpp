#include "lensdeevolver.h"
#include "lensgaindividual.h"
#include "eaevolverfunctions.h"

using namespace errut;
using namespace std;

namespace grale
{

LensDEEvolver::LensDEEvolver(
		const shared_ptr<eatk::RandomNumberGenerator> &rng,
		const shared_ptr<eatk::DifferentialEvolutionMutation> &mut,
		double F,
		const shared_ptr<eatk::DifferentialEvolutionCrossover> &cross,
		double CR,
		const shared_ptr<eatk::FitnessComparison> &fitComp, int objectiveNumber, size_t numObjectives,
		const shared_ptr<eatk::NonDominatedSetCreator> &ndCreator,
		bool needStrictlyBetter)
	: eatk::DifferentialEvolutionEvolver(rng, mut, F, cross, CR, fitComp, objectiveNumber, numObjectives, ndCreator, needStrictlyBetter)
{
}

bool_t LensDEEvolver::check(const std::shared_ptr<eatk::Population> &population)
{
	bool_t r = evolverCheck(*population);
	if (!r)
		return r;
	return eatk::DifferentialEvolutionEvolver::check(population);
}

// We need to override this function to copy the calculated scale factors to the genomes
bool_t LensDEEvolver::createNewPopulation(size_t generation, shared_ptr<eatk::Population> &population, size_t targetPopulationSize)
{
	copyScaleFactorFromFitnessToGenome(*population);

	m_popDump.checkDumpLoad(generation, *population);

	return eatk::DifferentialEvolutionEvolver::createNewPopulation(generation, population, targetPopulationSize);
}

LensJADEEvolver::LensJADEEvolver(
		const shared_ptr<eatk::RandomNumberGenerator> &rng,
		const shared_ptr<eatk::DifferentialEvolutionMutation> &mut,
		const shared_ptr<eatk::DifferentialEvolutionCrossover> &cross,
		const shared_ptr<eatk::FitnessComparison> &fitComp,
		int objectiveNumber, // negative means multi-objective
		double p, double c,
		bool useArchive,
		double initMuF,
		double initMuCR,
		size_t numObjectives,
		const shared_ptr<eatk::NonDominatedSetCreator> &ndCreator,
		bool needStrictlyBetter
	) 
	: eatk::JADEEvolver(rng, mut, cross, fitComp, objectiveNumber, p, c, useArchive, initMuF, initMuCR, numObjectives, ndCreator, needStrictlyBetter)
{
}

bool_t LensJADEEvolver::check(const std::shared_ptr<eatk::Population> &population)
{
	bool_t r = evolverCheck(*population);
	if (!r)
		return r;
	return eatk::JADEEvolver::check(population);
}

// We need to override this function to copy the calculated scale factors to the genomes
bool_t LensJADEEvolver::createNewPopulation(size_t generation, shared_ptr<eatk::Population> &population, size_t targetPopulationSize)
{
	copyScaleFactorFromFitnessToGenome(*population);

	m_popDump.checkDumpLoad(generation, *population);

	return eatk::JADEEvolver::createNewPopulation(generation, population, targetPopulationSize);
}

void LensJADEEvolver::onMutationCrossoverSettings(double muF, double muCR) const
{
	//cerr << "DEBUG: muF = " << muF << ", muCR = " << muCR << endl;
}

} // end namespace
