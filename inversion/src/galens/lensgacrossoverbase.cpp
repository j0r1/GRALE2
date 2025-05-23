#include "lensgacrossoverbase.h"
#include "lensgaindividual.h"
#include "eaevolverfunctions.h"
#include "utils.h"
#include <cmath>
#include <limits>
#include <cstdlib>

using namespace std;
using namespace errut;

namespace grale
{

LensGACrossoverBase::LensGACrossoverBase(double beta, bool elitism, bool includeBest, double crossoverRate,
			const shared_ptr<eatk::RandomNumberGenerator> &rng,
			bool allowNegative,
			const shared_ptr<eatk::GenomeMutation> &mutation)
	: m_beta(beta), m_bestWithoutMutation(elitism), m_bestWithMutation(includeBest), m_crossoverRate(crossoverRate),
		m_rng(rng), m_cross(rng, allowNegative),
		m_mutation(mutation)
{
}

bool_t LensGACrossoverBase::check(const shared_ptr<eatk::Population> &population)
{
	bool_t r;

	if (!(r = evolverCheck(*population, true)))
		return r;

	vector<shared_ptr<eatk::Genome>> testParents = { population->individual(0)->genome(), population->individual(0)->genome() };
	if (!(r = m_cross.check(testParents)))
		return "Error in genome crossover check: " + r.getErrorString();

	if (!(r = m_mutation->check(population->individual(0)->genomeRef())))
		return "Error checking mutation: " + r.getErrorString();

	return true;
}

bool_t LensGACrossoverBase::createNewPopulation(size_t generation, shared_ptr<eatk::Population> &population, size_t targetPopulationSize)
{
	// Scale factor is calculated and stored in fitness, copy it back to genome
	// Note: do this before the checkDumpLoad call, so that that call can also check
	//       for NaN (before scale factor is calculated/transferred it is also set
	//       to NaN)
	copyScaleFactorFromFitnessToGenome(*population);

	m_popDump.checkDumpLoad(generation, *population);

	bool_t r;
	if (generation == 0)
	{
		if (!(r = sortCheck(population)))
			return r;
	}
	if (population->size() != targetPopulationSize)
		return "Expecting size " + to_string(targetPopulationSize) + " but got " + to_string(population->size());
	
	if (!(r = sort(population, targetPopulationSize)))
		return "Error sorting population: " + r.getErrorString();

	copyPopulationIndex(population);

	shared_ptr<eatk::Population> newPop = make_shared<eatk::Population>();

	size_t mutOffset = elitism(population, newPop);

	if (!(r = crossover(generation, population, newPop)))
		return r;

	if (!(r = mutation(mutOffset, newPop)))
		return r;

	swap(population, newPop);
	
	return true;
}

void LensGACrossoverBase::copyPopulationIndex(const std::shared_ptr<eatk::Population> &population)
{
	for (size_t i = 0 ; i < population->size() ; i++)
	{
		LensGAIndividual *pInd = static_cast<LensGAIndividual *>(population->individual(i).get());
		assert(pInd);

		pInd->m_ownIndex = (int)i;
	}
}

inline bool areParentsParentsDifferent(const LensGAIndividual &ind1, const LensGAIndividual &ind2)
{
	const LensGAGenome &parent1 = static_cast<const LensGAGenome &>(ind1.genomeRef());
	const LensGAGenome &parent2 = static_cast<const LensGAGenome &>(ind2.genomeRef());
	int a1 = parent1.m_parent1;
	int a2 = parent1.m_parent2;
	int b1 = parent2.m_parent1;
	int b2 = parent2.m_parent2;
	
	if (a1 < 0 || b1 < 0)
		return true;

	if (!(a1 == b1 || a1 == b2 || b1 == a2 || ((a2 >= 0) && a2 == b2)))
		return true;

	return false;
}

void LensGACrossoverBase::pickParentsNoInbreed(const shared_ptr<eatk::Population> &population, 
												 LensGAIndividual **pParent1, LensGAIndividual **pParent2)
{
	bool ok = false;
	int count = 0;
	// In original version, an attempt was made to prevent inbreeding

	do 
	{
		pickParentsRaw(population, pParent1, pParent2);
		// cout << "Trying parents " << (*pParent1)->fitness()->toString() << " and " << (*pParent2)->fitness()->toString() << endl;

		// prevent inbreeding
		ok = areParentsParentsDifferent(**pParent1, **pParent2);
		count++;
	} while (count < 10 && !ok);
	// cout << "Accepted after " << count << " tries" << endl;
}

errut::bool_t LensGACrossoverBase::crossover(size_t generation, shared_ptr<eatk::Population> &population, shared_ptr<eatk::Population> &newPop)
{
	vector<shared_ptr<eatk::Genome>> parents(2);
	vector<shared_ptr<eatk::Genome>> offspring;
	bool_t r;

	while (newPop->size() < population->size())
	{
		if (m_rng->getRandomDouble() < m_crossoverRate)
		{
			LensGAIndividual *pParent1 = nullptr, *pParent2 = nullptr;
			pickParentsNoInbreed(population, &pParent1, &pParent2);
			assert(pParent1 && pParent2);

			parents[0] = pParent1->genome();
			parents[1] = pParent2->genome();
			
			if (!(r = m_cross.generateOffspring(parents, offspring)))
				return "Error in crossover: " + r.getErrorString();

			for (auto &g : offspring)
			{
				auto ind = make_shared<LensGAIndividual>(g, population->individual(0)->fitness()->createCopy(false), generation);
				LensGAGenome &genome = static_cast<LensGAGenome&>(ind->genomeRef());
				genome.m_parent1 = pParent1->m_ownIndex;
				genome.m_parent2 = pParent2->m_ownIndex;
				// Own offset will be set later
				newPop->append(ind);
			}
		}
		else // clone
		{
			LensGAIndividual *pParent = pickParent(population);
			assert(pParent);

			auto ind = pParent->createCopy();

			LensGAIndividual &lensInd = static_cast<LensGAIndividual&>(*ind);
			LensGAGenome &g = static_cast<LensGAGenome&>(lensInd.genomeRef());
			g.m_parent1 = pParent->m_ownIndex;
			g.m_parent2 = -1;
			// Own offset will be set later
			newPop->append(ind);
		}
	}
	return true;
}

bool_t LensGACrossoverBase::mutation(size_t mutOffset, shared_ptr<eatk::Population> &newPop)
{
	bool_t r;

	for (size_t i = mutOffset ; i < newPop->size() ; i++)
	{
		bool isChanged = false;

		auto &ind = newPop->individual(i);

		if (!(r = m_mutation->mutate(ind->genomeRef(), isChanged)))
			return "Error mutating genome: " + r.getErrorString();
		if (isChanged)
			ind->fitness()->setCalculated(false);
	}
	return true;
}

}
