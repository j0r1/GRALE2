#include "lensgacrossoverbase.h"
#include "lensgaindividual.h"
#include "utils.h"
#include <cmath>
#include <limits>
#include <cstdlib>
#include <serut/fileserializer.h>

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

	auto getDumpInfo = [](const string &keyGen, const string &keyFn, size_t &genNr, string &filename)
	{
		genNr = numeric_limits<size_t>::max();

		if (!std::getenv(keyGen.c_str()))
			return;

		bool_t r;
		int gen;

		if (!(r = getenv(keyGen, gen, 0)))
			cerr << "WARNING: " << keyGen << " value should be positive value: " << r.getErrorString() << endl;
		else
		{
			genNr = (size_t)gen;

			getenv(keyFn, filename);
			if (filename.empty())
				cerr << "WARNING: expecting " << keyFn << " to be set" << endl;
		}
	};

	getDumpInfo("GRALE_DUMPPOP_GENERATION", "GRALE_DUMPPOP_FILENAME", m_dumpPopulationGeneration, m_dumpPopulationFilename);
	getDumpInfo("GRALE_LOADPOP_GENERATION", "GRALE_LOADPOP_FILENAME", m_loadPopulationGeneration, m_loadPopulationFilename);
}

bool_t LensGACrossoverBase::check(const shared_ptr<eatk::Population> &population)
{
	if (population->size() == 0)
		return "Empty population";
	auto &i = population->individual(0);
	if (!dynamic_cast<const LensGAGenome *>(i->genomePtr()))
		return "Genome is of wrong type";
	if (!dynamic_cast<const LensGAFitness *>(i->fitnessPtr()))
		return "Fitness is of wrong type";

	bool_t r;
	vector<shared_ptr<eatk::Genome>> testParents = { population->individual(0)->genome(), population->individual(0)->genome() };
	if (!(r = m_cross.check(testParents)))
		return "Error in genome crossover check: " + r.getErrorString();

	if (!(r = m_mutation->check(population->individual(0)->genomeRef())))
		return "Error checking mutation: " + r.getErrorString();

	return true;
}

void LensGACrossoverBase::dumpPopulation(const eatk::Population &population, const std::string &filename)
{
	serut::FileSerializer fSer;

	if (!fSer.open(filename, serut::FileSerializer::WriteOnly))
	{
		cerr << "WARNING: can't open dump population file '" << m_dumpPopulationFilename << ": " << fSer.getErrorString() << endl;
		return;
	}

	cerr << "DEBUG: writing current population to " << filename << endl;

	fSer.writeInt32(population.size());

	for (auto &ind : population.individuals())
	{
		auto pInd = dynamic_cast<const LensGAIndividual *>(ind.get());
		if (!pInd)
		{
			cerr << "WARNING: individual is of incorrect type, can't dump to file" << endl;
			return;
		}
		pInd->write(fSer);
	}
}

void LensGACrossoverBase::loadPopulation(eatk::Population &population, const std::string &filename)
{
	serut::FileSerializer fSer;

	if (!fSer.open(filename, serut::FileSerializer::ReadOnly))
	{
		cerr << "WARNING: can't open dump population file '" << m_dumpPopulationFilename << ": " << fSer.getErrorString() << endl;
		return;
	}

	cerr << "DEBUG: loading current population from " << filename << endl;

	int32_t num = 0;
	fSer.readInt32(&num);
	cerr << "DEBUG: reading " << num << " individuals" << endl;

	std::shared_ptr<eatk::Individual> refInd = population.individual(0);
	std::vector<shared_ptr<eatk::Individual>> newPop;

	bool_t r;

	for (int32_t i = 0 ; i < num ; i++)
	{
		auto newInd = refInd->createCopy();
		auto pInd = dynamic_cast<LensGAIndividual *>(newInd.get());
		if (!(r = pInd->read(fSer)))
		{
			cerr << "WARNING: can't read an individual: " << r.getErrorString() << endl;
			return;
		}
		newPop.push_back(newInd);
	}

	// swap pops
	population.clear();
	for (auto &i : newPop)
		population.append(i);
}

bool_t LensGACrossoverBase::createNewPopulation(size_t generation, shared_ptr<eatk::Population> &population, size_t targetPopulationSize)
{
	if (generation == m_loadPopulationGeneration && !m_loadPopulationFilename.empty())
		loadPopulation(*population, m_loadPopulationFilename);

	if (generation == m_dumpPopulationGeneration && !m_dumpPopulationFilename.empty())
		dumpPopulation(*population, m_dumpPopulationFilename);

	bool_t r;
	if (generation == 0)
	{
		if (!(r = sortCheck(population)))
			return r;
	}
	if (population->size() != targetPopulationSize)
		return "Expecting size " + to_string(targetPopulationSize) + " but got " + to_string(population->size());

	// TODO: do this in another function?
	// Scale factor is calculated and stored in fitness, copy it back to genome
	copyScaleFactorFromFitnessToGenome(population);
	
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

void LensGACrossoverBase::copyScaleFactorFromFitnessToGenome(const shared_ptr<eatk::Population> &population)
{
	for (auto &i : population->individuals())
	{
		auto &g = static_cast<LensGAGenome &>(i->genomeRef());
		auto &f = static_cast<LensGAFitness &>(i->fitnessRef());
		g.m_scaleFactor = f.m_scaleFactor;
	}
}

inline bool areParentsParentsDifferent(const LensGAIndividual &parent1, const LensGAIndividual &parent2)
{
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
				ind->m_parent1 = pParent1->m_ownIndex;
				ind->m_parent2 = pParent2->m_ownIndex;
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
			lensInd.m_parent1 = pParent->m_ownIndex;
			lensInd.m_parent2 = -1;
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
