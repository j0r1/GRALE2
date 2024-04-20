#pragma once

#include <errut/booltype.h>
#include <eatk/individual.h>
#include <vector>

class PreferredIndividualSelector
{
public:
	PreferredIndividualSelector() { }
	virtual ~PreferredIndividualSelector() { }
	virtual errut::bool_t select(const std::vector<std::shared_ptr<eatk::Individual>> &best,
								 std::shared_ptr<eatk::Individual> &selected)
	{
		return "Not implemented";
	}
};

class SubsequentBestIndividualSelector : public PreferredIndividualSelector
{
public:
	SubsequentBestIndividualSelector(size_t numObjectives,
									 const std::shared_ptr<eatk::FitnessComparison> &fitComp)
		: m_numObjectives(numObjectives), m_cmp(fitComp) { }
	~SubsequentBestIndividualSelector() { }
	
	errut::bool_t select(const std::vector<std::shared_ptr<eatk::Individual>> &bestIndividuals,
								 std::shared_ptr<eatk::Individual> &selected) override
	{
		std::vector<std::shared_ptr<eatk::Individual>> genomes = bestIndividuals;
		std::vector<std::shared_ptr<eatk::Individual>> genomes2;
		std::shared_ptr<eatk::Individual> best;

		if (genomes.size() == 0)
			return "No best individuals to select one from";

		for (size_t comp = 0 ; comp < m_numObjectives ; comp++)
		{
			best = genomes[0];

			for (auto &i : genomes)
			{
				if (m_cmp->isFitterThan(i->fitnessRef(), best->fitnessRef(), comp))
					best = i;
			}

			genomes2.clear();

			// Ok, now we know a genome that has the lowest 'comp' fitness, let's see it there
			// are others which perfom equally well
			for (auto &i : genomes)
			{
				// if 'best' is not fitter than 'i' it must have the same fitness with respect to this component
				if (!m_cmp->isFitterThan(best->fitnessRef(), i->fitnessRef(), comp))
					genomes2.push_back(i);
			}

			swap(genomes, genomes2);
			genomes2.clear();
		}
		selected = best;
		return true;
	}
private:
	size_t m_numObjectives;
	std::shared_ptr<eatk::FitnessComparison> m_cmp;
};


