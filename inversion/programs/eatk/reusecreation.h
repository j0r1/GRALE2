#pragma once

#include <eatk/individual.h>
#include <eatk/population.h>
#include <list>

class ReuseCreation : public eatk::IndividualCreation
{
public:
	ReuseCreation(const eatk::Population &pop)
	{
		if (pop.size() == 0)
			return;

		m_referenceIndividual = pop.individual(0)->createCopy();
		m_referenceGenome = m_referenceIndividual->genomePtr()->createCopy();
		m_referenceFitness = m_referenceIndividual->fitnessPtr()->createCopy();

		for (auto &i : pop.individuals())
			m_genomePool.push_back(i->genomePtr()->createCopy());
	}

	std::shared_ptr<eatk::Genome> createInitializedGenome() override
	{
		if (m_genomePool.size() == 0)
			return nullptr;
		auto genome = m_genomePool.front();
		m_genomePool.pop_front();
		return genome;
	}

	std::shared_ptr<eatk::Genome> createUnInitializedGenome() override
	{ 
		if (!m_referenceGenome.get())
			return nullptr;
		return m_referenceGenome->createCopy(false);
	}

	std::shared_ptr<eatk::Fitness> createEmptyFitness() override
	{
		if (!m_referenceFitness.get())
			return nullptr;
		return m_referenceFitness->createCopy(false);
	}

	std::shared_ptr<eatk::Individual> createReferenceIndividual() override { return m_referenceIndividual; }
private:
	std::list<std::shared_ptr<eatk::Genome>> m_genomePool;
	std::shared_ptr<eatk::Individual> m_referenceIndividual;
	std::shared_ptr<eatk::Genome> m_referenceGenome;
	std::shared_ptr<eatk::Fitness> m_referenceFitness;
};


