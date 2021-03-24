#include "lensgamultiobjectivecrossover.h"
#include "lensgaindividual.h"
#include "lensgafitnesscomparison.h"
#include <mogal2/basicnondominatedsetcreator.h>
#include <mogal2/fitnessbasedduplicateremoval.h>

using namespace std;
using namespace errut;

namespace grale
{

LensGAMultiObjectiveCrossover::LensGAMultiObjectiveCrossover(double beta, bool elitism, bool includeBest, double crossoverRate,
            const shared_ptr<mogal2::RandomNumberGenerator> &rng,
            bool allowNegative,
            const shared_ptr<mogal2::GenomeMutation> &mutation,
            size_t numObjectives)
    : LensGACrossoverBase(beta, elitism, includeBest, crossoverRate, rng, allowNegative, mutation),
      m_sortedPop(make_shared<mogal2::BasicNonDominatedSetCreator>(make_shared<LensGAFitnessComparison>(), numObjectives),
                  make_shared<mogal2::FitnessBasedDuplicateRemoval>(make_shared<LensGAFitnessComparison>(), numObjectives) )
{
}

bool_t LensGAMultiObjectiveCrossover::sortCheck(const std::shared_ptr<mogal2::Population> &population)
{
    bool_t r;
    if (!(r = m_sortedPop.check(*population)))
        return "Error in ND sort check: " + r.getErrorString();
    return true;
}

bool_t LensGAMultiObjectiveCrossover::sort(std::shared_ptr<mogal2::Population> &population, size_t targetPopulationSize)
{
    bool_t r;
    if (!(r = m_sortedPop.processPopulation(population, targetPopulationSize)))
        return "Error in ND sort: " + r.getErrorString();
    return true;
}

size_t LensGAMultiObjectiveCrossover::elitism(shared_ptr<mogal2::Population> &population, shared_ptr<mogal2::Population> &newPop)
{
    // cout << "Population is: " << endl;
    // for (auto &i : population->individuals())
    // {
    //     auto *pInd = static_cast<LensGAIndividual*>(i.get());
    //     cout << "  " << pInd->fitnessRef().toString() << " " << pInd->m_parent1 << " " << pInd->m_parent2 << " " << pInd->m_ownIndex << endl;
    // }

    size_t populationSize = population->size();
	size_t eliteCount = (size_t)(((double)(populationSize)*0.005)+0.5); // TODO: just copied from old algorithm
	auto &ndSet = getBestIndividuals();
    size_t nonDominatedSetSize = ndSet.size();
	
	if (eliteCount < 1)
	 	eliteCount = 1;
	if (eliteCount > nonDominatedSetSize)
	 	eliteCount = nonDominatedSetSize;

	size_t mutOffset = 0;

	if ((m_bestWithoutMutation || m_bestWithMutation) && eliteCount > 0)
	{
	 	vector<int> bestGenomeIndices(nonDominatedSetSize);

        auto common = [&bestGenomeIndices, nonDominatedSetSize,&newPop, &ndSet, eliteCount, this](auto extraCode, const string &name)
        {
	 		for (size_t i = 0 ; i < nonDominatedSetSize ; i++)
	 			bestGenomeIndices[i] = i;

	 		for (size_t i = 0 ; i < eliteCount ; i++)
	 		{
				size_t randomIndex = pickRandomNumber(nonDominatedSetSize-i);
				size_t idx = bestGenomeIndices[randomIndex];

				bestGenomeIndices[randomIndex] = bestGenomeIndices[nonDominatedSetSize-i-1];
			
				// TODO: this 0 as the parent index is not correct. perhaps a better
				//       system should be thought of?
                auto ind = ndSet[idx]->createCopy();
                LensGAIndividual &ind2 = static_cast<LensGAIndividual&>(*ind);
                ind2.m_parent1 = 0;
                ind2.m_parent2 = -1;
                // std::cout << name << ind2.fitnessRef().toString() << std::endl;

                newPop->append(ind);
                
                extraCode();
	 		}
        };

	 	if (m_bestWithoutMutation)
	 	{
             common([&mutOffset](){
                mutOffset++;
             },"Elitism I of ");
	 	}

	 	if (m_bestWithMutation)
            common([](){},"Elitism II of ");
	}

    return mutOffset;
}

LensGAIndividual *LensGAMultiObjectiveCrossover::pickParent(const shared_ptr<mogal2::Population> &population)
{
    size_t s = pickBetaDistIndex(m_sortedPop.getNumberOfSets());
    size_t i = pickRandomNumber(m_sortedPop.getSetSize(s));
    return static_cast<LensGAIndividual*>(m_sortedPop.getIndividual(s, i).get());
}

void LensGAMultiObjectiveCrossover::pickParentsRaw(const std::shared_ptr<mogal2::Population> &population, LensGAIndividual **pParent1, LensGAIndividual **pParent2)
{
    size_t s1 = pickBetaDistIndex(m_sortedPop.getNumberOfSets());
    size_t s2 = pickBetaDistIndex(m_sortedPop.getNumberOfSets());
    
    size_t ss1 = m_sortedPop.getSetSize(s1);
    size_t i1 = 0;
    if (ss1 > 1)
        i1 = pickRandomNumber(ss1);
    
    size_t ss2 = m_sortedPop.getSetSize(s2);
    size_t i2 = 0;
    if (ss2 > 1)
        i2 = pickRandomNumber(ss2);

    // cout << "s1 = " << s1 << " i1 = " << i1 << " s2 = " << s2 << " i2 = " << i2 << endl;
    *pParent1 = static_cast<LensGAIndividual*>(m_sortedPop.getIndividual(s1, i1).get());
    *pParent2 = static_cast<LensGAIndividual*>(m_sortedPop.getIndividual(s2, i2).get());
}

}