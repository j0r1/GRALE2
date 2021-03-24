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
    // TODO
    return 0;
}

LensGAIndividual *LensGAMultiObjectiveCrossover::pickParent(const shared_ptr<mogal2::Population> &population)
{
    size_t s = pickBetaDistIndex(m_sortedPop.getNumberOfSets());
    size_t setSize = m_sortedPop.getSetSize(s);
    assert(setSize > 0);
    size_t i = (size_t)(m_rng->getRandomDouble() * (double)setSize);
    if (i >= setSize)
        i = setSize - 1;
    assert(i >= 0 && i < setSize);
    return static_cast<LensGAIndividual*>(m_sortedPop.getIndividual(s, i).get());
}

}