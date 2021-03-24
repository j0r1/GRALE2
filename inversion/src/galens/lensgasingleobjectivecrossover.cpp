#include "lensgasingleobjectivecrossover.h"
#include "lensgaindividual.h"
#include "lensgafitnesscomparison.h"

using namespace std;
using namespace errut;

namespace grale
{

LensGASingleObjectiveCrossover::LensGASingleObjectiveCrossover(double beta, bool elitism, bool includeBest, double crossoverRate,
            const shared_ptr<mogal2::RandomNumberGenerator> &rng,
            bool allowNegative,
            const shared_ptr<mogal2::GenomeMutation> &mutation)
    : LensGACrossoverBase(beta, elitism, includeBest, crossoverRate, rng, allowNegative, mutation),
        m_sortedPop(make_shared<LensGAFitnessComparison>())
{
}

bool_t LensGASingleObjectiveCrossover::sortCheck(const std::shared_ptr<mogal2::Population> &population)
{
    bool_t r;
    if (!(r = m_sortedPop.check(*population)))
        return "Error checking sorter: " + r.getErrorString();
    return true;
}

bool_t LensGASingleObjectiveCrossover::sort(std::shared_ptr<mogal2::Population> &population, size_t targetPopulationSize)
{
    bool_t r;

    if (!(r = m_sortedPop.processPopulation(population, targetPopulationSize)))
        return r;

    return true;
}

size_t LensGASingleObjectiveCrossover::elitism(shared_ptr<mogal2::Population> &population, shared_ptr<mogal2::Population> &newPop)
{
    auto appendBest = [&newPop, &population]()
    {
        auto ind = population->individual(0)->createCopy();
        LensGAIndividual &i = static_cast<LensGAIndividual&>(*ind);
        i.m_parent1 = 0;
        i.m_parent2 = -1;
        newPop->append(ind);
    };

    size_t mutOffset = 0;
    if (m_bestWithoutMutation)
    {
        mutOffset = 1;
        appendBest();
    }
    if (m_bestWithMutation)
        appendBest();

    return mutOffset;
}

LensGAIndividual *LensGASingleObjectiveCrossover::pickParent(const shared_ptr<mogal2::Population> &population)
{
    size_t r = pickBetaDistIndex(population->size());

    return static_cast<LensGAIndividual*>(population->individual(r).get());
}


}