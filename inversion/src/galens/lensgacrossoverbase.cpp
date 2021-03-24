#include "lensgacrossoverbase.h"
#include "lensgaindividual.h"
#include <cmath>

using namespace std;
using namespace errut;

namespace grale
{

LensGACrossoverBase::LensGACrossoverBase(double beta, bool elitism, bool includeBest, double crossoverRate,
            const shared_ptr<mogal2::RandomNumberGenerator> &rng,
            bool allowNegative,
            const shared_ptr<mogal2::GenomeMutation> &mutation)
    : m_beta(beta), m_bestWithoutMutation(elitism), m_bestWithMutation(includeBest), m_crossoverRate(crossoverRate),
        m_rng(rng), m_cross(rng, allowNegative),
        m_mutation(mutation)
{

}

bool_t LensGACrossoverBase::check(const shared_ptr<mogal2::Population> &population)
{
    if (population->size() == 0)
        return "Empty population";
    auto &i = population->individual(0);
    if (!dynamic_cast<const LensGAGenome *>(i->genomePtr()))
        return "Genome is of wrong type";
    if (!dynamic_cast<const LensGAFitness *>(i->fitnessPtr()))
        return "Fitness is of wrong type";

    bool_t r;
    vector<shared_ptr<mogal2::Genome>> testParents = { population->individual(0)->genome(), population->individual(0)->genome() };
    if (!(r = m_cross.check(testParents)))
        return "Error in genome crossover check: " + r.getErrorString();

    if (!(r = m_mutation->check(population->individual(0)->genomeRef())))
        return "Error checking mutation: " + r.getErrorString();

    return true;
}

bool_t LensGACrossoverBase::createNewPopulation(size_t generation, shared_ptr<mogal2::Population> &population, size_t targetPopulationSize)
{
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

    shared_ptr<mogal2::Population> newPop = make_shared<mogal2::Population>();

    size_t mutOffset = elitism(population, newPop);

    if (!(r = crossover(generation, population, newPop)))
        return r;

    if (!(r = mutation(mutOffset, newPop)))
        return r;

    swap(population, newPop);
    
    return true;
}

void LensGACrossoverBase::copyPopulationIndex(const std::shared_ptr<mogal2::Population> &population)
{
    for (size_t i = 0 ; i < population->size() ; i++)
    {
        LensGAIndividual *pInd = static_cast<LensGAIndividual *>(population->individual(i).get());
        assert(pInd);

        pInd->m_ownIndex = (int)i;
    }
}

void LensGACrossoverBase::copyScaleFactorFromFitnessToGenome(const shared_ptr<mogal2::Population> &population)
{
    for (auto &i : population->individuals())
    {
        auto &g = static_cast<LensGAGenome &>(i->genomeRef());
        auto &f = static_cast<LensGAFitness &>(i->fitnessRef());
        g.m_scaleFactor = f.m_scaleFactor;
    }
}

void LensGACrossoverBase::pickParents(const shared_ptr<mogal2::Population> &population, 
                                                 LensGAIndividual **pParent1, LensGAIndividual **pParent2)
{
    bool ok;
    int count = 0;
    // In original version, an attempt was made to prevent inbreeding

    do 
    {
        ok = false;

        *pParent1 = pickParent(population);
        *pParent2 = pickParent(population);

        // prevent inbreeding

        int a1 = (*pParent1)->m_parent1;
        int a2 = (*pParent1)->m_parent2;
        int b1 = (*pParent2)->m_parent1;
        int b2 = (*pParent2)->m_parent2;

        if (a1 < 0 || b1 < 0) // one of them is a brand new genome
            ok = true;
        else
        {
            if (!(a1 == b1 || a1 == b2 || b1 == a2 || ((a2 >= 0) && a2 == b2)))
                ok = true;
        }
        count++;
    } while (count < 10 && !ok);

}

errut::bool_t LensGACrossoverBase::crossover(size_t generation, shared_ptr<mogal2::Population> &population, shared_ptr<mogal2::Population> &newPop)
{
    vector<shared_ptr<mogal2::Genome>> parents(2);
    vector<shared_ptr<mogal2::Genome>> offspring;
    bool_t r;

    while (newPop->size() < population->size())
    {
        if (m_rng->getRandomDouble() < m_crossoverRate)
        {
            LensGAIndividual *pParent1 = nullptr, *pParent2 = nullptr;
            pickParents(population, &pParent1, &pParent2);
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

bool_t LensGACrossoverBase::mutation(size_t mutOffset, shared_ptr<mogal2::Population> &newPop)
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
