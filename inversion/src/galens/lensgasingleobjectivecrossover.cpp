#include "lensgasingleobjectivecrossover.h"
#include "lensgaindividual.h"
#include "lensgafitnesscomparison.h"
#include <cmath>

using namespace std;
using namespace errut;

namespace grale
{

LensGASingleFitnessCrossover::LensGASingleFitnessCrossover(double beta, bool elitism, bool includeBest, double crossoverRate,
            const shared_ptr<mogal2::RandomNumberGenerator> &rng,
            bool allowNegative,
            const shared_ptr<mogal2::GenomeMutation> &mutation)
    : m_beta(beta), m_bestWithoutMutation(elitism), m_bestWithMutation(includeBest), m_crossoverRate(crossoverRate),
        m_rng(rng), m_sortedPop(make_shared<LensGAFitnessComparison>()), m_cross(rng, allowNegative),
        m_mutation(mutation)
{

}

bool_t LensGASingleFitnessCrossover::check(const shared_ptr<mogal2::Population> &population)
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

bool_t LensGASingleFitnessCrossover::createNewPopulation(size_t generation, shared_ptr<mogal2::Population> &population, size_t targetPopulationSize)
{
    bool_t r;
    if (generation == 0)
    {
        if (!(r = m_sortedPop.check(*population)))
            return "Error checking sorter: " + r.getErrorString();
    }
    if (population->size() != targetPopulationSize)
        return "Expecting size " + to_string(targetPopulationSize) + " but got " + to_string(population->size());

    // TODO: do this in another function?
    // Scale factor is calculated and stored in fitness, copy it back to genome
    copyScaleFactorFromFitnessToGenome(population);
    
    if (!(r = m_sortedPop.processPopulation(population, targetPopulationSize)))
        return "Error sorting population: " + r.getErrorString();

    shared_ptr<mogal2::Population> newPop = make_shared<mogal2::Population>();

    size_t mutOffset = elitism(population, newPop);

    if (!(r = crossover(generation, population, newPop)))
        return r;

    if (!(r = mutation(mutOffset, newPop)))
        return r;

    swap(population, newPop);
    
    return true;
}

void LensGASingleFitnessCrossover::copyScaleFactorFromFitnessToGenome(shared_ptr<mogal2::Population> &population)
{
    for (auto &i : population->individuals())
    {
        auto &g = static_cast<LensGAGenome &>(i->genomeRef());
        auto &f = static_cast<LensGAFitness &>(i->fitnessRef());
        g.m_scaleFactor = f.m_scaleFactor;
    }
}

size_t LensGASingleFitnessCrossover::elitism(shared_ptr<mogal2::Population> &population, shared_ptr<mogal2::Population> &newPop)
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

errut::bool_t LensGASingleFitnessCrossover::crossover(size_t generation, shared_ptr<mogal2::Population> &population, shared_ptr<mogal2::Population> &newPop)
{
    auto pickParent = [&population, this]()
    {
        double x = m_rng->getRandomDouble();
        double val = (1.0-pow(x, 1.0/(1.0+m_beta)))*((double)population->size());
        int r = (int)val;
        if (r >= (int)population->size())
                r = (int)population->size() - 1;
        return r;
    };

    auto getLensIndividual = [&population](size_t i) -> const LensGAIndividual &
    {
        return static_cast<LensGAIndividual&>(*population->individual(i));
    };

    // In original version, an attempt was made to prevent inbreeding
    auto pickParentIndices = [&pickParent, &getLensIndividual](int &index1, int &index2)
    {
        bool ok;
        int count = 0;

        do 
        {
            ok = false;

            index1 = pickParent();
            index2 = pickParent();

            // prevent inbreeding

            int a1 = getLensIndividual(index1).m_parent1;
            int a2 = getLensIndividual(index1).m_parent2;
            int b1 = getLensIndividual(index2).m_parent1;
            int b2 = getLensIndividual(index2).m_parent2;

            if (a1 < 0 || b1 < 0) // one of them is a brand new genome
                ok = true;
            else
            {
                if (!(a1 == b1 || a1 == b2 || b1 == a2 || ((a2 >= 0) && a2 == b2)))
                    ok = true;
            }
            count++;
        } while (count < 10 && !ok);
    };

    vector<shared_ptr<mogal2::Genome>> parents(2);
    vector<shared_ptr<mogal2::Genome>> offspring;
    bool_t r;

    while (newPop->size() < population->size())
    {
        if (m_rng->getRandomDouble() < m_crossoverRate)
        {
            int index1 = -1, index2 = -1;
            pickParentIndices(index1, index2);
            parents[0] = population->individual(index1)->genome();
            parents[1] = population->individual(index2)->genome();
            
            if (!(r = m_cross.generateOffspring(parents, offspring)))
                return "Error in crossover: " + r.getErrorString();

            for (auto &g : offspring)
            {
                auto ind = make_shared<LensGAIndividual>(g, population->individual(0)->fitness()->createCopy(false), generation);
                ind->m_parent1 = index1;
                ind->m_parent2 = index2;
                newPop->append(ind);
            }
        }
        else // clone
        {
            int index = pickParent();
            auto ind = population->individual(index)->createCopy();

            LensGAIndividual &lensInd = static_cast<LensGAIndividual&>(*ind);
            lensInd.m_parent1 = index;
            lensInd.m_parent2 = -1;
            newPop->append(ind);
        }
    }
    return true;
}

bool_t LensGASingleFitnessCrossover::mutation(size_t mutOffset, shared_ptr<mogal2::Population> &newPop)
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