#include "lensgagenomecrossover.h"
#include "lensgaindividual.h"

using namespace std;
using namespace errut;

namespace grale
{

LensGAGenomeCrossover::LensGAGenomeCrossover(const shared_ptr<eatk::RandomNumberGenerator> &rng, bool allowNegative)
    : m_rng(rng), m_allowNegative(allowNegative)
{
}

bool_t LensGAGenomeCrossover::check(const vector<shared_ptr<eatk::Genome>> &parents)
{
    if (parents.size() != 2)
        return "Expecting two parents";
    if (dynamic_cast<LensGAGenome*>(parents[0].get()) == 0)
        return "Genome is of wrong type";
    return true;
}

bool_t LensGAGenomeCrossover::generateOffspring(const vector<shared_ptr<eatk::Genome>> &parents,
                                        vector<shared_ptr<eatk::Genome>> &generatedOffspring)
{
    assert(parents.size() == 2);
    LensGAGenome *pParents[2] = {
        static_cast<LensGAGenome*>(parents[1].get()),
        static_cast<LensGAGenome*>(parents[0].get())
    };

    size_t numBasisFunctions = pParents[0]->m_weights.size();

    generatedOffspring.clear();
    auto offspring = parents[0]->createCopy(false);
    generatedOffspring.push_back(offspring);
    LensGAGenome *pOff = static_cast<LensGAGenome*>(offspring.get());

    vector<float> &newBasisFunctionsWeights = pOff->m_weights;

    auto pickParent = [pParents, this]()
    {
        float x = m_rng->getRandomFloat();
        const LensGAGenome *pParent = pParents[(x < 0.5f)?1:0];
        return pParent;
    };

    auto genomeUniformCrossover = [&newBasisFunctionsWeights, numBasisFunctions, &pickParent](auto check)
    {
        for (int i = 0 ; i < numBasisFunctions ; i++)
        {
            const LensGAGenome *pParent = pickParent();
            newBasisFunctionsWeights[i] = pParent->m_scaleFactor * pParent->m_weights[i];
            check(newBasisFunctionsWeights[i]);
        }
    };

    auto noChange = [](float &x) { };
    auto clamp = [](float &x) { if (x < 0) x = 0; };

    if (m_allowNegative)
        genomeUniformCrossover(noChange);
    else
        genomeUniformCrossover(clamp);

    vector<float> &newSheetValues = pOff->m_sheets;
    for (size_t i = 0 ; i < newSheetValues.size() ; i++)
        newSheetValues[i] = pickParent()->m_sheets[i];

    return true;
}

}