#include "lensgagenomemutation.h"
#include "lensgaindividual.h"
#include "constants.h"
#include "mathfunctions.h"

using namespace std;
using namespace errut;

namespace grale
{

LensGAGenomeMutation::LensGAGenomeMutation(const shared_ptr<eatk::RandomNumberGenerator> &rng, float chanceMultiplier,
                    bool allowNegativeValues, float mutationAmplitude, bool absoluteMutation)
    : m_rng(rng), m_chanceMultiplier(chanceMultiplier), m_allowNegative(allowNegativeValues),
      m_mutationAmplitude(mutationAmplitude), m_absoluteMutation(absoluteMutation)
{
}

bool_t LensGAGenomeMutation::check(const eatk::Genome &genome)
{
    if (!dynamic_cast<const LensGAGenome*>(&genome))
        return "Genome is of wrong type";

    return true;
}

bool_t LensGAGenomeMutation::mutate(eatk::Genome &genome, bool &isChanged)
{
    LensGAGenome &g = static_cast<LensGAGenome&>(genome);

    size_t numBasisFunctions = g.m_weights.size();
    float chance = m_chanceMultiplier/((float)numBasisFunctions);

    auto absVal = [](auto x) { return ABS(x); };
    auto noChange = [](auto x) { return x; };
    auto getMaxVal = [&g, numBasisFunctions](auto fn)
    {
        float maxVal = 0;

        for (int i = 0 ; i < numBasisFunctions ; i++)
        {
            float x = fn(g.m_weights[i]);
            if (maxVal < x)
                maxVal = x;
        }
        return maxVal;
    };

    float maxVal = (m_allowNegative)?getMaxVal(absVal):getMaxVal(noChange);

    // In the past, this 'maxVal' value was rescaled to 0.5, the unit range then 
    // corresponds to twice this. This means that the unit based mutation amplitude
    // needs to be scaled by 2*maxVal
    float rescale = 2.0f*maxVal;
    float mutationAmplitude = m_mutationAmplitude * rescale;

    auto chanceSetUniform = [&isChanged, chance, this, rescale](float &x, float mult, float offset)
    {
        if (m_rng->getRandomFloat() < chance)
        {
            x = (m_rng->getRandomFloat()*mult - offset)*rescale;
            isChanged = true;
        }
    };

    auto chanceSetUniformAllScalableBasisFunctions = [this, &g, chance, numBasisFunctions, chanceSetUniform](float mult, float offset)
    {
        for (int i = 0 ; i < numBasisFunctions ; i++)
            chanceSetUniform(g.m_weights[i], mult, offset);
    };

    auto chanceSetSmallDiff = [&isChanged, chance, this, mutationAmplitude](float &target, float yMin, float yMax)
    {
        if (m_rng->getRandomFloat() < chance)
        {
            // allow larger mutations with smaller probablility
            // p(x) = (2/Pi)*1/(x^2+1)
            // cfr anomalous diffusion
            float p = m_rng->getRandomFloat()*2.0f-1.0f;
            float x = TAN(p*(float)(CONST_PI/4.0))*mutationAmplitude;
            float y = x+target;

            if (y < yMin)
                y = yMin;
            else if (y > yMax)
                y = yMax;
            
            target = y;

            isChanged = true;
        }
    };

    auto chanceSetSmallDiffAllScalableBasisFunctions = [this, chance, &g, numBasisFunctions, chanceSetSmallDiff](float yMin, float yMax)
    {
        for (int i = 0 ; i < numBasisFunctions ; i++)
            chanceSetSmallDiff(g.m_weights[i], yMin, yMax);
    };

    if (m_absoluteMutation)
    {
        if (m_allowNegative)
            chanceSetUniformAllScalableBasisFunctions(2.0f, 1.0f);
        else
            chanceSetUniformAllScalableBasisFunctions(1.0f, 0.0f);
    }
    else
    {
        if (m_allowNegative)
            chanceSetSmallDiffAllScalableBasisFunctions(-rescale, rescale);
        else
            chanceSetSmallDiffAllScalableBasisFunctions(0.0f, rescale);
    }

    auto chanceSetSimpleUniform = [&isChanged, chance, this](float &x)
    {
        if (m_rng->getRandomFloat() < chance)
        {
            x = m_rng->getRandomFloat();
            isChanged = true;
        }
    };

    if (m_absoluteMutation)
    {
        for (auto &v : g.m_sheets)
            chanceSetSimpleUniform(v);
    }
    else
    {
        for (auto &v : g.m_sheets)
            chanceSetSmallDiff(v, 0.0f, 1.0f);
    }

    return true;
}

}
