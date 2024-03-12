#include "lensdemutationcrossover.h"
#include "lensgaindividual.h"
#include <limits>

using namespace std;
using namespace errut;

namespace grale
{

bool_t LensDEMutation::check(const eatk::Genome &g)
{
	if (!dynamic_cast<const LensGAGenome*>(&g))
		return "Genome is of wrong type";
	return true;
}

shared_ptr<eatk::Genome> LensDEMutation::mutate(const vector<const eatk::Genome*> &genomes, const vector<double> &weights)
{
	assert(genomes.size() > 0);
	assert(genomes.size() == weights.size());

	shared_ptr<eatk::Genome> result = genomes[0]->createCopy();
	LensGAGenome &g0 = static_cast<LensGAGenome&>(*result);

	// Initialize everything to zero, so we can safely add things
	g0.m_weights.assign(g0.m_weights.size(), 0);
	g0.m_sheets.assign(g0.m_sheets.size(), 0);

	// We can't set this to NaN, since it will be used in the next crossover step
	// But we've incorporated the scalefactors, so this should just be 1
	g0.m_scaleFactor = 1.0;

	for (size_t j = 0 ; j < genomes.size() ; j++)
	{
		float f = (float)weights[j];
		const LensGAGenome &g1 = static_cast<const LensGAGenome&>(*genomes[j]);

		// Add the weights, need to take the scale factor into account as well!
		assert(g0.m_weights.size() == g1.m_weights.size());
		for (size_t i = 0 ; i < g0.m_weights.size() ; i++)
			g0.m_weights[i] += f * g1.m_weights[i] * g1.m_scaleFactor;

		// Add the sheet values, no scale factor here
		assert(g0.m_sheets.size() == g1.m_sheets.size());
		for (size_t i = 0 ; i < g0.m_sheets.size() ; i++)
			g0.m_sheets[i] += f * g1.m_sheets[i];
	}

	return result;
}

LensDECrossover::LensDECrossover(const shared_ptr<eatk::RandomNumberGenerator> &rng, bool allowNeg)
	: m_rng(rng), m_allowNegative(allowNeg) 
{
}

errut::bool_t LensDECrossover::check(const eatk::Genome &g)
{
	if (!dynamic_cast<const LensGAGenome*>(&g))
		return "Genome is of wrong type";
	return true;
}

errut::bool_t LensDECrossover::crossover(double CR, eatk::Genome &mutantDest, const eatk::Genome &origVector)
{
	LensGAGenome &g0 = static_cast<LensGAGenome&>(mutantDest);
	const LensGAGenome &g1 = static_cast<const LensGAGenome&>(origVector);

	size_t numWeigths = g0.m_weights.size();
	size_t numSheets = g0.m_sheets.size();
	size_t totalSize = numWeigths + numSheets;
	size_t rndIdx = ((size_t)m_rng->getRandomUint32())%(totalSize);

	auto getValue = [numWeigths,numSheets](const LensGAGenome &g, size_t i)
	{
		if (i < numWeigths)
			return g.m_weights[i] * g.m_scaleFactor; // Take scale factor into account!
		
		i -= numWeigths;
		assert(i < numSheets);
		return g.m_sheets[i]; // no scale factor here
	};

	// Note that no scale factor will be used here
	auto setValue = [numWeigths,numSheets](LensGAGenome &g, size_t i, float value)
	{
		if (i < numWeigths)
			g.m_weights[i] = value;
		else
		{
			i -= numWeigths;
			assert(i < numSheets);
			g.m_sheets[i] = value;
		}
	};

	for (size_t i = 0 ; i < totalSize ; i++)
	{
		float val = numeric_limits<float>::quiet_NaN();

		// We make sure to get the value in either case, so we get something that incorporated
		// the scale factor; we can then store the value again
		if (i != rndIdx && m_rng->getRandomDouble() > CR)
			val = getValue(g1, i);
		else
			val = getValue(g0, i);

		setValue(g0, i, val);
	}

	// Reset the scale factor
	g0.m_scaleFactor = numeric_limits<float>::quiet_NaN();

	if (!m_allowNegative) // Enforce bounds on the weights
	{
		// g1 should already be within bounds, we won't check this (only in assert)!
		for (size_t i = 0 ; i < numWeigths ; i++)
		{
			if (g0.m_weights[i] < 0)
				g0.m_weights[i] = g1.m_weights[i]/2.0f;

			assert(g0.m_weights[i] >= 0);
			assert(g1.m_weights[i] >= 0);
		}
	}

	// Always enforce bounds on the sheets - TODO: what is done in the normal GA?
	for (size_t i = 0 ; i < numSheets ; i++)
	{
		// g1 should already be within bounds, we won't check this (only in assert)!
		if (g0.m_sheets[i] < 0)
			g0.m_sheets[i] = g1.m_sheets[i]/2.0f;

		assert(g0.m_sheets[i] >= 0);
		assert(g1.m_sheets[i] >= 0);
	}

	return true;
}

} // end namespace
