#include "populationdump.h"
#include "lensgaindividual.h"
#include "utils.h"
#include <serut/fileserializer.h>
#include <eatk/vectorgenomefitness.h>
#include <memory>
#include <iostream>
#include <cmath>

using namespace errut;
using namespace std;

namespace grale
{

PopulationDump::PopulationDump(bool freeform) : m_freeform(freeform)
{
	auto getDumpInfo = [](const string &keyGen, const string &keyFn, size_t &genNr, string &filename)
	{
		genNr = numeric_limits<size_t>::max();

		if (!std::getenv(keyGen.c_str()))
			return;

		bool_t r;
		int gen;

		if (!(r = getenv(keyGen, gen, 0)))
			cerr << "WARNING: " << keyGen << " value should be positive value: " << r.getErrorString() << endl;
		else
		{
			genNr = (size_t)gen;

			getenv(keyFn, filename);
			if (filename.empty())
				cerr << "WARNING: expecting " << keyFn << " to be set" << endl;
		}
	};

	getDumpInfo("GRALE_DUMPPOP_GENERATION", "GRALE_DUMPPOP_FILENAME", m_dumpPopulationGeneration, m_dumpPopulationFilename);
	getDumpInfo("GRALE_LOADPOP_GENERATION", "GRALE_LOADPOP_FILENAME", m_loadPopulationGeneration, m_loadPopulationFilename);
}

PopulationDump::~PopulationDump()
{
}

void PopulationDump::checkDumpLoad(size_t generation, eatk::Population &population)
{
	if (generation == m_loadPopulationGeneration && !m_loadPopulationFilename.empty())
		loadPopulation(population, m_loadPopulationFilename);

	if (generation == m_dumpPopulationGeneration && !m_dumpPopulationFilename.empty())
		dumpPopulation(population, m_dumpPopulationFilename);
}

void PopulationDump::dumpPopulation(const eatk::Population &population, const std::string &filename)
{
	serut::FileSerializer fSer;

	if (!fSer.open(filename, serut::FileSerializer::WriteOnly))
	{
		cerr << "WARNING: can't open dump population file '" << m_dumpPopulationFilename << ": " << fSer.getErrorString() << endl;
		return;
	}

	cerr << "DEBUG: writing current population to " << filename << endl;

	fSer.writeInt32(population.size());

	if (m_freeform)
	{
		for (size_t idx = 0 ; idx < population.size() ; idx++)
		{
			auto pInd = dynamic_cast<const LensGAIndividual *>(population.individual(idx).get());
			if (!pInd)
			{
				cerr << "WARNING: individual is of incorrect type, can't dump to file" << endl;
				return;
			}
			pInd->write(fSer);

			auto pGenome = dynamic_cast<const LensGAGenome *>(pInd->genomePtr());
			auto pFitness = dynamic_cast<const LensGAFitness *>(pInd->fitnessPtr());
			if (!pGenome || !pFitness)
				cerr << "WARNING: Genome or fitness is not of the correct type for a LensGAIndividual" << endl;
			else
			{
				pGenome->checkNaN(idx);
				pFitness->checkNaN(idx);
			}
		}
	}
	else
	{
		for (size_t idx = 0 ; idx < population.size() ; idx++)
		{
			auto pInd = population.individual(idx).get();
			auto pGenome = dynamic_cast<const eatk::FloatVectorGenome *>(pInd->genomePtr());
			auto pFitness = dynamic_cast<const eatk::FloatVectorFitness *>(pInd->fitnessPtr());

			if (!pGenome || !pFitness)
				cerr << "WARNING: Genome or fitness is not of the correct type (float vector)" << endl;
			else
			{
				const vector<float> &gv = pGenome->getValues();
				const vector<float> &fv = pFitness->getValues();
				array<int32_t,2> sizes = { (int32_t)gv.size(), (int32_t)fv.size() };

				fSer.writeInt32s(sizes.data(), sizes.size());
				fSer.writeFloats(gv);
				fSer.writeFloats(fv);

				bool genomeNaN = false;
				bool fitnessNaN = false;
				for (auto x : gv)
					if (std::isnan(x))
						genomeNaN = true;
				for (auto x : fv)
					if (std::isnan(x))
						fitnessNaN = true;

				if (genomeNaN)
					cerr << "WARNING: detected NaN in genome of individual " << idx << endl;
				if (fitnessNaN)
					cerr << "WARNING: detected NaN in fitness of individual " << idx << endl;
			}
		}
	}
}

void PopulationDump::loadPopulation(eatk::Population &population, const std::string &filename)
{
	serut::FileSerializer fSer;

	if (!fSer.open(filename, serut::FileSerializer::ReadOnly))
	{
		cerr << "WARNING: can't open dump population file '" << m_dumpPopulationFilename << ": " << fSer.getErrorString() << endl;
		return;
	}

	cerr << "DEBUG: loading current population from " << filename << endl;

	int32_t num = 0;
	fSer.readInt32(&num);
	cerr << "DEBUG: reading " << num << " individuals" << endl;

	std::shared_ptr<eatk::Individual> refInd = population.individual(0);
	std::vector<shared_ptr<eatk::Individual>> newPop;

	bool_t r;

	if (m_freeform)
	{
		for (int32_t i = 0 ; i < num ; i++)
		{
			auto newInd = refInd->createCopy();
			auto pInd = dynamic_cast<LensGAIndividual *>(newInd.get());
			if (!(r = pInd->read(fSer)))
			{
				cerr << "WARNING: can't read an individual: " << r.getErrorString() << endl;
				return;
			}
			newPop.push_back(newInd);
		}
	}
	else
	{
		for (int32_t i = 0 ; i < num ; i++)
		{
			auto newInd = refInd->createCopy();
			auto pGenome = dynamic_cast<eatk::FloatVectorGenome *>(newInd->genomePtr());
			auto pFitness = dynamic_cast<eatk::FloatVectorFitness *>(newInd->fitnessPtr());

			if (!pGenome || !pFitness)
				cerr << "WARNING: Genome or fitness is not of the correct type to load into (float vector)" << endl;
			else
			{
				array<int32_t,2> sizes;
				if (!fSer.readInt32s(sizes.data(), sizes.size()))
					cerr << "WARNING: couldn't read float vector genome/fitness sizes" << endl;
				else
				{
					vector<float> &gv = pGenome->getValues();
					vector<float> &fv = pFitness->getValues();

					if ((size_t)sizes[0] != gv.size() ||
						(size_t)sizes[1] != fv.size())
						cerr << "WARNING: stored genome/fitness sizes don't match in memory ones" << endl;
					else
					{
						if (!fSer.readFloats(gv) || !fSer.readFloats(fv))
							cerr << "WARNING: can't read stored float values for genome or fitness" << endl;
						else
							newPop.push_back(newInd);
					}
				}
			}
		}
	}

	// swap pops
	population.clear();
	for (auto &i : newPop)
		population.append(i);
}

}

