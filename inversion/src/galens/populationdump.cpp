#include "populationdump.h"
#include "lensgaindividual.h"
#include "utils.h"
#include <serut/fileserializer.h>
#include <memory>
#include <iostream>

using namespace errut;
using namespace std;

namespace grale
{

PopulationDump::PopulationDump()
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

	for (auto &ind : population.individuals())
	{
		auto pInd = dynamic_cast<const LensGAIndividual *>(ind.get());
		if (!pInd)
		{
			cerr << "WARNING: individual is of incorrect type, can't dump to file" << endl;
			return;
		}
		pInd->write(fSer);
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

	// swap pops
	population.clear();
	for (auto &i : newPop)
		population.append(i);
}

}

