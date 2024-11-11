#pragma once

#include "graleconfig.h"
#include <eatk/population.h>

namespace grale
{

class PopulationDump
{
public:
	PopulationDump(bool m_freeform = true);
	~PopulationDump();

	void checkDumpLoad(size_t generation, eatk::Population &population);
private:
	void dumpPopulation(const eatk::Population &population, const std::string &filename);
	void loadPopulation(eatk::Population &population, const std::string &filename);

	size_t m_dumpPopulationGeneration, m_loadPopulationGeneration;
	std::string m_dumpPopulationFilename, m_loadPopulationFilename;
	bool m_freeform;
};

}
