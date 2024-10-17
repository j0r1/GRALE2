#pragma once

#include "graleconfig.h"
#include "lensinversionparametersbase.h"
#include "lensgaindividual.h"
#include "gravitationallens.h"
#include "lensfitnessobject.h"
#include <eatk/genomefitness.h>
#include <eatk/calculation.h>

#include <iostream>

namespace grale
{

class LensGAGenomeCalculatorLogger
{
public:
	virtual void log(const std::string &s) const = 0;
};

// TODO: perhaps this should be renamed?
class LensGAGenomeCalculator : public eatk::GenomeFitnessCalculation
{
public:
	void setLogger(const std::shared_ptr<LensGAGenomeCalculatorLogger> &lgr) { m_logger = lgr; }

	virtual errut::bool_t init(const LensInversionParametersBase &params) = 0;
	virtual errut::bool_t createLens(const eatk::Genome &genome, std::unique_ptr<GravitationalLens> &lens) const = 0;

	virtual size_t getNumberOfObjectives() const = 0;
	// virtual bool allowNegativeValues() const = 0;
	// virtual size_t getNumberOfBasisFunctions() const = 0;
	// virtual size_t getNumberOfSheets() const = 0;

	void log(const std::string &s) const
	{
		if (m_logger.get()) 
			m_logger->log(s);
		else
			std::cerr << "NO_LOGGER_SET:" << s << std::endl;
	}
private:
	std::shared_ptr<LensGAGenomeCalculatorLogger> m_logger;
};

}
