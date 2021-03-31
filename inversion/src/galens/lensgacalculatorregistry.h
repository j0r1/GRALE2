#pragma once

#include "graleconfig.h"
#include "lensinversionparametersbase.h"
#include "lensgagenomecalculator.h"
#include "lensfitnessobject.h"
#include <errut/booltype.h>
#include <memory>
#include <map>

namespace grale
{

class LensGACalculatorFactory
{
public:
	LensGACalculatorFactory() { }
	~LensGACalculatorFactory() { }

	virtual std::unique_ptr<LensInversionParametersBase> createParametersInstance() = 0;
	virtual std::unique_ptr<LensGAGenomeCalculator> createCalculatorInstance(std::unique_ptr<grale::LensFitnessObject> fitObj) = 0;
};

class LensGACalculatorRegistry
{
public:
	static LensGACalculatorRegistry &instance();

	~LensGACalculatorRegistry();

	errut::bool_t registerCalculatorFactory(const std::string &name, std::unique_ptr<LensGACalculatorFactory> factory);
	LensGACalculatorFactory *getFactory(const std::string &name);
private:
	LensGACalculatorRegistry();

	std::map<std::string, std::unique_ptr<LensGACalculatorFactory>> m_registry;
	static std::unique_ptr<LensGACalculatorRegistry> s_instance;
};

}