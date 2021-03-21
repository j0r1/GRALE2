#include "log.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "lensinversiongafactoryparamssingleplanecpu.h"
#include "lensinversiongafactorycommon.h"
#include "inversioncommunicator.h"
#include "randomnumbergenerator.h"
#include "constants.h"
#include "lensgaindividual.h"
#include "lensgagenomemutation.h"
#include "lensgagenomecrossover.h"
#include "lensgafitnesscomparison.h"
#include "lensgasingleobjectivecrossover.h"
#include "galensmodule.h"
#include "lensgastopcriterion.h"
#include <serut/memoryserializer.h>
#include <mogal/geneticalgorithm.h>
#include <mogal2/geneticalgorithm.h>
#include <mogal2/randomnumbergenerator.h>
#include <mogal2/singlethreadedpopulationfitnesscalculation.h>
#include <mogal2/multithreadedpopulationfitnesscalculation.h>
#include <mogal2/stopcriterion.h>
#include <mogal2/simplesortedpopulation.h>
#include <serut/vectorserializer.h>
#ifndef WIN32
#include <fcntl.h>
#endif // !WIN32

#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <limits>

using namespace std;
using namespace serut;
using namespace errut;

// Wrap the GSL based RNG for now
class RNG : public mogal2::RandomNumberGenerator
{
public:
	RNG(const mogal::RandomNumberGenerator *pRng) : m_pRng(pRng) { }
	~RNG() { }
    double getRandomDouble() override { return m_pRng->pickRandomNumber(); }
    float getRandomFloat() override { return (float)m_pRng->pickRandomNumber(); }
	uint32_t getRandomUint32() override
	{
		cerr << "getRandomUint32 NOT IMPLEMENTED" << endl;
		exit(-1);
		return 0;
	}
private:
	const mogal::RandomNumberGenerator *m_pRng;
};


class LensFitnessCalculation : public mogal2::GenomeFitnessCalculation
{
public:
	LensFitnessCalculation(grale::LensInversionGAFactoryCommon &factory) : m_pFactory(&factory) { }
	LensFitnessCalculation(const shared_ptr<grale::LensInversionGAFactoryCommon> &factory)
	{
		m_spFactory = factory; // keep a reference to the lib
		m_pFactory = factory.get();
	}

	~LensFitnessCalculation() { }

	errut::bool_t calculate(const mogal2::Genome &genome, mogal2::Fitness &fitness)
	{
		const grale::LensGAGenome &g = dynamic_cast<const grale::LensGAGenome &>(genome);
		grale::LensGAFitness &f = dynamic_cast<grale::LensGAFitness &>(fitness);
		
		// Anything that needs to get communicated back needs to be in the fitness, so
		// we'll store the scalefactor in the fitness and transfer it to the genome in a
		// later stage
		if (!m_pFactory->calculateFitness(g.m_weights, g.m_sheets, f.m_scaleFactor, f.m_fitnesses.data()))
			return "Error calculating fitness: " + m_pFactory->getErrorString();
		return true;
	}
private:
	grale::LensInversionGAFactoryCommon *m_pFactory;
	shared_ptr<grale::LensInversionGAFactoryCommon> m_spFactory;
};

class MyGA : public mogal2::GeneticAlgorithm
{
public:
	MyGA() { }
	~MyGA() { }

	// errut::bool_t onBeforeFitnessCalculation(size_t generation, std::shared_ptr<mogal2::Population> &population)
	// {
	// 	cout << "# Generation " << generation << ", before calculation: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
	// 	return true;
	// }

    // errut::bool_t onFitnessCalculated(size_t generation, std::shared_ptr<mogal2::Population> &population)
	// {
	// 	cout << "# Generation " << generation << ", calculated: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
	// 	return true;
	// }
};

// TODO: rename this, is from copy-paste
class SingleCoreCommunicator : public InversionCommunicator
{
public:
	SingleCoreCommunicator() { }
	~SingleCoreCommunicator() { }
protected:
	string getVersionInfo() const override { return "Single core genetic algorithm engine"; }

	bool_t runGA(int popSize, mogal::GAFactory &factory, mogal::GeneticAlgorithmParams &params,
	             const std::string &moduleDir, const std::string &moduleFile,
	             const std::vector<uint8_t> &factoryParamBytes) override
	{
		bool_t r;
		grale::LensInversionGAFactoryCommon &gaFactory = dynamic_cast<grale::LensInversionGAFactoryCommon&>(factory);
		if (gaFactory.getNumberOfFitnessComponents() != 1)
			return "Currently only for one fitness component";

		shared_ptr<RNG> rng = make_shared<RNG>(gaFactory.getRandomNumberGenerator());
		MyGA ga;

		auto mutation = make_shared<grale::LensGAGenomeMutation>(rng, 
					   gaFactory.getChanceMultiplier(),
					   gaFactory.allowNegativeValues(),
					   gaFactory.getMutationAmplitude(),
					   gaFactory.useAbsoluteMutation());

		grale::LensGAIndividualCreation creation(rng, gaFactory.getNumberOfBasisFunctions(), gaFactory.getNumberOfSheets(),
		                  gaFactory.allowNegativeValues(), gaFactory.getNumberOfFitnessComponents());

		grale::LensGASingleFitnessCrossover cross(params.getBeta(),
					      params.useElitism(),
						  params.alwaysIncludeBestGenome(),
						  params.getCrossOverRate(),
						  rng,
						  gaFactory.allowNegativeValues(),
						  mutation);

		size_t numThreads = 1;
		if (getenv("NUMTHREADS"))
			numThreads = stoi(getenv("NUMTHREADS"));
		cerr << "Using " << numThreads << " threads " << endl;

		vector<shared_ptr<mogal2::GenomeFitnessCalculation>> genomeFitnessCalculators = { make_shared<LensFitnessCalculation>(gaFactory) };
		grale::GALensModule module;

		if (!module.open(moduleDir, moduleFile))
			return "Couldn't open module: " + module.getErrorString();

		for (size_t i = 2 ; i < numThreads ; i++)
		{
			auto newFac = shared_ptr<grale::LensInversionGAFactoryCommon>(dynamic_cast<grale::LensInversionGAFactoryCommon*>(module.createFactoryInstance()));
			auto facParams = shared_ptr<mogal::GAFactoryParams>(newFac->createParamsInstance());

			serut::VectorSerializer ser(factoryParamBytes);
			if (!facParams->read(ser))
				return "Coudln't deserialize factory parameters: " + facParams->getErrorString();

			if (!newFac->init(facParams.get()))
				return "Unable to init factory copy: " + newFac->getErrorString();
			
			genomeFitnessCalculators.push_back(make_shared<LensFitnessCalculation>(newFac));
		}
		
		shared_ptr<mogal2::PopulationFitnessCalculation> calc;
		if (numThreads <= 1)
			calc = make_shared<mogal2::SingleThreadedPopulationFitnessCalculation>(genomeFitnessCalculators[0]);
		else
		{
			auto mpCalc = make_shared<mogal2::MultiThreadedPopulationFitnessCalculation>();
			if (!(r = mpCalc->initThreadPool(genomeFitnessCalculators)))
				return "Unable to initialize threads: " + r.getErrorString();
			calc = mpCalc;
		}
		
		grale::LensGAStropCriterion stop(gaFactory.getMaximumNumberOfGenerations(), mutation);
		if (!(r = stop.initialize(gaFactory.getFitnessObject())))
			return "Error initializing convergence checker: " + r.getErrorString();

		if (!(r = ga.run(creation, cross, *calc, stop, popSize)))
			return "Error running GA: " + r.getErrorString();

		auto &bestSolutions = cross.getBestIndividuals();
		cout << "Best: " << endl;
		for (auto &b: bestSolutions)
			cout << b->fitness()->toString() << endl;

		return true;
	}
};

int main(int argc, char *argv[])
{
	grale::LOG.init(argv[0]);
#ifndef WIN32
	// This was added to be able to launch instance in gdb, in which case
	// stdin/stdout can no longer be used for communication with the
	// executable (needed to control gdb).
	if (argc == 3)
	{
		string inputPipeName(argv[1]);
		string outputPipeName(argv[2]);
		cerr << "Opening " << outputPipeName << " for communication" << endl;
		// Nasty: override the file desc for communication
		stdOutFileDescriptor = open(outputPipeName.c_str(), O_WRONLY); 
		cerr << "Opening " << inputPipeName << " for communication" << endl;
		// Nasty: override the file desc for communication
		stdInFileDescriptor = open(inputPipeName.c_str(), O_RDWR); 
	}
#endif // WIN32

	SingleCoreCommunicator comm;

	bool_t r = comm.run();
	if (!r)
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}

	return 0;
}
