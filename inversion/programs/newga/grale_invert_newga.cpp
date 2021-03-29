#include "newgacommunicatorbase.h"
#include "lensgacalculatorregistry.h"
#include "lensinversiongafactorysingleplanecpu.h"
#include "lensinversiongafactoryparamssingleplanecpu.h"
#include "lensfitnessobject.h"
#include "lensfitnessgeneral.h"
#include "lensgacalculatorregistry.h"
#include <serut/vectorserializer.h>
#include <mogal/gafactorymultiobjective.h>
#ifndef WIN32
#include <fcntl.h>
#include <unistd.h>
#endif // !WIN32
#include <memory>

using namespace std;
using namespace serut;
using namespace errut;

class NewGACommunicatorThreads : public NewGACommunicatorBase
{
public:
	NewGACommunicatorThreads(size_t numThreads) : m_numThreads(numThreads)
	{
		cerr << "Using " << numThreads << " threads " << endl;
	}

	~NewGACommunicatorThreads() { }
protected:
	string getVersionInfo() const override { return "MOGAL2 Thread based algorithm, " + to_string(m_numThreads) + " threads"; }

	bool_t getCalculator(const std::string &lensFitnessObjectType, 
									grale::LensGACalculatorFactory &calcFactory, 
									const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
									const std::vector<uint8_t> &factoryParamBytes,
									grale::LensGAIndividualCreation &creation,
									std::shared_ptr<mogal2::PopulationFitnessCalculation> &calc) override
	{
		bool_t r;
		vector<shared_ptr<mogal2::GenomeFitnessCalculation>> genomeFitnessCalculators = { genomeCalculator };

		for (size_t i = 2 ; i < m_numThreads ; i++)
		{
			unique_ptr<grale::LensFitnessObject> fitObj = grale::LensFitnessObjectRegistry::instance().createFitnessObject(lensFitnessObjectType);
			if (!fitObj.get())
				return "No fitness object with name '" + lensFitnessObjectType + "' is known";

			auto calculatorParams = calcFactory.createParametersInstance();
			if (!(r = InversionCommunicator::loadFromBytes(*calculatorParams, factoryParamBytes)))
				return "Can't load calculator parameters from received data: " + r.getErrorString();
		
	 		shared_ptr<grale::LensGAGenomeCalculator> calculatorInstance = calcFactory.createCalculatorInstance(move(fitObj));
			if (!(r = calculatorInstance->init(*calculatorParams)))
			 	return "Unable to initialize calculator: " + r.getErrorString();

			genomeFitnessCalculators.push_back(calculatorInstance);
		}
		
		if (m_numThreads <= 1)
			calc = make_shared<mogal2::SingleThreadedPopulationFitnessCalculation>(genomeFitnessCalculators[0]);
		else
		{
			auto mpCalc = make_shared<mogal2::MultiThreadedPopulationFitnessCalculation>();
			if (!(r = mpCalc->initThreadPool(genomeFitnessCalculators)))
				return "Unable to initialize threads: " + r.getErrorString();
			calc = mpCalc;
		}
		return true;
	}
private:
	size_t m_numThreads;
};

class GAFactoryHelper : public grale::LensInversionGAFactorySinglePlaneCPU, public mogal::GAFactoryMultiObjective
{
public:
	GAFactoryHelper(unique_ptr<grale::LensFitnessObject> fitObj)
		: m_fitObj(move(fitObj)) { }
	GAFactoryHelper() { }

	// Fr now, we're not going to create new one, just use the previously set one
	grale::LensFitnessObject *createFitnessObject() override
	{
		cout << "GAFactoryHelpen: createFitnessObject " << (void*)m_fitObj.get() << endl;
		if (!m_fitObj.get())
		{
			setErrorString("LensFitnessObject was already retrieved");
			return nullptr;
		}
		auto obj = m_fitObj.release(); // The caller will take care of this
		return obj;
	}

	bool subInit(grale::LensFitnessObject *pFitnessObject)
	{
		int numFitness = pFitnessObject->getNumberOfFitnessComponents();
		setNumberOfFitnessComponents(numFitness);
		return true;
	}

	void onGeneticAlgorithmStep(int generation,  
			                    bool *generationInfoChanged, bool *stopAlgorithm)
	{
		sendMessage("Should not be called");
		*stopAlgorithm = true;
	}

	float getChanceMultiplier() { return 1.0f; }
	bool useAbsoluteMutation() { return true; }
	float getMutationAmplitude() { return 0.0f; }

	mogal::Genome *selectPreferredGenome(const std::list<mogal::Genome *> &nonDominatedSet) const
	{
		return nullptr;
	}
private:
	unique_ptr<grale::LensFitnessObject> m_fitObj;
};

class GAFactoryWrapperLensGAGenomeCalculator: public grale::LensGAGenomeCalculator
{
public:
	GAFactoryWrapperLensGAGenomeCalculator(unique_ptr<GAFactoryHelper> helper)
		: m_helperFactory(move(helper)) { }

	bool_t init(const grale::LensInversionParametersBase &params) override
	{
		// Convert the parameters to the right class
		serut::VectorSerializer ser;
		if (!params.write(ser))
			return "Error serializing parameters: " + params.getErrorString();
		
		grale::LensInversionGAFactoryParamsSinglePlaneCPU gaParams;
		if (!gaParams.read(ser))
			return "Error re-reading parameters: " + gaParams.getErrorString();

		if (!m_helperFactory->init(&gaParams))
			return "Can't init helper factory: " + m_helperFactory->getErrorString();

		return true;
	}

	bool_t createLens(const grale::LensGAGenome &genome, std::unique_ptr<grale::GravitationalLens> &lens) const override
	{
		string errStr = "unknown error";

		lens = move(unique_ptr<grale::GravitationalLens>(m_helperFactory->createLens(genome.m_weights, genome.m_sheets, genome.m_scaleFactor, errStr)));
		if (!lens.get())
			return errStr;
		return true;
	}

	size_t getNumberOfObjectives() const override { return m_helperFactory->getNumberOfFitnessComponents(); }
	bool allowNegativeValues() const override { return m_helperFactory->allowNegativeValues(); }
	size_t getNumberOfBasisFunctions() const override { return m_helperFactory->getNumberOfBasisFunctions(); }
    size_t getNumberOfSheets() const override { return m_helperFactory->getNumberOfSheets(); }
	size_t getMaximumNumberOfGenerations() const override { return m_helperFactory->getMaximumNumberOfGenerations(); }
	const grale::LensFitnessObject &getLensFitnessObject() const override { return m_helperFactory->getFitnessObject(); }

	bool_t calculate(const mogal2::Genome &genome, mogal2::Fitness &fitness) override
	{
		const grale::LensGAGenome &g = static_cast<const grale::LensGAGenome&>(genome);
		grale::LensGAFitness &f = static_cast<grale::LensGAFitness&>(fitness);
		if (!m_helperFactory->calculateFitness(g.m_weights, g.m_sheets, f.m_scaleFactor, f.m_fitnesses.data()))
			return "Unable to calculate fitness: " + m_helperFactory->getErrorString();
		return true;		
	}
private:
	unique_ptr<GAFactoryHelper> m_helperFactory;
};

class GeneralFactory : public grale::LensGACalculatorFactory
{
public:
    std::unique_ptr<grale::LensInversionParametersBase> createParametersInstance() override
	{
		return make_unique<grale::LensInversionParametersSinglePlaneCPU>();
	}

    std::unique_ptr<grale::LensGAGenomeCalculator> createCalculatorInstance(unique_ptr<grale::LensFitnessObject> fitObj) override
	{
		auto helper = make_unique<GAFactoryHelper>(move(fitObj));
		auto wrapper = make_unique<GAFactoryWrapperLensGAGenomeCalculator>(move(helper));
		return wrapper;
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

	size_t numThreads = 1;
	if (getenv("NUMTHREADS"))
		numThreads = stoi(getenv("NUMTHREADS"));

	cout << "Registering name..." << endl;
	grale::LensGACalculatorRegistry::instance().registerCalculatorFactory("cpu",
		make_unique<GeneralFactory>());

	grale::LensFitnessObjectRegistry::instance().registerLensFitnessObject<grale::LensFitnessGeneral>("general");

	NewGACommunicatorThreads comm(numThreads);

	bool_t r = comm.run();
	if (!r)
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}

	return 0;
}
