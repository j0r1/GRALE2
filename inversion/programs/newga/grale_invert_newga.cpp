#include "log.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "lensinversiongafactoryparamssingleplanecpu.h"
#include "lensinversiongafactorycommon.h"
#include "inversioncommunicator.h"
#include "randomnumbergenerator.h"
#include <errut/booltype.h>
#include <serut/memoryserializer.h>
#include <mogal/geneticalgorithm.h>
#include <mogal2/geneticalgorithm.h>
#include <mogal2/randomnumbergenerator.h>
#include <mogal2/vectorgenomefitness.h>
#include <mogal2/singlethreadedpopulationfitnesscalculation.h>
#include <mogal2/stopcriterion.h>
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

class LensFitness : public mogal2::Fitness
{
public:
	LensFitness(size_t numObjectives)
	{
		m_fitnesses.resize(numObjectives);
		m_scaleFactor = numeric_limits<float>::quiet_NaN();
	}

	std::shared_ptr<Fitness> createCopy(bool copyContents = true) const
	{
		auto c = make_shared<LensFitness>(m_fitnesses.size());
		if (copyContents)
		{
			for (size_t i = 0 ; i < m_fitnesses.size() ; i++)
				c->m_fitnesses[i] = m_fitnesses[i];
			c->m_scaleFactor = m_scaleFactor;
		}
		return c;
	}

	std::string toString() const
	{
		if (!isCalculated())
			return "?";

		stringstream ss;
		ss << "[";
		for (auto x : m_fitnesses)
			ss << " " << x;
		
		ss << " * " << m_scaleFactor << " ]";
		return ss.str();
	}

	errut::bool_t MPI_BroadcastLayout(int root, MPI_Comm communicator)
	{
		int s = (int)m_fitnesses.size();
		MPI_Bcast(&s, 1, MPI_INT, root, communicator);
		m_fitnesses.resize(s);
		return true;
	}

	errut::bool_t MPI_Send(int dest, int tag, MPI_Comm communicator,
	                               std::vector<MPI_Request> &requests) const
	{
		requests.resize(2);
		MPI_Isend(m_fitnesses.data(), m_fitnesses.size(), MPI_FLOAT, dest, tag, communicator, &requests[0]);
		MPI_Isend(&m_scaleFactor, 1, MPI_FLOAT, dest, tag, communicator, &requests[1]);
		return true;
	}
	errut::bool_t MPI_Recv(int src, int tag, MPI_Comm communicator,
								   std::vector<MPI_Request> &requests)
	{
		requests.resize(2);
		MPI_Irecv(m_fitnesses.data(), m_fitnesses.size(), MPI_FLOAT, src, tag, communicator, &requests[0]);
		MPI_Irecv(&m_scaleFactor, 1, MPI_FLOAT, src, tag, communicator, &requests[1]);
		return true;
	}

	vector<float> m_fitnesses;
	float m_scaleFactor;
};

class LensGenome : public mogal2::Genome
{
public:
	LensGenome(size_t numWeights, size_t numSheetValues)
	{
		m_weights.resize(numWeights);
		m_sheets.resize(numSheetValues);
		m_scaleFactor = numeric_limits<float>::quiet_NaN(); // Will be returned by the fitness calculation
	}

	std::shared_ptr<Genome> createCopy(bool copyContents = true) const
	{
		auto c = make_shared<LensGenome>(m_weights.size(), m_sheets.size());
		if (copyContents)
		{
			for (size_t i = 0 ; i < m_weights.size() ; i++)
				c->m_weights[i] = m_weights[i];
			for (size_t i = 0 ; i < m_sheets.size() ; i++)
				c->m_sheets[i] = m_sheets[i];
			c->m_scaleFactor = m_scaleFactor;
		}
		return c;
	}
	std::string toString() const
	{
		stringstream ss;
		ss << "[";
		for (auto x : m_weights)
			ss << " " << x;
		
		if (m_sheets.size() > 0)
		{
			ss << " |";
			for (auto x : m_sheets)
				ss << " " << x;
		}
		ss << " * " << m_scaleFactor << " ]";
		return ss.str();
	}

	errut::bool_t MPI_BroadcastLayout(int root, MPI_Comm communicator)
	{
		int sizes[2] = { (int)m_weights.size(), (int)m_sheets.size() };
		MPI_Bcast(sizes, 2, MPI_INT, root, communicator);
		m_weights.resize(sizes[0]);
		m_sheets.resize(sizes[1]);
		return true;
	}
	
	errut::bool_t MPI_Send(int dest, int tag, MPI_Comm communicator,
	                               std::vector<MPI_Request> &requests) const
	{
		size_t n = (m_sheets.empty())?1:2;
		requests.resize(n);

		MPI_Isend(m_weights.data(), m_weights.size(), MPI_FLOAT, dest, tag, communicator, &requests[0]);
		if (!m_sheets.empty())
			MPI_Isend(m_sheets.data(), m_sheets.size(), MPI_FLOAT, dest, tag, communicator, &requests[1]);
		return true;
	}
	
	errut::bool_t MPI_Recv(int src, int tag, MPI_Comm communicator,
								   std::vector<MPI_Request> &requests)
	{
		size_t n = (m_sheets.empty())?1:2;
		requests.resize(n);

		MPI_Irecv(m_weights.data(), m_weights.size(), MPI_FLOAT, src, tag, communicator, &requests[0]);
		if (!m_sheets.empty())
			MPI_Irecv(m_sheets.data(), m_sheets.size(), MPI_FLOAT, src, tag, communicator, &requests[1]);
		return true;		
	}

	vector<float> m_weights;
	vector<float> m_sheets;
	float m_scaleFactor;
};

class Creation : public mogal2::IndividualCreation
{
public:
	Creation(grale::LensInversionGAFactoryCommon &factory) : m_pFactory(&factory) { }
	~Creation() { }

    std::shared_ptr<mogal2::Genome> createInitializedGenome() override
	{
		unique_ptr<grale::LensInversionGenome> g(static_cast<grale::LensInversionGenome*>(m_pFactory->createNewGenome()));

		auto &basisWeights = g->getBasisFunctionWeights();
		auto &sheetValues = g->getSheetValues();

		shared_ptr<LensGenome> genome = make_shared<LensGenome>(basisWeights.size(), sheetValues.size());
		genome->m_weights = basisWeights;
		genome->m_sheets = sheetValues;
		return genome;
	}

    std::shared_ptr<mogal2::Fitness> createEmptyFitness() override
	{
		return make_shared<LensFitness>(m_pFactory->getNumberOfFitnessComponents());
	}
private:
	grale::LensInversionGAFactoryCommon *m_pFactory;
};

// Wrap the GSL based RNG for now
class RNG : public mogal2::RandomNumberGenerator
{
public:
	RNG(const mogal::RandomNumberGenerator *pRng) : m_pRng(pRng) { }
	~RNG() { }
    double getRandomDouble() override { return m_pRng->pickRandomNumber(); }
    float getRandomFloat() override
	{
		cerr << "getRandomFloat NOT IMPLEMENTED" << endl;
		exit(-1);
		return 0;
	}
	uint32_t getRandomUint32() override
	{
		cerr << "getRandomUint32 NOT IMPLEMENTED" << endl;
		exit(-1);
		return 0;
	}
private:
	const mogal::RandomNumberGenerator *m_pRng;
};

class MyCrossover : public mogal2::PopulationCrossover
{
public:
	errut::bool_t check(const std::shared_ptr<mogal2::Population> &population)
	{
		if (population->size() == 0)
			return "Empty population";
		auto &i = population->individual(0);
		if (!dynamic_cast<const LensGenome *>(i->genomePtr()))
			return "Genome is of wrong type";
		if (!dynamic_cast<const LensFitness *>(i->fitnessPtr()))
			return "Fitness is of wrong type";
		return true;
	}

	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<mogal2::Population> &population, size_t targetPopulationSize)
	{
		return "TODO: MyCrossover";
	}
};

class LensFitnessCalculation : public mogal2::GenomeFitnessCalculation
{
public:
	LensFitnessCalculation(grale::LensInversionGAFactoryCommon &factory) : m_pFactory(&factory) { }
	~LensFitnessCalculation() { }

	errut::bool_t calculate(const mogal2::Genome &genome, mogal2::Fitness &fitness)
	{
		const LensGenome &g = dynamic_cast<const LensGenome &>(genome);
		LensFitness &f = dynamic_cast<LensFitness &>(fitness);
		
		// Anything that needs to get communicated back needs to be in the fitness, so
		// we'll store the scalefactor in the fitness and transfer it to the genome in a
		// later stage
		if (!m_pFactory->calculateFitness(g.m_weights, g.m_sheets, f.m_scaleFactor, f.m_fitnesses.data()))
			return "Error calculating fitness: " + m_pFactory->getErrorString();
		return true;
	}
private:
	grale::LensInversionGAFactoryCommon *m_pFactory;
};

class MyGA : public mogal2::GeneticAlgorithm
{
public:
	MyGA() { }
	~MyGA() { }

	errut::bool_t onBeforeFitnessCalculation(size_t generation, std::shared_ptr<mogal2::Population> &population)
	{
		cout << "# Generation " << generation << ", before calculation: " << endl;
		for (auto &i : population->individuals())
			cout << i->fitness()->toString() << endl;
		return true;
	}

    errut::bool_t onFitnessCalculated(size_t generation, std::shared_ptr<mogal2::Population> &population)
	{
		cout << "# Generation " << generation << ", calculated: " << endl;
		for (auto &i : population->individuals())
			cout << i->fitness()->toString() << endl;
		return true;
	}
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
		shared_ptr<RNG> rng = make_shared<RNG>(gaFactory.getRandomNumberGenerator());
		MyGA ga;

		Creation creation(gaFactory);
		MyCrossover cross;
		mogal2::SingleThreadedPopulationFitnessCalculation calc(make_shared<LensFitnessCalculation>(gaFactory));
		mogal2::FixedGenerationsStopCriterion stop(10); // for testing

		if (!(r = ga.run(creation, cross, calc, stop, popSize)))
			return "Error running GA: " + r.getErrorString();

		return true;
	}

//	const mogal::GeneticAlgorithm *getGeneticAlgorithm() const override { return &m_ga; }
//private:
//	CommonGA<GeneticAlgorithm> m_ga;
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
