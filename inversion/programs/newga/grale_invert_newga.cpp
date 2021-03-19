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
#include <serut/memoryserializer.h>
#include <mogal/geneticalgorithm.h>
#include <mogal2/geneticalgorithm.h>
#include <mogal2/randomnumbergenerator.h>
#include <mogal2/singlethreadedpopulationfitnesscalculation.h>
#include <mogal2/stopcriterion.h>
#include <mogal2/simplesortedpopulation.h>
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

class Creation : public mogal2::IndividualCreation
{
public:
	Creation(const std::shared_ptr<mogal2::RandomNumberGenerator> &rng,
		     size_t numBasisFunctions, size_t numSheets, bool allowNegative,
			 size_t numObjectives)
	 : m_rng(rng), m_numBasisFunctions(numBasisFunctions), m_numSheets(numSheets), m_allowNegative(allowNegative),
	   m_numObjectives(numObjectives)
	{
	}

	~Creation() { }

    std::shared_ptr<mogal2::Genome> createInitializedGenome() override
	{
		shared_ptr<grale::LensGAGenome> genome = make_shared<grale::LensGAGenome>(m_numBasisFunctions, m_numSheets);

		for (auto &x : genome->m_weights)
			x = (float)m_rng->getRandomDouble();

		if (m_allowNegative)
		{
			for (auto &x : genome->m_weights)
				x -= 0.5f;
		}

		for (auto &s : genome->m_sheets)
			s = (float)m_rng->getRandomDouble();

		return genome;
	}

    std::shared_ptr<mogal2::Fitness> createEmptyFitness() override
	{
		return make_shared<grale::LensGAFitness>(m_numObjectives);
	}

	std::shared_ptr<mogal2::Individual> createReferenceIndividual() override
	{
		return std::make_shared<grale::LensGAIndividual>(nullptr, nullptr);
	}
private:
	std::shared_ptr<mogal2::RandomNumberGenerator> m_rng;
	const size_t m_numBasisFunctions, m_numSheets, m_numObjectives;
	const bool m_allowNegative;
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

// TODO: this is for a single fitness measure
class MyCrossover : public mogal2::PopulationCrossover
{
public:
	MyCrossover(double beta, bool elitism, bool includeBest, double crossoverRate,
				const shared_ptr<mogal2::RandomNumberGenerator> &rng,
				bool allowNegative,
				const shared_ptr<mogal2::GenomeMutation> &mutation)
		: m_beta(beta), m_bestWithoutMutation(elitism), m_bestWithMutation(includeBest), m_crossoverRate(crossoverRate),
		  m_rng(rng), m_sortedPop(make_shared<grale::LensGAFitnessComparison>()), m_cross(rng, allowNegative),
		  m_mutation(mutation)
	{

	}

	errut::bool_t check(const std::shared_ptr<mogal2::Population> &population)
	{
		if (population->size() == 0)
			return "Empty population";
		auto &i = population->individual(0);
		if (!dynamic_cast<const grale::LensGAGenome *>(i->genomePtr()))
			return "Genome is of wrong type";
		if (!dynamic_cast<const grale::LensGAFitness *>(i->fitnessPtr()))
			return "Fitness is of wrong type";

		bool_t r;
		vector<shared_ptr<mogal2::Genome>> testParents = { population->individual(0)->genome(), population->individual(0)->genome() };
		if (!(r = m_cross.check(testParents)))
			return "Error in genome crossover check: " + r.getErrorString();

		if (!(r = m_mutation->check(population->individual(0)->genomeRef())))
			return "Error checking mutation: " + r.getErrorString();

		return true;
	}

	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<mogal2::Population> &population, size_t targetPopulationSize)
	{
		bool_t r;
		if (generation == 0)
		{
			if (!(r = m_sortedPop.check(*population)))
				return "Error checking sorter: " + r.getErrorString();
		}

		// TODO: do this in another function?
		// Scale factor is calculated and stored in fitness, copy it back to genome
		for (auto &i : population->individuals())
		{
			auto &g = static_cast<grale::LensGAGenome &>(i->genomeRef());
			auto &f = static_cast<grale::LensGAFitness &>(i->fitnessRef());
			g.m_scaleFactor = f.m_scaleFactor;
		}

		if (population->size() != targetPopulationSize)
			return "Expecting size " + to_string(targetPopulationSize) + " but got " + to_string(population->size());
		
		if (!(r = m_sortedPop.processPopulation(population, targetPopulationSize)))
			return "Error sorting population: " + r.getErrorString();

		shared_ptr<mogal2::Population> newPop = make_shared<mogal2::Population>();

		auto appendBest = [&newPop, &population]()
		{
			auto ind = population->individual(0)->createCopy();
			grale::LensGAIndividual &i = static_cast<grale::LensGAIndividual&>(*ind);
			i.m_parent1 = 0;
			i.m_parent2 = -1;
			newPop->append(ind);
		};

		size_t mutOffset = 0;
		if (m_bestWithoutMutation)
		{
			mutOffset = 1;
			appendBest();
		}
		if (m_bestWithMutation)
			appendBest();
		
		auto pickParent = [&population, this]()
		{
			double x = m_rng->getRandomDouble();
			double val = (1.0-std::pow(x, 1.0/(1.0+m_beta)))*((double)population->size());
			int r = (int)val;
			if (r >= (int)population->size())
					r = (int)population->size() - 1;
			return r;
		};

		auto getLensIndividual = [&population](size_t i) -> const grale::LensGAIndividual &
		{
			return static_cast<grale::LensGAIndividual&>(*population->individual(i));
		};

		// In original version, an attempt was made to prevent inbreeding
		auto pickParentIndices = [&pickParent, &getLensIndividual](int &index1, int &index2)
		{
			bool ok;
			int count = 0;

			do 
			{
				ok = false;

				index1 = pickParent();
				index2 = pickParent();

				// prevent inbreeding

				int a1 = getLensIndividual(index1).m_parent1;
				int a2 = getLensIndividual(index1).m_parent2;
				int b1 = getLensIndividual(index2).m_parent1;
				int b2 = getLensIndividual(index2).m_parent2;

				if (a1 < 0 || b1 < 0) // one of them is a brand new genome
					ok = true;
				else
				{
					if (!(a1 == b1 || a1 == b2 || b1 == a2 || ((a2 >= 0) && a2 == b2)))
						ok = true;
				}
				count++;
			} while (count < 10 && !ok);
		};

		vector<shared_ptr<mogal2::Genome>> parents(2);
		vector<shared_ptr<mogal2::Genome>> offspring;
		while (newPop->size() < targetPopulationSize)
		{
			if (m_rng->getRandomDouble() < m_crossoverRate)
			{
				int index1 = -1, index2 = -1;
				pickParentIndices(index1, index2);
				parents[0] = population->individual(index1)->genome();
				parents[1] = population->individual(index2)->genome();
				
				if (!(r = m_cross.generateOffspring(parents, offspring)))
					return "Error in crossover: " + r.getErrorString();

				for (auto &g : offspring)
				{
					auto ind = make_shared<grale::LensGAIndividual>(g, population->individual(0)->fitness()->createCopy(false), generation);
					ind->m_parent1 = index1;
					ind->m_parent2 = index2;
					newPop->append(ind);
				}
			}
			else // clone
			{
				int index = pickParent();
				auto ind = population->individual(index)->createCopy();

				grale::LensGAIndividual &lensInd = static_cast<grale::LensGAIndividual&>(*ind);
				lensInd.m_parent1 = index;
				lensInd.m_parent2 = -1;
				newPop->append(ind);
			}
		}

		for (size_t i = mutOffset ; i < newPop->size() ; i++)
		{
			bool isChanged = false;

			auto &ind = newPop->individual(i);

			if (!(r = m_mutation->mutate(ind->genomeRef(), isChanged)))
				return "Error mutating genome: " + r.getErrorString();
			if (isChanged)
				ind->fitness()->setCalculated(false);
		}
		
		swap(population, newPop);
		
		return true;
	}

	shared_ptr<mogal2::RandomNumberGenerator> m_rng;
	mogal2::SimpleSortedPopulation m_sortedPop;
	grale::LensGAGenomeCrossover m_cross;
	double m_beta, m_bestWithMutation, m_bestWithoutMutation, m_crossoverRate;
	shared_ptr<mogal2::GenomeMutation> m_mutation;
};

class LensFitnessCalculation : public mogal2::GenomeFitnessCalculation
{
public:
	LensFitnessCalculation(grale::LensInversionGAFactoryCommon &factory) : m_pFactory(&factory) { }
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
		if (gaFactory.getNumberOfFitnessComponents() != 1)
			return "Currently only for one fitness component";

		shared_ptr<RNG> rng = make_shared<RNG>(gaFactory.getRandomNumberGenerator());
		MyGA ga;

		auto mutation = make_shared<grale::LensGAGenomeMutation>(rng, 
					   gaFactory.getChanceMultiplier(),
					   gaFactory.allowNegativeValues(),
					   gaFactory.getMutationAmplitude(),
					   gaFactory.useAbsoluteMutation());

		Creation creation(rng, gaFactory.getNumberOfBasisFunctions(), gaFactory.getNumberOfSheets(),
		                  gaFactory.allowNegativeValues(), gaFactory.getNumberOfFitnessComponents());

		MyCrossover cross(params.getBeta(),
					      params.useElitism(),
						  params.alwaysIncludeBestGenome(),
						  params.getCrossOverRate(),
						  rng,
						  gaFactory.allowNegativeValues(),
						  mutation);
		mogal2::SingleThreadedPopulationFitnessCalculation calc(make_shared<LensFitnessCalculation>(gaFactory));
		mogal2::FixedGenerationsStopCriterion stop(11); // for testing

		if (!(r = ga.run(creation, cross, calc, stop, popSize)))
			return "Error running GA: " + r.getErrorString();

		auto &bestSolutions = cross.m_sortedPop.getBestIndividuals();
		cout << "Best: " << endl;
		for (auto &b: bestSolutions)
			cout << b->fitness()->toString() << endl;

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
