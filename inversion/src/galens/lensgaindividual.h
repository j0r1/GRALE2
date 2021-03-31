#pragma once

#include "graleconfig.h"
#include <eatk/population.h>
#include <eatk/randomnumbergenerator.h>
#include <limits>

namespace grale
{

class LensGAFitness : public eatk::Fitness
{
public:
	LensGAFitness(size_t numObjectives);

	std::shared_ptr<Fitness> createCopy(bool copyContents = true) const override;
	std::string toString() const override;

	errut::bool_t MPI_BroadcastLayout(int root, MPI_Comm communicator) override;
	errut::bool_t MPI_Send(int dest, int tag, MPI_Comm communicator, std::vector<MPI_Request> &requests) const override;
	errut::bool_t MPI_Recv(int src, int tag, MPI_Comm communicator, std::vector<MPI_Request> &requests) override;

	std::vector<float> m_fitnesses;
	float m_scaleFactor;
};

class LensGAGenome : public eatk::Genome
{
public:
	LensGAGenome(size_t numWeights, size_t numSheetValues);

	std::shared_ptr<Genome> createCopy(bool copyContents = true) const override;
	std::string toString() const override;

	errut::bool_t MPI_BroadcastLayout(int root, MPI_Comm communicator) override;
	errut::bool_t MPI_Send(int dest, int tag, MPI_Comm communicator, std::vector<MPI_Request> &requests) const override;
	errut::bool_t MPI_Recv(int src, int tag, MPI_Comm communicator,
								   std::vector<MPI_Request> &requests) override;

	std::vector<float> m_weights;
	std::vector<float> m_sheets;
	float m_scaleFactor;
};

class LensGAIndividual : public eatk::Individual
{
public:
	LensGAIndividual(const std::shared_ptr<eatk::Genome> &genome, const std::shared_ptr<eatk::Fitness> &fitness,
			   size_t introducedInGeneration = std::numeric_limits<size_t>::max());

	std::shared_ptr<eatk::Individual> createNew(const std::shared_ptr<eatk::Genome> &genome, const std::shared_ptr<eatk::Fitness> &fitness,
			   size_t introducedInGeneration = std::numeric_limits<size_t>::max()) const override;

	int m_parent1, m_parent2;
	int m_ownIndex;
};

class LensGAIndividualCreation : public eatk::IndividualCreation
{
public:
	LensGAIndividualCreation(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
		     size_t numBasisFunctions, size_t numSheets, bool allowNegative,
			 size_t numObjectives);
	~LensGAIndividualCreation();

    std::shared_ptr<eatk::Genome> createUnInitializedGenome() override;
	std::shared_ptr<eatk::Genome> createInitializedGenome() override;
    std::shared_ptr<eatk::Fitness> createEmptyFitness() override;
	std::shared_ptr<eatk::Individual> createReferenceIndividual() override;
private:
	std::shared_ptr<eatk::RandomNumberGenerator> m_rng;
	const size_t m_numBasisFunctions, m_numSheets, m_numObjectives;
	const bool m_allowNegative;
};

}