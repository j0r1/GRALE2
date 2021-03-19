#pragma once

#include "graleconfig.h"
#include <mogal2/population.h>
#include <limits>

namespace grale
{

class LensGAFitness : public mogal2::Fitness
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

class LensGAGenome : public mogal2::Genome
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

class LensGAIndividual : public mogal2::Individual
{
public:
	LensGAIndividual(std::shared_ptr<mogal2::Genome> genome, std::shared_ptr<mogal2::Fitness> fitness,
			   size_t introducedInGeneration = std::numeric_limits<size_t>::max());

	std::shared_ptr<mogal2::Individual> createNew(std::shared_ptr<mogal2::Genome> genome, std::shared_ptr<mogal2::Fitness> fitness,
			   size_t introducedInGeneration = std::numeric_limits<size_t>::max()) const override;

	int m_parent1, m_parent2;
};

}