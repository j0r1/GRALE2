#include "lensgaindividual.h"
#include <sstream>

using namespace eatk;
using namespace std;
using namespace errut;

namespace grale
{

LensGAFitness::LensGAFitness(size_t numObjectives)
{
	m_fitnesses.resize(numObjectives);
	m_scaleFactor = numeric_limits<float>::quiet_NaN();
}

shared_ptr<Fitness> LensGAFitness::createCopy(bool copyContents) const
{
	auto c = make_shared<LensGAFitness>(m_fitnesses.size());
	if (copyContents)
	{
		for (size_t i = 0 ; i < m_fitnesses.size() ; i++)
			c->m_fitnesses[i] = m_fitnesses[i];
		c->m_scaleFactor = m_scaleFactor;
		if (isCalculated())
			c->setCalculated();
	}
	return c;
}

string LensGAFitness::toString() const
{
	if (!isCalculated())
		return "?";

	stringstream ss;
	// ss << "[ ";
	for (auto x : m_fitnesses)
		ss << x << " ";
	
	// ss << "* " << m_scaleFactor << " ]";
	return ss.str();
}

bool_t LensGAFitness::read(serut::SerializationInterface &si)
{
	int32_t num;
	if (!si.readInt32(&num))
		return si.getErrorString();

	m_fitnesses.resize(num);
	if (!si.readFloats(m_fitnesses) || !si.readFloat(&m_scaleFactor))
		return si.getErrorString();
	return true;
}

bool_t LensGAFitness::write(serut::SerializationInterface &si) const
{
	int32_t num = (int32_t)m_fitnesses.size();
	if (!si.writeInt32(num) || !si.writeFloats(m_fitnesses) || !si.writeFloat(m_scaleFactor))
		return si.getErrorString();
	return true;
}

#ifdef EATKCONFIG_MPISUPPORT
bool_t LensGAFitness::MPI_BroadcastLayout(int root, MPI_Comm communicator)
{
	int s = (int)m_fitnesses.size();
	MPI_Bcast(&s, 1, MPI_INT, root, communicator);
	m_fitnesses.resize(s);
	return true;
}

bool_t LensGAFitness::MPI_Send(int dest, int tag, MPI_Comm communicator,
							   vector<MPI_Request> &requests) const
{
	requests.resize(2);
	MPI_Isend(m_fitnesses.data(), m_fitnesses.size(), MPI_FLOAT, dest, tag, communicator, &requests[0]);
	MPI_Isend(&m_scaleFactor, 1, MPI_FLOAT, dest, tag, communicator, &requests[1]);
	return true;
}

bool_t LensGAFitness::MPI_Recv(int src, int tag, MPI_Comm communicator,
							   vector<MPI_Request> &requests)
{
	requests.resize(2);
	MPI_Irecv(m_fitnesses.data(), m_fitnesses.size(), MPI_FLOAT, src, tag, communicator, &requests[0]);
	MPI_Irecv(&m_scaleFactor, 1, MPI_FLOAT, src, tag, communicator, &requests[1]);
	return true;
}
#endif // EATKCONFIG_MPISUPPORT

LensGAGenome::LensGAGenome(size_t numWeights, size_t numSheetValues)
{
	m_weights.resize(numWeights);
	m_sheets.resize(numSheetValues);
	m_scaleFactor = numeric_limits<float>::quiet_NaN(); // Will be returned by the fitness calculation
}

shared_ptr<eatk::Genome> LensGAGenome::createCopy(bool copyContents) const
{
	auto c = make_shared<LensGAGenome>(m_weights.size(), m_sheets.size());
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

string LensGAGenome::toString() const
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

bool_t LensGAGenome::read(serut::SerializationInterface &si)
{
	int32_t nums[2];
	if (!si.readInt32s(nums, 2))
		return si.getErrorString();

	m_weights.resize(nums[0]);
	m_sheets.resize(nums[1]);
	if (!si.readFloats(m_weights) || !si.readFloats(m_sheets) || !si.readFloat(&m_scaleFactor))
		return si.getErrorString();
	return true;
}

bool_t LensGAGenome::write(serut::SerializationInterface &si) const
{
	int32_t nums[2] = { (int32_t)m_weights.size(), (int32_t)m_sheets.size() };
	if (!si.writeInt32s(nums, 2) || !si.writeFloats(m_weights) || !si.writeFloats(m_sheets) || !si.writeFloat(m_scaleFactor))
		return si.getErrorString();
	return true;
}

#ifdef EATKCONFIG_MPISUPPORT
bool_t LensGAGenome::MPI_BroadcastLayout(int root, MPI_Comm communicator)
{
	int sizes[2] = { (int)m_weights.size(), (int)m_sheets.size() };
	MPI_Bcast(sizes, 2, MPI_INT, root, communicator);
	m_weights.resize(sizes[0]);
	m_sheets.resize(sizes[1]);
	return true;
}

bool_t LensGAGenome::MPI_Send(int dest, int tag, MPI_Comm communicator,
							  vector<MPI_Request> &requests) const
{
	size_t n = (m_sheets.empty())?1:2;
	requests.resize(n);

	MPI_Isend(m_weights.data(), m_weights.size(), MPI_FLOAT, dest, tag, communicator, &requests[0]);
	if (!m_sheets.empty())
		MPI_Isend(m_sheets.data(), m_sheets.size(), MPI_FLOAT, dest, tag, communicator, &requests[1]);
	return true;
}

bool_t LensGAGenome::MPI_Recv(int src, int tag, MPI_Comm communicator,
							  vector<MPI_Request> &requests)
{
	size_t n = (m_sheets.empty())?1:2;
	requests.resize(n);

	MPI_Irecv(m_weights.data(), m_weights.size(), MPI_FLOAT, src, tag, communicator, &requests[0]);
	if (!m_sheets.empty())
		MPI_Irecv(m_sheets.data(), m_sheets.size(), MPI_FLOAT, src, tag, communicator, &requests[1]);
	return true;		
}
#endif //EATKCONFIG_MPISUPPORT

LensGAIndividual::LensGAIndividual(const shared_ptr<Genome> &genome, const shared_ptr<Fitness> &fitness,
			size_t introducedInGeneration)
	: Individual(genome, fitness, introducedInGeneration), m_parent1(-1), m_parent2(-1), m_ownIndex(-1)
{
}

shared_ptr<Individual> LensGAIndividual::createNew(const shared_ptr<eatk::Genome> &genome, const shared_ptr<eatk::Fitness> &fitness,
			size_t introducedInGeneration) const
{
	return make_shared<LensGAIndividual>(genome, fitness, introducedInGeneration);
}

bool_t LensGAIndividual::read(serut::SerializationInterface &si)
{
	auto pGenome = dynamic_cast<LensGAGenome *>(genomePtr());
	auto pFitness = dynamic_cast<LensGAFitness *>(fitnessPtr());
	if (!pGenome || !pFitness)
		return "Genome or fitness is not of the correct type for a LensGAIndividual";

	// TODO: m_lastMutationGeneration? m_introducedInGeneration?

	if (!si.readInt32(&m_parent1) || !si.readInt32(&m_parent2) || !si.readInt32(&m_ownIndex))
		return si.getErrorString();

	bool_t r;
	if (!(r = pGenome->read(si)))
		return r;
	if (!(r = pFitness->read(si)))
		return r;
	return true;
}

bool_t LensGAIndividual::write(serut::SerializationInterface &si) const
{
	auto pGenome = dynamic_cast<const LensGAGenome *>(genomePtr());
	auto pFitness = dynamic_cast<const LensGAFitness *>(fitnessPtr());
	if (!pGenome || !pFitness)
		return "Genome or fitness is not of the correct type for a LensGAIndividual";

	// TODO: m_lastMutationGeneration? m_introducedInGeneration?
	
	if (!si.writeInt32(m_parent1) || !si.writeInt32(m_parent2) || !si.writeInt32(m_ownIndex))
		return si.getErrorString();

	bool_t r;
	if (!(r = pGenome->write(si)))
		return r;
	if (!(r = pFitness->write(si)))
		return r;

	return true;
}

LensGAIndividualCreation::LensGAIndividualCreation(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
			 size_t numBasisFunctions, size_t numSheets, bool allowNegative,
			 size_t numObjectives)
	 : m_rng(rng), m_numBasisFunctions(numBasisFunctions), m_numSheets(numSheets), m_allowNegative(allowNegative),
	   m_numObjectives(numObjectives)
{
}

LensGAIndividualCreation::~LensGAIndividualCreation()
{
}

shared_ptr<eatk::Genome> LensGAIndividualCreation::createUnInitializedGenome()
{
	shared_ptr<LensGAGenome> genome = make_shared<LensGAGenome>(m_numBasisFunctions, m_numSheets);
	return genome;
}

shared_ptr<eatk::Genome> LensGAIndividualCreation::createInitializedGenome()
{
	shared_ptr<LensGAGenome> genome = make_shared<LensGAGenome>(m_numBasisFunctions, m_numSheets);

	for (auto &x : genome->m_weights)
		x = m_rng->getRandomFloat();

	if (m_allowNegative)
	{
		for (auto &x : genome->m_weights)
			x -= 0.5f;
	}

	for (auto &s : genome->m_sheets)
		s = m_rng->getRandomFloat();

	return genome;
}

shared_ptr<eatk::Fitness> LensGAIndividualCreation::createEmptyFitness()
{
	return make_shared<LensGAFitness>(m_numObjectives);
}

shared_ptr<eatk::Individual> LensGAIndividualCreation::createReferenceIndividual()
{
	return std::make_shared<LensGAIndividual>(nullptr, nullptr);
}

}
