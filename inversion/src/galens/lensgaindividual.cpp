#include "lensgaindividual.h"
#include <sstream>

using namespace mogal2;
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

LensGAGenome::LensGAGenome(size_t numWeights, size_t numSheetValues)
{
    m_weights.resize(numWeights);
    m_sheets.resize(numSheetValues);
    m_scaleFactor = numeric_limits<float>::quiet_NaN(); // Will be returned by the fitness calculation
}

shared_ptr<mogal2::Genome> LensGAGenome::createCopy(bool copyContents) const
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

LensGAIndividual::LensGAIndividual(shared_ptr<Genome> genome, shared_ptr<Fitness> fitness,
            size_t introducedInGeneration)
    : Individual(genome, fitness, introducedInGeneration), m_parent1(-1), m_parent2(-1)
{
}

shared_ptr<Individual> LensGAIndividual::createNew(shared_ptr<mogal2::Genome> genome, shared_ptr<mogal2::Fitness> fitness,
            size_t introducedInGeneration) const
{
    return make_shared<LensGAIndividual>(genome, fitness, introducedInGeneration);
}

LensGAIndividualCreation::LensGAIndividualCreation(const std::shared_ptr<mogal2::RandomNumberGenerator> &rng,
		     size_t numBasisFunctions, size_t numSheets, bool allowNegative,
			 size_t numObjectives)
	 : m_rng(rng), m_numBasisFunctions(numBasisFunctions), m_numSheets(numSheets), m_allowNegative(allowNegative),
	   m_numObjectives(numObjectives)
{
}

LensGAIndividualCreation::~LensGAIndividualCreation()
{
}

shared_ptr<mogal2::Genome> LensGAIndividualCreation::createInitializedGenome()
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

shared_ptr<mogal2::Fitness> LensGAIndividualCreation::createEmptyFitness()
{
    return make_shared<LensGAFitness>(m_numObjectives);
}

shared_ptr<mogal2::Individual> LensGAIndividualCreation::createReferenceIndividual()
{
    return std::make_shared<LensGAIndividual>(nullptr, nullptr);
}

}