#include "openclsingleplanedeflection.h"

using namespace std;
using namespace errut;

namespace grale
{

using namespace oclutils;

OpenCLSinglePlaneDeflection::OpenCLSinglePlaneDeflection()
{
}

OpenCLSinglePlaneDeflection::~OpenCLSinglePlaneDeflection()
{
}

bool_t OpenCLSinglePlaneDeflection::init(const std::vector<Vector2Df> &thetas, // already transformed into the correct units
					   const std::vector<int> &templateIntParameters, // these cannot change
					   const std::vector<float> &templateFloatParameters, // only floating point params can change
					   const std::vector<size_t> changeableParameterIndices,
					   size_t numParamSets, // number of genomes for example
					   const std::string &deflectionKernelCode,
                       size_t devIdx
					   ) // TODO: calculate betas from this as well?
{
    if (m_init)
        return "Already initialized";
    
    auto cl = make_unique<OpenCLMultiKernel<3>>();

    string libName = cl->getLibraryName();
    if (!cl->loadLibrary(libName))
        return "Can't load library '" + libName + "': " + cl->getErrorString();
    
    if (!cl->init(devIdx))
        return "Can't init device " + to_string(devIdx) + ": " + cl->getErrorString();

    auto cleanup = [this, &cl](string r)
    {
        m_clThetas.dealloc(*cl);
        m_clIntParams.dealloc(*cl);
        m_clFloatParams.dealloc(*cl);
        return r;
    };

    auto ctx = cl->getContext();
    auto queue = cl->getCommandQueue();

    // Copy the thetas to GPU    
    vector<cl_float> thetasFlat;
    for (auto t: thetas)
    {
        thetasFlat.push_back(t.getX());
        thetasFlat.push_back(t.getY());
    }
    
    bool_t r;
    if (!(r = m_clThetas.realloc(*cl, ctx, sizeof(cl_float)*thetasFlat.size())) ||
        !(r = m_clThetas.enqueueWriteBuffer(*cl, queue, thetasFlat, true)))
        return "Can't copy thetas to GPU: " + cleanup(r.getErrorString());

    vector<cl_int> clIntParams;
    for (auto x : templateIntParameters)
        clIntParams.push_back(x);
    clIntParams.push_back(-45678); // sentinel, and make sure that something is present
    
    // For the floating point parameters, we need a copy for each
    // genome, so it can change the values for that genome
    vector<cl_float> clFloatParams;
    for (size_t i = 0 ; i < numParamSets ; i++)
        for (auto x : templateFloatParameters)
            clFloatParams.push_back(x);
    clFloatParams.push_back(-45678); // sentinel, and make sure that something is present


    if (!(r = m_clIntParams.realloc(*cl, ctx, sizeof(cl_int)*clIntParams.size())) ||
        !(r = m_clIntParams.enqueueWriteBuffer(*cl, queue, clIntParams, true)))
        return "Can't copy integer parameters to GPU: " + cleanup(r.getErrorString());

    if (!(r = m_clFloatParams.realloc(*cl, ctx, sizeof(cl_float)*clFloatParams.size())) ||
        !(r = m_clFloatParams.enqueueWriteBuffer(*cl, queue, clFloatParams, true)))
        return "Can't copy floating point parameters to GPU: " + cleanup(r.getErrorString());

    // Sanity checks
    if (changeableParameterIndices.size() > templateFloatParameters.size())
        return cleanup("Too many changeable indices");

    vector<bool> indexCheck(templateFloatParameters.size(), false);
    for (auto x : changeableParameterIndices)
    {
        if (x >= templateFloatParameters.size())
            return cleanup("Invalid index " + to_string(x) + " in changeable parameter indices");
        if (indexCheck[x])
            return cleanup("Trying to change index " + to_string(x) + " multiple times");
        indexCheck[x] = true;
    }

    // TODO: compile kernel(s)
    
    swap(m_allFloatParams, clFloatParams);
    m_numFloatParams = templateFloatParameters.size();
    m_changeableParameterIndices = changeableParameterIndices;
    m_cl = move(cl);
    m_init = true;

    return true;
}

void OpenCLSinglePlaneDeflection::destroy()
{
    if (!m_init)
        return;

    m_clThetas.dealloc(*m_cl);
    m_cl = nullptr;
    m_init = false;
}

bool_t OpenCLSinglePlaneDeflection::calculateDeflection(const std::vector<float> &parameters) // should have numChangebleParams * numParamSets length
{
    return true;
}

}
