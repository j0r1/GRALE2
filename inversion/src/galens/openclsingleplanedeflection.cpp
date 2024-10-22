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
					   const std::string &deflectionKernelCode, const std::string &lensRoutineName,
                       bool uploadFullParameters,
                       int devIdx
					   ) // TODO: calculate betas from this as well?
{
    if (m_init)
        return "Already initialized";
    
    auto cl = make_unique<OpenCLMultiKernel<3>>();

    string libName = cl->getLibraryName();
    if (!cl->loadLibrary(libName))
        return "Can't load library '" + libName + "': " + cl->getErrorString();

    if (devIdx < 0)
    {
        devIdx = cl->getRotatedDeviceIndex();
        if (devIdx < 0) // something went wrong
            return "Can't get rotated device index: " + cl->getErrorString();
    }
    
    if (!cl->init(devIdx))
        return "Can't init device " + to_string(devIdx) + ": " + cl->getErrorString();

    auto cleanup = [this, &cl](string r)
    {
        m_clThetas.dealloc(*cl);
        m_clIntParams.dealloc(*cl);
        m_clFloatParams.dealloc(*cl);
        m_clAllResults.dealloc(*cl);
        m_clChangedParamsBuffer.dealloc(*cl);
        m_clChangeableParamIndices.dealloc(*cl);
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

    if (thetasFlat.size() == 0)
        return cleanup("No theta positions");

    bool_t r;
    if (!(r = m_clThetas.realloc(*cl, ctx, sizeof(cl_float)*thetasFlat.size())) ||
        !(r = m_clThetas.enqueueWriteBuffer(*cl, queue, thetasFlat, true)))
        return "Can't copy thetas to GPU: " + cleanup(r.getErrorString());

    vector<cl_int> clIntParams;
    for (auto x : templateIntParameters)
        clIntParams.push_back(x);
    clIntParams.push_back(-45678); // sentinel, and make sure that something is present
    
    vector<cl_float> clFloatParams;
    for (auto x : templateFloatParameters)
        clFloatParams.push_back(x);

    if (!(r = m_clIntParams.realloc(*cl, ctx, sizeof(cl_int)*clIntParams.size())) ||
        !(r = m_clIntParams.enqueueWriteBuffer(*cl, queue, clIntParams, true)))
        return "Can't copy integer parameters to GPU: " + cleanup(r.getErrorString());

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

    // compile kernel(s)

    string deflectionKernel = deflectionKernelCode;
    deflectionKernel += R"XYZ(

__kernel void calculateDeflectionAngles(int numPoints, int numParamSets, int numFloatParams,
                                    __global const int *pIntParams,
                                    __global const float *pFloatParamsBase,
                                    __global const float *pThetas,
                                    __global float *pAllResults
                                    )
{
    const int i = get_global_id(0);
    if (i >= numPoints)
        return;
    const int paramSet = get_global_id(1);
    if (paramSet >= numParamSets)
        return;

    float2 theta;
    int thetaOffset = i*2;
    theta.x = pThetas[thetaOffset+0];
    theta.y = pThetas[thetaOffset+1];

    __global const float *pFloatParams = pFloatParamsBase + paramSet*numFloatParams;
    LensQuantities r = )XYZ" + lensRoutineName + R"XYZ((theta, pIntParams, pFloatParams);

    __global float *pResults = pAllResults + 6*numPoints*paramSet; // for each genome, 6 values
    int resultOffset = 6*i; // 6 resulting values per point

    pResults[resultOffset + 0] = r.alphaX;
    pResults[resultOffset + 1] = r.alphaY;
    pResults[resultOffset + 2] = r.axx;
    pResults[resultOffset + 3] = r.ayy;
    pResults[resultOffset + 4] = r.axy;
    pResults[resultOffset + 5] = r.potential;
}
)XYZ";

    // cerr << "Compiling OpenCL program:" << endl;
    // cerr << deflectionKernel << endl;

    string faillog;
    if (!cl->loadKernel(deflectionKernel, "calculateDeflectionAngles", faillog, 0))
    {
        cerr << faillog << endl;
        return cleanup(cl->getErrorString());
    }

    // Upload full parameters, or only the changed ones (we'll need an extra kernel)
    if (!uploadFullParameters)
    {
        vector<cl_int> clChangeableParams;
        for (auto x : changeableParameterIndices)
            clChangeableParams.push_back((cl_int)x);
        clChangeableParams.push_back(-67890); // sentinel and avoid empty vector

        if (!(r = m_clChangeableParamIndices.realloc(*cl, ctx, sizeof(cl_int)*clChangeableParams.size())) ||
            !(r = m_clChangeableParamIndices.enqueueWriteBuffer(*cl, queue, clChangeableParams, true)))
            return "Can't copy changeable parameter indices to GPU: " + cleanup(r.getErrorString());

        // And create kernel to incorporate changed parameters into full parameters

        string src = R"XYZ(

__kernel void fillInChangedParameters(__global const float *pChangedParamsBase,
                                      __global float *pFloatParamsBase,
                                      __global const int *pIndices,
                                      int numFloatParams,
                                      int numChangedParams,
                                      int numParamSets)
{
    const int i = get_global_id(0);
    if (i >= numChangedParams)
        return;
    const int paramSet = get_global_id(1);
    if (paramSet >= numParamSets)
        return;

    __global const float *pChangedParams = pChangedParamsBase + paramSet*numChangedParams;
    __global float *pFloatParams = pFloatParamsBase + paramSet*numFloatParams;
    int destIdx = pIndices[i];
    pFloatParams[destIdx] = pChangedParams[i];
}

)XYZ";

        if (!cl->loadKernel(src, "fillInChangedParameters", faillog, 1))
        {
            cerr << faillog << endl;
            return cleanup(cl->getErrorString());
        }
    }

    // Ok

    m_uploadFullParameters = uploadFullParameters;
    m_floatParamsCopy = clFloatParams;
    m_numPoints = thetas.size();
    m_numFloatParams = templateFloatParameters.size();
    m_currentNumParamSets = 0;
    m_maxNumParamSets = 0;
    m_changeableParameterIndices = changeableParameterIndices;
    m_cl = move(cl);
    m_init = true;
    m_devIdx = devIdx;

    return true;
}

void OpenCLSinglePlaneDeflection::destroy()
{
    if (!m_init)
        return;

    m_clThetas.dealloc(*m_cl);
    m_clIntParams.dealloc(*m_cl);
    m_clFloatParams.dealloc(*m_cl);
    m_clAllResults.dealloc(*m_cl);
    m_clChangedParamsBuffer.dealloc(*m_cl);
    m_clChangeableParamIndices.dealloc(*m_cl);

    m_cl = nullptr;
    m_init = false;
}

bool_t OpenCLSinglePlaneDeflection::calculateDeflection(const std::vector<float> &parameters,  // should have numChangebleParams * numParamSets length
									  std::vector<Vector2Df> &allAlphas,
									  std::vector<float> &allAxx,
									  std::vector<float> &allAyy,
									  std::vector<float> &allAxy,
									  std::vector<float> &allPotentials)
{
    if (!m_init)
        return "Not initialized";

    size_t changeableSize = m_changeableParameterIndices.size();
    
    if (parameters.size()%changeableSize != 0)
        return "Invalid number of parameters";
    
    size_t numParamSets = parameters.size()/changeableSize;
    if (numParamSets == 0)
        return "No parameters specified";

    bool_t r;
    auto ctx = m_cl->getContext();
    auto queue = m_cl->getCommandQueue();

    // Need to do this here, in the following 'if' we'll be using the size
    m_allResultsBuffer.resize(numParamSets*m_numPoints*6); // 6 properties per point

    if (numParamSets > m_maxNumParamSets) // We need to prepare/upload some things
    {
        m_maxNumParamSets = numParamSets;

        vector<cl_float> clFloatParams;
        for (size_t i = 0 ; i < numParamSets ; i++)
            for (auto x : m_floatParamsCopy)
                clFloatParams.push_back(x);
        clFloatParams.push_back(-45678); // sentinel, and make sure that something is present

        if (!(r = m_clFloatParams.realloc(*m_cl, ctx, sizeof(cl_float)*clFloatParams.size())) ||
            !(r = m_clFloatParams.enqueueWriteBuffer(*m_cl, queue, clFloatParams, true)))
            return "Can't copy floating point parameters to GPU: " + r.getErrorString();

        swap(m_allFloatParams, clFloatParams);
        
        // allocate for outputs
        if (!(r = m_clAllResults.realloc(*m_cl, ctx, sizeof(cl_float)*m_allResultsBuffer.size())))
            return "Can't allocate results buffer: " + r.getErrorString();

        if (!m_uploadFullParameters)
        {
            if (!(r = m_clChangedParamsBuffer.realloc(*m_cl, ctx, sizeof(cl_float)*(numParamSets*changeableSize + 1))))
                return "Can't allocate changed parameters buffer: " + r.getErrorString();
        }
    }
    m_currentNumParamSets = numParamSets;
    
    cl_int clNumParamSets = (cl_int)numParamSets;
    cl_int clNumFloatParams = (cl_int)m_numFloatParams;

    if (m_uploadFullParameters)
    {
        // Fill in the new parameters in the total parameters on the CPU, and upload this
        if (parameters.size() > 0)
        {
            for (size_t i = 0 ; i < numParamSets ; i++)
            {
                float *pFloatParams = m_allFloatParams.data() + m_numFloatParams*i;
                const float *pParamSet = parameters.data() + i*changeableSize;
                for (size_t j = 0 ; j < changeableSize ; j++)
                {
                    size_t destIdx = m_changeableParameterIndices[j];
                    assert(destIdx < m_numFloatParams);

                    pFloatParams[destIdx] = pParamSet[j];
                }
            }
        }
        // Upload this to the GPU
        if (!(r = m_clFloatParams.enqueueWriteBuffer(*m_cl, queue, m_allFloatParams, true)))
            return "Can't copy floating point parameters to GPU: " + r.getErrorString();
    }    
    else // Upload the changed parameters to the GPU and fill in using a kernel
    {
        if (!(r = m_clChangedParamsBuffer.enqueueWriteBuffer(*m_cl, queue, parameters, true)))
            return "Can't copy changed parameters to GPU: " + r.getErrorString();

        // trigger a kernel to fill in the changed parameters into the full ones
        auto kernel = m_cl->getKernel(1);
        cl_int clNumChangeable = (cl_int)changeableSize;
        cl_int err = 0;

        auto setKernArg = [this, &err, &kernel](cl_uint idx, size_t size, const void *ptr)
        {
            err |= m_cl->clSetKernelArg(kernel, idx, size, ptr);
        };

        setKernArg(0, sizeof(cl_mem), (void*)&m_clChangedParamsBuffer.m_pMem);
        setKernArg(1, sizeof(cl_mem), (void*)&m_clFloatParams.m_pMem);
        setKernArg(2, sizeof(cl_mem), (void*)&m_clChangeableParamIndices.m_pMem);
        setKernArg(3, sizeof(cl_int), (void*)&clNumFloatParams);
        setKernArg(4, sizeof(cl_int), (void*)&clNumChangeable);
        setKernArg(5, sizeof(cl_int), (void*)&clNumParamSets);
        if (err != CL_SUCCESS)
            return "Error setting kernel arguments";

        size_t workSize[2] = { changeableSize, numParamSets };
        err = m_cl->clEnqueueNDRangeKernel(queue, kernel, 2, nullptr, workSize, nullptr, 0, nullptr, nullptr);
        if (err != CL_SUCCESS)
            return "Error enqueing parameter fill in kernel: " + to_string(err);
    }

    auto kernel = m_cl->getKernel(0);
    cl_int clNumPoints = (cl_int)m_numPoints;
    cl_int err = 0;

    auto setKernArg = [this, &err, &kernel](cl_uint idx, size_t size, const void *ptr)
    {
        err |= m_cl->clSetKernelArg(kernel, idx, size, ptr);
    };

    // TODO: I guess we don't need to do this every time? Perhaps for the integer stuff
    //       we do, so that the pointer is valid
    setKernArg(0, sizeof(cl_int), (void*)&clNumPoints);
    setKernArg(1, sizeof(cl_int), (void*)&clNumParamSets);
    setKernArg(2, sizeof(cl_int), (void*)&clNumFloatParams);
    setKernArg(3, sizeof(cl_mem), (void*)&m_clIntParams.m_pMem);
    setKernArg(4, sizeof(cl_mem), (void*)&m_clFloatParams.m_pMem);
    setKernArg(5, sizeof(cl_mem), (void*)&m_clThetas.m_pMem);
    setKernArg(6, sizeof(cl_mem), (void*)&m_clAllResults.m_pMem);
    if (err != CL_SUCCESS)
        return "Error setting kernel arguments";

    size_t workSize[2] = { m_numPoints, numParamSets };
    err = m_cl->clEnqueueNDRangeKernel(queue, kernel, 2, nullptr, workSize, nullptr, 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
        return "Error enqueing calculation kernel: " + to_string(err);

    if (!(r = m_clAllResults.enqueueReadBuffer(*m_cl, queue, m_allResultsBuffer, nullptr, nullptr, true)))
        return "Error reading calculation results: " + r.getErrorString();

    // Copy the results to the arrays in the parameters
    size_t totalPoints = m_numPoints*numParamSets;
    allAlphas.resize(totalPoints);
    allAxx.resize(totalPoints);
    allAyy.resize(totalPoints);
    allAxy.resize(totalPoints);
    allPotentials.resize(totalPoints);

    for (size_t i = 0, j = 0 ; i < totalPoints ; i++, j += 6)
    {
        assert(j+5 < m_allResultsBuffer.size());
        allAlphas[i] = Vector2Df(m_allResultsBuffer[j+0], m_allResultsBuffer[j+1]);
        allAxx[i] = m_allResultsBuffer[j+2];
        allAyy[i] = m_allResultsBuffer[j+3];
        allAxy[i] = m_allResultsBuffer[j+4];
        allPotentials[i] = m_allResultsBuffer[j+5];
    }

    return true;
}

std::unique_ptr<OpenCLSinglePlaneDeflectionInstance> OpenCLSinglePlaneDeflectionInstance::s_instance;
std::mutex OpenCLSinglePlaneDeflectionInstance::s_instanceMutex;
std::set<uint64_t> OpenCLSinglePlaneDeflectionInstance::s_users;
bool OpenCLSinglePlaneDeflectionInstance::s_initTried = false;

errut::bool_t OpenCLSinglePlaneDeflectionInstance::initInstance(uint64_t userId,const std::vector<Vector2Df> &thetas, // already transformed into the correct units
					   const std::vector<int> &templateIntParameters, // these cannot change
					   const std::vector<float> &templateFloatParameters, // only floating point params can change
					   const std::vector<size_t> changeableParameterIndices,
					   const std::string &deflectionKernelCode, const std::string &lensRoutineName,
					   bool uploadFullParameters,
					   int devIdx)
{
    lock_guard<mutex> guard(s_instanceMutex);

    if (s_users.find(userId) != s_users.end())
        return "An OpenCLSinglePlaneDeflectionInstance is already registered for this user";

    if (devIdx < 0)
        devIdx = -1;

    if (s_instance.get()) // Already initialized
    {
        if (devIdx != s_instance->getRequestedDeviceIndex())
            return "An OpenCLSinglePlaneDeflectionInstance already exists, but for a different device index (" + to_string(s_instance->getDeviceIndex()) + ") than requested (" + to_string(devIdx) + ")";

        cerr << "INFO: registered user in OpenCLSinglePlaneDeflectionInstance (2): " << userId << endl;
        s_users.insert(userId);
        return true;
    }

    if (s_initTried)
        return "GPU initialization failed before";
    s_initTried = true;

    unique_ptr<OpenCLSinglePlaneDeflectionInstance> oclCalc = make_unique<OpenCLSinglePlaneDeflectionInstance>();
    bool_t r = oclCalc->init(thetas, templateIntParameters, templateFloatParameters, changeableParameterIndices,
                             deflectionKernelCode, lensRoutineName, uploadFullParameters, devIdx);
    oclCalc->m_requestedDevIdx = devIdx; // TODO: store this in a better way?

    if (!r)
        return "Can't init GPU: " + r.getErrorString();

    s_users.insert(userId);
    cerr << "INFO: registered user in OpenCLSinglePlaneDeflectionInstance: " << userId << endl;
    cerr << "INFO: requested GPU " << oclCalc->getRequestedDeviceIndex() << ", got " << oclCalc->getDeviceIndex() << endl;

    s_instance = move(oclCalc);
    return true;
}

void OpenCLSinglePlaneDeflectionInstance::releaseInstance(uint64_t userId)
{
    lock_guard<mutex> guard(s_instanceMutex);

    auto it = s_users.find(userId);
    if (it == s_users.end())
    {
        cerr << "WARNING: OpenCLSinglePlaneDeflectionInstance::releaseInstance: can't find user " << userId << endl;
        return;
    }
    s_users.erase(it);

    if (s_users.empty())
    {
        cerr << "INFO: no more users of OpenCLSinglePlaneDeflectionInstance, removing static instance" << endl;
        s_initTried = false;
        s_instance = nullptr;
    }
}

OpenCLSinglePlaneDeflectionInstance &OpenCLSinglePlaneDeflectionInstance::instance()
{
    auto pInst = s_instance.get();
    assert(pInst);
    return *pInst;
}

OpenCLSinglePlaneDeflectionInstance::OpenCLSinglePlaneDeflectionInstance()
{

}

OpenCLSinglePlaneDeflectionInstance::~OpenCLSinglePlaneDeflectionInstance()
{

}

// This is called on a new iteration, so that we know how much genomes are
// coming in before we're going to start the GPU calculation
void OpenCLSinglePlaneDeflectionInstance::setTotalGenomesToCalculate(size_t num)
{
    lock_guard<mutex> guard(m_mutex);
    assert(num > 0);
    assert(m_totalGenomesToCalculate == 0 || m_totalGenomesToCalculate == num);

    m_totalGenomesToCalculate = num;
}

bool_t OpenCLSinglePlaneDeflectionInstance::scheduleCalculation(const eatk::FloatVectorGenome &genome)
{
    lock_guard<mutex> guard(m_mutex);

    assert(!m_calculationDone);

    const auto *pGenome = &genome;
    assert(m_genomeOffsets.find(pGenome) == m_genomeOffsets.end()); // Make sure we don't have an entry yet

    assert(genome.getValues().size() == m_changeableParameterIndices.size());
    size_t offset = m_floatBuffer.size()/m_changeableParameterIndices.size();
    for (auto x : genome.getValues())
        m_floatBuffer.push_back(x);

    m_genomeOffsets[pGenome] = offset; // Keep track where the results for this genome will be stored

    if (offset+1 == m_totalGenomesToCalculate)
    {
        //cerr << "Got " << m_totalGenomesToCalculate << " genomes, need to calculate!" << endl;
        bool_t r = calculateDeflection(m_floatBuffer, m_allAlphas, m_allAxx, m_allAyy, m_allAxy, m_allPotentials);
        if (!r)
            return r;
        //cerr << "Delections calculated successfully" << endl;
        m_calculationDone = true;

        // cerr << "GPU all alphas: " << endl;
        // for (auto a : m_allAlphas)
        //     cerr << a.getX() << " " << a.getY() << endl;
    }
    return true;
}

bool OpenCLSinglePlaneDeflectionInstance::getResultsForGenome(const eatk::FloatVectorGenome &genome,
	                         vector<Vector2Df> &alphas, vector<float> &axx,
							 vector<float> &ayy, vector<float> &axy,
							 vector<float> &potential)
{
    lock_guard<mutex> guard(m_mutex);
    if (!m_calculationDone)
        return false;

    const auto *pGenome = &genome;
    auto it = m_genomeOffsets.find(pGenome);
    assert(it != m_genomeOffsets.end());

    int offset = it->second;
    int pointOffsetStart = offset*m_numPoints;
    int pointOffsetEnd = pointOffsetStart + m_numPoints;
    // cerr << "offset for genome " << (void*)pGenome << " is " << offset 
    //      << ", numPoints = " << m_numPoints << endl;

    assert(pointOffsetStart < m_allAlphas.size());
    assert(pointOffsetEnd <= m_allAlphas.size());

    alphas.resize(m_numPoints);
    axx.resize(m_numPoints);
    ayy.resize(m_numPoints);
    axy.resize(m_numPoints);
    potential.resize(m_numPoints);

    assert(sizeof(Vector2Df) == sizeof(float)*2);
    memcpy(alphas.data(), m_allAlphas.data()+pointOffsetStart, m_numPoints*2*sizeof(float));
    memcpy(axx.data(), m_allAxx.data()+pointOffsetStart, m_numPoints*sizeof(float));
    memcpy(ayy.data(), m_allAyy.data()+pointOffsetStart, m_numPoints*sizeof(float));
    memcpy(axy.data(), m_allAxy.data()+pointOffsetStart, m_numPoints*sizeof(float));
    memcpy(potential.data(), m_allPotentials.data()+pointOffsetStart, m_numPoints*sizeof(float));

    m_genomeOffsets.erase(it);
    if (m_genomeOffsets.empty())
    {
        //cerr << "All results retrieved, ready to start again" << endl;
        m_calculationDone = false;
        m_totalGenomesToCalculate = 0;
        m_floatBuffer.resize(0);
    }
    return true;
}

} // end namespace