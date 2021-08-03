#include "openclcalculator.h"
#include "constants.h"
#include "compositelens.h"
#include "utils.h"
#include <sstream>
#include <sys/time.h>

using namespace errut;
using namespace std;

namespace grale
{

unique_ptr<OpenCLCalculator> OpenCLCalculator::s_pInstance;
mutex OpenCLCalculator::s_instanceMutex;
bool OpenCLCalculator::s_initTried = false;
set<uint64_t> OpenCLCalculator::s_users;

bool_t OpenCLCalculator::initInstance(int devIdx, const vector<ImagesDataExtended *> &allImages,
                                        const vector<ImagesDataExtended *> &shortImages,
                                        const vector<float> &zds,
                                        const Cosmology &cosm,
                                        const vector<vector<shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
                                        const vector<shared_ptr<GravitationalLens>> &unscaledLensesPerPlane,
                                        uint64_t userId
                                        )
{
    lock_guard<mutex> guard(s_instanceMutex);

    if (s_users.find(userId) != s_users.end())
        return "An OpenCLCalculator instance is already registered for this user";

    if (devIdx < 0)
        devIdx = -1;

    if (s_pInstance.get()) // Already initialized
    {
        if (devIdx != s_pInstance->getDeviceIndex())
            return "An OpenCLCalculator instance already exists, but for a different device index (" + to_string(s_pInstance->getDeviceIndex()) + ") than requested (" + to_string(devIdx) + ")";

        cerr << "INFO: registered user in OpenCLCalculator (2): " << userId << endl;
        s_users.insert(userId);
        return true;
    }

    if (s_initTried)
        return "GPU initialization failed before";
    s_initTried = true;

    unique_ptr<OpenCLCalculator> oclCalc = make_unique<OpenCLCalculator>();
    bool_t r = oclCalc->initAll(devIdx, allImages, shortImages, zds, cosm, planeBasisLenses, unscaledLensesPerPlane);
    if (!r)
        return "Can't init GPU: " + r.getErrorString();

    s_users.insert(userId);
    cerr << "INFO: registered user in OpenCLCalculator: " << userId << endl;

    s_pInstance = move(oclCalc);
    return true;
}

void OpenCLCalculator::releaseInstance(uint64_t userId)
{
    lock_guard<mutex> guard(s_instanceMutex);

    auto it = s_users.find(userId);
    if (it == s_users.end())
    {
        cerr << "WARNING: OpenCLCalculator::releaseInstance: can't find user " << userId << endl;
        return;
    }
    s_users.erase(it);

    if (s_users.empty())
    {
        cerr << "INFO: no more users of OpenCLCalculator, removing static instance" << endl;
        s_initTried = false;
        s_pInstance = nullptr;
    }
}

OpenCLCalculator &OpenCLCalculator::instance()
{
    OpenCLCalculator *pInst = s_pInstance.get();
    assert(pInst);
    return *pInst;
}

void OpenCLCalculator::setGenomesToCalculate(size_t s)
{
    lock_guard<mutex> guard(m_mutex);
    m_numGenomesToCalculate = s;
}

#define CHECKERRORSTATE do { if (m_errorState) return "Previous error was encountered"; } while(0)

bool_t OpenCLCalculator::startNewBackprojection(const LensGAGenome &g)
{
    lock_guard<mutex> guard(m_mutex);

    CHECKERRORSTATE;

    if (m_hasCalculated) // In this case we're starting a new population's fitness calculation
    {
        m_hasCalculated = false;
        m_uploadedWeights = false;
        m_states.clear();
        m_nextGenomeIndex = 0;
    }

    size_t idx = m_nextGenomeIndex;
    m_states[&g] = State(m_nextGenomeIndex);
    m_nextGenomeIndex++;

    size_t numBfWeights = g.m_weights.size();
    size_t numSheetWeights = g.m_sheets.size();
    assert(numBfWeights + numSheetWeights == m_numWeights);

    size_t startBfIdx = idx * m_numWeights;
    size_t endBfIdx = startBfIdx + m_numWeights;

    if (endBfIdx > m_allBasisFunctionWeights.size())
        m_allBasisFunctionWeights.resize(endBfIdx);

    if (numSheetWeights > 0)
    {
        // Store the sheet weights at the correct position
        size_t genomeOffset = 0;
        assert(m_numPlanes == g.m_sheets.size());
        for (size_t i = 0 ; i < m_numPlanes ; i++)
        {
            assert(i < m_planeWeightOffsets.size() && (i+1) < m_planeWeightOffsets.size());
            size_t start = m_planeWeightOffsets[i];
            size_t end = m_planeWeightOffsets[i+1];

            assert(end > start && end > 0);
            size_t num = end-start;

            assert(genomeOffset + num - 1 <= g.m_weights.size());
            memcpy(m_allBasisFunctionWeights.data() + startBfIdx + start, g.m_weights.data() + genomeOffset, sizeof(float)*(num-1));
            genomeOffset += (num-1);

            assert(i < g.m_sheets.size());
            m_allBasisFunctionWeights[startBfIdx + end-1] = g.m_sheets[i];
        }
        assert(genomeOffset == g.m_weights.size());
    }
    else
        memcpy(m_allBasisFunctionWeights.data() + startBfIdx, g.m_weights.data(), sizeof(float)*g.m_weights.size());

    return true;
}

bool_t OpenCLCalculator::scheduleUploadAndCalculation(const LensGAGenome &g, int &calculationIdentifier,
                                                      const std::vector<std::pair<float,float>> &scaleFactors, bool useShort)
{
    bool_t r;

    {
        lock_guard<mutex> guard(m_mutex);

        CHECKERRORSTATE;
        
        const FullOrShortClMem &fullOrShort = (useShort)?m_short:m_full;
        size_t numSteps = scaleFactors.size();
        if (!m_beingScheduled.get())
        {
            auto getContextMemory = [this]() -> unique_ptr<ContextMemory>
            {
                if (m_recycledContextMemory.size() > 0)
                {
                    auto r = move(m_recycledContextMemory.back());
                    m_recycledContextMemory.pop_back();
                    r->resizeCPUBuffers();
                    return r;
                }
                return make_unique<ContextMemory>();
            };

            m_beingScheduled = make_unique<CalculationContext>(m_nextCalculationContextIdentifier++, numSteps, m_numWeights,
                                                               m_common, fullOrShort, getContextMemory());
        }
        else
        {
            if (m_beingScheduled->m_numSteps != numSteps || m_beingScheduled->m_fullOrShort.m_numPoints != fullOrShort.m_numPoints)
            {
                calculationIdentifier = -1; // signal that we can't schedule yet
                return true; // This is not an error
            }
        }

        // Ok, at this point we have a context that's compatible
        auto it = m_states.find(&g);
        if (it == m_states.end())
        {
            m_errorState = true;
            return "Can't find genome state";
        }
        State &state = it->second;

        // TODO: check for some memory limit?
        // Note that the lock is still acquired!
        if (!(r = m_beingScheduled->schedule(*this, g, state.m_genomeIndex, scaleFactors)))
        {
            m_errorState = true;
            return "Can't schedule calculation: " + r.getErrorString();
        }

        calculationIdentifier = m_beingScheduled->m_identifier;
    } // releases the lock
    
    // If the GPU isn't busy, it can start calculating the next piece of work
    if (!(r = checkCalculateScheduledContext())) // sets errorstate itself
        return r;

    return true;
}

bool_t OpenCLCalculator::checkCalculateScheduledContext()
{
    lock_guard<mutex> guard(m_mutex);

    CHECKERRORSTATE;

    if (m_beingCalculated.get()) // Already something being calculated
        return true; // Nothing to do

    if (!m_beingScheduled.get()) // Nothing to do
        return true;

    if (m_numGenomesToCalculate != m_states.size()) // Don't have all genome weights yet, bail for now
        return true;

    // Perhaps the weights still need to be uploaded!
    if (!m_uploadedWeights)
    {
        m_uploadedWeights = true;

        bool_t r;
        assert(m_allBasisFunctionWeights.size() > 0);
        if (!(r = m_common.m_devAllWeights.realloc(*this, m_allBasisFunctionWeights)) ||
            !(r = m_common.m_devAllWeights.enqueueWriteBuffer(*this, m_allBasisFunctionWeights)) )
        {
            m_errorState = true;
            return "Error reallocating/uploading GPU buffer for basis function weights" + r.getErrorString();
        }
    }

    cl_event evt;
    bool_t r;
    if (!(r = m_beingScheduled->calculate(*this, &evt)))
    {
        m_errorState = true;
        return r;
    }

    // Do this after the 'calculate' in case an error is returned. That would cause the destructor
    // to keep waiting for the calculation to finish
    swap(m_beingCalculated, m_beingScheduled);

    cl_int err = clSetEventCallback(evt, CL_COMPLETE, staticEventNotify, this);
    if (err != CL_SUCCESS)
    {
        m_errorState = true;
        return "Error setting event callback: " + to_string(err);
    }

    return true;
}

bool_t OpenCLCalculator::CalculationContext::calculate(OpenCLKernel &cl, cl_event *pEvt)
{
    bool_t r;
    if (!(r = m_mem->m_devFactors.realloc(cl, m_mem->m_allFactors)) ||
        !(r = m_mem->m_devFactors.enqueueWriteBuffer(cl, m_mem->m_allFactors)))
        return "Can't send scale factors to GPU: " + r.getErrorString();

    if (!(r = m_mem->m_devBetas.realloc(cl, m_mem->m_allBetas)))
        return "Error reserving GPU memory for backprojected points: " + r.getErrorString();

    if (!(r = m_mem->m_devGenomeIndexForBetaIndex.realloc(cl, m_mem->m_genomeIndexForBetaIndex)) ||
        !(r = m_mem->m_devGenomeIndexForBetaIndex.enqueueWriteBuffer(cl, m_mem->m_genomeIndexForBetaIndex)))
        return "Error uploading mapping from genomes in beta buffer to weight indices: " + r.getErrorString();

    cl_kernel kernel = cl.getKernel();
    cl_int clNumPoints = (cl_int)m_fullOrShort.m_numPoints;
    cl_int clNumScaleFactors = (cl_int)m_numSteps;
    cl_int clNumGenomes = (cl_int)m_mem->m_genomeIndexForBetaIndex.size();
    cl_int clGenomeSize = m_numWeights;
    cl_int err = 0;

    err |= cl.clSetKernelArg(kernel, 0, sizeof(cl_int), (void*)&clNumPoints);
    err |= cl.clSetKernelArg(kernel, 1, sizeof(cl_int), (void*)&clNumScaleFactors);
    err |= cl.clSetKernelArg(kernel, 2, sizeof(cl_int), (void*)&clNumGenomes);
    err |= cl.clSetKernelArg(kernel, 3, sizeof(cl_int), (void*)&clGenomeSize); // Doesn't need to be reset every time ; TODO: inporate in kernel as constant?
    err |= cl.clSetKernelArg(kernel, 4, sizeof(cl_mem), (void*)&m_fullOrShort.m_pDevImages);
    err |= cl.clSetKernelArg(kernel, 5, sizeof(cl_mem), (void*)&m_mem->m_devBetas.m_pMem);
    err |= cl.clSetKernelArg(kernel, 6, sizeof(cl_mem), (void*)&m_fullOrShort.m_pDevUsedPlanes);
    err |= cl.clSetKernelArg(kernel, 7, sizeof(cl_mem), (void*)&m_fullOrShort.m_pDevDpoints);
    err |= cl.clSetKernelArg(kernel, 8, sizeof(cl_mem), (void*)&m_common.m_pDevDmatrix);
    err |= cl.clSetKernelArg(kernel, 9, sizeof(cl_mem), (void*)&m_common.m_pDevAllIntParams);
    err |= cl.clSetKernelArg(kernel, 10, sizeof(cl_mem), (void*)&m_common.m_pDevAllFloatParams);
    err |= cl.clSetKernelArg(kernel, 11, sizeof(cl_mem), (void*)&m_common.m_devAllWeights.m_pMem);
    err |= cl.clSetKernelArg(kernel, 12, sizeof(cl_mem), (void*)&m_common.m_pDevCenters);
    err |= cl.clSetKernelArg(kernel, 13, sizeof(cl_mem), (void*)&m_common.m_pDevPlaneIntParamOffsets);
    err |= cl.clSetKernelArg(kernel, 14, sizeof(cl_mem), (void*)&m_common.m_pDevPlaneFloatParamOffsets);
    err |= cl.clSetKernelArg(kernel, 15, sizeof(cl_mem), (void*)&m_common.m_pDevPlaneWeightOffsets);
    err |= cl.clSetKernelArg(kernel, 16, sizeof(cl_mem), (void*)&m_mem->m_devFactors.m_pMem);
    err |= cl.clSetKernelArg(kernel, 17, sizeof(cl_mem), (void*)&m_mem->m_devGenomeIndexForBetaIndex.m_pMem);
    if (err != CL_SUCCESS)
        return "Error setting kernel arguments";

    size_t workSize[3] = { (size_t)clNumPoints, (size_t)clNumGenomes, (size_t)clNumScaleFactors };
    err = cl.clEnqueueNDRangeKernel(cl.getCommandQueue(), kernel, 3, nullptr, workSize, nullptr, 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
        return "Error enqueuing kernel: " + to_string(err);
    
    // cerr << "Enqueued kernel for calculation " << m_identifier << endl;

    // Also queue reading the results to CPU memory
    if (!(r = m_mem->m_devBetas.enqueueReadBuffer(cl, m_mem->m_allBetas, pEvt)))
        return "Error enqueuing read buffer";
    
    // cerr << "Enqueued read for calculation " << m_identifier << endl;

    return true;
}

bool_t OpenCLCalculator::CalculationContext::schedule(OpenCLKernel &cl, const LensGAGenome &g, size_t genomeWeightsIndex,
                                                      const vector<pair<float,float>> &scaleFactors)
{
    bool_t r;

    assert(m_betaIndexForGenome.find(&g) == m_betaIndexForGenome.end()); // Shouldn't exist yet
    assert(m_numSteps == scaleFactors.size());

    size_t betaIndex = m_mem->m_genomeIndexForBetaIndex.size();
    m_mem->m_genomeIndexForBetaIndex.push_back(genomeWeightsIndex);
    m_betaIndexForGenome[&g] = betaIndex;

    size_t startIdx = betaIndex * m_numSteps;
    size_t endIdx = startIdx + m_numSteps;

    if (endIdx > m_mem->m_allFactors.size())
        m_mem->m_allFactors.resize(endIdx);

    for (size_t i = 0, j = startIdx ; i < scaleFactors.size() ; i++, j++)
    {
        assert(j < endIdx && j < m_mem->m_allFactors.size());
        m_mem->m_allFactors[j] = scaleFactors[i].second; // .second is the real value to be used
    }

    size_t numGenomesScheduled = m_mem->m_genomeIndexForBetaIndex.size();
    size_t floatsCurrentlyNeeded = m_fullOrShort.m_numPoints*m_numSteps*numGenomesScheduled * 2;
    if (floatsCurrentlyNeeded != m_mem->m_allBetas.size())
        m_mem->m_allBetas.resize(floatsCurrentlyNeeded);

    return true;
}

bool_t OpenCLCalculator::isCalculationDone(const LensGAGenome &g, int calculationIdentifier, size_t *pGenomeIndex,
                                           const float **ppAllBetas, size_t *pNumGenomes, size_t *pNumSteps,
                                           size_t *pNumPoints, bool *pDone)
{
    *pDone = false;

    {
        lock_guard<mutex> lock(m_mutex);

        CHECKERRORSTATE;

        // auto logIt = [this,calculationIdentifier](const string &pref)
        // {
        //     cerr << pref << ", not done for thread " << std::this_thread::get_id() << ", calculationIdentifier = " << calculationIdentifier << endl;            
        //     if (m_beingScheduled.get())
        //         cerr << "  m_beingScheduled: " << m_beingScheduled->m_identifier << endl;
        //     if (m_beingCalculated.get())
        //         cerr << "  m_beingCalculated: " << m_beingCalculated->m_identifier << " " << ((m_beingCalculated->m_calculated)?"calculated":"not calculated") << endl;
        //     if (m_doneCalculating.get())
        //         cerr << "  m_doneCalculating: " << m_doneCalculating->m_identifier << endl;
        // };

        // First check if we can advance m_beingCalculated to m_doneCalculating
        if (!m_doneCalculating.get() && m_beingCalculated.get() && m_beingCalculated->m_calculated)
            swap(m_beingCalculated, m_doneCalculating);

        if (!m_doneCalculating.get()) // Nothing done yet
        {
            // logIt("Nothing done yet");
            return true;
        }

        if (calculationIdentifier != m_doneCalculating->m_identifier) // Is of different context, skip for now
        {
            // logIt("Mismatch identifier");
            return true;
        }
        
        assert(m_doneCalculating->m_betaIndexForGenome.find(&g) != m_doneCalculating->m_betaIndexForGenome.end());
        *pGenomeIndex = m_doneCalculating->m_betaIndexForGenome[&g];
        *ppAllBetas = m_doneCalculating->m_mem->m_allBetas.data();
        *pNumGenomes = m_doneCalculating->m_mem->m_genomeIndexForBetaIndex.size();
        *pNumSteps = m_doneCalculating->m_numSteps;
        *pNumPoints = m_doneCalculating->m_fullOrShort.m_numPoints;
        *pDone = true;
        assert(m_doneCalculating->m_mem->m_allBetas.size() == (*pNumSteps)*(*pNumGenomes)*(*pNumPoints)*2);
    }

    bool_t r;
    if (!(r = checkCalculateScheduledContext())) // sets error state itself
        return r;

    return true;
}

bool_t OpenCLCalculator::setCalculationProcessed(const LensGAGenome &g, int calculationIdentifier)
{
    lock_guard<mutex> lock(m_mutex);
    
    if (!m_doneCalculating.get() || m_doneCalculating->m_identifier != calculationIdentifier)
    {
        m_errorState = true;
        return "m_doneCalculating doesn't match the requested calculation identifier";
    }
    
    auto it = m_doneCalculating->m_betaIndexForGenome.find(&g);
    if (it == m_doneCalculating->m_betaIndexForGenome.end())
    {
        m_errorState = true;
        return "Couldn't find genome in map of calculated genomes";
    }

    m_doneCalculating->m_betaIndexForGenome.erase(it);
    
    // If everyone's results were processed, we can recycle the memory, and release this
    if (m_doneCalculating->m_betaIndexForGenome.size() == 0)
    {
        m_recycledContextMemory.push_back(move(m_doneCalculating->m_mem));
        //cleanupContextMemory(*(m_doneCalculating->m_mem));
        m_doneCalculating = nullptr;
    }

    return true;
}

void OpenCLCalculator::cleanupContextMemory(ContextMemory &mem)
{
    mem.m_devFactors.dealloc(*this);
    mem.m_devBetas.dealloc(*this);
    mem.m_devGenomeIndexForBetaIndex.dealloc(*this);
}

void OpenCLCalculator::staticEventNotify(cl_event event, cl_int eventCommandStatus, void *userData)
{
    OpenCLCalculator *pInst = reinterpret_cast<OpenCLCalculator*>(userData);
    assert(pInst);
    pInst->eventNotify(event, eventCommandStatus);
}

//#define DUMPLASTBETAS

void OpenCLCalculator::eventNotify(cl_event event, cl_int eventCommandStatus)
{
    lock_guard<mutex> lock(m_mutex);
    if (eventCommandStatus != CL_COMPLETE)
    {
        cerr << "Unexpected result from OpenCL, bailing!" << endl;
        exit(-1);
    }

    // cout << "Calculation Done!" << endl;
    // cout.flush();

    m_hasCalculated = true; // Signal that we can accept new genome weights at a later point

    assert(m_beingCalculated.get());
    // Here we just set a flag; we'll do further processing based on this in isCalculationDone
    m_beingCalculated->m_calculated = true;

#ifdef DUMPLASTBETAS
    size_t numGenomesInBetas = m_beingCalculated->m_mem->m_genomeIndexForBetaIndex.size();
    size_t numStepsInBetas = m_beingCalculated->m_numSteps;
    size_t numPointsInBetas = m_beingCalculated->m_fullOrShort.m_numPoints;
    assert(2*numGenomesInBetas*numStepsInBetas*numPointsInBetas == m_beingCalculated->m_mem->m_allBetas.size());
    cout << std::dec;
    cout << "# genomes = " << numGenomesInBetas << endl;
    cout << "# numSteps = " << numStepsInBetas << endl;
    cout << "# numPoints = " << numPointsInBetas << endl;
    for (int g = 0 ; g < numGenomesInBetas ; g++)
    {
        int genomeOffset = (numStepsInBetas*numPointsInBetas*2)*g;
        for (int f = 0 ; f < numStepsInBetas ; f++)
        {
            int scaleOffset = genomeOffset + (numPointsInBetas*2)*f;
            for (int i = 0 ; i < numPointsInBetas ; i++)
                cout << "\t" << m_beingCalculated->m_mem->m_allBetas[scaleOffset + i*2];
            cout << endl;
            for (int i = 0 ; i < numPointsInBetas ; i++)
                cout << "\t" << m_beingCalculated->m_mem->m_allBetas[scaleOffset + i*2 + 1];
            cout << endl;
            cout << endl;
        }
        cout << endl;
    }
    cout.flush();

    exit(-1);
#endif // DUMPLASTBETAS
}

OpenCLCalculator::OpenCLCalculator()
    : m_devIdx(-1), m_errorState(false)
{
}

OpenCLCalculator::~OpenCLCalculator()
{
    // TODO: wait a maximum amount of time
    
    // Wait till GPU is finished calculating and release
    cerr << "INFO: OpenCLCalculator destructor start" << endl;

    auto cleanup = [this](unique_ptr<CalculationContext> &ctx)
    {
        cleanupContextMemory(*(ctx->m_mem));
        ctx = nullptr;
    };

    bool firstIteration = true;
    bool done = false;
    while (!done)
    {
        done = true;

        lock_guard<mutex> guard(m_mutex);
        if (m_beingScheduled.get())
            cleanup(m_beingScheduled);

        if (m_beingCalculated.get())
        {
            if (!m_beingCalculated->m_calculated)
            {
                done = false;
                if (firstIteration)
                    cerr << "WARNING: waiting for destruction till GPU is done" << endl;
            }
            else
                cleanup(m_beingCalculated);
        }

        if (m_doneCalculating.get())
            cleanup(m_doneCalculating);

        firstIteration = false;
    }

    lock_guard<mutex> guard(m_mutex);

    m_full.dealloc(*this);
    if (m_shortImagesAreAllImages)
        m_short.zeroAll();
    else
        m_short.dealloc(*this);

    for (auto &ctx : m_recycledContextMemory)
        cleanupContextMemory(*ctx);

	m_perNodeCounter.reset();

    cerr << "INFO: OpenCLCalculator destructor end" << endl;
}

bool_t OpenCLCalculator::initAll(int devIdx, const vector<ImagesDataExtended *> &allImages,
                        const vector<ImagesDataExtended *> &shortImages,
                        const vector<float> &zds,
                        const Cosmology &cosm,
                        const vector<vector<shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
                        const vector<shared_ptr<GravitationalLens>> &unscaledLensesPerPlane)
{
    bool_t r;

    if (devIdx < 0)
        devIdx = -1;

    if (!(r = initGPU(devIdx)))
        return "Error initializing GPU: " + r.getErrorString();

    m_devIdx = devIdx;

    m_angularScale = ANGLE_ARCSEC; // TODO: calculate something better?
    m_potentialScale = ANGLE_ARCSEC*ANGLE_ARCSEC; // TODO

    if (!(r = setupImages(allImages, shortImages)))
        return "Can't setup images for GPU backprojection: " + r.getErrorString();

    if (!(r = setupMultiPlaneDistanceMatrix(cosm, zds)))
        return "Can't setup multi-plane distance matrix: " + r.getErrorString();

    if (!(r = setupAngularDiameterDistances(cosm, zds, allImages, m_full.m_pDevDpoints, m_full.m_pDevUsedPlanes)))
        return "Can't setup angular diameter distances for all points: " + r.getErrorString();
    if (!m_shortImagesAreAllImages)
    {
        if (!(r = setupAngularDiameterDistances(cosm, zds, shortImages, m_short.m_pDevDpoints, m_short.m_pDevUsedPlanes)))
            return "Can't setup angular diameter distances for short points: " + r.getErrorString();
    }
    else
    {
        m_short.m_pDevUsedPlanes = m_full.m_pDevUsedPlanes;
        m_short.m_pDevDpoints = m_full.m_pDevDpoints;
    }
    
    if (!(r = setupBasisFunctions(planeBasisLenses, unscaledLensesPerPlane)))
        return "Can't setup GPU code for basis functions: " + r.getErrorString();

    m_hasCalculated = true;
    m_nextCalculationContextIdentifier = 0;

    return true;
}

bool_t OpenCLCalculator::analyzeCompositeLenses(const vector<vector<shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
                                const vector<shared_ptr<GravitationalLens>> &unscaledLensesPerPlane,
                                map<string,string> &subRoutineCode,
                                vector<string> &compLensSubRoutineNames,
                                int &compLensRecursion
                                )
{
    vector<shared_ptr<GravitationalLens>> allLensModels;
    for (auto &plane : planeBasisLenses)
        for (auto &bl : plane)
            allLensModels.push_back(bl->m_pLens);
    for (auto &l : unscaledLensesPerPlane)
        allLensModels.push_back(l);
    
    map<string,string> totalSubCode;
    vector<string> totalSubNames;
    int totalMaxRecurs = 0;

    for (auto &l : allLensModels)
    {
        if (l->getLensType() != GravitationalLens::Composite)
            continue;
        const CompositeLens &compLens = static_cast<const CompositeLens&>(*l);

        map<string,string> subCode;
        vector<string> subNames;
        int recurs = compLens.findCLSubroutines(subCode, subNames, false, false);

        if (recurs > totalMaxRecurs)
            totalMaxRecurs = recurs;

        for (auto &it : subCode)
            totalSubCode[it.first] = it.second;
        
        if (totalSubNames.size() == 0)
            totalSubNames = subNames;
        else
        {
            assert(totalSubNames.size() == subNames.size());
            for (size_t i = 0 ; i < subNames.size() ; i++)
            {
                if (subNames[i].length() > 0)
                    totalSubNames[i] = subNames[i];
            }
        }
    }

    swap(subRoutineCode, totalSubCode);
    swap(compLensSubRoutineNames, totalSubNames);
    compLensRecursion = totalMaxRecurs;

    return true;
}

bool_t OpenCLCalculator::getAlphaCodeForPlane(const string &functionName,
                            const vector<shared_ptr<LensInversionBasisLensInfo>> &basisLenses,
                            const GravitationalLens *pUnscaledLens,
                            map<string,string> &subRoutineCode,
                            int &intParamCount, int &floatParamCount, int &numWeights,
                            vector<float> &centers,
                            vector<cl_int> &intParams,
                            vector<float> &floatParams,
                            string &generatedCode)
{
    int totalIntParamCount = 0, totalFloatParamCount = 0, weightCount = 0;

    auto addCodeForLens = [](const string &fname, int num, int iCnt, int fCnt, bool usePlaneScale) -> string
    {
        stringstream ss;
        assert(num >= 0 && iCnt >= 0 && fCnt >= 0);

        if (num > 1)
            ss << "    for (int i = 0 ; i < " << num << " ; i++)\n";

        ss << "    {\n";
        ss << "        const float2 center = (float2)(pCenters[0], pCenters[1]);\n";
        ss << "        LensQuantities l = " << fname << "(theta-center, pIntParams, pFloatParams);\n";
        ss << "        const float w = *pWeights;\n";
        ss << "\n";
        ss << "        pWeights++;\n";
        ss << "        pCenters += 2;\n";
        if (usePlaneScale)
        {
            ss << "        alpha.x += w*l.alphaX*scalableFunctionScale;\n";
            ss << "        alpha.y += w*l.alphaY*scalableFunctionScale;\n";
        }
        else
        {
            ss << "        alpha.x += w*l.alphaX;\n";
            ss << "        alpha.y += w*l.alphaY;\n";
        }
        
        if (iCnt > 0)
            ss << "        pIntParams += " << iCnt << ";\n";
        if (fCnt > 0)
            ss << "        pFloatParams += " << fCnt << ";\n";
        ss << "    }\n";

        return ss.str();
    };

    stringstream codeStream;
    codeStream << "float2 " << functionName << "(float2 theta, __global const int *pIntParams, __global const float *pFloatParams, __global const float *pWeights,\n"
                << "                              __global const float *pCenters, float scalableFunctionScale)\n"
                << "{\n"
                << "    float2 alpha = (float2)(0, 0);\n";

    auto saveCode = [&subRoutineCode](const GravitationalLens &l) -> string
    {
        string fnName;
        string fnCode = l.getCLProgram(fnName, false, false);
        if (l.getLensType() == GravitationalLens::Composite) // Code for compositelens will be added later
            fnCode = "";
        
        subRoutineCode[fnName] = fnCode; // May overwrite, shouldn't matter, should be the same code
        return fnName;
    };

    auto addCenter = [&centers,this](Vector2Dd pos)
    {
        centers.push_back(pos.getX()/m_angularScale);
        centers.push_back(pos.getY()/m_angularScale);
    };

    auto processParameters = [&totalIntParamCount, &totalFloatParamCount, &weightCount,
                                &intParams, &floatParams, this](const GravitationalLens &l, int &iCnt, int &fCnt) -> bool_t
    {
        if (!l.getCLParameterCounts(&iCnt, &fCnt))
            return "Can't get OpenCL parameters count for lens: " + l.getErrorString();
        assert(iCnt >= 0 && fCnt >= 0);
        totalIntParamCount += iCnt;
        totalFloatParamCount += fCnt;
        weightCount++;

        size_t io = intParams.size(), fo = floatParams.size();
        intParams.resize(io + (size_t)iCnt);
        floatParams.resize(fo + (size_t)fCnt);
        if (!l.getCLParameters(m_angularScale, m_potentialScale, intParams.data() + io, floatParams.data() + fo))
            return "Can't get OpenCL parameters for lens: " + l.getErrorString();
        return true;
    };

    bool_t r;
    string prevBfName;
    int prevBfType = -1, prevBfCount = 0, prevIntCount = -1, prevFloatCount = -1;
    for (auto &bl : basisLenses)
    {
        string fnName = saveCode(*(bl->m_pLens));
        addCenter(bl->m_center);

        int iCnt = -1, fCnt = -1;
        if (!(r = processParameters(*(bl->m_pLens), iCnt, fCnt)))
            return r;
        
        int t = bl->m_pLens->getLensType();
        if (prevBfType != t || prevIntCount != iCnt || prevFloatCount != fCnt)
        {
            if (prevBfCount)
                codeStream << addCodeForLens(prevBfName, prevBfCount, prevIntCount, prevFloatCount, true);
            prevBfType = t;
            prevBfCount = 1;
            prevIntCount = iCnt;
            prevFloatCount = fCnt;
            prevBfName = fnName;
        }
        else
            prevBfCount++;
    }

    if (prevBfCount)
        codeStream << addCodeForLens(prevBfName, prevBfCount, prevIntCount, prevFloatCount, true);

    if (pUnscaledLens)
    {
        string fnName = saveCode(*pUnscaledLens);
        addCenter({0,0});
        
        int iCnt = -1, fCnt = -1;
        if (!(r = processParameters(*pUnscaledLens, iCnt, fCnt)))
            return r;
        codeStream << addCodeForLens(fnName, 1, iCnt, fCnt, false);
    }

    codeStream << "    return alpha;\n";
    codeStream << "}\n";

    intParamCount += totalIntParamCount;
    floatParamCount += totalFloatParamCount;
    numWeights += weightCount;

    generatedCode = codeStream.str();
    return true;
}

bool_t OpenCLCalculator::getMultiPlaneTraceCode(const vector<vector<shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
                                const vector<shared_ptr<GravitationalLens>> &unscaledLensesPerPlane,
                                string &resultingCode)
{
    map<string, string> subRoutineCode;
    vector<string> compLensSubNames;
    int compLensRecursion;
    bool_t r;

    // If there's one or more composite lenses among the basis lenses, the OpenCL code for that
    // lens may itself need other simple lens models code'. Analyze this.
    if (!(r = analyzeCompositeLenses(planeBasisLenses, unscaledLensesPerPlane, subRoutineCode, compLensSubNames, compLensRecursion)))
        return "Error analyzing basis functions for composite lenses: " + r.getErrorString();

    bool haveUnscaledLenses = (unscaledLensesPerPlane.size() > 0)?true:false;
    assert((haveUnscaledLenses && (planeBasisLenses.size() == unscaledLensesPerPlane.size())) || (!haveUnscaledLenses && (unscaledLensesPerPlane.size() == 0)));

    string alphaCode;
    stringstream code;

    code << "\n";
    code << "#define MAXPLANES " << planeBasisLenses.size() << "\n";
    code << "\n";
    code << "float2 getAlpha(int lpIdx, float2 theta, __global const int *pAllIntParams, __global const float *pAllFloatParams,\n";
    code << "                __global const float *pAllWeights, __global const float *pAllCenters, __global const int *pPlaneIntParamOffsets,\n";
    code << "                __global const int *pPlaneFloatParamOffsets, __global const int *pPlaneWeightOffsets,\n";
    code << "                float scalableFunctionScale)\n";
    code << "{\n";
    
    int intParamCount = 0, floatParamCount = 0, numPlaneWeights = 0;
    vector<cl_int> planeIntParamOffsets, planeFloatParamOffsets, planeWeightOffsets;
    vector<float> basisFunctionCenters;
    vector<cl_int> allIntParams;
    vector<float> allFloatParams;

    for (size_t i = 0 ; i < planeBasisLenses.size() ; i++)
    {
        planeIntParamOffsets.push_back(intParamCount);
        planeFloatParamOffsets.push_back(floatParamCount);
        planeWeightOffsets.push_back(numPlaneWeights);

        const GravitationalLens *pUnscaledLens = (haveUnscaledLenses)?unscaledLensesPerPlane[i].get():nullptr;

        string alphaFnName = "getAlpha_" + to_string(i);
        string generatedCode;

        if (!(r = getAlphaCodeForPlane(alphaFnName, planeBasisLenses[i], pUnscaledLens, subRoutineCode, intParamCount, floatParamCount, 
                                            numPlaneWeights, basisFunctionCenters, allIntParams, allFloatParams, generatedCode)))
            return "Can't get alpha code for plane: " + r.getErrorString();

        alphaCode += generatedCode;

        code << "    if (lpIdx == " << i << ")";
        code << "        return " << alphaFnName << "(theta, pAllIntParams + pPlaneIntParamOffsets[" << i << "], pAllFloatParams + pPlaneFloatParamOffsets[" << i << "],\n";
        code << "                                     pAllWeights + pPlaneWeightOffsets[" << i << "], pAllCenters + pPlaneWeightOffsets[" << i << "]*2, scalableFunctionScale);\n";
    }
    code << "    return (float2)(0.0f/0.0f, 0.0f/0.0f);\n";
    code << "}\n";

    m_numWeights = (size_t)numPlaneWeights;
    m_planeWeightOffsets.clear();
    for (auto x : planeWeightOffsets)
        m_planeWeightOffsets.push_back(x);
    m_planeWeightOffsets.push_back(m_numWeights); // allows us to use planeIdx+1 to refer to the start of the next plane

    allIntParams.push_back(-12345); // Add a sentinel and avoid a length zero array
    allFloatParams.push_back(-12345);

    // cout << "Total number of weights: " << m_numWeights << endl;
    // cout << "Plane weight offsets:" << endl;
    // for (auto o : planeWeightOffsets)
    //     cout << "    " << o << endl;
    // cout << "Plane int offsets:" << endl;
    // for (auto o : planeIntParamOffsets)
    //     cout << "    " << o << endl;
    // cout << "Plane float offsets:" << endl;
    // for (auto o : planeFloatParamOffsets)
    //     cout << "    " << o << endl;
    // cout << "Basis function centers: " << endl;
    // for (size_t i = 0 ; i < basisFunctionCenters.size() ; i += 2)
    //     cout << "    " << basisFunctionCenters[i] << "," << basisFunctionCenters[i+1] << endl;
    // cout << "All integer parameters:" << endl;
    // for (auto p : allIntParams)
    //     cout << "    " << p << endl;
    // cout << "All float parameters:" << endl;
    // for (auto p : allFloatParams)
    //     cout << "    " << p << endl;

    auto uploadOffsets = [this](const vector<cl_int> &offsets, const string &msg, cl_mem &dest) -> bool_t
    {
        cl_int err;
        dest = clCreateBuffer(getContext(), CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*offsets.size(), (void*)offsets.data(), &err);
        if (err != CL_SUCCESS || dest == nullptr)
            return "Error uploading " + msg + "to GPU";
        return true;
    };

    if (!(r = uploadOffsets(planeWeightOffsets, "plane weight offsets", m_common.m_pDevPlaneWeightOffsets)) ||
        !(r = uploadOffsets(planeIntParamOffsets, "plane int param offsets", m_common.m_pDevPlaneIntParamOffsets)) ||
        !(r = uploadOffsets(planeFloatParamOffsets, "plane float param offsets", m_common.m_pDevPlaneFloatParamOffsets)) )
        return r;
    
    cl_int err;
    m_common.m_pDevCenters = clCreateBuffer(getContext(), CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float)*basisFunctionCenters.size(), basisFunctionCenters.data(), &err);
    if (err != CL_SUCCESS || m_common.m_pDevCenters == nullptr)
        return "Error uploading basis function centers to GPU";
    m_common.m_pDevAllIntParams = clCreateBuffer(getContext(), CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(cl_int)*allIntParams.size(), allIntParams.data(), &err);
    if (err != CL_SUCCESS || m_common.m_pDevAllIntParams == nullptr)
        return "Error uploading integer parameters to GPU";
    m_common.m_pDevAllFloatParams = clCreateBuffer(getContext(), CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float)*allFloatParams.size(), allFloatParams.data(), &err);
    if (err != CL_SUCCESS || m_common.m_pDevAllFloatParams == nullptr)
        return "Error uploading float parameters to GPU";
    
    code << R"XYZ(
float2 multiPlaneTrace(float2 theta, int numPlanes, __global const float *Dsrc, __global const float *Dmatrix,
            __global const int *pAllIntParams, __global const float *pAllFloatParams, __global const float *pAllWeights,
            __global const float *pAllCenters, __global const int *pPlaneIntParamOffsets,
            __global const int *pPlaneFloatParamOffsets, __global const int *pPlaneWeightOffsets,
            float scalableFunctionScale)
{
    float2 T[MAXPLANES];
    float2 alphas[MAXPLANES];

    for (int j = 0 ; j < numPlanes ; j++)
    {
        T[j] = theta;
        for (int i = 0 ; i < j-1 ; i++)
            T[j] -= (Dmatrix[(i+1)*MAXPLANES + j]/Dmatrix[0 + j]) * alphas[i];
        const int i = j-1;
        if (i >= 0)
        {
            const float2 alpha = getAlpha(i, T[i], pAllIntParams, pAllFloatParams, pAllWeights, pAllCenters, pPlaneIntParamOffsets, pPlaneFloatParamOffsets, pPlaneWeightOffsets, scalableFunctionScale);
            alphas[i] = alpha;
            const float dfrac = Dmatrix[(i+1)*MAXPLANES + j]/Dmatrix[0 + j];

            T[j] -= dfrac * alpha;
        }
    }

    float2 beta = theta;
    for (int i = 0 ; i < numPlanes-1 ; i++)
        beta -= (Dsrc[i+1]/Dsrc[0]) * alphas[i];

    const int i = numPlanes-1;
    if (i >= 0)
    {
        const float2 alpha = getAlpha(i, T[i], pAllIntParams, pAllFloatParams, pAllWeights, pAllCenters, pPlaneIntParamOffsets, pPlaneFloatParamOffsets, pPlaneWeightOffsets, scalableFunctionScale);
        const float dfrac = Dsrc[i+1]/Dsrc[0];
        beta -= dfrac * alpha;
    }
    return beta;
}
)XYZ";

    string allCode = planeBasisLenses[0][0]->m_pLens->getCLLensQuantitiesStructure();
    
    // Add the code for all lens models
    for (auto &it : subRoutineCode)
        allCode += it.second;

    // If there's a composite lens somewhere, we need to add code for this as well
    if (compLensSubNames.size() > 0)
    {
        string fnName;
        allCode += CompositeLens::getCLProgram(fnName, compLensSubNames, compLensRecursion, false, false);
    }

    allCode += alphaCode;
    allCode += code.str();
    resultingCode = allCode;

    return true;
}

bool_t OpenCLCalculator::setupBasisFunctions(const vector<vector<shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
                            const vector<shared_ptr<GravitationalLens>> &unscaledLensesPerPlane)
{
    bool_t r;
    string oclTraceCode;
    
    // This also uploads parameters and centers (and some offsets)
    if (!(r = getMultiPlaneTraceCode(planeBasisLenses, unscaledLensesPerPlane, oclTraceCode)))
        return "Unable to het multiplane OpenCL code: " + r.getErrorString();

    string kernel = oclTraceCode + R"XYZ(
__kernel void calculateBetas(const int numPoints, const int numScaleFactors, const int numGenomesInOutput, const int numGenomeWeights,
                            __global const float *pThetas, __global float *pBetas, 
                            __global const int *pNumPlanes, __global const float *DsrcAll, __global const float *Dmatrix,
                            __global const int *pAllIntParams, __global const float *pAllFloatParams,
                            __global const float *pAllGenomeWeights, __global const float *pAllCenters,
                            __global const int *pPlaneIntParamOffsets,
                            __global const int *pPlaneFloatParamOffsets, __global const int *pPlaneWeightOffsets,
                            __global const float *pScalableFunctionScales, __global const int *pWeightIndices)
{
    const int i = get_global_id(0);
    if (i >= numPoints)
        return;

    const int outputIdx = get_global_id(1);
    if (outputIdx >= numGenomesInOutput)
        return;

    const int genomeWeightIdx = pWeightIndices[outputIdx];

    const int j = get_global_id(2);
    if (j >= numScaleFactors)
        return;

    const int numPlanesForPoint = pNumPlanes[i];
    __global const float *pAllWeights = pAllGenomeWeights + genomeWeightIdx*numGenomeWeights;

    // For each point, a number of scales will be used
    const float2 theta = (float2)(pThetas[i*2+0], pThetas[i*2+1]);

    // For these weights (a genome), and for this point, we'll calculate the results for the requested scale factor
    const float scalableFunctionScale = pScalableFunctionScales[outputIdx*numScaleFactors + j];

    // Each Dsrc is vector of MAXPLANES+1 length
    __global const float *Dsrc = DsrcAll + (MAXPLANES+1)*i;
    const float2 beta = multiPlaneTrace(theta, numPlanesForPoint, Dsrc, Dmatrix,
                                        pAllIntParams, pAllFloatParams, pAllWeights, pAllCenters,
                                        pPlaneIntParamOffsets, pPlaneFloatParamOffsets, pPlaneWeightOffsets,
                                        scalableFunctionScale);

    const int offset = (outputIdx*numScaleFactors*numPoints + j*numPoints + i)*2;
    pBetas[offset + 0] = beta.x;
    pBetas[offset + 1] = beta.y;

    //TODO: for testing    
    //pBetas[offset + 0] = genomeWeightIdx;
    //pBetas[offset + 1] = j;
}
)XYZ";

    //cout << kernel << endl;

    string faillog;
    if (!loadKernel(kernel, "calculateBetas", faillog))
    {
        cerr << faillog << endl;
        return "Error compiling kernel: " + getErrorString();
    }

    return true;
}

bool_t OpenCLCalculator::initGPU(int devIdx)
{
    string library = getLibraryName();
    if (!loadLibrary(library))
        return "Can't load OpenCL library: " + getErrorString();

	if (devIdx < 0) // Means rotate over the available devices
	{
		string fileName = "/dev/shm/grale_mpopencl_nextdevice.dat";
		getenv("GRALE_OPENCL_AUTODEVICEFILE", fileName); // Doesn't change file name if envvar not set

		m_perNodeCounter = make_unique<PerNodeCounter>(fileName);

		int idx = m_perNodeCounter->getCount();
		if (idx < 0)
		{
			m_perNodeCounter.reset();
			return "Couldn't read per-node device index from file '" + fileName + "': " + m_perNodeCounter->getErrorString();
		}

		int numDevices = getDeviceCount();
		if (numDevices < 0)
			return "Error getting device count: " + getErrorString();
		if (numDevices == 0)
			return "Unexpectedly got zero GPU devices";

		devIdx = idx%numDevices;

		auto GetTimeStamp = []() {
			struct timeval tv;
			gettimeofday(&tv,NULL);
			return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
		};

		cerr << "DEBUG (" << GetTimeStamp() << "): Automatically using new device index " << devIdx << endl;
	}

    if (!OpenCLKernel::init(devIdx))
        return "Can't init specified GPU device: " + getErrorString();
    return true;
}

bool_t OpenCLCalculator::setupMultiPlaneDistanceMatrix(const Cosmology &cosm, const vector<float> &zds)
{
    size_t numPlanes = zds.size();
    size_t cols = numPlanes;
    size_t rows = numPlanes+1;
    
    vector<float> Dij(rows*cols, numeric_limits<float>::quiet_NaN());
    for (size_t j = 0 ; j < numPlanes ; j++)
        Dij[0 + j] = (float)(cosm.getAngularDiameterDistance(zds[j])/DIST_MPC);

    for (size_t i = 1 ; i <= numPlanes ; i++)
    {
        for (size_t j = i ; j < numPlanes ; j++)
            Dij[i*cols + j] = (float)(cosm.getAngularDiameterDistance(zds[i-1], zds[j])/DIST_MPC);
    }
    
    cl_int err;
    m_common.m_pDevDmatrix = clCreateBuffer(getContext(), CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float)*Dij.size(), Dij.data(), &err);
    if (err != CL_SUCCESS || m_common.m_pDevDmatrix == nullptr)
        return "Error uploading lens plane distance matrix to GPU";

    // cout << "Distance matrix uploaded: " << endl;
    // for (size_t i = 0 ; i <= numPlanes ; i++)
    // {
    //     for (size_t j = 0 ; j < numPlanes ; j++)
    //         cout << "\t" << Dij[cols*i+j];
    //     cout << endl;
    // }

    m_numPlanes = numPlanes;
    return true;
}

bool_t OpenCLCalculator::setupAngularDiameterDistances(const Cosmology &cosm, const vector<float> &zds,
                                        const vector<ImagesDataExtended*> &images,
                                        cl_mem &pDevDpoints, cl_mem &pDevUsedPlanes)
{
    // Build the image points redshift vector
    // TODO: For now, we're associating a zs to every single point, not just to the entire source
    vector<float> zss;
    for (auto img : images)
    {
        double zs = numeric_limits<double>::quiet_NaN();
        img->getExtraParameter("z", zs);

        int numImages = img->getNumberOfImages();
        for (int i = 0 ; i < numImages ; i++)
        {
            int numPoints = img->getNumberOfImagePoints(i);
            for (int p = 0 ; p < numPoints ; p++)
                zss.push_back((float)zs);
        }
    }

    size_t numPlanes = zds.size();
    size_t numPoints = zss.size();

    size_t numCols = numPlanes+1;
    size_t numRows = numPoints;
    vector<float> Dsources(numRows*numCols, numeric_limits<float>::quiet_NaN());
    vector<cl_int> usedPlanes(numPoints);

    for (size_t p = 0 ; p < numPoints ; p++)
    {
        float zs = zss[p];

        Dsources[p*numCols + 0] = (float)(cosm.getAngularDiameterDistance(zs)/DIST_MPC);
        size_t useCount = 0;
        for (size_t i = 0 ; i < numPlanes ; i++)
        {
            if (zs > zds[i])
            {
                useCount++;
                Dsources[p*numCols + (i+1)] = (float)(cosm.getAngularDiameterDistance(zds[i], zs)/DIST_MPC);
            }
        }
        usedPlanes[p] = useCount;
    }

    // cout << "Point distance info:" << endl;
    // for (size_t p = 0 ; p < numPoints ; p++)
    // {
    //     cout << p << ") " << usedPlanes[p] << "|\t";
    //     for (size_t i = 0 ; i < numCols ; i++)
    //         cout << "\t" << Dsources[p*numCols + i];
    //     cout << endl;
    // }
    
    // Upload info to GPU
    cl_int err;
    pDevDpoints = clCreateBuffer(getContext(), CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float)*Dsources.size(), Dsources.data(), &err);
    if (err != CL_SUCCESS || pDevDpoints == nullptr)
        return "Error uploading distances from lens planes to GPU";

    pDevUsedPlanes = clCreateBuffer(getContext(), CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float)*usedPlanes.size(), usedPlanes.data(), &err);
    if (err != CL_SUCCESS || pDevUsedPlanes == nullptr)
        return "Error uploading number of used lens planes to GPU";
    
    return true;
}

bool_t OpenCLCalculator::allocateAndUploadImages(const vector<ImagesDataExtended *> &images, cl_mem &pDevImages, size_t &numPoints)
{
    vector<float> allCoordinates;
    for (auto img : images)
    {
        int numImages = img->getNumberOfImages();
        for (int i = 0 ; i < numImages ; i++)
        {
            int numPoints = img->getNumberOfImagePoints(i);
            for (int p = 0 ; p < numPoints ; p++)
            {
                Vector2Dd pos = img->getImagePointPosition(i, p);
                pos /= m_angularScale;
                allCoordinates.push_back((float)pos.getX());
                allCoordinates.push_back((float)pos.getY());
            }
        }
    }

    cl_int err;
    pDevImages = clCreateBuffer(getContext(), CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float)*allCoordinates.size(), allCoordinates.data(), &err);
    if (err != CL_SUCCESS || pDevImages == nullptr)
        return "Error uploading image data to GPU";
    
    numPoints = allCoordinates.size()/2;
    // cout << "Uploaded coordinates for " << numPoints << " points to GPU" << endl;

    return true;
}

bool_t OpenCLCalculator::setupImages(const vector<ImagesDataExtended *> &allImages,
                    const vector<ImagesDataExtended *> &shortImages)
{
    bool_t r;
    if (!(r = allocateAndUploadImages(allImages, m_full.m_pDevImages, m_full.m_numPoints)))
        return r;

    if (shortImages.size() == 0) // indicates that all images should be used
    {
        m_shortImagesAreAllImages = true;
        m_short.m_pDevImages = m_full.m_pDevImages;
        m_short.m_numPoints = m_full.m_numPoints;
    }
    else
    {
        m_shortImagesAreAllImages = false;
        if (!(r = allocateAndUploadImages(shortImages, m_short.m_pDevImages, m_short.m_numPoints)))
            return r;
    }

    return true;
}

void OpenCLCalculator::CLMem::dealloc(OpenCLKernel &cl)
{
    if (!m_pMem)
        return;
    cl.clReleaseMemObject(m_pMem);
    m_pMem = nullptr;
    m_size = 0;
}

bool_t OpenCLCalculator::CLMem::realloc(OpenCLKernel &cl, size_t s) // Only reallocates if more memory is requested
{
    if (s <= m_size)
        return true;

    cl_int err;

    cl_mem pBuf = cl.clCreateBuffer(cl.getContext(), CL_MEM_READ_WRITE, s, nullptr, &err);
    if (pBuf == nullptr || err != CL_SUCCESS)
        return "Can't create buffer of size " + to_string(s) + " on GPU: code " + to_string(err);

    dealloc(cl);
    m_pMem = pBuf;
    m_size = s;
    return true;
}

bool_t OpenCLCalculator::CLMem::enqueueWriteBuffer(OpenCLKernel &cl, const void *pData, size_t s)
{
    if (s > m_size)
        return "Size exceeds GPU buffer size";
    cl_int err = cl.clEnqueueWriteBuffer(cl.getCommandQueue(), m_pMem, false, 0, s, pData, 0, nullptr, nullptr);
    if (err != CL_SUCCESS)
        return "Error enqueuing write buffer: code " + to_string(err);
    return true;
}

bool_t OpenCLCalculator::CLMem::enqueueReadBuffer(OpenCLKernel &cl, void *pData, size_t s, cl_event *pEvt)
{
    if (s > m_size)
        return "Reading beyond CPU buffer size";
    cl_int err = cl.clEnqueueReadBuffer(cl.getCommandQueue(), m_pMem, false, 0, s, pData, 0, nullptr, pEvt);
    if (err != CL_SUCCESS)
        return "Error enqueuing read buffer: code " + to_string(err);
    return true;
}

void OpenCLCalculator::CommonClMem::dealloc(OpenCLKernel &cl)
{
    auto release = [&cl](cl_mem x)
    {
        if (x)
            cl.clReleaseMemObject(x);
    };
    release(m_pDevDmatrix);
    release(m_pDevPlaneWeightOffsets);
    release(m_pDevPlaneIntParamOffsets);
    release(m_pDevPlaneFloatParamOffsets);
    release(m_pDevCenters);
    release(m_pDevAllIntParams);
    release(m_pDevAllFloatParams);
    zeroAll();
}

void OpenCLCalculator::FullOrShortClMem::dealloc(OpenCLKernel &cl)
{
    auto release = [&cl](cl_mem x)
    {
        if (x)
            cl.clReleaseMemObject(x);
    };
    release(m_pDevImages);
    release(m_pDevDpoints);
    release(m_pDevUsedPlanes);
    zeroAll();
}

} // end namespace