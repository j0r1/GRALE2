#pragma once

#include "graleconfig.h"
#include "openclkernel.h"
#include "imagesdataextended.h"
#include "cosmology.h"
#include "lensinversionbasislensinfo.h"
#include "lensgaindividual.h"
#include "pernodecounter.h"
#include "oclutils.h"
#include <errut/booltype.h>
#include <vector>
#include <memory>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <chrono>

namespace grale
{

class OpenCLMultiPlaneCalculator : public OpenCLKernel
{
public:
	static errut::bool_t initInstance(int devIdx, const std::vector<ImagesDataExtended *> &allImages,
										   const std::vector<ImagesDataExtended *> &shortImages,
										   const std::vector<float> &zds,
										   const Cosmology &cosm,
										   const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
										   const std::vector<std::shared_ptr<GravitationalLens>> &unscaledLensesPerPlane,
										   uint64_t userId,
										   const std::vector<std::shared_ptr<GravitationalLens>> &baseLensesPerPlane
										   );
	static void releaseInstance(uint64_t userId);
	static OpenCLMultiPlaneCalculator &instance();

	void setGenomesToCalculate(size_t s);
	errut::bool_t startNewBackprojection(const LensGAGenome &g);
	errut::bool_t scheduleUploadAndCalculation(const LensGAGenome &g, int &calculationIdentifier,
											   const std::vector<std::pair<float,float>> &scaleFactors, bool useShort);
	
	errut::bool_t isCalculationDone(const LensGAGenome &g, int calculationIdentifier, size_t *pGenomeIndex,
									const float **ppAllBetas, size_t *pNumGenomes, size_t *pNumSteps, size_t *pNumPoints, bool *pDone);
	errut::bool_t setCalculationProcessed(const LensGAGenome &g, int calculationIdentifier);

	double getAngularScale() const { return m_angularScale; }

	OpenCLMultiPlaneCalculator();
	~OpenCLMultiPlaneCalculator();
private:
	int getDeviceIndex() const { return m_devIdx; }
	errut::bool_t initAll(int devIdx, const std::vector<ImagesDataExtended *> &allImages,
							const std::vector<ImagesDataExtended *> &shortImages,
							const std::vector<float> &zds,
							const Cosmology &cosm,
							const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
							const std::vector<std::shared_ptr<GravitationalLens>> &unscaledLensesPerPlane,
							const std::vector<std::shared_ptr<GravitationalLens>> &baseLensesPerPlane);
	errut::bool_t analyzeCompositeLenses(const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
								  const std::vector<std::shared_ptr<GravitationalLens>> &unscaledLensesPerPlane,
								  const std::vector<std::shared_ptr<GravitationalLens>> &baseLensesPerPlane,
								  std::map<std::string,std::string> &subRoutineCode,
								  std::vector<std::string> &compLensSubRoutineNames,
								  int &compLensRecursion
								  );
	errut::bool_t getAlphaCodeForPlane(const std::string &functionName,
								const std::vector<std::shared_ptr<LensInversionBasisLensInfo>> &basisLenses,
								const GravitationalLens *pUnscaledLens,
								const GravitationalLens *pBaseLens,
								std::map<std::string,std::string> &subRoutineCode,
								int &intParamCount, int &floatParamCount, int &numWeights,
								std::vector<float> &centers,
								std::vector<cl_int> &intParams,
								std::vector<float> &floatParams,
								std::string &generatedCode);
	errut::bool_t getMultiPlaneTraceCode(const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
								  const std::vector<std::shared_ptr<GravitationalLens>> &unscaledLensesPerPlane,
								  const std::vector<std::shared_ptr<GravitationalLens>> &baseLensesPerPlane,
								  std::string &resultingCode);
	errut::bool_t setupBasisFunctions(const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
							   const std::vector<std::shared_ptr<GravitationalLens>> &unscaledLensesPerPlane,
							   const std::vector<std::shared_ptr<GravitationalLens>> &baseLensesPerPlane);
	errut::bool_t initGPU(int devIdx);
	errut::bool_t setupMultiPlaneDistanceMatrix(const Cosmology &cosm, const std::vector<float> &zds);
	errut::bool_t setupAngularDiameterDistances(const Cosmology &cosm, const std::vector<float> &zds,
										 const std::vector<ImagesDataExtended*> &images,
										 cl_mem &pDevDpoints, cl_mem &pDevUsedPlanes);
	errut::bool_t allocateAndUploadImages(const std::vector<ImagesDataExtended *> &images, cl_mem &pDevImages, size_t &numPoints);
	errut::bool_t setupImages(const std::vector<ImagesDataExtended *> &allImages,
					   const std::vector<ImagesDataExtended *> &shortImages);

	errut::bool_t checkCalculateScheduledContext();
	errut::bool_t checkAddFakeSource(std::vector<ImagesDataExtended *> &allImages, std::vector<ImagesDataExtended *> &shortImages);


	// static void staticEventNotify(cl_event event, cl_int event_command_status, void *user_data);
	// void eventNotify(cl_event event, cl_int event_command_status);

	oclutils::CLMem m_debug;

	class State
	{
	public:
		State(size_t genomeIndex = std::numeric_limits<size_t>::max()) : m_genomeIndex(genomeIndex) { }
		size_t m_genomeIndex;
	};

	struct CommonClMem
	{
		CommonClMem()
		{
			zeroAll();
		}

		void zeroAll()
		{
			m_pDevDmatrix = nullptr;
			m_pDevPlaneWeightOffsets = nullptr;
			m_pDevPlaneIntParamOffsets = nullptr;
			m_pDevPlaneFloatParamOffsets = nullptr;
			m_pDevCenters = nullptr;
			m_pDevAllIntParams = nullptr;
			m_pDevAllFloatParams = nullptr;
		}

		void dealloc(OpenCLKernel &cl);

		cl_mem m_pDevDmatrix;
		cl_mem m_pDevPlaneWeightOffsets, m_pDevPlaneIntParamOffsets, m_pDevPlaneFloatParamOffsets;
		cl_mem m_pDevCenters, m_pDevAllIntParams, m_pDevAllFloatParams;
		oclutils::CLMem m_devAllWeights;
	};

	struct FullOrShortClMem
	{
		FullOrShortClMem()
		{
			zeroAll();
		}

		void zeroAll()
		{
			m_pDevImages = nullptr;
			m_pDevDpoints = nullptr;
			m_pDevUsedPlanes = nullptr;
		}

		void dealloc(OpenCLKernel &cl);
		
		cl_mem m_pDevImages;
		cl_mem m_pDevDpoints;
		cl_mem m_pDevUsedPlanes;
		size_t m_numPoints;
	};

	class ContextMemory
	{
	public:
		void resizeCPUBuffers()
		{
			m_genomeIndexForBetaIndex.resize(0);
			m_allFactors.resize(0);
			m_allBetas.resize(0);
		}
		std::vector<cl_int> m_genomeIndexForBetaIndex;
		std::vector<float> m_allFactors, m_allBetas;
		oclutils::CLMem m_devFactors, m_devBetas, m_devGenomeIndexForBetaIndex;
	};

	std::vector<std::unique_ptr<ContextMemory>> m_recycledContextMemory;
	void cleanupContextMemory(ContextMemory &mem);

	class CalculationContext
	{
	public:
		CalculationContext(int ident, size_t numSteps, size_t numWeights, const CommonClMem &common,
						   const FullOrShortClMem &fullOrShort, std::unique_ptr<ContextMemory> ctxMem)
			: m_identifier(ident), m_numSteps(numSteps), m_numWeights(numWeights),
			  m_common(common), m_fullOrShort(fullOrShort), m_calculated(false), m_calcEvt(nullptr), m_evt(nullptr),
			  m_mem(std::move(ctxMem))
			{ }

		errut::bool_t schedule(OpenCLKernel &cl, const LensGAGenome &g, size_t genomeWeightsIndex,
							   const std::vector<std::pair<float,float>> &scaleFactors);
		errut::bool_t calculate(OpenCLKernel &cl, cl_event *pEvt);

		const int m_identifier;
		const size_t m_numSteps, m_numWeights;
		const CommonClMem &m_common;
		const FullOrShortClMem &m_fullOrShort;
		bool m_calculated;
		std::chrono::time_point<std::chrono::steady_clock> m_calcQueueTime;

		cl_event m_calcEvt, m_evt;

		std::unordered_map<const LensGAGenome *, size_t> m_betaIndexForGenome;
		std::unique_ptr<ContextMemory> m_mem;
	};

	int m_devIdx;
	bool m_errorState;

	std::vector<std::unique_ptr<ImagesDataExtended>> m_fakeImages;

	std::unique_ptr<CalculationContext> m_beingScheduled;
	std::unique_ptr<CalculationContext> m_beingCalculated;
	std::unique_ptr<CalculationContext> m_doneCalculating;
	int m_nextCalculationContextIdentifier;

	size_t m_numGenomesToCalculate;
	std::unordered_map<const LensGAGenome *, State> m_states;
	std::mutex m_mutex;
	std::vector<float> m_allBasisFunctionWeights;
	size_t m_nextGenomeIndex;
	bool m_hasCalculated;
	bool m_uploadedWeights;

	double m_angularScale = 0, m_potentialScale = 0;

	bool m_shortImagesAreAllImages;
	size_t m_numWeights, m_numPlanes;
	std::vector<size_t> m_planeWeightOffsets;

	CommonClMem m_common;
	FullOrShortClMem m_full, m_short;

	std::unique_ptr<PerNodeCounter> m_perNodeCounter;

	uint64_t m_timeoutCheckCounter;

	static std::unique_ptr<OpenCLMultiPlaneCalculator> s_pInstance;
	static std::mutex s_instanceMutex;
	static std::set<uint64_t> s_users;
	static bool s_initTried;
};

} // end namespace
