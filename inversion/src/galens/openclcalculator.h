#pragma once

#include "graleconfig.h"
#include "openclkernel.h"
#include "imagesdataextended.h"
#include "cosmology.h"
#include "lensinversionbasislensinfo.h"
#include "lensgaindividual.h"
#include <errut/booltype.h>
#include <vector>
#include <memory>
#include <mutex>
#include <map>

namespace grale
{

class OpenCLCalculator : public OpenCLKernel
{
public:
    static errut::bool_t initInstance(int devIdx, const std::vector<ImagesDataExtended *> &allImages,
	                                       const std::vector<ImagesDataExtended *> &shortImages,
										   const std::vector<float> &zds,
										   const Cosmology &cosm,
										   const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
										   const std::vector<std::shared_ptr<GravitationalLens>> &unscaledLensesPerPlane
										   );
	static OpenCLCalculator &instance();

	void setGenomesToCalculate(size_t s);
	errut::bool_t startNewBackprojection(const LensGAGenome &g, size_t &genomeIndex);

	enum DeviceStatus { Ready, UploadingWeights, Calculating, CalculationDone, Unknown };
    bool isCalculationDone();

	errut::bool_t scheduleUploadAndCalculation(size_t genomeIndex, const std::vector<std::pair<float,float>> &scaleFactors, bool useShort);
	errut::bool_t getGenomeIndex(const LensGAGenome &g, size_t &genomeIndex) const;

	OpenCLCalculator();
	~OpenCLCalculator();
private:
	errut::bool_t initAll(int devIdx, const std::vector<ImagesDataExtended *> &allImages,
	                        const std::vector<ImagesDataExtended *> &shortImages,
							const std::vector<float> &zds,
							const Cosmology &cosm,
							const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
							const std::vector<std::shared_ptr<GravitationalLens>> &unscaledLensesPerPlane);
	errut::bool_t analyzeCompositeLenses(const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
							      const std::vector<std::shared_ptr<GravitationalLens>> &unscaledLensesPerPlane,
								  std::map<std::string,std::string> &subRoutineCode,
								  std::vector<std::string> &compLensSubRoutineNames,
								  int &compLensRecursion
								  );
	errut::bool_t getAlphaCodeForPlane(const std::string &functionName,
		                        const std::vector<std::shared_ptr<LensInversionBasisLensInfo>> &basisLenses,
								const GravitationalLens *pUnscaledLens,
								std::map<std::string,std::string> &subRoutineCode,
								int &intParamCount, int &floatParamCount, int &numWeights,
								std::vector<float> &centers,
								std::vector<cl_int> &intParams,
								std::vector<float> &floatParams,
								std::string &generatedCode);
	errut::bool_t getMultiPlaneTraceCode(const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
							      const std::vector<std::shared_ptr<GravitationalLens>> &unscaledLensesPerPlane,
								  std::string &resultingCode);
	errut::bool_t setupBasisFunctions(const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &planeBasisLenses,
							   const std::vector<std::shared_ptr<GravitationalLens>> &unscaledLensesPerPlane);
	errut::bool_t initGPU(int devIdx);
	errut::bool_t setupMultiPlaneDistanceMatrix(const Cosmology &cosm, const std::vector<float> &zds);
	errut::bool_t setupAngularDiameterDistances(const Cosmology &cosm, const std::vector<float> &zds,
	                                     const std::vector<ImagesDataExtended*> &images,
										 cl_mem &pDevDpoints, cl_mem &pDevUsedPlanes);
	errut::bool_t allocateAndUploadImages(const std::vector<ImagesDataExtended *> &images, cl_mem &pDevImages, size_t &numPoints);
	errut::bool_t setupImages(const std::vector<ImagesDataExtended *> &allImages,
	                   const std::vector<ImagesDataExtended *> &shortImages);

    errut::bool_t uploadScaleFactorsAndBackproject(bool useShort, int numScaleFactors);

    static void staticEventNotify(cl_event event, cl_int event_command_status, void *user_data);
    void eventNotify(cl_event event, cl_int event_command_status);

	class CLMem
	{
	public:
		CLMem() : m_pMem(nullptr), m_size(0) { }
		void dealloc(OpenCLKernel &cl);
		errut::bool_t realloc(OpenCLKernel &cl, size_t s); // Only reallocates if more memory is requested
		template<class T> errut::bool_t realloc(OpenCLKernel &cl, const std::vector<T> &buffer)	{ return realloc(cl, buffer.size()*sizeof(T)); }

		errut::bool_t enqueueWriteBuffer(OpenCLKernel &cl, const void *pData, size_t s);
		template<class T> errut::bool_t enqueueWriteBuffer(OpenCLKernel &cl, const std::vector<T> &data) { return enqueueWriteBuffer(cl, data.data(), data.size()*sizeof(T)); }

        errut::bool_t enqueueReadBuffer(OpenCLKernel &cl, void *pData, size_t s, cl_event *pEvt);
        template<class T> errut::bool_t enqueueReadBuffer(OpenCLKernel &cl, std::vector<T> &data, cl_event *pEvt) { return enqueueReadBuffer(cl, data.data(), data.size()*sizeof(T), pEvt); }

		cl_mem m_pMem;
		size_t m_size;
	};

	class State // TODO: do we need more than a genome index?
	{
	public:
		State(size_t genomeIndex) : m_genomeIndex(genomeIndex) { }
		const size_t m_genomeIndex;
	};

	std::map<const LensGAGenome *, std::unique_ptr<State>> m_states; // TODO: probably don't need a pointer here
	size_t m_nextGenomeIndex;
	bool m_hasCalculated;
	std::mutex m_mutex;
	std::vector<float> m_allBasisFunctionWeights;
	std::vector<float> m_allFactors;
	size_t m_uploadInfoCount;
	DeviceStatus m_devStatus;

	double m_angularScale, m_potentialScale;

	bool m_shortImagesAreAllImages;
	cl_mem m_pDevAllImages, m_pDevShortImages;
	size_t m_numAllImagePoints, m_numShortImagePoints;
	size_t m_numWeights;

	cl_mem m_pDevDmatrix;
	cl_mem m_pDevDpointsAll, m_pDevDpointsShort;
	cl_mem m_pDevUsedPlanesAll, m_pDevUsedPlanesShort;
	cl_mem m_pDevPlaneWeightOffsets, m_pDevPlaneIntParamOffsets, m_pDevPlaneFloatParamOffsets;
	cl_mem m_pDevCenters, m_pDevAllIntParams, m_pDevAllFloatParams;

	CLMem m_devAllWeights, m_devBetas, m_devFactors;
	std::vector<float> m_allBetas;

	size_t m_genomesLeftToCalculate;

	static std::unique_ptr<OpenCLCalculator> s_pInstance;
	static std::mutex s_instanceMutex;
	static bool s_initTried;
};

} // end namespace