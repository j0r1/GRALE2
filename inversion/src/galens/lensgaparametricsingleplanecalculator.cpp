#include "lensgaparametricsingleplanecalculator.h"
#include "openclsingleplanedeflection.h"
#include "utils.h"
#include "constants.h"
#include <eatk/vectorgenomefitness.h>
#include <eatk/valuefitness.h>
#include <iostream>
#include <iomanip>
#include <list>
#include <cassert>
#include <sstream>
#include <mutex>
#include <fstream>

using namespace errut;
using namespace std;

namespace grale
{

class AccumulatedBetaStats
{
public:
	AccumulatedBetaStats(double deflScale) : m_deflScale(deflScale) { }

	void accumulate(const vector<float> &limits, const vector<size_t> &counts, size_t otherCounts, size_t penaltyCount)
	{
		if (m_limits.size() == 0)
		{
			assert(m_counts.size() == 0);
			assert(m_otherCounts == 0);
			m_limits = limits;
			m_counts = counts;
			m_otherCounts = otherCounts;
			m_penaltyCount = penaltyCount;
		}
		else
		{
			assert(m_limits.size() == limits.size());
			assert(m_counts.size() == counts.size());
			for (size_t i = 0 ; i < counts.size() ; i++)
				m_counts[i] += counts[i];
			m_otherCounts += otherCounts;
			m_penaltyCount += penaltyCount;
		}
	}

	~AccumulatedBetaStats()
	{
		stringstream ss;

		ss << "DEBUG: Beta Size Stats: ";
		for (size_t i = 0 ; i < m_limits.size() ; i++)
			ss << m_counts[i] << " < " << std::setprecision(8) << std::fixed << ((m_limits[i]*m_deflScale)/ANGLE_ARCSEC) << " | ";
		ss << " other: " << m_otherCounts << " ; penalty count = " << m_penaltyCount << endl;
		cerr << ss.str();
	}

	double m_deflScale = 0;
	size_t m_regCount = 0;

	vector<float> m_limits;
	vector<size_t> m_counts;
	size_t m_otherCounts = 0;
	size_t m_penaltyCount = 0;
};

class BetaSizeStats
{
public:
	BetaSizeStats(double deflScale)
	{
		for (auto s : vector<double> { 0.000001, 
				                       0.00001,
									   0.0001,
									   0.001,
									   0.01,
									   0.1 })
		{
			s *= ANGLE_ARCSEC;
			s /= deflScale;

			m_limits.push_back(s);
			m_counts.push_back(0);
		}

		lock_guard<mutex> guard(s_statsMutex);
		if (!s_accStats)
			s_accStats = make_unique<AccumulatedBetaStats>(deflScale);
		s_accStats->m_regCount++;
	}

	~BetaSizeStats()
	{
		lock_guard<mutex> guard(s_statsMutex);
		s_accStats->accumulate(m_limits, m_counts, m_otherCounts, m_penaltyCount);
		s_accStats->m_regCount--;
		if (s_accStats->m_regCount == 0)
			s_accStats = nullptr;
	}

	void count(float v)
	{
		for (size_t i = 0 ; i < m_limits.size() ; i++)
		{
			if (v < m_limits[i])
			{
				m_counts[i]++;
				return;
			}
		}
		m_otherCounts++;
	}

	vector<float> m_limits;
	vector<size_t> m_counts;
	size_t m_otherCounts = 0;
	size_t m_penaltyCount = 0;

	static mutex s_statsMutex;
	static unique_ptr<AccumulatedBetaStats> s_accStats;
};

mutex BetaSizeStats::s_statsMutex;
unique_ptr<AccumulatedBetaStats> BetaSizeStats::s_accStats;

LensGAParametricSinglePlaneCalculator::LensGAParametricSinglePlaneCalculator(std::unique_ptr<LensFitnessObject> fitObj)
	: m_fitObj(move(fitObj))
{
}

LensGAParametricSinglePlaneCalculator::~LensGAParametricSinglePlaneCalculator()
{
	OpenCLSinglePlaneDeflectionInstance::releaseInstance((uint64_t)this);
}

bool_t LensGAParametricSinglePlaneCalculator::init(const LensInversionParametersBase &params0)
{
	if (!dynamic_cast<const LensInversionParametersParametricSinglePlane*>(&params0))
		return "Parameters are not of the correct type";
	const LensInversionParametersParametricSinglePlane &params = static_cast<const LensInversionParametersParametricSinglePlane &>(params0);

	if (m_init)
		return "Already initialized";

	if (!m_fitObj.get())
		return "No fitness object was set";

	if (params.getImages().size() == 0)
		return "No images were specified";

	size_t numChangeAbleParams = params.getOffsets().size();
	if (numChangeAbleParams == 0)
		return "No changeable parameters";

	m_angularScale = params.getDeflectionScale();
	m_potScale = params.getPotentialScale();
	m_changeableParamIdx.clear();
	for (auto i : params.getOffsets())
	{
		assert(i >= 0);
		m_changeableParamIdx.push_back((size_t)i);
	}
	m_devIdx = params.getDeviceIndex();
	m_templateLens = params.getTemplateLens().createCopy();

	bool_t r;
	vector<ImagesDataExtended *> imagesPtrs;
	list<ImagesDataExtended *> imagesPtrsList;
	for (auto &img : params.getImages())
	{
		imagesPtrs.push_back(img.get());
		imagesPtrsList.push_back(img.get());
	}
	
	m_oclBp = make_unique<OclCalculatedBackProjector>();
	if (!(r = m_oclBp->init(imagesPtrs, m_angularScale)))
		return "Can't init backprojector: " + r.getErrorString();

	m_thetas.clear();
	m_distFrac.clear();
	m_pointMap.clear();

	bool havePosUncerts = false;
	// Use this flag to see if we should randomize the inputs, the
	// position uncertainty alone will not suffice: that (in the
	// future) should be able to indicate uncertainty for a bayesian
	// approach, where the input locations are fixed and the uncertainty
	// is used in a log prob
	if (params.getRandomizeInputPositions())
	{
		// Check
		for (auto &img : params.getImages())
			if (img->hasProperty(ImagesData::PositionUncertainty))
				havePosUncerts = true;

		if (!havePosUncerts)
			return "Randomization of input positions requested, but no position uncertainties were specified";
	}

	vector<float> posUncertainties;
	vector<float> retraceDistanceFractions;
	const vector<shared_ptr<ImagesDataExtended>> &images = params.getImages();
	const vector<bool> &retraceImage = params.shouldRetraceImages();
	vector<pair<int, float>> recalcThetaInfo;

	if (images.size() != retraceImage.size())
		return "Different number of images and retrace flags";

	bool anyRetrace = false;

	m_tracedSourcesFlags = make_shared<vector<bool>>();
	m_tracedSourcesPoints = make_shared<vector<vector<Vector2Df>>>();
	m_tracedSourcesConvergedFlags = make_shared<vector<vector<int>>>();

	for (size_t srcIdx = 0 ; srcIdx < images.size() ; srcIdx++)
	{
		auto &img = images[srcIdx];
		bool retrace = retraceImage[srcIdx];
		anyRetrace |= retrace;

		m_tracedSourcesFlags->push_back(retrace);
		m_tracedSourcesPoints->push_back(vector<Vector2Df>());
		m_tracedSourcesConvergedFlags->push_back(vector<int>());

		int numImages = img->getNumberOfImages();
		double frac = img->getDds()/img->getDs();
		assert(!isnan(frac));

		for (int i = 0 ; i < numImages ; i++)
		{
			int numPoints = img->getNumberOfImagePoints(i);
			for (int p = 0 ; p < numPoints ; p++)
			{
				Vector2Dd pos = img->getImagePointPosition(i, p);
				pos /= m_angularScale;

				Vector2Df floatPos((float)pos.getX(), (float)pos.getY());

				float uncert = 0;
				if (img->hasProperty(ImagesData::PositionUncertainty))
					uncert = (float)(img->getImagePointProperty(ImagesData::PositionUncertainty, i, p) / m_angularScale);

				auto newPointCallback = [this, srcIdx, havePosUncerts, uncert, retrace, frac,
					                     &posUncertainties, &recalcThetaInfo](Vector2Df pt, bool forceUnique, size_t thetaIndex) -> bool_t
				{
					if (havePosUncerts)
						posUncertainties.push_back(uncert);

					pair<int, float> traceInf = {-1, numeric_limits<float>::quiet_NaN() };
					if (retrace)
					{
						traceInf = { srcIdx, frac };
						m_bpPointInfo.push_back({ srcIdx, m_tracedSourcesPoints->back().size() });
						m_tracedSourcesPoints->back().push_back({ numeric_limits<float>::quiet_NaN(), numeric_limits<float>::quiet_NaN() });
						m_tracedSourcesConvergedFlags->back().push_back(false);
					}
					recalcThetaInfo.push_back(traceInf);
					return true;
				};

				auto existingPointCallback = [havePosUncerts, retrace, uncert, &posUncertainties, &recalcThetaInfo](Vector2Df pt, size_t thetaIndex) -> bool_t
				{
					size_t ptIdx = thetaIndex;
					if (havePosUncerts)
					{
						assert(ptIdx < posUncertainties.size());
						if (posUncertainties[ptIdx] != uncert)
							return "Points with same positions have different uncertainties";
					}

					if (retrace)
						return "Can't retrace points that are already known!";

					assert(ptIdx < recalcThetaInfo.size());
					if (recalcThetaInfo[ptIdx].first >= 0)
						return "Different retrace settings for identical points";

					return true;
				};

				bool forceNew = false;
				if (retrace)
					forceNew = true;
				if (!(r = m_pointMap.addPoint(floatPos, forceNew, newPointCallback, existingPointCallback)))
					return r;
	
				m_distFrac.push_back(frac);
			}
		}
	}

	m_thetas = m_pointMap.getPoints();
	if (!anyRetrace) // Don't do extra calculations if no retracing is needed
	{
		recalcThetaInfo.clear();
		m_tracedSourcesPoints = nullptr;
		m_tracedSourcesFlags = nullptr;
		m_tracedSourcesConvergedFlags = nullptr;
		m_bpPointInfo.clear();
	}
	else
	{
		// Keep some stats for now
		m_stats = make_unique<BetaSizeStats>(m_angularScale);

		// Set threshold to accept retrace as converged
		m_betaThresHold = (float)(params.getSourcePlaneDistanceThreshold()/m_angularScale);
		if (m_betaThresHold < 0)
			return "Negative convergence threshold for retracing is not allowed";
	}

	cerr << "INFO: total points " << m_pointMap.getPointMapping().size() << " =? " << m_distFrac.size() << ", unique points = " << m_thetas.size() << endl;
	// cerr << "Point map is:" << endl;
	// for (size_t i = 0 ; i < m_pointMap.getPointMapping().size() ; i++)
	// 	cerr << "    " << i << "\t -> " << m_pointMap.getPointMapping()[i] << endl;

	int numIntParams = 0, numFloatParams = 0;
	if (!m_templateLens->getCLParameterCounts(&numIntParams, &numFloatParams))
		return "Can't get floating point parameters: " + m_templateLens->getErrorString();
	
	m_intParams.resize(numIntParams);
	m_floatParams.resize(numFloatParams);
	m_intParams.push_back(-12345);
	m_floatParams.push_back(-12345);
	if (!m_templateLens->getCLParameters(m_angularScale, m_potScale, m_intParams.data(), m_floatParams.data()))
		return "Can't get float or in parameters: " + m_templateLens->getErrorString();

	m_kernelCode = m_templateLens->getCLLensProgram(m_angularScale, m_potScale, m_kernelName);
	if (m_kernelCode.length() == 0)
		return "Couldn't get OpenCL kernel code: " + m_templateLens->getErrorString();

	uint64_t posUncertSeed = params.getInitialPositionUncertaintySeed();
	if (posUncertainties.size() > 0)
	{
		m_randomizingInputPositions = true;
		cerr << "INFO: enabling EXPERIMENTAL positional uncertainty with seed " << posUncertSeed << endl;
	}
	else
	{
		m_randomizingInputPositions = false;
		cerr << "INFO: NOT enabling EXPERIMENTAL positional uncertainties" << endl;
	}

	vector<pair<size_t, string>> originParameterMapping = params.getOriginParameterMapping();
	m_numOriginParams = params.getNumberOfOriginParameters();

	m_extraClPriorCode = params.getOpenCLPriorCode();

	const TraceParameters &origTraceParams = params.getRetraceParameters();
	unique_ptr<TraceParameters> retraceParams = origTraceParams.createScaledCopy(m_angularScale, m_potScale);

	if (!(r = OpenCLSinglePlaneDeflectionInstance::initInstance((uint64_t)this, 
																	m_thetas, posUncertainties, m_intParams,
																	m_floatParams, m_changeableParamIdx,
																	m_kernelCode, m_kernelName,
																	m_extraClPriorCode,
																	m_devIdx,
																	posUncertSeed,
																	originParameterMapping,
																	m_numOriginParams,
																	recalcThetaInfo,
																	*retraceParams
																	)))
		return "Couldn't init OpenCLSinglePlaneDeflectionInstance: " + r.getErrorString();

	if (anyRetrace)
	{
		cerr << "INFO: retrace info: " << origTraceParams.getRetraceDescription() << endl;
		cerr << "INFO: scaled retrace info: " << retraceParams->getRetraceDescription() << endl;
	}

	list<ImagesDataExtended *> empty;

	if (!m_fitObj->init(params.getZd(), imagesPtrsList, empty, &params.getFitnessObjectParameters()))
		return "Unable to initialize fitness object: " + m_fitObj->getErrorString();

	vector<bool> needDefl, needDeriv, needPot, needSecondDeriv;
	m_fitObj->getTotalCalcFlags(needDefl, needDeriv, needPot, needSecondDeriv);

	auto any = [](auto &v) {
		for (auto x : v)
			if (x)
				return true;
		return false;
	};

	m_needDerivs = any(needDeriv);
	m_needPotentials = any(needPot);
	if (m_needPotentials)
	{
		// The projected images interface expects the potential scale to be the
		// angular scale squared, so we'll need a conversion factor here
		// TODO: use a different m_distFrac version that has this scale incorporated?
		m_potScaleConversion = (float)(m_potScale/(m_angularScale*m_angularScale));
	}

	if (any(needSecondDeriv))
		return "Fitness object needs second derivatives of deflection; not suppored";

	m_numObjectives = m_fitObj->getNumberOfFitnessComponents();

	m_initMin = params.getInitMin();
	m_initMax = params.getInitMax();
	m_hardMin = params.getHardMin();
	m_hardMax = params.getHardMax();

	if (m_initMin.size() != m_initMax.size() || m_initMin.size() != m_hardMin.size() ||
	    m_hardMin.size() != m_hardMax.size())
		return "Parameter ranges should all have the same size";

	{
		for (size_t i = 0 ; i < m_initMin.size() ; i++)
		{
			if (params.allowEqualValuesInInitialRange())
			{
				if (!(m_hardMin[i] <= m_initMin[i] && m_initMin[i] <= m_initMax[i] && m_initMax[i] <= m_hardMax[i]))
					return "Error in boundary check for parameter " + std::to_string(i) + ": hardMin = " + std::to_string(m_hardMin[i])
						+ " initMin = " + std::to_string(m_initMin[i]) + " initMax = " + std::to_string(m_initMax[i])
						+ " hardMin = " + std::to_string(m_hardMax[i]);
			}
			else
			{
				if (!(m_hardMin[i] <= m_initMin[i] && m_initMin[i] < m_initMax[i] && m_initMax[i] <= m_hardMax[i]))
					return "Error in boundary check for parameter " + std::to_string(i) + ": hardMin = " + std::to_string(m_hardMin[i])
						+ " initMin = " + std::to_string(m_initMin[i]) + " initMax = " + std::to_string(m_initMax[i])
						+ " hardMin = " + std::to_string(m_hardMax[i]);
			}
		}
	}

	if (m_numOriginParams == 0)
	{
		if (m_initMin.size() != m_changeableParamIdx.size())
			return "Incompatible sizes for bounds and number of changeable parameters";
	}
	else
	{
		if (m_initMin.size() != m_numOriginParams)
			return "Incompatible sizes for bounds and number of origin parameters";
	}

	// Prior stuff

	m_fitnessToAddPriorTo = -1;
	if (m_extraClPriorCode.length() > 0)
	{
		for (size_t i = 0 ; i < m_fitObj->getNumberOfFitnessComponents() ; i++)
		{
			if (m_fitObj->isNegativeLogProb_Overall(i))
			{
				if (m_fitnessToAddPriorTo >= 0)
					return "Found more than one fitness component to add prior info to: had " + to_string(m_fitnessToAddPriorTo) + ", but now also " + to_string(i);
				m_fitnessToAddPriorTo = i;
			}
		}

		if (m_fitnessToAddPriorTo >= 0)
			cerr << "INFO: adding prior information to fitness component " << m_fitnessToAddPriorTo << endl;
		else
			cerr << "INFO: couldn't find any fitness component to add prior information to" << endl;
	}
	else
		cerr << "INFO: no prior information detected" << endl;

	if (m_extraClPriorCode.length() > 0 && m_fitnessToAddPriorTo < 0)
	{
		if (!params.shouldAllowUnusedPriors())
			return "Detected prior information on parameters, but there's no fitness component to add this to (use 'allowUnusedPriors' flag to continue anyway)";
			
		cerr << "INFO: Detected prior information on parameters, but there's no fitness component to add this to - continuing because 'allowUnusedPriors' is set" << endl;
	}
 
	m_genomesToCalculateFitnessFor = params.getGenomesToCalculateFitness();
	for (auto &v : m_genomesToCalculateFitnessFor)
	{
		if (v.size() != m_initMin.size())
			return "Expecting " + to_string(m_initMin.size()) + " parameters in genomes to calculate, but found " + to_string(v.size());
	}

	m_infOnBoundsViolation = params.infinityOnBoundsViolation();
	m_init = true;

	std::cerr << "INFO: m_infOnBoundsViolation = " << ((m_infOnBoundsViolation)?1:0) << std::endl;

	return true;
}

std::shared_ptr<eatk::FitnessComparison> LensGAParametricSinglePlaneCalculator::getFitnessComparison() const
{
	return make_shared<eatk::VectorFitnessComparison<float>>();
}

bool_t LensGAParametricSinglePlaneCalculator::getFullFloatParamsFromOptimizableParameters(const std::vector<float> &optParams, std::vector<float> &fullFloatParams) const
{
	// Bounds have same dimension as parameters
	if (optParams.size() != m_initMin.size())
		return "Expecting " + to_string(m_initMin.size()) + " values, but got " + to_string(optParams.size());

	fullFloatParams = m_floatParams;

	auto fillChangedParamsInFull = [this,&fullFloatParams](const vector<float> &values)
	{
		assert(values.size() == m_changeableParamIdx.size());
		for (size_t i = 0 ; i < m_changeableParamIdx.size() ; i++)
		{
			size_t dst = m_changeableParamIdx[i];
			assert(dst < fullFloatParams.size());

			fullFloatParams[dst] = values[i];
		}
	};

	if (m_numOriginParams == 0)
	{
		// Fill in the solution from the changed parameters
		fillChangedParamsInFull(optParams);
	}
	else // transform the genome values to the actual parameters
	{
		vector<float> changeableParams;

		auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
		bool_t r = cl.getChangeableParametersFromOriginParameters(optParams, changeableParams);
		if (!r)
			return "Can't convert origin parameters to actual ones: " + r.getErrorString();

		fillChangedParamsInFull(changeableParams);
	}

	return true;
}

bool_t LensGAParametricSinglePlaneCalculator::createLens(const eatk::Genome &genome0, std::unique_ptr<GravitationalLens> &resultLens) const
{
	if (!m_init)
		return "Not initialized";

	unique_ptr<GravitationalLens> lens;
	bool_t r;

	assert(dynamic_cast<const eatk::FloatVectorGenome *>(&genome0));
	auto genome = static_cast<const eatk::FloatVectorGenome&>(genome0);
	vector<float> &values = genome.getValues();

	//cerr << "Creating lens from best genome" << genome0.toString() << endl;

	vector<float> fullFloatParams;
	if (!(r = getFullFloatParamsFromOptimizableParameters(values, fullFloatParams)))
		return r;


	//cerr << "Creating lens from Float params:" << endl;
	//for (auto f : fullFloatParams)
	//	cerr << "  " << f << endl;

	// And create a new lens from these values
	lens = m_templateLens->createLensFromCLFloatParams(m_angularScale, m_potScale, fullFloatParams.data());
	if (!lens)
		throw runtime_error("Unexpected: can't create lens from CL float params: " + m_templateLens->getErrorString());

	resultLens = move(lens);
	return true;
}

errut::bool_t LensGAParametricSinglePlaneCalculator::onNewCalculationStart(size_t iteration, size_t genomesForThisCalculator, size_t genomesForPopulationCalculator)
{
	if (!m_init)
		return "Not initialized";
	
	auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
	cl.setTotalGenomesToCalculate(iteration, genomesForPopulationCalculator);
	m_firstCalculationForNewGeneration = true;
	
	return true;
}

errut::bool_t LensGAParametricSinglePlaneCalculator::startNewCalculation(const eatk::Genome &genome0)
{
	const eatk::FloatVectorGenome &genome = static_cast<const eatk::FloatVectorGenome &>(genome0);

	auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
	bool_t r;
	if (!(r = cl.scheduleCalculation(genome)))
		return "Can't start new calculation for genome: " + r.getErrorString();

	return true;
}

errut::bool_t LensGAParametricSinglePlaneCalculator::pollCalculate(const eatk::Genome &genome0, eatk::Fitness &fitness0)
{
	const eatk::FloatVectorGenome &genome = static_cast<const eatk::FloatVectorGenome &>(genome0);
	assert(dynamic_cast<eatk::FloatVectorFitness *>(&fitness0));
	eatk::FloatVectorFitness &fitness = static_cast<eatk::FloatVectorFitness &>(fitness0);

	assert(m_numObjectives == fitness.getValues().size());

	auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
	if (m_firstCalculationForNewGeneration) // See if we can set some new randomized positions
	{
		if (!cl.getAdjustedThetas(m_adjustedThetas))
			return true; // No error, but not ready yet

		m_firstCalculationForNewGeneration = false;
		if (m_adjustedThetas.size() > 0) // Randomization was done
		{
			// Need to use point map to really get the full thetas
			const vector<size_t> &pointMap = m_pointMap.getPointMapping();

			m_fullAdjustedThetas.resize(pointMap.size());
			for (size_t i = 0 ; i < m_fullAdjustedThetas.size() ; i++)
			{
				assert(i < pointMap.size());
				size_t idxForPoint = pointMap[i];

				assert(idxForPoint < m_adjustedThetas.size());
				m_fullAdjustedThetas[i] = m_adjustedThetas[idxForPoint];
			}
			m_oclBp->setAdjustedThetas(m_fullAdjustedThetas);
		}
	}
#ifndef NDEBUG
	else
	{
		// Note: this check won't work anymore if getAdjustedThetas does a swap
		//       instead of a copy!
		vector<Vector2Df> check;
		if (cl.getAdjustedThetas(check))
		{
			assert(check == m_adjustedThetas);
		}
	}
#endif // !NDEBUG

	// This routine can reset the flag that indicates the calculation isn't finished, that's why the
	// getAdjustedThetas code needs to be first
	float negLogClPriorProb = 0;
	if (!cl.getResultsForGenome(genome, m_alphas, m_axx, m_ayy, m_axy, m_potential, m_tracedThetas, m_tracedBetaDiffs, negLogClPriorProb))
		return true; // No error, but not ready yet

	// Check if the bounds were violated (if requested), is for mcmc,
	// In GA or JADE the algorithm should protect agains going out of bounds
	if (m_infOnBoundsViolation)
	{
		assert(dynamic_cast<const eatk::FloatVectorGenome *>(&genome0));
		const vector<float> &values = genome.getValues();

		assert(values.size() == m_hardMin.size());
		bool violated = false;
		for (size_t i = 0 ; i < values.size() ; i++)
		{
			float x = values[i];
			if (x < m_hardMin[i] || x > m_hardMax[i])
			{
				violated = true;
				break;
			}
		}

		if (violated) // Set fitness to infinity, and return
		{
			vector<float> &fitnessValues = fitness.getValues();
			assert(fitnessValues.size() == 1);
			fitnessValues[0] = std::numeric_limits<float>::infinity();
			fitness.setCalculated();
			return true;
		}

		// Ok, not violated, continue the calculation
	}

	/*
	auto lens = createLensFromGenome(genome);
	vector<Vector2Df> alphasCPU;
	for (auto t : m_thetas)
	{
		Vector2Dd thetaDouble(t.getX()*m_angularScale, t.getY()*m_angularScale);
		Vector2Dd alphaCPU;
		if (!lens->getAlphaVector(thetaDouble, &alphaCPU))
			throw runtime_error("Can't get alpha on CPU: " + lens->getErrorString());
		
		alphaCPU /= m_angularScale;
		alphasCPU.push_back({(float)alphaCPU.getX(), (float)alphaCPU.getY()});
	}

	// cerr << "Alpha comparison:" << endl;
	assert(alphasCPU.size() == m_alphas.size());

	// for (size_t i = 0, j = 0 ; i < m_alphas.size() ; i++)
	//  	cerr << m_alphas[i].getX() << "," << m_alphas[i].getY() << " vs CPU "
	//  	     << alphasCPU[i].getX() << "," << alphasCPU[i].getY() << endl;
	float maxDiff = 0;
	for (size_t i = 0, j = 0 ; i < m_alphas.size() ; i++)
	{
		float dx = std::abs(m_alphas[i].getX()-alphasCPU[i].getX());
		float dy = std::abs(m_alphas[i].getY()-alphasCPU[i].getY());
		if (dx > maxDiff)
			maxDiff = dx;
		if (dy > maxDiff)
			maxDiff =  dy;
	}
	cerr << "Max diff: " << maxDiff << endl;
	*/

	const vector<size_t> &pointMap = m_pointMap.getPointMapping();
	
	m_betas.resize(pointMap.size()*2);
	m_scaledAlphas.resize(m_betas.size());
	assert(pointMap.size() == m_distFrac.size());
	for (size_t i = 0, j = 0 ; j < m_betas.size() ; i++, j += 2)
	{
		assert(i < pointMap.size());
		size_t idxForPoint = pointMap[i];
		float f = m_distFrac[i];

		assert(idxForPoint < m_thetas.size() && idxForPoint < m_alphas.size());
		Vector2Df scaledAlpha = m_alphas[idxForPoint];
		scaledAlpha *= f;

		m_betas[j] = m_thetas[idxForPoint].getX() - scaledAlpha.getX();
		m_betas[j+1] = m_thetas[idxForPoint].getY() - scaledAlpha.getY();
		m_scaledAlphas[j] = scaledAlpha.getX();
		m_scaledAlphas[j+1] = scaledAlpha.getY();
	}

	m_oclBp->setBetaBuffer(m_betas.data(), m_betas.size());
	m_oclBp->setAlphaBuffer(m_scaledAlphas.data(), m_alphas.size());

	if (m_needDerivs)
	{
		m_scaledAxx.resize(pointMap.size());
		m_scaledAyy.resize(m_scaledAxx.size());
		m_scaledAxy.resize(m_scaledAxx.size());

		for (size_t i = 0 ; i < m_scaledAxx.size() ; i++)
		{
			assert(i < pointMap.size());
			size_t idxForPoint = pointMap[i];
			float f = m_distFrac[i];

			assert(idxForPoint < m_axx.size() && idxForPoint < m_ayy.size() && idxForPoint < m_axy.size());
			m_scaledAxx[i] = f * m_axx[idxForPoint];
			m_scaledAyy[i] = f * m_ayy[idxForPoint];
			m_scaledAxy[i] = f * m_axy[idxForPoint];
		}

		m_oclBp->setDerivBuffers(m_scaledAxx.data(), m_scaledAyy.data(), m_scaledAxy.data(), m_scaledAxx.size());
	}

	if (m_needPotentials)
	{
		m_scaledPotentials.resize(pointMap.size());

		for (size_t i = 0 ; i < m_scaledPotentials.size() ; i++)
		{
			assert(i < pointMap.size());
			size_t idxForPoint = pointMap[i];
			float f = m_distFrac[i] * m_potScaleConversion; // TODO: use a new array in which the rescaling is included?

			assert(idxForPoint < m_potential.size());
			m_scaledPotentials[i] = f * m_potential[idxForPoint];
		}

		m_oclBp->setPotentialBuffers(m_scaledPotentials.data(), m_scaledPotentials.size());
	}

	if (m_tracedThetas.size() > 0)
	{
		// TODO: also use m_tracedBetaDiffs to see if the convergence made sense? Perhaps just for statistics?
		assert(m_tracedThetas.size() == m_bpPointInfo.size());
		assert(m_tracedBetaDiffs.size() == m_bpPointInfo.size());
		for (size_t i = 0 ; i < m_tracedThetas.size() ; i++)
		{
			auto [ srcIdx, thetaIdx ] = m_bpPointInfo[i];
			assert(m_tracedSourcesPoints.get());
			assert(m_tracedSourcesConvergedFlags.get());
			assert(srcIdx < m_tracedSourcesPoints->size());
			assert(srcIdx < m_tracedSourcesConvergedFlags->size());
			assert(thetaIdx < (*m_tracedSourcesPoints)[srcIdx].size());
			assert(thetaIdx < (*m_tracedSourcesConvergedFlags)[srcIdx].size());

			bool hasConverged = (m_tracedBetaDiffs[i] < m_betaThresHold)?true:false;
			(*m_tracedSourcesPoints)[srcIdx][thetaIdx] = m_tracedThetas[i];
			(*m_tracedSourcesConvergedFlags)[srcIdx][thetaIdx] = hasConverged;

			if (!hasConverged)
				m_stats->m_penaltyCount++; // Keep track for the stats
		}

#ifndef NDEBUG
		for (auto &src : (*m_tracedSourcesPoints))
		{
			for (auto &imgPt : src)
			{
				assert(!std::isnan(imgPt.getX()));
				assert(!std::isnan(imgPt.getY()));
			}
		}
#endif // !NDEBUG

		m_oclBp->setRetraceInfo(m_tracedSourcesFlags, m_tracedSourcesPoints, m_tracedSourcesConvergedFlags);

		// TODO: always do this? Or only in debug mode? Or flag during initialization?
		for (auto x : m_tracedBetaDiffs)
			m_stats->count(x);
	}

	vector<float> &fitnessValues = fitness.getValues();
	if (!m_fitObj->calculateOverallFitness(*m_oclBp, fitnessValues.data()))
		return "Can't calculate fitness: " + m_fitObj->getErrorString();

	if (m_fitnessToAddPriorTo >= 0)
	{
		const vector<float> &values = genome.getValues();

		// Add what was calculated in OpenCL
		fitnessValues[m_fitnessToAddPriorTo] += negLogClPriorProb;
		//cerr << "Added CL prior " << negLogClPriorProb << endl;
	}

	fitness.setCalculated();

#if 0
	ofstream f("/tmp/bpinfo.dat");
	if (f.is_open())
	{
		f << std::setprecision(17);
		f << "#scale " << m_oclBp->getAngularScale() << endl;
		f << "#sources " << m_oclBp->getNumberOfSources() << endl;
		for (int s = 0 ; s < m_oclBp->getNumberOfSources() ; s++)
		{
			f << "#images " << m_oclBp->getNumberOfImages(s) << endl;
			f << "#" << ((m_oclBp->hasRetracedThetas(s))?"retraced":"not_retraced") << endl;
			for (int i = 0 ; i < m_oclBp->getNumberOfImages(s) ; i++)
			{
				f << "#points " << m_oclBp->getNumberOfImagePoints(s, i) << endl;
				for (int p = 0 ; p < m_oclBp->getNumberOfImagePoints(s, i) ; p++)
				{
					f << m_oclBp->getBetas(s, i)[p].getX() << " " << m_oclBp->getBetas(s, i)[p].getY();
					f << " " << m_oclBp->getThetas(s, i)[p].getX() << " " << m_oclBp->getThetas(s, i)[p].getY();
					if (m_oclBp->hasRetracedThetas(s))
						f << " " << m_oclBp->getRetracedThetas(s, i)[p].getX() << " " << m_oclBp->getRetracedThetas(s, i)[p].getY();
					f << endl;
				}
			}
		}
	}
#endif
	return true;
}


} // end namespace
