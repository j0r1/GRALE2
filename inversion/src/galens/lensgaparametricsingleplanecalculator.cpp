#include "lensgaparametricsingleplanecalculator.h"
#include "openclsingleplanedeflection.h"
#include "utils.h"
#include <eatk/vectorgenomefitness.h>
#include <eatk/valuefitness.h>
#include <iostream>
#include <list>
#include <cassert>

using namespace errut;
using namespace std;

namespace grale
{

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

	if (numChangeAbleParams != params.getInitMin().size() ||
	    numChangeAbleParams != params.getInitMax().size() ||
		numChangeAbleParams != params.getHardMin().size() ||
		numChangeAbleParams != params.getHardMax().size())
		return "Internal error: Incompatibility in changeable parameter sizes";

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
	bool uploadFullParams = params.alwaysUploadFullParameters();

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

	for (auto &img : params.getImages())
	{
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

				int ptIdx = m_pointMap.addPoint(floatPos);
				if (ptIdx < 0) // returns -1 for a new point, point index for an existing one
				{
					m_thetas.push_back(floatPos);
					if (havePosUncerts)
						posUncertainties.push_back(uncert);
				}
				else
				{
					if (havePosUncerts)
					{
						assert(ptIdx < posUncertainties.size());
						if (posUncertainties[ptIdx] != uncert)
							return "Points with same positions have different uncertainties";
					}
				}
	
				m_distFrac.push_back(frac);
			}
		}
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

	bool uploadFullParameters = params.alwaysUploadFullParameters();
	if (!(r = OpenCLSinglePlaneDeflectionInstance::initInstance((uint64_t)this, 
																	m_thetas, posUncertainties, m_intParams,
																	m_floatParams, m_changeableParamIdx,
																	m_kernelCode, m_kernelName,
																	uploadFullParameters, m_devIdx,
																	posUncertSeed
																	)))
		return "Couldn't init OpenCLSinglePlaneDeflectionInstance: " + r.getErrorString();

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

	m_infOnBoundsViolation = params.infinityOnBoundsViolation();
	m_init = true;

	std::cerr << "INFO: m_infOnBoundsViolation = " << ((m_infOnBoundsViolation)?1:0) << std::endl;

	return true;
}

std::shared_ptr<eatk::FitnessComparison> LensGAParametricSinglePlaneCalculator::getFitnessComparison() const
{
	return make_shared<eatk::VectorFitnessComparison<float>>(); // TODO: replace with multi-objective?
}

bool_t LensGAParametricSinglePlaneCalculator::createLens(const eatk::Genome &genome0, std::unique_ptr<GravitationalLens> &resultLens) const
{
	if (!m_init)
		return "Not initialized";

	unique_ptr<GravitationalLens> lens;

	assert(dynamic_cast<const eatk::FloatVectorGenome *>(&genome0));
	auto genome = static_cast<const eatk::FloatVectorGenome&>(genome0);
	vector<float> &values = genome.getValues();

	assert(values.size() == m_changeableParamIdx.size());

	// Fill in the solution
	vector<float> fullFloatParams = m_floatParams;
	for (size_t i = 0 ; i < m_changeableParamIdx.size() ; i++)
	{
		size_t dst = m_changeableParamIdx[i];
		assert(dst < fullFloatParams.size());

		fullFloatParams[dst] = values[i];
	}

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
	
	return true;
}

errut::bool_t LensGAParametricSinglePlaneCalculator::startNewCalculation(const eatk::Genome &genome0)
{
	const eatk::FloatVectorGenome &genome = static_cast<const eatk::FloatVectorGenome &>(genome0);

	auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
	cl.scheduleCalculation(genome);

	return true;
}

errut::bool_t LensGAParametricSinglePlaneCalculator::pollCalculate(const eatk::Genome &genome0, eatk::Fitness &fitness0)
{
	const eatk::FloatVectorGenome &genome = static_cast<const eatk::FloatVectorGenome &>(genome0);
	assert(dynamic_cast<eatk::FloatVectorFitness *>(&fitness0));
	eatk::FloatVectorFitness &fitness = static_cast<eatk::FloatVectorFitness &>(fitness0);

	assert(m_numObjectives == fitness.getValues().size());

	auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
	if (!cl.getResultsForGenome(genome, m_alphas, m_axx, m_ayy, m_axy, m_potential))
		return true; // No error, but not ready yet

	// Check if the bounds were violated (if requested), is for mcmc
	if (m_infOnBoundsViolation)
	{
		assert(dynamic_cast<const eatk::FloatVectorGenome *>(&genome0));
		const eatk::FloatVectorGenome &genome = static_cast<const eatk::FloatVectorGenome&>(genome0);
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

	vector<float> &fitnessValues = fitness.getValues();
	if (!m_fitObj->calculateOverallFitness(*m_oclBp, fitnessValues.data()))
		return "Can't calculate fitness: " + m_fitObj->getErrorString();
	
	fitness.setCalculated();
	
	return true;
}


} // end namespace
