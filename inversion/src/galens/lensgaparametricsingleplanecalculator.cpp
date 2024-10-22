#include "lensgaparametricsingleplanecalculator.h"
#include "openclsingleplanedeflection.h"
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
	// TODO: we'll need Dds/Ds as well later
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

				if (m_pointMap.addPoint(floatPos)) // returns true for a new point, false if already present
					m_thetas.push_back(floatPos);
	
				m_distFrac.push_back(frac); // for now, we'll use one fraction per point
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

	bool uploadFullParameters = params.alwaysUploadFullParameters();
	if (!(r = OpenCLSinglePlaneDeflectionInstance::initInstance((uint64_t)this, 
																	m_thetas, m_intParams,
																	m_floatParams, m_changeableParamIdx,
																	m_kernelCode, m_kernelName,
																	true, m_devIdx)))
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

	if (any(needPot))
		return "Fitness object needs potential calculations; currently not supported";
	if (any(needSecondDeriv))
		return "Fitness object needs second derivatives of deflection; not suppored";

	m_numObjectives = m_fitObj->getNumberOfFitnessComponents();

	m_initMin = params.getInitMin();
	m_initMax = params.getInitMax();
	m_hardMin = params.getHardMin();
	m_hardMax = params.getHardMax();

	m_init = true;

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

errut::bool_t LensGAParametricSinglePlaneCalculator::onNewCalculationStart(size_t genomesForThisCalculator, size_t genomesForPopulationCalculator)
{
	if (!m_init)
		return "Not initialized";
	
	auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
	cl.setTotalGenomesToCalculate(genomesForPopulationCalculator);
	
	return true;
}

errut::bool_t LensGAParametricSinglePlaneCalculator::startNewCalculation(const eatk::Genome &genome0)
{
	const eatk::FloatVectorGenome &genome = static_cast<const eatk::FloatVectorGenome &>(genome0);

	auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
	cl.scheduleCalculation(genome);

	return true;
}

// TODO: use FloatVectorFitness, for multi-objective option!

errut::bool_t LensGAParametricSinglePlaneCalculator::pollCalculate(const eatk::Genome &genome0, eatk::Fitness &fitness0)
{
	const eatk::FloatVectorGenome &genome = static_cast<const eatk::FloatVectorGenome &>(genome0);
	assert(dynamic_cast<eatk::FloatVectorFitness *>(&fitness0));
	eatk::FloatVectorFitness &fitness = static_cast<eatk::FloatVectorFitness &>(fitness0);

	assert(m_numObjectives == fitness.getValues().size());

	auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
	if (!cl.getResultsForGenome(genome, m_alphas, m_axx, m_ayy, m_axy, m_potential))
		return true; // No error, but not ready yet

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

			if (idxForPoint >= m_axx.size() || idxForPoint >= m_ayy.size() || idxForPoint >= m_axy.size())
				cerr << "idx = " << idxForPoint << " " << m_axx.size() << " " << m_ayy.size() << " " << m_axy.size() << endl;
			assert(idxForPoint < m_thetas.size() && idxForPoint < m_axx.size() && idxForPoint < m_ayy.size() && idxForPoint < m_axy.size());
			m_scaledAxx[i] = f * m_axx[idxForPoint];
			m_scaledAyy[i] = f * m_ayy[idxForPoint];
			m_scaledAxy[i] = f * m_axy[idxForPoint];
		}

		m_oclBp->setDerivBuffers(m_scaledAxx.data(), m_scaledAyy.data(), m_scaledAxy.data(), m_scaledAxx.size());
	}

	vector<float> &fitnessValues = fitness.getValues();
	if (!m_fitObj->calculateOverallFitness(*m_oclBp, fitnessValues.data()))
		return "Can't calculate fitness: " + m_fitObj->getErrorString();
	
	fitness.setCalculated();
	
	return true;
}


} // end namespace
