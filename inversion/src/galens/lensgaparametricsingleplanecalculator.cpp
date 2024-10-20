#include "lensgaparametricsingleplanecalculator.h"
#include "openclsingleplanedeflection.h"
#include <eatk/vectorgenomefitness.h>
#include <eatk/valuefitness.h>
#include <iostream>
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
	for (auto &img : params.getImages())
		imagesPtrs.push_back(img.get());
	
	m_oclBp = make_unique<OclCalculatedBackProjector>();
	if (!(r = m_oclBp->init(imagesPtrs, m_angularScale)))
		return "Can't init backprojector: " + r.getErrorString();

	m_thetas.clear();
	m_distFrac.clear();
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
				m_thetas.push_back({(float)pos.getX(), (float)pos.getY()});
	
				m_distFrac.push_back(frac); // for now, we'll use one fraction per point
			}
		}
	}

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

	bool uploadFullParameters = false; // TODO: make this configurable
	if (!(r = OpenCLSinglePlaneDeflectionInstance::initInstance((uint64_t)this, 
																	m_thetas, m_intParams,
																	m_floatParams, m_changeableParamIdx,
																	m_kernelCode, m_kernelName,
																	true, m_devIdx)))
		return "Couldn't init OpenCLSinglePlaneDeflectionInstance: " + r.getErrorString();

	m_init = true;

	return true;
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
	assert(dynamic_cast<eatk::ValueFitness<float> *>(&fitness0));
	eatk::ValueFitness<float> &fitness = static_cast<eatk::ValueFitness<float> &>(fitness0);

	auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
	if (!cl.getResultsForGenome(genome, m_alphas, m_axx, m_ayy, m_axy, m_potential)) // Not ready yet
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

	// cerr << "Calculated betas: " << endl;
	m_betas.resize(m_alphas.size()*2);
	for (size_t i = 0, j = 0 ; i < m_alphas.size() ; i++)
	{
		m_betas[j++] = m_thetas[i].getX()-m_distFrac[i]*m_alphas[i].getX();
		m_betas[j++] = m_thetas[i].getY()-m_distFrac[i]*m_alphas[i].getY();
		// m_betas[j++] = m_distFrac[i]*alphasCPU[i].getX();
		// m_betas[j++] = m_distFrac[i]*alphasCPU[i].getY();

		// cerr << m_alphas[i].getX() << "," << m_alphas[i].getY() << " -> "
		//      << m_betas[j-2] << "," << m_betas[j-1] << endl;
	}

	m_oclBp->setBetaBuffer(m_betas.data(), m_betas.size());

	array<float,16> fitnessValues;
	float sentinel = -98765.0f;
	for (auto &x : fitnessValues)
		x = sentinel; // For now, set a sentinel
	if (!m_fitObj->calculateOverallFitness(*m_oclBp, fitnessValues.data()))
		return "Can't calculate fitness: " + m_fitObj->getErrorString();
	fitness.setValue(fitnessValues[0]);

	assert(fitnessValues[1] == sentinel);
	
	fitness.setCalculated();
	
	return true;
}


} // end namespace