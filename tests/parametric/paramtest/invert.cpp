#include "gravitationallens.h"
#include "randomnumbergenerator.h"
#include "mathfunctions.h"
#include "lensgastopcriterion.h"
#include "imagesdataextended.h"
#include "constants.h"
#include "openclsingleplanedeflection.h"
#include "oclcalculatedbackprojector.h"
#include "lensfitnessgeneral.h"
#include <eatk/evolutionaryalgorithm.h>
#include <eatk/vectordifferentialevolution.h>
#include <eatk/singlethreadedpopulationfitnesscalculation.h>
#include <eatk/multithreadedpopulationfitnesscalculation.h>
#include <eatk/jadeevolver.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <chrono>

using namespace grale;
using namespace std;
using namespace errut;

class Timer {
public:
    Timer() {
        start();
    }

    void start() {
        beg = std::chrono::steady_clock::now();
    }

    void stop() {
        end = std::chrono::steady_clock::now();
    }

    double duration() {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(end - beg).count();
    }
private:
    std::chrono::time_point<std::chrono::steady_clock> beg;
    std::chrono::time_point<std::chrono::steady_clock> end;
};

class MyGA : public eatk::EvolutionaryAlgorithm
{
public:
	MyGA() { }
	~MyGA()
	{
		double dtSum = 0;
		double dt2Sum = 0;
		for (auto dt : m_intervals)
		{
			dtSum += dt;
			dt2Sum += dt*dt;
		}

		double avg = dtSum/m_intervals.size();
		double stddev = SQRT(dt2Sum/m_intervals.size() - avg*avg);
		std::cerr << "Avg = " << avg/1e6 << " ms, stddev = " << stddev/1e6 << " ms" << std::endl;
	}

	Timer m_timer;
	std::vector <double> m_intervals;

	errut::bool_t onBeforeFitnessCalculation(size_t generation, const std::shared_ptr<eatk::Population> &population) override
	{
		m_timer.start();
	// 	cout << "# Generation " << generation << ", before calculation: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
		return true;
	}

    errut::bool_t onFitnessCalculated(size_t generation, const std::shared_ptr<eatk::Population> &population) override
	{
		m_timer.stop();
		m_intervals.push_back(m_timer.duration());
	// 	cout << "# Generation " << generation << ", calculated: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
		return true;
	}

	errut::bool_t onBeforeFitnessCalculation(size_t generation, const std::vector<std::shared_ptr<eatk::Population>> &populations) override
	{
		m_timer.start();
		return true;
	}

    errut::bool_t onFitnessCalculated(size_t generation, const std::vector<std::shared_ptr<eatk::Population>> &populations) override
	{
		m_timer.stop();
		m_intervals.push_back(m_timer.duration());
		return true;
	}
};




struct InversionParameters
{
	double m_zd = 0;
	unique_ptr<GravitationalLens> m_lens;
	double m_angularScale = 0, m_potentialScale = 0;
	size_t m_numVarParams = 0;
	vector<size_t> m_offsets;
	vector<float> m_initMin, m_initMax;
	vector<float> m_hardMin, m_hardMax;
	vector<ImagesDataExtended> m_images;
};

InversionParameters readInversionParameters(const string &templateLens, const string &paramsFile,
                                            const string &imgFile)
{
	unique_ptr<GravitationalLens> lens;
	string errStr;
	if (!GravitationalLens::load(templateLens, lens, errStr))
		throw runtime_error(errStr);

	cerr << "Lens loaded" << endl;

	ifstream config(paramsFile);
	if (!config.is_open())
		throw runtime_error("Can't read config file");

	double angularScale, potentialScale;
	config >> angularScale;
	config >> potentialScale;

	cerr << "angularScale = " << angularScale << endl;
	cerr << "potentialScale = " << potentialScale << endl;

	size_t numVariableParams;
	config >> numVariableParams;
	cerr << "numVariableParams = " << numVariableParams << endl;

	vector<size_t> offsets;
	vector<float> initMin, initMax, hardMin, hardMax;
	auto readVector = [&config, numVariableParams](auto &v)
	{
		using ElementType = typename std::decay_t<decltype(v)>::value_type;
		for (size_t i = 0 ; i < numVariableParams ; i++)
		{
			ElementType x;
			string s;
			config >> s;

			if (s == "inf")
				x = numeric_limits<ElementType>::infinity();
			else if (s == "-inf")
				x = -numeric_limits<ElementType>::infinity();
			else
			{
				stringstream ss(s);
				ss >> x;
			}
			v.push_back(x);
		}
	};

	readVector(offsets);
	readVector(initMin);
	readVector(initMax);
	readVector(hardMin);
	readVector(hardMax);

	auto printVector = [](const string &name, const auto &v)
	{
		cerr << name << " =";
		for (auto x : v)
			cerr << " " << x;
		cerr << endl;
	};

	printVector("initMin", initMin);
	printVector("initMax", initMax);
	printVector("hardMin", hardMin);
	printVector("hardMax", hardMax);

	ifstream imageFile(imgFile);
	if (!imageFile.is_open())
		throw runtime_error("Can't open images file " + imgFile);
	
	string line;
	vector<ImagesDataExtended> images;
	ImagesDataExtended *pCurImg = nullptr;

	getline(imageFile, line); // first line is redshift
	double zd = stod(line);
	cerr << "zd = " << zd << endl;

	while (getline(imageFile, line))
	{
		if (line.length() == 0)
		{
			pCurImg = nullptr;
			continue;
		}

		stringstream ss(line);
		double x, y, frac;
		ss >> x;
		ss >> y;
		ss >> frac;

		if (!pCurImg)
		{
			images.push_back(ImagesDataExtended(1.0, frac));
			pCurImg = &images.back();
			pCurImg->create(0, {});
			pCurImg->setExtraParameter("type", string("pointimages"));
		}

		int idx = pCurImg->addImage();
		if (idx < 0)
			throw runtime_error(pCurImg->getErrorString());
		int pt = pCurImg->addPoint(idx, Vector2Dd(x*ANGLE_ARCSEC, y*ANGLE_ARCSEC));
		if (pt < 0)
			throw runtime_error(pCurImg->getErrorString());
	}

	cerr << "Loaded info about " << images.size() << " sources" << endl;

	InversionParameters p;
	p.m_zd = zd;
	p.m_lens = move(lens);
	p.m_angularScale = angularScale;
	p.m_potentialScale = potentialScale;
	p.m_numVarParams = numVariableParams;
	p.m_offsets = offsets;
	p.m_initMin = initMin;
	p.m_initMax = initMax;
	p.m_hardMin = hardMin;
	p.m_hardMax = hardMax;
	std::swap(images, p.m_images);

	return p;
}


class Stop : public grale::LensGAStopCriterion
{
public:
	Stop()
		: grale::LensGAStopCriterion(0) { }
protected:
	void onReport(const std::string &s)	const override
	{
		cerr << s << endl;
	}
};

class MyFitnessCalc : public eatk::GenomeFitnessCalculation
{
public:
	MyFitnessCalc() = delete;
	MyFitnessCalc(unique_ptr<LensFitnessObject> initializedFitObj)
	{
		m_fitObj = move(initializedFitObj);
	}

	~MyFitnessCalc()
	{
		OpenCLSinglePlaneDeflectionInstance::releaseInstance((uint64_t)this);
	}

	bool_t init(double angularScale, double potentialScale, vector<ImagesDataExtended> &images,
	            GravitationalLens &templateLens, const vector<size_t> &changeableParamIdx,
				int devIdx)
	{
		if (m_init)
			return "Already initialized";

		if (changeableParamIdx.size() == 0)
			return "No changeable parameters";

		bool_t r;
		vector<ImagesDataExtended *> imagesPtrs;
		for (auto &img : images)
			imagesPtrs.push_back(&img);
		
		if (!(r = m_oclBp.init(imagesPtrs, angularScale)))
			return "Can't init backprojector: " + r.getErrorString();

		m_thetas.clear();
		m_distFrac.clear();
		// TODO: we'll need Dds/Ds as well later
		for (auto &img : images)
		{
			int numImages = img.getNumberOfImages();
			double frac = img.getDds()/img.getDs();
			assert(!isnan(frac));

			for (int i = 0 ; i < numImages ; i++)
			{
				int numPoints = img.getNumberOfImagePoints(i);
				for (int p = 0 ; p < numPoints ; p++)
				{
					Vector2Dd pos = img.getImagePointPosition(i, p);
					pos /= angularScale;
					m_thetas.push_back({(float)pos.getX(), (float)pos.getY()});
		
					m_distFrac.push_back(frac); // for now, we'll use one fraction per point
				}
			}
		}

		int numIntParams = 0, numFloatParams = 0;
		if (!templateLens.getCLParameterCounts(&numIntParams, &numFloatParams))
			return "Can't get floating point parameters: " + templateLens.getErrorString();
		
		m_intParams.resize(numIntParams);
		m_floatParams.resize(numFloatParams);
		m_intParams.push_back(-12345);
		m_floatParams.push_back(-12345);
		if (!templateLens.getCLParameters(angularScale, potentialScale, m_intParams.data(), m_floatParams.data()))
			return "Can't get float or in parameters: " + templateLens.getErrorString();

		m_kernelCode = templateLens.getCLLensProgram(angularScale, potentialScale, m_kernelName);
		if (m_kernelCode.length() == 0)
			return "Couldn't get OpenCL kernel code: " + templateLens.getErrorString();

		bool uploadFullParameters = false; // TODO: make this configurable
		if (!(r = OpenCLSinglePlaneDeflectionInstance::initInstance((uint64_t)this, 
																	 m_thetas, m_intParams,
																	 m_floatParams, changeableParamIdx,
																	 m_kernelCode, m_kernelName,
																	 true, devIdx)))
			return "Couldn't init OpenCLSinglePlaneDeflectionInstance: " + r.getErrorString();

		m_angularScale = angularScale;
		m_potScale = potentialScale;
		m_changeableParamIdx = changeableParamIdx;
		m_devIdx = devIdx;
		m_templateLens = templateLens.createCopy();
		m_init = true;

		return true;
	}

	unique_ptr<GravitationalLens> createLensFromGenome(const eatk::Genome &genome0)
	{
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

		return lens;
	}

	OclCalculatedBackProjector m_oclBp;
	bool m_init = false;
	std::unique_ptr<OpenCLSinglePlaneDeflection> m_clDef;
	double m_angularScale = 0;
	double m_potScale = 0;
	int m_devIdx = 0;
	vector<Vector2Df> m_thetas;
	string m_kernelCode;
	string m_kernelName;
	size_t m_numGenomesToCalculate = 0;
	vector<size_t> m_changeableParamIdx;
	vector<int> m_intParams;
	vector<float> m_floatParams;
	std::unique_ptr<GravitationalLens> m_templateLens;

	// This will be called for every generation
	// On the first one we can initialize the OpenCLSinglePlaneDeflection instance
	// (since we then know the number of genomes for this population calculator)
	// In subsequent generations we should check that this hasn't changed
	errut::bool_t onNewCalculationStart(size_t genomesForThisCalculator, size_t genomesForPopulationCalculator) override
	{
		//cerr << "onNewCalculationStart: " << genomesForThisCalculator << "," << genomesForPopulationCalculator << endl;
		if (!m_init)
			return "Not initialized";
		
		auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
		cl.setTotalGenomesToCalculate(genomesForPopulationCalculator);
		
		return true;
	}

	errut::bool_t startNewCalculation(const eatk::Genome &genome0) override
	{
		const eatk::FloatVectorGenome &genome = static_cast<const eatk::FloatVectorGenome &>(genome0);
		//cerr << "Need to calculate: " << genome.toString() << endl;

		auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
		cl.scheduleCalculation(genome);

		return true;
	}
	
	errut::bool_t pollCalculate(const eatk::Genome &genome0, eatk::Fitness &fitness0) override
	{
		const eatk::FloatVectorGenome &genome = static_cast<const eatk::FloatVectorGenome &>(genome0);
		assert(dynamic_cast<eatk::ValueFitness<float> *>(&fitness0));
		eatk::ValueFitness<float> &fitness = static_cast<eatk::ValueFitness<float> &>(fitness0);

		auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
		if (!cl.getResultsForGenome(genome, m_alphas, m_axx, m_ayy, m_axy, m_potential)) // Not ready yet
			return true; // No error, but not ready yet

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

		// cerr << "Calculated betas: " << endl;
		m_betas.resize(m_alphas.size()*2);
		for (size_t i = 0, j = 0 ; i < m_alphas.size() ; i++)
		{
			m_betas[j++] = m_distFrac[i]*m_alphas[i].getX();
			m_betas[j++] = m_distFrac[i]*m_alphas[i].getY();
			// m_betas[j++] = m_distFrac[i]*alphasCPU[i].getX();
			// m_betas[j++] = m_distFrac[i]*alphasCPU[i].getY();

			// cerr << m_alphas[i].getX() << "," << m_alphas[i].getY() << " -> "
			//      << m_betas[j-2] << "," << m_betas[j-1] << endl;
		}

		m_oclBp.setBetaBuffer(m_betas.data(), m_betas.size());

		array<float,16> fitnessValues;
		float sentinel = -98765.0f;
		for (auto &x : fitnessValues)
			x = sentinel; // For now, set a sentinel
		if (!m_fitObj->calculateOverallFitness(m_oclBp, fitnessValues.data()))
			return "Can't calculate fitness: " + m_fitObj->getErrorString();
		fitness.setValue(fitnessValues[0]);

		assert(fitnessValues[1] == sentinel);
		
		fitness.setCalculated();
		
		return true;
	}

	vector<Vector2Df> m_alphas;
	vector<float> m_axx, m_ayy, m_axy, m_potential;
	vector<float> m_betas;
	vector<double> m_distFrac;
	unique_ptr<LensFitnessObject> m_fitObj;
};

int main0(void)
{
	InversionParameters params = readInversionParameters("templatelens.lensdata", "inversionparams.txt", "images.txt");

	std::shared_ptr<grale::RandomNumberGenerator> rng = std::make_shared<grale::RandomNumberGenerator>();

	MyGA ea;

	eatk::VectorDifferentialEvolutionIndividualCreation<float,float> creation(params.m_initMin, params.m_initMax, rng);
	
	vector<shared_ptr<eatk::GenomeFitnessCalculation>> calculators;
	size_t numCalculators = 1;
	bool_t r;

	list<ImagesDataExtended*> imagesPtrs, emptyList;
	for (auto &x : params.m_images)
		imagesPtrs.push_back(&x);

	for (size_t i = 0 ; i < numCalculators ; i++)
	{
		unique_ptr<LensFitnessObject> fitObj = make_unique<LensFitnessGeneral>();
		auto fitParams = fitObj->getDefaultParametersInstance();

		fitParams->setParameter("fitness_pointimageoverlap_scaletype", string("MAD"));

		if (!fitObj->init(params.m_zd, imagesPtrs, emptyList, fitParams.get()))
			throw runtime_error("Can't init fitness object: " + fitObj->getErrorString());

		auto fitCalc = make_shared<MyFitnessCalc>(move(fitObj));
		if (!(r = fitCalc->init(params.m_angularScale, params.m_potentialScale,
	    	               params.m_images, *params.m_lens, params.m_offsets, 0)))
			throw runtime_error("Can't init fitness calculator: " + r.getErrorString());

		calculators.push_back(fitCalc);
	}

	shared_ptr<eatk::PopulationFitnessCalculation> popCalc;
	if (numCalculators == 1)
		popCalc = make_shared<eatk::SingleThreadedPopulationFitnessCalculation>(calculators[0]);
	else if (numCalculators > 1)
	{
		auto calc = make_shared<eatk::MultiThreadedPopulationFitnessCalculation>();
		calc->initThreadPool(calculators);
		popCalc = calc;
	}

	Stop stop;
	size_t popSize = 64;

	LensGAConvergenceParameters convParams;
	convParams.setConvergenceFactor(0.01);
	convParams.setHistorySize(500);
	stop.initialize(1, convParams); // TODO: single objective for now

	auto mut = make_shared<eatk::VectorDifferentialEvolutionMutation<float>>();
	auto cross = make_shared<eatk::VectorDifferentialEvolutionCrossover<float>>(rng, params.m_hardMin, params.m_hardMax);
	auto fitComp = make_shared<eatk::ValueFitnessComparison<float>>(); // TODO: replace with multi-objective?
	eatk::JADEEvolver evolver(rng, mut, cross, fitComp); // TODO: single objective for now

	r = ea.run(creation, evolver, *popCalc, stop, popSize, popSize, popSize*2);
	if (!r)
		throw runtime_error("Can't run EA: " + r.getErrorString());

	auto best = evolver.getBestIndividuals();
	auto calc = static_cast<MyFitnessCalc*>(calculators[0].get());
	auto solution = calc->createLensFromGenome(best[0]->genomeRef());

	cerr << "Writing solution" << endl;
	solution->save("solution.lensdata");
	return 0;
}

int main(void)
{
	try
	{
		main0();
	}
	catch(runtime_error &e)
	{
		cerr << "Got error: " << e.what() << endl;
		return -1;
	}
	return 0;
}