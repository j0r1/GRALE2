#include "gravitationallens.h"
#include "randomnumbergenerator.h"
#include "mathfunctions.h"
#include "lensgastopcriterion.h"
#include "imagesdataextended.h"
#include "constants.h"
#include "openclsingleplanedeflection.h"
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
	MyFitnessCalc() { }

	~MyFitnessCalc()
	{
		OpenCLSinglePlaneDeflectionInstance::releaseInstance((uint64_t)this);
	}

	bool_t init(double angularScale, double potentialScale, const vector<ImagesDataExtended> &images,
	            GravitationalLens &templateLens, const vector<size_t> &changeableParamIdx,
				int devIdx)
	{
		if (m_init)
			return "Already initialized";

		if (changeableParamIdx.size() == 0)
			return "No changeable parameters";

		m_thetas.clear();
		// TODO: we'll need Dds/Ds as well later
		for (auto &img : images)
		{
			int numImages = img.getNumberOfImages();
			for (int i = 0 ; i < numImages ; i++)
			{
				int numPoints = img.getNumberOfImagePoints(i);
				for (int p = 0 ; p < numPoints ; p++)
				{
					Vector2Dd pos = img.getImagePointPosition(i, p);
					pos /= angularScale;
					m_thetas.push_back({(float)pos.getX(), (float)pos.getY()});
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

		bool_t r;
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

		// TODO: make these member variable, better for memory!
		vector<Vector2Df> alphas;
		vector<float> axx, ayy, axy, potential;

		auto &cl = OpenCLSinglePlaneDeflectionInstance::instance();
		if (!cl.getResultsForGenome(genome, alphas, axx, ayy, axy, potential)) // Not ready yet
			return true; // No error, but not ready yet

		//cerr << "Got results for genome " << (void*)&genome << endl;
		
		this_thread::sleep_for(10ms);
		// TODO: set actual value!
		fitness.setCalculated();

		return true;
	}
};

int main0(void)
{
	InversionParameters params = readInversionParameters("templatelens.lensdata", "inversionparams.txt", "images.txt");

	std::shared_ptr<grale::RandomNumberGenerator> rng = std::make_shared<grale::RandomNumberGenerator>();

	MyGA ea;

	eatk::VectorDifferentialEvolutionIndividualCreation<float,float> creation(params.m_initMin, params.m_initMax, rng);
	
	vector<shared_ptr<eatk::GenomeFitnessCalculation>> calculators;
	size_t numCalculators = 4;
	bool_t r;

	for (size_t i = 0 ; i < numCalculators ; i++)
	{
		auto fitCalc = make_shared<MyFitnessCalc>();
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
	stop.initialize(1, convParams); // TODO: single objective for now

	auto mut = make_shared<eatk::VectorDifferentialEvolutionMutation<float>>();
	auto cross = make_shared<eatk::VectorDifferentialEvolutionCrossover<float>>(rng, params.m_hardMin, params.m_hardMax);
	auto fitComp = make_shared<eatk::ValueFitnessComparison<float>>(); // TODO: replace with multi-objective?
	eatk::JADEEvolver evolver(rng, mut, cross, fitComp); // TODO: single objective for now

	r = ea.run(creation, evolver, *popCalc, stop, popSize, popSize, popSize*2);
	if (!r)
		throw runtime_error("Can't run EA: " + r.getErrorString());
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