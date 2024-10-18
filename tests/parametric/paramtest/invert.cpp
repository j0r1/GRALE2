#include "gravitationallens.h"
#include "randomnumbergenerator.h"
#include "mathfunctions.h"
#include "lensgastopcriterion.h"
#include "imagesdataextended.h"
#include "constants.h"
#include <eatk/evolutionaryalgorithm.h>
#include <eatk/vectordifferentialevolution.h>
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
		}

		int idx = pCurImg->addImage();
		pCurImg->addPoint(idx, Vector2Dd(x*ANGLE_ARCSEC, y*ANGLE_ARCSEC));
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

class ParametricFitnessCalculation : public eatk::PopulationFitnessCalculation
{
	// TODO:
};

int main(void)
{
	InversionParameters params = readInversionParameters("templatelens.lensdata", "inversionparams.txt", "images.txt");
	
	std::shared_ptr<grale::RandomNumberGenerator> rng = std::make_shared<grale::RandomNumberGenerator>();

	MyGA ea;

	eatk::VectorDifferentialEvolutionIndividualCreation<float,float> creation(params.m_initMin, params.m_initMax, rng);
	ParametricFitnessCalculation fitFalc;
	Stop stop;
	size_t popSize = 512;

	auto mut = make_shared<eatk::VectorDifferentialEvolutionMutation<float>>();
	auto cross = make_shared<eatk::VectorDifferentialEvolutionCrossover<float>>(rng, params.m_hardMin, params.m_hardMax);
	auto fitComp = make_shared<eatk::VectorFitnessComparison<float>>(); // TODO: replace with multi-objective?
	eatk::JADEEvolver evolver(rng, mut, cross, fitComp); // TODO: single objective for now

	bool_t r = ea.run(creation, evolver, fitFalc, stop, popSize, popSize, popSize*2);
	if (!r)
		throw runtime_error("Can't run EA: " + r.getErrorString());
	return 0;
}
