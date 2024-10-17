#include "gravitationallens.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>

using namespace grale;
using namespace std;

struct InversionParameters
{
	unique_ptr<GravitationalLens> m_lens;
	double m_angularScale = 0, m_potentialScale = 0;
	size_t m_numVarParams = 0;
	vector<size_t> m_offsets;
	vector<float> m_initMin, m_initMax;
	vector<float> m_hardMin, m_hardMax;
};

InversionParameters readInversionParameters(const string &templateLens, const string &paramsFile)
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

	return p;
}

int main(void)
{
	InversionParameters params = readInversionParameters("templatelens.lensdata", "inversionparams.txt");

	return 0;
}
