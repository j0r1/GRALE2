#pragma once

#include "timer.h"
#include "mathfunctions.h"
#include <eatk/evolutionaryalgorithm.h>
#include <iostream>

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


