#include "multicontourfinder.h"
#include <thread>

using namespace std;
using namespace grale;

MultiContourFinder::MultiContourFinder(const std::vector<double> &values, Vector2Dd bottomLeft, Vector2Dd topRight,
									   int numX, int numY) : grale::ContourFinder(values, bottomLeft, topRight, numX, numY)
{
}

MultiContourFinder::~MultiContourFinder()
{
}

void MultiContourFinder::findContourThread(MultiContourFinder *instance, const std::vector<double> *levels, const std::vector<int> *levelIdx)
{
	for (int idx : *levelIdx)
		instance->m_contours[idx] = instance->findContour((*levels)[idx]);
}

bool MultiContourFinder::findContours(const std::vector<double> &levels, int numThreads)
{
	if (levels.size() == 0)
	{
		setErrorString("No levels specified");
		return false;
	}

	if (numThreads <= 0)
	{
		setErrorString("The number of threads must be positive");
		return false;
	}
	
	m_contours.resize(levels.size());

	vector<vector<int>> levelIdx(numThreads);
	for (size_t i = 0 ; i < levels.size() ; i++)
		levelIdx[i%numThreads].push_back(i);

	vector<thread*> threads(numThreads);
	for (int i = 0 ; i < numThreads ; i++)
		threads[i] = new thread(findContourThread, this, &levels, &(levelIdx[i]));

	for (int i = 0 ; i < numThreads ; i++)
	{
		threads[i]->join();
		delete threads[i];
	}

	return true;
}
