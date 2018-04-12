#ifndef MULTICONTOURFINDER_H

#define MULTICONTOURFINDER_H

#include <errut/errorbase.h>
#include <grale/contourfinder.h>
#include <vector>

class MultiContourFinder : public grale::ContourFinder, public errut::ErrorBase
{
public:
	MultiContourFinder(const std::vector<double> &values, grale::Vector2Dd bottomLeft, grale::Vector2Dd topRight,
			           int numX, int numY);
	~MultiContourFinder();

	bool findContours(const std::vector<double> &levels, int numThreads);
	const std::vector<std::vector<std::vector<grale::Vector2Dd>>> &getContours() const	{ return m_contours; }
private:
	static void findContourThread(MultiContourFinder *instance, const std::vector<double> *levels, const std::vector<int> *levelIdx);

	std::vector<std::vector<std::vector<grale::Vector2Dd>>> m_contours;
};

#endif // MULTICONTOURFINDER_H
