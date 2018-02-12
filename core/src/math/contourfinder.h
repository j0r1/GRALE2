#ifndef GRALE_CONTOURFINDER_H

#define GRALE_CONTOURFINDER_H

#include "graleconfig.h"
#include "vector2d.h"
#include <vector>

namespace grale
{

class GRALE_IMPORTEXPORT ContourFinder
{
public:
	ContourFinder(const std::vector<double> &values, Vector2Dd bottomLeft, Vector2Dd topRight,
                  int numX, int numY);
	virtual ~ContourFinder();
	
	std::vector<std::vector<Vector2Dd>> findContour(double level);

	Vector2Dd getBottomLeft() const										{ return m_bottomLeft; }
	Vector2Dd getTopRight() const										{ return m_topRight; }
	int getNumX() const													{ return m_numX; }
	int getNumY() const													{ return m_numY; }
	double getXStep() const												{ return m_xstep; }
	double getYStep() const												{ return m_ystep; }
private:
	const std::vector<double> m_values;
	const Vector2Dd m_bottomLeft, m_topRight;
	const int m_numX, m_numY;
	double m_xstep, m_ystep;
};

} // end namespace

#endif // GRALE_CONTOURFINDER_H
