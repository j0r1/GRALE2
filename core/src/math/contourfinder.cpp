#include "contourfinder.h"
#include <list>

using namespace std;

namespace grale
{

typedef pair<Vector2Dd,Vector2Dd> VecPair;
typedef pair<int, int> XY;

ContourFinder::ContourFinder(const vector<double> &values, Vector2Dd bottomLeft, Vector2Dd topRight,
                             int numX, int numY)
	: m_values(values),
	  m_bottomLeft(bottomLeft),
	  m_topRight(topRight),
	  m_numX(numX),
	  m_numY(numY)
{
	m_xstep = (m_topRight.getX() - m_bottomLeft.getX())/(m_numX-1);
	m_ystep = (m_topRight.getY() - m_bottomLeft.getY())/(m_numY-1);
}

ContourFinder::~ContourFinder()
{
}

inline Vector2Dd lineIntersection(double v0, double v1, Vector2Dd p0, Vector2Dd p1)
{
	Vector2Dd diff = p1;
	diff -= p0;

	double frac = (0.0-v0)/(v1-v0);
	Vector2Dd result = p0;
	result += frac*diff;
	return result;
}

inline bool hasIntersection(double v0, double v1, double v2, 
									Vector2Dd p0, Vector2Dd p1, Vector2Dd p2,
		                            Vector2Dd &pt1, Vector2Dd &pt2)
{
	// At this point, we're assuming that the length of the sides is 1
	vector<Vector2Dd> points;

	if (v0*v1 < 0)
		points.push_back(lineIntersection(v0, v1, p0, p1));
	if (v0*v2 < 0)
		points.push_back(lineIntersection(v0, v2, p0, p2));
	if (v1*v2 < 0)
		points.push_back(lineIntersection(v1, v2, p1, p2));
		
	if (points.size() >= 2)
	{
		pt1 = points[0];
		pt2 = points[1];
		return true;
	}
	return false;
}

inline bool touches(Vector2Dd p0, Vector2Dd p1)
{
	p0 -= p1;
	if (p0.getLength() < 1e-6)
		return true;
	return false;
}

int growCurrentLine(list<Vector2Dd> &currentLine, list<VecPair> &segments)
{
	int count = 0;

	auto it = segments.begin();
	while (it != segments.end())
	{
		VecPair line = *it;
		bool segmentAccepted = false;

		if (currentLine.empty())
		{
			currentLine.push_back(line.first);
			currentLine.push_back(line.second);
			segmentAccepted = true;
			count++;
		}
		else
		{
			if (touches(currentLine.front(), line.first))
			{
				currentLine.push_front(line.second);
				segmentAccepted = true;
				count++;
			}
			else if (touches(currentLine.front(), line.second))
			{
				currentLine.push_front(line.first);
				segmentAccepted = true;
				count++;
			}
			else if (touches(currentLine.back(), line.first))
			{
				currentLine.push_back(line.second);
				segmentAccepted = true;
				count++;
			}
			else if (touches(currentLine.back(), line.second))
			{
				currentLine.push_back(line.first);
				segmentAccepted = true;
				count++;
			}
		}

		if (!segmentAccepted)
			++it;
		else
			it = segments.erase(it);
	}

	return count;
}

vector<vector<Vector2Dd>> ContourFinder::findContour(double level)
{
	int numTotal = m_numX*m_numY;

	vector<list<VecPair>> segmentsPerPixel(numTotal);

	// Find the intersections for all triangles
	for (int y = 0 ; y < m_numY-1 ; y++)
	{
		for (int x = 0 ; x < m_numX-1 ; x++)
		{
			int idx = x+y*m_numX;
			double v00 = m_values[idx] - level;
			double v01 = m_values[idx+1] - level;
			double v10 = m_values[idx+m_numX] - level;
			double v11 = m_values[idx+m_numX+1] - level;

			// consider the v00, v01, v11 triangle
			{
				Vector2Dd pt1, pt2;
				if (hasIntersection(v00, v11, v01,
							        Vector2Dd(x, y), Vector2Dd(x+1, y+1), Vector2Dd(x+1, y),
								    pt1, pt2))
					segmentsPerPixel[idx].push_back(VecPair(pt1, pt2));
			}

			// consider v00, v11, v10 triangle
			{
				Vector2Dd pt1, pt2;
				if (hasIntersection(v00, v11, v10, 
							        Vector2Dd(x, y), Vector2Dd(x+1, y+1), Vector2Dd(x, y+1),
							        pt1, pt2))
					segmentsPerPixel[idx].push_back(VecPair(pt1, pt2));
			}
		}
	}

	vector<vector<Vector2Dd>> contours;

	// Link the segments together
	bool done = false;
	while (!done)
	{
		int startIdx = -1;
		list<XY> pixelsToConsider;

		// Look for a starting pixel
		for (int y = 0 ; startIdx < 0 && y < m_numY-1 ; y++)
		{
			for (int x = 0 ; startIdx < 0 && x < m_numX-1 ; x++)
			{
				int idx = x+y*m_numX;
				if (segmentsPerPixel[idx].size() > 0)
				{
					startIdx = idx;
					pixelsToConsider.push_back(XY(x,y));
				}
			}
		}

		if (startIdx < 0)
			done = true;
		else
		{
			list<Vector2Dd> currentLine;
			
			while (pixelsToConsider.size() > 0)
			{
				XY xy = pixelsToConsider.front();
				int x = xy.first;
				int y = xy.second;

				if (segmentsPerPixel[x+y*m_numX].size() == 0) // nothing there, just remove this pixel to consider
					pixelsToConsider.pop_front();
				else
				{
					auto &segments = segmentsPerPixel[x+y*m_numX];
					int numGrown = growCurrentLine(currentLine, segments);

					if (numGrown == 0)
						pixelsToConsider.pop_front();
					else
					{
						if (segments.size() == 0) // remove this pixel to consider
							pixelsToConsider.pop_front();

						// Add the neighbours
						if (x < m_numX-1)
							pixelsToConsider.push_back(XY(x+1, y));
						if (x > 0)
							pixelsToConsider.push_back(XY(x-1, y));
						if (y < m_numY-1)
							pixelsToConsider.push_back(XY(x, y+1));
						if (y > 0)
							pixelsToConsider.push_back(XY(x, y-1));
					}
				}
			}

			if (currentLine.size() > 0) // should always be the case
			{
				vector<Vector2Dd> l;
				
				for (auto p : currentLine)
					l.push_back(p);

				contours.push_back(l);
			}
		}
	}

	// Change coordinates based on the rect
	for (auto &contour : contours)
	{
		for (auto &pt : contour)
		{
			double x = pt.getX();
			double y = pt.getY();

			Vector2Dd newCoord(x * m_xstep, y * m_ystep);
			newCoord += m_bottomLeft;

			pt = newCoord;
		}
	}

	return contours;
}

} // end namespace
