#include "fitnessutil.h"
#include "fitnesscomponent.h"
#include "projectedimagesinterface.h"
#include "polygon2d.h"
#include "triangle2d.h"
#include <assert.h>
#include <limits>

using namespace std;

namespace grale
{

float calculateCausticPenaltyFitness(const ProjectedImagesInterface &iface,
		const vector<int> &sourceIndices,
		const vector<int> &gridIndices,
		const vector<vector<vector<pair<int,int>>>> &lineSegments,
		const vector<vector<vector<TriangleIndices>>> &critTriangles,
		vector<vector<vector<bool>>> &lineSegmentFlags,
		vector<vector<vector<Vector2D<float>>>> &lineSegmentIntersections,
		FitnessComponentCache *pCache	
		)
{
	assert(sourceIndices.size() == gridIndices.size());
	assert(lineSegments.size() == sourceIndices.size());
	assert(pCache);

	float causticFitness = 0;

	for (size_t sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int s = sourceIndices[sIdx];
		const int g = gridIndices[sIdx];

		for (int i = 0 ; i < lineSegments[sIdx].size() ; i++)
		{
			assert(lineSegments[sIdx].size() == iface.getNumberOfImages(g));

			// consider each 'image' separately

			// first, we calculate the intersections of each line segment with the inversemagnification == 0 plane

			for (int j = 0 ; j < lineSegments[sIdx][i].size() ; j++)
			{
				int idx1 = lineSegments[sIdx][i][j].first;
				int idx2 = lineSegments[sIdx][i][j].second;

				assert(idx1 >= 0 && idx1 < iface.getNumberOfImagePoints(g, i));
				assert(idx2 >= 0 && idx2 < iface.getNumberOfImagePoints(g, i));
				assert(idx1 != idx2);

				float mag1 = iface.getInverseMagnifications(g, i)[idx1];
				float mag2 = iface.getInverseMagnifications(g, i)[idx2];

				if (mag1*mag2 <= 0) // crosses a critical line
				{
					lineSegmentFlags[sIdx][i][j] = true;

					Vector2D<float> p1 = iface.getBetas(g, i)[idx1]; // work directly (approximately) in the source plane
					Vector2D<float> p2 = iface.getBetas(g, i)[idx2];
					float lambda = mag1/(mag1-mag2);
					
					float x = (1.0f-lambda)*p1.getX()+lambda*p2.getX();
					float y = (1.0f-lambda)*p1.getY()+lambda*p2.getY();

					lineSegmentIntersections[sIdx][i][j] = Vector2D<float>(x,y);

					// TODO: for debugging
		//			std::cerr << x << " " << y << std::endl;

				}
				else
					lineSegmentFlags[sIdx][i][j] = false;
			}

			// then, we check the triangles and see which have line segments with a critical line crossing
			vector<pair<Vector2D<float>, Vector2D<float>>> causticParts;

			for (int i = 0 ; i < critTriangles[sIdx].size() ; i++)
			{
				for (int j = 0 ; j < critTriangles[sIdx][i].size() ; j++)
				{
					// here, the triangle indices refer to line segments!

					TriangleIndices t = critTriangles[sIdx][i][j];
					Vector2D<float> intersections[3];
					int count = 0;
					
					for (int k = 0 ; k < 3 ; k++)
					{
						if (lineSegmentFlags[sIdx][i][t.getIndex(k)])
						{
							intersections[count] = lineSegmentIntersections[sIdx][i][t.getIndex(k)];
							count++;
						}
					}

					if (count == 2) // TODO: don't know how to handle other cases
					{
						causticParts.push_back(pair<Vector2D<float>, Vector2D<float>>(intersections[0], intersections[1]));
						// TODO: for debugging
	//							std::cerr << intersections[0].getX() << " " << intersections[0].getY() << std::endl;
	//							std::cerr << intersections[1].getX() << " " << intersections[1].getY() << std::endl;
	//							std::cerr << std::endl;
					}
				}
			}

			// now we have (approximate) line segments for the caustics, check if they intersect
			// the estimated source shape

			Polygon2D<float> bigpoly; 
			float minx = numeric_limits<float>::max();
			float maxx = -numeric_limits<float>::max();
			float miny = numeric_limits<float>::max();
			float maxy = -numeric_limits<float>::max();
			float scale = 0;

			getEstimatedSourceShape(s, iface, pCache, bigpoly, minx, maxx, miny, maxy);
			bool foundScale = pCache->getEstimatedSourceScale(s, scale);
			assert(foundScale);

			for (auto &part : causticParts)
			{
				Vector2D<float> p1 = part.first;
				Vector2D<float> p2 = part.second;

				if (p1.getX() < minx && p2.getX() < minx)
					continue;

				if (p1.getX() > maxx && p2.getX() > maxx)
					continue;

				if (p1.getY() < miny && p2.getY() < miny)
					continue;

				if (p1.getY() > maxy && p2.getY() > maxy)
					continue;

				if (bigpoly.isInside(p1))
				{
					if (bigpoly.isInside(p2)) // both are inside, count total length
					{
						Vector2D<float> diff = (p1-p2)/scale;
						causticFitness += diff.getLength();

						// TODO: for debugging
					//	DrawLine(p1,p2);
					}
					else // p2 is outside
					{
						// calculate intersection, there definitely is one

						const Vector2D<float> *pPolyPoints = bigpoly.getPoints();
						int numPoints = bigpoly.getNumberOfPoints();
						bool found = false;
						Line2D<float> caustLine(p1, p2-p1);

						for (int i = 0 ; i < numPoints && !found ; i++)
						{
							Vector2D<float> direction = pPolyPoints[(i+1)%numPoints] - pPolyPoints[i];
							Line2D<float> polyLine(pPolyPoints[i], direction);
							
							float factor = polyLine.getIntersectionFactor(caustLine);

							if (factor >= 0.0f && factor <= 1.0f)
							{
								Vector2D<float> intersectionPoint = caustLine.getIntersection(polyLine, factor);

								if (factor >= 0.0f && factor <= 1.0f)
								{
									Vector2D<float> diff = (intersectionPoint - p1)/scale;
									causticFitness += diff.getLength();
									found = true;

									// TODO: for debugging
					//				DrawLine(p1, intersectionPoint);
								}
							}
						}
					}
				}
				else // point1 is outside
				{
					if (bigpoly.isInside(p2))
					{
						// calculate intersection, there definitely is one

						const Vector2D<float> *pPolyPoints = bigpoly.getPoints();
						int numPoints = bigpoly.getNumberOfPoints();
						bool found = false;
						Line2D<float> caustLine(p1, p2-p1);

						for (int i = 0 ; i < numPoints && !found ; i++)
						{
							Vector2D<float> direction = pPolyPoints[(i+1)%numPoints] - pPolyPoints[i];
							Line2D<float> polyLine(pPolyPoints[i], direction);
							
							float factor = polyLine.getIntersectionFactor(caustLine);

							if (factor >= 0.0f && factor <= 1.0f)
							{
								Vector2D<float> intersectionPoint = caustLine.getIntersection(polyLine, factor);

								if (factor >= 0.0f && factor <= 1.0f)
								{
									Vector2D<float> diff = (intersectionPoint - p2)/scale;
									causticFitness += diff.getLength();
									found = true;

									// TODO: for debugging
						//			DrawLine(p2, intersectionPoint);
								}
							}
						}
					}
					else // both points outside
					{
						// calculate possible intersections
						const Vector2D<float> *pPolyPoints = bigpoly.getPoints();
						int numPoints = bigpoly.getNumberOfPoints();
						Line2D<float> caustLine(p1, p2-p1);
						Vector2D<float> intersections[2];
						int count = 0;

						for (int i = 0 ; i < numPoints && count < 2 ; i++)
						{
							Vector2D<float> direction = pPolyPoints[(i+1)%numPoints] - pPolyPoints[i];
							Line2D<float> polyLine(pPolyPoints[i], direction);
							
							float factor = polyLine.getIntersectionFactor(caustLine);

							if (factor >= 0.0f && factor <= 1.0f)
							{
								Vector2D<float> intersectionPoint = caustLine.getIntersection(polyLine, factor);

								if (factor >= 0.0f && factor <= 1.0f)
								{
									intersections[count] = intersectionPoint;
									count++;
								}
							}
						}

						if (count == 2)
						{
							Vector2D<float> diff = (intersections[0] - intersections[1])/scale;

							causticFitness += diff.getLength();

							// TODO: for debugging
						//	DrawLine(intersections[0], intersections[1]);
						}
					}
				}
			}
			// TODO: for debugging
			//std::cerr << std::endl;
		}
	}

	assert(!isnan(causticFitness));
	return causticFitness;
}

} // end namespace
