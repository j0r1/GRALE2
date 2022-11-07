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

// Just take the backprojected image and count the number of backprojected
// triangles that enclose it (as an estimate of the number of images this lens
// model will generate)

int getSingleImageNullSpaceHitCount(Vector2D<float> beta, const vector<TriangleIndices> &triangles,
									const Vector2D<float> *nullBetas)
{
	const float x = beta.getX();
	const float y = beta.getY();
	int hitCount = 0;
	
	for (int idx = 0 ; idx < triangles.size() ; idx++)
	{
		TriangleIndices t = triangles[idx];
		Vector2D<float> p[3];

		p[0] = nullBetas[t.getIndex(0)];
		p[1] = nullBetas[t.getIndex(1)];
		p[2] = nullBetas[t.getIndex(2)];

		if (p[0].getX() < x && p[1].getX() < x && p[2].getX() < x)
			continue;

		if (p[0].getX() > x && p[1].getX() > x && p[2].getX() > x)
			continue;

		if (p[0].getY() < y && p[1].getY() < y && p[2].getY() < y)
			continue;

		if (p[0].getY() > y && p[1].getY() > y && p[2].getY() > y)
			continue;

		Triangle2D<float> bpTriang(p[0],p[1],p[2]);
		if (bpTriang.isInside(beta))
			hitCount++;
	}

	if (hitCount < 1) // in other code this needs to be >= 1
		hitCount = 1;

	return hitCount;
}

float calculateNullFitness_PointImages(const PointGroupStorage &pointGroups,
		                               const ProjectedImagesInterface &interface, 
                                       const vector<int> &sourceIndices,
									   const vector<int> &nullIndices,
									   const vector<vector<TriangleIndices> > &nullTriangles,
									   const vector<float> &nullWeights
									   )
{
	assert(sourceIndices.size() == nullIndices.size());
	assert(nullIndices.size() == nullTriangles.size());
	assert(nullIndices.size() == nullWeights.size());
	assert(pointGroups.getNumberOfSources() == sourceIndices.size());

	float nullFitness = 0;

	for (int idx = 0 ; idx < sourceIndices.size() ; idx++)
	{
		const int s = sourceIndices[idx];
		const int np = nullIndices[idx];
		const int tp = idx;
		const int weigthIdx = idx;

		assert(np != s);
		assert(s >= 0 && s < interface.getNumberOfSources());
		assert(np >= 0 && np < interface.getNumberOfSources());
		assert(tp >= 0 && tp < nullTriangles.size());
		assert(weigthIdx >= 0 && weigthIdx < nullWeights.size());

		const vector<TriangleIndices> &triangles = nullTriangles[tp];
		if (triangles.size() == 0) // If no null space is used for this source, there should be no triangles
			continue;

		const Vector2D<float> *nullBetas = interface.getBetas(np); // these are the null space points
		int numImages = interface.getNumberOfImages(s);
		float nullenergy = 0;

		if (numImages > 1) 
		{
			// In this case, there are multiple images. We're going to take
			// the convex null of the backprojected images. We're using
			// point groups here, so we can also use extended images if
			// they have corresponding points marked

			assert(idx < pointGroups.getNumberOfSources());
			int numGrps = pointGroups.getNumberOfGroups(idx);

			for (int g = 0 ; g < numGrps ; g++)
			{
				auto getGroupPoint = [idx, s, g, &pointGroups, &interface](int ptIdx) -> Vector2Df 
				{
					const int srcIdx = idx;
					assert(srcIdx >= 0 && srcIdx < pointGroups.getNumberOfSources());
					assert(g >= 0 && g < pointGroups.getNumberOfGroups(srcIdx));
					assert(ptIdx >= 0 && ptIdx < pointGroups.getNumberOfGroupPoints(srcIdx, g));
					int imgNum, ptNum;

					pointGroups.getGroupPointIndices(srcIdx, g, ptIdx, &imgNum, &ptNum);
					assert(imgNum >= 0 && imgNum < interface.getNumberOfImages(s));
					assert(ptNum >= 0 && ptNum < interface.getNumberOfImagePoints(s, imgNum));
					return interface.getBetas(s, imgNum)[ptNum];
				};

				vector<Vector2D<float> > allpoints; // these will contain the backprojected image points
				float maxx = -numeric_limits<float>::max();
				float maxy = -numeric_limits<float>::max();
				float minx = numeric_limits<float>::max();
				float miny = numeric_limits<float>::max();

				auto updateMinMaxAllPoints = [&allpoints, &maxx, &maxy, &minx, &miny](Vector2Df v) {
					maxx = MAX(maxx, v.getX());
					maxy = MAX(maxy, v.getY());
					minx = MIN(minx, v.getX());
					miny = MIN(miny, v.getY());
					allpoints.push_back(v);
				};

				int numPts = pointGroups.getNumberOfGroupPoints(idx, g);
				if (numPts > 0)
					updateMinMaxAllPoints(getGroupPoint(0));

				for (int p = 1 ; p < numPts ; p++)
					updateMinMaxAllPoints(getGroupPoint(p));

				// If we only have two points, add a few points so that we can use the rest
				// of the code
				if (numPts == 2)
				{
					float centerX = (maxx+minx)/2.0;
					float centerY = (maxy+miny)/2.0;
					float dx = maxx-minx;
					float dy = maxy-miny;
					float r = SQRT(dx*dx+dy*dy)/50.0;

					float angles[2] = { (float)(45.0*CONST_PI/180.0), (float)(135.0*CONST_PI/180.0) };
					float extraPoints[2][2] = { { centerX+r*COS(angles[0]), centerY+r*SIN(angles[0]) },
												{ centerX+r*COS(angles[1]), centerY+r*SIN(angles[1]) } };

					for (int extra = 0 ; extra < 2 ; extra++)
					{
						float px = extraPoints[extra][0];
						float py = extraPoints[extra][1];
						allpoints.push_back(Vector2D<float>(px, py));
						maxx = MAX(maxx,px);
						maxy = MAX(maxy,py);
						minx = MIN(minx,px);
						miny = MIN(miny,py);
					}
				}

				if (allpoints.size() == 0) // no point, nothing to do
				{
				}
				else if (allpoints.size() == 1) // one point, use it as a constraint
				{
					nullenergy += getSingleImageNullSpaceHitCount(allpoints[0], nullTriangles[tp], nullBetas);
				}
				else // use convex hull and calculate overlap
				{
					Polygon2D<float> bigpoly;
					bigpoly.init(allpoints, true); // calculate convexhull

					// Now search for triangles overlapping the source
					
					for (int idx = 0 ; idx < triangles.size() ; idx++)
					{
						TriangleIndices t = nullTriangles[tp][idx];
						Vector2D<float> p[3];

						p[0] = nullBetas[t.getIndex(0)];
						p[1] = nullBetas[t.getIndex(1)];
						p[2] = nullBetas[t.getIndex(2)];

						if (p[0].getX() < minx && p[1].getX() < minx && p[2].getX() < minx)
							continue;

						if (p[0].getX() > maxx && p[1].getX() > maxx && p[2].getX() > maxx)
							continue;

						if (p[0].getY() < miny && p[1].getY() < miny && p[2].getY() < miny)
							continue;

						if (p[0].getY() > maxy && p[1].getY() > maxy && p[2].getY() > maxy)
							continue;

						// Calculate overlapping area

						Triangle2D<float> nulltriangle(p[0],p[1],p[2]);
						const Vector2D<float> *hullpoints = bigpoly.getPoints();

						int num = bigpoly.getNumberOfPoints() - 2;
						float totalCount = 0;

						for (int i = 0 ; i < num ; i++)
						{
							Triangle2D<float> polypart(hullpoints[0],hullpoints[1+i],hullpoints[2+i]);

							if (polypart.getOverlapArea(nulltriangle) > 0)
							{
								// Just estimate the number of null space triangles that overlap the source area
								// As soon as one of the polygon parts overlaps with the triangle, we should
								// count the triangle
								totalCount += 1.0f;
								break;
							}
						}

						nullenergy += totalCount; 
					}
				}
			}
		}
		else if (numImages == 1) 
		{
			// use every single point in this image as a separate constraint

			const int numPoints = interface.getNumberOfImagePoints(s, 0);
			const Vector2D<float> *pBetas = interface.getBetas(s, 0);
			int hitCount = 0;
			for (int p = 0 ; p < numPoints ; p++)
			{
				auto beta = pBetas[p];

				assert(tp >= 0 && tp < nullTriangles.size());
				hitCount += getSingleImageNullSpaceHitCount(beta, nullTriangles[tp], nullBetas);
			}

			nullenergy = (float)hitCount;
		}

		assert(weigthIdx < nullWeights.size());

		nullFitness += nullenergy * nullWeights[weigthIdx];
	}

	assert(!isnan(nullFitness));
	return nullFitness;
}

void getEstimatedSourceShape(int s, const ProjectedImagesInterface &iface, FitnessComponentCache *pCache,
		Polygon2D<float> &estimatedShape, float &minX, float &maxX, float &minY, float &maxY)
{
	if (pCache != 0)
	{
		if (pCache->getEstimatedSource(s, estimatedShape, minX, maxX, minY, maxY))
			return; // ok, found it in the cache
	}

	// Not found in cache, calculate the convex hull of the back-projected images
	
	int numimages = iface.getNumberOfImages(s); 
	vector<Vector2D<float>> allpoints; 

	float maxx = -numeric_limits<float>::max();
	float maxy = -numeric_limits<float>::max();
	float minx = numeric_limits<float>::max();
	float miny = numeric_limits<float>::max();

	auto updateMinMax = [&maxx, &maxy, &minx, &miny](Vector2Df p)
	{
		maxx = MAX(maxx, p.getX()); 
		maxy = MAX(maxy, p.getY()); 
		minx = MIN(minx, p.getX()); 
		miny = MIN(miny, p.getY()); 
	};

	// first calculate the convex hull for each image
	for (int i = 0 ; i < numimages ; i++) 
	{ 
		int numpoints = iface.getNumberOfImagePoints(s,i); 
		const Vector2D<float> *sourceBetas = iface.getBetas(s,i); 
		
		if (numpoints > 2)
		{
			Polygon2D<float> poly;

			poly.init(sourceBetas, numpoints, true);
			int numPolyPoints = poly.getNumberOfPoints(); 
			const Vector2D<float> *points = poly.getPoints(); 
			
			for (int p = 0 ; p < numPolyPoints ; p++) 
			{
				updateMinMax(points[p]);
				allpoints.push_back(points[p]);
			}
		}
		else
		{
			for (int i = 0 ; i < numpoints ; i++)
			{
				updateMinMax(sourceBetas[i]);
				allpoints.push_back(sourceBetas[i]);
			}
		}
	} 

	// Then calculate the convex hull of these convex hulls as the
	// current estimate of the source	
	assert(allpoints.size() > 1); // Should be more than one image with at least one point

	if (allpoints.size() < 3)
	{
		// add a small circle
		Vector2Df center = allpoints[0];
		float radius = 1e-6; // TODO: this should never actually be used

		if (allpoints.size() == 2)
		{
			center += allpoints[1];
			center *= 0.5;

			Vector2Df diff = allpoints[0];
			diff -= center;

			radius = diff.getLength();
		}

		int num = 8;
		for (int i = 0 ; i < num ; i++)
		{
			float x = center.getX() + radius*COS(2.0f*(float)CONST_PI*(float)i/(float)num);
			float y = center.getY() + radius*SIN(2.0f*(float)CONST_PI*(float)i/(float)num);
			Vector2Df v(x, y);

			allpoints.push_back(v);
			updateMinMax(v);
		}
	}
	
	Polygon2D<float> bigpoly; 
	bigpoly.init(allpoints, true); 

	// Store the final result in the cache (if possible)
	if (pCache)
		pCache->setEstimatedSource(s, bigpoly, minx, maxx, miny, maxy);

	// Make sure the caller gets the answer 
	estimatedShape = bigpoly;
	minX = minx;
	maxX = maxx;
	minY = miny;
	maxY = maxy;

	assert(!isnan(minX));
	assert(!isnan(maxX));
	assert(!isnan(minY));
	assert(!isnan(maxY));
}

float calculateNullFitness_ExtendedImages(const ProjectedImagesInterface &iface, 
                                       const vector<int> &sourceIndices,
									   const vector<int> &nullIndices,
									   const vector<vector<TriangleIndices> > &nullTriangles,
									   const vector<vector<double> > &nullTriangleAreas,
									   const vector<float> &nullWeights,
									   FitnessComponentCache *pCache
									   )
{
	assert(sourceIndices.size() == nullIndices.size());
	assert(nullIndices.size() == nullTriangles.size());
	assert(nullIndices.size() == nullWeights.size());

	float nullfitness = 0;

	for (int idx = 0 ; idx < sourceIndices.size() ; idx++)
	{
		const int s = sourceIndices[idx];
		const int np = nullIndices[idx];
		const int tp = idx;
		const int weigthIdx = idx;
	
		assert(np != s);
		assert(s >= 0 && s < iface.getNumberOfSources());
		assert(np >= 0 && np < iface.getNumberOfSources());
		assert(tp >= 0 && tp < nullTriangles.size() && tp < nullTriangleAreas.size());
		assert(weigthIdx >= 0 && weigthIdx < nullWeights.size());

		const vector<TriangleIndices> &triangles = nullTriangles[tp]; 
		if (triangles.size() == 0) // If no null space is used for this source, there should be no triangles
			continue;

		int numimages = iface.getNumberOfImages(s); 
		const Vector2D<float> *nullBetas = iface.getBetas(np); 
		float nullenergy = 0; 
		
		if (numimages > 1)
		{
			float maxx = -numeric_limits<float>::max();
			float maxy = -numeric_limits<float>::max();
			float minx = numeric_limits<float>::max();
			float miny = numeric_limits<float>::max();
			Polygon2D<float> bigpoly; 

			getEstimatedSourceShape(s, iface, pCache, bigpoly, minx, maxx, miny, maxy);
			
			// Now search for triangles overlapping the source 

			list<double>::const_iterator it2; 
			
			for (int idx = 0 ; idx < triangles.size() ; idx++) 
			{ 
				TriangleIndices t = triangles[idx]; 
				Vector2D<float> p[3]; 
				
				p[0] = nullBetas[t.getIndex(0)]; 
				p[1] = nullBetas[t.getIndex(1)]; 
				p[2] = nullBetas[t.getIndex(2)]; 
				
				//DrawLine(p[0],p[1]); 
				//DrawLine(p[1],p[2]); 
				//DrawLine(p[2],p[0]); 

				if (p[0].getX() < minx && p[1].getX() < minx && p[2].getX() < minx) 
					continue; 
				if (p[0].getX() > maxx && p[1].getX() > maxx && p[2].getX() > maxx) 
					continue; 
				if (p[0].getY() < miny && p[1].getY() < miny && p[2].getY() < miny) 
					continue; 
				if (p[0].getY() > maxy && p[1].getY() > maxy && p[2].getY() > maxy) 
					continue; 
				
				// Calculate overlapping area 
				Triangle2D<float> nulltriangle(p[0],p[1],p[2]); 
				const Vector2D<float> *hullpoints = bigpoly.getPoints(); 
				int num = bigpoly.getNumberOfPoints() - 2; 
				float totalarea = 0; 
				
				for (int i = 0 ; i < num ; i++) 
				{ 
					Triangle2D<float> polypart(hullpoints[0],hullpoints[1+i],hullpoints[2+i]); 
					
					totalarea += polypart.getOverlapArea(nulltriangle); 
				} 
				
				float scaledoriginalarea = nullTriangleAreas[tp][idx]; 
				float denom = nulltriangle.getArea();
				if (denom == 0.0f)
				{
					const float epsilon = 1e-6;
					// Avoid division by zero, leading to NaN
					denom = epsilon;
					cerr << "WARNING: avoiding division by zero!" << endl;
				}
				float fraction = totalarea/denom;
				float penalty = fraction*scaledoriginalarea; 
				
				//std::cerr << penalty << std::endl; 

				nullenergy += penalty; 
				assert(!isnan(penalty));
			} 

				// TODO: for debugging
	//			const Vector2D<float> *pPolyPoints = bigpoly.getPoints();
	//			int numPoints = bigpoly.getNumberOfPoints();
	//
	//			for (int i = 0 ; i < numPoints ; i++)
	//				DrawLine(pPolyPoints[(i+1)%numPoints],pPolyPoints[i]);
	//			std::cerr << std::endl;
		}
		else if (numimages == 1)
		{
			// Just take the backprojected image and count the number of backprojected
			// triangles that enclose it (as an estimate of the number of images this lens
			// model will generate)

			assert(iface.getNumberOfImagePoints(s, 0) == 1); 
			const Vector2D<float> beta = iface.getBetas(s, 0)[0];

			assert(tp >= 0 && tp < nullTriangles.size());
			int hitCount = getSingleImageNullSpaceHitCount(beta, nullTriangles[tp], nullBetas);

			nullenergy = (float)(hitCount-1);
			if (nullenergy < 0)
				nullenergy = 0;
		}

		assert(weigthIdx < nullWeights.size());

		nullfitness += nullenergy * nullWeights[weigthIdx];
		assert(!isnan(nullenergy));
	}

	assert(!isnan(nullfitness));
	return nullfitness;
}

} // end namespace
