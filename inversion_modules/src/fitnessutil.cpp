#include "fitnessutil.h"
#include "pointgroupstorage.h"
#include "fitnesscomponent.h"
#include <grale/projectedimagesinterface.h>
#include <grale/polygon2d.h>
#include <grale/triangle2d.h>
#include <assert.h>
#include <limits>
#include <tuple>

using namespace std;

namespace grale
{

float getScaleFactor_PointImages(const ProjectedImagesInterface &interface,
		                         const vector<int> &sourceIndices,
								 const vector<float> &sourceDistanceFractions)
{
	float scale = (interface.getAngularScale()*interface.getAngularScale())/(ANGLE_ARCSEC*ANGLE_ARCSEC);

	assert(sourceIndices.size() == sourceDistanceFractions.size());

	if (sourceIndices.size() > 1) // We're assuming that different source indices are present
	{
		float minx = numeric_limits<float>::max();
		float maxx = -numeric_limits<float>::max();
		float miny = numeric_limits<float>::max();
		float maxy = -numeric_limits<float>::max();

		for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
		{
			const int s = sourceIndices[sIdx];
			const float distFrac = sourceDistanceFractions[sIdx];

			assert(s < interface.getNumberOfSources());
			const int numImages = interface.getNumberOfImages(s);
			const float sqrtDistFrac = SQRT(distFrac);

			if (numImages > 1)
			{
				for (int i = 0 ; i < numImages ; i++)
				{
					Vector2D<float> beta = interface.getBetas(s, i)[0];
					float x = beta.getX()/sqrtDistFrac;
					float y = beta.getY()/sqrtDistFrac;
					
					maxx = MAX(maxx, x);
					minx = MIN(minx, x);
					maxy = MAX(maxy, y);
					miny = MIN(miny, y);
				}
			}
		}

		float dx = maxx-minx;
		float dy = maxy-miny;

		if (dx > 0 && dy > 0)
			scale = 1.0f/(dx*dx + dy*dy);
	}
	assert(!isnan(scale));
	return scale;
}

float calculateOverlapFitness_PointImages(const ProjectedImagesInterface &interface, 
		                                  const vector<int> &sourceIndices,
								          const vector<float> &sourceDistanceFractions,
										  float scale)
{
	float fitness = 0;
	int sourceCount = 0;

	assert(sourceIndices.size() == sourceDistanceFractions.size());

	for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int s = sourceIndices[sIdx];
		const float distFrac = sourceDistanceFractions[sIdx];

		assert(s < interface.getNumberOfSources());
		const int numImages = interface.getNumberOfImages(s);

		if (numImages > 1)
		{
			float sourceFitness = 0;
			
			for (int i = 0 ; i < numImages ; i++)
			{
				Vector2D<float> beta1 = interface.getBetas(s, i)[0];

				for (int j = i+1 ; j < numImages ; j++)
				{
					Vector2D<float> beta2 = interface.getBetas(s, j)[0];
					Vector2D<float> diff = beta2-beta1;

					sourceFitness += diff.getLengthSquared()/distFrac;
				}
			}

			int numSprings = (numImages*(numImages-1))/2;
			
			sourceFitness /= (float)numSprings;
			fitness += sourceFitness;
			sourceCount++; // Only count the sources that are actually used
		}
	}

	fitness *= scale;
	if (sourceCount > 0)
		fitness /= (float)sourceCount;

	assert(!isnan(fitness));
	return fitness;
}

float calculateOverlapFitness_Extended(const PointGroupStorage &pointGroups, const ProjectedImagesInterface &iface,
		                               const vector<int> &sourceIndices,
									   const std::vector<bool> &rectFlags,
									   const std::vector<bool> &groupFlags,
									   FitnessComponentCache *pCache
									   )
{
	float posfitness = 0;
	int sourceCount = 0;

	assert(pointGroups.getNumberOfSources() == sourceIndices.size());
	assert(rectFlags.size() == sourceIndices.size());
	assert(groupFlags.size() == sourceIndices.size());

	for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int s = sourceIndices[sIdx];
		//std::cerr << "source = " << s << endl;

		assert(s < iface.getNumberOfSources());
		int numimages = iface.getNumberOfImages(s);
		if (numimages < 2) // Sources with a single image may be present, they can be used in null space
			continue;

		float sourcefitness = 0;
		vector<Vector2D<float> > bottomleftpoints(numimages);
		vector<Vector2D<float> > toprightpoints(numimages);
		float scale = 0;

		for (int i = 0 ; i < numimages ; i++)
		{
			int numpoints = iface.getNumberOfImagePoints(s,i);
			const Vector2D<float> *betas = iface.getBetas(s,i);

			float maxx = betas[0].getX();
			float minx = betas[0].getX();
			float maxy = betas[0].getY();
			float miny = betas[0].getY();

			//cerr << "x = " << betas[0].getX() << " " << "y = " << betas[0].getY() << endl;
			
			for (int p = 1 ; p < numpoints ; p++)
			{
				maxx = MAX(maxx,betas[p].getX());
				maxy = MAX(maxy,betas[p].getY());
				minx = MIN(minx,betas[p].getX());
				miny = MIN(miny,betas[p].getY());
				
				//cerr << "x = " << betas[p].getX() << " " << "y = " << betas[p].getY() << endl;
			}

			bottomleftpoints[i] = Vector2D<float>(minx,miny);
			toprightpoints[i] = Vector2D<float>(maxx,maxy);
			float imgscale = SQRT((maxx-minx)*(maxx-minx)+(maxy-miny)*(maxy-miny));

			scale += imgscale;

			//std::cerr << "maxx = " << maxx << " maxy = " << maxy << " minx = " << minx << " miny = " << miny << std::endl;
		}

		scale /= (float)numimages;

		// Cache the scale for use by the caustic penalty
		if (pCache)
			pCache->setEstimatedSourceScale(s, scale);

		int numsprings = 0;

		// Surrounding square
		if (rectFlags[sIdx])
		{
			for (int i = 0 ; i < numimages ; i++)
			{
				Vector2D<float> corneri1(bottomleftpoints[i].getX(),toprightpoints[i].getY());
				Vector2D<float> corneri2(toprightpoints[i].getX(),toprightpoints[i].getY());
				Vector2D<float> corneri3(toprightpoints[i].getX(),bottomleftpoints[i].getY());
				Vector2D<float> corneri4(bottomleftpoints[i].getX(),bottomleftpoints[i].getY());
					
				for (int j = i+1 ; j < numimages ; j++)
				{
					Vector2D<float> cornerj1(bottomleftpoints[j].getX(),toprightpoints[j].getY());
					Vector2D<float> cornerj2(toprightpoints[j].getX(),toprightpoints[j].getY());
					Vector2D<float> cornerj3(toprightpoints[j].getX(),bottomleftpoints[j].getY());
					Vector2D<float> cornerj4(bottomleftpoints[j].getX(),bottomleftpoints[j].getY());
					
					Vector2D<float> diff1 = (corneri1-cornerj1)/scale;
					Vector2D<float> diff2 = (corneri2-cornerj2)/scale;
					Vector2D<float> diff3 = (corneri3-cornerj3)/scale;
					Vector2D<float> diff4 = (corneri4-cornerj4)/scale;

					sourcefitness += diff1.getLengthSquared()+diff2.getLengthSquared()+diff3.getLengthSquared()+diff4.getLengthSquared();
				}
			}

			numsprings += 4*(numimages*(numimages-1))/2;
		}
	
		// Point groups (if any)
		if (groupFlags[sIdx])
		{
			assert(sIdx < pointGroups.getNumberOfSources());
			int ng = pointGroups.getNumberOfGroups(sIdx);

			for (int g = 0 ; g < ng ; g++)
			{
				int np = pointGroups.getNumberOfGroupPoints(sIdx,g);
				
				for (int p = 0 ; p < np ; p++)
				{
					int img,point;

					pointGroups.getGroupPointIndices(sIdx,g,p,&img,&point);

					assert(s >= 0 && s < iface.getNumberOfSources());
					assert(img >= 0 && img < iface.getNumberOfImages(s));
					assert(point >= 0 && point < iface.getNumberOfImagePoints(s, img));
					Vector2D<float> point1 = iface.getBetas(s,img)[point];
					
					for (int q = p+1 ; q < np ; q++)
					{
						pointGroups.getGroupPointIndices(sIdx,g,q,&img,&point);
						assert(img >= 0 && img < iface.getNumberOfImages(s));
						assert(point >= 0 && point < iface.getNumberOfImagePoints(s, img));
						Vector2D<float> point2 = iface.getBetas(s,img)[point];

						Vector2D<float> diff = (point2-point1)/scale;

						sourcefitness += diff.getLengthSquared();
					}
				}
				
				numsprings += (np*(np-1))/2;
			}
		}

		if (numsprings != 0)
			sourcefitness /= (float)numsprings;

		posfitness += sourcefitness;
		sourceCount++;

		//std::cerr << "sourcefitness = " << sourcefitness << std::endl;	
	}

	//std::cerr << std::endl;

	if (sourceCount > 0)
		posfitness /= sourceCount;
	
	assert(!isnan(posfitness));
	return posfitness;
}

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

float calculateNullFitness_PointImages(const ProjectedImagesInterface &interface, 
                                       const vector<int> &sourceIndices,
									   const vector<int> &nullIndices,
									   const vector<vector<TriangleIndices> > &nullTriangles,
									   const vector<float> &nullWeights
									   )
{
	assert(sourceIndices.size() == nullIndices.size());
	assert(nullIndices.size() == nullTriangles.size());
	assert(nullIndices.size() == nullWeights.size());

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
			vector<Vector2D<float> > allpoints; // these will contain the backprojected image points
			float maxx = -numeric_limits<float>::max();
			float maxy = -numeric_limits<float>::max();
			float minx = numeric_limits<float>::max();
			float miny = numeric_limits<float>::max();

			// In this case, there are multiple images. We're going to take
			// the convex null of the backprojected images

			for (int i = 0 ; i < numImages ; i++)
			{
				assert(interface.getNumberOfImagePoints(s, i) == 1); 

				const Vector2D<float> *points = interface.getBetas(s, i);
				allpoints.push_back(points[0]);
				maxx = MAX(maxx,points[0].getX());
				maxy = MAX(maxy,points[0].getY());
				minx = MIN(minx,points[0].getX());
				miny = MIN(miny,points[0].getY());
			}

			// If we only have two points, add a few points so that we can use the rest
			// of the code
			if (numImages == 2)
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

			//nullenergy -= (float)numImages;
			if (nullenergy < 0)
				nullenergy = 0;

			// TODO: for debugging
	/*		{
				const Vector2D<float> *pPolyPoints = bigpoly.getPoints();
				int numPoints = bigpoly.getNumberOfPoints();
				
				cerr << bigpoly.getNumberOfPoints() << endl;

				for (int i = 0 ; i < numPoints ; i++)
					DrawLine(pPolyPoints[(i+1)%numPoints],pPolyPoints[i]);
				cerr << endl << endl;
			}
	*/
		}
		else if (numImages == 1)
		{
			const int numPoints = interface.getNumberOfImagePoints(s, 0);
			const Vector2D<float> *pBetas = interface.getBetas(s, 0);
			int hitCount = 0;
			for (int p = 0 ; p < numPoints ; p++)
			{
				auto beta = pBetas[p];

				assert(tp >= 0 && tp < nullTriangles.size());
				hitCount += getSingleImageNullSpaceHitCount(beta, nullTriangles[tp], nullBetas);
			}

			assert(numPoints <= hitCount);
			nullenergy = (float)(hitCount-numPoints);

			if (nullenergy < 0)
				nullenergy = 0;
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
	vector<Polygon2D<float>> polygons(numimages); 

	// first calculate the convex hull for each image
	for (int i = 0 ; i < numimages ; i++) 
	{ 
		int numpoints = iface.getNumberOfImagePoints(s,i); 
		const Vector2D<float> *sourceBetas = iface.getBetas(s,i); 
		
		polygons[i].init(sourceBetas, numpoints, true); 
	} 

	// Then calculate the convex hull of these convex hulls as the
	// current estimate of the source
	Polygon2D<float> bigpoly; 
	vector<Vector2D<float>> allpoints; 

	float maxx = -numeric_limits<float>::max();
	float maxy = -numeric_limits<float>::max();
	float minx = numeric_limits<float>::max();
	float miny = numeric_limits<float>::max();
	
	for (int i = 0 ; i < numimages ; i++) 
	{ 
		int numPolyPoints = polygons[i].getNumberOfPoints(); 
		const Vector2D<float> *points = polygons[i].getPoints(); 
		
		for (int p = 0 ; p < numPolyPoints ; p++) 
		{ 
			allpoints.push_back(points[p]); 
			
			maxx = MAX(maxx,points[p].getX()); 
			maxy = MAX(maxy,points[p].getY()); 
			minx = MIN(minx,points[p].getX()); 
			miny = MIN(miny,points[p].getY()); 
		} 
	} 
	
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
				float fraction = totalarea/nulltriangle.getArea(); 
				float penalty = fraction*scaledoriginalarea; 
				
				//std::cerr << penalty << std::endl; 

				nullenergy += penalty; 
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
	}

	assert(!isnan(nullfitness));
	return nullfitness;
}

float calculateWeakLensingFitness(const ProjectedImagesInterface &interface, const vector<int> &weakIndices,
								  const vector<bool> &reduced, const vector<double> &oneMinusKappaThreshold)
{
	assert(weakIndices.size() == reduced.size());
	assert(weakIndices.size() == oneMinusKappaThreshold.size());

	const float epsilon = 1e-6; // to avoid division by zero
	float shearFitness = 0;
	int usedPoints = 0;

	for (int sIdx = 0 ; sIdx < weakIndices.size() ; sIdx++)
	{
		const int s = weakIndices[sIdx];

		assert(s >= 0 && s < interface.getNumberOfSources());
		int numPoints = interface.getNumberOfImagePoints(s);
		const float *pRealReducedShear1 = interface.getOriginalShearComponent1s(s);
		const float *pRealReducedShear2 = interface.getOriginalShearComponent2s(s);
		const float *pCalcShear1 = interface.getShearComponents1(s);
		const float *pCalcShear2 = interface.getShearComponents2(s);
		const float *pConvergence = interface.getConvergence(s);

		assert(pRealReducedShear1 && pRealReducedShear2 && pCalcShear1 && pCalcShear2 && pConvergence);

		bool useReduced = reduced[sIdx];
		double threshold = oneMinusKappaThreshold[sIdx];
		
		for (int i = 0 ; i < numPoints ; i++)
		{
			float g1 = pCalcShear1[i];
			float g2 = pCalcShear2[i];
			float kappa = pConvergence[i];

			// Note: we're also using the threshold here in case regular shear is used.
			//       This probably doesn't make much sense but I leave it to the user to
			//       specify that the oneMinusKappaThreshold should be zero in that case.
			if (ABS(1.0-kappa) >= threshold)
			{
				if (!useReduced)
				{
					// In this case, despite the naming, pRealReducedShear is supposed to
					// hold the 'normal' shear and not the reduced shear
					float d1 = (g1-pRealReducedShear1[i]);
					float d2 = (g2-pRealReducedShear2[i]);

					shearFitness += d1*d1 + d2*d2;
				}
				else
				{
					// use reduced shear
					float factor = 1.0f-kappa;
					float g1 = pCalcShear1[i];
					float g2 = pCalcShear2[i];

					float d1 = (g1-factor*pRealReducedShear1[i]);
					float d2 = (g2-factor*pRealReducedShear2[i]);

					// This should yield the likelihood for a constant noise on the
					// measured shear parameters. The only deviation is the
					// epsilon: if the error bars would become infinitely large,
					// we replace them by something really large.
					shearFitness += (d1*d1 + d2*d2)/(factor*factor + epsilon);
				}

				usedPoints++;
			}
		}
	}

	if (usedPoints == 0)
		shearFitness = 1e30; // Avoid creating a solution that dominates this fitness because there are no points
	else
		shearFitness /= usedPoints; // this also covers all the weak lensing data sets

	assert(!isnan(shearFitness));
	return shearFitness;
}

float calculateTimeDelayFitness(const ProjectedImagesInterface &iface, const vector<int> &sourceIndices)
{
	float timeDelayFitness = 0;

	for (int s : sourceIndices)
	{
		assert(iface.getOriginalNumberOfTimeDelays(s) > 0);

		int numTimeDelays = iface.getOriginalNumberOfTimeDelays(s);
		float tdfit = 0;
		int count = 0;

		for (int i = 0 ; i < numTimeDelays ; i++)
		{
			int img1, point1;
			float originalDelay1;

			iface.getOriginalTimeDelay(s, i, &img1, &point1, &originalDelay1);

			if (originalDelay1 < -0.0001f) // negative values can be used to simply indicate a source position
				continue;

			for (int j = 0 ; j < numTimeDelays ; j++)
			{
				if (j == i)
					continue;

				int img2, point2;
				float originalDelay2;

				iface.getOriginalTimeDelay(s, j, &img2, &point2, &originalDelay2);

				if (originalDelay2 < -0.0001f)
					continue;

				float realTimeDelayDifference = originalDelay2 - originalDelay1;

				// We have a time delay difference. Now, we'll see how well we can fit this
				// for each backprojected image position

				float squaredSum = 0;

				for (int k = 0 ; k < numTimeDelays ; k++)
				{
					int srcImg, srcPoint;
					float dummyValue;
					Vector2D<float> beta1;

					iface.getOriginalTimeDelay(s, k, &srcImg, &srcPoint, &dummyValue);
					beta1 = iface.getBetas(s, srcImg)[srcPoint];

					for (int l = 0 ; l < numTimeDelays ; l++)
					{
						Vector2D<float> beta2;

						iface.getOriginalTimeDelay(s, l, &srcImg, &srcPoint, &dummyValue);
						beta2 = iface.getBetas(s, srcImg)[srcPoint];

						float delay1 = iface.getTimeDelay(s, img1, point1, beta1);
						float delay2 = iface.getTimeDelay(s, img2, point2, beta2);
						float calculatedDifference = delay2 - delay1;

						float relativeDiff = (realTimeDelayDifference - calculatedDifference)/realTimeDelayDifference;

						squaredSum += relativeDiff*relativeDiff;
					}
				}

				squaredSum /= (float)(numTimeDelays*numTimeDelays);


				tdfit += squaredSum;
				count++;
			}
		}
		
		if (count)
			tdfit /= (int)count;

		timeDelayFitness += tdfit;
	}

	assert(!isnan(timeDelayFitness));
	return timeDelayFitness;
}

float calculateTimeDelayFitnessExperimental(const ProjectedImagesInterface &iface, const vector<int> &sourceIndices)
{
	float timeDelayFitness = 0;

	double D_d = iface.getLensDistance();
	double z_d = iface.getLensRedshift();
	double dFactor = ((D_d*(1.0+z_d)/SPEED_C) * iface.getAngularScale()*iface.getAngularScale())/(60*60*24);
	float baseFactor = (float)dFactor;

	for (int s : sourceIndices)
	{
		assert(iface.getOriginalNumberOfTimeDelays(s) > 0);

		int numTimeDelays = iface.getOriginalNumberOfTimeDelays(s);
		float tdfit = 0;
		int count = 0;
	
		float factor = baseFactor/iface.getDistanceFraction(s);

		for (int i = 0 ; i < numTimeDelays ; i++)
		{
			int img1, point1;
			float originalDelay1;

			iface.getOriginalTimeDelay(s, i, &img1, &point1, &originalDelay1);

			if (originalDelay1 < -0.0001f) // negative values can be used to simply indicate a source position
				continue;

			Vector2Df alpha1 = iface.getAlphas(s, img1)[point1];
			float phi1 = iface.getLensPotential(s, img1)[point1];

			for (int j = i+1 ; j < numTimeDelays ; j++)
			{
				int img2, point2;
				float originalDelay2;

				iface.getOriginalTimeDelay(s, j, &img2, &point2, &originalDelay2);

				if (originalDelay2 < -0.0001f)
					continue;

				float realTimeDelayDifference = originalDelay2 - originalDelay1;

				Vector2Df alpha2 = iface.getAlphas(s, img2)[point2];
				float phi2 = iface.getLensPotential(s, img2)[point2];

				float phiDiff = (phi2-phi1);
				float alphaDiff = 0.5*(alpha2.getLengthSquared() - alpha1.getLengthSquared());

				float calculatedDifference = (alphaDiff - phiDiff)*factor;
				float relativeDiff = (realTimeDelayDifference - calculatedDifference)/realTimeDelayDifference;
					
				tdfit += relativeDiff*relativeDiff;
				count++;
			}
		}
		
		if (count)
			tdfit /= (int)count;

		timeDelayFitness += tdfit;
	}

	assert(!isnan(timeDelayFitness));
	return timeDelayFitness;
}


float calculateTimeDelayFitnessExperimental2(const ProjectedImagesInterface &iface, const vector<int> &sourceIndices)
{
	float timeDelayFitness = 0;

	double D_d = iface.getLensDistance();
	double z_d = iface.getLensRedshift();
	double dFactor = ((D_d*(1.0+z_d)/SPEED_C) * iface.getAngularScale()*iface.getAngularScale())/(60*60*24);
	float baseFactor = (float)dFactor;

	for (int s : sourceIndices)
	{
		assert(iface.getOriginalNumberOfTimeDelays(s) > 0);

		int numTimeDelays = iface.getOriginalNumberOfTimeDelays(s);
		float tdfit = 0;
		int count = 0;
	
		float factor = baseFactor/iface.getDistanceFraction(s);

		for (int i = 0 ; i < numTimeDelays ; i++)
		{
			int img1, point1;
			float originalDelay1;

			iface.getOriginalTimeDelay(s, i, &img1, &point1, &originalDelay1);

			if (originalDelay1 < -0.0001f) // negative values can be used to simply indicate a source position
				continue;

			Vector2Df alpha1 = iface.getAlphas(s, img1)[point1];
			Vector2Df theta1 = iface.getThetas(s, img1)[point1];
			float phi1 = iface.getLensPotential(s, img1)[point1];

			for (int j = i+1 ; j < numTimeDelays ; j++)
			{
				int img2, point2;
				float originalDelay2;

				iface.getOriginalTimeDelay(s, j, &img2, &point2, &originalDelay2);

				if (originalDelay2 < -0.0001f)
					continue;

				float realTimeDelayDifference = originalDelay2 - originalDelay1;

				Vector2Df alpha2 = iface.getAlphas(s, img2)[point2];
				Vector2Df theta2 = iface.getThetas(s, img2)[point2];
				float phi2 = iface.getLensPotential(s, img2)[point2];

				float phiDiff = (phi2-phi1);
				//float alphaDiff = 0.5*(alpha2.getLengthSquared() - alpha1.getLengthSquared());
				Vector2Df alphaSum = alpha2;
				alphaSum += alpha1;

				Vector2Df thetaDiff = theta2;
				thetaDiff -= theta1;

				float calculatedDifference = (
						0.5*(
							(thetaDiff.getX()*alphaSum.getX()) + (thetaDiff.getY()*alphaSum.getY())
						) 
						-
						phiDiff)*factor;

				float relativeDiff = (realTimeDelayDifference - calculatedDifference)/realTimeDelayDifference;
					
				tdfit += relativeDiff*relativeDiff;
				count++;
			}
		}
		
		if (count)
			tdfit /= (int)count;

		timeDelayFitness += tdfit;
	}

	assert(!isnan(timeDelayFitness));
	return timeDelayFitness;
}

void tdHelper_checkRefPoints(const ProjectedImagesInterface &iface, 
							 const std::vector<int> &sourceIndices,
							 std::vector<std::pair<int,int>> &referencePoints)
{
	if (referencePoints.size() != sourceIndices.size()) // We need to calculate the references
	{
		for (int s : sourceIndices)
		{
			assert(iface.getOriginalNumberOfTimeDelays(s) > 0);

			int numTimeDelays = iface.getOriginalNumberOfTimeDelays(s);
			double maxDelay = 0;
			pair<int,int> refPoints;

			for (int i = 0 ; i < numTimeDelays ; i++)
			{
				int img1, point1;
				float originalDelay1;
				iface.getOriginalTimeDelay(s, i, &img1, &point1, &originalDelay1);

				for (int j = i+1 ; j < numTimeDelays ; j++)
				{
					int img2, point2;
					float originalDelay2;
					iface.getOriginalTimeDelay(s, j, &img2, &point2, &originalDelay2);

					double td = std::abs(originalDelay2-originalDelay1);
					if (td > maxDelay)
					{
						maxDelay = td;
						refPoints = { i, j };
					}
				}
			}

			//cerr << "DEBUG: ref " << refPoints.first << "," << refPoints.second << ": " << maxDelay << endl;
			referencePoints.push_back(refPoints);
		}
	}
}

float calculateTimeDelayFitness_Relative(const ProjectedImagesInterface &iface, const vector<int> &sourceIndices,
		                                 vector<pair<int,int>> &referencePoints)
{
	tdHelper_checkRefPoints(iface, sourceIndices, referencePoints);

	float timeDelayFitness = 0;

	auto getAvgCalculatedDifference = [&iface](int s, int img1, int img2, int point1, int point2)
	{
		float avgDelay = 0;
		int count = 0;
		int numTimeDelays = iface.getOriginalNumberOfTimeDelays(s);

		for (int k = 0 ; k < numTimeDelays ; k++)
		{
			int srcImg, srcPoint;
			float dummyValue;
			Vector2D<float> beta1;

			iface.getOriginalTimeDelay(s, k, &srcImg, &srcPoint, &dummyValue);
			beta1 = iface.getBetas(s, srcImg)[srcPoint];

			for (int l = 0 ; l < numTimeDelays ; l++)
			{
				Vector2D<float> beta2;

				iface.getOriginalTimeDelay(s, l, &srcImg, &srcPoint, &dummyValue);
				beta2 = iface.getBetas(s, srcImg)[srcPoint];

				float delay1 = iface.getTimeDelay(s, img1, point1, beta1);
				float delay2 = iface.getTimeDelay(s, img2, point2, beta2);
				float calculatedDifference = delay2 - delay1;

				avgDelay += calculatedDifference;
				count++;
			}
		}
		return avgDelay/count;
	};

	auto getTDScales = [&iface,&getAvgCalculatedDifference](int s, int i1, int i2)
	{
		int img1, img2, point1, point2;
		float originalDelay1, originalDelay2;

		iface.getOriginalTimeDelay(s, i1, &img1, &point1, &originalDelay1);
		iface.getOriginalTimeDelay(s, i2, &img2, &point2, &originalDelay2);
		float scaleReal = std::abs(originalDelay2-originalDelay1);

		float scaleCalc = std::abs(getAvgCalculatedDifference(s, img1, img2, point1, point2));
		return pair<float,float>(scaleReal, scaleCalc);
	};

	for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int s = sourceIndices[sIdx];
		assert(iface.getOriginalNumberOfTimeDelays(s) > 0);

		int numTimeDelays = iface.getOriginalNumberOfTimeDelays(s);
		float tdfit = 0;
		int count = 0;

		int refIdx1, refIdx2;
		tie(refIdx1, refIdx2) = referencePoints[sIdx];
		float tdScaleReal, tdScaleCalc;
		tie(tdScaleReal, tdScaleCalc) = getTDScales(s, refIdx1, refIdx2);
		//cerr << "refIdx1 = " << refIdx1 << " refIdx2 = " << refIdx2 << " tdScaleReal = " << tdScaleReal << " tdScaleCalc = " << tdScaleCalc << endl;

		for (int i = 0 ; i < numTimeDelays ; i++)
		{
			int img1, point1;
			float originalDelay1;

			iface.getOriginalTimeDelay(s, i, &img1, &point1, &originalDelay1);

			if (originalDelay1 < -0.0001f) // negative values can be used to simply indicate a source position
				continue;

			for (int j = 0 ; j < numTimeDelays ; j++)
			{
				if (j == i)
					continue;

				int img2, point2;
				float originalDelay2;

				iface.getOriginalTimeDelay(s, j, &img2, &point2, &originalDelay2);

				if (originalDelay2 < -0.0001f)
					continue;

				float realTimeDelayDifference = originalDelay2 - originalDelay1;
				realTimeDelayDifference /= tdScaleReal;

				// We have a time delay difference. Now, we'll see how well we can fit this
				// for each backprojected image position

				float squaredSum = 0;

				for (int k = 0 ; k < numTimeDelays ; k++)
				{
					int srcImg, srcPoint;
					float dummyValue;
					Vector2D<float> beta1;

					iface.getOriginalTimeDelay(s, k, &srcImg, &srcPoint, &dummyValue);
					beta1 = iface.getBetas(s, srcImg)[srcPoint];

					for (int l = 0 ; l < numTimeDelays ; l++)
					{
						Vector2D<float> beta2;

						iface.getOriginalTimeDelay(s, l, &srcImg, &srcPoint, &dummyValue);
						beta2 = iface.getBetas(s, srcImg)[srcPoint];

						float delay1 = iface.getTimeDelay(s, img1, point1, beta1);
						float delay2 = iface.getTimeDelay(s, img2, point2, beta2);
						float calculatedDifference = delay2 - delay1;
						calculatedDifference /= tdScaleCalc;

						float relativeDiff = (realTimeDelayDifference - calculatedDifference)/realTimeDelayDifference;

						squaredSum += relativeDiff*relativeDiff;
					}
				}

				squaredSum /= (float)(numTimeDelays*numTimeDelays);

				tdfit += squaredSum;
				count++;
			}
		}
		
		if (count)
			tdfit /= (int)count;

		timeDelayFitness += tdfit;
	}

	assert(!isnan(timeDelayFitness));
	return timeDelayFitness;
}

float calculateTimeDelayFitnessExperimental2_Relative(const ProjectedImagesInterface &iface, 
		                                              const std::vector<int> &sourceIndices,
													  std::vector<std::pair<int,int>> &referencePoints)
{
	tdHelper_checkRefPoints(iface, sourceIndices, referencePoints);

	float timeDelayFitness = 0;

	double D_d = iface.getLensDistance();
	double z_d = iface.getLensRedshift();
	double dFactor = ((D_d*(1.0+z_d)/SPEED_C) * iface.getAngularScale()*iface.getAngularScale())/(60*60*24);
	float baseFactor = (float)dFactor;

	auto getCalculatedDifference = [](Vector2Df theta1, Vector2Df theta2, Vector2Df alpha1, 
			                          Vector2Df alpha2, float phi1, float phi2, float factor)
	{
		float phiDiff = (phi2-phi1);
		Vector2Df alphaSum = alpha2;
		alphaSum += alpha1;

		Vector2Df thetaDiff = theta2;
		thetaDiff -= theta1;

		float calculatedDifference = (
				0.5*(
					(thetaDiff.getX()*alphaSum.getX()) + (thetaDiff.getY()*alphaSum.getY())
				) 
				-
				phiDiff)*factor;
		return calculatedDifference;
	};

	auto getTDScales = [&iface,&getCalculatedDifference](int s, int i1, int i2, float factor)
	{
		int img1, img2, point1, point2;
		float originalDelay1, originalDelay2;

		iface.getOriginalTimeDelay(s, i1, &img1, &point1, &originalDelay1);
		iface.getOriginalTimeDelay(s, i2, &img2, &point2, &originalDelay2);
		float scaleReal = std::abs(originalDelay2-originalDelay1);

		Vector2Df alpha1 = iface.getAlphas(s, img1)[point1];
		Vector2Df theta1 = iface.getThetas(s, img1)[point1];
		float phi1 = iface.getLensPotential(s, img1)[point1];

		Vector2Df alpha2 = iface.getAlphas(s, img2)[point2];
		Vector2Df theta2 = iface.getThetas(s, img2)[point2];
		float phi2 = iface.getLensPotential(s, img2)[point2];

		float scaleCalc = std::abs(getCalculatedDifference(theta1, theta2, alpha1, alpha2, phi1, phi2, factor));
		return pair<float,float>(scaleReal, scaleCalc);
	};

	for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++ )
	{
		const int s = sourceIndices[sIdx];
		assert(iface.getOriginalNumberOfTimeDelays(s) > 0);

		int numTimeDelays = iface.getOriginalNumberOfTimeDelays(s);
		float tdfit = 0;
		int count = 0;
	
		float factor = baseFactor/iface.getDistanceFraction(s);

		int refIdx1, refIdx2;
		tie(refIdx1, refIdx2) = referencePoints[sIdx];
		float tdScaleReal, tdScaleCalc;
		tie(tdScaleReal, tdScaleCalc) = getTDScales(s, refIdx1, refIdx2, factor);
		//cerr << "refIdx1 = " << refIdx1 << " refIdx2 = " << refIdx2 << " tdScaleReal = " << tdScaleReal << " tdScaleCalc = " << tdScaleCalc << endl;

		for (int i = 0 ; i < numTimeDelays ; i++)
		{
			int img1, point1;
			float originalDelay1;

			iface.getOriginalTimeDelay(s, i, &img1, &point1, &originalDelay1);

			if (originalDelay1 < -0.0001f) // negative values can be used to simply indicate a source position
				continue;

			Vector2Df alpha1 = iface.getAlphas(s, img1)[point1];
			Vector2Df theta1 = iface.getThetas(s, img1)[point1];
			float phi1 = iface.getLensPotential(s, img1)[point1];

			for (int j = i+1 ; j < numTimeDelays ; j++)
			{
				int img2, point2;
				float originalDelay2;

				iface.getOriginalTimeDelay(s, j, &img2, &point2, &originalDelay2);

				if (originalDelay2 < -0.0001f)
					continue;

				float scaleFactor = 1e5;
				float realTimeDelayDifference = originalDelay2 - originalDelay1;
				realTimeDelayDifference /= tdScaleReal;
				realTimeDelayDifference *= scaleFactor;

				Vector2Df alpha2 = iface.getAlphas(s, img2)[point2];
				Vector2Df theta2 = iface.getThetas(s, img2)[point2];
				float phi2 = iface.getLensPotential(s, img2)[point2];

				float calculatedDifference = getCalculatedDifference(theta1, theta2, alpha1, alpha2, phi1, phi2, factor);
				calculatedDifference /= tdScaleCalc;
				calculatedDifference *= scaleFactor;

				float relativeDiff = (realTimeDelayDifference - calculatedDifference)/realTimeDelayDifference;
					
				tdfit += relativeDiff*relativeDiff;
				count++;
			}
		}
		
		if (count)
			tdfit /= (int)count;

		timeDelayFitness += tdfit;
	}

	assert(!isnan(timeDelayFitness));
	return timeDelayFitness;
}

float calculateKappaThresholdFitness(const ProjectedImagesInterface &iface, const vector<int> &sourceIndices,
		                             const vector<float> &kappaThresholds)
{
	float fitness = 0;

	assert(sourceIndices.size() == kappaThresholds.size());

	for (size_t sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int s = sourceIndices[sIdx];
		const float threshold = kappaThresholds[sIdx];

		assert(s >= 0 && s < iface.getNumberOfSources());
		assert(iface.getNumberOfImages(s) == 1);
		assert(iface.getNumberOfImagePoints(s, 0) > 0);

		const int numPoints = iface.getNumberOfImagePoints(s, 0);
		const float *pKappas = iface.getConvergence(s, 0);

		for (int i = 0 ; i < numPoints ; i++)
		{
			if (pKappas[i] > threshold)
				fitness += (pKappas[i]-threshold);
		}
	}
	assert(!isnan(fitness));
	return fitness;
}

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
