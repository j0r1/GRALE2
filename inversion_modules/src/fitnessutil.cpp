#include "fitnessutil.h"
#include "pointgroupstorage.h"
#include "fitnesscomponent.h"
#include <grale/projectedimagesinterface.h>
#include <grale/polygon2d.h>
#include <grale/triangle2d.h>
#include <assert.h>
#include <limits>
#include <tuple>
#include <complex>

using namespace std;

namespace grale
{

float getScaleFactor_PointImages(const ProjectedImagesInterface &interface,
		                         const std::vector<int> &sourceIndices)
{
	vector<int> sourceGroup(sourceIndices.size(), 0);
	vector<float> scaleFactors(1, 0);

	ScaleFactorWorkspace ws;
	ws.m_floats.resize(4);

	getScaleFactors_PointImages(interface, sourceIndices, sourceGroup, scaleFactors, ws, MinMax);
	float scale = scaleFactors[0];
	assert(scale != 0);
	return scale;
}

void getScaleFactors_PointImages(const ProjectedImagesInterface &interface,
		                         const vector<int> &sourceIndices,
								 const vector<int> &sourceGroup,
								 vector<float> &scaleFactors, // one for each group
								 ScaleFactorWorkspace &ws,
								 PointImageScaleType scaleType)
{
	float defaultScale = (interface.getAngularScale()*interface.getAngularScale())/(ANGLE_ARCSEC*ANGLE_ARCSEC);
	for (auto &s : scaleFactors)
		s = defaultScale;

	assert(scaleFactors.size() > 0);
	assert(sourceGroup.size() == sourceIndices.size());

	if (sourceIndices.size() > 1) // We're assuming that different source indices are present
	{
		if (scaleType == MinMax)
		{
			vector<float> &workSpace = ws.m_floats;

			workSpace.resize(scaleFactors.size()*4);
			for (int i = 0 ; i < scaleFactors.size() ; i++)
			{
				float minx = numeric_limits<float>::max();
				float maxx = -numeric_limits<float>::max();
				float miny = numeric_limits<float>::max();
				float maxy = -numeric_limits<float>::max();

				workSpace[i*4+0] = minx;
				workSpace[i*4+1] = maxx;
				workSpace[i*4+2] = miny;
				workSpace[i*4+3] = maxy;
			}

			for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
			{
				const int s = sourceIndices[sIdx];
				const int group = sourceGroup[sIdx];

				assert(group >= 0 && group < scaleFactors.size());
				assert(s < interface.getNumberOfSources());
				const int numImages = interface.getNumberOfImages(s);

				if (numImages > 1)
				{
					for (int i = 0 ; i < numImages ; i++)
					{
						Vector2D<float> beta = interface.getBetas(s, i)[0];
						float x = beta.getX();
						float y = beta.getY();
						
						float &minx = workSpace[group*4+0];
						float &maxx = workSpace[group*4+1];
						float &miny = workSpace[group*4+2];
						float &maxy = workSpace[group*4+3];

						maxx = MAX(maxx, x);
						minx = MIN(minx, x);
						maxy = MAX(maxy, y);
						miny = MIN(miny, y);
					}
				}
			}

			for (int group = 0 ; group < scaleFactors.size() ; group++)
			{
				float minx = workSpace[group*4+0];
				float maxx = workSpace[group*4+1];
				float miny = workSpace[group*4+2];
				float maxy = workSpace[group*4+3];

				float dx = maxx-minx;
				float dy = maxy-miny;

				if (dx > 0 && dy > 0)
				{
					float scale = 1.0f/(dx*dx + dy*dy);
					assert(!isnan(scale));

					scaleFactors[group] = scale;
				}
			}
		}
		else // MAD
		{
			assert(scaleType == MAD);

			vector<vector<float>> &vecWorkspace = ws.m_vecFloat;
			const int numGroups = scaleFactors.size();

			vecWorkspace.resize(numGroups*2); // one for x components, one for y

			for (int group = 0 ; group < scaleFactors.size() ; group++)
			{
				auto &xPoints = vecWorkspace[group*2+0];
				auto &yPoints = vecWorkspace[group*2+1];
				xPoints.resize(0);
				yPoints.resize(0);
			}

			for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
			{
				const int s = sourceIndices[sIdx];
				const int group = sourceGroup[sIdx];

				assert(group >= 0 && group < scaleFactors.size());
				assert(s < interface.getNumberOfSources());
				const int numImages = interface.getNumberOfImages(s);

				auto &xPoints = vecWorkspace[group*2+0];
				auto &yPoints = vecWorkspace[group*2+1];

				if (numImages > 1)
				{
					for (int i = 0 ; i < numImages ; i++)
					{
						Vector2D<float> beta = interface.getBetas(s, i)[0];
						float x = beta.getX();
						float y = beta.getY();
						
						xPoints.push_back(x);
						yPoints.push_back(y);
					}
				}
			}

			for (int group = 0 ; group < scaleFactors.size() ; group++)
			{
				auto &xPoints = vecWorkspace[group*2+0];
				auto &yPoints = vecWorkspace[group*2+1];
				size_t len = xPoints.size();
				size_t len2 = len/2;

				float xMed = 0, yMed = 0;

				for (int k = 0 ; k < 2 ; k++)
				{
					sort(xPoints.begin(), xPoints.end());
					sort(yPoints.begin(), yPoints.end());

					if (len % 2 == 1) // odd
					{
						xMed = xPoints[len2];
						yMed = yPoints[len2];
					}
					else // even
					{
						xMed = 0.5f*(xPoints[len2] + xPoints[len2-1]);
						yMed = 0.5f*(yPoints[len2] + yPoints[len2-1]);
					}

					if (k == 0)
					{
						for (size_t i = 0 ; i < len ; i++)
						{
							xPoints[i] = std::abs(xPoints[i]-xMed);
							yPoints[i] = std::abs(yPoints[i]-yMed);
						}
					}
				}

				// Actually the square of the scale is used in the other routine
				float scaleSquared = 1.0f/(xMed*xMed + yMed*yMed);
				assert(!isnan(scaleSquared));

				scaleFactors[group] = scaleSquared;
			}
		}
	}
}

float calculateOverlapFitness_PointImages(const ProjectedImagesInterface &interface, 
		                                  const std::vector<int> &sourceIndices,
										  float scale)
{
	vector<int> sourceGroup(sourceIndices.size(), 0);
	vector<float> scaleFactors(1, scale);

	return calculateOverlapFitness_PointImages(interface, sourceIndices, sourceGroup, scaleFactors);
}

float calculateOverlapFitness_PointImages(const ProjectedImagesInterface &interface, 
		                                  const vector<int> &sourceIndices,
										  const vector<int> &sourceGroups,
										  const vector<float> &scaleFactors)
{
	float fitness = 0;
	int sourceCount = 0;

	assert(sourceGroups.size() == sourceIndices.size());
	assert(scaleFactors.size() > 0);

	for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int s = sourceIndices[sIdx];
		const int group = sourceGroups[sIdx];

		assert(group >= 0 && group < scaleFactors.size());
		assert(s < interface.getNumberOfSources());
		const int numImages = interface.getNumberOfImages(s);
		const float scale = scaleFactors[group];

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

					sourceFitness += diff.getLengthSquared();
				}
			}

			int numSprings = (numImages*(numImages-1))/2;
			
			sourceFitness /= (float)numSprings;
			fitness += sourceFitness * scale;
			sourceCount++; // Only count the sources that are actually used
		}
	}

	if (sourceCount > 0)
		fitness /= (float)sourceCount;

	assert(!isnan(fitness));
	return fitness;
}

float calculateOverlapFitness_PointGroups(const PointGroupStorage &pointGroups,
		                                  const ProjectedImagesInterface &interface,
										  const std::vector<int> &sourceIndices,
										  PointGroupRMSType t)
{
	assert(sourceIndices.size() == pointGroups.getNumberOfSources());
	
	int numContribs = 0;
	float fitness = 0;

	vector<Vector2Df> betas;

	for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int src = sourceIndices[sIdx];

		const int groups = pointGroups.getNumberOfGroups(sIdx);
		for (int g = 0 ; g < groups ; g++)
		{
			int numPoints = pointGroups.getNumberOfGroupPoints(sIdx, g);
			if (numPoints <= 1) // can't do anything with only one point
				continue;

			float groupRms = 0;

			auto getGroupPointBetaAndDerivs = [src, sIdx, g, &pointGroups, &interface](int idx,  
					                                       float &axx, float &ayy, float &axy,
														   bool withDerivs = true) -> Vector2Df
			{
				int imgNum = 0, ptNum = 0;
				pointGroups.getGroupPointIndices(sIdx, g, idx, &imgNum, &ptNum);

				assert(imgNum >= 0 && imgNum < interface.getNumberOfImages(src));
				assert(ptNum >= 0 && ptNum < interface.getNumberOfImagePoints(src, imgNum));

				if (withDerivs)
				{
					axx = interface.getDerivativesXX(src, imgNum)[ptNum];
					ayy = interface.getDerivativesYY(src, imgNum)[ptNum];
					axy = interface.getDerivativesXY(src, imgNum)[ptNum];
				}

				return interface.getBetas(src, imgNum)[ptNum];
			};

			auto getGroupPointBeta = [getGroupPointBetaAndDerivs](int idx) -> Vector2Df
			{
				float dummy;
				return getGroupPointBetaAndDerivs(idx, dummy, dummy, dummy, false);
			};

			auto getBetas_All = [numPoints, getGroupPointBeta](vector<Vector2Df> &betas)
			{
				betas.clear();
				for (int p = 0 ; p < numPoints ; p++)
					betas.push_back(getGroupPointBeta(p));
			};

			auto getBetas_Average = [numPoints, getGroupPointBeta](vector<Vector2Df> &betas)
			{
				betas.clear();
				Vector2Df avgBeta { 0, 0 };
				for (int p = 0 ; p < numPoints ; p++)
					avgBeta += getGroupPointBeta(p);
				avgBeta /= numPoints;
				betas.push_back(avgBeta);
			};

			if (t == AllBetas)
				getBetas_All(betas);
			else
			{
				assert(t == AverageBeta);
				getBetas_Average(betas);
			}

			int ptCount = 0;
			for (int p1 = 0 ; p1 < numPoints ; p1++)
			{
				// The point theta0 corresponds to a source plane position beta0
				float axx0, ayy0, axy0;
				Vector2Df beta0 = getGroupPointBetaAndDerivs(p1, axx0, ayy0, axy0);

				float det = (1.0f-axx0)*(1.0f-ayy0)-axy0*axy0;
				const float epsilon = 1e-10; // to avoid division by zero
				if (ABS(det) < epsilon)
					det = copysign(epsilon, det);
				const float invDet = 1.0f/det;

				// inverse of magnification matrix
				const float AI[2][2] = { { invDet*(1.0f-ayy0), invDet*axy0 },
					                     { invDet*axy0, (1.0f-axx0)*invDet } };

				for (auto dBeta : betas)
				{
					ptCount++;

					dBeta -= beta0;

					// The beta corresponding to this point (if close to beta0),
					// would correspond to a theta that differs by this much from
					// theta0
					Vector2Df dTheta { AI[0][0]*dBeta.getX() + AI[0][1]*dBeta.getY(),
						               AI[1][0]*dBeta.getX() + AI[1][1]*dBeta.getY() };

					// TODO: optionally use something else here?
					groupRms += dTheta.getLengthSquared();
				}
			}

			groupRms /= ptCount;

			fitness += groupRms;
			numContribs++;
		}
	}

	if (numContribs > 0)
		fitness = SQRT(fitness/numContribs);
	
	// Let's express this in arcsec, to keep things understandable
	float s = (float)(interface.getAngularScale()/ANGLE_ARCSEC);
	return fitness*s;
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

	float avgImgScale = 0;
	int imgScaleCount = 0;

	assert(pointGroups.getNumberOfSources() == sourceIndices.size());
	assert(rectFlags.size() == sourceIndices.size());
	assert(groupFlags.size() == sourceIndices.size());

	vector<bool> processedSource(sourceIndices.size(), false);

	// In the first iteration, we'll handle the extended images
	// In the second iteration we'll process the point images in a similar
	// manner, but with the source scale based on the extended images
	for (int it = 0 ; it < 2 ; it++)
	{
		for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
		{
			if (processedSource[sIdx])
				continue;

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

			bool pointImage = false;
			for (int i = 0 ; i < numimages ; i++)
			{
				int numpoints = iface.getNumberOfImagePoints(s,i);
				if (numpoints == 1)
					pointImage = true;

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

			}

			scale /= (float)numimages;
			//cerr << "scale for " << s << " is " << scale << endl;

			if (it == 0)
			{
				if (pointImage) // handle point images when the scale has been set
				{
					//cerr << "is point image, skipping for now" << endl;
					continue;
				}

				avgImgScale += scale;
				imgScaleCount++;
			}
			else
			{
				assert(!pointImage); // in the second iteraction we should only process pointimages
				scale = avgImgScale;
				//cerr << "forcing scale to " << scale << endl;

				if (scale == 0)
					return std::numeric_limits<float>::quiet_NaN();
			}

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

						//cerr << "between " << i << " and " << j << ": " << diff1.getLengthSquared() << " "
						//	<< diff2.getLengthSquared() << " "
						//	<< diff3.getLengthSquared() << " "
						//	<< diff4.getLengthSquared() << endl;
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

			//cerr << "sourcefitness: " << sourcefitness << endl;
			posfitness += sourcefitness;
			sourceCount++;

			processedSource[sIdx] = true;

			if (it == 0)
				avgImgScale /= imgScaleCount;

			//std::cerr << "sourcefitness = " << sourcefitness << std::endl;	
		}
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

// For the bayesian version, we assume that Dds/Ds of the extended images
// was set to one; the actual Dds/Ds should be stored in the shear weights
// entries. If set to zero, the actual redshift is unknown and a weighted
// average will be used based on distanceFractionWeights (assumed to be
// normalized)
float calculateWeakLensingFitness_Bayes(const ProjectedImagesInterface &interface, const vector<int> &weakIndices,
								  const vector<pair<float,float>> &unknownDistFracWeightsNormed)
{
	const float epsilon = 1e-6; // to avoid division by zero
	float shearFitness = 0;
	int usedPoints = 0;

	auto calculateEllipticityProbability = [epsilon](float f1, float f2, 
	                                          float axx, float ayy, float axy, float distFrac) -> float
	{
		float gamma1 = 0.5f*(axx-ayy)*distFrac;
		float gamma2 = axy*distFrac;
		float kappa = 0.5f*(axx+ayy)*distFrac;
		
		float oneMinusKappa = 1.0f-kappa;
		if (ABS(oneMinusKappa) < epsilon) // avoid division by zero
				oneMinusKappa = copysign(epsilon, oneMinusKappa);

		float g1 = gamma1/oneMinusKappa;
		float g2 = gamma2/oneMinusKappa;

		// This is based on a Bayesian calculation, assuming the ellipticities
		// are known without error. Further assuming a uniform distribution 
		// for the b/a ratio of the elliptic source shapes, and a uniform 
		// prior on the basis function weights.
		// Then, only the jacobians of the ellipSrc to ellipImg transforms
		// weighted by possibly unknown distance fractions

		// TODO: this fitness can (and will) get negative. The algorithm does
		//       seem to keep working, but for now I'm not really sure if there
		//       are unintended consequences of this.
		
		// weights are ignored currently
		float gSq = g1*g1 + g2*g2;
		float fSq = f1*f1 + f2*f2;
		complex<float> eImg = { f1, f2 };
		complex<float> g = { g1, g2 };
		float F = 0;
		float jacRoot = 0;
		
		// A factor 2 has been omitted, and logprob sign reversed so that we'll be
		// looking for a minimum
		if (gSq <= 1)
		{
			F = abs( (eImg - g)/(1.0f - conj(g)*eImg) );
			jacRoot = (gSq-1.0f)/(fSq*gSq - 2.0f*f1*g1 - 2.0f*f2*g2 + 1.0f);
		}
		else
		{
			F = abs( (1.0f - g*conj(eImg))/(conj(eImg) - conj(g)) );
			jacRoot = (gSq-1.0f)/(fSq + gSq - 2.0f*f1*g2 - 2.0f*f2*g2);
		}

		float jac = jacRoot*jacRoot;
		// Note that this is actualy still a proportionality, depending on the b/a distribution
		// With this proportionality, every b/a, from 0 to 1 is possible with equal probability
		float F1 = F+1.0f;
		float prob = 1.0f/(float(CONST_PI) * (F + epsilon) * F1*F1)*jac; // avoid div by zero
		return prob;
	};

	for (int sIdx = 0 ; sIdx < weakIndices.size() ; sIdx++)
	{
		const int s = weakIndices[sIdx];

		assert(s >= 0 && s < interface.getNumberOfSources());
		int numPoints = interface.getNumberOfImagePoints(s);
		const float *pEll1 = interface.getOriginalShearComponent1s(s);
		const float *pEll2 = interface.getOriginalShearComponent2s(s);
		const float *pDistFrac = interface.getShearWeights(s);
		const float *pAxx = interface.getDerivativesXX(s);
		const float *pAyy = interface.getDerivativesYY(s);
		const float *pAxy = interface.getDerivativesXY(s);

		assert(pEll1 && pEll2 && pAxx && pAyy && pAxy);
		
		for (int i = 0 ; i < numPoints ; i++)
		{
			float axx = pAxx[i];
			float ayy = pAyy[i];
			float axy = pAxy[i];
			float distFrac = pDistFrac[i];
			// we're abusing the weights to store distance fraction, which should be set
			// explicitly (the default value is 1)
			assert(distFrac != 1);

			float f1 = pEll1[i]; // entries represent measured galaxy ellipticities
			float f2 = pEll2[i];

			float elliptProb = 0;
			if (distFrac != 0)
				elliptProb = calculateEllipticityProbability(f1, f2, axx, ayy, axy, distFrac);
			else
			{
				// Unknown redshift/distance fraction, use weighted average according to some distribution
				for (auto fracAndProb : unknownDistFracWeightsNormed)
					elliptProb += calculateEllipticityProbability(f1, f2, axx, ayy, axy, fracAndProb.first) * fracAndProb.second;
			}

			if (elliptProb <= 0) // avoid problems with log
				elliptProb = epsilon; 
			
			float logElliptProb = std::log(elliptProb);				
			shearFitness += -logElliptProb; // use negative to search for a minimum
			usedPoints++;
		}
	}

	if (usedPoints == 0)
		shearFitness = 1e30; // Avoid creating a solution that dominates this fitness because there are no points
	else
		shearFitness /= usedPoints; // this also covers all the weak lensing data sets

	assert(!isnan(shearFitness));
	return shearFitness;
}

float calculateWeakLensingFitness(const ProjectedImagesInterface &interface, const vector<int> &weakIndices,
								  WeakLensingType type, const vector<float> &oneMinusKappaThreshold)
{
	assert(weakIndices.size() == oneMinusKappaThreshold.size());

	const float epsilon = 1e-6; // to avoid division by zero
	float shearFitness = 0;
	int usedPoints = 0;

	for (int sIdx = 0 ; sIdx < weakIndices.size() ; sIdx++)
	{
		const int s = weakIndices[sIdx];

		assert(s >= 0 && s < interface.getNumberOfSources());
		int numPoints = interface.getNumberOfImagePoints(s);
		const float *pStoredShear1 = interface.getOriginalShearComponent1s(s);
		const float *pStoredShear2 = interface.getOriginalShearComponent2s(s);
		const float *pWeights = interface.getShearWeights(s);
		const float *pCalcShear1 = interface.getShearComponents1(s);
		const float *pCalcShear2 = interface.getShearComponents2(s);
		const float *pConvergence = interface.getConvergence(s);

		assert(pStoredShear1 && pStoredShear2 && pCalcShear1 && pCalcShear2 && pConvergence && pWeights);
		float threshold = oneMinusKappaThreshold[sIdx];
		
		for (int i = 0 ; i < numPoints ; i++)
		{
			float gamma1 = pCalcShear1[i]; 
			float gamma2 = pCalcShear2[i];
			float kappa = pConvergence[i];
			float oneMinusKappa = 1.0f-kappa;
			float weight = pWeights[i];

			// Note: we're also using the threshold here in case regular shear is used.
			//       This probably doesn't make much sense but I leave it to the user to
			//       specify that the oneMinusKappaThreshold should be zero in that case.
			if (ABS(oneMinusKappa) >= threshold)
			{
				if (ABS(oneMinusKappa) < epsilon) // avoid division by zero
					oneMinusKappa = copysign(epsilon, oneMinusKappa);

				float d1 = 0, d2 = 0;
				if (type == RealShear) // Real shear, for testing
				{
					// In this case, despite the naming, pStoredShear is supposed to
					// hold the 'normal' shear and not the reduced shear
					d1 = (gamma1-pStoredShear1[i]);
					d2 = (gamma2-pStoredShear2[i]);
				}
				else // Reduced
				{
					float g1 = gamma1/oneMinusKappa;
					float g2 = gamma2/oneMinusKappa;

					if (type == RealReducedShear) // use just reduced shear (for testing)
					{
						d1 = (g1-pStoredShear1[i]);
						d2 = (g2-pStoredShear2[i]);

						shearFitness += (d1*d1 + d2*d2)*weight;
					}
					else if (type == AveragedEllipticities) // Ellipticities
					{
						assert(type == AveragedEllipticities);
						complex<float> Ee { pStoredShear1[i], pStoredShear2[i] }; // entry represents averaged ellipticities
						complex<float> g { g1, g2 };

						// See eg https://arxiv.org/pdf/astro-ph/0509252.pdf
						// E(e) = g if    |g| <= 1
						// E(e) = 1/g* if |g| > 1
						
						if (std::abs(g) > 1.0f)
							g = 1.0f/std::conj(g);

						d1 = (g.real()-Ee.real());
						d2 = (g.imag()-Ee.imag());
						
						shearFitness += (d1*d1 + d2*d2)*weight;
					}
					else // BayesianEllipticities is covered somewhere else, signal error
					{
						return numeric_limits<float>::quiet_NaN();
					}
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

float calculateTimeDelayFitnessPaper2009(const ProjectedImagesInterface &iface, const vector<int> &sourceIndices,
		                        const vector<float> &tdScaleFactors)
{
	assert(tdScaleFactors.size() == 0 || tdScaleFactors.size() == sourceIndices.size());

	float timeDelayFitness = 0;

	for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int s = sourceIndices[sIdx];
		const float scaleFactor = (tdScaleFactors.size() == 0)?0:tdScaleFactors[sIdx];

		assert(iface.getOriginalNumberOfTimeDelays(s) > 0);
		assert(scaleFactor >= 0);

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

						float denom = (scaleFactor == 0)?realTimeDelayDifference:scaleFactor;
						float relativeDiff = (realTimeDelayDifference - calculatedDifference)/denom;

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

float calculateTimeDelayFitnessNoSrc(const ProjectedImagesInterface &iface, const vector<int> &sourceIndices,
                                             const vector<float> &tdScaleFactors)
{
	assert(tdScaleFactors.size() == 0 || tdScaleFactors.size() == sourceIndices.size());

	float timeDelayFitness = 0;

	double D_d = iface.getLensDistance();
	double z_d = iface.getLensRedshift();
	double dFactor = ((D_d*(1.0+z_d)/SPEED_C) * iface.getAngularScale()*iface.getAngularScale())/(60*60*24);
	float baseFactor = (float)dFactor;

	for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int s = sourceIndices[sIdx];
		const float scaleFactor = (tdScaleFactors.size() == 0)?0:tdScaleFactors[sIdx];

		assert(iface.getOriginalNumberOfTimeDelays(s) > 0);
		assert(scaleFactor >= 0);

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

				float denom = (scaleFactor == 0)?realTimeDelayDifference:scaleFactor;
				float relativeDiff = (realTimeDelayDifference - calculatedDifference)/denom;
					
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

// Create one point group for the point images
shared_ptr<ImagesDataExtended> addGroupsToPointImages(const ImagesDataExtended &imgDat, string &errStr)
{
	shared_ptr<ImagesDataExtended> Bad;
	shared_ptr<ImagesDataExtended> grpDat0 { new ImagesDataExtended(imgDat) };
	ImagesData &grpDat = *grpDat0.get();
	bool hasPt = false;

	if (grpDat.getNumberOfGroups() != 0)
	{
		errStr = "Image already contains point groups";
		return Bad;
	}

	int newGrp = grpDat.addGroup();

	for (int img = 0 ; img < grpDat.getNumberOfImages() ; img++)
	{
		int numPts = grpDat.getNumberOfImagePoints(img);
		if (numPts != 1)
		{
			errStr = "Image id " + to_string(img) + " has " + to_string(numPts) + " points (should be exactly 1)";
			return Bad;
		}

		grpDat.addGroupPoint(newGrp, img, 0);
		hasPt = true;
	}

	if (!hasPt)
	{
		errStr = "No points were present, no groups could be added";
		return Bad;
	}
	return grpDat0;
}

} // end namespace
